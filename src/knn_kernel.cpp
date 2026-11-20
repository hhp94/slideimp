// The K-NN kernel: input validation, the masked copy, the neighbor search,
// the weights, the result layout and the accumulation, with the exported
// entry point impute_knn_brute() at the end. Everything but the entry point
// is local to this file.

#include <Rcpp.h>
#include <RcppThread.h> // RcppThread::ProgressBar, checkUserInterrupt
#include "loc_timer.h"
#include "matrix_checks.h" // col_ptr, stop_on_bad_col_index, stop_on_inf
#include "par_for.h"       // parallelFor with a serial path at one thread

#include <algorithm> // std::max, std::min, std::stable_sort
#include <cmath>     // std::abs, std::isfinite, std::isnan, std::pow, std::sqrt
#include <cstddef>   // std::size_t
#include <cstdint>   // std::uint8_t
#include <limits>    // std::numeric_limits
#include <memory>    // std::unique_ptr
#include <utility>   // std::swap
#include <vector>

// Indices that come from R (rows, columns, k) are int, bounded by nrow() and
// ncol(); a column's offset into a buffer is taken in size_t by col_ptr()
// (src/matrix_checks.h).

namespace {

// single source of truth for the mask storage type
using MaskByte = std::uint8_t;

// GroupLayout: single source of truth for group boundaries.
// Groups 1+2 (imputed + missing-no-impute) share identical kernel treatment
// and live together in obj_masked / nmiss_masked at local cols [0, n_masked()).
// n_imp is retained so initialize_result knows which columns actually get
// imputed.
//
// complete_start() is a "virtual" index: there is no physical column at that
// position in obj_masked. It's used as a tag in nn_columns to mark
// "this neighbor came from group 3, offset (local - complete_start()) into
// grp_complete".
//
// The groups are disjoint and in range (validate_knn_inputs), so every count
// here is at most ncol().
struct GroupLayout {
  int n_imp;
  int n_mni;
  int n_complete;
  int mni_start() const { return n_imp; }
  int n_masked() const { return n_imp + n_mni; }
  int complete_start() const { return n_imp + n_mni; }
  int n_working() const { return n_imp + n_mni + n_complete; }
};

void validate_knn_inputs(const Rcpp::NumericMatrix &obj, const int k,
                         const Rcpp::IntegerVector &grp_impute,
                         const Rcpp::IntegerVector &grp_miss_no_imp,
                         const Rcpp::IntegerVector &grp_complete,
                         const int method, const double dist_pow) {
  const int n_rows = obj.nrow();
  const int n_cols = obj.ncol();

  if (n_rows == 0) {
    Rcpp::stop("'obj' must have at least one row");
  }

  if (n_cols < 2) {
    Rcpp::stop("'obj' must have at least two columns");
  }

  // NA_integer_ is INT_MIN, so an NA k fails here too.
  if (k < 1) {
    Rcpp::stop("'k' must be >= 1");
  }

  if (method != 0 && method != 1) {
    Rcpp::stop("'method' must be 0 (Euclidean) or 1 (Manhattan)");
  }

  if (!std::isfinite(dist_pow) || dist_pow < 0.0) {
    Rcpp::stop("'dist_pow' must be finite and >= 0");
  }

  stop_on_bad_col_index(grp_impute, n_cols, "grp_impute");
  stop_on_bad_col_index(grp_miss_no_imp, n_cols, "grp_miss_no_imp");
  stop_on_bad_col_index(grp_complete, n_cols, "grp_complete");

  // The three groups partition the working set. The kernel relies on that
  // when it tags a group-3 neighbor by its offset into grp_complete, and a
  // column listed twice would appear as its own neighbor at distance zero
  // and take the whole weight, so the partition is checked, not assumed.
  std::vector<MaskByte> seen(n_cols, 0);
  auto claim_group = [&](const Rcpp::IntegerVector &g, const char *name) {
    for (R_xlen_t i = 0; i < g.size(); ++i) {
      if (seen[g[i]]) {
        Rcpp::stop("'%s' contains column index %d, which already belongs to "
                   "another group. The three groups must be disjoint",
                   name, g[i]);
      }
      seen[g[i]] = 1;
    }
  };

  claim_group(grp_impute, "grp_impute");
  claim_group(grp_miss_no_imp, "grp_miss_no_imp");
  claim_group(grp_complete, "grp_complete");

  // Group 3 is read straight out of obj with no mask
  // (impute_column_values), because its columns are declared complete. A
  // NaN there is never excluded from the weighted sum, so it produces a
  // missing imputed value with no error. This is where "complete" is
  // enforced; it costs one pass over the group-3 block, against the
  // k * n_complete distance passes the kernel is about to do.
  const double *x = obj.begin();
  for (R_xlen_t i = 0; i < grp_complete.size(); ++i) {
    const double *col = col_ptr(x, grp_complete[i], n_rows);
    for (int r = 0; r < n_rows; ++r) {
      if (std::isnan(col[r])) {
        Rcpp::stop("'grp_complete' column %d has a missing value at row %d. "
                   "'grp_complete' columns must be fully observed",
                   grp_complete[i] + 1, r + 1);
      }
    }
  }

  // disjoint and in range, so at most n_cols
  const int n_working = static_cast<int>(
      grp_impute.size() + grp_miss_no_imp.size() + grp_complete.size());

  if (grp_impute.size() > 0) {
    if (n_working < 2) {
      Rcpp::stop("need at least one candidate neighbor column");
    }

    if (k > n_working - 1) {
      Rcpp::stop("'k' is %d, more than the %d candidate neighbor columns", k,
                 n_working - 1);
    }
  }
}

// -----------------------------------------------------------------------------
// metric definitions, add new metrics by defining a struct with an accumulate()
// that maps a coordinate diff to its contribution to the unnormalized distance.
// -----------------------------------------------------------------------------
struct EuclideanMetric {
  static inline double accumulate(double diff) { return diff * diff; }
};

struct ManhattanMetric {
  static inline double accumulate(double diff) { return std::abs(diff); }
};

// GRAIN size for early exit
constexpr int GRAIN = 16;

constexpr double DIST_INF = std::numeric_limits<double>::infinity();

// -----------------------------------------------------------------------------
// calc_distance_raw - Groups 1 & 2 (both target and other may have missing)
// -----------------------------------------------------------------------------
// Bound math (conservative but correct):
//   final_dist >= partial_dist (all terms non-negative)
//   final_n_valid <= partial_n_valid + remaining_rows
//   => if partial_dist > worst * (partial_n_valid + remaining), prune.
// The check sits under `if constexpr (Bound)`, so the fill-phase callers, which
// instantiate Bound = false, do not run a vacuous comparison - the branch is
// not compiled into them at all. The `worst_dist = +inf` default only covers a
// Bound = true caller that has no worst distance yet.
// -----------------------------------------------------------------------------
template <typename Metric, bool Bound>
inline double
calc_distance_raw(const double *target_ptr, const MaskByte *target_nmiss,
                  const double *other_ptr, const MaskByte *other_nmiss,
                  const int n_rows, const double worst_dist = DIST_INF) {
  double dist = 0.0;
  double n_valid = 0.0;
  const int n_full = n_rows / GRAIN;
  int r = 0;

  for (int chunk = 0; chunk < n_full; ++chunk) {
    for (int i = 0; i < GRAIN; ++i) {
      const int rr = r + i;
      const double valid = target_nmiss[rr] & other_nmiss[rr]; // implicit
      const double diff = target_ptr[rr] - other_ptr[rr];
      dist += valid * Metric::accumulate(diff);
      n_valid += valid;
    }
    r += GRAIN;
    if constexpr (Bound) {
      if (dist > worst_dist * (n_valid + static_cast<double>(n_rows - r))) {
        return DIST_INF;
      }
    }
  }
  for (; r < n_rows; ++r) {
    const double valid = target_nmiss[r] & other_nmiss[r];
    const double diff = target_ptr[r] - other_ptr[r];
    dist += valid * Metric::accumulate(diff);
    n_valid += valid;
  }

  if (n_valid == 0.0) {
    return DIST_INF;
  }

  return dist / n_valid;
}

// -----------------------------------------------------------------------------
// calc_distance_raw_complete - Group 3 (other side fully observed)
// -----------------------------------------------------------------------------
// Bound math (exact): n_valid is known up front, so if
//   partial_dist > worst_dist * n_valid
// then final_dist / n_valid > worst_dist, and we can prune.
// -----------------------------------------------------------------------------
template <typename Metric, bool Bound>
inline double calc_distance_raw_complete(const double *target_ptr,
                                         const MaskByte *target_nmiss,
                                         const double *other_ptr,
                                         const int n_rows, const double n_valid,
                                         const double worst_dist = DIST_INF) {
  double dist = 0.0;
  [[maybe_unused]] double unnorm_bound = 0.0;
  if constexpr (Bound) {
    unnorm_bound = worst_dist * n_valid;
  }

  const int n_full = n_rows / GRAIN;
  int r = 0;

  for (int chunk = 0; chunk < n_full; ++chunk) {
    for (int i = 0; i < GRAIN; ++i) {
      const int rr = r + i;
      const double diff = target_ptr[rr] - other_ptr[rr];
      dist += target_nmiss[rr] * Metric::accumulate(diff);
    }
    r += GRAIN;
    if constexpr (Bound) {
      if (dist > unnorm_bound) {
        return DIST_INF;
      }
    }
  }
  for (; r < n_rows; ++r) {
    const double diff = target_ptr[r] - other_ptr[r];
    dist += target_nmiss[r] * Metric::accumulate(diff);
  }
  // n_valid is the target's own observed count, guaranteed positive by the
  // guard at the top of distance_vector_impl().
  return dist / n_valid;
}

// -----------------------------------------------------------------------------
// Neighbor tracking
// -----------------------------------------------------------------------------
struct NeighborInfo {
  double distance;
  int index;
  NeighborInfo(double d, int i) : distance(d), index(i) {}
};

// fill phase: we just need any k neighbors - no ordering yet.
inline void insert_before_k(std::vector<NeighborInfo> &top_k, double dist,
                            int idx) {
  top_k.emplace_back(dist, idx);
}

// replacement phase: replace worst if new is better, then bubble into place.
inline void insert_if_better_than_worst(std::vector<NeighborInfo> &top_k,
                                        double dist, int idx) {
  if (dist >= top_k.back().distance) {
    return;
  }

  top_k.back() = NeighborInfo(dist, idx);

  for (std::size_t i = top_k.size() - 1;
       i > 0 && top_k[i].distance < top_k[i - 1].distance; --i) {
    std::swap(top_k[i], top_k[i - 1]);
  }
}

// -----------------------------------------------------------------------------
// distance_vector_impl - templated body, fully monomorphized per Metric.
// -----------------------------------------------------------------------------
// Two phases:
//  1. Fill: collect up to k candidates with no pruning (Bound=false).
//  2. Replacement: keep scanning, pruning against the current worst (Bound=true).
// Each counter (p1, p2, c) is advanced by the fill phase and resumed by the
// replacement phase, so every column is visited exactly once.
// A candidate that shares no observed row with the target has distance +inf.
// It is not a neighbor: such entries are dropped from the tail of the sorted
// result before it is returned, so the caller may receive fewer than k.
// Runs in a worker, so it reads plain pointers only.
// -----------------------------------------------------------------------------
template <typename Metric>
std::vector<NeighborInfo>
distance_vector_impl(const double *obj_masked, const MaskByte *nmiss_masked,
                     const int n_rows, const GroupLayout &layout,
                     const int index, const int k,
                     const std::vector<double> &n_valid_vec, const double *obj,
                     const int *grp_complete) {
  const int n_masked = layout.n_masked();
  const int complete_start = layout.complete_start();
  const std::size_t k_size = static_cast<std::size_t>(k);
  const double *target_ptr = col_ptr(obj_masked, index, n_rows);
  const MaskByte *target_nmiss_ptr = col_ptr(nmiss_masked, index, n_rows);
  const double target_n_valid = n_valid_vec[index];

  // A target with no observed row shares no row with any candidate, so by
  // the rule above it has no neighbors at all. calc_distance_raw() reaches
  // that answer on its own, but calc_distance_raw_complete() divides by
  // this count and would return 0/0 = NaN, which is not ordered against
  // anything and would corrupt the selection. Settling it once per target
  // column costs one branch and keeps the two distance functions agreeing.
  if (target_n_valid == 0.0) {
    return {};
  }

  std::vector<NeighborInfo> top_k;
  top_k.reserve(k_size);

  auto masked_dist_unbounded = [&](int p) -> double {
    return calc_distance_raw<Metric, false>(
        target_ptr, target_nmiss_ptr, col_ptr(obj_masked, p, n_rows),
        col_ptr(nmiss_masked, p, n_rows), n_rows);
  };

  auto masked_dist_bounded = [&](int p, double worst) -> double {
    return calc_distance_raw<Metric, true>(
        target_ptr, target_nmiss_ptr, col_ptr(obj_masked, p, n_rows),
        col_ptr(nmiss_masked, p, n_rows), n_rows, worst);
  };

  auto complete_dist_unbounded = [&](int c) -> double {
    return calc_distance_raw_complete<Metric, false>(
        target_ptr, target_nmiss_ptr, col_ptr(obj, grp_complete[c], n_rows),
        n_rows, target_n_valid);
  };

  auto complete_dist_bounded = [&](int c, double worst) -> double {
    return calc_distance_raw_complete<Metric, true>(
        target_ptr, target_nmiss_ptr, col_ptr(obj, grp_complete[c], n_rows),
        n_rows, target_n_valid, worst);
  };

  int p1 = 0;
  int p2 = index + 1;
  int c = 0;

  // fill to k nearest neighbors
  for (; c < layout.n_complete && top_k.size() < k_size; ++c) {
    insert_before_k(top_k, complete_dist_unbounded(c), complete_start + c);
  }
  for (; p1 < index && top_k.size() < k_size; ++p1) {
    insert_before_k(top_k, masked_dist_unbounded(p1), p1);
  }
  for (; p2 < n_masked && top_k.size() < k_size; ++p2) {
    insert_before_k(top_k, masked_dist_unbounded(p2), p2);
  }

  // stable_sort, not sort: the fill is in scan order and the replacement
  // phase evicts top_k.back(), so the order of equal distances decides which
  // tied neighbor is dropped. Keeping it makes the first-encountered
  // candidate win, as the strict < comparisons below already do.
  std::stable_sort(top_k.begin(), top_k.end(),
                   [](const NeighborInfo &a, const NeighborInfo &b) {
                     return a.distance < b.distance;
                   });

  // replacing worst distance
  for (; c < layout.n_complete; ++c) {
    insert_if_better_than_worst(top_k,
                                complete_dist_bounded(c, top_k.back().distance),
                                complete_start + c);
  }
  for (; p1 < index; ++p1) {
    insert_if_better_than_worst(
        top_k, masked_dist_bounded(p1, top_k.back().distance), p1);
  }
  for (; p2 < n_masked; ++p2) {
    insert_if_better_than_worst(
        top_k, masked_dist_bounded(p2, top_k.back().distance), p2);
  }

  // top_k is sorted ascending, so zero-overlap candidates (distance +inf)
  // sit at the tail. Drop them: they carry no information about the target
  // and must not enter the weighted average with a nonzero weight.
  while (!top_k.empty() && !std::isfinite(top_k.back().distance)) {
    top_k.pop_back();
  }

  return top_k;
}

// dispatch once per target column, then the whole body is monomorphized.
std::vector<NeighborInfo>
distance_vector(const double *obj_masked, const MaskByte *nmiss_masked,
                const int n_rows, const GroupLayout &layout, const int index,
                const int k, const std::vector<double> &n_valid_vec,
                const int method, const double *obj, const int *grp_complete) {
  if (method == 0) {
    return distance_vector_impl<EuclideanMetric>(
        obj_masked, nmiss_masked, n_rows, layout, index, k, n_valid_vec, obj,
        grp_complete);
  } else {
    return distance_vector_impl<ManhattanMetric>(
        obj_masked, nmiss_masked, n_rows, layout, index, k, n_valid_vec, obj,
        grp_complete);
  }
}

// inverse-distance KNN weights, rescaled so the largest weight is 1
std::vector<double> knn_weights(const std::vector<NeighborInfo> &top_k,
                                double dist_pow) {
  const std::size_t n = top_k.size();
  std::vector<double> w(n, 1.0);

  if (n == 0 || dist_pow == 0.0) {
    return w; // unweighted KNN
  }

  double min_pos = DIST_INF;
  bool has_zero = false;

  for (std::size_t j = 0; j < n; ++j) {
    const double d = top_k[j].distance;
    if (d <= 0.0) {
      has_zero = true;
    } else if (std::isfinite(d)) {
      min_pos = std::min(min_pos, d);
    }
  }

  // distance_vector() has already dropped infinite distances, so this means
  // every neighbor is an exact duplicate (all distances zero): equal weights.
  if (!std::isfinite(min_pos)) {
    return w;
  }

  if (has_zero) {
    // exact-duplicate neighbors (d == 0) dominate with weight 1;
    // positive-distance neighbors get tiny but nonzero weights.
    const double floor =
        min_pos * std::sqrt(std::numeric_limits<double>::epsilon());
    for (std::size_t j = 0; j < n; ++j) {
      const double d = std::max(top_k[j].distance, floor);
      w[j] = std::pow(floor / d, dist_pow);
    }
  } else {
    for (std::size_t j = 0; j < n; ++j) {
      w[j] = std::pow(min_pos / top_k[j].distance, dist_pow);
    }
  }

  return w;
}

// The result has one entry per missing cell of group 1 (local columns
// [0, n_imp) of the masked copy), column by column and down each column:
// `loc` holds its 1-based row and original 1-based column, `value` holds NA
// until impute_column_values() writes it. col_offsets[i] is where column i's
// entries start, and col_offsets[n_imp] is their total.
// n_col_valid holds the observed count of each masked column, taken during
// the copy-with-mask pass, so the missing counts need no second pass over
// the mask. Only entries [0, layout.n_imp) are read here.
void initialize_result(const std::vector<MaskByte> &nmiss_masked,
                       const int n_rows, const Rcpp::IntegerVector &grp_impute,
                       const GroupLayout &layout,
                       const std::vector<int> &n_col_valid,
                       std::vector<std::size_t> &col_offsets,
                       Rcpp::IntegerMatrix &loc, Rcpp::NumericVector &value) {
  const int n_imp = layout.n_imp;

  col_offsets.clear();
  loc = Rcpp::IntegerMatrix(0, 2);
  value = Rcpp::NumericVector(0);

  if (n_imp == 0) {
    return;
  }

  col_offsets.assign(static_cast<std::size_t>(n_imp) + 1, 0);

  // at most n_imp * n_rows, so a size_t sum cannot overflow
  std::size_t sum_missing = 0;

  for (int i = 0; i < n_imp; ++i) {
    const int n_miss_i = n_rows - n_col_valid[i];
    sum_missing += static_cast<std::size_t>(n_miss_i);
    col_offsets[i + 1] = sum_missing;
  }

  if (sum_missing == 0) {
    return;
  }

  // `loc` is an R matrix, whose dims are int.
  if (sum_missing > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("the columns to impute hold %d missing values, but the result "
               "has a row for each and an R matrix holds at most %d rows. "
               "Impute fewer columns at a time with 'subset'",
               sum_missing, std::numeric_limits<int>::max());
  }

  const int n_total = static_cast<int>(sum_missing);
  loc = Rcpp::IntegerMatrix(Rcpp::no_init_matrix(n_total, 2));
  value = Rcpp::NumericVector(n_total, NA_REAL);
  int *loc_row = loc.begin();
  int *loc_col = loc.begin() + n_total;

  for (int i = 0; i < n_imp; ++i) {
    const MaskByte *mask_col = col_ptr(nmiss_masked.data(), i, n_rows);
    const int orig_col_1based = grp_impute[i] + 1;
    std::size_t w = col_offsets[i];

    // copy_with_mask counted the ones of this mask column into n_col_valid,
    // so it holds exactly as many zeros as its slice of the result has rows.
    for (int r = 0; r < n_rows; ++r) {
      if (mask_col[r] == 0) {
        loc_row[w] = r + 1;           // 1-based row
        loc_col[w] = orig_col_1based; // original 1-based column
        ++w;
      }
    }
  }
}

// Runs in a worker: plain pointers and std containers only, never an Rcpp
// object. `out` is this column's slice of the result values, and loc_rows
// the same slice of the first column of `loc`: the n_miss 1-based rows the
// values stand for.
//
// nn_columns holds tagged positions: values < layout.complete_start() index
// directly into obj_masked/nmiss_masked; values >= layout.complete_start()
// must be de-tagged via (nn_columns[j] - complete_start()) and looked up
// through grp_complete into obj.
//
// Neighbors are split into masked (groups 1+2) and unmasked (group 3) buckets
// before accumulation. Each bucket gets its own clean inner loop.
//
// Unmasked simplification: every row is valid, so wt_ptr[r] += w for every r,
// every unmasked neighbor collapses to a single "total_unmasked_w" added to
// every row once at the end. The ws_ptr accumulation still varies per row
// (different obj values).
void impute_column_values(double *out, const double *obj_masked,
                          const MaskByte *nmiss_masked, const int n_rows,
                          const GroupLayout &layout,
                          const std::vector<int> &nn_columns,
                          const std::vector<double> &nn_weights,
                          const int *loc_rows, const int n_miss,
                          const double *obj, const int *grp_complete) {
  if (n_miss == 0) {
    return;
  }

  std::vector<double> weighted_sums(n_miss, 0.0);
  std::vector<double> weight_totals(n_miss, 0.0);

  double *ws_ptr = weighted_sums.data();
  double *wt_ptr = weight_totals.data();

  const int complete_start = layout.complete_start();

  std::vector<std::size_t> masked_j;
  std::vector<std::size_t> complete_j;
  masked_j.reserve(nn_columns.size());
  complete_j.reserve(nn_columns.size());

  for (std::size_t j = 0; j < nn_columns.size(); ++j) {
    if (nn_columns[j] < complete_start) {
      masked_j.push_back(j);
    } else {
      complete_j.push_back(j);
    }
  }

  // masked neighbors (groups 1+2) ----
  for (std::size_t j : masked_j) {
    const double w = nn_weights[j];
    const int local = nn_columns[j];
    const double *obj_col = col_ptr(obj_masked, local, n_rows);
    const MaskByte *nmiss_col = col_ptr(nmiss_masked, local, n_rows);

    for (int r = 0; r < n_miss; ++r) {
      const int row_idx = loc_rows[r] - 1; // loc holds 1-based rows
      double valid = nmiss_col[row_idx];   // uint8 -> double: exact 0.0/1.0
      double weight = w * valid;
      ws_ptr[r] += weight * obj_col[row_idx];
      wt_ptr[r] += weight;
    }
  }

  // complete neighbors (group 3) ----
  // valid == 1 for every row, so weight == w uniformly across rows.
  double total_complete_w = 0.0;

  for (std::size_t j : complete_j) {
    const double w = nn_weights[j];
    const int c = nn_columns[j] - complete_start;
    const double *obj_col = col_ptr(obj, grp_complete[c], n_rows);

    for (int r = 0; r < n_miss; ++r) {
      ws_ptr[r] += w * obj_col[loc_rows[r] - 1];
    }
    total_complete_w += w;
  }

  if (total_complete_w > 0.0) {
    for (int r = 0; r < n_miss; ++r) {
      wt_ptr[r] += total_complete_w;
    }
  }

  // a row where none of the chosen neighbors is observed keeps the NA
  // initialize_result() wrote.
  for (int r = 0; r < n_miss; ++r) {
    if (wt_ptr[r] > 0.0) {
      out[r] = ws_ptr[r] / wt_ptr[r];
    }
  }
}

} // namespace

// -----------------------------------------------------------------------------
// Entry point
// -----------------------------------------------------------------------------
// Returns list(loc, value): `loc` an integer matrix with a row per missing
// cell of grp_impute (1-based row, 1-based column of obj), `value` the
// imputed value of each, NA where no neighbor is observed. R writes them
// with obj[loc] <- value.
// [[Rcpp::export]]
Rcpp::List impute_knn_brute(const Rcpp::NumericMatrix &obj, const int k,
                            const Rcpp::IntegerVector &grp_impute,
                            const Rcpp::IntegerVector &grp_miss_no_imp,
                            const Rcpp::IntegerVector &grp_complete,
                            const int method, const double dist_pow,
                            int cores = 1, const bool pb = false) {
  validate_knn_inputs(obj, k, grp_impute, grp_miss_no_imp, grp_complete, method,
                      dist_pow);
  stop_on_inf(obj);
  const GroupLayout layout{static_cast<int>(grp_impute.size()),
                           static_cast<int>(grp_miss_no_imp.size()),
                           static_cast<int>(grp_complete.size())};
  const int n_rows = obj.nrow();
  const int n_masked = layout.n_masked();
  const double *obj_ptr = obj.begin();
  const int *grp_complete_ptr = grp_complete.begin();

  const std::size_t n_cells =
      static_cast<std::size_t>(n_rows) * static_cast<std::size_t>(n_masked);
  std::vector<double> obj_masked(n_cells);
  std::vector<MaskByte> nmiss_masked(n_cells);
  std::vector<int> n_col_valid(n_masked, 0);

  auto copy_with_mask = [&](int local_pos, int orig_pos) {
    const double *src = col_ptr(obj_ptr, orig_pos, n_rows);
    double *dst = col_ptr(obj_masked.data(), local_pos, n_rows);
    MaskByte *mask_dst = col_ptr(nmiss_masked.data(), local_pos, n_rows);
    int valid_count = 0;
    for (int r = 0; r < n_rows; ++r) {
      double v = src[r];
      // `finite` is 1 when observed, 0 when NaN. Single isnan call, used
      // directly as both the mask byte and the count increment. The
      // ternary zero out the NaN for kernel math.
      MaskByte finite = !std::isnan(v);
      dst[r] = finite ? v : 0.0;
      mask_dst[r] = finite;
      valid_count += finite;
    }
    n_col_valid[local_pos] = valid_count;
  };

  // group 1
  for (int i = 0; i < layout.n_imp; ++i) {
    copy_with_mask(i, grp_impute[i]);
  }

  // early exit after group 1. Shouldn't really reach this from R.
  std::vector<std::size_t> col_offsets;
  Rcpp::IntegerMatrix loc;
  Rcpp::NumericVector value;
  initialize_result(nmiss_masked, n_rows, grp_impute, layout, n_col_valid,
                    col_offsets, loc, value);

  if (value.size() == 0) {
    return Rcpp::List::create(Rcpp::Named("loc") = loc,
                              Rcpp::Named("value") = value);
  }

  // group 2
  for (int p = 0; p < layout.n_mni; ++p) {
    copy_with_mask(layout.mni_start() + p, grp_miss_no_imp[p]);
  }

  // group 3: read directly from obj via grp_complete - no copy.
  // n_valid_vec built directly from the group-1 portion of n_col_valid
  std::vector<double> n_valid_vec(layout.n_imp);
  for (int i = 0; i < layout.n_imp; ++i) {
    n_valid_vec[i] = static_cast<double>(n_col_valid[i]);
  }

  // The workers get plain pointers and std containers, captured by name, so
  // none can reach an Rcpp object.
  const double *obj_masked_ptr = obj_masked.data();
  const MaskByte *nmiss_masked_ptr = nmiss_masked.data();
  double *value_ptr = value.begin();
  const int *loc_row_ptr = loc.begin();

  // parallelFor setup
  cores = std::max(1, cores);
  const std::size_t n_threads = static_cast<std::size_t>(cores);
  const std::size_t n_batches = static_cast<std::size_t>(cores);
  std::unique_ptr<RcppThread::ProgressBar> bar;
  if (pb) {
    bar = std::make_unique<RcppThread::ProgressBar>(layout.n_imp, 2);
  }
  LOC_TIMER_OBJ(knn_tm);
  LOC_TIC(knn_tm, "impute_total");

  par_for(
      0, static_cast<std::size_t>(layout.n_imp),
      [obj_masked_ptr, nmiss_masked_ptr, n_rows, &layout, k, &n_valid_vec,
       method, obj_ptr, grp_complete_ptr, dist_pow, value_ptr, loc_row_ptr,
       &col_offsets, pb, &bar](std::size_t i) {
        if (i % 5 == 0) {
          RcppThread::checkUserInterrupt();
        }
        // LOC_TIMER_SCOPED(knn_tm, "neighbor_search");
        std::vector<NeighborInfo> top_k =
            distance_vector(obj_masked_ptr, nmiss_masked_ptr, n_rows, layout,
                            static_cast<int>(i), k, n_valid_vec, method,
                            obj_ptr, grp_complete_ptr);
        const std::size_t n_neighbors = top_k.size();
        // With no neighbor (no candidate shares an observed row with this
        // column), its result values keep the NA they were initialized
        // with; the column still counts towards the progress bar.
        if (n_neighbors > 0) {
          std::vector<int> nn_columns(n_neighbors);
          for (std::size_t jj = 0; jj < n_neighbors; ++jj) {
            nn_columns[jj] = top_k[jj].index;
          }
          const std::vector<double> weights = knn_weights(top_k, dist_pow);
          const std::size_t offset = col_offsets[i];
          impute_column_values(value_ptr + offset, obj_masked_ptr,
                               nmiss_masked_ptr, n_rows, layout, nn_columns,
                               weights, loc_row_ptr + offset,
                               static_cast<int>(col_offsets[i + 1] - offset),
                               obj_ptr, grp_complete_ptr);
        }
        if (pb) {
          ++(*bar);
        }
      },
      n_threads, n_batches);

  LOC_TOC(knn_tm, "impute_total");
  LOC_TIMER_DUMP_RAW(knn_tm, "knn_tm_raw");
  return Rcpp::List::create(Rcpp::Named("loc") = loc,
                            Rcpp::Named("value") = value);
}
