#include "imputed_value.h"
#include <limits>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#if defined(_OPENMP)
#include <omp.h>
#endif

// =============================================================================
// Strict lower-triangular cache for pairwise distances
// =============================================================================
struct StrictLowerTriangularMatrix
{
  size_t n;
  std::vector<double> data;
  std::vector<size_t> row_offsets;

  explicit StrictLowerTriangularMatrix(size_t size)
      : n(size),
        data(size * (size - 1) / 2, 0.0),
        row_offsets(size + 1, 0)
  {
    if (size == 0)
    {
      throw std::invalid_argument("Size must be at least 1");
    }
    // precompute once: row_offsets[r] = r*(r-1)/2. Eliminate calculations in hot path
    for (size_t r = 1; r <= size; ++r)
    {
      row_offsets[r] = row_offsets[r - 1] + (r - 1);
    }
  }

  double &operator()(size_t i, size_t j) noexcept
  {
    return data[row_offsets[i] + j];
  }

  const double &operator()(size_t i, size_t j) const noexcept
  {
    return data[row_offsets[i] + j];
  }
};

// ======================== TEMPLATING README (Part 1) =========================
// We use non-type template parameters so the compiler generates separate
// functions at compile time, eliminating runtime branching in the hot path:
//  - int Method :  0 = Euclidean, 1 = Manhattan
//  - bool Bound :  when true, every GRAIN rows we check whether the partial
//                  distance already exceeds the current worst neighbor.
//                  If so we return Inf immediately, skipping remaining rows.
//                  When false, the original SIMD-friendly loop is emitted.
// =============================================================================

constexpr arma::uword GRAIN = 16; // best for hundreds to thousands of rows

// =============================================================================
// calc_distance_raw — called for Groups 1 & 2 (both target and other have missing)
// =============================================================================
// Bound math (conservative but correct):
//   final_dist >= partial_dist (all terms non-negative)
//   final_n_valid <= partial_n_valid + remaining_rows
//   => final_dist / final_n_valid >= partial_dist / (partial_n_valid + remaining)
//   So if partial_dist > worst_dist * (partial_n_valid + remaining), prune.
// =============================================================================
template <int Method, bool Bound>
inline double calc_distance_raw(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_nmiss,
    const arma::uword n_rows,
    // default is infinity so lambdas can be written more cleanly
    const double worst_dist = std::numeric_limits<double>::infinity())
{
  double dist = 0.0;
  double n_valid = 0.0;

  if constexpr (!Bound)
  {
    // unbounded, full vectorization
#if defined(_OPENMP)
#pragma omp simd reduction(+ : dist, n_valid)
#endif
    for (arma::uword r = 0; r < n_rows; ++r)
    {
      double valid = target_nmiss[r] * other_nmiss[r];
      double diff = target_ptr[r] - other_ptr[r];
      if constexpr (Method == 0)
      {
        dist += valid * diff * diff;
      }
      else
      {
        dist += valid * std::abs(diff);
      }
      n_valid += valid;
    }
  }
  else
  {
    // bounded, precompute chunk count and check bound between chunks
    const arma::uword n_full = n_rows / GRAIN;
    arma::uword r = 0;
    for (arma::uword chunk = 0; chunk < n_full; ++chunk)
    {
#if defined(_OPENMP)
#pragma omp simd reduction(+ : dist, n_valid)
#endif
      for (arma::uword i = 0; i < GRAIN; ++i)
      {
        arma::uword rr = r + i;
        double valid = target_nmiss[rr] * other_nmiss[rr];
        double diff = target_ptr[rr] - other_ptr[rr];
        if constexpr (Method == 0)
        {
          dist += valid * diff * diff;
        }
        else
        {
          dist += valid * std::abs(diff);
        }
        n_valid += valid;
      }
      r += GRAIN;
      // loose bound: partial_dist > worst * (n_valid_so_far + remaining_rows).
      // i.e., we assume all remaining rows are valid
      if (dist > worst_dist * (n_valid + static_cast<double>(n_rows - r)))
      {
        return arma::datum::inf;
      }
    }
    // remainder rows
    for (; r < n_rows; ++r)
    {
      double valid = target_nmiss[r] * other_nmiss[r];
      double diff = target_ptr[r] - other_ptr[r];
      if constexpr (Method == 0)
      {
        dist += valid * diff * diff;
      }
      else
      {
        dist += valid * std::abs(diff);
      }
      n_valid += valid;
    }
  }
  if (n_valid == 0.0)
  {
    return arma::datum::inf;
  }
  return dist / n_valid;
}

template <int Method>
inline double calc_distance(
    const arma::mat &obj,
    const arma::mat &nmiss,
    arma::uword idx1,
    arma::uword idx2)
{
  return calc_distance_raw<Method, false>(
      obj.colptr(idx1), nmiss.colptr(idx1),
      obj.colptr(idx2), nmiss.colptr(idx2),
      obj.n_rows);
}

// =============================================================================
// calc_distance_raw_complete — called for Group 3 (other side fully observed)
// =============================================================================
// Bound math (exact):
//   n_valid is known up front (precomputed target_n_valid).
//   final_dist >= partial_dist, so:
//   if partial_dist > worst_dist * n_valid, then final_dist/n_valid > worst_dist.
// =============================================================================
template <int Method, bool Bound>
inline double calc_distance_raw_complete(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const arma::uword n_rows,
    const double n_valid,
    const double worst_dist = std::numeric_limits<double>::infinity())
{
  double dist = 0.0;

  if constexpr (!Bound)
  {
    // unbounded path
#if defined(_OPENMP)
#pragma omp simd reduction(+ : dist)
#endif
    for (arma::uword r = 0; r < n_rows; ++r)
    {
      double diff = target_ptr[r] - other_ptr[r];
      if constexpr (Method == 0)
      {
        dist += target_nmiss[r] * diff * diff;
      }
      else
      {
        dist += target_nmiss[r] * std::abs(diff);
      }
    }
  }
  else
  {
    // bounded path. Updated bound with new worst_dist
    const double unnorm_bound = worst_dist * n_valid;
    const arma::uword n_full = n_rows / GRAIN;
    arma::uword r = 0;
    for (arma::uword chunk = 0; chunk < n_full; ++chunk)
    {
#if defined(_OPENMP)
#pragma omp simd reduction(+ : dist)
#endif
      for (arma::uword i = 0; i < GRAIN; ++i)
      {
        arma::uword rr = r + i;
        double diff = target_ptr[rr] - other_ptr[rr];
        if constexpr (Method == 0)
        {
          dist += target_nmiss[rr] * diff * diff;
        }
        else
        {
          dist += target_nmiss[rr] * std::abs(diff);
        }
      }
      r += GRAIN;
      if (dist > unnorm_bound)
      {
        return arma::datum::inf;
      }
    }
    for (; r < n_rows; ++r)
    {
      double diff = target_ptr[r] - other_ptr[r];
      if constexpr (Method == 0)
      {
        dist += target_nmiss[r] * diff * diff;
      }
      else
      {
        dist += target_nmiss[r] * std::abs(diff);
      }
    }
  }
  return dist / n_valid;
}

// =============================================================================
// Neighbor tracking
// =============================================================================
struct NeighborInfo
{
  double distance;
  arma::uword index;
  NeighborInfo(double d, arma::uword i) : distance(d), index(i) {}
};

// fill phase: we just need any k neighbors — no ordering yet.
// emplace_back is the absolute fastest way to collect them.
inline void insert_before_k(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
  top_k.emplace_back(dist, idx);
}

// compare the new distance with the worst. This converges when enough neighbor
// candidates are iterated through with the early return
inline void insert_if_better_than_worst(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
  if (dist >= top_k.back().distance)
  {
    return;
  }
  top_k.back() = NeighborInfo(dist, idx);
  for (size_t i = top_k.size() - 1; i > 0 && top_k[i].distance < top_k[i - 1].distance; --i)
  {
    std::swap(top_k[i], top_k[i - 1]);
  }
}

// =============================================================================
// Unified distance_vector. Meat of the algorithm, called once per imputed column
//   Group 1a: grp_impute[0 .. index-1]       (cached)
//   Group 1b: grp_impute[index+1 .. end]     (cached)
//   Group 2:  grp_miss_no_imp                (on-the-fly)
//   Group 3:  grp_complete                   (on-the-fly, optimized kernel)
// =============================================================================
// ======================== TEMPLATING README (Part 2) =========================
// This function is templated on two non-type parameters:
//   1. int Method  -> chooses Euclidean or Manhattan at compile time
//   2. bool UseCache -> decides whether to read from the cache or compute
// =============================================================================
template <int Method, bool UseCache>
std::vector<NeighborInfo> distance_vector(
    const arma::mat &obj,
    const arma::mat &nmiss,
    const arma::uword index,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const arma::uword k,
    const arma::vec &n_valid_vec,
    const StrictLowerTriangularMatrix *cache_ptr = nullptr)
{
  const arma::uword n_rows = obj.n_rows;
  const arma::uword target_col = grp_impute(index);
  const double *target_ptr = obj.colptr(target_col);
  const double *target_nmiss_ptr = nmiss.colptr(target_col);
  const double target_n_valid = n_valid_vec(index);

  std::vector<NeighborInfo> top_k;
  if constexpr (UseCache)
  {
    top_k.reserve(std::max(k, grp_impute.n_elem - 1));
  }
  else
  {
    top_k.reserve(k);
  }

  // ---- fill phase lambdas (Bound=false) ----
  // Group 1: other grp_impute columns (cached when UseCache)
  auto impute_dist_before = [&](arma::uword p) -> double
  {
    if constexpr (UseCache)
    {
      return (*cache_ptr)(index, p);
    }
    else
    {
      return calc_distance_raw<Method, false>(
          target_ptr, target_nmiss_ptr,
          obj.colptr(grp_impute(p)), nmiss.colptr(grp_impute(p)),
          n_rows);
    }
  };
  auto impute_dist_after = [&](arma::uword p) -> double
  {
    if constexpr (UseCache)
    {
      return (*cache_ptr)(p, index);
    }
    else
    {
      return calc_distance_raw<Method, false>(
          target_ptr, target_nmiss_ptr,
          obj.colptr(grp_impute(p)), nmiss.colptr(grp_impute(p)),
          n_rows);
    }
  };
  // Group 2: grp_miss_no_imp (on-the-fly)
  auto miss_no_imp_dist = [&](arma::uword p) -> double
  {
    return calc_distance_raw<Method, false>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_miss_no_imp(p)), nmiss.colptr(grp_miss_no_imp(p)),
        n_rows);
  };
  // Group 3: grp_complete (on-the-fly, optimized kernel)
  auto complete_dist = [&](arma::uword p) -> double
  {
    return calc_distance_raw_complete<Method, false>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_complete(p)),
        n_rows, target_n_valid);
  };

  // ---- replacement phase lambdas (Bound=true) ----
  // Group 1 bounded: only used in !UseCache replacement path
  auto impute_dist_before_bounded = [&](arma::uword p, double worst) -> double
  {
    return calc_distance_raw<Method, true>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_impute(p)), nmiss.colptr(grp_impute(p)),
        n_rows, worst);
  };
  auto impute_dist_after_bounded = [&](arma::uword p, double worst) -> double
  {
    return calc_distance_raw<Method, true>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_impute(p)), nmiss.colptr(grp_impute(p)),
        n_rows, worst);
  };
  // Group 2 bounded
  auto miss_no_imp_dist_bounded = [&](arma::uword p, double worst) -> double
  {
    return calc_distance_raw<Method, true>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_miss_no_imp(p)), nmiss.colptr(grp_miss_no_imp(p)),
        n_rows, worst);
  };
  // Group 3 bounded
  auto complete_dist_bounded = [&](arma::uword p, double worst) -> double
  {
    return calc_distance_raw_complete<Method, true>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_complete(p)),
        n_rows, target_n_valid, worst);
  };

  arma::uword i = 0;
  arma::uword j = index + 1;
  arma::uword m = 0;
  arma::uword c = 0;

  // ==========================================================================
  // Fill phase: collect first k neighbors (Bound=false).
  //  UseCache:  Scan ALL Group 1 (free), keep best k via nth_element,
  //             then top up from Group 3 -> Group 2 if needed.
  //  !UseCache: Group 3 (fastest kernel) -> Group 1 -> Group 2
  // ==========================================================================
  if constexpr (UseCache)
  {
    // scan ALL Group 1 — every distance is a free cache lookup.
    // selecting the best k seeds a much tighter bound for replacement.
    // ---- Group 1a: impute columns before self (cached, free) ----
    for (; i < index; ++i)
    {
      insert_before_k(top_k, impute_dist_before(i), grp_impute(i));
    }
    // ---- Group 1b: impute columns after self (cached, free) ----
    for (; j < grp_impute.n_elem; ++j)
    {
      insert_before_k(top_k, impute_dist_after(j), grp_impute(j));
    }
    // keep only the best k via partial sort
    if (top_k.size() > k)
    {
      std::nth_element(top_k.begin(), top_k.begin() + k, top_k.end(),
                       [](const NeighborInfo &a, const NeighborInfo &b)
                       { return a.distance < b.distance; });
      top_k.erase(top_k.begin() + k, top_k.end());
    }
    // else if Group 1 alone didn't fill k, top up from Group 3 then Group 2
    if (top_k.size() < k)
    {
      // ---- Group 3: complete (fast kernel) ----
      for (; c < grp_complete.n_elem && top_k.size() < k; ++c)
      {
        insert_before_k(top_k, complete_dist(c), grp_complete(c));
      }
      // ---- Group 2: miss_no_imp (most expensive) ----
      for (; m < grp_miss_no_imp.n_elem && top_k.size() < k; ++m)
      {
        insert_before_k(top_k, miss_no_imp_dist(m), grp_miss_no_imp(m));
      }
    }
  }
  else
  {
    // ---- Group 3: complete (fastest kernel) ----
    for (; c < grp_complete.n_elem && top_k.size() < k; ++c)
    {
      insert_before_k(top_k, complete_dist(c), grp_complete(c));
    }
    // ---- Group 1a: impute columns before self ----
    for (; i < index && top_k.size() < k; ++i)
    {
      insert_before_k(top_k, impute_dist_before(i), grp_impute(i));
    }
    // ---- Group 1b: impute columns after self ----
    for (; j < grp_impute.n_elem && top_k.size() < k; ++j)
    {
      insert_before_k(top_k, impute_dist_after(j), grp_impute(j));
    }
    // ---- Group 2: miss_no_imp (most expensive) ----
    for (; m < grp_miss_no_imp.n_elem && top_k.size() < k; ++m)
    {
      insert_before_k(top_k, miss_no_imp_dist(m), grp_miss_no_imp(m));
    }
  }

  // sort once between fill and replacement phases.
  std::sort(top_k.begin(), top_k.end(),
            [](const NeighborInfo &a, const NeighborInfo &b)
            { return a.distance < b.distance; });

  // ==========================================================================
  // Replacement phase — order chosen to tighten the bound fastest.
  //  UseCache:  Group 1 fully exhausted in fill -> Group 3 -> Group 2
  //  !UseCache: Group 3 (fastest, exact bound) -> Group 1 -> Group 2
  // ==========================================================================
  if constexpr (UseCache)
  {
    // Group 1 already fully consumed in fill phase — go straight to Group 3.
    // ---- Group 3: complete (fast kernel, exact bound) ----
    for (; c < grp_complete.n_elem; ++c)
    {
      insert_if_better_than_worst(
          top_k,
          complete_dist_bounded(c, top_k.back().distance),
          grp_complete(c));
    }
    // ---- Group 2: miss_no_imp (most expensive, loosest bound) ----
    for (; m < grp_miss_no_imp.n_elem; ++m)
    {
      insert_if_better_than_worst(
          top_k,
          miss_no_imp_dist_bounded(m, top_k.back().distance),
          grp_miss_no_imp(m));
    }
  }
  else
  {
    // ---- Group 3: complete (fastest kernel, exact bound -> tighten first) ----
    for (; c < grp_complete.n_elem; ++c)
    {
      insert_if_better_than_worst(
          top_k,
          complete_dist_bounded(c, top_k.back().distance),
          grp_complete(c));
    }
    // ---- Group 1a: impute before self ----
    for (; i < index; ++i)
    {
      insert_if_better_than_worst(
          top_k,
          impute_dist_before_bounded(i, top_k.back().distance),
          grp_impute(i));
    }
    // ---- Group 1b: impute after self ----
    for (; j < grp_impute.n_elem; ++j)
    {
      insert_if_better_than_worst(
          top_k,
          impute_dist_after_bounded(j, top_k.back().distance),
          grp_impute(j));
    }
    // ---- Group 2: miss_no_imp (most expensive, loosest bound) ----
    for (; m < grp_miss_no_imp.n_elem; ++m)
    {
      insert_if_better_than_worst(
          top_k,
          miss_no_imp_dist_bounded(m, top_k.back().distance),
          grp_miss_no_imp(m));
    }
  }
  return top_k;
}

// =============================================================================
// Main imputation function
// =============================================================================
template <int Method, bool UseCache>
void impute_knn_brute_impl(
    arma::mat &result,
    const arma::mat &obj,
    const arma::mat &nmiss,
    const arma::uword k,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const arma::uvec &col_offsets,
    const std::vector<arma::uvec> &rows_to_impute_vec,
    const double dist_pow,
    int cores,
    const StrictLowerTriangularMatrix *cache_ptr)
{
  // vectorized n_valid for all impute columns instead of doing in distance_vector
  arma::vec n_valid_vec(grp_impute.n_elem);
  for (arma::uword i = 0; i < grp_impute.n_elem; ++i)
  {
    n_valid_vec(i) = arma::accu(nmiss.col(grp_impute(i)));
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic)
#endif
  for (arma::uword i = 0; i < grp_impute.n_elem; ++i)
  {
    std::vector<NeighborInfo> top_k = distance_vector<Method, UseCache>(
        obj, nmiss, i, grp_impute, grp_miss_no_imp, grp_complete, k,
        n_valid_vec, cache_ptr);

    arma::uword n_neighbors = top_k.size();
    if (n_neighbors == 0)
    {
      continue;
    }

    arma::uvec nn_columns(n_neighbors);
    arma::vec weights(n_neighbors);
    for (arma::uword j = 0; j < n_neighbors; ++j)
    {
      nn_columns(j) = top_k[j].index;
      weights(j) = 1.0 / std::pow(top_k[j].distance + epsilon, dist_pow);
    }

    impute_column_values(
        result, obj, nmiss,
        col_offsets(i), grp_impute(i),
        nn_columns, weights,
        rows_to_impute_vec[i]);
  }
}

// =============================================================================
// Cache setup + dispatch
// =============================================================================
// ======================== TEMPLATING README (Part 3) =========================
// This is the runtime -> compile-time bridge.
//  * The user passes `cache = true/false` (runtime bool)
//  * We decide here whether to build the cache and then call the *templated*
//    version with UseCache = true or false.
//  * Once inside the templated `impute_knn_brute_impl<Method, true/false>`,
//    everything is compile-time again -> zero runtime branching in the hot path.
// =============================================================================
template <int Method>
void dispatch_cache(
    arma::mat &result,
    const arma::mat &obj,
    const arma::mat &nmiss,
    const arma::uword k,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const arma::uvec &col_offsets,
    const std::vector<arma::uvec> &rows_to_impute_vec,
    const double dist_pow,
    const bool use_cache,
    int cores)
{
  if (use_cache)
  {
    StrictLowerTriangularMatrix cache(grp_impute.n_elem);

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic)
#endif
    for (arma::uword row = 1; row < grp_impute.n_elem; ++row)
    {
      for (arma::uword col = 0; col < row; ++col)
      {
        cache(row, col) = calc_distance<Method>(
            obj, nmiss, grp_impute(row), grp_impute(col));
      }
    }

    impute_knn_brute_impl<Method, true>(
        result, obj, nmiss, k, grp_impute, grp_miss_no_imp, grp_complete,
        col_offsets, rows_to_impute_vec, dist_pow, cores, &cache);
  }
  else
  {
    impute_knn_brute_impl<Method, false>(
        result, obj, nmiss, k, grp_impute, grp_miss_no_imp, grp_complete,
        col_offsets, rows_to_impute_vec, dist_pow, cores, nullptr);
  }
}

// [[Rcpp::export]]
arma::mat impute_knn_brute(
    const arma::mat &obj,
    const arma::mat &nmiss,
    const arma::uword k,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const int method,
    const double dist_pow,
    const bool cache = true,
    int cores = 1)
{
  arma::uvec col_offsets;
  std::vector<arma::uvec> rows_to_impute_vec;
  arma::mat result = initialize_result_matrix(nmiss, grp_impute, col_offsets, rows_to_impute_vec);
  if (result.n_rows == 0)
  {
    return result;
  }

  switch (method)
  {
  case 0:
    dispatch_cache<0>(result, obj, nmiss, k, grp_impute, grp_miss_no_imp,
                      grp_complete, col_offsets, rows_to_impute_vec, dist_pow, cache, cores);
    break;
  case 1:
    dispatch_cache<1>(result, obj, nmiss, k, grp_impute, grp_miss_no_imp,
                      grp_complete, col_offsets, rows_to_impute_vec, dist_pow, cache, cores);
    break;
  default:
    throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");
  }

  return result;
}
