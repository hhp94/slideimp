#include "imputed_value.h"
#include "loc_timer.h"
#include <stdexcept>
#include <cmath>
#include <RcppThread.h>
#include <algorithm>

// -----------------------------------------------------------------------------
// metric definitions, add new metrics by defining a struct with an accumulate()
// that maps a coordinate diff to its contribution to the unnormalized distance.
// -----------------------------------------------------------------------------
struct EuclideanMetric
{
    static inline double accumulate(double diff) { return diff * diff; }
};
struct ManhattanMetric
{
    static inline double accumulate(double diff) { return std::abs(diff); }
};

// GRAIN size for early exit
constexpr arma::uword GRAIN = 16;

// -----------------------------------------------------------------------------
// calc_distance_raw — Groups 1 & 2 (both target and other may have missing)
// -----------------------------------------------------------------------------
// Bound math (conservative but correct):
//   final_dist >= partial_dist (all terms non-negative)
//   final_n_valid <= partial_n_valid + remaining_rows
//   => if partial_dist > worst * (partial_n_valid + remaining), prune.
// With worst_dist = +inf (the default), the bound check is always false, so the
// fill-phase callers effectively run the unbounded computation.
// -----------------------------------------------------------------------------
template <typename Metric, bool Bound>
inline double calc_distance_raw(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_nmiss,
    const arma::uword n_rows,
    const double worst_dist = std::numeric_limits<double>::infinity())
{
    double dist = 0.0;
    double n_valid = 0.0;
    const arma::uword n_full = n_rows / GRAIN;
    arma::uword r = 0;

    for (arma::uword chunk = 0; chunk < n_full; ++chunk)
    {
        for (arma::uword i = 0; i < GRAIN; ++i)
        {
            const arma::uword rr = r + i;
            const double valid = target_nmiss[rr] * other_nmiss[rr];
            const double diff = target_ptr[rr] - other_ptr[rr];
            dist += valid * Metric::accumulate(diff);
            n_valid += valid;
        }
        r += GRAIN;
        if constexpr (Bound)
        {
            if (dist > worst_dist * (n_valid + static_cast<double>(n_rows - r)))
            {
                return arma::datum::inf;
            }
        }
    }
    for (; r < n_rows; ++r)
    {
        const double valid = target_nmiss[r] * other_nmiss[r];
        const double diff = target_ptr[r] - other_ptr[r];
        dist += valid * Metric::accumulate(diff);
        n_valid += valid;
    }

    if (n_valid == 0.0)
    {
        return arma::datum::inf;
    }

    return dist / n_valid;
}

// =============================================================================
// calc_distance_raw_complete — Group 3 (other side fully observed)
// -----------------------------------------------------------------------------
// Bound math (exact): n_valid is known up front, so if
//   partial_dist > worst_dist * n_valid
// then final_dist / n_valid > worst_dist, and we can prune.
// =============================================================================
template <typename Metric, bool Bound>
inline double calc_distance_raw_complete(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const arma::uword n_rows,
    const double n_valid,
    const double worst_dist = std::numeric_limits<double>::infinity())
{
    double dist = 0.0;
    [[maybe_unused]] double unnorm_bound = 0.0;
    if constexpr (Bound)
    {
        unnorm_bound = worst_dist * n_valid;
    }

    const arma::uword n_full = n_rows / GRAIN;
    arma::uword r = 0;

    for (arma::uword chunk = 0; chunk < n_full; ++chunk)
    {
        for (arma::uword i = 0; i < GRAIN; ++i)
        {
            const arma::uword rr = r + i;
            const double diff = target_ptr[rr] - other_ptr[rr];
            dist += target_nmiss[rr] * Metric::accumulate(diff);
        }
        r += GRAIN;
        if constexpr (Bound)
        {
            if (dist > unnorm_bound)
            {
                return arma::datum::inf;
            }
        }
    }
    for (; r < n_rows; ++r)
    {
        const double diff = target_ptr[r] - other_ptr[r];
        dist += target_nmiss[r] * Metric::accumulate(diff);
    }
    return dist / n_valid;
}

// -----------------------------------------------------------------------------
// Neighbor tracking
// -----------------------------------------------------------------------------
struct NeighborInfo
{
    double distance;
    arma::uword index;
    NeighborInfo(double d, arma::uword i) : distance(d), index(i) {}
};

// fill phase: we just need any k neighbors — no ordering yet.
inline void insert_before_k(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
    top_k.emplace_back(dist, idx);
}

// replacement phase: replace worst if new is better, then bubble into place.
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

// -----------------------------------------------------------------------------
// distance_vector_impl — templated body, fully monomorphized per Metric.
// -----------------------------------------------------------------------------
// Two phases:
//  1. Fill: collect up to k candidates with no pruning (Bound=false).
//  2. Replacement: keep scanning, pruning against the current worst (Bound=true).
// Each counter (p1, p2, c) is advanced by the fill phase and resumed by the
// replacement phase, so every column is visited exactly once.
// -----------------------------------------------------------------------------
template <typename Metric>
std::vector<NeighborInfo> distance_vector_impl(
    const arma::mat &obj_masked,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword index,
    const arma::uword k,
    const arma::vec &n_valid_vec,
    const arma::mat &obj,
    const arma::uvec &grp_complete)
{
    const arma::uword n_rows = obj_masked.n_rows;
    const arma::uword n_masked = layout.n_masked();
    const arma::uword complete_start = layout.complete_start();
    const double *target_ptr = obj_masked.colptr(index);
    const double *target_nmiss_ptr = nmiss_masked.colptr(index);
    const double target_n_valid = n_valid_vec(index);

    std::vector<NeighborInfo> top_k;
    top_k.reserve(k);

    auto masked_dist_unbounded = [&](arma::uword p) -> double
    {
        return calc_distance_raw<Metric, false>(
            target_ptr, target_nmiss_ptr,
            obj_masked.colptr(p), nmiss_masked.colptr(p),
            n_rows);
    };

    auto masked_dist_bounded = [&](arma::uword p, double worst) -> double
    {
        return calc_distance_raw<Metric, true>(
            target_ptr, target_nmiss_ptr,
            obj_masked.colptr(p), nmiss_masked.colptr(p),
            n_rows, worst);
    };

    auto complete_dist_unbounded = [&](arma::uword c) -> double
    {
        return calc_distance_raw_complete<Metric, false>(
            target_ptr, target_nmiss_ptr,
            obj.colptr(grp_complete(c)),
            n_rows, target_n_valid);
    };

    auto complete_dist_bounded = [&](arma::uword c, double worst) -> double
    {
        return calc_distance_raw_complete<Metric, true>(
            target_ptr, target_nmiss_ptr,
            obj.colptr(grp_complete(c)),
            n_rows, target_n_valid, worst);
    };

    arma::uword p1 = 0;
    arma::uword p2 = index + 1;
    arma::uword c = 0;

    // fill to k nearest neighbors
    for (; c < layout.n_complete && top_k.size() < k; ++c)
    {
        insert_before_k(top_k, complete_dist_unbounded(c), complete_start + c);
    }
    for (; p1 < index && top_k.size() < k; ++p1)
    {
        insert_before_k(top_k, masked_dist_unbounded(p1), p1);
    }
    for (; p2 < n_masked && top_k.size() < k; ++p2)
    {
        insert_before_k(top_k, masked_dist_unbounded(p2), p2);
    }

    std::sort(top_k.begin(), top_k.end(),
              [](const NeighborInfo &a, const NeighborInfo &b)
              { return a.distance < b.distance; });

    // replacing worst distance
    for (; c < layout.n_complete; ++c)
    {
        insert_if_better_than_worst(
            top_k, complete_dist_bounded(c, top_k.back().distance), complete_start + c);
    }
    for (; p1 < index; ++p1)
    {
        insert_if_better_than_worst(
            top_k, masked_dist_bounded(p1, top_k.back().distance), p1);
    }
    for (; p2 < n_masked; ++p2)
    {
        insert_if_better_than_worst(
            top_k, masked_dist_bounded(p2, top_k.back().distance), p2);
    }

    return top_k;
}

// dispatch once per target column, then the whole body is monomorphized.
std::vector<NeighborInfo> distance_vector(
    const arma::mat &obj_masked,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword index,
    const arma::uword k,
    const arma::vec &n_valid_vec,
    const int method,
    const arma::mat &obj,
    const arma::uvec &grp_complete)
{
    if (method == 0)
    {
        return distance_vector_impl<EuclideanMetric>(
            obj_masked, nmiss_masked, layout, index, k, n_valid_vec,
            obj, grp_complete);
    }
    else
    {
        return distance_vector_impl<ManhattanMetric>(
            obj_masked, nmiss_masked, layout, index, k, n_valid_vec,
            obj, grp_complete);
    }
}

// -----------------------------------------------------------------------------
// Entry point
// -----------------------------------------------------------------------------
// Build reordered working matrices then run the imputation directly.
// Layout of obj_masked (columns):
//  [ 0 -> n_imp )                           grp_impute
//  [ n_imp -> (n_imp + n_mni) )             grp_miss_no_imp
//
// `nmiss_masked` covers only the first two regions. For groups 1 and 2, NaN
// entries in obj_masked are replaced with 0.0 (required for correctness of
// the masked kernel) and the corresponding `nmiss_masked` entry is 0.0;
// everything else is 1.0.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat impute_knn_brute(
    const arma::mat &obj,
    const arma::uword k,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const int method,
    const double dist_pow,
    int cores = 1)
{
    if (method != 0 && method != 1)
    {
        throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");
    }

    GroupLayout layout{grp_impute.n_elem, grp_miss_no_imp.n_elem, grp_complete.n_elem};
    const arma::uword n_rows = obj.n_rows;
    const arma::uword n_masked = layout.n_masked();

    arma::mat obj_masked(n_rows, n_masked);
    arma::mat nmiss_masked(n_rows, n_masked);

    auto copy_with_mask = [&](arma::uword local_pos, arma::uword orig_pos)
    {
        const double *src = obj.colptr(orig_pos);
        double *dst = obj_masked.colptr(local_pos);
        double *mask_dst = nmiss_masked.colptr(local_pos);
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            double v = src[r];
            bool is_nan = std::isnan(v);
            dst[r] = is_nan ? 0.0 : v;
            mask_dst[r] = is_nan ? 0.0 : 1.0;
        }
    };

    // group 1
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        copy_with_mask(i, grp_impute(i));
    }

    // early exit after group 1. Shouldn't really reach this from R.
    arma::uvec col_offsets;
    std::vector<arma::uvec> rows_to_impute_vec;
    arma::mat result = initialize_result_matrix(
        nmiss_masked, grp_impute, layout, col_offsets, rows_to_impute_vec);

    if (result.n_rows == 0)
    {
        return result;
    }

    // group 2
    for (arma::uword p = 0; p < layout.n_mni; ++p)
    {
        copy_with_mask(layout.mni_start() + p, grp_miss_no_imp(p));
    }

    // group 3: read directly from obj via grp_complete — no copy.
    arma::vec n_valid_vec(layout.n_imp);
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        n_valid_vec(i) = arma::accu(nmiss_masked.col(i));
    }

    // parallelFor setup
    cores = std::max(1, cores);
    const size_t n_threads = static_cast<size_t>(cores);
    const size_t n_batches = static_cast<size_t>(cores);

    LOC_TIMER_OBJ(knn_tm);
    LOC_TIC(knn_tm, "impute_total");

    RcppThread::parallelFor(
        0,
        layout.n_imp,
        [&](arma::uword i)
        {
            // LOC_TIMER_SCOPED(knn_tm, "neighbor_search");
            std::vector<NeighborInfo> top_k = distance_vector(
                obj_masked, nmiss_masked, layout, i, k, n_valid_vec, method,
                obj, grp_complete);
            const arma::uword n_neighbors = top_k.size();
            if (n_neighbors == 0)
            {
                return; // was `continue`
            }
            arma::uvec nn_columns(n_neighbors);
            arma::vec weights(n_neighbors);
            for (arma::uword jj = 0; jj < n_neighbors; ++jj)
            {
                nn_columns(jj) = top_k[jj].index;
                weights(jj) = 1.0 / std::pow(top_k[jj].distance + epsilon, dist_pow);
            }
            impute_column_values(
                result, obj_masked, nmiss_masked, layout,
                col_offsets(i), nn_columns, weights, rows_to_impute_vec[i],
                obj, grp_complete);
        },
        n_threads, n_batches);

    LOC_TOC(knn_tm, "impute_total");
    LOC_TIMER_DUMP_RAW(knn_tm, "knn_tm_raw");
    return result;
}
