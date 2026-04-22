#include "imputed_value.h"
#include "loc_timer.h"
#include <limits>
#include <stdexcept>
#include <cmath>
#include <cstring>
#include <algorithm>
#if defined(_OPENMP)
#include <omp.h>
#endif

constexpr arma::uword GRAIN = 16;

// =============================================================================
// calc_distance_raw — Groups 1 & 2 (both target and other may have missing)
// -----------------------------------------------------------------------------
// Bound math (conservative but correct):
//   final_dist >= partial_dist (all terms non-negative)
//   final_n_valid <= partial_n_valid + remaining_rows
//   => if partial_dist > worst * (partial_n_valid + remaining), prune.
// With worst_dist = +inf (the default), the bound check is always false, so the
// fill-phase callers effectively run the unbounded computation.
// =============================================================================
inline double calc_distance_raw(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_nmiss,
    const arma::uword n_rows,
    const int method,
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
            if (method == 0)
                dist += valid * diff * diff;
            else
                dist += valid * std::abs(diff);
            n_valid += valid;
        }
        r += GRAIN;
        if (dist > worst_dist * (n_valid + static_cast<double>(n_rows - r)))
            return arma::datum::inf;
    }
    for (; r < n_rows; ++r)
    {
        const double valid = target_nmiss[r] * other_nmiss[r];
        const double diff = target_ptr[r] - other_ptr[r];
        if (method == 0)
            dist += valid * diff * diff;
        else
            dist += valid * std::abs(diff);
        n_valid += valid;
    }

    if (n_valid == 0.0)
        return arma::datum::inf;
    return dist / n_valid;
}

// =============================================================================
// calc_distance_raw_complete — Group 3 (other side fully observed)
// -----------------------------------------------------------------------------
// Bound math (exact): n_valid is known up front, so if
//   partial_dist > worst_dist * n_valid
// then final_dist / n_valid > worst_dist, and we can prune.
// =============================================================================
inline double calc_distance_raw_complete(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const arma::uword n_rows,
    const double n_valid,
    const int method,
    const double worst_dist = std::numeric_limits<double>::infinity())
{
    double dist = 0.0;
    const double unnorm_bound = worst_dist * n_valid;
    const arma::uword n_full = n_rows / GRAIN;
    arma::uword r = 0;

    for (arma::uword chunk = 0; chunk < n_full; ++chunk)
    {
        for (arma::uword i = 0; i < GRAIN; ++i)
        {
            const arma::uword rr = r + i;
            const double diff = target_ptr[rr] - other_ptr[rr];
            if (method == 0)
                dist += target_nmiss[rr] * diff * diff;
            else
                dist += target_nmiss[rr] * std::abs(diff);
        }
        r += GRAIN;
        if (dist > unnorm_bound)
            return arma::datum::inf;
    }
    for (; r < n_rows; ++r)
    {
        const double diff = target_ptr[r] - other_ptr[r];
        if (method == 0)
            dist += target_nmiss[r] * diff * diff;
        else
            dist += target_nmiss[r] * std::abs(diff);
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

// Fill phase: we just need any k neighbors — no ordering yet.
inline void insert_before_k(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
    top_k.emplace_back(dist, idx);
}

// Replacement phase: replace worst if new is better, then bubble into place.
inline void insert_if_better_than_worst(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
    if (dist >= top_k.back().distance)
        return;
    top_k.back() = NeighborInfo(dist, idx);
    for (size_t i = top_k.size() - 1; i > 0 && top_k[i].distance < top_k[i - 1].distance; --i)
        std::swap(top_k[i], top_k[i - 1]);
}

// =============================================================================
// distance_vector
// -----------------------------------------------------------------------------
// Two phases:
//  1. Fill: collect up to k candidates with no pruning (worst = +inf default).
//  2. Replacement: keep scanning, pruning against the current worst neighbor.
// Each counter (i, j, m, c) is advanced by the fill phase and resumed by the
// replacement phase, so every column is visited exactly once.
// =============================================================================
std::vector<NeighborInfo> distance_vector(
    const arma::mat &obj_reordered,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword index,
    const arma::uword k,
    const arma::vec &n_valid_vec,
    const int method)
{
    const arma::uword n_rows = obj_reordered.n_rows;
    const arma::uword n_imp = layout.n_imp;
    const arma::uword mni_start = layout.mni_start();
    const arma::uword complete_start = layout.complete_start();

    const double *target_ptr = obj_reordered.colptr(index);
    const double *target_nmiss_ptr = nmiss_masked.colptr(index);
    const double target_n_valid = n_valid_vec(index);

    std::vector<NeighborInfo> top_k;
    top_k.reserve(k);

    auto impute_dist = [&](arma::uword p,
                           double worst = std::numeric_limits<double>::infinity()) -> double
    {
        return calc_distance_raw(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(p), nmiss_masked.colptr(p),
            n_rows, method, worst);
    };
    auto miss_no_imp_dist = [&](arma::uword p,
                                double worst = std::numeric_limits<double>::infinity()) -> double
    {
        const arma::uword local = mni_start + p;
        return calc_distance_raw(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(local), nmiss_masked.colptr(local),
            n_rows, method, worst);
    };
    auto complete_dist = [&](arma::uword p,
                             double worst = std::numeric_limits<double>::infinity()) -> double
    {
        return calc_distance_raw_complete(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(complete_start + p),
            n_rows, target_n_valid, method, worst);
    };

    arma::uword i = 0;
    arma::uword j = index + 1;
    arma::uword m = 0;
    arma::uword c = 0;

    // fill phase ----
    for (; c < layout.n_complete && top_k.size() < k; ++c)
        insert_before_k(top_k, complete_dist(c), complete_start + c);
    for (; i < index && top_k.size() < k; ++i)
        insert_before_k(top_k, impute_dist(i), i);
    for (; j < n_imp && top_k.size() < k; ++j)
        insert_before_k(top_k, impute_dist(j), j);
    for (; m < layout.n_mni && top_k.size() < k; ++m)
        insert_before_k(top_k, miss_no_imp_dist(m), mni_start + m);

    std::sort(top_k.begin(), top_k.end(),
              [](const NeighborInfo &a, const NeighborInfo &b)
              { return a.distance < b.distance; });

    // replacement phase ----
    for (; c < layout.n_complete; ++c)
        insert_if_better_than_worst(
            top_k, complete_dist(c, top_k.back().distance), complete_start + c);
    for (; i < index; ++i)
        insert_if_better_than_worst(
            top_k, impute_dist(i, top_k.back().distance), i);
    for (; j < n_imp; ++j)
        insert_if_better_than_worst(
            top_k, impute_dist(j, top_k.back().distance), j);
    for (; m < layout.n_mni; ++m)
        insert_if_better_than_worst(
            top_k, miss_no_imp_dist(m, top_k.back().distance), mni_start + m);

    return top_k;
}

// =============================================================================
// Entry point
// -----------------------------------------------------------------------------
// Build reordered working matrices then run the imputation directly.
// Layout of obj_reordered (columns):
//  [ 0 -> n_imp )                           grp_impute
//  [ n_imp -> (n_imp + n_mni) )             grp_miss_no_imp
//  [ (n_imp + n_mni) -> n_working )         grp_complete
//
// `nmiss_masked` covers only the first two regions. For groups 1 and 2, NaN
// entries in obj_reordered are replaced with 0.0 (required for correctness of
// the masked kernel) and the corresponding nmiss_masked entry is 0.0;
// everything else is 1.0.
// =============================================================================
//
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
        throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");

    GroupLayout layout{grp_impute.n_elem, grp_miss_no_imp.n_elem, grp_complete.n_elem};
    const arma::uword n_rows = obj.n_rows;

    arma::mat obj_reordered(n_rows, layout.n_working());
    arma::mat nmiss_masked(n_rows, layout.n_imp + layout.n_mni);

    // copy a group 1/2 column: NaN -> 0 in obj_reordered, mask = {0,1} in nmiss_masked.
    auto copy_with_mask = [&](arma::uword local_pos, arma::uword orig_pos)
    {
        const double *src = obj.colptr(orig_pos);
        double *dst = obj_reordered.colptr(local_pos);
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
        copy_with_mask(i, grp_impute(i));

    // early exit: initialize_result_matrix only touches the first n_imp cols
    // of nmiss_masked, which is why group 1 is filled first.
    arma::uvec col_offsets;
    std::vector<arma::uvec> rows_to_impute_vec;
    arma::mat result = initialize_result_matrix(
        nmiss_masked, grp_impute, layout, col_offsets, rows_to_impute_vec);
    if (result.n_rows == 0)
        return result;

    // group 2
    for (arma::uword p = 0; p < layout.n_mni; ++p)
        copy_with_mask(layout.mni_start() + p, grp_miss_no_imp(p));

    // group 3: NaN-free so we can just memcpy.
    for (arma::uword p = 0; p < layout.n_complete; ++p)
        std::memcpy(
            obj_reordered.colptr(layout.complete_start() + p),
            obj.colptr(grp_complete(p)),
            n_rows * sizeof(double));

    arma::vec n_valid_vec(layout.n_imp);
    for (arma::uword i = 0; i < layout.n_imp; ++i)
        n_valid_vec(i) = arma::accu(nmiss_masked.col(i));

    LOC_TIMER_OBJ(knn_tm);
    LOC_TIC(knn_tm, "impute_total");

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic)
#endif
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        LOC_TIMER_SCOPED(knn_tm, "neighbor_search");

        std::vector<NeighborInfo> top_k = distance_vector(
            obj_reordered, nmiss_masked, layout, i, k, n_valid_vec, method);

        const arma::uword n_neighbors = top_k.size();
        if (n_neighbors == 0)
            continue;

        arma::uvec nn_columns(n_neighbors);
        arma::vec weights(n_neighbors);
        for (arma::uword jj = 0; jj < n_neighbors; ++jj)
        {
            nn_columns(jj) = top_k[jj].index;
            weights(jj) = 1.0 / std::pow(top_k[jj].distance + epsilon, dist_pow);
        }

        impute_column_values(
            result, obj_reordered, nmiss_masked, layout,
            col_offsets(i), nn_columns, weights, rows_to_impute_vec[i]);
    }

    LOC_TOC(knn_tm, "impute_total");
    LOC_TIMER_DUMP_RAW(knn_tm, "knn_tm_raw");
    return result;
}
