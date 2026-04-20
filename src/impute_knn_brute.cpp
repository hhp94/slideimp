#include "imputed_value.h"
#include <limits>
#include <stdexcept>
#include <cmath>
#include <cstring>
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

constexpr arma::uword GRAIN = 32; // best for hundreds to thousands of rows

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

template <int Method, bool UseCache>
std::vector<NeighborInfo> distance_vector(
    const arma::mat &obj_reordered,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword index, // local position in [0, n_imp)
    const arma::uword k,
    const arma::vec &n_valid_vec,
    const StrictLowerTriangularMatrix *cache_ptr = nullptr)
{
    const arma::uword n_rows = obj_reordered.n_rows;
    const arma::uword n_imp = layout.n_imp;
    const arma::uword mni_start = layout.mni_start();
    const arma::uword complete_start = layout.complete_start();

    const double *target_ptr = obj_reordered.colptr(index);
    const double *target_nmiss_ptr = nmiss_masked.colptr(index);
    const double target_n_valid = n_valid_vec(index);

    std::vector<NeighborInfo> top_k;
    if constexpr (UseCache)
    {
        top_k.reserve(std::max(k, n_imp - 1));
    }
    else
    {
        top_k.reserve(k);
    }

    // ---- fill-phase lambdas (Bound = false) ----
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
                obj_reordered.colptr(p), nmiss_masked.colptr(p),
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
                obj_reordered.colptr(p), nmiss_masked.colptr(p),
                n_rows);
        }
    };
    auto miss_no_imp_dist = [&](arma::uword p) -> double
    {
        const arma::uword local = mni_start + p;
        return calc_distance_raw<Method, false>(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(local), nmiss_masked.colptr(local),
            n_rows);
    };
    auto complete_dist = [&](arma::uword p) -> double
    {
        return calc_distance_raw_complete<Method, false>(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(complete_start + p),
            n_rows, target_n_valid);
    };

    // ---- replacement-phase lambdas (Bound = true) ----
    auto impute_dist_before_bounded = [&](arma::uword p, double worst) -> double
    {
        return calc_distance_raw<Method, true>(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(p), nmiss_masked.colptr(p),
            n_rows, worst);
    };
    auto impute_dist_after_bounded = [&](arma::uword p, double worst) -> double
    {
        return calc_distance_raw<Method, true>(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(p), nmiss_masked.colptr(p),
            n_rows, worst);
    };
    auto miss_no_imp_dist_bounded = [&](arma::uword p, double worst) -> double
    {
        const arma::uword local = mni_start + p;
        return calc_distance_raw<Method, true>(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(local), nmiss_masked.colptr(local),
            n_rows, worst);
    };
    auto complete_dist_bounded = [&](arma::uword p, double worst) -> double
    {
        return calc_distance_raw_complete<Method, true>(
            target_ptr, target_nmiss_ptr,
            obj_reordered.colptr(complete_start + p),
            n_rows, target_n_valid, worst);
    };

    arma::uword i = 0;
    arma::uword j = index + 1;
    arma::uword m = 0;
    arma::uword c = 0;

    // ==========================================================================
    // Fill phase
    // ==========================================================================
    if constexpr (UseCache)
    {
        // scan ALL of group 1 (cache hits are free)
        for (; i < index; ++i)
        {
            insert_before_k(top_k, impute_dist_before(i), /*local=*/i);
        }
        for (; j < n_imp; ++j)
        {
            insert_before_k(top_k, impute_dist_after(j), /*local=*/j);
        }
        if (top_k.size() > k)
        {
            std::nth_element(top_k.begin(), top_k.begin() + k, top_k.end(),
                             [](const NeighborInfo &a, const NeighborInfo &b)
                             { return a.distance < b.distance; });
            top_k.erase(top_k.begin() + k, top_k.end());
        }
        if (top_k.size() < k)
        {
            for (; c < layout.n_complete && top_k.size() < k; ++c)
            {
                insert_before_k(top_k, complete_dist(c), complete_start + c);
            }
            for (; m < layout.n_mni && top_k.size() < k; ++m)
            {
                insert_before_k(top_k, miss_no_imp_dist(m), mni_start + m);
            }
        }
    }
    else
    {
        for (; c < layout.n_complete && top_k.size() < k; ++c)
        {
            insert_before_k(top_k, complete_dist(c), complete_start + c);
        }
        for (; i < index && top_k.size() < k; ++i)
        {
            insert_before_k(top_k, impute_dist_before(i), i);
        }
        for (; j < n_imp && top_k.size() < k; ++j)
        {
            insert_before_k(top_k, impute_dist_after(j), j);
        }
        for (; m < layout.n_mni && top_k.size() < k; ++m)
        {
            insert_before_k(top_k, miss_no_imp_dist(m), mni_start + m);
        }
    }

    std::sort(top_k.begin(), top_k.end(),
              [](const NeighborInfo &a, const NeighborInfo &b)
              { return a.distance < b.distance; });

    // ==========================================================================
    // Replacement phase
    // ==========================================================================
    if constexpr (UseCache)
    {
        for (; c < layout.n_complete; ++c)
        {
            insert_if_better_than_worst(
                top_k,
                complete_dist_bounded(c, top_k.back().distance),
                complete_start + c);
        }
        for (; m < layout.n_mni; ++m)
        {
            insert_if_better_than_worst(
                top_k,
                miss_no_imp_dist_bounded(m, top_k.back().distance),
                mni_start + m);
        }
    }
    else
    {
        for (; c < layout.n_complete; ++c)
        {
            insert_if_better_than_worst(
                top_k,
                complete_dist_bounded(c, top_k.back().distance),
                complete_start + c);
        }
        for (; i < index; ++i)
        {
            insert_if_better_than_worst(
                top_k,
                impute_dist_before_bounded(i, top_k.back().distance),
                i);
        }
        for (; j < n_imp; ++j)
        {
            insert_if_better_than_worst(
                top_k,
                impute_dist_after_bounded(j, top_k.back().distance),
                j);
        }
        for (; m < layout.n_mni; ++m)
        {
            insert_if_better_than_worst(
                top_k,
                miss_no_imp_dist_bounded(m, top_k.back().distance),
                mni_start + m);
        }
    }
    return top_k;
}

// =============================================================================
// impute_knn_brute_impl
// =============================================================================
template <int Method, bool UseCache>
void impute_knn_brute_impl(
    arma::mat &result,
    const arma::mat &obj_reordered,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword k,
    const arma::uvec &col_offsets,
    const std::vector<arma::uvec> &rows_to_impute_vec,
    const double dist_pow,
    int cores,
    const StrictLowerTriangularMatrix *cache_ptr)
{
    arma::vec n_valid_vec(layout.n_imp);
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        n_valid_vec(i) = arma::accu(nmiss_masked.col(i));
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic)
#endif
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        std::vector<NeighborInfo> top_k = distance_vector<Method, UseCache>(
            obj_reordered, nmiss_masked, layout, i, k, n_valid_vec, cache_ptr);

        arma::uword n_neighbors = top_k.size();
        if (n_neighbors == 0)
        {
            continue;
        }

        arma::uvec nn_columns(n_neighbors);
        arma::vec weights(n_neighbors);
        for (arma::uword jj = 0; jj < n_neighbors; ++jj)
        {
            nn_columns(jj) = top_k[jj].index; // local position in obj_reordered
            weights(jj) = 1.0 / std::pow(top_k[jj].distance + epsilon, dist_pow);
        }

        impute_column_values(
            result, obj_reordered, nmiss_masked, layout,
            col_offsets(i),
            nn_columns, weights,
            rows_to_impute_vec[i]);
    }
}

// =============================================================================
// dispatch_cache
// =============================================================================
template <int Method>
void dispatch_cache(
    arma::mat &result,
    const arma::mat &obj_reordered,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword k,
    const arma::uvec &col_offsets,
    const std::vector<arma::uvec> &rows_to_impute_vec,
    const double dist_pow,
    const bool use_cache,
    int cores)
{
    if (use_cache)
    {
        StrictLowerTriangularMatrix cache(layout.n_imp);
// Fill cache upfront
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic)
#endif
        for (arma::uword row = 1; row < layout.n_imp; ++row)
        {
            for (arma::uword col = 0; col < row; ++col)
            {
                cache(row, col) = calc_distance<Method>(
                    obj_reordered, nmiss_masked, row, col);
            }
        }

        impute_knn_brute_impl<Method, true>(
            result, obj_reordered, nmiss_masked, layout, k,
            col_offsets, rows_to_impute_vec, dist_pow, cores, &cache);
    }
    else
    {
        impute_knn_brute_impl<Method, false>(
            result, obj_reordered, nmiss_masked, layout, k,
            col_offsets, rows_to_impute_vec, dist_pow, cores, nullptr);
    }
}

// =============================================================================
// Entry point
// -----------------------------------------------------------------------------
// First, we build the reordered working matrices, then dispatche the cache
// Layout of obj_reordered (columns):
//  [ 0 -> n_imp ) grp_impute
//  [ n_imp -> (n_imp + n_mni)) grp_miss_no_imp
//  [ (n_imp + n_mni) -> (n_imp + n_mni + n_complete)) grp_complete
//
// `nmiss_masked` covers only the first two regions (groups that can have NaNs).
// For groups 1 and 2, NaN entries in `obj_reordered` are replaced with 0.0
// (required for SIMD correctness in the kernel) and the corresponding
// nmiss_masked entry is 0.0; all other entries are 1.0.
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
    const bool cache = true,
    int cores = 1)
{
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

    // group 1 and 2: copy and fill in zeros
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        copy_with_mask(i, grp_impute(i));
    }
    // early exit here. Shouldn't happen unless called directly. Note:
    // nmiss_masked is size nrows * (n_imp + n_mni) but the initialize_result_matrix
    // only touches n_imp columns.
    arma::uvec col_offsets;
    std::vector<arma::uvec> rows_to_impute_vec;
    arma::mat result = initialize_result_matrix(
        nmiss_masked, grp_impute, layout, col_offsets, rows_to_impute_vec);
    if (result.n_rows == 0)
    {
        return result;
    }
    // continue filling group 2
    for (arma::uword p = 0; p < layout.n_mni; ++p)
    {
        copy_with_mask(layout.mni_start() + p, grp_miss_no_imp(p));
    }
    // group 3: NaN-free so we can just copy
    for (arma::uword p = 0; p < layout.n_complete; ++p)
    {
        std::memcpy(
            obj_reordered.colptr(layout.complete_start() + p),
            obj.colptr(grp_complete(p)),
            n_rows * sizeof(double));
    }

    switch (method)
    {
    case 0:
        dispatch_cache<0>(result, obj_reordered, nmiss_masked, layout, k,
                          col_offsets, rows_to_impute_vec, dist_pow, cache, cores);
        break;
    case 1:
        dispatch_cache<1>(result, obj_reordered, nmiss_masked, layout, k,
                          col_offsets, rows_to_impute_vec, dist_pow, cache, cores);
        break;
    default:
        throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");
    }

    return result;
}
