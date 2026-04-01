// [[Rcpp::depends(RcppArmadillo)]]

#include "imputed_value.h"
#include <RcppArmadillo.h>
#include <vector>
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

  explicit StrictLowerTriangularMatrix(size_t size)
      : n(size), data(size * (size - 1) / 2, 0.0)
  {
    if (size == 0)
      throw std::invalid_argument("Size must be at least 1");
  }

  double &operator()(size_t i, size_t j) noexcept
  {
    return data[i * (i - 1) / 2 + j];
  }

  const double &operator()(size_t i, size_t j) const noexcept
  {
    return data[i * (i - 1) / 2 + j];
  }
};

// ======================== TEMPLATING README (Part 1) =========================
// We use a non-type template parameter `Method` (0 = Euclidean, 1 = Manhattan)
// so the compiler can generate two completely separate functions at compile time.
// This eliminates the runtime `switch` / `if` that would otherwise sit inside
// the hottest loop (distance_vector). The specializations below are what actually
// call the correct distance routine.

// These are functions that are called in group 1 and 2
template <int Method>
inline double calc_distance_raw(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_nmiss,
    const arma::uword n_rows)
{
  double dist = 0.0;
  double n_valid = 0.0;
#if defined(_OPENMP)
#pragma omp simd reduction(+:dist, n_valid)
#endif
  for (arma::uword r = 0; r < n_rows; ++r)
  {
    double valid = target_nmiss[r] * other_nmiss[r];
    double diff = target_ptr[r] - other_ptr[r];
    if constexpr (Method == 0)
    {
      dist += valid * diff * diff;
    }
    // Add other distances here
    else
    {
      dist += valid * std::abs(diff);
    }
    n_valid += valid;
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
  return calc_distance_raw<Method>(
      obj.colptr(idx1), nmiss.colptr(idx1),
      obj.colptr(idx2), nmiss.colptr(idx2),
      obj.n_rows);
}

// This is an optimized kernel called in group 3 only since all other_nmiss is 1
template <int Method>
inline double calc_distance_raw_complete(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_nmiss,
    const double *__restrict__ other_ptr,
    const arma::uword n_rows,
    const double n_valid)
{
  double dist = 0.0;
#if defined(_OPENMP)
#pragma omp simd reduction(+:dist)
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

inline void insert_before_k(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
  auto it = std::upper_bound(
      top_k.begin(),
      top_k.end(),
      dist,
      [](double d, const NeighborInfo &ni)
      {
        return d < ni.distance;
      });
  top_k.insert(it, NeighborInfo(dist, idx));
}

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
// Unified distance_vector
//   Group 1a: grp_impute[0 .. index-1]       (cached)
//   Group 1b: grp_impute[index+1 .. end]     (cached)
//   Group 2:  grp_miss_no_imp                (on-the-fly)
//   Group 3:  grp_complete                   (on-the-fly, optimized calc_distance_)
// =============================================================================
// ======================== TEMPLATING README (Part 2) =========================
// This function is templated on TWO non-type parameters:
//   1. int Method -> chooses Euclidean or Manhattan at compile time
//   2. bool UseCache -> decides whether to read from the cache or compute on-the-fly
//
// Because both are known at compile time, the compiler can:
//   * inline the correct distance function (no virtual call, no branch)
//   * completely eliminate the `if (UseCache)` test inside the hot lambdas
//   * generate two completely separate code paths (cached vs. non-cached)
//
// The `constexpr if` below is the key: when UseCache = false the whole
// cache-lookup branch is compiled away.
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
    const StrictLowerTriangularMatrix *cache_ptr = nullptr)
{
  const arma::uword n_rows = obj.n_rows;
  const arma::uword target_col = grp_impute(index);
  const double *target_ptr = obj.colptr(target_col);
  const double *target_nmiss_ptr = nmiss.colptr(target_col);
  const double target_n_valid = arma::accu(nmiss.col(target_col));

  std::vector<NeighborInfo> top_k;
  top_k.reserve(k);
  arma::uword remaining = k;

  // Group 1: other grp_impute columns (cached when UseCache)
  auto impute_dist_before = [&](arma::uword p) -> double
  {
    if constexpr (UseCache)
    {
      return (*cache_ptr)(index, p);
    }
    else
    {
      return calc_distance_raw<Method>(
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
      return calc_distance_raw<Method>(
          target_ptr, target_nmiss_ptr,
          obj.colptr(grp_impute(p)), nmiss.colptr(grp_impute(p)),
          n_rows);
    }
  };

  // Group 2: grp_miss_no_imp (on-the-fly)
  auto miss_no_imp_dist = [&](arma::uword p) -> double
  {
    return calc_distance_raw<Method>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_miss_no_imp(p)), nmiss.colptr(grp_miss_no_imp(p)),
        n_rows);
  };

  // Group 3: grp_complete (on-the-fly, optimized calc_distance_)
  auto complete_dist = [&](arma::uword p) -> double
  {
    return calc_distance_raw_complete<Method>(
        target_ptr, target_nmiss_ptr,
        obj.colptr(grp_complete(p)),
        n_rows, target_n_valid);
  };

  // ---- Group 1a: impute columns before self, fill ----
  arma::uword i = 0;
  for (; i < index && remaining > 0; ++i)
  {
    insert_before_k(top_k, impute_dist_before(i), grp_impute(i));
    --remaining;
  }
  // ---- Group 1b: impute columns after self, fill ----
  arma::uword j = index + 1;
  for (; j < grp_impute.n_elem && remaining > 0; ++j)
  {
    insert_before_k(top_k, impute_dist_after(j), grp_impute(j));
    --remaining;
  }
  // ---- Group 2: miss_no_imp, fill ----
  arma::uword m = 0;
  for (; m < grp_miss_no_imp.n_elem && remaining > 0; ++m)
  {
    insert_before_k(top_k, miss_no_imp_dist(m), grp_miss_no_imp(m));
    --remaining;
  }
  // ---- Group 3: complete, fill ----
  arma::uword c = 0;
  for (; c < grp_complete.n_elem && remaining > 0; ++c)
  {
    insert_before_k(top_k, complete_dist(c), grp_complete(c));
    --remaining;
  }
  // ---- Group 1a: impute columns before self, replacement ----
  for (; i < index; ++i)
  {
    insert_if_better_than_worst(top_k, impute_dist_before(i), grp_impute(i));
  }
  // ---- Group 1b: impute columns after self, replacement ----
  for (; j < grp_impute.n_elem; ++j)
  {
    insert_if_better_than_worst(top_k, impute_dist_after(j), grp_impute(j));
  }
  // ---- Group 2: miss_no_imp, replacement ----
  for (; m < grp_miss_no_imp.n_elem; ++m)
  {
    insert_if_better_than_worst(top_k, miss_no_imp_dist(m), grp_miss_no_imp(m));
  }
  // ---- Group 3: complete — replacement ----
  for (; c < grp_complete.n_elem; ++c)
  {
    insert_if_better_than_worst(top_k, complete_dist(c), grp_complete(c));
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
    const double dist_pow,
    int cores,
    const StrictLowerTriangularMatrix *cache_ptr)
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < grp_impute.n_elem; ++i)
  {
    std::vector<NeighborInfo> top_k = distance_vector<Method, UseCache>(
        obj, nmiss, i, grp_impute, grp_miss_no_imp, grp_complete, k, cache_ptr);

    arma::uword n_neighbors = top_k.size();
    if (n_neighbors == 0)
    {
      continue;
    }

    arma::uvec nn_columns(n_neighbors);
    arma::vec nn_dists(n_neighbors);
    for (arma::uword j = 0; j < n_neighbors; ++j)
    {
      nn_columns(j) = top_k[j].index;
      nn_dists(j) = top_k[j].distance;
    }

    arma::vec weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);

    arma::umat nn_columns_mat(n_neighbors, 1);
    nn_columns_mat.col(0) = nn_columns;
    impute_column_values(result, obj, nmiss,
                         col_offsets(i), grp_impute(i),
                         nn_columns_mat, weights);
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
    const double dist_pow,
    const bool use_cache,
    int cores)
{
  if (use_cache)
  {
    StrictLowerTriangularMatrix cache(grp_impute.n_elem);

#ifdef _OPENMP
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
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
        col_offsets, dist_pow, cores, &cache);
  }
  else
  {
    impute_knn_brute_impl<Method, false>(
        result, obj, nmiss, k, grp_impute, grp_miss_no_imp, grp_complete,
        col_offsets, dist_pow, cores, nullptr);
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
  arma::uvec n_miss_impute(grp_impute.n_elem);
  for (arma::uword i = 0; i < grp_impute.n_elem; ++i)
  {
    n_miss_impute(i) = obj.n_rows - arma::accu(nmiss.col(grp_impute(i)));
  }

  arma::uvec col_offsets;
  arma::mat result = initialize_result_matrix(nmiss, grp_impute, n_miss_impute, col_offsets);
  if (result.n_rows == 0)
  {
    return result;
  }

  switch (method)
  {
  case 0:
    dispatch_cache<0>(result, obj, nmiss, k, grp_impute, grp_miss_no_imp,
                      grp_complete, col_offsets, dist_pow, cache, cores);
    break;
  case 1:
    dispatch_cache<1>(result, obj, nmiss, k, grp_impute, grp_miss_no_imp,
                      grp_complete, col_offsets, dist_pow, cache, cores);
    break;
  default:
    throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");
  }

  return result;
}
