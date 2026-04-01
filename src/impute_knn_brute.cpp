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

// ============================================================================
// Strict lower-triangular cache for pairwise distances
// ============================================================================
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

  size_t total_bytes() const noexcept
  {
    return data.size() * sizeof(double);
  }
};

// ============================================================================
// Distance functions
// ============================================================================

// Euclidean (method 0): called from distance_vector where target is hoisted outside
inline double calc_distance_euclidean(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_miss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_miss,
    const arma::uword n_rows)
{
  double dist = 0.0;
  double n_valid = 0.0;
  for (arma::uword r = 0; r < n_rows; ++r)
  {
    double valid = (1.0 - target_miss[r]) * (1.0 - other_miss[r]);
    double diff = target_ptr[r] - other_ptr[r];
    dist += valid * diff * diff;
    n_valid += valid;
  }
  if (static_cast<int>(n_valid) == 0) return arma::datum::inf;
  return dist / n_valid;
}

// Manhattan (method 1)
inline double calc_distance_manhattan(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_miss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_miss,
    const arma::uword n_rows)
{
  double dist = 0.0;
  double n_valid = 0.0;
  for (arma::uword r = 0; r < n_rows; ++r)
  {
    double valid = (1.0 - target_miss[r]) * (1.0 - other_miss[r]);
    double diff = target_ptr[r] - other_ptr[r];
    dist += valid * std::abs(diff);
    n_valid += valid;
  }
  if (static_cast<int>(n_valid) == 0) return arma::datum::inf;
  return dist / n_valid;
}

// ============================================================================
// Armadillo-index versions for cache-fill (both sides looked up here)
// ============================================================================
template <int Method>
inline double calc_distance(
    const arma::mat &obj,
    const arma::mat &miss,
    arma::uword idx1,
    arma::uword idx2);

template <>
inline double calc_distance<0>(
    const arma::mat &obj,
    const arma::mat &miss,
    arma::uword idx1,
    arma::uword idx2)
{
  return calc_distance_euclidean(
      obj.colptr(idx1), miss.colptr(idx1),
      obj.colptr(idx2), miss.colptr(idx2),
      obj.n_rows);
}

template <>
inline double calc_distance<1>(
    const arma::mat &obj,
    const arma::mat &miss,
    arma::uword idx1,
    arma::uword idx2)
{
  return calc_distance_manhattan(
      obj.colptr(idx1), miss.colptr(idx1),
      obj.colptr(idx2), miss.colptr(idx2),
      obj.n_rows);
}

// ============================================================================
// Method dispatch helper
// ============================================================================
template <int Method>
inline double calc_distance_raw(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_miss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_miss,
    const arma::uword n_rows);

template <>
inline double calc_distance_raw<0>(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_miss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_miss,
    const arma::uword n_rows)
{
  return calc_distance_euclidean(target_ptr, target_miss, other_ptr, other_miss, n_rows);
}

template <>
inline double calc_distance_raw<1>(
    const double *__restrict__ target_ptr,
    const double *__restrict__ target_miss,
    const double *__restrict__ other_ptr,
    const double *__restrict__ other_miss,
    const arma::uword n_rows)
{
  return calc_distance_manhattan(target_ptr, target_miss, other_ptr, other_miss, n_rows);
}

// ============================================================================
// Neighbor tracking
// ============================================================================
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

// ============================================================================
// distance_vector — on the fly
// ============================================================================
template <int Method>
std::vector<NeighborInfo> distance_vector_brute(
    const arma::mat &obj,
    const arma::mat &miss,
    const arma::uword target_col,
    const arma::uword n_cols,
    const arma::uword k)
{
  const arma::uword n_rows = obj.n_rows;
  const double *target_ptr = obj.colptr(target_col);
  const double *target_miss = miss.colptr(target_col);

  std::vector<NeighborInfo> top_k;
  top_k.reserve(k);
  arma::uword remaining = k;

  auto process = [&](arma::uword col)
  {
    double dist = calc_distance_raw<Method>(
        target_ptr, target_miss,
        obj.colptr(col), miss.colptr(col),
        n_rows);
    if (remaining > 0)
    {
      insert_before_k(top_k, dist, col);
      --remaining;
    }
    else
    {
      insert_if_better_than_worst(top_k, dist, col);
    }
  };
  for (arma::uword i = 0; i < target_col; ++i)
  {
    process(i);
  }
  for (arma::uword i = target_col + 1; i < n_cols; ++i)
  {
    process(i);
  }
  return top_k;
}

// ============================================================================
// distance_vector — cached path (target pointers hoisted for non-miss lookups)
// ============================================================================
template <int Method>
std::vector<NeighborInfo> distance_vector_cached(
    const arma::mat &obj,
    const arma::mat &miss,
    const arma::uword index,
    const arma::uvec &col_index_miss,
    const arma::uvec &col_index_non_miss,
    const StrictLowerTriangularMatrix &cache,
    const arma::uword k)
{
  const arma::uword n_rows = obj.n_rows;
  const arma::uword target_col = col_index_miss(index);
  const double *target_ptr = obj.colptr(target_col);
  const double *target_miss_ptr = miss.colptr(target_col);

  std::vector<NeighborInfo> top_k;
  top_k.reserve(k);
  arma::uword remaining = k;

  // ---- Fill phase (first k neighbors) ----
  arma::uword i = 0;
  for (; i < index && remaining > 0; ++i)
  {
    double dist = cache(index, i);
    insert_before_k(top_k, dist, col_index_miss(i));
    --remaining;
  }
  arma::uword j = index + 1;
  for (; j < col_index_miss.n_elem && remaining > 0; ++j)
  {
    double dist = cache(j, index);
    insert_before_k(top_k, dist, col_index_miss(j));
    --remaining;
  }
  arma::uword m = 0;
  for (; m < col_index_non_miss.n_elem && remaining > 0; ++m)
  {
    double dist = calc_distance_raw<Method>(
        target_ptr, target_miss_ptr,
        obj.colptr(col_index_non_miss(m)), miss.colptr(col_index_non_miss(m)),
        n_rows);
    insert_before_k(top_k, dist, col_index_non_miss(m));
    --remaining;
  }
  // ---- Replacement phase (remaining neighbors) ----
  for (; i < index; ++i)
  {
    double dist = cache(index, i);
    insert_if_better_than_worst(top_k, dist, col_index_miss(i));
  }
  for (; j < col_index_miss.n_elem; ++j)
  {
    double dist = cache(j, index);
    insert_if_better_than_worst(top_k, dist, col_index_miss(j));
  }
  for (; m < col_index_non_miss.n_elem; ++m)
  {
    double dist = calc_distance_raw<Method>(
        target_ptr, target_miss_ptr,
        obj.colptr(col_index_non_miss(m)), miss.colptr(col_index_non_miss(m)),
        n_rows);
    insert_if_better_than_worst(top_k, dist, col_index_non_miss(m));
  }
  return top_k;
}

// ============================================================================
// Main imputation function
// ============================================================================
template <int Method, bool cache>
void impute_knn_brute_impl(
    arma::mat &result,
    const arma::mat &obj,
    const arma::mat &miss,
    const arma::uword k,
    const arma::uvec &col_index_miss,
    const arma::uvec &col_offsets,
    const double dist_pow,
    int cores,
    const arma::uvec *col_index_non_miss_ptr,
    const StrictLowerTriangularMatrix *cache_ptr)
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    std::vector<NeighborInfo> top_k;
    if constexpr (cache)
    {
      top_k = distance_vector_cached<Method>(
          obj, miss, i, col_index_miss,
          *col_index_non_miss_ptr, *cache_ptr, k);
    }
    else
    {
      top_k = distance_vector_brute<Method>(
          obj, miss, col_index_miss(i), obj.n_cols, k);
    }

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
    impute_column_values(result, obj, miss,
                         col_offsets(i), col_index_miss(i),
                         nn_columns_mat, weights);
  }
}

// ============================================================================
// Cache setup + dispatch
// ============================================================================
template <int Method>
void dispatch_cache(
    arma::mat &result,
    const arma::mat &obj,
    const arma::mat &miss,
    const arma::uword k,
    const arma::uvec &col_index_miss,
    const arma::uvec &col_offsets,
    const arma::uvec &n_col_miss,
    const double dist_pow,
    const bool use_cache,
    int cores)
{
  if (use_cache)
  {
    arma::uvec col_index_non_miss = arma::find(n_col_miss == 0);
    StrictLowerTriangularMatrix cache(col_index_miss.n_elem);

#ifdef _OPENMP
    omp_set_num_threads(cores);
#pragma omp parallel for schedule(dynamic)
#endif
    for (arma::uword row = 1; row < col_index_miss.n_elem; ++row)
    {
      for (arma::uword col = 0; col < row; ++col)
      {
        cache(row, col) = calc_distance<Method>(
            obj, miss, col_index_miss(row), col_index_miss(col));
      }
    }

    impute_knn_brute_impl<Method, true>(
        result, obj, miss, k, col_index_miss, col_offsets, dist_pow, cores,
        &col_index_non_miss, &cache);
  }
  else
  {
    impute_knn_brute_impl<Method, false>(
        result, obj, miss, k, col_index_miss, col_offsets, dist_pow, cores,
        nullptr, nullptr);
  }
}

// [[Rcpp::export]]
arma::mat impute_knn_brute(
    const arma::mat &obj,
    const arma::mat &miss,
    const arma::uword k,
    const arma::uvec &n_col_miss,
    const int method,
    const double dist_pow,
    const bool cache = true,
    int cores = 1)
{
  arma::uvec col_index_miss = arma::find(n_col_miss > 0);
  arma::uvec col_offsets;
  arma::mat result = initialize_result_matrix(miss, col_index_miss, n_col_miss, col_offsets);
  if (result.n_rows == 0)
  {
    return result;
  }
  switch (method)
  {
  case 0:
    dispatch_cache<0>(result, obj, miss, k, col_index_miss, col_offsets,
                      n_col_miss, dist_pow, cache, cores);
    break;
  case 1:
    dispatch_cache<1>(result, obj, miss, k, col_index_miss, col_offsets,
                      n_col_miss, dist_pow, cache, cores);
    break;
  default:
    throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");
  }

  return result;
}
