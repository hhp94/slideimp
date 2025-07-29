// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace Rcpp;

struct StrictLowerTriangularMatrix
{
  size_t n;                       // size of triangular matrix (4 for a 4x4 matrix)
  std::vector<double> data;       // is a dynamic vector double underneath
  std::vector<size_t> row_starts; // cache row offsets positions

  explicit StrictLowerTriangularMatrix(size_t size) : n(size), data(), row_starts()
  {
    if (size == 0)
    {
      throw std::invalid_argument("Size must be at least 1");
    }
    // Overflow check for n*(n-1)/2
    size_t max_val = std::numeric_limits<size_t>::max();
    if (size > 1 && (size - 1) > max_val / size)
    {
      throw std::overflow_error("Matrix size too large. Reduce n");
    }

    size_t num_elements = (size * (size - 1)) / 2;
    data.assign(num_elements, 0.0);

    row_starts.resize(size);
    size_t offset = 0;
    // pre-calculate row offsets at initialization once to make operation() faster
    for (size_t i = 0; i < size; ++i)
    {
      row_starts[i] = offset;
      offset += i;
    }
  }

  double &operator()(size_t i, size_t j) noexcept
  {
    return data[row_starts[i] + j];
  }

  const double &operator()(size_t i, size_t j) const noexcept
  {
    return data[row_starts[i] + j];
  }

  size_t size() const noexcept { return n; }

  void print() const
  {
    for (size_t i = 0; i < n; ++i)
    {
      for (size_t j = 0; j < n; ++j)
      {
        if (j < i)
        {
          Rcout << (*this)(i, j) << " ";
        }
        else
        {
          Rcout << "0 ";
        }
      }
      Rcout << std::endl;
    }
  }
};

// Helper function for squared Euclidean (method 0)
double calc_distance_sqeuclid(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2,
    const double total_rows)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);

  double dist = 0.0;
  arma::uword n_valid = 0;
#if defined(_OPENMP)
#pragma omp simd reduction(+ : dist) reduction(+ : n_valid)
#endif
  for (arma::uword r = 0; r < obj.n_rows; ++r)
  {
    if (miss1_ptr[r] == 0 && miss2_ptr[r] == 0)
    {
      double diff = col1_ptr[r] - col2_ptr[r];
      dist += diff * diff;
      ++n_valid;
    }
  }

  if (n_valid == 0)
  {
    return arma::datum::inf;
  }

  dist *= (total_rows / n_valid);

  return dist;
}

// Helper function for Manhattan (method 1)
double calc_distance_manhattan(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2,
    const double total_rows)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);

  double dist = 0.0;
  arma::uword n_valid = 0;
#if defined(_OPENMP)
#pragma omp simd reduction(+ : dist) reduction(+ : n_valid)
#endif
  for (arma::uword r = 0; r < obj.n_rows; ++r)
  {
    if (miss1_ptr[r] == 0 && miss2_ptr[r] == 0)
    {
      double diff = col1_ptr[r] - col2_ptr[r];
      dist += std::abs(diff);
      ++n_valid;
    }
  }

  if (n_valid == 0)
  {
    return arma::datum::inf;
  }

  dist *= (total_rows / n_valid);

  return dist;
}

// Helper function for impute.knn (method 2)
double calc_distance_knn(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2,
    const double total_rows)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);

  double dist = 0.0;
  arma::uword n_valid = 0;
#if defined(_OPENMP)
#pragma omp simd reduction(+ : dist) reduction(+ : n_valid)
#endif
  for (arma::uword r = 0; r < obj.n_rows; ++r)
  {
    if (miss1_ptr[r] == 0 && miss2_ptr[r] == 0)
    {
      double diff = col1_ptr[r] - col2_ptr[r];
      dist += diff * diff;
      ++n_valid;
    }
  }
  if (n_valid == 0)
  {
    return arma::datum::inf;
  }

  dist /= n_valid;
  return dist;
}

using dist_func_t = double (*)(
    const arma::mat &,
    const arma::umat &,
    const arma::uword,
    const arma::uword,
    const double);

arma::vec distance_vector(
    const arma::mat &obj,                     // data with missing
    const arma::umat &miss,                   // miss matrix correspond to full data
    const arma::uword &index,                 // index of `index_miss`, not the obj
    const arma::uvec &index_miss,             // actual index of all missing columns, refers to the obj
    const arma::uvec &index_not_miss,         // actual index of all non missing columns, refers to the obj
    const StrictLowerTriangularMatrix &cache, // cache of value using lower triangular matrix
    const dist_func_t calc_dist)
{
  // Each missing column needs a vector of distance to all other columns. We iterate
  // column wise instead of calculating all the distances upfront because the distance
  // from missing to non-missing columns can be discarded after each iteration to
  // lower memory usage. This keeps memory usage to minimum.
  arma::vec dist_vec(obj.n_cols, arma::fill::none);
  const double total_rows = static_cast<double>(obj.n_rows);

  // dist_vec is arranged in two parts:: all the missing columns, then the non missing columns

  // First iterate through all the missing column. Missing columns has 3 parts: 1) i < index
  // if the current column is < index, then the value has already been calculated, fetch from cache.
  // This will run from 0 to right before index
  for (arma::uword i = 0; i < index; ++i)
  {
    dist_vec(i) = cache(index, i);
  }
  // 2) then the next position is the index, which is the distance between the column and itself,
  // which is inf for the sake of sorting
  dist_vec(index) = arma::datum::inf;
  // 3) then from index to index_miss.n_elem, we calculate and cache the values. Automatically
  // skiped for last missing column where all values are fetched from cache
  for (arma::uword i = index + 1; i < index_miss.n_elem; ++i)
  {
    dist_vec(i) = cache(i, index);
  }
  // Compute the remaining part of the vector (distance from missing to non missing columns).
  // If all columns have missing values (i.e., the matrix is square with
  // obj.n_cols == index_miss.n_elem), this loop is automatically skipped.
  for (arma::uword i = index_miss.n_elem; i < obj.n_cols; ++i)
  {
    // Offset i to index into index_not_miss from 0
    dist_vec(i) = calc_dist(obj, miss, index_miss(index), index_not_miss(i - index_miss.n_elem), total_rows);
  }

  return dist_vec;
}

arma::uvec find_knn_indices_arma(const arma::vec &distances, arma::uword k)
{
  arma::uvec indices = arma::regspace<arma::uvec>(0, distances.n_elem - 1);
  // Partition the indices vector such that the k-th smallest element is in its
  // head k position. Elements before it are smaller or equal.
  std::nth_element(
      indices.begin(),
      // nth_element is 0-indexed, so we use k-1 to get the k-th element
      indices.begin() + (k - 1),
      indices.end(),
      // Custom comparator lambda function that compares the distances at the given indices.
      [&](arma::uword a, arma::uword b)
      {
        return distances(a) < distances(b);
      });

  // Return the first k indices, which now correspond to the k smallest distances.
  return indices.head(k);
}

//' Impute missing values in a matrix using k-nearest neighbors (k-NN)
//'
//' This function imputes missing values in a matrix using a k-nearest neighbors
//' approach based on the specified distance metric. It processes only columns
//' with missing values, calculating distances to other columns to find the
//' k-nearest neighbors, and imputes missing values by averaging non-missing
//' values from these neighbors.
//'
//' @param obj Numeric matrix with missing values represented as NA (NaN).
//' @param miss Logical matrix (0/1) indicating missing values (1 = missing).
//' @param k Number of nearest neighbors to use for imputation.
//' @param n_col_miss Integer vector specifying the count of missing values per column.
//' @param method Integer specifying the distance metric: 0 = Euclidean, 1 = Manhattan, 2 = impute.knn.
//' @param cores Number of CPU cores to use for parallel processing (default = 1).
//' @return A matrix with imputed values where missing values were present.
// [[Rcpp::export]]
arma::mat impute_knn_naive(
    const arma::mat &obj,        // data with NA (Nan) not prefilled
    const arma::umat &miss,      // missing data matrix
    const arma::uword k,         // n neighbors
    const arma::uvec n_col_miss, // vector of missing per column
    const int method,
    int cores = 1)
{
  // pick the distance method.
  dist_func_t calc_dist = nullptr;

  switch (method)
  {
  case 0:
    calc_dist = calc_distance_sqeuclid;
    break;
  case 1:
    calc_dist = calc_distance_manhattan;
    break;
  case 2:
    calc_dist = calc_distance_knn;
    break;
  default:
    throw std::invalid_argument("0 = Euclid; 1 = Manhattan; 2 = impute.knn.");
  }

  // Only compute distance of cols with any missing and do it once
  arma::uvec col_index_miss = arma::find(n_col_miss);
  if (col_index_miss.n_elem == 0)
  {
    return obj;
  }
  arma::uvec col_index_non_miss = arma::find(n_col_miss == 0);
  arma::mat imputed = obj; // Copy full object.
  // To use the `distance_vector` function, we have to rearrange the neighbor into
  // missing then non-missing columns. This index is the source of truth needed
  // to get back the actual column after nearest neighbor calculation
  arma::uvec neighbor_index = arma::join_vert(col_index_miss, col_index_non_miss);
  // Initialize the cache with size of col_index_miss
  StrictLowerTriangularMatrix cache(col_index_miss.n_elem);
  double total_rows = static_cast<double>(obj.n_rows);
#if defined(_OPENMP)
  omp_set_num_threads(cores);
#else
  cores = 1;
#endif
  // Fill the cache upfront with pairwise distances between missing columns
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for (arma::uword row = 1; row < col_index_miss.n_elem; ++row)
  {
    for (arma::uword col = 0; col < row; ++col)
    {
      cache(row, col) = calc_dist(obj, miss, col_index_miss(row), col_index_miss(col), total_rows);
    }
  }
  // Loop through only missing columns
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    arma::vec dist_vec = distance_vector(
        obj, miss, i, col_index_miss, col_index_non_miss, cache, calc_dist);
    arma::uvec candidate_indices = arma::find_finite(dist_vec);
    if (candidate_indices.n_elem == 0)
    {
      continue;
    }
    // Caps `k` to the number of available neighbors.
    arma::uword n_neighbors = std::min(k, candidate_indices.n_elem);
    // index of nn relative to dist_vec
    arma::uvec nn_dist_indices = find_knn_indices_arma(dist_vec, n_neighbors);
    // actual index of nn columns relative to obj columns
    arma::uvec nn_columns = neighbor_index.elem(nn_dist_indices);

    // now we move on to nn imputation. Find the rows to impute for this column
    arma::uword target_col_idx = col_index_miss(i);
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));

    // man, range-for-loop is trippy. kind of not used to this but it works
    for (arma::uword row_idx : rows_to_impute)
    {
      double sum = 0.0;
      arma::uword count = 0;
#if defined(_OPENMP)
#pragma omp simd reduction(+ : sum) reduction(+ : count)
#endif
      for (arma::uword neighbor_col_idx : nn_columns)
      {
        if (miss(row_idx, neighbor_col_idx) == 0)
        {
          sum += obj(row_idx, neighbor_col_idx);
          ++count;
        }
      }
      if (count > 0)
      {
        imputed(row_idx, target_col_idx) = sum / count;
      }
    }
  }

  return imputed;
}
