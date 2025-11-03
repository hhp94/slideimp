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
          Rcpp::Rcout << (*this)(i, j) << " ";
        }
        else
        {
          Rcpp::Rcout << "0 ";
        }
      }
      Rcpp::Rcout << std::endl;
    }
  }
};

// Helper function for Euclidean (method 0)
double calc_distance_euclid(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);
  double dist = 0.0;
  arma::uword n_valid = 0;
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

  dist /= static_cast<double>(n_valid);

  return dist;
}

// Helper function for Manhattan (method 1)
double calc_distance_manhattan(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword idx1,
    const arma::uword idx2)
{
  const double *col1_ptr = obj.colptr(idx1);
  const double *col2_ptr = obj.colptr(idx2);
  const arma::uword *miss1_ptr = miss.colptr(idx1);
  const arma::uword *miss2_ptr = miss.colptr(idx2);
  double dist = 0.0;
  arma::uword n_valid = 0;
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

  dist /= static_cast<double>(n_valid);

  return dist;
}

// struct hold the distance and index relative to distance_vector for each neighbor
struct NeighborInfo
{
  double distance;   // calculated distance to this neighbor
  arma::uword index; // index of this neighbor relative to distance_vector
  NeighborInfo(double d, arma::uword i) : distance(d), index(i) {}
};

// always insert when we have fewer than k neighbors
inline void insert_before_k(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
  // inf will always sort to the end. So 1 branch only
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

// only insert if better than worst when we have exactly k neighbors
inline void insert_if_better_than_worst(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx)
{
  // inf automatically return. Only 1 branch left
  if (dist >= top_k.back().distance)
    return;
  // replace the worst neighbor
  top_k.back() = NeighborInfo(dist, idx);
  // Bubble up to maintain sorted order
  for (size_t i = top_k.size() - 1; i > 0 && top_k[i].distance < top_k[i - 1].distance; --i)
  {
    std::swap(top_k[i], top_k[i - 1]);
  }
}

using dist_func_t = double (*)(
    const arma::mat &,
    const arma::umat &,
    const arma::uword,
    const arma::uword);

std::vector<NeighborInfo> distance_vector(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword &index,
    const arma::uvec &index_miss,
    const arma::uvec &index_not_miss,
    const StrictLowerTriangularMatrix &cache,
    const dist_func_t calc_dist,
    const arma::uword k)
{
  std::vector<NeighborInfo> top_k_neighbors;
  top_k_neighbors.reserve(k);
  // how many neighbors we need to process before changing strategy
  arma::uword remaining_to_fill = k;

  // From 1 to k: fill up to k neighbors
  // Same strategy as before of splitting from 1 -> index -> beyond
  arma::uword i = 0;
  for (; i < index && remaining_to_fill > 0; ++i)
  {
    double dist = cache(index, i);
    insert_before_k(top_k_neighbors, dist, i);
    --remaining_to_fill;
  }
  // columns after current index
  arma::uword j = index + 1;
  for (; j < index_miss.n_elem && remaining_to_fill > 0; ++j)
  {
    double dist = cache(j, index);
    insert_before_k(top_k_neighbors, dist, j);
    --remaining_to_fill;
  }
  // non-missing columns
  arma::uword m = 0;
  for (; m < index_not_miss.n_elem && remaining_to_fill > 0; ++m)
  {
    double dist = calc_dist(obj, miss, index_miss(index), index_not_miss(m));
    arma::uword global_idx = index_miss.n_elem + m;
    insert_before_k(top_k_neighbors, dist, global_idx);
    --remaining_to_fill;
  }
  // From k to n: replacement phase - only insert if better than worst and no k left
  for (; i < index; ++i)
  {
    double dist = cache(index, i);
    insert_if_better_than_worst(top_k_neighbors, dist, i);
  }
  // columns after current index
  for (; j < index_miss.n_elem; ++j)
  {
    double dist = cache(j, index);
    insert_if_better_than_worst(top_k_neighbors, dist, j);
  }
  // non-missing columns (if still building)
  for (; m < index_not_miss.n_elem; ++m)
  {
    double dist = calc_dist(obj, miss, index_miss(index), index_not_miss(m));
    arma::uword global_idx = index_miss.n_elem + m;
    insert_if_better_than_worst(top_k_neighbors, dist, global_idx);
  }
  return top_k_neighbors;
}

//' Impute missing values in a matrix using k-nearest neighbors (k-NN) with brute-force
//'
//' This function imputes missing values in a matrix using a k-nearest neighbors
//' approach based on the specified distance metric. It processes only columns
//' with missing values, calculating distances to other columns to find the
//' k-nearest neighbors, and imputes missing values using either a simple average
//' or a weighted average (inverse distance weighting) of non-missing values from
//' these neighbors, depending on the `weighted` parameter.
//'
//'
//' @param obj Numeric matrix with missing values represented as NA (NaN).
//' @param miss Logical matrix (0/1) indicating missing values (1 = missing).
//' @param k Number of nearest neighbors to use for imputation.
//' @param n_col_miss Integer vector specifying the count of missing values per column.
//' @param method Integer specifying the distance metric: 0 = Euclidean, 1 = Manhattan.
//' @param dist_pow A positive double that controls the penalty for larger distances in
//' the weighted mean imputation. Must be greater than zero: values between 0 and 1 apply a softer penalty,
//' 1 is linear (default), and values greater than 1 apply a harsher penalty.
//' @param cores Number of CPU cores to use for parallel processing (default = 1).
//' @return A matrix where the first column is the 1-based row index, the second column is the 1-based column index,
//' and the subsequent `n_imp` columns contain the imputed values (one column per bootstrap replicate if `n_imp > 1`).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat impute_knn_brute(
    const arma::mat &obj,         // Input data matrix with missing values as NaN
    const arma::umat &miss,       // Matrix of same size as obj, with 1 for missing (NA) and 0 for present
    const arma::uword k,          // Number of nearest neighbors to use
    const arma::uvec &n_col_miss, // Vector containing the count of missing values per column
    const int method,             // Distance metric: 0=Euclidean, 1=Manhattan method.
    const double dist_pow,        // Power for distance penalty in weighted average
    int cores = 1)                // Number of cores for parallel processing
{
  // Select the distance calculation function based on the method
  dist_func_t calc_dist = nullptr;
  switch (method)
  {
  case 0:
    calc_dist = calc_distance_euclid;
    break;
  case 1:
    calc_dist = calc_distance_manhattan;
    break;
  default:
    throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");
  }
  // Find columns that contain missing values
  arma::uvec col_index_miss = arma::find(n_col_miss > 0);

  // Initialize result matrix and get column offsets using helper function
  arma::uvec col_offsets;
  arma::mat result = initialize_result_matrix(miss, col_index_miss, n_col_miss, col_offsets);
  if (result.n_rows == 0)
  {
    return result; // No missing values
  }
  // For finding neighbors, we need all columns. Separate them into missing and non-missing.
  arma::uvec col_index_non_miss = arma::find(n_col_miss == 0);
  // Create a combined index vector to look up original column indices after finding neighbors.
  arma::uvec neighbor_index = arma::join_vert(col_index_miss, col_index_non_miss);
  // Initialize a cache to store pairwise distances between columns that have missing data.
  StrictLowerTriangularMatrix cache(col_index_miss.n_elem);
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  // Pre-fill the cache with distances in parallel. A bit unoptimal because it can be
  // that k neighbors
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (arma::uword row = 1; row < col_index_miss.n_elem; ++row)
  {
    for (arma::uword col = 0; col < row; ++col)
    {
      cache(row, col) = calc_dist(obj, miss, col_index_miss(row), col_index_miss(col));
    }
  }
  // Main imputation loop
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    // Get top-k neighbors directly
    std::vector<NeighborInfo> top_k = distance_vector(
        obj, miss, i, col_index_miss, col_index_non_miss, cache, calc_dist, k);
    arma::uword n_neighbors = top_k.size();
    if (n_neighbors == 0)
    {
      continue; // Skip if no valid neighbors are found
    }
    // Extract neighbor indices (already sorted by distance)
    arma::uvec nn_columns(n_neighbors);
    arma::vec nn_dists(n_neighbors);
    for (arma::uword j = 0; j < n_neighbors; ++j)
    {
      nn_columns(j) = neighbor_index(top_k[j].index);
      nn_dists(j) = top_k[j].distance;
    }
    // Calculate weights once (same for all methods)
    arma::vec weights(n_neighbors);
    weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
    arma::uword target_col_idx = col_index_miss(i);
    // Choose imputation method based on n_pmm

    // Single deterministic imputation
    arma::umat nn_columns_mat(n_neighbors, 1);
    nn_columns_mat.col(0) = nn_columns;

    impute_column_values(
        result, obj, miss,
        col_offsets(i), target_col_idx,
        nn_columns_mat, weights);
  }
  return result;
}

//' @title Find K-Nearest Neighbors for Columns with Missing Values
//'
//' @description
//' This function returns the information of nearest neighbors found.
//' Used for benchmarking and unit testing logic of \code{impute_knn_brute} only.
//'
//' @param obj R matrix
//' @param miss is.na(obj)
//' @param k n neighbor
//' @param n_col_miss Integer vector specifying the count of missing values per column.
//' @param n_col_name Character vector of same length as n_col_miss containing column names.
//' @param method distance metric
//' @param cores n cores
//'
//' @return
//' A named list where each element corresponds to one column with missing values
//' (i.e., where n_col_miss > 0). Each sub-list contains:
//'
//' \itemize{
//'   \item \code{indices}: An integer vector of 1-based indices of the nearest neighbor columns.
//'   \item \code{distances}: A numeric vector of distances to the nearest neighbors.
//'   \item \code{n_neighbors}: The number of valid neighbors found.
//' }
//'
//' If no columns have missing values, an empty list is returned.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List find_knn_brute(
    const arma::mat &obj,                    // Input data matrix with missing values as NaN
    const arma::umat &miss,                  // Matrix of same size as obj, with 1 for missing (NA) and 0 for present
    const arma::uword k,                     // Number of nearest neighbors to use
    const arma::uvec &n_col_miss,            // Vector containing the count of missing values per column
    const Rcpp::CharacterVector &n_col_name, // Vector of column names (same length as n_col_miss)
    const int method,                        // Distance metric: 0=Euclidean, 1=Manhattan method
    int cores = 1)                           // Number of cores for parallel processing
{
  // Select the distance calculation function based on the method
  dist_func_t calc_dist = nullptr;
  switch (method)
  {
  case 0:
    calc_dist = calc_distance_euclid;
    break;
  case 1:
    calc_dist = calc_distance_manhattan;
    break;
  default:
    throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan, 2=impute.knn");
  }

  arma::uvec col_index_miss = arma::find(n_col_miss > 0);
  if (col_index_miss.n_elem == 0)
  {
    return Rcpp::List();
  }
  arma::uvec col_index_non_miss = arma::find(n_col_miss == 0);
  arma::uvec neighbor_index = arma::join_vert(col_index_miss, col_index_non_miss);
  StrictLowerTriangularMatrix cache(col_index_miss.n_elem);
#if defined(_OPENMP)
  omp_set_num_threads(cores);
#endif
  // Pre-fill the cache with distances in parallel
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
  for (arma::uword row = 1; row < col_index_miss.n_elem; ++row)
  {
    for (arma::uword col = 0; col < row; ++col)
    {
      cache(row, col) = calc_dist(obj, miss, col_index_miss(row), col_index_miss(col));
    }
  }
  // Initialize result list for all columns with missing values
  Rcpp::List result(col_index_miss.n_elem);

  // Main loop: iterate through all columns with missing values
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    // Get top-k neighbors directly
    std::vector<NeighborInfo> top_k = distance_vector(
        obj, miss, i, col_index_miss, col_index_non_miss, cache, calc_dist, k);
    arma::uword n_neighbors = top_k.size();
    Rcpp::List neighbor_info;
    if (n_neighbors == 0)
    {
      // No valid neighbors found
      neighbor_info = Rcpp::List::create(
          Rcpp::Named("indices") = Rcpp::IntegerVector(),
          Rcpp::Named("distances") = Rcpp::NumericVector(),
          Rcpp::Named("n_neighbors") = 0);
    }
    else
    {
      // Extract neighbor indices and distances
      Rcpp::IntegerVector nn_indices(n_neighbors);
      Rcpp::NumericVector nn_distances(n_neighbors);

      for (arma::uword j = 0; j < n_neighbors; ++j)
      {
        // Convert back to R 1-based indexing
        nn_indices[j] = neighbor_index(top_k[j].index) + 1;
        nn_distances[j] = top_k[j].distance;
      }

      neighbor_info = Rcpp::List::create(
          Rcpp::Named("indices") = nn_indices,
          Rcpp::Named("distances") = nn_distances,
          Rcpp::Named("n_neighbors") = n_neighbors);
    }

    // Store in result list
    result[i] = neighbor_info;
  }

  // Set names for the list elements from n_col_name
  Rcpp::CharacterVector result_names(col_index_miss.n_elem);
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    result_names[i] = n_col_name[col_index_miss(i)];
  }
  result.names() = result_names;

  return result;
}
