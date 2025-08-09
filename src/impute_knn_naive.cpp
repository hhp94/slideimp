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
  double weight = static_cast<double>(obj.n_rows) / n_valid;
  return std::sqrt(weight * dist);
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
  double weight = static_cast<double>(obj.n_rows) / n_valid;
  return weight * dist;
}

// Helper function for impute.knn (method 2)
double calc_distance_knn(
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

  dist /= n_valid;

  return dist;
}

using dist_func_t = double (*)(
    const arma::mat &,
    const arma::umat &,
    const arma::uword,
    const arma::uword);

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
  // 3) then from index to index_miss.n_elem, the cache is in the symmetric position
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
    dist_vec(i) = calc_dist(obj, miss, index_miss(index), index_not_miss(i - index_miss.n_elem));
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

constexpr double epsilon = 1e-10;

//' Impute missing values in a matrix using k-nearest neighbors (k-NN)
//'
//' This function imputes missing values in a matrix using a k-nearest neighbors
//' approach based on the specified distance metric. It processes only columns
//' with missing values, calculating distances to other columns to find the
//' k-nearest neighbors, and imputes missing values using either a simple average
//' or a weighted average (inverse distance weighting) of non-missing values from
//' these neighbors, depending on the `weighted` parameter.
//'
//' @param obj Numeric matrix with missing values represented as NA (NaN).
//' @param miss Logical matrix (0/1) indicating missing values (1 = missing).
//' @param k Number of nearest neighbors to use for imputation.
//' @param n_col_miss Integer vector specifying the count of missing values per column.
//' @param method Integer specifying the distance metric: 0 = Euclidean, 1 = Manhattan.
//' @param weighted Boolean controls for the imputed value to be a simple mean or weighted mean by inverse distance.
//' @param dist_pow A positive double that controls the penalty for larger distances in
//' the weighted mean imputation. Must be greater than zero: values between 0 and 1 apply a softer penalty,
//' 1 is linear (default), and values greater than 1 apply a harsher penalty.
//' @param cores Number of CPU cores to use for parallel processing (default = 1).
//' @return A matrix where the first column is the index of missing value as if calculated in R and column > 2 are imputed values.
//'
//' @export
// [[Rcpp::export]]
arma::mat impute_knn_naive(
    const arma::mat &obj,         // Input data matrix with missing values as NaN
    const arma::umat &miss,       // Matrix of same size as obj, with 1 for missing (NA) and 0 for present
    const arma::uword k,          // Number of nearest neighbors to use
    const arma::uvec &n_col_miss, // Vector containing the count of missing values per column
    const int method,             // Distance metric: 0=Euclidean, 1=Manhattan, 2=impute.knn's method
    const bool weighted,          // TRUE for a weighted average, FALSE for a simple average
    const double dist_pow,        // Power for distance penalty in weighted average
    int cores)                    // Number of cores for parallel processing
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
  case 2:
    calc_dist = calc_distance_knn;
    break;
  default:
    throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan, 2=impute.knn");
  }

  // Find columns that contain missing values
  arma::uvec col_index_miss = arma::find(n_col_miss > 0);
  if (col_index_miss.n_elem == 0)
  {
    // Return empty matrix if no missing values exist
    return arma::mat(0, 2);
  }

  // For finding neighbors, we need all columns. Separate them into missing/non-missing.
  arma::uvec col_index_non_miss = arma::find(n_col_miss == 0);

  // Create a combined index vector to look up original column indices after finding neighbors.
  // This groups columns with missing data first, which simplifies distance calculations.
  arma::uvec neighbor_index = arma::join_vert(col_index_miss, col_index_non_miss);

  // Initialize a cache to store pairwise distances between columns that have missing data.
  // This avoids redundant distance calculations.
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

  // Initialize the result matrix: [1-based linear index, imputed value]
  // This allows for easy merging back into a matrix in R.
  arma::uword sum_missing = arma::accu(n_col_miss);
  arma::mat result(sum_missing, 2);
  result.fill(arma::datum::nan);

  // To efficiently fill the 'result' matrix, we calculate offsets.
  // 'col_offsets(i)' will be the starting row in 'result' for the i-th column with missing data.
  arma::uvec miss_counts = n_col_miss.elem(col_index_miss);
  arma::uvec col_offsets(miss_counts.n_elem + 1);
  col_offsets.fill(arma::fill::zeros);
  col_offsets.subvec(1, miss_counts.n_elem) = arma::cumsum(miss_counts);

  // Pre-fill result with linear indices to avoid "continue" skipping the indices
  arma::uword offset = 0;
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    const arma::uword target_col_idx = col_index_miss(i);
    const arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
      const arma::uword row_idx = rows_to_impute(r);
      const arma::uword res_row = offset + r;
      result(res_row, 0) = target_col_idx * obj.n_rows + row_idx + 1;
    }
    offset += rows_to_impute.n_elem;
  }

  // Main imputation loop: iterate through only the columns that need imputation
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    // Calculate distances from the current target column to all other potential neighbor columns
    arma::vec dist_vec = distance_vector(
        obj, miss, i, col_index_miss, col_index_non_miss, cache, calc_dist);

    // Find available neighbors (those with a finite distance)
    arma::uvec candidate_indices = arma::find_finite(dist_vec);
    arma::uword n_neighbors = std::min(k, candidate_indices.n_elem);

    if (n_neighbors == 0)
    {
      continue; // Skip if no valid neighbors are found
    }

    // Find the indices of the 'k' nearest neighbors relative to 'dist_vec'
    arma::uvec nn_dist_indices = find_knn_indices_arma(dist_vec, n_neighbors);

    // Get the original column indices (in 'obj') of these neighbors
    arma::uvec nn_columns = neighbor_index.elem(nn_dist_indices);
    // Get the distances to these neighbors
    arma::vec nn_dists = dist_vec.elem(nn_dist_indices);

    // Precompute weights for each neighbor (these are the same for all missing rows in this column)
    arma::vec nn_weights(nn_columns.n_elem);
    if (weighted)
    {
      nn_weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
    }
    else
    {
      nn_weights.fill(1.0); // Simple average means all weights are 1
    }

    // Get the original index of the column we are currently imputing
    arma::uword target_col_idx = col_index_miss(i);
    // Find which rows are missing in this specific column
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));

    // For each missing cell in this column, calculate the imputed value
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
      arma::uword row_idx = rows_to_impute(r); // The actual row index of the missing value
      double weighted_sum = 0.0;
      double weight_total = 0.0;

      // Aggregate values from neighbors
      for (arma::uword j = 0; j < nn_columns.n_elem; ++j)
      {
        arma::uword neighbor_col_idx = nn_columns(j);
        // A neighbor can only contribute if its value in the same row is NOT missing
        if (miss(row_idx, neighbor_col_idx) == 0)
        {
          double weight = nn_weights(j);
          weighted_sum += weight * obj(row_idx, neighbor_col_idx);
          weight_total += weight;
        }
      }

      double imputed_value = (weight_total > 0.0) ? (weighted_sum / weight_total) : arma::datum::nan;

      // Find the correct row in the result matrix to write to
      const arma::uword res_row = col_offsets(i) + r;
      result(res_row, 1) = imputed_value;
    }
  }
  return result;
}
