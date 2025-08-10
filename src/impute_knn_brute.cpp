// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <random>
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

// struct hold the distance and index relative to distance_vector for each neighbor
struct NeighborInfo
{
  double distance;   // The calculated distance to this neighbor
  arma::uword index; // The index of this neighbor relative to distance_vector
  NeighborInfo(double d, arma::uword i) : distance(d), index(i) {}
};

void insert_if_better(std::vector<NeighborInfo> &top_k, double dist, arma::uword idx, arma::uword k)
{
  // Step 1: Reject infinite distances immediately
  if (!std::isfinite(dist))
    return;

  // Step 2: If we haven't filled k neighbors yet, always insert
  if (top_k.size() < k)
  {
    // Find the correct position to maintain sorted order (smallest distances first)
    auto it = std::upper_bound(
        top_k.begin(),
        top_k.end(),
        dist,
        [](double d, const NeighborInfo &ni)
        {
          return d < ni.distance;
        });
    // Insert the new neighbor at the correct position
    top_k.insert(it, NeighborInfo(dist, idx));
  }
  // Step 3: If we already have k neighbors, only insert if this one is better
  else if (dist < top_k.back().distance)
  {
    // Replace the worst neighbor (last in sorted list)
    top_k.back() = NeighborInfo(dist, idx);

    // Bubble the new element to its correct position
    // Since we only changed the last element, we just need to bubble it forward
    for (size_t i = top_k.size() - 1; i > 0 && top_k[i].distance < top_k[i - 1].distance; --i)
    {
      std::swap(top_k[i], top_k[i - 1]);
    }
  }
  // Step 4: If dist >= top_k.back().distance, do nothing (this neighbor is worse)
}

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

  // Process missing columns before current index (0 to index-1)
  for (arma::uword i = 0; i < index; ++i)
  {
    double dist = cache(index, i);
    insert_if_better(top_k_neighbors, dist, i, k);
  }

  // Process missing columns after current index (index+1 to n_elem-1)
  for (arma::uword i = index + 1; i < index_miss.n_elem; ++i)
  {
    double dist = cache(i, index);
    insert_if_better(top_k_neighbors, dist, i, k);
  }

  // Process non-missing columns (unchanged)
  for (arma::uword i = 0; i < index_not_miss.n_elem; ++i)
  {
    double dist = calc_dist(obj, miss, index_miss(index), index_not_miss(i));
    arma::uword global_idx = index_miss.n_elem + i; // Offset to match original indexing
    insert_if_better(top_k_neighbors, dist, global_idx, k);
  }

  return top_k_neighbors;
}

constexpr double epsilon = 1e-10;

//' Impute missing values in a matrix using k-nearest neighbors (k-NN) with brute-force
//'
//' This function imputes missing values in a matrix using a k-nearest neighbors
//' approach based on the specified distance metric. It processes only columns
//' with missing values, calculating distances to other columns to find the
//' k-nearest neighbors, and imputes missing values using either a simple average
//' or a weighted average (inverse distance weighting) of non-missing values from
//' these neighbors, depending on the `weighted` parameter.
//'
//' When `nboot > 1`, bootstrapping is enabled: for each missing value, `nboot` imputed values are generated
//' by resampling the k nearest neighbors with replacement and using simple averages (weighted is forced to FALSE).
//' This provides variability estimates for imputation uncertainty. Random number generation is seeded for reproducibility,
//' with per-column offsets to ensure independence.
//'
//' @param obj Numeric matrix with missing values represented as NA (NaN).
//' @param miss Logical matrix (0/1) indicating missing values (1 = missing).
//' @param k Number of nearest neighbors to use for imputation.
//' @param n_col_miss Integer vector specifying the count of missing values per column.
//' @param method Integer specifying the distance metric: 0 = Euclidean, 1 = Manhattan, 2 = impute.knn method.
//' @param weighted Boolean controls for the imputed value to be a simple mean or weighted mean by inverse distance.
//'   Note: Forced to FALSE when `nboot > 1`.
//' @param dist_pow A positive double that controls the penalty for larger distances in
//' the weighted mean imputation. Must be greater than zero: values between 0 and 1 apply a softer penalty,
//' 1 is linear (default), and values greater than 1 apply a harsher penalty.
//' @param nboot Integer specifying the number of bootstrap replicates for imputation (default = 1). If > 1, enables bootstrapping.
//' @param seed Integer seed for random number generation during bootstrapping (default = 42). Only used when `nboot > 1`.
//' @param cores Number of CPU cores to use for parallel processing (default = 1).
//' @return A matrix where the first column is the 1-based row index, the second column is the 1-based column index,
//' and the subsequent `nboot` columns contain the imputed values (one column per bootstrap replicate if `nboot > 1`).
//'
//' @export
// [[Rcpp::export]]
arma::mat impute_knn_brute(
    const arma::mat &obj,         // Input data matrix with missing values as NaN
    const arma::umat &miss,       // Matrix of same size as obj, with 1 for missing (NA) and 0 for present
    const arma::uword k,          // Number of nearest neighbors to use
    const arma::uvec &n_col_miss, // Vector containing the count of missing values per column
    const int method,             // Distance metric: 0=Euclidean, 1=Manhattan, 2=impute.knn's method
    bool weighted,                // TRUE for a weighted average, FALSE for a simple average
    const double dist_pow,        // Power for distance penalty in weighted average
    const arma::uword nboot = 1,  // Number of boot strap ( > 1)
    const arma::uword seed = 42,  // seed
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
  case 2:
    calc_dist = calc_distance_knn;
    break;
  default:
    throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan, 2=impute.knn");
  }
  // column number = 2 (row index, column index) + nboot (each boot is an additional column)
  const arma::uword n_col_result = 2 + nboot;
  // Find columns that contain missing values
  arma::uvec col_index_miss = arma::find(n_col_miss > 0);
  if (col_index_miss.n_elem == 0)
  {
    // Return empty matrix if no missing values exist
    return arma::mat(0, n_col_result);
  }
  // For finding neighbors, we need all columns. Separate them into missing and non-missing.
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
  // Initialize the result matrix: [row index, column index, imputed values...]
  // This allows for easy merging back into a matrix in R.
  arma::uword sum_missing = arma::accu(n_col_miss);
  arma::mat result(sum_missing, n_col_result);
  result.fill(arma::datum::nan);
  // To efficiently fill the 'result' matrix, we calculate offsets.
  // 'col_offsets(i)' will be the starting row in 'result' for the i-th column with missing data.
  arma::uvec miss_counts = n_col_miss.elem(col_index_miss);
  arma::uvec col_offsets(miss_counts.n_elem + 1);
  col_offsets.fill(arma::fill::zeros);
  col_offsets.subvec(1, miss_counts.n_elem) = arma::cumsum(miss_counts);
  // Pre-fill result with row and column indices to avoid "continue" skipping the indices
  arma::uword offset = 0;
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    const arma::uword target_col_idx = col_index_miss(i);
    const arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
      const arma::uword row_idx = rows_to_impute(r);
      const arma::uword res_row = offset + r;
      result(res_row, 0) = row_idx + 1;        // R Row index (1-based)
      result(res_row, 1) = target_col_idx + 1; // R Column index (1-based)
    }
    offset += rows_to_impute.n_elem;
  }
  // Setup bootstrap stuff
  if (nboot > 1)
  {
    weighted = false;
  }
  // Main imputation loop: iterate through only the columns that need imputation
#if defined(_OPENMP)
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
    for (arma::uword j = 0; j < n_neighbors; ++j)
    {
      nn_columns(j) = neighbor_index(top_k[j].index);
    }
    // Init the boot strap neighbors and distance matrix
    arma::umat nn_columns_mat(n_neighbors, nboot);
    arma::mat nn_weights_mat(n_neighbors, nboot);
    arma::uword target_col_idx = col_index_miss(i);

    if (nboot > 1)
    {
      std::mt19937 gen(seed + target_col_idx);
      std::uniform_int_distribution<arma::uword> dist(0, n_neighbors - 1);
      for (arma::uword b = 0; b < nboot; ++b)
      {
        arma::uvec resampled(n_neighbors);
        for (arma::uword j = 0; j < n_neighbors; ++j)
        {
          arma::uword idx = dist(gen);
          resampled(j) = nn_columns(idx);
        }
        nn_columns_mat.col(b) = resampled;
      }
      nn_weights_mat.fill(1.0);
    }
    else
    {
      nn_columns_mat.col(0) = nn_columns;
      if (weighted)
      {
        // Extract distances for weight calculation
        arma::vec nn_dists(n_neighbors);
        for (arma::uword j = 0; j < n_neighbors; ++j)
        {
          nn_dists(j) = top_k[j].distance;
        }
        nn_weights_mat.col(0) = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
      }
      else
      {
        nn_weights_mat.col(0).fill(1.0);
      }
    }
    // Get the original index of the column we are currently imputing
    // Find which rows are missing in this specific column
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
    // For each missing cell in this column, calculate the imputed value
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
      arma::uword row_idx = rows_to_impute(r); // The actual row index of the missing value
      // Find the correct row in the result matrix to write to
      const arma::uword res_row = col_offsets(i) + r;
      for (arma::uword b = 0; b < nboot; ++b)
      {
        double weighted_sum = 0.0;
        double weight_total = 0.0;
        // Aggregate values from neighbors
        for (arma::uword j = 0; j < n_neighbors; ++j)
        {
          arma::uword neighbor_col_idx = nn_columns_mat(j, b);
          // A neighbor can only contribute if its value in the same row is NOT missing
          if (miss(row_idx, neighbor_col_idx) == 0)
          {
            double weight = nn_weights_mat(j, b);
            weighted_sum += weight * obj(row_idx, neighbor_col_idx);
            weight_total += weight;
          }
        }
        double imputed_value = (weight_total > 0.0) ? (weighted_sum / weight_total) : arma::datum::nan;
        result(res_row, 2 + b) = imputed_value; // boot starts at column 2
      }
    }
  }
  return result;
}
