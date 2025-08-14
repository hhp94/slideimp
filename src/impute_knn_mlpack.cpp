// [[Rcpp::depends(mlpack, RcppArmadillo)]]

#include "imputed_value.h"
#include <mlpack.h>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <RcppArmadillo.h>
#include <stdexcept>
#if defined(_OPENMP)
#include <omp.h>
#endif

//' Impute missing values in a matrix using treed k-nearest neighbors (k-NN)
//'
//' k-NN using KDTree or BallTree with optional bootstrap support for uncertainty estimation.
//'
//' When `n_imp > 1`, bootstrapping is enabled: for each missing value, `n_imp` imputed values are generated
//' by resampling the k nearest neighbors with replacement and using simple averages (weighted is forced to FALSE).
//' This provides variability estimates for imputation uncertainty.
//'
//' @param obj Numeric matrix with missing values pre-filled with colMeans.
//' @param miss Logical matrix (0/1) indicating missing values (1 = missing).
//' @param k Number of nearest neighbors to use for imputation.
//' @param n_col_miss Integer vector specifying the count of missing values per column.
//' @param method Integer specifying the distance metric: 0 = Euclidean, 1 = Manhattan.
//' @param tree Which type of tree? "kd" or "ball".
//' @param weighted Boolean controls for the imputed value to be a simple mean or weighted mean by inverse distance.
//'   Note: Forced to FALSE when `n_imp > 1`.
//' @param dist_pow A positive double that controls the penalty for larger distances in
//' the weighted mean imputation. Must be greater than zero: values between 0 and 1 apply a softer penalty,
//' 1 is linear (default), and values greater than 1 apply a harsher penalty.
//' @param n_imp Integer specifying the number of bootstrap replicates for imputation (default = 1). If > 1, enables bootstrapping.
//' @param n_pmm Integer specifying the number of donors for pmm.
//' @param seed Integer seed for random number generation during bootstrapping (default = 42). Only used when `n_imp > 1`.
//' @param cores Number of CPU cores to use for parallel processing (default = 1).
//' @return A matrix where the first column is the 1-based row index, the second column is the 1-based column index,
//' and the subsequent `n_imp` columns contain the imputed values (one column per bootstrap replicate if `n_imp > 1`).
//'
//' @export
// [[Rcpp::export]]
arma::mat impute_knn_mlpack(
    const arma::mat &obj,         // data with NA pre-filled. So there's no NA
    const arma::umat &miss,       // missing data matrix
    const arma::uword k,          // n neighbors
    const arma::uvec &n_col_miss, // vector of missing per column
    const int method,             // 0 = "euclidean" or 1 = "manhattan"
    const std::string tree,       // "kd" or "ball"
    bool weighted,                // weighted average for imputation or not
    const double dist_pow,        // controls distance penalty for weighted average
    const arma::uword n_imp = 1,  // Number of imputation
    arma::uword n_pmm = 0,        // Number of pmm donors. If n_imp = 1, then n_pmm = 0 else if n_pmm > 0, then use pmm, else use bootstrap
    const arma::uword seed = 42,  // seed for RNG
    const int cores = 1)          // Number of cores for parallel processing
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  // Find columns that contain missing values
  arma::uvec col_index_miss = arma::find(n_col_miss > 0);
  // Initialize result matrix and get column offsets using helper function
  arma::uvec col_offsets;
  arma::mat result = initialize_result_matrix(miss, col_index_miss, n_col_miss, n_imp, col_offsets);
  if (result.n_rows == 0)
  {
    return result; // No missing values
  }
  arma::mat query_mat = obj.cols(col_index_miss);
  // Matrices to store output
  arma::umat resultingNeighbors;
  arma::mat resultingDistances;
  // Perform k-NN search based on tree and method
  if (tree == "kd" && method == 0)
  {
    using KNNType = mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::KDTree>;
    KNNType nn(obj);
    nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
  }
  else if (tree == "kd" && method == 1)
  {
    using KNNType = mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::ManhattanDistance, arma::mat, mlpack::KDTree>;
    KNNType nn(obj);
    nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
  }
  else if (tree == "ball" && method == 0)
  {
    using KNNType = mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::BallTree>;
    KNNType nn(obj);
    nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
  }
  else if (tree == "ball" && method == 1)
  {
    using KNNType = mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::ManhattanDistance, arma::mat, mlpack::BallTree>;
    KNNType nn(obj);
    nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
  }
  else
  {
    throw std::invalid_argument("Invalid `tree` or `method`. Use 'kd' or 'ball' for `tree`, and 0 for 'euclidean' or 1 for 'manhattan' for `method`.");
  }
  if (n_imp > 1 && n_pmm == 0)
  {
    // if n_imp > 1, and n_pmm == 0, then it's bootstrapping. In which case
    // we fix weighted to false
    weighted = false;
  }
  // Main imputation loop
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
  {
    // Get the indices of the k nearest neighbors (skip the 0-th neighbor, which is self)
    arma::uvec nn_columns = resultingNeighbors(arma::span(1, k), i);
    // Get the corresponding distances for these neighbors
    arma::vec nn_dists = resultingDistances(arma::span(1, k), i);
    // Copy pasted from impute_knn_brute
    const arma::uword n_neighbors = nn_columns.n_elem;
    // Calculate weights once (same for all methods)
    arma::vec weights(n_neighbors);
    if (weighted)
    {
      weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
    }
    else
    {
      weights.fill(1.0);
    }
    arma::uword target_col_idx = col_index_miss(i);
    // Single Imputation
    if (n_imp == 1)
    {
      arma::umat nn_columns_mat(n_neighbors, 1);
      nn_columns_mat.col(0) = nn_columns;

      impute_column_values(
          result, obj, miss,
          col_offsets(i), target_col_idx,
          nn_columns_mat, weights,
          n_imp);
    }
    else if (n_imp > 1 && n_pmm > 0)
    {
      // PMM
      impute_column_values_pmm(
          result, obj, miss,
          col_offsets(i), target_col_idx,
          nn_columns, weights,
          n_imp, n_pmm, seed);
    }
    else if (n_imp > 1 && n_pmm == 0)
    {
      // Neighbor Boot strap
      arma::umat nn_columns_mat(n_neighbors, n_imp);
      for (arma::uword b = 0; b < n_imp; ++b)
      {
        nn_columns_mat.col(b) = nn_columns;
      }
      resample_neighbor(nn_columns_mat, seed, target_col_idx);
      impute_column_values(
          result, obj, miss,
          col_offsets(i), target_col_idx,
          nn_columns_mat, weights,
          n_imp);
    }
  }

  return result;
}
