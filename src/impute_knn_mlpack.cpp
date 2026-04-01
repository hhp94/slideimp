#include "imputed_value.h"
#include <mlpack.h>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <RcppArmadillo.h>
#include <stdexcept>
#if defined(_OPENMP)
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat impute_knn_mlpack(
    const arma::mat &obj,         // data with NA pre-filled. So there's no NA
    const arma::mat &nmiss,       // missing data matrix
    const arma::uword k,          // n neighbors
    const arma::uvec &n_col_miss, // vector of missing per column
    const int method,             // 0 = "euclidean" or 1 = "manhattan"
    const std::string tree,       // "kd" or "ball"
    const double dist_pow,        // controls distance penalty for weighted average
    const int cores = 1)          // Number of cores for parallel processing
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  // Find columns that contain missing values
  arma::uvec col_index_miss = arma::find(n_col_miss > 0);
  // Initialize result matrix and get column offsets using helper function
  arma::uvec col_offsets;
  arma::mat result = initialize_result_matrix(nmiss, col_index_miss, n_col_miss, col_offsets);
  if (result.n_rows == 0)
  {
    return result;
  }
  arma::mat query_mat = obj.cols(col_index_miss);
  // Matrices to store output
  arma::umat resultingNeighbors;
  arma::mat resultingDistances;
  // Perform K-NN search based on tree and method
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
    weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
    arma::uword target_col_idx = col_index_miss(i);
    // Choose imputation method based on n_pmm
    // Single deterministic imputation
    arma::umat nn_columns_mat(n_neighbors, 1);
    nn_columns_mat.col(0) = nn_columns;

    impute_column_values(
        result, obj, nmiss,
        col_offsets(i), target_col_idx,
        nn_columns_mat, weights);
  }
  return result;
}
