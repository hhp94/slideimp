#include "imputed_value.h"
#include <mlpack.h>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <stdexcept>
#if defined(_OPENMP)
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat impute_knn_mlpack(
    const arma::mat &obj,         // data with NA pre-filled (no NA values)
    const arma::mat &nmiss,       // missing data matrix
    const arma::uword k,          // number of neighbors
    const arma::uvec &grp_impute, // 0-based indices of columns to impute
    const int method,             // 0 = euclidean, 1 = manhattan
    const double dist_pow,        // controls distance penalty for weighted average
    const int cores = 1)          // number of cores for parallel processing
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

  arma::uvec col_offsets;
  std::vector<arma::uvec> rows_to_impute_vec;
  arma::mat result = initialize_result_matrix(nmiss, grp_impute, col_offsets, rows_to_impute_vec);
  if (result.n_rows == 0)
  {
    return result;
  }

  arma::mat query_mat = obj.cols(grp_impute);

  // Matrices to store output
  arma::umat resultingNeighbors;
  arma::mat resultingDistances;

  if (method == 0)
  {
    using KNNType = mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::EuclideanDistance, arma::mat, mlpack::BallTree>;
    KNNType nn(obj);
    nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
  }
  else if (method == 1)
  {
    using KNNType = mlpack::NeighborSearch<mlpack::NearestNeighborSort, mlpack::ManhattanDistance, arma::mat, mlpack::BallTree>;
    KNNType nn(obj);
    nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
  }
  else
  {
    throw std::invalid_argument("Invalid `method`. Use 0 for 'euclidean' or 1 for 'manhattan'.");
  }

  // Main imputation loop
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (arma::uword i = 0; i < grp_impute.n_elem; ++i)
  {
    // Get the indices of the k nearest neighbors (skip the 0-th neighbor, which is self)
    arma::uvec nn_columns = resultingNeighbors(arma::span(1, k), i);

    // Get the corresponding distances for these neighbors
    arma::vec nn_dists = resultingDistances(arma::span(1, k), i);
    const arma::uword n_neighbors = nn_columns.n_elem;

    // Calculate weights once (same for all methods)
    arma::vec weights(n_neighbors);
    weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);

    arma::uword target_col_idx = grp_impute(i);

    // Single deterministic imputation
    impute_column_values(
        result, obj, nmiss,
        col_offsets(i), target_col_idx,
        nn_columns, weights,
        rows_to_impute_vec[i]);
  }

  return result;
}
