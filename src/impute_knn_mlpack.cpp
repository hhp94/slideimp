// [[Rcpp::depends(mlpack, RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <stdexcept>
#include <cmath>
#if defined(_OPENMP)
#include <omp.h>
#endif

constexpr double epsilon = 1e-10;


//' Impute missing values in a matrix using treed k-nearest neighbors (k-NN)
//'
//' k-NN using KDTree or BallTree.
//'
//' @param obj Numeric matrix with missing values pre-filled with colMeans.
//' @param miss Logical matrix (0/1) indicating missing values (1 = missing).
//' @param k Number of nearest neighbors to use for imputation.
//' @param n_col_miss Integer vector specifying the count of missing values per column.
//' @param method Integer specifying the distance metric: 0 = Euclidean, 1 = Manhattan.
//' @param tree Which type of tree? "kd" or "ball".
//' @param weighted Boolean controls for the imputed value to be a simple mean or weighted mean by inverse distance.
//' @param dist_pow A positive double that controls the penalty for larger distances in
//' the weighted mean imputation. Must be greater than zero: values between 0 and 1 apply a softer penalty,
//' 1 is linear (default), and values greater than 1 apply a harsher penalty.
//' @param cores Number of CPU cores to use for parallel processing (default = 1).
//' @return A matrix with imputed values where missing values were present.
//'
//' @export
// [[Rcpp::export]]
arma::mat impute_knn_mlpack(
    const arma::mat &obj,        // data with NA pre-filled. So there's no NA
    const arma::umat &miss,      // missing data matrix
    const arma::uword k,         // n neighbors
    const arma::uvec n_col_miss, // vector of missing per column
    const int method,            // 0 = "euclidean" or 1 = "manhattan"
    const std::string tree,      // "kd" or "ball"
    const bool weighted,         // weighted average for imputation or not
    const double dist_pow,       // controls distance penalty for weighted average
    const int cores)
{
#if defined(_OPENMP)
  omp_set_num_threads(cores);
#else
  cores = 1;
#endif
    // col_index_miss holds columns that has any missing
    arma::uvec col_index_miss = arma::find(n_col_miss);
    if (col_index_miss.n_elem == 0)
    {
        return obj;
    }
    arma::mat query_mat = obj.cols(col_index_miss);
    // Matrices to store output
    arma::umat resultingNeighbors;
    arma::mat resultingDistances;
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

    // Nearest neighbors and distance are now stored. We proceed to imputation column wise.
    arma::mat imputed = obj;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
    {
        // Get the indices of the k nearest neighbors for the i-th query point (skip the 0-th neighbor, which is the point itself)
        arma::uvec nn_columns = resultingNeighbors(arma::span(1, k), i);
        // Get the corresponding distances for these neighbors
        arma::vec nn_dists = resultingDistances(arma::span(1, k), i);
        arma::vec nn_weights(nn_columns.n_elem);
        if (weighted)
        {
            nn_weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
        }
        else
        {
            nn_weights.fill(1.0); // weight = 1 is just simple average
        }
        // now we move on to nn imputation. Find the rows to impute for this column
        arma::uword target_col_idx = col_index_miss(i);
        arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
        // man, range-for-loop is trippy. kind of not used to this but it works
        for (arma::uword row_idx : rows_to_impute)
        {
            double weighted_sum = 0.0;
            double weight_total = 0.0;
            for (arma::uword j = 0; j < nn_columns.n_elem; ++j)
            {
                arma::uword neighbor_col_idx = nn_columns(j);
                if (miss(row_idx, neighbor_col_idx) == 0)
                {
                    double weight = nn_weights(j); // Use precomputed weight
                    weighted_sum += weight * obj(row_idx, neighbor_col_idx);
                    weight_total += weight;
                }
            }
            if (weight_total > 0.0)
            {
                imputed(row_idx, target_col_idx) = weighted_sum / weight_total;
            }
        }
    }
    return imputed;
}