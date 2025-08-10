// [[Rcpp::depends(mlpack, RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <stdexcept>
#include <cmath>
#include <random>
#if defined(_OPENMP)
#include <omp.h>
#endif

constexpr double epsilon = 1e-10;

//' Impute missing values in a matrix using treed k-nearest neighbors (k-NN)
//'
//' k-NN using KDTree or BallTree with optional bootstrap support for uncertainty estimation.
//'
//' When `nboot > 1`, bootstrapping is enabled: for each missing value, `nboot` imputed values are generated
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
//'   Note: Forced to FALSE when `nboot > 1`.
//' @param dist_pow A positive double that controls the penalty for larger distances in
//' the weighted mean imputation. Must be greater than zero: values between 0 and 1 apply a softer penalty,
//' 1 is linear (default), and values greater than 1 apply a harsher penalty.
//' @param nboot Integer specifying the number of bootstrap replicates for imputation (default = 1). If > 1, enables bootstrapping.
//' @param seed Integer seed for random number generation during bootstrapping (default = 42). Only used when `nboot > 1`.
//' @param cores Number of CPU cores to use for parallel processing (default = 1).
//' @return A matrix where the first column is the 1-based linear index of the missing value (as calculated in R),
//' and the subsequent `nboot` columns contain the imputed values (one column per bootstrap replicate if `nboot > 1`).
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
    const arma::uword nboot = 1,  // Number of bootstrap replicates
    const arma::uword seed = 42,  // seed for RNG
    const int cores = 1)          // Number of cores for parallel processing
{
#if defined(_OPENMP)
    omp_set_num_threads(cores);
#endif
    // col_index_miss holds columns that has any missing
    const arma::uword n_col_result = 1 + nboot;
    // Find columns that contain missing values
    arma::uvec col_index_miss = arma::find(n_col_miss > 0);
    if (col_index_miss.n_elem == 0)
    {
        // Return empty matrix if no missing values exist
        return arma::mat(0, n_col_result);
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
    // Setup bootstrap: force weighted to false if bootstrapping
    if (nboot > 1)
    {
        weighted = false;
    }
    // Nearest neighbors and distance are now stored. We proceed to imputation column wise.
    arma::uword sum_missing = arma::accu(n_col_miss);
    arma::mat result(sum_missing, n_col_result);
    result.fill(arma::datum::nan);

    arma::uvec miss_counts = n_col_miss.elem(col_index_miss);
    arma::uvec col_offsets(miss_counts.n_elem + 1);
    col_offsets.fill(arma::fill::zeros);
    col_offsets.subvec(1, miss_counts.n_elem) = arma::cumsum(miss_counts);

    // Pre-fill result with linear indices
    for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
    {
        const arma::uword target_col_idx = col_index_miss(i);
        const arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
        for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
        {
            const arma::uword row_idx = rows_to_impute(r);
            const arma::uword res_row = col_offsets(i) + r;
            result(res_row, 0) = target_col_idx * obj.n_rows + row_idx + 1;
        }
    }

#if defined(_OPENMP)
#pragma omp parallel for
#endif
    for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
    {
        // Get the indices of the k nearest neighbors for the i-th query point (skip the 0-th neighbor, which is self)
        arma::uvec nn_columns = resultingNeighbors(arma::span(1, k), i);
        // Get the corresponding distances for these neighbors
        arma::vec nn_dists = resultingDistances(arma::span(1, k), i);

        // Prepare storage for potentially bootstrapped neighbors and weights
        arma::umat nn_columns_mat(nn_columns.n_elem, nboot);
        arma::mat nn_weights_mat(nn_columns.n_elem, nboot);

        arma::uword target_col_idx = col_index_miss(i);

        if (nboot > 1)
        {
            // Setup thread-safe random number generator
            std::mt19937 gen(seed + target_col_idx);
            std::uniform_int_distribution<arma::uword> dist(0, nn_columns.n_elem - 1);
            for (arma::uword b = 0; b < nboot; ++b)
            {
                arma::uvec resampled(nn_columns.n_elem);
                for (arma::uword j = 0; j < nn_columns.n_elem; ++j)
                {
                    arma::uword idx = dist(gen);
                    resampled(j) = nn_columns(idx);
                }
                nn_columns_mat.col(b) = resampled;
            }
            // Since weighted = false for bootstrap, weights are all 1
            nn_weights_mat.fill(1.0);
        }
        else
        {
            nn_columns_mat.col(0) = nn_columns;
            // Precompute weights for each neighbor
            if (weighted)
            {
                nn_weights_mat.col(0) = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
            }
            else
            {
                nn_weights_mat.col(0).fill(1.0); // weight = 1 is just simple average
            }
        }
        // Find the rows to impute for this column
        arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
        for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
        {
            arma::uword row_idx = rows_to_impute(r);
            const arma::uword res_row = col_offsets(i) + r;

            for (arma::uword b = 0; b < nboot; ++b)
            {
                double weighted_sum = 0.0;
                double weight_total = 0.0;

                for (arma::uword j = 0; j < nn_columns_mat.n_rows; ++j)
                {
                    arma::uword neighbor_col_idx = nn_columns_mat(j, b);
                    if (miss(row_idx, neighbor_col_idx) == 0)
                    {
                        double weight = nn_weights_mat(j, b);
                        weighted_sum += weight * obj(row_idx, neighbor_col_idx);
                        weight_total += weight;
                    }
                }
                double imputed_value = (weight_total > 0.0) ? (weighted_sum / weight_total) : arma::datum::nan;
                result(res_row, 1 + b) = imputed_value;
            }
        }
    }
    return result;
}
