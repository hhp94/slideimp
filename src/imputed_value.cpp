#include "imputed_value.h"
#include <random>
#include <cmath>

arma::mat initialize_result_matrix(
    const arma::umat &miss,
    const arma::uvec &col_index_miss,
    const arma::uvec &n_col_miss,
    const arma::uword n_imp,
    arma::uvec &col_offsets)
{
    // 2 columns for ij indices + n_imp columns for values
    const arma::uword n_col_result = 2 + n_imp;

    // Find columns that contain missing values
    if (col_index_miss.n_elem == 0)
    {
        // Return empty matrix if no missing values exist
        return arma::mat(0, n_col_result);
    }

    // Calculate total number of missing values
    arma::uword sum_missing = arma::accu(n_col_miss);
    arma::mat result(sum_missing, n_col_result);
    result.fill(arma::datum::nan);

    // Calculate offsets for efficient filling
    arma::uvec miss_counts = n_col_miss.elem(col_index_miss);
    col_offsets.set_size(miss_counts.n_elem + 1);
    col_offsets.fill(arma::fill::zeros);
    col_offsets.subvec(1, miss_counts.n_elem) = arma::cumsum(miss_counts);

    // Pre-fill result with row and column indices (1-based for R)
    for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
    {
        const arma::uword target_col_idx = col_index_miss(i);
        const arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));

        for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
        {
            const arma::uword row_idx = rows_to_impute(r);
            const arma::uword res_row = col_offsets(i) + r;
            result(res_row, 0) = row_idx + 1;        // R Row index (1-based)
            result(res_row, 1) = target_col_idx + 1; // R Column index (1-based)
        }
    }

    return result;
}

void impute_column_values(
    arma::mat &result,
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword col_offset,
    const arma::uword target_col_idx,
    const arma::umat &nn_columns_mat,
    const arma::mat &nn_weights_mat,
    const arma::uword n_imp)
{
    // Find which rows are missing in this specific column
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));

    // For each missing cell in this column, calculate the imputed value
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
        arma::uword row_idx = rows_to_impute(r);
        const arma::uword res_row = col_offset + r;

        for (arma::uword b = 0; b < n_imp; ++b)
        {
            double weighted_sum = 0.0;
            double weight_total = 0.0;

            // Aggregate values from neighbors
            for (arma::uword j = 0; j < nn_columns_mat.n_rows; ++j)
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
            result(res_row, 2 + b) = imputed_value; // Bootstrap values start at column 2
        }
    }
}

void prepare_bootstrap_neighbors(
    arma::umat &nn_columns_mat,
    arma::mat &nn_weights_mat,
    const arma::uvec &nn_columns,
    const arma::vec &nn_dists,
    const arma::uword n_imp,
    const arma::uword seed,
    const arma::uword target_col_idx,
    const bool weighted,
    const double dist_pow)
{
    arma::uword n_neighbors = nn_columns.n_elem;

    // Initialize matrices
    nn_columns_mat.set_size(n_neighbors, n_imp);
    nn_weights_mat.set_size(n_neighbors, n_imp);

    if (n_imp > 1)
    {
        // Setup thread-safe random number generator for bootstrapping
        std::mt19937 gen(seed + target_col_idx);
        std::uniform_int_distribution<arma::uword> dist(0, n_neighbors - 1);

        for (arma::uword b = 0; b < n_imp; ++b)
        {
            arma::uvec resampled(n_neighbors);
            for (arma::uword j = 0; j < n_neighbors; ++j)
            {
                arma::uword idx = dist(gen);
                resampled(j) = nn_columns(idx);
            }
            nn_columns_mat.col(b) = resampled;
        }
        // For bootstrap, weights are all 1 (not weighted)
        nn_weights_mat.fill(1.0);
    }
    else
    {
        nn_columns_mat.col(0) = nn_columns;

        // Compute weights based on distances
        if (weighted)
        {
            nn_weights_mat.col(0) = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);
        }
        else
        {
            nn_weights_mat.col(0).fill(1.0); // Simple average
        }
    }
}
