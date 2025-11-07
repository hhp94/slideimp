#include "imputed_value.h"
#include <random>
#include <cmath>
#include <algorithm>

// initialize matrix to store the result. Also modify col_offsets
arma::mat initialize_result_matrix(
    const arma::umat &miss,
    const arma::uvec &col_index_miss,
    const arma::uvec &n_col_miss,
    arma::uvec &col_offsets)
{
    // 2 columns for ij indices + n_imp columns for values. There's only one imputation
    const arma::uword n_col_result = 3;

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
    // col 0 offset = 0, col 1 -> N offset = cumsum number of missing in each col
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

// basically the imputation for loop down each column. Modify the result matrix.
// group the bootstrap and single imputation together because for these we only
// have to work on rows with missing.
void impute_column_values(
    arma::mat &result,
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword col_offset,
    const arma::uword target_col_idx,
    const arma::umat &nn_columns_mat,
    const arma::vec &nn_weights)
{
    // Find which rows are missing in this specific column
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
    // For each missing cell in this column, calculate the imputed value
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
        const arma::uword row_idx = rows_to_impute(r);
        const arma::uword res_row = col_offset + r;
        double weighted_sum = 0.0;
        double weight_total = 0.0;
        // Aggregate values from neighbors
        for (arma::uword j = 0; j < nn_columns_mat.n_rows; ++j)
        {
            arma::uword neighbor_col_idx = nn_columns_mat(j, 0);
            // A neighbor can only contribute if its value in the same row is NOT missing
            if (miss(row_idx, neighbor_col_idx) == 0)
            {
                double weight = nn_weights(j);
                weighted_sum += weight * obj(row_idx, neighbor_col_idx);
                weight_total += weight;
            }
        }
        double imputed_value = (weight_total > 0.0) ? (weighted_sum / weight_total) : arma::datum::nan;
        result(res_row, 2) = imputed_value;
    }
}

//' @title Weighted Row Mean
//'
//' @description Calculate weighted row means for specified columns, accounting for missing values.
//'
//' @param obj A numeric matrix containing the data.
//' @param miss An unsigned integer matrix indicating missing values (0 for observed, 1 for missing).
//' @param nn_columns An unsigned integer vector of column indices for the neighbors.
//' @param nn_weights A numeric vector of weights corresponding to the neighbors.
//'
//' @return A column vector containing the weighted row means, with NaN where computation is not possible.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec weighted_row_means(
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uvec &nn_columns,
    const arma::vec &nn_weights)
{
    const arma::uword n_rows = obj.n_rows;
    const arma::uword n_neighbors = nn_columns.n_elem;
    arma::vec predicted(n_rows);
    // Pre-fetch column pointers to avoid repeated indexing
    std::vector<const double *> data_cols(n_neighbors);
    std::vector<const arma::uword *> miss_cols(n_neighbors);
    for (arma::uword j = 0; j < n_neighbors; ++j)
    {
        data_cols[j] = obj.colptr(nn_columns(j));
        miss_cols[j] = miss.colptr(nn_columns(j));
    }
    const double *weight_ptr = nn_weights.memptr();
    // Process each row
    for (arma::uword r = 0; r < n_rows; ++r)
    {
        double numerator = 0.0;
        double denominator = 0.0;
        for (arma::uword j = 0; j < n_neighbors; ++j)
        {
            if (miss_cols[j][r] == 0)
            {
                double weight = weight_ptr[j];
                numerator += data_cols[j][r] * weight;
                denominator += weight;
            }
        }
        // If denominator is zero then return NaN
        predicted(r) = (denominator > 0.0) ? (numerator / denominator) : arma::datum::nan;
    }

    return predicted;
}
