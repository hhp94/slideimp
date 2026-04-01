#include "imputed_value.h"
#include <random>
#include <cmath>
#include <algorithm>

// initialize matrix to store the result. Also modify col_offsets
arma::mat initialize_result_matrix(
    const arma::Mat<unsigned char> &miss,
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
            result(res_row, 0) = row_idx + 1; // R Row index (1-based)
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
    const arma::Mat<unsigned char> &miss,
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
