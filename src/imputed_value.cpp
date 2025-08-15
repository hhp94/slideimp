#include "imputed_value.h"
#include <random>
#include <cmath>
#include <algorithm>

// initialize matrix to store the result. Also modify col_offsets
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
    const arma::vec &nn_weights,
    const arma::uword n_imp)
{
    // Find which rows are missing in this specific column
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));
    // For each missing cell in this column, calculate the imputed value
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
        const arma::uword row_idx = rows_to_impute(r);
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
                    double weight = nn_weights(j);
                    weighted_sum += weight * obj(row_idx, neighbor_col_idx);
                    weight_total += weight;
                }
            }
            double imputed_value = (weight_total > 0.0) ? (weighted_sum / weight_total) : arma::datum::nan;
            result(res_row, 2 + b) = imputed_value; // Bootstrap values start at column 2
        }
    }
}

void resample_neighbor(
    arma::umat &nn_columns_mat,
    const arma::uword seed,
    const arma::uword target_col_idx)
{
    // if (nn_columns_mat.n_cols <= 1)
    //     return;
    std::mt19937 gen(seed + target_col_idx);

    const arma::uword n_neighbors = nn_columns_mat.n_rows;
    const arma::uvec orig_columns = nn_columns_mat.col(0);

    std::uniform_int_distribution<arma::uword> dist(0, n_neighbors - 1);
    // sampling with replacement the neighbor matrix
    for (arma::uword b = 0; b < nn_columns_mat.n_cols; ++b)
    {
        for (arma::uword j = 0; j < n_neighbors; ++j)
        {
            const arma::uword rnd_idx = dist(gen);
            nn_columns_mat(j, b) = orig_columns(rnd_idx);
        }
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
    std::vector<const double*> data_cols(n_neighbors);
    std::vector<const arma::uword*> miss_cols(n_neighbors);
    for (arma::uword j = 0; j < n_neighbors; ++j) {
        data_cols[j] = obj.colptr(nn_columns(j));
        miss_cols[j] = miss.colptr(nn_columns(j));
    }
    const double* weight_ptr = nn_weights.memptr();
    // Process each row
    for (arma::uword r = 0; r < n_rows; ++r) {
        double numerator = 0.0;
        double denominator = 0.0;
        for (arma::uword j = 0; j < n_neighbors; ++j) {
            if (miss_cols[j][r] == 0) {
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

//' @title k-NN Impute using Predictive Mean Matching
//'
//' @description Rcpp implementation of predictive mean matching (PMM) for imputation. Exported for testing only.
//'
//' @param result Result matrix where column 2:2+n_imp contains imputed values (modified in place).
//' @param obj R matrix
//' @param miss is.na(obj).
//' @param col_offset Basically `cumsum(colSums(is.na(obj)))`.
//' @param target_col_idx The index of the target column to impute in obj.
//' @param nn_columns Column idx of nearest neighbors.
//' @param nn_weights Weights corresponding to nn_columns.
//' @param n_imp The number of imputations to perform for each missing value.
//' @param n_pmm The number of closest donors to consider for PMM.
//' @param seed RNG seed.
//'
//' @return void (modifies the result matrix in place).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
void impute_column_values_pmm(
    arma::mat &result,
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword col_offset,
    const arma::uword target_col_idx,
    const arma::uvec &nn_columns,
    const arma::vec &nn_weights,
    const arma::uword n_imp,
    const arma::uword n_pmm,
    const arma::uword seed)
{
    // Find rows with observed values (donors)
    arma::uvec donor_rows = arma::find(miss.col(target_col_idx) == 0);

    // Calculate predicted values for all rows
    arma::vec predicted = weighted_row_means(obj, miss, nn_columns, nn_weights);

    // Find rows that need imputation in the target column
    arma::uvec rows_to_impute = arma::find(miss.col(target_col_idx));

    // Setup random number generator
    std::mt19937 gen(seed + target_col_idx);

    // For each missing value
    for (arma::uword r = 0; r < rows_to_impute.n_elem; ++r)
    {
        const arma::uword row_idx = rows_to_impute(r);
        const arma::uword res_row = col_offset + r;
        // target_pred is the target value to calculate distance from
        const double target_pred = predicted(row_idx);
        // Skip if predicted value is NaN
        if (!std::isfinite(target_pred))
        {
            continue;
        }
        // Calculate distances between target predicted value and donor predicted values
        arma::vec donor_distances(donor_rows.n_elem);
        for (arma::uword d = 0; d < donor_rows.n_elem; ++d)
        {
            donor_distances(d) = std::abs(predicted(donor_rows(d)) - target_pred);
        }
        // Find the n_pmm closest donors
        arma::uword n_closest = std::min(n_pmm, donor_rows.n_elem);
        // Create indices for sorting
        arma::uvec indices = arma::linspace<arma::uvec>(0, donor_rows.n_elem - 1, donor_rows.n_elem);
        // Use nth_element to partially sort and find the n_closest donors
        if (n_closest < donor_rows.n_elem)
        {
            std::nth_element(
                indices.begin(),
                indices.begin() + n_closest,
                indices.end(),
                [&donor_distances](arma::uword i, arma::uword j)
                {
                    return donor_distances(i) < donor_distances(j);
                });
        }
        // Then for each imputation
        for (arma::uword b = 0; b < n_imp; ++b)
        {
            // Randomly select from the n_closest donors
            std::uniform_int_distribution<arma::uword> dist(0, n_closest - 1);
            arma::uword selected_idx = indices(dist(gen));
            arma::uword selected_donor_row = donor_rows(selected_idx);
            // Use the observed value from the selected donor
            result(res_row, 2 + b) = obj(selected_donor_row, target_col_idx);
        }
    }
}
