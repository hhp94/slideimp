#ifndef IMPUTED_VALUE_H
#define IMPUTED_VALUE_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

constexpr double epsilon = 1e-10;

// Helper function to initialize the result matrix with row and column indices.
// Returns the result matrix and populates col_offsets and rows_to_impute_vec.
arma::mat initialize_result_matrix(
    const arma::mat &nmiss,
    const arma::uvec &col_index_miss,
    arma::uvec &col_offsets,
    std::vector<arma::uvec> &rows_to_impute_vec);

// Helper function to impute missing values for a single column
void impute_column_values(
    arma::mat &result,
    const arma::mat &obj,
    const arma::mat &nmiss,
    const arma::uword col_offset,
    const arma::uword target_col_idx,
    const arma::uvec &nn_columns_vec,
    const arma::vec &nn_weights,
    const arma::uvec &rows_to_impute);

#endif
