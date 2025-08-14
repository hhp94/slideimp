#ifndef IMPUTED_VALUE_H
#define IMPUTED_VALUE_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

constexpr double epsilon = 1e-10;

// Helper function to initialize the result matrix with row and column indices
// Returns the result matrix and populates col_offsets for later use
arma::mat initialize_result_matrix(
    const arma::umat &miss,
    const arma::uvec &col_index_miss,
    const arma::uvec &n_col_miss,
    const arma::uword n_imp,
    arma::uvec &col_offsets);

// Helper function to impute missing values for a single column
void impute_column_values(
    arma::mat &result,
    const arma::mat &obj,
    const arma::umat &miss,
    const arma::uword col_offset,
    const arma::uword target_col_idx,
    const arma::umat &nn_columns_mat,
    const arma::mat &nn_weights_mat,
    const arma::uword n_imp);

// Helper function to prepare neighbor columns and weights for bootstrapping
void prepare_bootstrap_neighbors(
    arma::umat &nn_columns_mat,
    arma::mat &nn_weights_mat,
    const arma::uvec &nn_columns,
    const arma::vec &nn_dists,
    const arma::uword n_imp,
    const arma::uword seed,
    const arma::uword target_col_idx,
    const bool weighted,
    const double dist_pow);

#endif
