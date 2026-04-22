#ifndef IMPUTED_VALUE_H
#define IMPUTED_VALUE_H

#include <RcppArmadillo.h>
#include <vector>

constexpr double epsilon = 1e-10;

// GroupLayout: single source of truth for group boundaries.
// Groups 1+2 (imputed + missing-no-impute) share identical kernel treatment
// and live together in obj_masked / nmiss_masked at local cols [0, n_masked()).
// n_imp is retained so initialize_result_matrix knows which columns actually
// get imputed.
//
// complete_start() is a "virtual" index: there is no physical column at that
// position in obj_masked. It's used as a tag in nn_columns to mark
// "this neighbor came from group 3, offset (local - complete_start()) into
// grp_complete".
struct GroupLayout
{
    arma::uword n_imp;
    arma::uword n_mni;
    arma::uword n_complete;
    arma::uword mni_start() const { return n_imp; }
    arma::uword n_masked() const { return n_imp + n_mni; }
    arma::uword complete_start() const { return n_imp + n_mni; }
    arma::uword n_working() const { return n_imp + n_mni + n_complete; }
};

arma::mat initialize_result_matrix(
    const arma::mat &nmiss_masked,
    const arma::uvec &grp_impute,
    const GroupLayout &layout,
    arma::uvec &col_offsets,
    std::vector<arma::uvec> &rows_to_impute_vec);

// nn_columns holds tagged positions: values < layout.complete_start() index
// directly into obj_masked/nmiss_masked; values >= layout.complete_start()
// must be de-tagged via (nn_columns(j) - complete_start()) and looked up
// through grp_complete into obj.
void impute_column_values(
    arma::mat &result,
    const arma::mat &obj_masked,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword col_offset,
    const arma::uvec &nn_columns,
    const arma::vec &nn_weights,
    const arma::uvec &rows_to_impute,
    const arma::mat &obj,
    const arma::uvec &grp_complete);

#endif
