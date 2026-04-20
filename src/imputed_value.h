#ifndef IMPUTED_VALUE_H
#define IMPUTED_VALUE_H

#include <RcppArmadillo.h>
#include <vector>

constexpr double epsilon = 1e-10;

// GroupLayout: single source of truth for group boundaries in the
// reordered working matrices. All local positions are indices into
// obj_reordered / nmiss_masked, NOT into the original obj.
struct GroupLayout
{
    arma::uword n_imp;
    arma::uword n_mni;
    arma::uword n_complete;

    arma::uword mni_start() const { return n_imp; }
    arma::uword complete_start() const { return n_imp + n_mni; }
    arma::uword n_working() const { return n_imp + n_mni + n_complete; }
};

// nmiss_masked covers only groups 1+2 (local cols [0, n_imp + n_mni)).
// grp_impute is needed only for writing original 1-based column indices
// into the result matrix — after this function returns, the original
// coordinate space is not touched again.
arma::mat initialize_result_matrix(
    const arma::mat &nmiss_masked,
    const arma::uvec &grp_impute,
    const GroupLayout &layout,
    arma::uvec &col_offsets,
    std::vector<arma::uvec> &rows_to_impute_vec);

// nn_columns holds LOCAL positions in obj_reordered. Neighbors with
// local_pos < layout.complete_start() need the mask; others don't.
void impute_column_values(
    arma::mat &result,
    const arma::mat &obj_reordered,
    const arma::mat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword col_offset,
    const arma::uvec &nn_columns,
    const arma::vec &nn_weights,
    const arma::uvec &rows_to_impute);

#endif
