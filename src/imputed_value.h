#ifndef IMPUTED_VALUE_H
#define IMPUTED_VALUE_H

#include <RcppArmadillo.h>
#include <cstdint>
#include <vector>
#include <cmath>           // std::isinf (used by stop_on_inf)

// single source of truth for the mask storage type
using mask_t = uint8_t;
using MaskMat = arma::Mat<mask_t>;

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

// n_col_valid: per-column valid-entry counts for the masked region
// (groups 1+2). Populated by the caller during the copy-with-mask pass,
// so initialize_result_matrix can derive missing counts without a second
// accu over the byte matrix. Only entries [0, layout.n_imp) are read here.
arma::mat initialize_result_matrix(
    const MaskMat &nmiss_masked,
    const arma::uvec &grp_impute,
    const GroupLayout &layout,
    const arma::uvec &n_col_valid,
    arma::uvec &col_offsets,
    std::vector<arma::uvec> &rows_to_impute_vec);

// nn_columns holds tagged positions: values < layout.complete_start() index
// directly into obj_masked/nmiss_masked; values >= layout.complete_start()
// must be de-tagged via (nn_columns(j) - complete_start()) and looked up
// through grp_complete into obj.
void impute_column_values(
    arma::mat &result,
    const arma::mat &obj_masked,
    const MaskMat &nmiss_masked,
    const GroupLayout &layout,
    const arma::uword col_offset,
    const arma::uvec &nn_columns,
    const arma::vec &nn_weights,
    const arma::uvec &rows_to_impute,
    const arma::mat &obj,
    const arma::uvec &grp_complete);
#endif

static inline void stop_on_inf(const arma::mat &obj)
{
    const arma::uword n_rows = obj.n_rows;
    const arma::uword n_cols = obj.n_cols;

    for (arma::uword c = 0; c < n_cols; ++c)
    {
        const double *col = obj.colptr(c);
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            if (std::isinf(col[r]))
            {
                Rcpp::stop("Infinite value found at row %u, column %u. "
                           "Infinite values are not supported.",
                           r + 1, c + 1);
            }
        }
    }
}
