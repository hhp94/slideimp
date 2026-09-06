#ifndef IMPUTED_VALUE_H
#define IMPUTED_VALUE_H

#include <RcppArmadillo.h>
#include <cstdint>
#include <vector>
#include <cmath>
#include <string>
#include "matrix_checks.h"

// single source of truth for the mask storage type
using mask_t = uint8_t;
using MaskMat = arma::Mat<mask_t>;

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

static inline void validate_knn_inputs(
    const arma::mat &obj,
    const arma::uword k,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const int method,
    const double dist_pow)
{
    if (obj.n_rows == 0)
    {
        Rcpp::stop("obj must have at least one row.");
    }

    if (obj.n_cols < 2)
    {
        Rcpp::stop("obj must have at least two columns.");
    }

    if (k == 0)
    {
        Rcpp::stop("k must be >= 1.");
    }

    if (method != 0 && method != 1)
    {
        Rcpp::stop("Invalid method: 0=Euclidean, 1=Manhattan.");
    }

    if (!std::isfinite(dist_pow) || dist_pow < 0.0)
    {
        Rcpp::stop("dist_pow must be finite and >= 0.");
    }

    auto check_group = [&](const arma::uvec &g, const char *name)
    {
        for (arma::uword i = 0; i < g.n_elem; ++i)
        {
            if (g(i) >= obj.n_cols)
            {
                Rcpp::stop(
                    std::string(name) +
                    " contains out-of-bounds column index " +
                    std::to_string(g(i)) +
                    ". obj has " +
                    std::to_string(obj.n_cols) +
                    " columns.");
            }
        }
    };

    check_group(grp_impute, "grp_impute");
    check_group(grp_miss_no_imp, "grp_miss_no_imp");
    check_group(grp_complete, "grp_complete");

    // The three groups partition the working set. The kernel relies on that
    // when it tags a group-3 neighbor by its offset into grp_complete, and a
    // column listed twice would appear as its own neighbor at distance zero
    // and take the whole weight, so the partition is checked, not assumed.
    std::vector<mask_t> seen(obj.n_cols, 0);
    auto claim_group = [&](const arma::uvec &g, const char *name)
    {
        for (arma::uword i = 0; i < g.n_elem; ++i)
        {
            if (seen[g(i)])
            {
                Rcpp::stop(
                    std::string(name) +
                    " contains column index " +
                    std::to_string(g(i)) +
                    ", which already belongs to another group. The three "
                    "groups must be disjoint.");
            }
            seen[g(i)] = 1;
        }
    };

    claim_group(grp_impute, "grp_impute");
    claim_group(grp_miss_no_imp, "grp_miss_no_imp");
    claim_group(grp_complete, "grp_complete");

    // Group 3 is read straight out of obj with no mask
    // (impute_column_values), because its columns are declared complete. A
    // NaN there is never excluded from the weighted sum, so it produces a
    // missing imputed value with no error. This is where "complete" is
    // enforced; it costs one pass over the group-3 block, against the
    // k * n_complete distance passes the kernel is about to do.
    for (arma::uword i = 0; i < grp_complete.n_elem; ++i)
    {
        const double *col = obj.colptr(grp_complete(i));
        for (arma::uword r = 0; r < obj.n_rows; ++r)
        {
            if (std::isnan(col[r]))
            {
                Rcpp::stop(
                    std::string("grp_complete column ") +
                    std::to_string(grp_complete(i) + 1) +
                    " has a missing value at row " +
                    std::to_string(r + 1) +
                    ". grp_complete columns must be fully observed.");
            }
        }
    }

    const arma::uword n_working =
        grp_impute.n_elem + grp_miss_no_imp.n_elem + grp_complete.n_elem;

    if (grp_impute.n_elem > 0)
    {
        if (n_working < 2)
        {
            Rcpp::stop("Need at least one candidate neighbor column.");
        }

        if (k > n_working - 1)
        {
            Rcpp::stop(
                std::string("k = ") +
                std::to_string(k) +
                " exceeds available candidate columns = " +
                std::to_string(n_working - 1) +
                ".");
        }
    }
}

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
