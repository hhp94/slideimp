#include "imputed_value.h"
#include <cmath>
#include <algorithm>

// =============================================================================
// initialize_result_matrix
// -----------------------------------------------------------------------------
// Only iterates group 1 (local positions [0, n_imp) in nmiss_masked).
// Populates result columns 0 and 1 (1-based row / original 1-based column).
// Column 2 is filled later by impute_column_values.
// =============================================================================
arma::mat initialize_result_matrix(
    const arma::mat &nmiss_masked,
    const arma::uvec &grp_impute,
    const GroupLayout &layout,
    arma::uvec &col_offsets,
    std::vector<arma::uvec> &rows_to_impute_vec)
{
    const arma::uword n_col_result = 3; // column 1, 2 are row, col. 3 is value
    const arma::uword n_imp = layout.n_imp;

    if (n_imp == 0)
    {
        return arma::mat(0, n_col_result);
    }

    // per-column missing counts, read directly from the masked matrix
    arma::uvec n_col_miss(n_imp);
    for (arma::uword i = 0; i < n_imp; ++i)
    {
        n_col_miss(i) = nmiss_masked.n_rows - arma::accu(nmiss_masked.col(i));
    }

    arma::uword sum_missing = arma::accu(n_col_miss);
    if (sum_missing == 0)
    {
        return arma::mat(0, n_col_result);
    }

    arma::mat result(sum_missing, n_col_result);
    result.fill(arma::datum::nan);

    col_offsets.set_size(n_imp + 1);
    col_offsets.fill(arma::fill::zeros);
    col_offsets.subvec(1, n_imp) = arma::cumsum(n_col_miss);

    rows_to_impute_vec.resize(n_imp);

    for (arma::uword i = 0; i < n_imp; ++i)
    {
        // Group 1 lives at local positions [0, n_imp) in nmiss_masked
        rows_to_impute_vec[i] = arma::find(nmiss_masked.col(i) == 0);
        const arma::uvec &rows = rows_to_impute_vec[i];
        const arma::uword orig_col_1based = grp_impute(i) + 1;

        for (arma::uword r = 0; r < rows.n_elem; ++r)
        {
            const arma::uword res_row = col_offsets(i) + r;
            result(res_row, 0) = rows(r) + 1;     // 1-based row
            result(res_row, 1) = orig_col_1based; // original 1-based column
        }
    }

    return result;
}

// =============================================================================
// impute_column_values
// -----------------------------------------------------------------------------
// Neighbors are split into masked (groups 1+2) and unmasked (group 3) buckets
// before accumulation. Each bucket gets its own clean inner loop.
//
// Unmasked simplification: every row is valid, so wt_ptr[r] += w for every r,
// every unmasked neighbor collapses to a single "total_unmasked_w" added to
// every row once at the end. The ws_ptr accumulation still varies per row
// (different obj values).
// =============================================================================
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
    const arma::uvec &grp_complete)
{
    const arma::uword n_rows = rows_to_impute.n_elem;
    if (n_rows == 0)
    {
        return;
    }

    arma::vec weighted_sums(n_rows, arma::fill::zeros);
    arma::vec weight_totals(n_rows, arma::fill::zeros);

    double *__restrict__ ws_ptr = weighted_sums.memptr();
    double *__restrict__ wt_ptr = weight_totals.memptr();

    const arma::uword complete_start = layout.complete_start();

    std::vector<arma::uword> masked_j;
    std::vector<arma::uword> complete_j;
    masked_j.reserve(nn_columns.n_elem);
    complete_j.reserve(nn_columns.n_elem);

    for (arma::uword j = 0; j < nn_columns.n_elem; ++j)
    {
        if (nn_columns(j) < complete_start)
        {
            masked_j.push_back(j);
        }
        else
        {
            complete_j.push_back(j);
        }
    }

    // masked neighbors (groups 1+2) ----
    for (arma::uword j : masked_j)
    {
        const double w = nn_weights(j);
        const arma::uword local = nn_columns(j);
        const double *obj_col = obj_masked.colptr(local);
        const double *nmiss_col = nmiss_masked.colptr(local);

        for (arma::uword r = 0; r < n_rows; ++r)
        {
            const arma::uword row_idx = rows_to_impute(r);
            double valid = nmiss_col[row_idx];
            double weight = w * valid;
            ws_ptr[r] += weight * obj_col[row_idx];
            wt_ptr[r] += weight;
        }
    }

    // complete neighbors (group 3) ----
    // valid == 1 for every row, so weight == w uniformly across rows.
    double total_complete_w = 0.0;

    for (arma::uword j : complete_j)
    {
        const double w = nn_weights(j);
        const arma::uword c = nn_columns(j) - complete_start;
        const double *obj_col = obj.colptr(grp_complete(c));

        for (arma::uword r = 0; r < n_rows; ++r)
        {
            ws_ptr[r] += w * obj_col[rows_to_impute(r)];
        }
        total_complete_w += w;
    }

    if (total_complete_w > 0.0)
    {
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            wt_ptr[r] += total_complete_w;
        }
    }

    for (arma::uword r = 0; r < n_rows; ++r)
    {
        result(col_offset + r, 2) =
            (wt_ptr[r] > 0.0) ? (ws_ptr[r] / wt_ptr[r]) : arma::datum::nan;
    }
}
