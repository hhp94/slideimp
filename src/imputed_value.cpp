#include "imputed_value.h"
#include <random>
#include <cmath>
#include <algorithm>

// initialize matrix to store the result. Also modify col_offsets and
// precompute the row indices that need imputation for each column.
arma::mat initialize_result_matrix(
    const arma::mat &nmiss,
    const arma::uvec &col_index_miss,
    arma::uvec &col_offsets,
    std::vector<arma::uvec> &rows_to_impute_vec)
{
    // col 1 = row index, col 2 = col index, col 3 = missing values
    const arma::uword n_col_result = 3;

    if (col_index_miss.n_elem == 0)
    {
        return arma::mat(0, n_col_result);
    }

    // compute per-column missing counts internally
    arma::uvec n_col_miss(col_index_miss.n_elem);
    for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
    {
        n_col_miss(i) = nmiss.n_rows - arma::accu(nmiss.col(col_index_miss(i)));
    }

    arma::uword sum_missing = arma::accu(n_col_miss);
    if (sum_missing == 0)
    {
        return arma::mat(0, n_col_result);
    }

    arma::mat result(sum_missing, n_col_result);
    result.fill(arma::datum::nan);

    col_offsets.set_size(n_col_miss.n_elem + 1);
    col_offsets.fill(arma::fill::zeros);
    col_offsets.subvec(1, n_col_miss.n_elem) = arma::cumsum(n_col_miss);

    // Precompute row indices and fill result with ij indices (1-based for R)
    rows_to_impute_vec.resize(col_index_miss.n_elem);
    for (arma::uword i = 0; i < col_index_miss.n_elem; ++i)
    {
        const arma::uword target_col_idx = col_index_miss(i);
        rows_to_impute_vec[i] = arma::find(nmiss.col(target_col_idx) == 0);
        const arma::uvec &rows_to_impute = rows_to_impute_vec[i];

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

// ======================== LOOP ORDER NOTE ====================================
// The intuitive approach is: for each missing row, loop over the k neighbors
// and accumulate a weighted average. That's easy to understand, but it means the
// inner loop touches k *different* column pointers on every iteration so the CPU
// can't prefetch any of them, and the compiler can't vectorize because consecutive
// iterations read from unrelated memory. This is a slight optimization
//
// We iterate NEIGHBORS instead in the outer loop and ROWS in the inner loop.
// This gives us:
//  1. The neighbor's weight `w` is loaded once and reused for every miss row.
//  2. The inner loop reads through exactly one base pointer (the neighbor's column)
//  plus varying row offsets.
// =============================================================================
void impute_column_values(
    arma::mat &result,
    const arma::mat &obj,
    const arma::mat &nmiss,
    const arma::uword col_offset,
    const arma::uword target_col_idx,
    const arma::uvec &nn_columns_vec,
    const arma::vec &nn_weights,
    const arma::uvec &rows_to_impute)
{
    const arma::uword n_rows = rows_to_impute.n_elem;
    if (n_rows == 0)
    {
        return;
    }
    // accumulation buffers: one entry per missing row in this column. This enables
    // vectorization. Previously, we had two double to hold the accumulators
    arma::vec weighted_sums(n_rows, arma::fill::zeros);
    arma::vec weight_totals(n_rows, arma::fill::zeros);
    double *__restrict__ ws_ptr = weighted_sums.memptr();
    double *__restrict__ wt_ptr = weight_totals.memptr();

    // pre-cache the column pointers for each neighbor because this doesn't change
    // in the inner loop.
    std::vector<const double *> nn_nmiss_ptrs(nn_columns_vec.n_elem);
    std::vector<const double *> nn_obj_ptrs(nn_columns_vec.n_elem);
    for (arma::uword j = 0; j < nn_columns_vec.n_elem; ++j)
    {
        nn_nmiss_ptrs[j] = nmiss.colptr(nn_columns_vec(j));
        nn_obj_ptrs[j] = obj.colptr(nn_columns_vec(j));
    }

    // ---- outer loop over neighbors ----
    for (arma::uword j = 0; j < nn_columns_vec.n_elem; ++j)
    {
        const double w = nn_weights(j);
        const double *obj_col = nn_obj_ptrs[j];
        const double *nmiss_col = nn_nmiss_ptrs[j];
        // ---- inner loop over missing rows ----
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            const arma::uword row_idx = rows_to_impute(r);
            double valid = nmiss_col[row_idx];
            double weight = w * valid;
            ws_ptr[r] += weight * obj_col[row_idx];
            wt_ptr[r] += weight;
        }
    }

    // ---- write imputed values to the result matrix ----
    for (arma::uword r = 0; r < n_rows; ++r)
    {
        result(col_offset + r, 2) = (wt_ptr[r] > 0.0) ? (ws_ptr[r] / wt_ptr[r]) : arma::datum::nan;
    }
}
