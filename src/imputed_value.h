#ifndef SLIDEIMP_IMPUTED_VALUE_H
#define SLIDEIMP_IMPUTED_VALUE_H

#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "matrix_checks.h"

// single source of truth for the mask storage type
using MaskByte = std::uint8_t;

// GroupLayout: single source of truth for group boundaries.
// Groups 1+2 (imputed + missing-no-impute) share identical kernel treatment
// and live together in obj_masked / nmiss_masked at local cols [0, n_masked()).
// n_imp is retained so initialize_result knows which columns actually get
// imputed.
//
// complete_start() is a "virtual" index: there is no physical column at that
// position in obj_masked. It's used as a tag in nn_columns to mark
// "this neighbor came from group 3, offset (local - complete_start()) into
// grp_complete".
//
// The groups are disjoint and in range (validate_knn_inputs), so every count
// here is at most ncol().
struct GroupLayout {
  int n_imp;
  int n_mni;
  int n_complete;
  int mni_start() const { return n_imp; }
  int n_masked() const { return n_imp + n_mni; }
  int complete_start() const { return n_imp + n_mni; }
  int n_working() const { return n_imp + n_mni + n_complete; }
};

static inline void
validate_knn_inputs(const Rcpp::NumericMatrix &obj, const int k,
                    const Rcpp::IntegerVector &grp_impute,
                    const Rcpp::IntegerVector &grp_miss_no_imp,
                    const Rcpp::IntegerVector &grp_complete, const int method,
                    const double dist_pow) {
  const int n_rows = obj.nrow();
  const int n_cols = obj.ncol();

  if (n_rows == 0) {
    Rcpp::stop("'obj' must have at least one row");
  }

  if (n_cols < 2) {
    Rcpp::stop("'obj' must have at least two columns");
  }

  // NA_integer_ is INT_MIN, so an NA k fails here too.
  if (k < 1) {
    Rcpp::stop("'k' must be >= 1");
  }

  if (method != 0 && method != 1) {
    Rcpp::stop("'method' must be 0 (Euclidean) or 1 (Manhattan)");
  }

  if (!std::isfinite(dist_pow) || dist_pow < 0.0) {
    Rcpp::stop("'dist_pow' must be finite and >= 0");
  }

  stop_on_bad_col_index(grp_impute, n_cols, "grp_impute");
  stop_on_bad_col_index(grp_miss_no_imp, n_cols, "grp_miss_no_imp");
  stop_on_bad_col_index(grp_complete, n_cols, "grp_complete");

  // The three groups partition the working set. The kernel relies on that
  // when it tags a group-3 neighbor by its offset into grp_complete, and a
  // column listed twice would appear as its own neighbor at distance zero
  // and take the whole weight, so the partition is checked, not assumed.
  std::vector<MaskByte> seen(n_cols, 0);
  auto claim_group = [&](const Rcpp::IntegerVector &g, const char *name) {
    for (R_xlen_t i = 0; i < g.size(); ++i) {
      if (seen[g[i]]) {
        Rcpp::stop("'%s' contains column index %d, which already belongs to "
                   "another group. The three groups must be disjoint",
                   name, g[i]);
      }
      seen[g[i]] = 1;
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
  const double *x = obj.begin();
  for (R_xlen_t i = 0; i < grp_complete.size(); ++i) {
    const double *col = col_ptr(x, grp_complete[i], n_rows);
    for (int r = 0; r < n_rows; ++r) {
      if (std::isnan(col[r])) {
        Rcpp::stop("'grp_complete' column %d has a missing value at row %d. "
                   "'grp_complete' columns must be fully observed",
                   grp_complete[i] + 1, r + 1);
      }
    }
  }

  // disjoint and in range, so at most n_cols
  const int n_working = static_cast<int>(
      grp_impute.size() + grp_miss_no_imp.size() + grp_complete.size());

  if (grp_impute.size() > 0) {
    if (n_working < 2) {
      Rcpp::stop("need at least one candidate neighbor column");
    }

    if (k > n_working - 1) {
      Rcpp::stop("'k' is %d, more than the %d candidate neighbor columns", k,
                 n_working - 1);
    }
  }
}

// The result has one entry per missing cell of group 1 (local columns
// [0, n_imp) of the masked copy), column by column and down each column:
// `loc` holds its 1-based row and original 1-based column, `value` holds NA
// until impute_column_values() writes it. col_offsets[i] is where column i's
// entries start, rows_to_impute[i] the 0-based rows they stand for.
// n_col_valid holds the observed count of each masked column, taken during
// the copy-with-mask pass, so the missing counts need no second pass over
// the mask. Only entries [0, layout.n_imp) are read here.
void initialize_result(const std::vector<MaskByte> &nmiss_masked,
                       const int n_rows, const Rcpp::IntegerVector &grp_impute,
                       const GroupLayout &layout,
                       const std::vector<int> &n_col_valid,
                       std::vector<std::size_t> &col_offsets,
                       std::vector<std::vector<int>> &rows_to_impute,
                       Rcpp::IntegerMatrix &loc, Rcpp::NumericVector &value);

// Runs in a worker: plain pointers and std containers only, never an Rcpp
// object. `out` is this column's slice of the result values.
//
// nn_columns holds tagged positions: values < layout.complete_start() index
// directly into obj_masked/nmiss_masked; values >= layout.complete_start()
// must be de-tagged via (nn_columns[j] - complete_start()) and looked up
// through grp_complete into obj.
void impute_column_values(double *out, const double *obj_masked,
                          const MaskByte *nmiss_masked, const int n_rows,
                          const GroupLayout &layout,
                          const std::vector<int> &nn_columns,
                          const std::vector<double> &nn_weights,
                          const std::vector<int> &rows_to_impute,
                          const double *obj, const int *grp_complete);

#endif
