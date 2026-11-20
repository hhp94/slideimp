#include "imputed_value.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// initialize_result
// -----------------------------------------------------------------------------
// Only iterates group 1 (local positions [0, n_imp) in nmiss_masked).
// Fills `loc` (1-based row / original 1-based column) and sets `value` to NA.
// The values are written later by impute_column_values.
//
// n_col_valid is computed once during the copy-with-mask pass so we never
// sum the byte mask.
// -----------------------------------------------------------------------------
void initialize_result(const std::vector<MaskByte> &nmiss_masked,
                       const int n_rows, const Rcpp::IntegerVector &grp_impute,
                       const GroupLayout &layout,
                       const std::vector<int> &n_col_valid,
                       std::vector<std::size_t> &col_offsets,
                       Rcpp::IntegerMatrix &loc, Rcpp::NumericVector &value) {
  const int n_imp = layout.n_imp;

  col_offsets.clear();
  loc = Rcpp::IntegerMatrix(0, 2);
  value = Rcpp::NumericVector(0);

  if (n_imp == 0) {
    return;
  }

  col_offsets.assign(static_cast<std::size_t>(n_imp) + 1, 0);

  // at most n_imp * n_rows, so a size_t sum cannot overflow
  std::size_t sum_missing = 0;

  for (int i = 0; i < n_imp; ++i) {
    const int n_miss_i = n_rows - n_col_valid[i];
    sum_missing += static_cast<std::size_t>(n_miss_i);
    col_offsets[i + 1] = sum_missing;
  }

  if (sum_missing == 0) {
    return;
  }

  // `loc` is an R matrix, whose dims are int.
  if (sum_missing > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("the columns to impute hold %d missing values, but the result "
               "has a row for each and an R matrix holds at most %d rows. "
               "Impute fewer columns at a time with 'subset'",
               sum_missing, std::numeric_limits<int>::max());
  }

  const int n_total = static_cast<int>(sum_missing);
  loc = Rcpp::IntegerMatrix(Rcpp::no_init_matrix(n_total, 2));
  value = Rcpp::NumericVector(n_total, NA_REAL);
  int *loc_row = loc.begin();
  int *loc_col = loc.begin() + n_total;

  for (int i = 0; i < n_imp; ++i) {
    const MaskByte *mask_col = col_ptr(nmiss_masked.data(), i, n_rows);
    const int orig_col_1based = grp_impute[i] + 1;
    std::size_t w = col_offsets[i];

    // copy_with_mask counted the ones of this mask column into n_col_valid,
    // so it holds exactly as many zeros as its slice of the result has rows.
    for (int r = 0; r < n_rows; ++r) {
      if (mask_col[r] == 0) {
        loc_row[w] = r + 1;           // 1-based row
        loc_col[w] = orig_col_1based; // original 1-based column
        ++w;
      }
    }
  }
}

// -----------------------------------------------------------------------------
// impute_column_values
// -----------------------------------------------------------------------------
// Neighbors are split into masked (groups 1+2) and unmasked (group 3) buckets
// before accumulation. Each bucket gets its own clean inner loop.
//
// Unmasked simplification: every row is valid, so wt_ptr[r] += w for every r,
// every unmasked neighbor collapses to a single "total_unmasked_w" added to
// every row once at the end. The ws_ptr accumulation still varies per row
// (different obj values).
// -----------------------------------------------------------------------------
void impute_column_values(double *out, const double *obj_masked,
                          const MaskByte *nmiss_masked, const int n_rows,
                          const GroupLayout &layout,
                          const std::vector<int> &nn_columns,
                          const std::vector<double> &nn_weights,
                          const int *loc_rows, const int n_miss,
                          const double *obj, const int *grp_complete) {
  if (n_miss == 0) {
    return;
  }

  std::vector<double> weighted_sums(n_miss, 0.0);
  std::vector<double> weight_totals(n_miss, 0.0);

  double *ws_ptr = weighted_sums.data();
  double *wt_ptr = weight_totals.data();

  const int complete_start = layout.complete_start();

  std::vector<std::size_t> masked_j;
  std::vector<std::size_t> complete_j;
  masked_j.reserve(nn_columns.size());
  complete_j.reserve(nn_columns.size());

  for (std::size_t j = 0; j < nn_columns.size(); ++j) {
    if (nn_columns[j] < complete_start) {
      masked_j.push_back(j);
    } else {
      complete_j.push_back(j);
    }
  }

  // masked neighbors (groups 1+2) ----
  for (std::size_t j : masked_j) {
    const double w = nn_weights[j];
    const int local = nn_columns[j];
    const double *obj_col = col_ptr(obj_masked, local, n_rows);
    const MaskByte *nmiss_col = col_ptr(nmiss_masked, local, n_rows);

    for (int r = 0; r < n_miss; ++r) {
      const int row_idx = loc_rows[r] - 1; // loc holds 1-based rows
      double valid = nmiss_col[row_idx];   // uint8 -> double: exact 0.0/1.0
      double weight = w * valid;
      ws_ptr[r] += weight * obj_col[row_idx];
      wt_ptr[r] += weight;
    }
  }

  // complete neighbors (group 3) ----
  // valid == 1 for every row, so weight == w uniformly across rows.
  double total_complete_w = 0.0;

  for (std::size_t j : complete_j) {
    const double w = nn_weights[j];
    const int c = nn_columns[j] - complete_start;
    const double *obj_col = col_ptr(obj, grp_complete[c], n_rows);

    for (int r = 0; r < n_miss; ++r) {
      ws_ptr[r] += w * obj_col[loc_rows[r] - 1];
    }
    total_complete_w += w;
  }

  if (total_complete_w > 0.0) {
    for (int r = 0; r < n_miss; ++r) {
      wt_ptr[r] += total_complete_w;
    }
  }

  // a row where none of the chosen neighbors is observed keeps the NA
  // initialize_result() wrote.
  for (int r = 0; r < n_miss; ++r) {
    if (wt_ptr[r] > 0.0) {
      out[r] = ws_ptr[r] / wt_ptr[r];
    }
  }
}
