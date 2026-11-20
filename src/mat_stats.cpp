#include <Rcpp.h>

#include "matrix_checks.h"

#include <cmath>
#include <cstddef>
#include <string>

// Column and row statistics over a numeric matrix, read in place. Each export
// takes the matrix as `obj` and returns the R type its value naturally has:
// counts as an integer vector, statistics as a double vector. Indices and
// offsets follow src/matrix_checks.h.

// Sample variance of each column, by Welford's update.
// [[Rcpp::export]]
Rcpp::NumericVector col_vars_internal(const Rcpp::NumericMatrix &obj) {
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  const double *x = obj.begin();
  Rcpp::NumericVector out(Rcpp::no_init(obj.ncol()));
  double *o = out.begin();

  for (std::size_t k = 0; k < n_cols; ++k) {
    const double *col = col_ptr(x, k, n_rows);

    double mean = 0.0;
    double mean2 = 0.0;
    std::size_t n = 0;

    for (std::size_t i = 0; i < n_rows; ++i) {
      const double v = col[i];

      // missingness is NA/NaN.
      // Inf/-Inf are not treated as missing; they enter the arithmetic.
      // This gives R-like behavior: variance involving Inf can become NaN.
      if (!std::isnan(v)) {
        ++n;

        const double delta = v - mean;
        mean += delta / static_cast<double>(n);

        const double delta2 = v - mean;
        mean2 += delta * delta2;
      }
    }

    // NA when fewer than two non-missing observations are available.
    o[k] = (n < 2) ? NA_REAL : mean2 / static_cast<double>(n - 1);
  }

  return out;
}

// Inf/-Inf only. For callers that reject infinities but legitimately accept a
// column with no observed value: `mean_imp_col()` leaves such columns
// untouched by design, and `sample_na_loc()` reports them through its own
// variance pre-check.
// [[Rcpp::export]]
void check_inf(const Rcpp::NumericMatrix &obj) { stop_on_inf(obj); }

// Inf/-Inf and all-NA/NaN columns. The gate for the imputation entry points.
// [[Rcpp::export]]
void check_finite(const Rcpp::NumericMatrix &obj) {
  const std::size_t nr = static_cast<std::size_t>(obj.nrow());
  const std::size_t nc = static_cast<std::size_t>(obj.ncol());
  const double *x = obj.begin();

  for (std::size_t col = 0; col < nc; ++col) {
    const double *cp = col_ptr(x, col, nr);

    bool has_finite = false;

    // one pass per column, two refusals:
    // Inf/-Inf is rejected outright, wherever it appears. A column that ends
    // the pass with no finite cell holds nothing but NA/NaN and is rejected
    // too. Inf is checked first, so the second message never has to mention it.
    for (std::size_t row = 0; row < nr; ++row) {
      const double v = cp[row];

      if (std::isinf(v)) {
        Rcpp::stop("'obj' holds an infinite value at row %d, column %d. "
                   "Infinite values are not supported",
                   row + 1, col + 1);
      }

      // Inf is already ruled out, so not-NaN means finite here.
      if (!std::isnan(v)) {
        has_finite = true;
      }
    }

    if (!has_finite) {
      Rcpp::stop("'obj' has a column with no observed value, only NA or NaN. "
                 "Remove that column before imputation");
    }
  }
}

// NA/NaN count of each column. Inf/-Inf are not counted as missing.
// [[Rcpp::export]]
Rcpp::IntegerVector col_miss_internal(const Rcpp::NumericMatrix &obj) {
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  const double *x = obj.begin();
  Rcpp::IntegerVector out(Rcpp::no_init(obj.ncol()));
  int *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j) {
    const double *p = col_ptr(x, j, n_rows);

    int count = 0;

    for (std::size_t i = 0; i < n_rows; ++i) {
      count += std::isnan(p[i]);
    }

    o[j] = count;
  }

  return out;
}

// NA/NaN count of each row. Inf/-Inf are not counted as missing.
// [[Rcpp::export]]
Rcpp::IntegerVector row_miss_internal(const Rcpp::NumericMatrix &obj) {
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  const double *x = obj.begin();
  Rcpp::IntegerVector out(obj.nrow()); // zero-filled
  int *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j) {
    const double *p = col_ptr(x, j, n_rows);

    for (std::size_t i = 0; i < n_rows; ++i) {
      o[i] += std::isnan(p[i]);
    }
  }

  return out;
}
