#include <Rcpp.h>

#include "matrix_checks.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

// Column and row statistics over a numeric matrix, read in place. Each export
// takes the matrix as `obj` and returns the R type its value naturally has:
// counts as an integer vector, statistics as a double vector. Every index and
// offset is std::size_t, so a matrix past 2^31 cells indexes correctly; a
// count is bounded by nrow() or ncol(), which R holds as int.

// Row 1 the column minima, row 2 the maxima. Comparisons skip NaN.
// [[Rcpp::export]]
Rcpp::NumericMatrix col_min_max(const Rcpp::NumericMatrix &obj)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  if (n_rows == 0)
  {
    Rcpp::stop("`obj` has no rows, so its columns have no minimum or maximum.");
  }

  const double *x = obj.begin();
  Rcpp::NumericMatrix out(2, obj.ncol());
  double *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j)
  {
    const double *col = x + j * n_rows;
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();

    for (std::size_t i = 0; i < n_rows; ++i)
    {
      const double v = col[i];
      if (v < lo)
        lo = v;
      if (v > hi)
        hi = v;
    }

    o[2 * j] = lo;
    o[2 * j + 1] = hi;
  }
  return out;
}

// Sample variance of each column, by Welford's update.
// [[Rcpp::export]]
Rcpp::NumericVector col_vars_internal(const Rcpp::NumericMatrix &obj)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  const double *x = obj.begin();
  Rcpp::NumericVector out(Rcpp::no_init(obj.ncol()));
  double *o = out.begin();

  for (std::size_t k = 0; k < n_cols; ++k)
  {
    const double *col = x + k * n_rows;

    double mean = 0.0;
    double mean2 = 0.0;
    std::size_t n = 0;

    for (std::size_t i = 0; i < n_rows; ++i)
    {
      const double v = col[i];

      // missingness is NA/NaN.
      // Inf/-Inf are not treated as missing; they enter the arithmetic.
      // This gives R-like behavior: variance involving Inf can become NaN.
      if (!std::isnan(v))
      {
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
void check_inf(const Rcpp::NumericMatrix &obj)
{
  stop_on_inf(obj);
}

// Inf/-Inf and all-NA/NaN columns. The gate for the imputation entry points.
// [[Rcpp::export]]
void check_finite(const Rcpp::NumericMatrix &obj)
{
  const std::size_t nr = static_cast<std::size_t>(obj.nrow());
  const std::size_t nc = static_cast<std::size_t>(obj.ncol());
  const double *x = obj.begin();

  for (std::size_t col = 0; col < nc; ++col)
  {
    const double *cp = x + col * nr;

    bool has_finite = false;

    // one pass per column, two refusals:
    // Inf/-Inf is rejected outright, wherever it appears. A column that ends
    // the pass with no finite cell holds nothing but NA/NaN and is rejected
    // too. Inf is checked first, so the second message never has to mention it.
    for (std::size_t row = 0; row < nr; ++row)
    {
      const double v = cp[row];

      if (std::isinf(v))
      {
        Rcpp::stop(
            std::string("Infinite value found at row ") +
            std::to_string(row + 1) +
            ", column " +
            std::to_string(col + 1) +
            ". Infinite values are not supported.");
      }

      // Inf is already ruled out, so not-NaN means finite here.
      if (!std::isnan(v))
      {
        has_finite = true;
      }
    }

    if (!has_finite)
    {
      Rcpp::stop("All NA/NaN column detected");
    }
  }
}

// NA/NaN count of each column. Inf/-Inf are not counted as missing.
// [[Rcpp::export]]
Rcpp::IntegerVector col_miss_internal(const Rcpp::NumericMatrix &obj)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  const double *x = obj.begin();
  Rcpp::IntegerVector out(Rcpp::no_init(obj.ncol()));
  int *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j)
  {
    const double *p = x + j * n_rows;

    int count = 0;

    for (std::size_t i = 0; i < n_rows; ++i)
    {
      count += std::isnan(p[i]);
    }

    o[j] = count;
  }

  return out;
}

// NA/NaN count of each row. Inf/-Inf are not counted as missing.
// [[Rcpp::export]]
Rcpp::IntegerVector row_miss_internal(const Rcpp::NumericMatrix &obj)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  const double *x = obj.begin();
  Rcpp::IntegerVector out(obj.nrow()); // zero-filled
  int *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j)
  {
    const double *p = x + j * n_rows;

    for (std::size_t i = 0; i < n_rows; ++i)
    {
      o[i] += std::isnan(p[i]);
    }
  }

  return out;
}
