// RcppArmadillo.h only for arma_uword_bytes(), which stays until Armadillo
// leaves the package (dev/eigen_conversion.md, Phase 3). It comes first:
// RcppArmadillo.h refuses to follow Rcpp.h, which matrix_checks.h includes.
#include <RcppArmadillo.h>
#include <RcppThread.h>

#include "matrix_checks.h"
#include "par_for.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// The exports take `const Rcpp::NumericMatrix &`. RcppArmadillo built the
// `const arma::mat &` they took before from that same Rcpp::Matrix
// (RcppArmadilloAs.h, ArmaMat_InputParameter), so the same inputs are
// accepted, coerced and refused, and a double matrix is read in place. Dims
// are int in R; every index and offset below is std::size_t, so a matrix
// past 2^31 cells indexes correctly. A worker reads raw pointers taken
// before the parallel region and nothing of R's.

// [[Rcpp::export]]
Rcpp::NumericMatrix col_min_max(const Rcpp::NumericMatrix &obj)
{
  // internal helper only, guaranteed no NA/Inf/NaN
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  // Armadillo's min() refused an empty column with this exception; a call
  // with no rows fails as it did.
  if (n_rows == 0 && n_cols > 0)
  {
    throw std::logic_error("min(): object has no elements");
  }

  const double *x = obj.begin();
  Rcpp::NumericMatrix out(2, obj.ncol()); // row 0 = min, row 1 = max
  double *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j)
  {
    const double *col = x + j * n_rows;

    // Armadillo's start values and comparisons, so NaN is skipped. It kept
    // one running value over even rows and one over odd rows, so when the
    // extreme is a zero held as both 0 and -0 it can return the other sign.
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

// [[Rcpp::export]]
Rcpp::NumericMatrix col_vars_internal(const Rcpp::NumericMatrix &obj,
                                      int cores = 1)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  // 1 x p, the shape Armadillo's rowvec came back to R in.
  Rcpp::NumericMatrix out(1, obj.ncol());

  const double *x = obj.begin();
  double *o = out.begin();
  const double na = NA_REAL;

  // parallelFor setup
  cores = std::max(1, cores);
  const size_t n_threads = static_cast<size_t>(cores);

  par_for(
      0,
      n_cols,
      [&](size_t k)
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

        // return NA when fewer than two non-missing observations are available.
        o[k] = (n < 2) ? na : mean2 / static_cast<double>(n - 1);
      },
      n_threads);

  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix mean_imp_col_internal(const Rcpp::NumericMatrix &obj,
                                          const Rcpp::NumericVector &col_idx,
                                          int cores = 1)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  // `col_idx` was an arma::uvec, which RcppArmadillo filled by reading the
  // argument as doubles and casting each to an unsigned index
  // (RcppArmadilloAs.h, ArmaVec_InputParameter). A value is a column when it
  // truncates to 0 .. n_cols - 1, so 2.7 is column 2 and -0.5 is column 0,
  // as the cast made them; anything else is out of bounds. The cast was
  // undefined for NA, Inf and values outside the unsigned range, and on x86
  // it took 1e20 and Inf as column 0; those are refused here too.
  std::vector<uint8_t> needs_imp(n_cols, 0);
  const R_xlen_t n_idx = col_idx.size();

  for (R_xlen_t i = 0; i < n_idx; ++i)
  {
    const double v = col_idx[i];
    if (!(v > -1.0 && v < static_cast<double>(n_cols)))
    {
      Rcpp::stop("`col_idx` out of bounds");
    }
    needs_imp[static_cast<std::size_t>(v)] = 1;
  }

  const double *x = obj.begin();

  // validate selected columns before parallel work.
  // NA/NaN are allowed because they are the values to impute.
  // Inf/-Inf are treated as invalid observed values.
  for (std::size_t c = 0; c < n_cols; ++c)
  {
    if (!needs_imp[c])
      continue;

    const double *ptr = x + c * n_rows;

    for (std::size_t i = 0; i < n_rows; ++i)
    {
      if (std::isinf(ptr[i]))
      {
        Rcpp::stop(
            "Infinite value found in column %d, row %d. "
            "Mean imputation requires finite observed values.",
            c + 1,
            i + 1);
      }
    }
  }

  Rcpp::NumericMatrix out(Rcpp::no_init_matrix(obj.nrow(), obj.ncol()));
  double *o = out.begin();
  const uint8_t *imp = needs_imp.data();

  cores = std::max(1, cores);
  const size_t n_threads = static_cast<size_t>(cores);
  const size_t col_bytes = n_rows * sizeof(double);

  par_for(
      0,
      n_cols,
      [&](size_t c)
      {
        const double *in_ptr = x + c * n_rows;
        double *out_ptr = o + c * n_rows;

        if (!imp[c])
        {
          std::memcpy(out_ptr, in_ptr, col_bytes);
          return;
        }

        double sum = 0.0;
        std::size_t n = 0;

        for (std::size_t i = 0; i < n_rows; ++i)
        {
          const double v = in_ptr[i];

          if (!std::isnan(v))
          {
            sum += v;
            ++n;
          }
        }

        if (n == 0 || n == n_rows)
        {
          std::memcpy(out_ptr, in_ptr, col_bytes);
          return;
        }

        const double mean_val = sum / static_cast<double>(n);

        for (std::size_t i = 0; i < n_rows; ++i)
        {
          const double v = in_ptr[i];
          out_ptr[i] = std::isnan(v) ? mean_val : v;
        }
      },
      n_threads);

  return out;
}

// Width of an Armadillo index, in bytes. Exists so the test suite can assert
// that -DARMA_64BIT_WORD=1 survived in src/Makevars: without it arma indices
// are 32 bits and a matrix past 2^32-1 elements wraps silently rather than
// erroring. Nothing else calls this.
// [[Rcpp::export]]
int arma_uword_bytes()
{
  return static_cast<int>(sizeof(arma::uword));
}

// Inf/-Inf only. For callers that reject infinities but legitimately accept a
// column with no observed value: `mean_imp_col()` leaves such columns
// untouched by design, and `sample_na_loc()` reports them through its own
// variance pre-check.
// [[Rcpp::export]]
void check_inf(const Rcpp::NumericMatrix &obj)
{
  stop_on_inf(obj.begin(),
              static_cast<std::size_t>(obj.nrow()),
              static_cast<std::size_t>(obj.ncol()));
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

// n_cols x 1 doubles, the shape and type Armadillo's uvec came back to R in.
// [[Rcpp::export]]
Rcpp::NumericMatrix col_miss_internal(const Rcpp::NumericMatrix &obj)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  Rcpp::NumericMatrix out(obj.ncol(), 1);
  const double *x = obj.begin();
  double *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j)
  {
    const double *p = x + j * n_rows;

    std::size_t count = 0;

    for (std::size_t i = 0; i < n_rows; ++i)
    {
      // count NA/NaN as missing.
      // Inf/-Inf are not counted as missing.
      count += std::isnan(p[i]);
    }

    o[j] = static_cast<double>(count);
  }

  return out;
}

// n_rows x 1 doubles, as col_miss_internal(). A count is a whole number
// below 2^31, so adding in double is exact.
// [[Rcpp::export]]
Rcpp::NumericMatrix row_miss_internal(const Rcpp::NumericMatrix &obj)
{
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());

  Rcpp::NumericMatrix out(obj.nrow(), 1); // zero-filled
  const double *x = obj.begin();
  double *o = out.begin();

  for (std::size_t j = 0; j < n_cols; ++j)
  {
    const double *p = x + j * n_rows;

    for (std::size_t i = 0; i < n_rows; ++i)
    {
      // count NA/NaN as missing.
      // Inf/-Inf are not counted as missing.
      o[i] += std::isnan(p[i]);
    }
  }

  return out;
}
