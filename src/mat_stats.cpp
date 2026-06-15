#include <RcppArmadillo.h>
#include <RcppThread.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>

// [[Rcpp::export]]
arma::mat col_min_max(const arma::mat &mat)
{
  // internal helper only, guaranteed no NA/Inf/NaN
  const arma::uword p = mat.n_cols;
  arma::mat out(2, p); // row 0 = min, row 1 = max
  for (arma::uword j = 0; j < p; ++j)
  {
    out(0, j) = mat.col(j).min();
    out(1, j) = mat.col(j).max();
  }
  return out;
}

// [[Rcpp::export]]
arma::rowvec col_vars_internal(const arma::mat &mat, int cores = 1)
{
  arma::rowvec out(mat.n_cols);

  // parallelFor setup
  cores = std::max(1, cores);
  const size_t n_threads = static_cast<size_t>(cores);
  const size_t n_cols = static_cast<size_t>(mat.n_cols);

  RcppThread::parallelFor(
      0,
      n_cols,
      [&](size_t kk)
      {
        const arma::uword k = static_cast<arma::uword>(kk);
        const double *col = mat.colptr(k);

        double mean = 0.0;
        double mean2 = 0.0;
        arma::uword n = 0;

        for (arma::uword i = 0; i < mat.n_rows; ++i)
        {
          const double x = col[i];

          // missingness is NA/NaN.
          // Inf/-Inf are not treated as missing; they enter the arithmetic.
          // This gives R-like behavior: variance involving Inf can become NaN.
          if (!std::isnan(x))
          {
            ++n;

            const double delta = x - mean;
            mean += delta / static_cast<double>(n);

            const double delta2 = x - mean;
            mean2 += delta * delta2;
          }
        }

        // return NaN when fewer than two non-missing observations are available.
        out(k) = (n < 2) ? NA_REAL : mean2 / static_cast<double>(n - 1);
      },
      n_threads);

  return out;
}

// [[Rcpp::export]]
arma::mat mean_imp_col_internal(const arma::mat &mat,
                                const arma::uvec &col_idx,
                                int cores = 1)
{
  if (col_idx.n_elem > 0 && col_idx.max() >= mat.n_cols)
  {
    Rcpp::stop("`col_idx` out of bounds");
  }

  std::vector<uint8_t> needs_imp(mat.n_cols, 0);

  for (arma::uword i = 0; i < col_idx.n_elem; ++i)
  {
    needs_imp[col_idx(i)] = 1;
  }

  const arma::uword n_rows = mat.n_rows;

  // validate selected columns before parallel work.
  // NA/NaN are allowed because they are the values to impute.
  // Inf/-Inf are treated as invalid observed values.
  for (arma::uword c = 0; c < mat.n_cols; ++c)
  {
    if (!needs_imp[c])
      continue;

    const double *ptr = mat.colptr(c);

    for (arma::uword i = 0; i < n_rows; ++i)
    {
      if (std::isinf(ptr[i]))
      {
        Rcpp::stop(
            "Infinite value found in column %d, row %d. "
            "Mean imputation requires finite observed values.",
            static_cast<int>(c + 1),
            static_cast<int>(i + 1));
      }
    }
  }

  arma::mat out(n_rows, mat.n_cols, arma::fill::none);

  cores = std::max(1, cores);
  const size_t n_threads = static_cast<size_t>(cores);
  const size_t n_cols = static_cast<size_t>(mat.n_cols);
  const size_t col_bytes = static_cast<size_t>(n_rows) * sizeof(double);

  RcppThread::parallelFor(
      0,
      n_cols,
      [&](size_t cc)
      {
        const arma::uword c = static_cast<arma::uword>(cc);

        const double *in_ptr = mat.colptr(c);
        double *out_ptr = out.colptr(c);

        if (!needs_imp[c])
        {
          std::memcpy(out_ptr, in_ptr, col_bytes);
          return;
        }

        double sum = 0.0;
        arma::uword n = 0;

        for (arma::uword i = 0; i < n_rows; ++i)
        {
          const double x = in_ptr[i];

          if (!std::isnan(x))
          {
            sum += x;
            ++n;
          }
        }

        if (n == 0 || n == n_rows)
        {
          std::memcpy(out_ptr, in_ptr, col_bytes);
          return;
        }

        const double mean_val = sum / static_cast<double>(n);

        for (arma::uword i = 0; i < n_rows; ++i)
        {
          const double x = in_ptr[i];
          out_ptr[i] = std::isnan(x) ? mean_val : x;
        }
      },
      n_threads);

  return out;
}

// [[Rcpp::export]]
void check_finite(const arma::mat &mat)
{
  const arma::uword nr = mat.n_rows;
  const arma::uword nc = mat.n_cols;

  for (arma::uword col = 0; col < nc; ++col)
  {
    const double *cp = mat.colptr(col);

    bool has_finite = false;

    for (arma::uword row = 0; row < nr; ++row)
    {
      if (std::isfinite(cp[row]))
      {
        has_finite = true;
        break;
      }
    }

    if (!has_finite)
    {
      Rcpp::stop("All NA/NaN/Inf column detected");
    }
  }
}

// [[Rcpp::export]]
arma::uvec col_miss_internal(const arma::mat &obj)
{
  const arma::uword n_rows = obj.n_rows;
  const arma::uword n_cols = obj.n_cols;

  arma::uvec out(n_cols);

  for (arma::uword j = 0; j < n_cols; ++j)
  {
    const double *p = obj.colptr(j);

    arma::uword count = 0;

    for (arma::uword i = 0; i < n_rows; ++i)
    {
      // count NA/NaN as missing.
      // Inf/-Inf are not counted as missing.
      count += std::isnan(p[i]);
    }

    out(j) = count;
  }

  return out;
}

// [[Rcpp::export]]
arma::uvec row_miss_internal(const arma::mat &obj)
{
  const arma::uword n_rows = obj.n_rows;
  const arma::uword n_cols = obj.n_cols;

  arma::uvec out(n_rows, arma::fill::zeros);
  for (arma::uword j = 0; j < n_cols; ++j)
  {
    const double *p = obj.colptr(j);

    for (arma::uword i = 0; i < n_rows; ++i)
    {
      // count NA/NaN as missing.
      // Inf/-Inf are not counted as missing.
      out(i) += std::isnan(p[i]);
    }
  }

  return out;
}
