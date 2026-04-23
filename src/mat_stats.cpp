#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <cmath>

// [[Rcpp::export]]
arma::mat col_min_max(const arma::mat &mat)
{
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

  RcppThread::parallelFor(
      0,
      mat.n_cols,
      [&](arma::uword k)
      {
        const double *col = mat.colptr(k);
        double mean = 0.0;
        double mean2 = 0.0;
        arma::uword n = 0;
        for (arma::uword i = 0; i < mat.n_rows; ++i)
        {
          if (std::isfinite(col[i]))
          {
            ++n;
            double delta = col[i] - mean;
            mean += delta / static_cast<double>(n);
            double delta2 = col[i] - mean;
            mean2 += delta * delta2;
          }
        }
        out(k) = (n < 2) ? arma::datum::nan : mean2 / (n - 1);
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
    Rcpp::stop("col_idx out of bounds");
  }

  std::vector<uint8_t> needs_imp(mat.n_cols, 0);
  for (arma::uword i = 0; i < col_idx.n_elem; ++i)
  {
    needs_imp[col_idx(i)] = 1;
  }
  // initialize out matrix
  const arma::uword n_rows = mat.n_rows;
  arma::mat out(n_rows, mat.n_cols, arma::fill::none);

  cores = std::max(1, cores);
  const size_t n_threads = static_cast<size_t>(cores);
  const size_t col_bytes = n_rows * sizeof(double);

  RcppThread::parallelFor(
      0, mat.n_cols,
      [&](arma::uword c)
      {
        const double *in_ptr = mat.colptr(c);
        double *out_ptr = out.colptr(c);

        if (!needs_imp[c])
        {
          std::memcpy(out_ptr, in_ptr, col_bytes);
          return;
        }

        // pass 1: sum + count
        double sum = 0.0;
        arma::uword n = 0;
        for (arma::uword i = 0; i < n_rows; ++i)
        {
          const double x = in_ptr[i];
          const bool fin = std::isfinite(x);
          sum += fin ? x : 0.0;
          n += fin ? 1u : 0u;
        }

        if (n == 0 || n == n_rows)
        {
          std::memcpy(out_ptr, in_ptr, col_bytes);
          return;
        }

        // pass 2: write.
        const double mean_val = sum / static_cast<double>(n);
        for (arma::uword i = 0; i < n_rows; ++i)
        {
          const double x = in_ptr[i];
          out_ptr[i] = std::isfinite(x) ? x : mean_val;
        }
      },
      n_threads);

  return out;
}

// [[Rcpp::export]]
void check_finite(const arma::mat &mat)
{
  const int nr = mat.n_rows;
  const int nc = mat.n_cols;

  for (int col = 0; col < nc; ++col)
  {
    const double *cp = mat.colptr(col);
    bool has_finite = false;
    for (int row = 0; row < nr; ++row)
    {
      if (std::isfinite(cp[row]))
      {
        has_finite = true;
        break;
      }
    }
    if (!has_finite)
    {
      Rcpp::stop("All NA/Inf column detected");
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
  // column-major scan keeps memory access sequential
  for (arma::uword j = 0; j < n_cols; ++j)
  {
    const double *p = obj.colptr(j);
    for (arma::uword i = 0; i < n_rows; ++i)
    {
      out(i) += std::isnan(p[i]);
    }
  }
  return out;
}
