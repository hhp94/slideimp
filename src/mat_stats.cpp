#include <RcppArmadillo.h>
#include <cmath>
#if defined(_OPENMP)
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat col_min_max(const arma::mat &mat) {
  const arma::uword p = mat.n_cols;
  arma::mat out(2, p); // row 0 = min, row 1 = max
  for (arma::uword j = 0; j < p; ++j) {
    out(0, j) = mat.col(j).min();
    out(1, j) = mat.col(j).max();
  }
  return out;
}

// [[Rcpp::export]]
arma::rowvec col_vars_internal(const arma::mat &mat, arma::uword cores = 1)
{
  arma::rowvec out(mat.n_cols);
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (arma::uword k = 0; k < mat.n_cols; ++k)
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
  }
  return out;
}

// [[Rcpp::export]]
arma::mat mean_imp_col_internal(const arma::mat &mat,
                                const arma::uvec &col_idx,
                                arma::uword cores = 1)
{
  if (col_idx.max() >= mat.n_cols)
  {
    Rcpp::stop("col_idx out of bounds");
  }

  arma::mat out = mat;

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (arma::uword k = 0; k < col_idx.n_elem; ++k)
  {
    const arma::uword c = col_idx(k);
    const double *col_ptr = mat.colptr(c);
    double *out_ptr = out.colptr(c);

    double sum = 0.0;
    arma::uword n = 0;

    // compute mean
    for (arma::uword i = 0; i < mat.n_rows; ++i)
    {
      if (std::isfinite(col_ptr[i]))
      {
        sum += col_ptr[i];
        ++n;
      }
    }

    // impute only if mean can be calculated
    if (n > 0)
    {
      double mean_val = sum / static_cast<double>(n);
      for (arma::uword i = 0; i < mat.n_rows; ++i)
      {
        if (!std::isfinite(col_ptr[i]))
        {
          out_ptr[i] = mean_val;
        }
      }
    }
  }
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
