#ifndef PCA_LINALG_UTILS_H
#define PCA_LINALG_UTILS_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstring>

inline constexpr double PCA_TOL = 1e-15;

namespace pca_detail
{
  // sklearn svd_flip convention: for each column, force the largest-magnitude
  // entry to be non-negative.
  static inline void canonicalize_signs(arma::mat &V)
  {
    const arma::uword n = V.n_rows, k = V.n_cols;
    for (arma::uword j = 0; j < k; ++j)
    {
      double *col = V.colptr(j);
      arma::uword imax = 0;
      double amax = std::abs(col[0]);
      for (arma::uword i = 1; i < n; ++i)
      {
        const double a = std::abs(col[i]);
        if (a > amax)
        {
          amax = a;
          imax = i;
        }
      }
      if (col[imax] < 0.0)
      {
        for (arma::uword i = 0; i < n; ++i)
          col[i] = -col[i];
      }
    }
  }

  static inline void scale_cols_inplace(double *M, const double *s,
                                        arma::uword nr, arma::uword nc)
  {
    for (arma::uword j = 0; j < nc; ++j)
    {
      const double sj = s[j];
      double *col = M + j * nr;
      for (arma::uword i = 0; i < nr; ++i)
      {
        col[i] *= sj;
      }
    }
  }

  static inline void scale_rows_inplace(double *M, const double *w,
                                        arma::uword nr, arma::uword nc)
  {
    for (arma::uword j = 0; j < nc; ++j)
    {
      double *col = M + j * nr;
      for (arma::uword i = 0; i < nr; ++i)
      {
        col[i] *= w[i];
      }
    }
  }

  static inline void copy_scale_rows(double *dst, const double *src,
                                     const double *w,
                                     arma::uword nr, arma::uword nc)
  {
    for (arma::uword j = 0; j < nc; ++j)
    {
      double *dcol = dst + j * nr;
      const double *scol = src + j * nr;
      for (arma::uword i = 0; i < nr; ++i)
      {
        dcol[i] = scol[i] * w[i];
      }
    }
  }

  static inline void mirror_upper_to_lower(arma::mat &A)
  {
    const arma::uword n = A.n_rows;
    double *M = A.memptr();
    for (arma::uword j = 0; j < n; ++j)
    {
      for (arma::uword i = 0; i < j; ++i)
      {
        M[j + i * n] = M[i + j * n];
      }
    }
  }

  static inline void gemm_nn(const arma::mat &A,
                             const arma::mat &B,
                             arma::mat &C)
  {
    C.set_size(A.n_rows, B.n_cols);

    if (A.n_rows == 0 || B.n_cols == 0)
      return;

    if (A.n_cols == 0)
    {
      C.zeros();
      return;
    }

    const char trN = 'N';
    const double alpha = 1.0;
    const double beta = 0.0;

    const arma::blas_int m = static_cast<arma::blas_int>(A.n_rows);
    const arma::blas_int n = static_cast<arma::blas_int>(B.n_cols);
    const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);

    const arma::blas_int lda = static_cast<arma::blas_int>(A.n_rows);
    const arma::blas_int ldb = static_cast<arma::blas_int>(B.n_rows);
    const arma::blas_int ldc = static_cast<arma::blas_int>(C.n_rows);

    arma::dgemm_(&trN, &trN,
                 &m, &n, &k,
                 &alpha,
                 A.memptr(), &lda,
                 B.memptr(), &ldb,
                 &beta,
                 C.memptr(), &ldc,
                 1, 1);
  }

  static inline void gemm_tn(const arma::mat &A,
                             const arma::mat &B,
                             arma::mat &C)
  {
    // C = A.t() * B
    C.set_size(A.n_cols, B.n_cols);

    if (A.n_cols == 0 || B.n_cols == 0)
      return;

    if (A.n_rows == 0)
    {
      C.zeros();
      return;
    }

    const char trT = 'T';
    const char trN = 'N';
    const double alpha = 1.0;
    const double beta = 0.0;

    const arma::blas_int m = static_cast<arma::blas_int>(A.n_cols);
    const arma::blas_int n = static_cast<arma::blas_int>(B.n_cols);
    const arma::blas_int k = static_cast<arma::blas_int>(A.n_rows);

    const arma::blas_int lda = static_cast<arma::blas_int>(A.n_rows);
    const arma::blas_int ldb = static_cast<arma::blas_int>(B.n_rows);
    const arma::blas_int ldc = static_cast<arma::blas_int>(C.n_rows);

    arma::dgemm_(&trT, &trN,
                 &m, &n, &k,
                 &alpha,
                 A.memptr(), &lda,
                 B.memptr(), &ldb,
                 &beta,
                 C.memptr(), &ldc,
                 1, 1);
  }

  static inline void gemm_nt(const arma::mat &A,
                             const arma::mat &B,
                             arma::mat &C)
  {
    // C = A * B.t()
    C.set_size(A.n_rows, B.n_rows);

    if (A.n_rows == 0 || B.n_rows == 0)
      return;

    if (A.n_cols == 0)
    {
      C.zeros();
      return;
    }

    const char trN = 'N';
    const char trT = 'T';
    const double alpha = 1.0;
    const double beta = 0.0;

    const arma::blas_int m = static_cast<arma::blas_int>(A.n_rows);
    const arma::blas_int n = static_cast<arma::blas_int>(B.n_rows);
    const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);

    const arma::blas_int lda = static_cast<arma::blas_int>(A.n_rows);
    const arma::blas_int ldb = static_cast<arma::blas_int>(B.n_rows);
    const arma::blas_int ldc = static_cast<arma::blas_int>(C.n_rows);

    arma::dgemm_(&trN, &trT,
                 &m, &n, &k,
                 &alpha,
                 A.memptr(), &lda,
                 B.memptr(), &ldb,
                 &beta,
                 C.memptr(), &ldc,
                 1, 1);
  }

} // namespace pca_detail

#endif // PCA_LINALG_UTILS_H
