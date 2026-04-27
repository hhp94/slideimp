#ifndef GRAM_OPS_H
#define GRAM_OPS_H

#include <RcppArmadillo.h>
#include <cstring>
#include "pca_linalg_utils.h"

inline constexpr double SYRK_ALPHA = 1.0;
inline constexpr double SYRK_BETA = 0.0;

// ---------------------------------------------------------------------------
// shared sizing/syrk workspace.
// ---------------------------------------------------------------------------
struct GramWorkspace
{
  char trans = 'N';
  arma::blas_int nr = 0;
  arma::blas_int nc = 0;
  arma::blas_int n = 0;
  arma::blas_int k = 0;
  arma::blas_int lda = 0;
  arma::blas_int ldc = 0;

  bool init(arma::uword nrows_A, arma::uword ncols_A, bool tall)
  {
    nr = static_cast<arma::blas_int>(nrows_A);
    nc = static_cast<arma::blas_int>(ncols_A);
    trans = tall ? 'T' : 'N';
    n = tall ? nc : nr;
    k = tall ? nr : nc;
    lda = nr;
    ldc = n;
    return true;
  }
};

// ---------------------------------------------------------------------------
// cached fixed-column contribution to the Gram.
// ---------------------------------------------------------------------------
struct GramCache
{
  arma::uword n_fixed = 0;
  bool tall = false;
  bool active = false;
  arma::mat Gram_fixed;
  arma::mat X_chg_scaled;
};

// ---------------------------------------------------------------------------
// LAPACK syrk_upper wrapper.
// ---------------------------------------------------------------------------
inline void syrk_upper(const arma::mat &A, arma::mat &C,
                       const GramWorkspace &ws)
{
  const char uplo = 'U';
  arma::dsyrk_(&uplo, &ws.trans, &ws.n, &ws.k,
               &SYRK_ALPHA, A.memptr(), &ws.lda,
               &SYRK_BETA, C.memptr(), &ws.ldc,
               1, 1);
}

// ---------------------------------------------------------------------------
// initialize the fixed-block cache.
// ---------------------------------------------------------------------------
inline bool gram_cache_init(GramCache &cache,
                            const arma::mat &Xhat_perm,
                            const double *sw,
                            const arma::uword nrX,
                            const arma::uword n_fixed,
                            const arma::uword n_mc,
                            const bool tall,
                            arma::mat &X_work)
{
  cache.n_fixed = n_fixed;
  cache.tall = tall;
  cache.active = (n_fixed > 0);
  if (n_fixed == 0)
  {
    return true;
  }

  if (tall)
  {
    pca_detail::copy_scale_rows(X_work.memptr(), Xhat_perm.memptr(), sw, nrX, n_fixed);

    // zero-fill so the lower triangle is well-defined. dsyrk only writes the
    // upper triangle
    cache.Gram_fixed.zeros(n_fixed, n_fixed);
    const arma::blas_int nn = static_cast<arma::blas_int>(n_fixed);
    const arma::blas_int kk = static_cast<arma::blas_int>(nrX);
    const arma::blas_int lda = static_cast<arma::blas_int>(nrX);
    const arma::blas_int ldc = nn;
    const char uplo = 'U';
    const char trans = 'T';
    arma::dsyrk_(&uplo, &trans, &nn, &kk,
                 &SYRK_ALPHA, X_work.memptr(), &lda,
                 &SYRK_BETA, cache.Gram_fixed.memptr(), &ldc,
                 1, 1);
  }
  else
  {
    arma::mat X_fixed_scaled(nrX, n_fixed);
    pca_detail::copy_scale_rows(X_fixed_scaled.memptr(), Xhat_perm.memptr(), sw,
                                nrX, n_fixed);

    // zero-fill so the lower triangle is well-defined.
    cache.Gram_fixed.zeros(nrX, nrX);
    const arma::blas_int nn = static_cast<arma::blas_int>(nrX);
    const arma::blas_int kk = static_cast<arma::blas_int>(n_fixed);
    const arma::blas_int lda = static_cast<arma::blas_int>(nrX);
    const arma::blas_int ldc = nn;
    const char uplo = 'U';
    const char trans = 'N';
    arma::dsyrk_(&uplo, &trans, &nn, &kk,
                 &SYRK_ALPHA, X_fixed_scaled.memptr(), &lda,
                 &SYRK_BETA, cache.Gram_fixed.memptr(), &ldc,
                 1, 1);

    cache.X_chg_scaled.set_size(nrX, n_mc);
  }
  return true;
}

// ---------------------------------------------------------------------------
// form the (upper-triangular) weighted Gram matrix into AA_NxN.
//
// tall (n = ncols): AA = Xhat^T diag(sw^2) Xhat
//  - X_work[:, n_fixed:] is filled with sw-scaled changing cols.
//  - When cache.active, Gram_fixed is copied into the top-left block and
//  dgemm/dsyrk fill the cross and changing-changing blocks.
// wide (n = nrows): AA = Xhat diag(sw^2) Xhat^T
//  - When cache.active, Gram_fixed is the base and dsyrk(beta=1) adds the
//  contribution of the changing columns. cache.X_chg_scaled is filled.
// ---------------------------------------------------------------------------
inline void form_weighted_gram(const arma::mat &Xhat,
                               const double *sw,
                               const bool tall,
                               arma::mat &X_work,
                               arma::mat &AA_NxN,
                               const GramWorkspace &gram_ws,
                               GramCache &gram_cache)
{
  const arma::uword nr = static_cast<arma::uword>(gram_ws.nr);
  const arma::uword nc = static_cast<arma::uword>(gram_ws.nc);
  const arma::uword n_fixed = gram_cache.n_fixed;
  const arma::uword n_mc = nc - n_fixed;
  const bool use_cache = gram_cache.active;
  if (tall)
  {
    // X_work already contains sw * Xhat.
    // fixed columns initialized in cache init.
    // changing columns maintained by impute_restandardize().
    if (use_cache)
    {
      const arma::blas_int ldAA = gram_ws.ldc;
      const arma::blas_int lda = gram_ws.lda;

      for (arma::uword j = 0; j < n_fixed; ++j)
      {
        std::memcpy(AA_NxN.colptr(j),
                    gram_cache.Gram_fixed.colptr(j),
                    (j + 1) * sizeof(double));
      }

      if (n_mc > 0)
      {
        // fixed-changing cross block
        {
          const arma::blas_int m = static_cast<arma::blas_int>(n_fixed);
          const arma::blas_int n_rhs = static_cast<arma::blas_int>(n_mc);
          const arma::blas_int kk = static_cast<arma::blas_int>(nr);
          const char trA = 'T';
          const char trB = 'N';

          arma::dgemm_(&trA, &trB, &m, &n_rhs, &kk,
                       &SYRK_ALPHA,
                       X_work.memptr(), &lda,
                       X_work.colptr(n_fixed), &lda,
                       &SYRK_BETA,
                       AA_NxN.colptr(n_fixed), &ldAA,
                       1, 1);
        }
        // changing-changing block
        {
          const arma::blas_int nn = static_cast<arma::blas_int>(n_mc);
          const arma::blas_int kk = static_cast<arma::blas_int>(nr);
          const char uplo = 'U';
          const char trans = 'T';
          double *C_sub = AA_NxN.colptr(n_fixed) + n_fixed;

          arma::dsyrk_(&uplo, &trans, &nn, &kk,
                       &SYRK_ALPHA,
                       X_work.colptr(n_fixed), &lda,
                       &SYRK_BETA,
                       C_sub, &ldAA,
                       1, 1);
        }
      }
    }
    else
    {
      syrk_upper(X_work, AA_NxN, gram_ws);
    }
  }
  else
  {
    if (use_cache)
    {
      // gram_cache.X_chg_scaled already contains sw * X_changing.
      for (arma::uword j = 0; j < nr; ++j)
      {
        std::memcpy(AA_NxN.colptr(j),
                    gram_cache.Gram_fixed.colptr(j),
                    (j + 1) * sizeof(double));
      }

      if (n_mc > 0)
      {
        const arma::blas_int nn = gram_ws.n;
        const arma::blas_int kk = static_cast<arma::blas_int>(n_mc);
        const arma::blas_int lda = static_cast<arma::blas_int>(nr);
        const arma::blas_int ldc = gram_ws.ldc;
        const char uplo = 'U';
        const char trans = 'N';
        const double beta_one = 1.0;

        arma::dsyrk_(&uplo, &trans, &nn, &kk,
                     &SYRK_ALPHA,
                     gram_cache.X_chg_scaled.memptr(), &lda,
                     &beta_one,
                     AA_NxN.memptr(), &ldc,
                     1, 1);
      }
    }
    else
    {
      syrk_upper(Xhat, AA_NxN, gram_ws);

      double *M = AA_NxN.memptr();
      const arma::uword ng = AA_NxN.n_rows;

      for (arma::uword j = 0; j < ng; ++j)
      {
        const double sj = sw[j];
        double *col = M + j * ng;

        for (arma::uword i = 0; i <= j; ++i)
        {
          col[i] *= sw[i] * sj;
        }
      }
    }
  }
}

#endif // GRAM_OPS_H
