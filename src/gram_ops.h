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

  void init(arma::uword nrows_A, arma::uword ncols_A, bool tall)
  {
    nr = static_cast<arma::blas_int>(nrows_A);
    nc = static_cast<arma::blas_int>(ncols_A);
    trans = tall ? 'T' : 'N';
    n = tall ? nc : nr;
    k = tall ? nr : nc;
    lda = nr;
    ldc = n;
  }
};

// ---------------------------------------------------------------------------
// cached fixed-column contribution to the Gram.
// ---------------------------------------------------------------------------
struct GramCache
{
  arma::uword n_fixed = 0;
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
inline void gram_cache_init(GramCache &cache,
                            const arma::mat &Xhat_perm,
                            const double *sw,
                            const arma::uword nrX,
                            const arma::uword n_fixed,
                            const arma::uword n_mc,
                            const bool tall,
                            arma::mat &X_work)
{
  cache.n_fixed = n_fixed;
  cache.active = (n_fixed > 0);
  if (n_fixed == 0)
  {
    return;
  }

  GramWorkspace fixed_ws;
  fixed_ws.init(nrX, n_fixed, tall);

  if (tall)
  {
    pca_detail::copy_scale_rows(X_work.memptr(), Xhat_perm.memptr(), sw, nrX, n_fixed);

    // zero-fill so the lower triangle is well-defined. dsyrk only writes the
    // upper triangle
    cache.Gram_fixed.zeros(n_fixed, n_fixed);
    syrk_upper(X_work, cache.Gram_fixed, fixed_ws);
  }
  else
  {
    arma::mat X_fixed_scaled(nrX, n_fixed);
    pca_detail::copy_scale_rows(X_fixed_scaled.memptr(), Xhat_perm.memptr(), sw,
                                nrX, n_fixed);

    // zero-fill so the lower triangle is well-defined.
    cache.Gram_fixed.zeros(nrX, nrX);
    syrk_upper(X_fixed_scaled, cache.Gram_fixed, fixed_ws);

    cache.X_chg_scaled.set_size(nrX, n_mc);
  }
}

// ---------------------------------------------------------------------------
// where the sw-scaled changing-column block lives, if anywhere. The caller
// must keep that block scaled by sw between form_weighted_gram() calls.
//   tall:            the tail columns of X_work
//   wide, cached:    the compact GramCache::X_chg_scaled block
//   wide, non-cached: nowhere (form_weighted_gram scales the Gram instead)
// ---------------------------------------------------------------------------
inline double *gram_changing_block(GramCache &cache,
                                   arma::mat &X_work,
                                   const arma::uword n_fixed,
                                   const arma::uword n_mc,
                                   const bool tall)
{
  if (n_mc == 0)
  {
    return nullptr;
  }
  if (tall)
  {
    return X_work.colptr(n_fixed);
  }
  if (cache.active)
  {
    return cache.X_chg_scaled.memptr();
  }
  return nullptr;
}

// ---------------------------------------------------------------------------
// copy the upper-triangular part (rows 0..j) of the leading ncols columns.
// ---------------------------------------------------------------------------
inline void copy_upper_cols(arma::mat &dst, const arma::mat &src,
                            const arma::uword ncols)
{
  for (arma::uword j = 0; j < ncols; ++j)
  {
    std::memcpy(dst.colptr(j), src.colptr(j), (j + 1) * sizeof(double));
  }
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

      copy_upper_cols(AA_NxN, gram_cache.Gram_fixed, n_fixed);

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
      copy_upper_cols(AA_NxN, gram_cache.Gram_fixed, nr);

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
