#ifndef SLIDEIMP_GRAM_OPS_H
#define SLIDEIMP_GRAM_OPS_H

// The weighted Gram of the PCA chain, formed with R's BLAS through
// src/si_lapack.h. Matrices are raw column-major buffers with std::size_t
// dimensions; nothing here uses an Armadillo, Eigen or R type. A function
// that calls BLAS returns the wrapper's status, si_la::ok on success, and
// raising the R error is the caller's.

#include "si_lapack.h"
#include "pca_linalg_utils.h"
#include <cmath>
#include <cstddef>
#include <cstring>
#include <vector>

inline constexpr double SYRK_ALPHA = 1.0;
inline constexpr double SYRK_BETA = 0.0;

// ---------------------------------------------------------------------------
// the dimensions of a Gram and of the dsyrk that forms it, no storage: A is
// nr x nc, stored with leading dimension nr, and the Gram is n x n, stored
// with leading dimension n. When tall, n = nc and the Gram is A'A;
// otherwise n = nr and it is AA'. k is the inner dimension of the product,
// and trans the dsyrk option that gives it.
// ---------------------------------------------------------------------------
struct GramDims {
  char trans = 'N';
  bool tall = false;
  std::size_t nr = 0;
  std::size_t nc = 0;
  std::size_t n = 0;
  std::size_t k = 0;

  void init(const std::size_t nrows_A, const std::size_t ncols_A,
            const bool tall_A) {
    nr = nrows_A;
    nc = ncols_A;
    tall = tall_A;
    trans = tall ? 'T' : 'N';
    n = tall ? nc : nr;
    k = tall ? nr : nc;
  }
};

// ---------------------------------------------------------------------------
// cached fixed-column contribution to the Gram.
//   Gram_fixed:   the fixed columns' Gram, upper triangle written, lower
//                 zero; n_fixed x n_fixed (tall) or nr x nr (wide)
//   X_chg_scaled: nr x n_mc, the sw-scaled changing columns; wide and
//                 cached only, empty otherwise
// Both are empty while the cache is inactive.
// ---------------------------------------------------------------------------
struct GramCache {
  std::size_t n_fixed = 0;
  bool active = false;
  std::vector<double> Gram_fixed;
  std::vector<double> X_chg_scaled;
};

// ---------------------------------------------------------------------------
// dsyrk into the upper triangle of C: C = A'A (tall) or AA' (wide), with the
// dimensions of dims.
// ---------------------------------------------------------------------------
[[nodiscard]] inline int syrk_upper(const double *A, double *C,
                                    const GramDims &dims) {
  return si_la::syrk('U', dims.trans, dims.n, dims.k, SYRK_ALPHA, A, dims.nr,
                     SYRK_BETA, C, dims.n);
}

// ---------------------------------------------------------------------------
// initialize the fixed-block cache. Xhat_perm is nrX x (n_fixed + n_mc),
// fixed columns first. In the tall case X_work (nrX x (n_fixed + n_mc))
// receives the sw-scaled fixed columns; the wide case does not touch it.
// ---------------------------------------------------------------------------
[[nodiscard]] inline int
gram_cache_init(GramCache &cache, const double *Xhat_perm, const double *sw,
                const std::size_t nrX, const std::size_t n_fixed,
                const std::size_t n_mc, const bool tall, double *X_work) {
  cache.n_fixed = n_fixed;
  cache.active = (n_fixed > 0);
  if (n_fixed == 0) {
    return si_la::ok;
  }

  GramDims fixed_dims;
  fixed_dims.init(nrX, n_fixed, tall);

  if (tall) {
    pca_detail::copy_scale_rows(X_work, Xhat_perm, sw, nrX, n_fixed);

    // zero-fill so the lower triangle is well-defined. dsyrk only writes the
    // upper triangle
    cache.Gram_fixed.assign(n_fixed * n_fixed, 0.0);
    return syrk_upper(X_work, cache.Gram_fixed.data(), fixed_dims);
  }

  std::vector<double> X_fixed_scaled(nrX * n_fixed);
  pca_detail::copy_scale_rows(X_fixed_scaled.data(), Xhat_perm, sw, nrX,
                              n_fixed);

  // zero-fill so the lower triangle is well-defined.
  cache.Gram_fixed.assign(nrX * nrX, 0.0);
  const int st =
      syrk_upper(X_fixed_scaled.data(), cache.Gram_fixed.data(), fixed_dims);

  cache.X_chg_scaled.assign(nrX * n_mc, 0.0);
  return st;
}

// ---------------------------------------------------------------------------
// where the sw-scaled changing-column block lives, if anywhere. The caller
// must keep that block scaled by sw between form_weighted_gram() calls.
//   tall:            the tail columns of X_work (nrX rows)
//   wide, cached:    the compact GramCache::X_chg_scaled block
//   wide, non-cached: nowhere (form_weighted_gram scales the Gram instead)
// ---------------------------------------------------------------------------
inline double *gram_changing_block(GramCache &cache, double *X_work,
                                   const std::size_t nrX,
                                   const std::size_t n_fixed,
                                   const std::size_t n_mc, const bool tall) {
  if (n_mc == 0) {
    return nullptr;
  }
  if (tall) {
    return X_work + n_fixed * nrX;
  }
  if (cache.active) {
    return cache.X_chg_scaled.data();
  }
  return nullptr;
}

// ---------------------------------------------------------------------------
// copy the upper-triangular part (rows 0..j) of the leading ncols columns of
// src, leading dimension ld_src, into dst, leading dimension ld_dst.
// ---------------------------------------------------------------------------
inline void copy_upper_cols(double *dst, const std::size_t ld_dst,
                            const double *src, const std::size_t ld_src,
                            const std::size_t ncols) {
  for (std::size_t j = 0; j < ncols; ++j) {
    std::memcpy(dst + j * ld_dst, src + j * ld_src, (j + 1) * sizeof(double));
  }
}

// ---------------------------------------------------------------------------
// form the (upper-triangular) weighted Gram matrix into AA_NxN, which is
// gram_dims.n x gram_dims.n.
//
// tall (n = ncols): AA = Xhat^T diag(sw^2) Xhat
//  - X_work[:, n_fixed:] is filled with sw-scaled changing cols.
//  - When cache.active, Gram_fixed is copied into the top-left block and
//  dgemm/dsyrk fill the cross and changing-changing blocks.
// wide (n = nrows): AA = Xhat diag(sw^2) Xhat^T
//  - When cache.active, Gram_fixed is the base and dsyrk(beta=1) adds the
//  contribution of the changing columns. cache.X_chg_scaled is filled.
//
// Xhat and X_work are gram_dims.nr x gram_dims.nc; X_work is read in the
// tall case only, Xhat in the wide non-cached case only.
// ---------------------------------------------------------------------------
[[nodiscard]] inline int
form_weighted_gram(const double *Xhat, const double *sw, const double *X_work,
                   double *AA_NxN, const GramDims &gram_dims,
                   const GramCache &gram_cache) {
  const std::size_t nr = gram_dims.nr;
  const std::size_t nc = gram_dims.nc;
  const std::size_t n = gram_dims.n;
  const std::size_t n_fixed = gram_cache.n_fixed;
  const std::size_t n_mc = nc - n_fixed;
  const bool use_cache = gram_cache.active;
  if (gram_dims.tall) {
    // X_work already contains sw * Xhat.
    // fixed columns initialized in cache init.
    // changing columns maintained by impute_restandardize().
    if (!use_cache) {
      return syrk_upper(X_work, AA_NxN, gram_dims);
    }

    copy_upper_cols(AA_NxN, n, gram_cache.Gram_fixed.data(), n_fixed, n_fixed);

    if (n_mc == 0) {
      return si_la::ok;
    }

    const double *X_chg = X_work + n_fixed * nr;
    double *AA_chg_cols = AA_NxN + n_fixed * n;

    // fixed-changing cross block
    const int st = si_la::gemm('T', 'N', n_fixed, n_mc, nr, SYRK_ALPHA, X_work,
                               nr, X_chg, nr, SYRK_BETA, AA_chg_cols, n);
    if (st != si_la::ok) {
      return st;
    }

    // changing-changing block
    return si_la::syrk('U', 'T', n_mc, nr, SYRK_ALPHA, X_chg, nr, SYRK_BETA,
                       AA_chg_cols + n_fixed, n);
  }

  if (use_cache) {
    // gram_cache.X_chg_scaled already contains sw * X_changing.
    copy_upper_cols(AA_NxN, n, gram_cache.Gram_fixed.data(), nr, nr);

    if (n_mc == 0) {
      return si_la::ok;
    }

    const double beta_one = 1.0;
    return si_la::syrk('U', 'N', n, n_mc, SYRK_ALPHA,
                       gram_cache.X_chg_scaled.data(), nr, beta_one, AA_NxN, n);
  }

  const int st = syrk_upper(Xhat, AA_NxN, gram_dims);
  if (st != si_la::ok) {
    return st;
  }

  for (std::size_t j = 0; j < n; ++j) {
    const double sj = sw[j];
    double *col = AA_NxN + j * n;

    for (std::size_t i = 0; i <= j; ++i) {
      col[i] *= sw[i] * sj;
    }
  }
  return si_la::ok;
}

// ---------------------------------------------------------------------------
// whether the upper triangle of the n x n Gram, diagonal included, holds only
// finite values: the triangle form_weighted_gram() writes. The input is
// finite, so a NaN or Inf here is a cross-product that overflowed.
// ---------------------------------------------------------------------------
[[nodiscard]] inline bool gram_upper_finite(const double *A,
                                            const std::size_t n) {
  for (std::size_t j = 0; j < n; ++j) {
    const double *col = A + j * n;
    for (std::size_t i = 0; i <= j; ++i) {
      if (!std::isfinite(col[i])) {
        return false;
      }
    }
  }
  return true;
}

#endif // SLIDEIMP_GRAM_OPS_H
