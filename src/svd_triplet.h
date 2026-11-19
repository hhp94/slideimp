#ifndef SVD_TRIPLET_H
#define SVD_TRIPLET_H

#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include "gram_ops.h"
#include "eig_sym_sel.h"
#include "hybrid_topk_eig.h"
#include "pca_linalg_utils.h"
#include "loc_timer.h"

// ---------------------------------------------------------------------------
// copy the leading ncp columns of the n-row eigvecs into dst while writing
// the d_inv-scaled (and, when row_scale is non-null, row-scaled) copy into
// Vn. dst and Vn are n x ncp.
// ---------------------------------------------------------------------------
static inline void split_eigvec_block(const double *eigvecs,
                                      const std::size_t n,
                                      double *dst,
                                      double *Vn,
                                      const double *d_inv,
                                      const double *row_scale,
                                      const std::size_t ncp)
{
  for (std::size_t j = 0; j < ncp; ++j)
  {
    const double dj = d_inv[j];
    const double *ecol = eigvecs + j * n;
    double *dcol = dst + j * n;
    double *ncol = Vn + j * n;
    if (row_scale)
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        dcol[i] = e;
        ncol[i] = e * dj * row_scale[i];
      }
    }
    else
    {
      for (std::size_t i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        dcol[i] = e;
        ncol[i] = e * dj;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// the sum of the diagonal of the n x n A, as arma::trace() took it
// (fn_trace.hpp): two partial sums, over the even and the odd cells, added
// at the end. dev/CPP.md bans such a reduction; this one is kept by the
// maintainer's ruling, because a plain sum moves results through the EM
// loop (dev/eigen_conversion.md, step 5 of Phase 2), and Phase 4 revisits
// it with the rest of Armadillo's choices.
// ---------------------------------------------------------------------------
static inline double gram_trace(const double *A, const std::size_t n)
{
  double val1 = 0.0;
  double val2 = 0.0;
  std::size_t i = 0;
  std::size_t j = 1;
  for (; j < n; i += 2, j += 2)
  {
    val1 += A[i + i * n];
    val2 += A[j + j * n];
  }
  if (i < n)
  {
    val1 += A[i + i * n];
  }
  return val1 + val2;
}

inline constexpr char GRAM_NOT_FINITE_MSG[] =
    "the cross-product of the data is not finite; values this large overflow "
    "when squared, rescale the input";

// ---------------------------------------------------------------------------
// SVD_triplet
//
// U is returned in sqrt-row-weight space; the caller removes the weighting
// (fused with its own column scaling of U).
//
// In the wide case V is written only for the changing columns, rows
// n_fixed .. n_elig-1; the caller fills the fixed rows once after the EM
// loop from the same Vn. In the tall case V is complete.
//
// Every buffer is the caller's and sized by it, with nr, nc and n those of
// gram_ws and topk that of hyb_ctx.eig_range: Xhat nr x nc; X_work nr x nc
// in the tall case, unread in the wide case; AA_NxN n x n; U nr x ncp; V
// nc x ncp; Vn n x ncp; eigvecs n x topk; eigvals and vs_top topk; d_inv
// ncp.
// ---------------------------------------------------------------------------
inline void SVD_triplet(const double *Xhat,
                        const double *sw,
                        const std::size_t ncp,
                        const bool tall,
                        double *vs_top,
                        double &trace_val,
                        double *U,
                        double *V,
                        const double *X_work,
                        double *AA_NxN,
                        double *d_inv,
                        double *Vn,
                        double *eigvals,
                        double *eigvecs,
                        const GramWorkspace &gram_ws,
                        const GramCache &gram_cache,
                        HybridEigContext &hyb_ctx,
                        int outer_iter
                            LOC_TIMER_PARAM(timer))
{
  const std::size_t nr = gram_ws.nr;
  const std::size_t nc = gram_ws.nc;
  const std::size_t n = gram_ws.n;

  LOC_TIC(timer, "form_gram");
  const int gram_st = form_weighted_gram(Xhat, sw, tall, X_work, AA_NxN,
                                         gram_ws, gram_cache);
  LOC_TOC(timer, "form_gram");
  if (gram_st != si_la::ok)
  {
    Rcpp::stop("internal error: dsyrk or dgemm refused its arguments "
               "(status %d)",
               gram_st);
  }

  // The one test for NaN and Inf in the chain's linear algebra, before
  // either solver sees the Gram.
  if (!gram_upper_finite(AA_NxN, n))
  {
    Rcpp::stop(GRAM_NOT_FINITE_MSG);
  }

  trace_val = gram_trace(AA_NxN, n);

  LOC_TIC(timer, "eig");
  const int eig_st =
      hybrid_topk_eig(eigvals, eigvecs, AA_NxN, n, n, hyb_ctx, outer_iter);
  LOC_TOC(timer, "eig");
  if (eig_st > 0)
  {
    Rcpp::stop("the symmetric eigensolver failed (LAPACK dsyevr, info = %d)",
               eig_st);
  }
  if (eig_st == si_la::err_nonfinite)
  {
    // The Gram passed the test above, so this is a value that became NaN or
    // Inf between it and dsyevr.
    Rcpp::stop(GRAM_NOT_FINITE_MSG);
  }
  if (eig_st == si_la::err_alloc)
  {
    Rcpp::stop("could not allocate the workspace of dsyevr");
  }
  if (eig_st != si_la::ok)
  {
    Rcpp::stop("internal error: the symmetric eigensolver returned status %d",
               eig_st);
  }

  const std::size_t ke = hyb_ctx.eig_range->topk;
  for (std::size_t i = 0; i < ke; ++i)
  {
    const double e = eigvals[i];
    vs_top[i] = (e > 0.0) ? std::sqrt(e) : 0.0;
  }
  for (std::size_t i = 0; i < ncp; ++i)
  {
    d_inv[i] = (vs_top[i] > PCA_TOL) ? (1.0 / vs_top[i]) : 0.0;
  }

  LOC_TIC(timer, "recover_factors");
  if (tall)
  {
    split_eigvec_block(eigvecs, n, V, Vn, d_inv, nullptr, ncp);
    const int st = pca_detail::gemm_nn(X_work, nr, nc, Vn, n, ncp, U);
    if (st != si_la::ok)
    {
      Rcpp::stop("internal error: dgemm refused its arguments (status %d)", st);
    }
  }
  else
  {
    split_eigvec_block(eigvecs, n, U, Vn, d_inv, sw, ncp);
    // inside the EM loop only the changing columns' rows of V are read, and
    // the permutation keeps those columns in the tail block of Xhat, so one
    // dgemm over that block fills rows n_fixed .. n_elig-1 in place. The
    // caller sized V once, before the EM loop, so the rows written earlier
    // stay where they are.
    const std::size_t n_fixed = gram_cache.n_fixed;
    const int st = pca_detail::gemm_tn_block(Xhat, nr, nc, n_fixed,
                                             nc - n_fixed, Vn, n, ncp, V);
    if (st != si_la::ok)
    {
      Rcpp::stop("internal error: dgemm refused its arguments (status %d)", st);
    }
  }
  LOC_TOC(timer, "recover_factors");
}

#endif // SVD_TRIPLET_H
