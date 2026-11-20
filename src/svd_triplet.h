#ifndef SLIDEIMP_SVD_TRIPLET_H
#define SLIDEIMP_SVD_TRIPLET_H

#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <vector>
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
inline void split_eigvec_block(const double *eigvecs, const std::size_t n,
                               double *dst, double *Vn, const double *d_inv,
                               const double *row_scale, const std::size_t ncp) {
  for (std::size_t j = 0; j < ncp; ++j) {
    const double dj = d_inv[j];
    const double *ecol = eigvecs + j * n;
    double *dcol = dst + j * n;
    double *ncol = Vn + j * n;
    if (row_scale) {
      for (std::size_t i = 0; i < n; ++i) {
        const double e = ecol[i];
        dcol[i] = e;
        ncol[i] = e * dj * row_scale[i];
      }
    } else {
      for (std::size_t i = 0; i < n; ++i) {
        const double e = ecol[i];
        dcol[i] = e;
        ncol[i] = e * dj;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// the sum of the diagonal of the n x n A, in order
// ---------------------------------------------------------------------------
inline double gram_trace(const double *A, const std::size_t n) {
  double trace = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    trace += A[i + i * n];
  }
  return trace;
}

inline constexpr char GRAM_NOT_FINITE_MSG[] =
    "the cross-product of the data is not finite; values this large overflow "
    "when squared, rescale the input";

// ---------------------------------------------------------------------------
// stops with an internal error unless st, the si_la status of a dgemm, is
// si_la::ok. Kernel setup refuses every dimension past int, so no call is
// expected to stop here.
// ---------------------------------------------------------------------------
inline void stop_if_gemm_refused(const int st) {
  if (st != si_la::ok) {
    Rcpp::stop("internal error: dgemm refused its arguments (status %d)", st);
  }
}

// ---------------------------------------------------------------------------
// the buffers SVD_triplet() reads and writes, sized once by the constructor,
// before the EM loop, and zero-filled; every call reuses them. With nr, nc
// and n those of the GramDims and topk the number of eigenpairs asked for:
// AA_NxN n x n; X_work nr x nc in the tall case, empty in the wide case; U
// nr x ncp; V nc x ncp; Vn n x ncp; eigvecs n x topk; eigvals and vs_top
// topk; d_inv ncp.
//
// In the tall case the caller keeps X_work holding the sw-scaled Xhat
// (gram_cache_init(), gram_changing_block()). After a call U, V, Vn, vs_top
// and eigvecs hold its results for the caller; eigvals and d_inv are
// SVD_triplet()'s own, and the caller does not read them.
// ---------------------------------------------------------------------------
struct SolveBuffers {
  std::vector<double> AA_NxN;
  std::vector<double> X_work;
  std::vector<double> U;
  std::vector<double> V;
  std::vector<double> Vn;
  std::vector<double> eigvecs;
  std::vector<double> eigvals;
  std::vector<double> vs_top;
  std::vector<double> d_inv;

  SolveBuffers(const GramDims &dims, const std::size_t ncp,
               const std::size_t topk)
      : AA_NxN(dims.n * dims.n), X_work(dims.tall ? dims.nr * dims.nc : 0),
        U(dims.nr * ncp), V(dims.nc * ncp), Vn(dims.n * ncp),
        eigvecs(dims.n * topk), eigvals(topk), vs_top(topk), d_inv(ncp) {}
};

// ---------------------------------------------------------------------------
// SVD_triplet
//
// U is returned in sqrt-row-weight space; the caller removes the weighting
// (fused with its own column scaling of U).
//
// In the wide case V is written only for the changing columns, rows
// n_fixed .. nc-1; the caller fills the fixed rows once after the EM
// loop from the same Vn. In the tall case V is complete.
//
// Xhat is gram_dims.nr x gram_dims.nc, and buf was sized for gram_dims, ncp
// and the topk of hyb_ctx.eig_range. Returns the trace of the Gram, taken
// before the eigensolver overwrites it.
// ---------------------------------------------------------------------------
inline double SVD_triplet(const double *Xhat, const double *sw,
                          const std::size_t ncp, SolveBuffers &buf,
                          const GramDims &gram_dims,
                          const GramCache &gram_cache,
                          HybridEigContext &hyb_ctx,
                          int outer_iter LOC_TIMER_PARAM(timer)) {
  const std::size_t nr = gram_dims.nr;
  const std::size_t nc = gram_dims.nc;
  const std::size_t n = gram_dims.n;
  double *AA_NxN = buf.AA_NxN.data();
  const double *X_work = buf.X_work.data();
  double *U = buf.U.data();
  double *V = buf.V.data();
  double *Vn = buf.Vn.data();
  double *eigvecs = buf.eigvecs.data();
  double *eigvals = buf.eigvals.data();
  double *vs_top = buf.vs_top.data();
  double *d_inv = buf.d_inv.data();

  LOC_TIC(timer, "form_gram");
  const int gram_st =
      form_weighted_gram(Xhat, sw, X_work, AA_NxN, gram_dims, gram_cache);
  LOC_TOC(timer, "form_gram");
  if (gram_st != si_la::ok) {
    Rcpp::stop("internal error: dsyrk or dgemm refused its arguments "
               "(status %d)",
               gram_st);
  }

  // The Gram is tested for NaN and Inf here, before either solver sees it.
  // On the exact path si_la::syevr_range tests the same triangle again.
  if (!gram_upper_finite(AA_NxN, n)) {
    Rcpp::stop(GRAM_NOT_FINITE_MSG);
  }

  const double trace_val = gram_trace(AA_NxN, n);

  LOC_TIC(timer, "eig");
  const int eig_st =
      hybrid_topk_eig(eigvals, eigvecs, AA_NxN, n, hyb_ctx, outer_iter);
  LOC_TOC(timer, "eig");
  if (eig_st > 0) {
    Rcpp::stop("the symmetric eigensolver failed (LAPACK dsyevr, info = %d)",
               eig_st);
  }
  if (eig_st == si_la::err_nonfinite) {
    // The Gram passed the test above, so this is a value that became NaN or
    // Inf between it and dsyevr.
    Rcpp::stop(GRAM_NOT_FINITE_MSG);
  }
  if (eig_st == si_la::err_alloc) {
    Rcpp::stop("could not allocate the workspace of dsyevr");
  }
  if (eig_st != si_la::ok) {
    Rcpp::stop("internal error: the symmetric eigensolver returned status %d",
               eig_st);
  }

  const std::size_t ke = hyb_ctx.eig_range->topk;
  for (std::size_t i = 0; i < ke; ++i) {
    const double e = eigvals[i];
    vs_top[i] = (e > 0.0) ? std::sqrt(e) : 0.0;
  }
  for (std::size_t i = 0; i < ncp; ++i) {
    d_inv[i] = (vs_top[i] > PCA_TOL) ? (1.0 / vs_top[i]) : 0.0;
  }

  LOC_TIC(timer, "recover_factors");
  if (gram_dims.tall) {
    split_eigvec_block(eigvecs, n, V, Vn, d_inv, nullptr, ncp);
    const int st = pca_detail::gemm_nn(X_work, nr, nc, Vn, n, ncp, U);
    stop_if_gemm_refused(st);
  } else {
    split_eigvec_block(eigvecs, n, U, Vn, d_inv, sw, ncp);
    // inside the EM loop only the changing columns' rows of V are read, and
    // the permutation keeps those columns in the tail block of Xhat, so one
    // dgemm over that block fills rows n_fixed .. nc-1 in place. The
    // caller sized V once, before the EM loop, so the rows written earlier
    // stay where they are.
    const std::size_t n_fixed = gram_cache.n_fixed;
    const int st = pca_detail::gemm_tn_block(Xhat, nr, nc, n_fixed,
                                             nc - n_fixed, Vn, n, ncp, V);
    stop_if_gemm_refused(st);
  }
  LOC_TOC(timer, "recover_factors");
  return trace_val;
}

#endif // SLIDEIMP_SVD_TRIPLET_H
