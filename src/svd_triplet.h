#ifndef SVD_TRIPLET_H
#define SVD_TRIPLET_H

#include <RcppArmadillo.h>
#include "gram_ops.h"
#include "eig_sym_sel.h"
#include "hybrid_topk_eig.h"
#include "pca_linalg_utils.h"
#include "loc_timer.h"

// ---------------------------------------------------------------------------
// copy the leading ncp eigenvector columns into dst while writing the
// d_inv-scaled (and, when row_scale is non-null, row-scaled) copy into Vn.
// ---------------------------------------------------------------------------
static inline void split_eigvec_block(const arma::mat &eigvecs,
                                      arma::mat &dst,
                                      arma::mat &Vn,
                                      const arma::vec &d_inv,
                                      const double *row_scale,
                                      const arma::uword ncp)
{
  const arma::uword n = eigvecs.n_rows;
  dst.set_size(n, ncp);
  Vn.set_size(n, ncp);
  const double *EV = eigvecs.memptr();
  double *Dp = dst.memptr();
  double *Vnp = Vn.memptr();
  const double *dp = d_inv.memptr();
  for (arma::uword j = 0; j < ncp; ++j)
  {
    const double dj = dp[j];
    const double *ecol = EV + j * n;
    double *dcol = Dp + j * n;
    double *ncol = Vnp + j * n;
    if (row_scale)
    {
      for (arma::uword i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        dcol[i] = e;
        ncol[i] = e * dj * row_scale[i];
      }
    }
    else
    {
      for (arma::uword i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        dcol[i] = e;
        ncol[i] = e * dj;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// SVD_triplet
//
// U is returned in sqrt-row-weight space; the caller removes the weighting
// (fused with its own column scaling of U).
// ---------------------------------------------------------------------------
inline void SVD_triplet(const arma::mat &Xhat,
                        const double *sw,
                        const arma::uword ncp,
                        const bool tall,
                        arma::vec &vs_top,
                        double &trace_val,
                        arma::mat &U,
                        arma::mat &V,
                        arma::mat &X_work,
                        arma::mat &AA_NxN,
                        arma::vec &d_inv,
                        arma::mat &Vn,
                        arma::vec &eigvals,
                        arma::mat &eigvecs,
                        GramWorkspace &gram_ws,
                        GramCache &gram_cache,
                        HybridEigContext &hyb_ctx,
                        int outer_iter
                            LOC_TIMER_PARAM(timer))
{
  LOC_TIC(timer, "form_gram");
  form_weighted_gram(Xhat, sw, tall, X_work, AA_NxN, gram_ws, gram_cache);
  LOC_TOC(timer, "form_gram");

  trace_val = arma::trace(AA_NxN);

  LOC_TIC(timer, "eig");
  if (!hybrid_topk_eig(eigvals, eigvecs, AA_NxN, hyb_ctx, outer_iter))
  {
    Rcpp::stop("top-k symmetric eigensolver failed to converge");
  }
  LOC_TOC(timer, "eig");

  const arma::uword ke = eigvals.n_elem;
  vs_top.set_size(ke);
  d_inv.set_size(ncp);
  {
    const double *ep = eigvals.memptr();
    double *vp = vs_top.memptr();
    double *dp = d_inv.memptr();
    for (arma::uword i = 0; i < ke; ++i)
    {
      const double e = ep[i];
      vp[i] = (e > 0.0) ? std::sqrt(e) : 0.0;
    }
    for (arma::uword i = 0; i < ncp; ++i)
    {
      dp[i] = (vp[i] > PCA_TOL) ? (1.0 / vp[i]) : 0.0;
    }
  }

  LOC_TIC(timer, "recover_factors");
  if (tall)
  {
    split_eigvec_block(eigvecs, V, Vn, d_inv, nullptr, ncp);
    pca_detail::gemm_nn(X_work, Vn, U);
  }
  else
  {
    split_eigvec_block(eigvecs, U, Vn, d_inv, sw, ncp);
    pca_detail::gemm_tn(Xhat, Vn, V);
  }
  LOC_TOC(timer, "recover_factors");
}

#endif // SVD_TRIPLET_H
