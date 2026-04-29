#ifndef SVD_TRIPLET_H
#define SVD_TRIPLET_H

#include <RcppArmadillo.h>
#include <chrono>
#include "gram_ops.h"
#include "eig_sym_sel.h"
#include "hybrid_topk_eig.h"
#include "pca_linalg_utils.h"
#include "loc_timer.h"

// ---------------------------------------------------------------------------
// SVD_triplet
// ---------------------------------------------------------------------------
inline void SVD_triplet(const arma::mat &Xhat,
                        const double *sw,
                        const double *isw,
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
                        EigSymWorkspace &eig_ws,
                        GramCache &gram_cache,
                        HybridEigContext *hyb_ctx,
                        int outer_iter
                            LOC_TIMER_PARAM(timer))
{
  LOC_TIC(timer, "form_gram");
  form_weighted_gram(Xhat, sw, tall, X_work, AA_NxN, gram_ws, gram_cache);
  LOC_TOC(timer, "form_gram");

  trace_val = arma::trace(AA_NxN);

  LOC_TIC(timer, "eig");

  bool ok = false;

  if (hyb_ctx)
  {
    if (hyb_ctx->auto_timing_active())
    {
      const auto t0 = std::chrono::steady_clock::now();

      ok = hybrid_topk_eig(eigvals, eigvecs, AA_NxN, *hyb_ctx, outer_iter);

      const auto t1 = std::chrono::steady_clock::now();

      const double seconds = std::chrono::duration<double>(t1 - t0).count();

      hyb_ctx->record_auto_probe_time(seconds);
    }
    else
    {
      ok = hybrid_topk_eig(eigvals, eigvecs, AA_NxN, *hyb_ctx, outer_iter);
    }
  }
  else
  {
    ok = eig_sym_sel(eigvals, eigvecs, AA_NxN, eig_ws);
  }

  if (!ok)
  {
    Rcpp::stop("`eig_sym_sel` failed to converge");
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
    const arma::uword n = eigvecs.n_rows;
    V.set_size(n, ncp);
    Vn.set_size(n, ncp);
    const double *EV = eigvecs.memptr();
    double *Vp = V.memptr();
    double *Vnp = Vn.memptr();
    const double *dp = d_inv.memptr();
    for (arma::uword j = 0; j < ncp; ++j)
    {
      const double dj = dp[j];
      const double *ecol = EV + j * n;
      double *vcol = Vp + j * n;
      double *ncol = Vnp + j * n;
      for (arma::uword i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        vcol[i] = e;
        ncol[i] = e * dj;
      }
    }
    pca_detail::gemm_nn(X_work, Vn, U);
  }
  else
  {
    const arma::uword n = eigvecs.n_rows;
    U.set_size(n, ncp);
    Vn.set_size(n, ncp);
    const double *EV = eigvecs.memptr();
    double *Up = U.memptr();
    double *Vnp = Vn.memptr();
    const double *dp = d_inv.memptr();
    for (arma::uword j = 0; j < ncp; ++j)
    {
      const double dj = dp[j];
      const double *ecol = EV + j * n;
      double *ucol = Up + j * n;
      double *ncol = Vnp + j * n;
      for (arma::uword i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        ucol[i] = e;
        ncol[i] = e * dj * sw[i];
      }
    }
    pca_detail::gemm_tn(Xhat, Vn, V);
  }
  pca_detail::scale_rows_inplace(U.memptr(), isw, U.n_rows, U.n_cols);
  LOC_TOC(timer, "recover_factors");
}

#endif // SVD_TRIPLET_H
