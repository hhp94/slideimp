#ifndef HYBRID_TOPK_EIG_H
#define HYBRID_TOPK_EIG_H
//
// Hybrid top-k symmetric eigensolver.
//

#include <RcppArmadillo.h>
#include "eig_sym_sel.h"
#include "lobpcg_warm.h"

inline bool make_temporal_P(arma::mat &P,
                            const arma::mat &X_curr,
                            const arma::mat &X_prev,
                            const LOBPCGOptions &opt)
{
  using namespace lobpcg_detail;

  if (X_prev.n_rows != X_curr.n_rows || X_prev.n_cols != X_curr.n_cols)
  {
    return false;
  }

  P = X_prev - X_curr * (X_curr.t() * X_prev);

  const double p_fro = arma::norm(P, "fro");
  const double min_temporal_norm = 1e-10 * std::sqrt(static_cast<double>(X_curr.n_cols));

  if (!(p_fro > min_temporal_norm))
  {
    P.reset();
    return false;
  }

  if (qr_orthonormalize(P, opt.qr_rel_tol) == 0 ||
      P.n_cols != X_curr.n_cols)
  {
    P.reset();
    return false;
  }

  P -= X_curr * (X_curr.t() * P);

  if (qr_orthonormalize(P, opt.qr_rel_tol) == 0 ||
      P.n_cols != X_curr.n_cols)
  {
    P.reset();
    return false;
  }

  return true;
}

enum HybridEigPath
{
  HYB_NOT_RUN = -1,
  HYB_DSYEVR_DIRECT = 0,
  HYB_LOBPCG_OK = 1,
  HYB_DSYEVR_FALLBACK = 2
};

enum HybridEigReason
{
  HYB_REASON_NOT_RUN = -100,
  HYB_REASON_WARMUP = -10,
  HYB_REASON_NO_SEED = -11,
  HYB_REASON_DISABLED = -12
};

struct HybridEigContext
{
  EigSymWorkspace *eig_ws = nullptr;
  LOBPCGState lobpcg_state;
  LOBPCGOptions lobpcg_opt;

  arma::mat X_prev;
  int warmup_iters = 10; // controlled by hyb_ctx context

  int n_exact = 0;
  int n_lobpcg_ok = 0;
  int n_lobpcg_bad = 0;
  std::vector<int> path_log; // 0=exact_direct, 1=lobpcg_ok, 2=exact_fallback
  std::vector<int> reason_log;
  std::vector<int> lobpcg_iter_log;
  std::vector<int> lobpcg_fail_iter_log;
  std::vector<double> lobpcg_max_rel_res_log;
};

inline bool hybrid_topk_eig(arma::vec &eigvals,
                            arma::mat &eigvecs,
                            arma::mat &A,
                            HybridEigContext &ctx,
                            int outer_iter)
{
  auto log_at = [&](int path, int reason, const LOBPCGResult *res = nullptr)
  {
    if (outer_iter < 0)
    {
      return;
    }

    const size_t oi = static_cast<size_t>(outer_iter);

    if (oi < ctx.path_log.size())
    {
      ctx.path_log[oi] = path;
    }
    if (oi < ctx.reason_log.size())
    {
      ctx.reason_log[oi] = reason;
    }
    if (oi < ctx.lobpcg_iter_log.size())
    {
      ctx.lobpcg_iter_log[oi] = res ? res->iterations : -1;
    }
    if (oi < ctx.lobpcg_fail_iter_log.size())
    {
      ctx.lobpcg_fail_iter_log[oi] = res ? res->fail_iter : -1;
    }
    if (oi < ctx.lobpcg_max_rel_res_log.size())
    {
      ctx.lobpcg_max_rel_res_log[oi] =
          res ? res->max_rel_res : std::numeric_limits<double>::quiet_NaN();
    }
  };

  const bool past_warmup = outer_iter >= ctx.warmup_iters;
  const bool has_seed = ctx.lobpcg_state.seeded();
  const bool lobpcg_enabled = ctx.lobpcg_opt.maxiter > 0;
  const bool try_lobpcg = past_warmup && has_seed && lobpcg_enabled;

  int direct_reason = HYB_REASON_WARMUP;
  if (!lobpcg_enabled)
  {
    direct_reason = HYB_REASON_DISABLED;
  }
  else if (!past_warmup)
  {
    direct_reason = HYB_REASON_WARMUP;
  }
  else if (!has_seed)
  {
    direct_reason = HYB_REASON_NO_SEED;
  }

  LOBPCGResult lob_res;
  bool attempted_lobpcg = false;

  if (try_lobpcg)
  {
    attempted_lobpcg = true;

    pca_detail::mirror_upper_to_lower(A);

    lob_res = lobpcg_solve(A, ctx.lobpcg_state, ctx.lobpcg_opt);

    if (lob_res.converged)
    {
      ctx.n_lobpcg_ok++;

      eigvals = std::move(lob_res.eigvals);
      eigvecs = ctx.lobpcg_state.X;

      if (!ctx.X_prev.is_empty())
      {
        arma::mat P_temporal;
        if (make_temporal_P(P_temporal, eigvecs, ctx.X_prev, ctx.lobpcg_opt))
        {
          ctx.lobpcg_state.P = std::move(P_temporal);
        }
      }

      ctx.X_prev = eigvecs;

      log_at(HYB_LOBPCG_OK, lob_res.status, &lob_res);
      return true;
    }

    ++ctx.n_lobpcg_bad;
  }

  if (!eig_sym_sel(eigvals, eigvecs, A, *ctx.eig_ws))
  {
    return false;
  }

  ++ctx.n_exact;

  if (attempted_lobpcg)
  {
    log_at(HYB_DSYEVR_FALLBACK, lob_res.status, &lob_res);
  }
  else
  {
    log_at(HYB_DSYEVR_DIRECT, direct_reason, nullptr);
  }

  if (!ctx.X_prev.is_empty())
  {
    seed_lobpcg_state(ctx.lobpcg_state, eigvecs, ctx.X_prev, ctx.lobpcg_opt);
  }

  ctx.X_prev = eigvecs;
  return true;
}

#endif // HYBRID_TOPK_EIG_H
