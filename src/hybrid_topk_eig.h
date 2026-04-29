#ifndef HYBRID_TOPK_EIG_H
#define HYBRID_TOPK_EIG_H
//
// Hybrid top-k symmetric eigensolver.
//

#include <RcppArmadillo.h>
#include <limits>
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
  HYB_REASON_DISABLED = -12,
  HYB_REASON_AUTO_EXACT_PROBE = -13,
  HYB_REASON_AUTO_CHOSE_EXACT = -14,
  HYB_REASON_FORCED_EXACT = -15
};

enum class HybridSolverMode : int
{
  Exact = 0,
  Lobpcg = 1,
  Auto = 2
};

enum class HybridAutoChoice : int
{
  Undecided = 0,
  Exact = 3,
  Lobpcg = 4
};

enum class HybridAutoProbe : int
{
  None = 0,
  Exact = 1,
  Lobpcg = 2
};

// v1 auto policy.
//
// auto currently uses a fixed number of exact and LOBPCG probe calls.
inline constexpr int HYB_AUTO_N_PROBE_ITER_DEFAULT = 10;
inline constexpr double HYB_AUTO_MARGIN_DEFAULT = 0.10;

struct HybridEigContext
{
  EigSymWorkspace *eig_ws = nullptr;
  LOBPCGState lobpcg_state;
  LOBPCGOptions lobpcg_opt;

  arma::mat X_prev;
  int warmup_iters = 10;

  HybridSolverMode solver_mode = HybridSolverMode::Lobpcg;

  HybridAutoChoice auto_choice = HybridAutoChoice::Undecided;
  HybridAutoProbe last_auto_probe = HybridAutoProbe::None;

  int auto_n_probe_iter = HYB_AUTO_N_PROBE_ITER_DEFAULT;
  double auto_margin = HYB_AUTO_MARGIN_DEFAULT;

  int auto_exact_n = 0;
  int auto_lobpcg_n = 0;
  double auto_exact_sec = 0.0;
  double auto_lobpcg_sec = 0.0;

  int n_exact = 0;
  int n_lobpcg_ok = 0;
  int n_lobpcg_bad = 0;

  std::vector<int> path_log; // 0=exact_direct, 1=lobpcg_ok, 2=exact_fallback
  std::vector<int> reason_log;
  std::vector<int> lobpcg_iter_log;
  std::vector<int> lobpcg_fail_iter_log;
  std::vector<double> lobpcg_max_rel_res_log;

  int auto_exact_target() const
  {
    return warmup_iters + auto_n_probe_iter;
  }

  bool auto_timing_active() const
  {
    return solver_mode == HybridSolverMode::Auto &&
           auto_choice == HybridAutoChoice::Undecided;
  }

  void choose_auto_exact()
  {
    auto_choice = HybridAutoChoice::Exact;
    // once auto chooses exact, discard LOBPCG state so exact mode has no
    // continuing seed-maintenance overhead.
    lobpcg_state.X.reset();
    lobpcg_state.P.reset();
    X_prev.reset();
  }

  void choose_auto_lobpcg()
  {
    auto_choice = HybridAutoChoice::Lobpcg;
  }

  void maybe_decide_auto()
  {
    if (solver_mode != HybridSolverMode::Auto)
    {
      return;
    }

    if (auto_choice != HybridAutoChoice::Undecided)
    {
      return;
    }

    if (auto_lobpcg_n >= auto_n_probe_iter && auto_exact_n > 0)
    {
      const double exact_mean =
          auto_exact_sec / static_cast<double>(auto_exact_n);

      const double lobpcg_mean =
          auto_lobpcg_sec / static_cast<double>(auto_lobpcg_n);

      if (exact_mean > 0.0 &&
          lobpcg_mean < (1.0 - auto_margin) * exact_mean)
      {
        choose_auto_lobpcg();
      }
      else
      {
        choose_auto_exact();
      }

      return;
    }
    // if the exact probe phase is complete but we still cannot build a usable
    // LOBPCG seed, auto cannot fairly probe LOBPCG. Choose exact.
    if (auto_exact_n >= auto_exact_target() && !lobpcg_state.seeded())
    {
      choose_auto_exact();
    }
  }

  void record_auto_probe_time(double seconds)
  {
    if (solver_mode == HybridSolverMode::Auto &&
        auto_choice == HybridAutoChoice::Undecided)
    {
      if (last_auto_probe == HybridAutoProbe::Exact)
      {
        auto_exact_sec += seconds;
        ++auto_exact_n;
      }
      else if (last_auto_probe == HybridAutoProbe::Lobpcg)
      {
        auto_lobpcg_sec += seconds;
        ++auto_lobpcg_n;
      }

      maybe_decide_auto();
    }

    last_auto_probe = HybridAutoProbe::None;
  }

  double auto_exact_mean_sec() const
  {
    return auto_exact_n > 0
               ? auto_exact_sec / static_cast<double>(auto_exact_n)
               : std::numeric_limits<double>::quiet_NaN();
  }

  double auto_lobpcg_mean_sec() const
  {
    return auto_lobpcg_n > 0
               ? auto_lobpcg_sec / static_cast<double>(auto_lobpcg_n)
               : std::numeric_limits<double>::quiet_NaN();
  }

  bool keep_lobpcg_seed_state() const
  {
    if (solver_mode == HybridSolverMode::Exact)
    {
      return false;
    }

    if (lobpcg_opt.maxiter <= 0)
    {
      return false;
    }

    if (solver_mode == HybridSolverMode::Auto &&
        auto_choice == HybridAutoChoice::Exact)
    {
      return false;
    }

    return true;
  }

  int solver_chosen_code() const
  {
    if (solver_mode == HybridSolverMode::Exact)
    {
      return 0;
    }

    if (solver_mode == HybridSolverMode::Lobpcg)
    {
      return 1;
    }

    if (auto_choice == HybridAutoChoice::Exact)
    {
      return 3;
    }

    if (auto_choice == HybridAutoChoice::Lobpcg)
    {
      return 4;
    }
    // solver = auto, but the run stopped before the auto probe completed.
    return 2;
  }
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

  ctx.last_auto_probe = HybridAutoProbe::None;

  const bool past_warmup = outer_iter >= ctx.warmup_iters;
  const bool has_seed = ctx.lobpcg_state.seeded();
  const bool lobpcg_enabled = ctx.lobpcg_opt.maxiter > 0;

  bool try_lobpcg = false;
  int direct_reason = HYB_REASON_WARMUP;

  if (ctx.solver_mode == HybridSolverMode::Exact)
  {
    direct_reason = HYB_REASON_FORCED_EXACT;
  }
  else if (ctx.solver_mode == HybridSolverMode::Lobpcg)
  {
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
    else
    {
      try_lobpcg = true;
    }
  }
  else // HybridSolverMode::Auto
  {
    if (!lobpcg_enabled)
    {
      direct_reason = HYB_REASON_DISABLED;
      ctx.choose_auto_exact();
    }
    else if (ctx.auto_choice == HybridAutoChoice::Exact)
    {
      direct_reason = HYB_REASON_AUTO_CHOSE_EXACT;
    }
    else if (ctx.auto_choice == HybridAutoChoice::Lobpcg)
    {
      if (has_seed)
      {
        try_lobpcg = true;
      }
      else
      {
        direct_reason = HYB_REASON_NO_SEED;
      }
    }
    else // auto undecided
    {
      if (ctx.auto_exact_n < ctx.auto_exact_target())
      {
        direct_reason = HYB_REASON_AUTO_EXACT_PROBE;
        ctx.last_auto_probe = HybridAutoProbe::Exact;
      }
      else if (!has_seed)
      {
        direct_reason = HYB_REASON_NO_SEED;
        ctx.choose_auto_exact();
      }
      else if (ctx.auto_lobpcg_n < ctx.auto_n_probe_iter)
      {
        try_lobpcg = true;
        ctx.last_auto_probe = HybridAutoProbe::Lobpcg;
      }
      else
      {
        ctx.maybe_decide_auto();

        if (ctx.auto_choice == HybridAutoChoice::Lobpcg && has_seed)
        {
          try_lobpcg = true;
        }
        else
        {
          direct_reason = HYB_REASON_AUTO_CHOSE_EXACT;
        }
      }
    }
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

  if (ctx.keep_lobpcg_seed_state())
  {
    if (!ctx.X_prev.is_empty())
    {
      seed_lobpcg_state(ctx.lobpcg_state, eigvecs, ctx.X_prev, ctx.lobpcg_opt);
    }

    ctx.X_prev = eigvecs;
  }

  return true;
}

#endif // HYBRID_TOPK_EIG_H
