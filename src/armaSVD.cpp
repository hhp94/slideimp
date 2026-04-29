#ifndef PCA_IMP_DIAGNOSTICS
#define PCA_IMP_DIAGNOSTICS 0
#endif

#include "svd_triplet.h"
#include "imputed_value.h"
#include "loc_timer.h"
#include <limits>

template <typename F>
static inline void for_each_missing_in_col(arma::uword cidx,
                                           const arma::uvec &offsets,
                                           const arma::uvec &rows_flat,
                                           F &&f)
{
  const arma::uword beg = offsets[cidx];
  const arma::uword end = offsets[cidx + 1];
  for (arma::uword p = beg; p < end; ++p)
  {
    f(rows_flat[p]);
  }
}

// ---------------------------------------------------------------------------
// impute_restandardize
// ---------------------------------------------------------------------------
template <bool Scale>
void impute_restandardize(arma::mat &Xhat,
                          const arma::mat &fittedX,
                          arma::rowvec &mean_p,
                          arma::rowvec &et,
                          const double *wptr,
                          const double *swptr,
                          double *X_chg_scaled,
                          const arma::uword nrX,
                          const arma::uvec &miss_cols,
                          const arma::uvec &miss_rows_flat,
                          const arma::uvec &miss_rows_offsets,
                          const double *obs_sum,
                          const double *obs_sumsq)
{
  double *xptr = Xhat.memptr();
  const double *fxp = fittedX.memptr();

  const arma::uword n_mc = miss_cols.n_elem;

  for (arma::uword cidx = 0; cidx < n_mc; ++cidx)
  {
    const arma::uword j = miss_cols[cidx];
    const arma::uword off = j * nrX;

    double *col = xptr + off;
    const double *fit_col = fxp + off;

    // Compact changing-column weighted buffer.
    //
    // Tall:
    //   points into X_work.colptr(n_fixed)
    //
    // Wide cached:
    //   points into gram_cache.X_chg_scaled
    //
    // Wide non-cached:
    //   nullptr
    double *scaled_col = X_chg_scaled ? (X_chg_scaled + cidx * nrX) : nullptr;

    const double old_mu = mean_p(j);

    double old_et;
    if constexpr (Scale)
    {
      old_et = et(j);
    }
    else
    {
      old_et = 1.0;
    }

    double imp_sum = 0.0;
    double imp_sumsq = 0.0;

    for_each_missing_in_col(cidx, miss_rows_offsets, miss_rows_flat,
                            [&](arma::uword i)
                            {
                              const double v_std = fit_col[i];
                              col[i] = v_std;

                              const double w = wptr[i];

                              if constexpr (Scale)
                              {
                                const double v = v_std * old_et + old_mu;
                                imp_sum += w * v;
                                imp_sumsq += w * v * v;
                              }
                              else
                              {
                                const double v = v_std + old_mu;
                                imp_sum += w * v;
                              }
                            });

    const double new_mu = obs_sum[j] + imp_sum;

    if constexpr (Scale)
    {
      const double new_var = obs_sumsq[j] + imp_sumsq - new_mu * new_mu;
      const double new_et = std::sqrt(std::max(0.0, new_var));
      const double new_et_eff = std::max(new_et, PCA_TOL);

      const double inv_new = 1.0 / new_et_eff;
      const double a = old_et * inv_new;
      const double b = (old_mu - new_mu) * inv_new;

      if (scaled_col)
      {
        for (arma::uword i = 0; i < nrX; ++i)
        {
          const double z = col[i] * a + b;
          col[i] = z;
          scaled_col[i] = z * swptr[i];
        }
      }
      else
      {
        for (arma::uword i = 0; i < nrX; ++i)
        {
          col[i] = col[i] * a + b;
        }
      }

      et(j) = new_et_eff;
    }
    else
    {
      const double b = old_mu - new_mu;

      if (scaled_col)
      {
        for (arma::uword i = 0; i < nrX; ++i)
        {
          const double z = col[i] + b;
          col[i] = z;
          scaled_col[i] = z * swptr[i];
        }
      }
      else
      {
        for (arma::uword i = 0; i < nrX; ++i)
        {
          col[i] += b;
        }
      }
    }

    mean_p(j) = new_mu;
  }
}

static inline void impute_restandardize_dispatch(
    const bool scale,
    arma::mat &Xhat,
    const arma::mat &fittedX,
    arma::rowvec &mean_p,
    arma::rowvec &et,
    const double *wptr,
    const double *swptr,
    double *X_chg_scaled,
    const arma::uword nrX,
    const arma::uvec &miss_cols,
    const arma::uvec &miss_rows_flat,
    const arma::uvec &miss_rows_offsets,
    const double *obs_sum,
    const double *obs_sumsq)
{
  if (scale)
  {
    impute_restandardize<true>(Xhat, fittedX, mean_p, et,
                               wptr, swptr, X_chg_scaled, nrX,
                               miss_cols, miss_rows_flat, miss_rows_offsets,
                               obs_sum, obs_sumsq);
  }
  else
  {
    impute_restandardize<false>(Xhat, fittedX, mean_p, et,
                                wptr, swptr, X_chg_scaled, nrX,
                                miss_cols, miss_rows_flat, miss_rows_offsets,
                                obs_sum, obs_sumsq);
  }
}

static inline void initialize_missing_start(arma::mat &Xhat,
                                            arma::mat &fittedX,
                                            const arma::uword init,
                                            const arma::uword nrX,
                                            const arma::uvec &miss_cols,
                                            const arma::uvec &miss_rows_flat,
                                            const arma::uvec &miss_rows_offsets)
{
  // Allocate the full fitted matrix once. It will be overwritten by DGEMM
  // in the hot loop. Before the first DGEMM, only missing positions matter.
  fittedX.zeros(nrX, Xhat.n_cols);

  const arma::uword n_mc = miss_cols.n_elem;

  if (n_mc == 0)
  {
    return;
  }

  if (init == 0)
  {
    for (arma::uword ci = 0; ci < n_mc; ++ci)
    {
      const arma::uword j = miss_cols[ci];

      double *col = Xhat.colptr(j);
      double *fit_col = fittedX.colptr(j);

      for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                              [&](arma::uword i)
                              {
                                col[i] = 0.0;
                                fit_col[i] = 0.0;
                              });
    }

    return;
  }

  const arma::uword n_missing = miss_rows_flat.n_elem;

  if (n_missing == 0)
  {
    return;
  }

  Rcpp::NumericVector rand_vals = Rcpp::rnorm(n_missing);

  for (arma::uword ci = 0; ci < n_mc; ++ci)
  {
    const arma::uword j = miss_cols[ci];

    double *col = Xhat.colptr(j);
    double *fit_col = fittedX.colptr(j);

    arma::uword p = miss_rows_offsets[ci];

    for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                            [&](arma::uword i)
                            {
                              const double z = rand_vals[p++];
                              col[i] = z;
                              fit_col[i] = z;
                            });
  }
}

static inline double missing_weighted_residual_ss(
    const arma::mat &Xhat,
    const arma::mat &fittedX,
    const double *wptr,
    const arma::uword nrX,
    const arma::uvec &miss_cols,
    const arma::uvec &miss_rows_flat,
    const arma::uvec &miss_rows_offsets)
{
  const double *xhp = Xhat.memptr();
  const double *fxp = fittedX.memptr();

  double out = 0.0;

  const arma::uword n_mc = miss_cols.n_elem;

  for (arma::uword ci = 0; ci < n_mc; ++ci)
  {
    const arma::uword j = miss_cols[ci];
    const arma::uword off = j * nrX;

    const double *xh = xhp + off;
    const double *fx = fxp + off;

    double s = 0.0;

    for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                            [&](arma::uword i)
                            {
                              const double d = xh[i] - fx[i];
                              s += wptr[i] * d * d;
                            });

    out += s;
  }

  return out;
}

// [[Rcpp::export]]
Rcpp::List pca_imp_internal_cpp(
    const arma::mat &obj,
    const arma::uvec &eligible_idx,
    const arma::uword ncp,
    const bool scale,
    const bool regularized,
    const double threshold,
    const arma::uword init,
    const arma::uword maxiter,
    const arma::uword miniter,
    arma::rowvec row_w,
    const double coeff_ridge,
    const int solver,
    const arma::uword warmup_iters,
    const double lobpcg_tol,
    const arma::uword lobpcg_maxiter)
{
  LOC_TIMER_OBJ(pca_imp_gram);
  LOC_TIC(pca_imp_gram, "pca_imp_internal_cpp_total");
  stop_on_inf(obj);

  // -------------------------------------------------------------------------
  // basic argument validation before we set up workspaces.
  // -------------------------------------------------------------------------
  const arma::uword nrX = obj.n_rows;
  const arma::uword n_elig = eligible_idx.n_elem;

  if (maxiter == 0)
  {
    Rcpp::stop("maxiter must be >= 1");
  }
  if (miniter > maxiter)
  {
    Rcpp::stop("miniter must be <= maxiter");
  }
  if (ncp == 0)
  {
    Rcpp::stop("ncp must be >= 1");
  }
  if (solver < 0 || solver > 2)
  {
    Rcpp::stop("solver must be 0 (exact), 1 (lobpcg), or 2 (auto)");
  }

  const bool solver_auto = (solver == 2);
  const bool lobpcg_requested = (solver != 0);

  if (lobpcg_requested && lobpcg_maxiter == 0)
  {
    Rcpp::stop("lobpcg_maxiter must be >= 1 when solver is 1 (lobpcg) or 2 (auto)");
  }

  if (n_elig == 0)
  {
    Rcpp::stop("eligible_idx is empty");
  }
  if (row_w.n_elem != nrX)
  {
    Rcpp::stop("row_w length must equal n_rows(obj)");
  }
  if (eligible_idx.max() >= obj.n_cols)
  {
    Rcpp::stop("eligible_idx contains out-of-bounds column indices");
  }

  const arma::uword min_problem_dim = std::min(nrX, n_elig);
  if (ncp >= min_problem_dim)
  {
    Rcpp::stop("ncp + 1 must be <= min(n_rows, n_eligible_cols)");
  }

  // prevent nrX * n_elig overflow later.
  if (n_elig != 0 && nrX > (std::numeric_limits<arma::uword>::max)() / n_elig)
  {
    Rcpp::stop("problem dimensions are too large");
  }

  // BLAS/LAPACK integer range guard
  {
    const arma::uword blas_max =
        static_cast<arma::uword>((std::numeric_limits<arma::blas_int>::max)());

    const arma::uword k_eig = regularized ? (ncp + 1) : ncp;
    const arma::uword n_gram_check = (nrX >= n_elig) ? n_elig : nrX;

    if (nrX > blas_max || n_elig > blas_max || ncp > blas_max || k_eig > blas_max || n_gram_check > blas_max)
    {
      Rcpp::stop("problem dimensions exceed BLAS/LAPACK integer range");
    }

    // EigSymWorkspace uses workspace minima 26*n and 10*n.
    if (n_gram_check > blas_max / 26)
    {
      Rcpp::stop("problem dimension too large for LAPACK dsyevr workspace");
    }

    // LOBPCG internally solves small Rayleigh-Ritz problems up to ~3*k_eig.
    // only apply this guard when the requested solver can actually use LOBPCG.
    if (lobpcg_requested && lobpcg_maxiter > 0)
    {
      const long double rr_dim = 3.0L * static_cast<long double>(k_eig);
      const long double dsyevd_lwork = 1.0L + 6.0L * rr_dim + 2.0L * rr_dim * rr_dim;
      const long double dsyevd_liwork = 3.0L + 5.0L * rr_dim;

      if (rr_dim > static_cast<long double>(blas_max) ||
          dsyevd_lwork > static_cast<long double>(blas_max) ||
          dsyevd_liwork > static_cast<long double>(blas_max))
      {
        Rcpp::stop("ncp too large for LOBPCG internal LAPACK workspace");
      }
    }

    const arma::uword int_max = static_cast<arma::uword>((std::numeric_limits<int>::max)());

    if (lobpcg_requested)
    {
      if (warmup_iters > int_max || lobpcg_maxiter > int_max)
      {
        Rcpp::stop("LOBPCG iteration controls exceed int range");
      }

      if (solver_auto)
      {
        const arma::uword auto_probe =
            static_cast<arma::uword>(HYB_AUTO_N_PROBE_ITER_DEFAULT);

        if (warmup_iters > int_max - auto_probe)
        {
          Rcpp::stop("warmup_iters + auto probe iterations exceed int range");
        }
      }
    }
  }

  if (!row_w.is_finite())
  {
    Rcpp::stop("row_w must be finite");
  }
  if (arma::any(row_w < 0.0))
  {
    Rcpp::stop("row_w must be non-negative");
  }

  const double row_w_sum = arma::accu(row_w);
  if (!(row_w_sum > 0.0))
  {
    Rcpp::stop("row_w must have positive sum");
  }

  const double nrX_d = static_cast<double>(nrX);
  const double n_elig_d = static_cast<double>(n_elig);
  const double ncp_d = static_cast<double>(ncp);

  const double denom_sigma = (n_elig_d - ncp_d) * (nrX_d - 1.0 - ncp_d);
  if (regularized && !(denom_sigma > 0.0))
  {
    Rcpp::stop("regularized=true requires n_elig > ncp and nrX - 1 > ncp");
  }

  double old = arma::datum::inf;

  row_w /= row_w_sum;
  const double *wptr = row_w.memptr();

  LOC_TIC(pca_imp_gram, "pass1_count_missing");
  arma::uvec col_nmiss(n_elig, arma::fill::zeros);
  for (arma::uword j_local = 0; j_local < n_elig; ++j_local)
  {
    const double *p = obj.colptr(eligible_idx[j_local]);
    arma::uword c = 0;
    for (arma::uword i = 0; i < nrX; ++i)
    {
      c += static_cast<arma::uword>(std::isnan(p[i]));
    }
    col_nmiss[j_local] = c;
  }
  LOC_TOC(pca_imp_gram, "pass1_count_missing");

  arma::uvec fixed_cols = arma::find(col_nmiss == 0);
  arma::uvec changing_cols = arma::find(col_nmiss > 0);
  const arma::uword n_fixed = fixed_cols.n_elem;
  const arma::uword n_mc = changing_cols.n_elem;

  arma::uvec perm(n_elig);
  if (n_fixed > 0)
  {
    perm.head(n_fixed) = fixed_cols;
  }
  if (n_mc > 0)
  {
    perm.tail(n_mc) = changing_cols;
  }

  arma::uvec eligible_idx_perm = eligible_idx.elem(perm);

  // only build a non-empty miss_cols_idx when we actually
  // have missing columns (the downstream loops gate on n_mc > 0 anyway)
  arma::uvec miss_cols_idx;
  if (n_mc > 0)
  {
    miss_cols_idx = arma::regspace<arma::uvec>(n_fixed, n_elig - 1);
  }

  arma::uvec miss_rows_offsets(n_mc + 1, arma::fill::zeros);
  if (n_mc > 0)
  {
    miss_rows_offsets.tail(n_mc) = arma::cumsum(col_nmiss.elem(changing_cols));
  }
  const arma::uword n_missing = miss_rows_offsets[n_mc];
  arma::uvec miss_rows_flat(n_missing);

  arma::mat Xhat(nrX, n_elig);
  arma::mat fittedX;

  arma::rowvec denom(n_elig, arma::fill::zeros);
  arma::rowvec weighted_sum(n_elig, arma::fill::zeros);
  arma::rowvec sum_sq(n_elig, arma::fill::zeros);

  LOC_TIC(pca_imp_gram, "pass2_scan");
  {
    for (arma::uword jp = 0; jp < n_elig; ++jp)
    {
      const arma::uword orig_col = eligible_idx_perm[jp];
      const double *src = obj.colptr(orig_col);
      double *dst = Xhat.colptr(jp);

      const bool has_miss = (jp >= n_fixed);
      const arma::uword cidx = has_miss ? (jp - n_fixed) : 0;
      arma::uword p = has_miss ? miss_rows_offsets[cidx] : 0;

      double d_acc = 0.0, ws_acc = 0.0, ss_acc = 0.0;
      for (arma::uword i = 0; i < nrX; ++i)
      {
        const double v = src[i];
        if (std::isnan(v))
        {
          dst[i] = 0.0;
          if (has_miss)
          {
            miss_rows_flat[p++] = i;
          }
        }
        else
        {
          dst[i] = v;
          const double w = wptr[i];
          d_acc += w;
          ws_acc += w * v;
          ss_acc += w * v * v;
        }
      }
      denom[jp] = d_acc;
      weighted_sum[jp] = ws_acc;
      sum_sq[jp] = ss_acc;
    }
  }
  LOC_TOC(pca_imp_gram, "pass2_scan");

  // every eligible column must have at least one observed row with positive
  // weight, else mean_p (and the downstream Gram) is NaN.
  for (arma::uword j = 0; j < n_elig; ++j)
  {
    if (!(denom[j] > 0.0))
    {
      Rcpp::stop("eligible column has zero observed weight (column %u)",
                 static_cast<unsigned>(eligible_idx_perm[j] + 1));
    }
  }

  const double *wsptr = weighted_sum.memptr();
  const double *ssptr = sum_sq.memptr();
  arma::rowvec mean_p = weighted_sum / denom;
  arma::rowvec et;
  if (scale)
  {
    // clamp tiny-negative variance from rounding before sqrt, and floor the
    // scale at PCA_TOL so 1/et is finite. Storing the floored value here is
    // consistent with impute_restandardize<true> committing
    // et(j) = max(new_et, PCA_TOL).
    arma::rowvec var = sum_sq / denom - mean_p % mean_p;
    var.transform([](double v)
                  { return v < 0.0 ? 0.0 : v; });
    et = arma::sqrt(var);
    et.transform([](double e)
                 { return e < PCA_TOL ? PCA_TOL : e; });
  }

  if (scale)
  {
    for (arma::uword j = 0; j < n_elig; ++j)
    {
      double *col = Xhat.colptr(j);
      const double mu = mean_p[j];
      const double inv_et = 1.0 / et[j];

      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] = (col[i] - mu) * inv_et;
      }
    }
  }
  else
  {
    for (arma::uword j = 0; j < n_elig; ++j)
    {
      double *col = Xhat.colptr(j);
      const double mu = mean_p[j];

      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] -= mu;
      }
    }
  }

  initialize_missing_start(Xhat, fittedX, init, nrX,
                           miss_cols_idx, miss_rows_flat, miss_rows_offsets);

  // ---------------------------------------------------------------------------
  // SVD / solve buffers.
  // ---------------------------------------------------------------------------
  double trace_val;
  double sigma2 = 0.0;
  const arma::colvec row_w_col = row_w.t();
  const arma::colvec sqrt_row_w = arma::sqrt(row_w_col);
  arma::colvec inv_sqrt_row_w(nrX);
  for (arma::uword i = 0; i < nrX; ++i)
  {
    const double s = sqrt_row_w[i];
    inv_sqrt_row_w[i] = (s > PCA_TOL) ? (1.0 / s) : (1.0 / PCA_TOL);
  }
  const double *sptr = sqrt_row_w.memptr();
  const double *isptr = inv_sqrt_row_w.memptr();

  const bool tall = (nrX >= n_elig);
  const arma::uword k = regularized ? ncp + 1 : ncp;
  const arma::uword n_gram = tall ? n_elig : nrX;

  arma::mat U, V, X_work, AA_NxN, Vn, eigvecs;
  arma::vec vs_top, d_inv, eigvals, lambda_shrinked;

  GramWorkspace gram_ws;
  if (!gram_ws.init(nrX, n_elig, tall))
  {
    Rcpp::stop("gram workspace init failed");
  }
  EigSymWorkspace eig_ws;
  if (!eig_ws.init(n_gram, k))
  {
    Rcpp::stop("eig workspace query failed");
  }

  AA_NxN.set_size(n_gram, n_gram);
  U.set_size(nrX, ncp);
  V.set_size(n_elig, ncp);
  Vn.set_size(tall ? n_elig : nrX, ncp);
  if (tall)
  {
    X_work.set_size(nrX, n_elig);
  }
  eigvecs.set_size(n_gram, k);
  eigvals.set_size(k);
  vs_top.set_size(k);
  d_inv.set_size(ncp);
  lambda_shrinked.set_size(ncp);

  LOC_TIC(pca_imp_gram, "gram_cache_init");
  GramCache gram_cache;
  if (!gram_cache_init(gram_cache, Xhat, sptr, nrX, n_fixed, n_mc, tall, X_work))
  {
    Rcpp::stop("gram_cache_init failed");
  }
  LOC_TOC(pca_imp_gram, "gram_cache_init");
  double *weighted_chg_ptr = nullptr;

  if (n_mc > 0)
  {
    if (tall)
    {
      // X_work is full nrX x n_elig in tall mode.
      // fixed columns were initialized in gram_cache_init().
      // changing columns will be maintained by impute_restandardize().
      weighted_chg_ptr = X_work.colptr(n_fixed);
    }
    else if (gram_cache.active)
    {
      // wide cached mode uses a compact nrX x n_mc changing block.
      weighted_chg_ptr = gram_cache.X_chg_scaled.memptr();
    }
    // wide non-cached n_fixed == 0 case intentionally remains nullptr.
  }

  HybridEigContext hyb_ctx;
  hyb_ctx.eig_ws = &eig_ws;
  hyb_ctx.solver_mode = static_cast<HybridSolverMode>(solver);
  hyb_ctx.warmup_iters = lobpcg_requested ? static_cast<int>(warmup_iters) : 0;
  hyb_ctx.lobpcg_opt.tol = lobpcg_tol;
  hyb_ctx.lobpcg_opt.maxiter = lobpcg_requested ? static_cast<int>(lobpcg_maxiter) : 0;
  hyb_ctx.auto_n_probe_iter = HYB_AUTO_N_PROBE_ITER_DEFAULT;
  hyb_ctx.auto_margin = HYB_AUTO_MARGIN_DEFAULT;
#if PCA_IMP_DIAGNOSTICS
  hyb_ctx.path_log.assign(maxiter, HYB_NOT_RUN);
  hyb_ctx.reason_log.assign(maxiter, HYB_REASON_NOT_RUN);
  hyb_ctx.lobpcg_iter_log.assign(maxiter, -1);
  hyb_ctx.lobpcg_fail_iter_log.assign(maxiter, -1);
  hyb_ctx.lobpcg_max_rel_res_log.assign(
      maxiter, std::numeric_limits<double>::quiet_NaN());
#endif

  const double min_dim = std::min(n_elig_d, nrX_d - 1.0);
  const double scale_factor = (nrX_d * n_elig_d) / min_dim;

  // ---------------------------------------------------------------------------
  // hot loop
  // ---------------------------------------------------------------------------
  double objective = 0.0;
  const double *xhp = Xhat.memptr();

#if PCA_IMP_DIAGNOSTICS
  arma::mat eigval_hist(ncp, maxiter, arma::fill::zeros);
  arma::vec obj_hist(maxiter, arma::fill::zeros);
  arma::vec subspace_cos_hist(maxiter, arma::fill::zeros);
  arma::mat V_prev;
  arma::uword iters_done = 0;
#endif

  for (arma::uword nb_iter = 1;; ++nb_iter)
  {
    if (nb_iter % 5 == 0)
    {
      Rcpp::checkUserInterrupt();
    }
    LOC_TIC(pca_imp_gram, "restandardize");
    impute_restandardize_dispatch(
        scale,
        Xhat,
        fittedX,
        mean_p,
        et,
        wptr,
        sptr,
        weighted_chg_ptr,
        nrX,
        miss_cols_idx,
        miss_rows_flat,
        miss_rows_offsets,
        wsptr,
        ssptr);
    LOC_TOC(pca_imp_gram, "restandardize");

    LOC_TIC(pca_imp_gram, "svd");
    SVD_triplet(Xhat, sptr, isptr, ncp, tall,
                vs_top, trace_val, U, V, X_work, AA_NxN,
                d_inv, Vn, eigvals, eigvecs,
                gram_ws, eig_ws, gram_cache,
                &hyb_ctx, static_cast<int>(nb_iter - 1) LOC_TIMER_ARG(pca_imp_gram));
    LOC_TOC(pca_imp_gram, "svd");

    LOC_TIC(pca_imp_gram, "post_svd");
    sigma2 = 0.0;
    const double sum_top_sq = arma::dot(vs_top.head(ncp), vs_top.head(ncp));
    // clamp tail energy once and use consistently in both sigma2 (was unclamped)
    // and the objective.
    const double tail = std::max(0.0, trace_val - sum_top_sq);
    if (regularized)
    {
      sigma2 = scale_factor * tail / denom_sigma;
      sigma2 = std::min(sigma2 * coeff_ridge, vs_top(ncp) * vs_top(ncp));
    }

    double top_corr = 0.0;
    {
      const double *vp = vs_top.memptr();
      double *lp = lambda_shrinked.memptr();
      for (arma::uword i = 0; i < ncp; ++i)
      {
        const double v = vp[i];
        const double lam = (v > PCA_TOL) ? (v - sigma2 / v) : 0.0;
        lp[i] = lam;
        const double d = v - lam;
        top_corr += d * d;
      }
    }
    const double objective_full = tail + top_corr;
    LOC_TOC(pca_imp_gram, "post_svd");

    LOC_TIC(pca_imp_gram, "reconstruct");

    // after this, U contains U %*% diag(lambda_shrinked).
    pca_detail::scale_cols_inplace(U.memptr(),
                                   lambda_shrinked.memptr(),
                                   U.n_rows,
                                   ncp);
    // we are reconstructing the full fittedX from U and V
    pca_detail::gemm_nt(U, V, fittedX);

    LOC_TOC(pca_imp_gram, "reconstruct");

    LOC_TIC(pca_imp_gram, "objective");
    {
      const double miss_contrib =
          missing_weighted_residual_ss(Xhat,
                                       fittedX,
                                       wptr,
                                       nrX,
                                       miss_cols_idx,
                                       miss_rows_flat,
                                       miss_rows_offsets);

      objective = objective_full - miss_contrib;
    }
    LOC_TOC(pca_imp_gram, "objective");

#if PCA_IMP_DIAGNOSTICS
    eigval_hist.col(nb_iter - 1) = vs_top.head(ncp);
    obj_hist(nb_iter - 1) = objective;
    {
      arma::mat V_new = eigvecs.head_cols(ncp);
      if (V_prev.n_elem > 0)
      {
        arma::vec s = arma::svd(V_prev.t() * V_new);
        subspace_cos_hist(nb_iter - 1) = s.min();
      }
      else
      {
        subspace_cos_hist(nb_iter - 1) = arma::datum::nan;
      }
      V_prev = std::move(V_new);
    }
    iters_done = nb_iter;
#endif

    // missMDA convergence test:
    //
    // criterion <- abs(1 - objective / old)
    // old <- objective
    //
    // missMDA increments nb.iter before testing `nb.iter > 5`.
    // Since `nb_iter` here is the number of completed iterations,
    // `nb_iter >= miniter` with default miniter = 5 reproduces that behavior.
    const double criterion = std::abs(1.0 - objective / old);
    old = objective;

    if (!std::isnan(criterion))
    {
      if (criterion < threshold && nb_iter >= miniter)
      {
        break;
      }
      if (objective < threshold && nb_iter >= miniter)
      {
        break;
      }
    }

    if (nb_iter >= maxiter)
    {
      Rcpp::warning("Stopped after " + std::to_string(maxiter) + " iterations");
      break;
    }
  }

  // ---------------------------------------------------------------------------
  // postprocessing.
  // ---------------------------------------------------------------------------
  LOC_TIC(pca_imp_gram, "postprocessing");
  // fittedX already contains the final full reconstruction from the last
  // hot-loop DGEMM. U was scaled by lambda_shrinked in that reconstruction step.
  const double *fxp = fittedX.memptr();

  double sse = 0.0;
  for (arma::uword j = 0; j < n_elig; ++j)
  {
    const double ej = scale ? et[j] : 1.0;
    const double ej2 = ej * ej;
    const arma::uword off = j * nrX;
    const double *xh = xhp + off;
    const double *fx = fxp + off;
    double s = 0.0;
    for (arma::uword i = 0; i < nrX; ++i)
    {
      const double d = xh[i] - fx[i];
      s += d * d;
    }
    sse += s * ej2;
  }
  for (arma::uword ci = 0; ci < n_mc; ++ci)
  {
    const arma::uword j_local = miss_cols_idx[ci];
    const double ej = scale ? et[j_local] : 1.0;
    const double ej2 = ej * ej;
    const arma::uword off = j_local * nrX;
    const double *xh = xhp + off;
    const double *fx = fxp + off;
    double s = 0.0;
    for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                            [&](arma::uword i)
                            {
                              const double d = xh[i] - fx[i];
                              s += d * d;
                            });
    sse -= s * ej2;
  }
  const arma::uword n_obs = nrX * n_elig - n_missing;
  const double mse = (n_obs > 0) ? sse / static_cast<double>(n_obs) : arma::datum::nan;

  arma::mat imputed_values(n_missing, 3);
  {
    arma::uword iv = 0;
    for (arma::uword ci = 0; ci < n_mc; ++ci)
    {
      const arma::uword j_local = miss_cols_idx[ci];
      const double ej = scale ? et[j_local] : 1.0;
      const double mj = mean_p[j_local];
      const double orig_col_1based = static_cast<double>(eligible_idx_perm[j_local] + 1);

      // missMDA::imputePCA() returns completeObs using the final Xhat at
      // missing positions, not the final fittedX. fittedX is still used above
      // for the nb.init selection MSE.
      const double *xh = Xhat.colptr(j_local);

      for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                              [&](arma::uword i)
                              {
                                imputed_values(iv, 0) = static_cast<double>(i + 1);
                                imputed_values(iv, 1) = orig_col_1based;
                                imputed_values(iv, 2) = xh[i] * ej + mj;
                                ++iv;
                              });
    }
  }

  LOC_TOC(pca_imp_gram, "postprocessing");
  LOC_TOC(pca_imp_gram, "pca_imp_internal_cpp_total");
#if PCA_IMP_DIAGNOSTICS
  eigval_hist.resize(ncp, iters_done);
  obj_hist.resize(iters_done);
  subspace_cos_hist.resize(iters_done);

  hyb_ctx.path_log.resize(iters_done);
  hyb_ctx.reason_log.resize(iters_done);
  hyb_ctx.lobpcg_iter_log.resize(iters_done);
  hyb_ctx.lobpcg_fail_iter_log.resize(iters_done);
  hyb_ctx.lobpcg_max_rel_res_log.resize(iters_done);

  return Rcpp::List::create(
      Rcpp::Named("imputed_values") = imputed_values,
      Rcpp::Named("mse") = mse,
      Rcpp::Named("solver_chosen") = hyb_ctx.solver_chosen_code(),
      Rcpp::Named("eigval_hist") = eigval_hist,
      Rcpp::Named("obj_hist") = obj_hist,
      Rcpp::Named("subspace_cos_hist") = subspace_cos_hist,
      Rcpp::Named("n_exact") = hyb_ctx.n_exact,
      Rcpp::Named("n_lobpcg_ok") = hyb_ctx.n_lobpcg_ok,
      Rcpp::Named("n_lobpcg_bad") = hyb_ctx.n_lobpcg_bad,
      Rcpp::Named("auto_n_probe_iter") = hyb_ctx.auto_n_probe_iter,
      Rcpp::Named("auto_exact_n") = hyb_ctx.auto_exact_n,
      Rcpp::Named("auto_lobpcg_n") = hyb_ctx.auto_lobpcg_n,
      Rcpp::Named("auto_exact_mean_sec") = hyb_ctx.auto_exact_mean_sec(),
      Rcpp::Named("auto_lobpcg_mean_sec") = hyb_ctx.auto_lobpcg_mean_sec(),
      Rcpp::Named("eig_path") = hyb_ctx.path_log,
      Rcpp::Named("eig_reason") = hyb_ctx.reason_log,
      Rcpp::Named("eig_lobpcg_iter") = hyb_ctx.lobpcg_iter_log,
      Rcpp::Named("eig_lobpcg_fail_iter") = hyb_ctx.lobpcg_fail_iter_log,
      Rcpp::Named("eig_lobpcg_max_rel_res") = hyb_ctx.lobpcg_max_rel_res_log);
#else
  return Rcpp::List::create(
      Rcpp::Named("imputed_values") = imputed_values,
      Rcpp::Named("mse") = mse,
      Rcpp::Named("solver_chosen") = hyb_ctx.solver_chosen_code(),
      Rcpp::Named("n_exact") = hyb_ctx.n_exact,
      Rcpp::Named("n_lobpcg_ok") = hyb_ctx.n_lobpcg_ok,
      Rcpp::Named("n_lobpcg_bad") = hyb_ctx.n_lobpcg_bad);
#endif
}
