#include "eig_sym_sel.h"
#include "imputed_value.h"
#include "loc_timer.h"

// small helper for more readable codes that deals with the CSC
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
// impute_restandardize: column indices in miss_cols_idx and miss_rows_flat
// are local positions into Xhat (NOT original obj indices)
// ---------------------------------------------------------------------------
template <bool Scale>
void impute_restandardize(arma::mat &Xhat,
                          const arma::mat &fittedX,
                          arma::rowvec &mean_p,
                          arma::rowvec &et,
                          const double *wptr,
                          const arma::uword nrX,
                          const arma::uvec &miss_cols,
                          const arma::uvec &miss_rows_flat,
                          const arma::uvec &miss_rows_offsets,
                          const double *obs_sum,
                          const double *obs_sumsq)
{
  double *xptr = Xhat.memptr();
  const double *fptr = fittedX.memptr();
  const arma::uword n_mc = miss_cols.n_elem;

  for (arma::uword cidx = 0; cidx < n_mc; ++cidx)
  {
    const arma::uword j = miss_cols[cidx];
    const arma::uword off = j * nrX;
    double *col = xptr + off;
    const double *fcol = fptr + off;

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
                              const double v_std = fcol[i];
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
      const double inv_new = 1.0 / std::max(new_et, PCA_TOL);
      const double a = old_et * inv_new;
      const double b = (old_mu - new_mu) * inv_new;
      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] = col[i] * a + b;
      }
      et(j) = new_et;
    }
    else
    {
      const double b = old_mu - new_mu;
      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] += b;
      }
    }
    mean_p(j) = new_mu;
  }
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
    const double coeff_ridge)
{
  LOC_TIMER_OBJ(pca_imp_gram);
  LOC_TIC(pca_imp_gram, "pca_imp_internal_cpp_total");
  stop_on_inf(obj);
  const arma::uword nrX = obj.n_rows;
  const arma::uword n_elig = eligible_idx.n_elem;
  const double nrX_d = static_cast<double>(nrX);
  const double n_elig_d = static_cast<double>(n_elig);
  double old = arma::datum::inf;

  row_w /= arma::accu(row_w);
  const double *wptr = row_w.memptr();

  LOC_TIC(pca_imp_gram, "pass1_count_missing");
  // ---------------------------------------------------------------------------
  // pass 1: count NaN per eligible column; build CSC offsets
  // ---------------------------------------------------------------------------
  arma::urowvec col_nmiss(n_elig, arma::fill::zeros);
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
  arma::uvec miss_cols_idx = arma::find(col_nmiss > 0); // sorted ascending
  const arma::uword n_mc = miss_cols_idx.n_elem;
  arma::uvec miss_rows_offsets(n_mc + 1, arma::fill::zeros);
  miss_rows_offsets.tail(n_mc) = arma::cumsum(col_nmiss.elem(miss_cols_idx));
  const arma::uword n_missing = miss_rows_offsets[n_mc];
  arma::uvec miss_rows_flat(n_missing);
  LOC_TOC(pca_imp_gram, "pass1_count_missing");

  // ---------------------------------------------------------------------------
  // pass 2: fill Xhat (NaN -> 0, else v), miss_rows_flat, and running sums
  // over observed rows. Single column-major sweep. Since miss_cols_idx is
  // sorted and we iterate j_local in order, cidx advances lazily to avoid a
  // per-column binary search.
  // ---------------------------------------------------------------------------
  arma::mat Xhat(nrX, n_elig);
  arma::rowvec denom(n_elig, arma::fill::zeros);
  arma::rowvec weighted_sum(n_elig, arma::fill::zeros);
  arma::rowvec sum_sq(n_elig, arma::fill::zeros);

  LOC_TIC(pca_imp_gram, "pass2_scan");
  {
    arma::uword cidx = 0;
    for (arma::uword j_local = 0; j_local < n_elig; ++j_local)
    {
      const double *src = obj.colptr(eligible_idx[j_local]);
      double *dst = Xhat.colptr(j_local);

      const bool has_miss = (cidx < n_mc) && (miss_cols_idx[cidx] == j_local);
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
      denom[j_local] = d_acc;
      weighted_sum[j_local] = ws_acc;
      sum_sq[j_local] = ss_acc;
      if (has_miss)
      {
        ++cidx;
      }
    }
  }
  LOC_TOC(pca_imp_gram, "pass2_scan");

  const double *wsptr = weighted_sum.memptr(); // reused inside restandardize
  const double *ssptr = sum_sq.memptr();
  arma::rowvec mean_p = weighted_sum / denom;
  arma::rowvec et;
  if (scale)
  {
    et = arma::sqrt(sum_sq / denom - mean_p % mean_p);
  }
  // center (and optionally scale) the entire Xhat. This contaminates missing
  // positions (currently 0) with -mean_p / et. We overwrite them immediately
  // below with either 0 or rnorm according to missMDA init behavior.
  Xhat.each_row() -= mean_p;
  if (scale)
  {
    Xhat.each_row() /= et;
  }
  if (init == 0)
  {
    for (arma::uword ci = 0; ci < n_mc; ++ci)
    {
      double *col = Xhat.colptr(miss_cols_idx[ci]);
      for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                              [&](arma::uword i)
                              { col[i] = 0.0; });
    }
  }
  else if (n_missing > 0)
  {
    Rcpp::NumericVector rand_vals = Rcpp::rnorm(n_missing);
    arma::uword idx = 0;
    for (arma::uword ci = 0; ci < n_mc; ++ci)
    {
      double *col = Xhat.colptr(miss_cols_idx[ci]);
      for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                              [&](arma::uword i)
                              { col[i] = rand_vals[idx++]; });
    }
  }

  // ---------------------------------------------------------------------------
  // SVD / solve buffers, sized on n_elig instead of ncX.
  // ---------------------------------------------------------------------------
  arma::mat fittedX = Xhat; // seed first restandardize with Xhat's current state
  double trace_val;
  double sigma2 = 0.0;
  const arma::colvec row_w_col = row_w.t();
  const arma::colvec sqrt_row_w = arma::sqrt(row_w_col);
  const arma::colvec inv_sqrt_row_w = 1.0 / sqrt_row_w;
  const double *sptr = sqrt_row_w.memptr();
  const double *isptr = inv_sqrt_row_w.memptr();
  const bool tall = (nrX >= n_elig);
  const arma::uword k = ncp + 1; // clamped in R

  arma::mat U, V, X_work, AA_NxN, Vn, eigvecs;
  arma::vec vs_top, d_inv, eigvals, lambda_shrinked;
  GramWorkspace gram_ws;
  if (!gram_ws.init(nrX, n_elig, tall, k))
  {
    Rcpp::stop("workspace query failed");
  }
  AA_NxN.set_size(gram_ws.n, gram_ws.n);
  U.set_size(nrX, ncp);
  V.set_size(n_elig, ncp);
  Vn.set_size(tall ? n_elig : nrX, ncp);
  // intentional no-op: fittedX was copy-constructed from Xhat above at
  // exactly (nrX, n_elig) and seeds the first impute_restandardize. In
  // current Armadillo, set_size with matching dims preserves values
  // (init_warm returns early). DO NOT change to fittedX.zeros(...): that
  // would corrupt the iter-1 seed and diverge from missMDA.
  fittedX.set_size(nrX, n_elig);
  fittedX.set_size(nrX, n_elig);
  if (tall)
  {
    X_work.set_size(nrX, n_elig);
  }
  eigvecs.set_size(tall ? n_elig : nrX, k);
  eigvals.set_size(k);
  vs_top.set_size(k);
  d_inv.set_size(ncp);
  lambda_shrinked.set_size(ncp);

  const double ncp_d = static_cast<double>(ncp);
  const double min_dim = std::min(n_elig_d, nrX_d - 1.0);
  const double denom_sigma = (n_elig_d - ncp_d) * (nrX_d - 1.0 - ncp_d);
  const double scale_factor = (nrX_d * n_elig_d) / min_dim;

  // ---------------------------------------------------------------------------
  // hot loop
  // ---------------------------------------------------------------------------
  double objective = 0.0;
  const double *xhp = Xhat.memptr(); // stable: Xhat modified in place
  double *fxp = fittedX.memptr();    // allow to refetch

  for (arma::uword nb_iter = 1;; ++nb_iter)
  {
    LOC_TIC(pca_imp_gram, "restandardize");
    if (scale)
    {
      impute_restandardize<true>(Xhat, fittedX, mean_p, et, wptr, nrX,
                                 miss_cols_idx, miss_rows_flat, miss_rows_offsets,
                                 wsptr, ssptr);
    }
    else
    {
      impute_restandardize<false>(Xhat, fittedX, mean_p, et, wptr, nrX,
                                  miss_cols_idx, miss_rows_flat, miss_rows_offsets,
                                  wsptr, ssptr);
    }
    LOC_TOC(pca_imp_gram, "restandardize");

    LOC_TIC(pca_imp_gram, "svd");
    SVD_triplet(Xhat, sptr, isptr, ncp, tall,
                vs_top, trace_val, U, V, X_work, AA_NxN,
                d_inv, Vn, eigvals, eigvecs, gram_ws LOC_TIMER_ARG(pca_imp_gram));
    LOC_TOC(pca_imp_gram, "svd");

    LOC_TIC(pca_imp_gram, "post_svd");
    sigma2 = 0.0;
    const double sum_top_sq = arma::dot(vs_top.head(ncp), vs_top.head(ncp));
    const double sum_tail_sq = trace_val - sum_top_sq;
    if (regularized)
    {
      sigma2 = scale_factor * sum_tail_sq / denom_sigma;
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
        const double d = v - lam; // = sigma2/v in the regularized path,
        top_corr += d * d;        //   = v when we clamped lam to 0.
      }
    }
    // full-matrix weighted SSE via SVD identity (O(ncp), not O(n*p)):
    const double objective_full = std::max(0.0, sum_tail_sq) + top_corr;
    LOC_TOC(pca_imp_gram, "post_svd");

    LOC_TIC(pca_imp_gram, "reconstruct");
    scale_cols_inplace(U.memptr(), lambda_shrinked.memptr(), U.n_rows, ncp);
    fittedX = U * V.t();
    // refetch: operator= reuses the buffer when dims match (they do here),
    // but this decouples fxp from that Armadillo contract.
    fxp = fittedX.memptr();
    LOC_TOC(pca_imp_gram, "reconstruct");

    LOC_TIC(pca_imp_gram, "objective");
    {
      double miss_contrib = 0.0;
      for (arma::uword ci = 0; ci < n_mc; ++ci)
      {
        const arma::uword j_local = miss_cols_idx[ci];
        const arma::uword off = j_local * nrX;
        const double *xh = xhp + off;
        const double *fx = fxp + off;
        double s = 0.0;
        for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                                [&](arma::uword i)
                                {
                                  const double d = xh[i] - fx[i];
                                  s += wptr[i] * d * d;
                                });
        miss_contrib += s;
      }
      objective = objective_full - miss_contrib;
    }
    LOC_TOC(pca_imp_gram, "objective");

    const double criterion = std::isinf(old) ? 1.0 : std::abs(1.0 - objective / old);
    old = objective;

    if (criterion < threshold && nb_iter >= miniter)
    {
      break;
    }
    if (nb_iter >= maxiter)
    {
      Rcpp::warning("Stopped after " + std::to_string(maxiter) + " iterations");
      break;
    }
  }

  // ---------------------------------------------------------------------------
  // Postprocessing.
  //
  // SSE in original (unstandardized) scale at observed positions uses the
  // identity that holds exactly at every iteration:
  //  Xhat[i,j] = (X[i,j] - mean_p[j]) / et[j] (scale = true, obs i)
  //  Xhat[i,j] =  X[i,j] - mean_p[j] (scale = false, obs i)
  //
  // Therefore at observed positions:
  //  (fittedX[i,j] * ej + mean_p[j] - X[i,j])^2 = (fittedX - Xhat)^2 * ej^2
  // with ej = et[j] if scale else 1.
  // Structure mirrors the objective loop: dense full sum per column, then
  // subtract missing contributions via CSC.
  // ---------------------------------------------------------------------------
  LOC_TIC(pca_imp_gram, "postprocessing");

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

  // emit imputed triples in original 1-based (row, col) coordinates.
  // j_local -> eligible_idx[j_local] + 1 (back to original column space).
  arma::mat imputed_values(n_missing, 3);
  {
    arma::uword iv = 0;
    for (arma::uword ci = 0; ci < n_mc; ++ci)
    {
      const arma::uword j_local = miss_cols_idx[ci];
      const double ej = scale ? et[j_local] : 1.0;
      const double mj = mean_p[j_local];
      const double orig_col_1based =
          static_cast<double>(eligible_idx[j_local] + 1);
      const arma::uword off = j_local * nrX;
      const double *xh = xhp + off;
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
  // // fittedX back to original scale. stays (nrX x n_elig)
  // if (scale)
  // {
  //   scale_cols_inplace(fittedX.memptr(), et.memptr(), nrX, n_elig);
  // }
  // fittedX.each_row() += mean_p;
  // const arma::uvec eligible_cols_R = eligible_idx + 1;

  LOC_TOC(pca_imp_gram, "postprocessing");
  LOC_TOC(pca_imp_gram, "pca_imp_internal_cpp_total");

  return Rcpp::List::create(
      Rcpp::Named("imputed_values") = imputed_values,
      // Rcpp::Named("fittedX") = fittedX,
      // Rcpp::Named("eligible_cols") = eligible_cols_R,
      Rcpp::Named("mse") = mse);
}
