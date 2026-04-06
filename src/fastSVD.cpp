#include "partial_eig.h"
// #include <rcpptimer.h>

constexpr double PCA_TOL = 1e-15;

// helpers to address the heavy overhead in svd copy creating functions like %=
static inline void scale_cols_inplace(double *M, const double *s,
                                      arma::uword nr, arma::uword nc)
{
  for (arma::uword j = 0; j < nc; ++j)
  {
    const double sj = s[j];
    double *col = M + j * nr;
    for (arma::uword i = 0; i < nr; ++i)
    {
      col[i] *= sj;
    }
  }
}

static inline void scale_rows_inplace(double *M, const double *w,
                                      arma::uword nr, arma::uword nc)
{
  for (arma::uword j = 0; j < nc; ++j)
  {
    double *col = M + j * nr;
    for (arma::uword i = 0; i < nr; ++i)
    {
      col[i] *= w[i];
    }
  }
}

// Copy source -> dest while scaling row i by w(i): dst(i,j) = src(i,j) * w(i).
static inline void copy_scale_rows(double *dst, const double *src,
                                   const double *w,
                                   arma::uword nr, arma::uword nc)
{
  for (arma::uword j = 0; j < nc; ++j)
  {
    double *dcol = dst + j * nr;
    const double *scol = src + j * nr;
    for (arma::uword i = 0; i < nr; ++i)
    {
      dcol[i] = scol[i] * w[i];
    }
  }
}

// ---------------------------------------------------------------------------
// SVD with symmetric eigs_sym on gram matrix
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
                        GramWorkspace &gram_ws)
{
  // Rcpp::Timer timer("SVD_triplet");
  const arma::uword nr = static_cast<arma::uword>(gram_ws.nr);
  const arma::uword nc = static_cast<arma::uword>(gram_ws.nc);

  // timer.tic("form_gram");
  if (tall)
  {
    // X_work = diag(sqrt_w) * Xhat, fused copy+scale (single pass)
    X_work.set_size(nr, nc);
    copy_scale_rows(X_work.memptr(), Xhat.memptr(), sw, nr, nc);
    syrk_upper(X_work, AA_NxN, gram_ws); // AA = X_work^T X_work
  }
  else
  {
    syrk_upper(Xhat, AA_NxN, gram_ws); // AA = Xhat Xhat^T (upper)
    // Apply D * AA * D on upper triangle only.
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
  // timer.toc("form_gram");

  trace_val = arma::trace(AA_NxN);

  // timer.tic("eig");
  if (!partial_eig_sym(eigvals, eigvecs, AA_NxN, gram_ws))
  {
    Rcpp::stop("`partial_eig_sym` failed to converge");
  }
  // timer.toc("eig");

  // Fused clamp+sqrt into vs_top, and build d_inv in the same pass.
  const arma::uword ke = eigvals.n_elem; // = gram_ws.topk
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

  // timer.tic("recover_factors");
  if (tall)
  {
    // V = eigvecs[:,0:ncp]; Vn(i,j) = V(i,j) * d_inv(j).
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
    U = X_work * Vn;
  }
  else
  {
    // U = eigvecs[:,0:ncp]; Vn(i,j) = U(i,j) * d_inv(j) * sqrt_w(i). Fused.
    const arma::uword n = eigvecs.n_rows; // == nr
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
    V = Xhat.t() * Vn;
  }
  // Undo row weighting on U in place.
  scale_rows_inplace(U.memptr(), isw, U.n_rows, U.n_cols);
  // timer.toc("recover_factors");
}

// ---------------------------------------------------------------------------
// eliminate the scale branch so inner loops auto-vectorize in the while hot loop
// ---------------------------------------------------------------------------
template <bool Scale>
void restandardize(arma::mat &Xhat,
                   arma::rowvec &mean_p,
                   arma::rowvec &et,
                   const double *wptr,
                   const arma::uword nrX,
                   const arma::uword ncX)
{
  double *xptr = Xhat.memptr();
  for (arma::uword j = 0; j < ncX; ++j)
  {
    double *col = xptr + j * nrX;
    const double mj = mean_p(j);
    double ej;
    if constexpr (Scale)
    {
      ej = et(j);
    }
    else
    {
      ej = 1.0;
    }
    // pass 1: unstandardize, accumulate weighted mean
    double new_mean = 0.0;
    double sum_sq = 0.0;
    for (arma::uword i = 0; i < nrX; ++i)
    {
      double val;
      if constexpr (Scale)
      {
        val = col[i] * ej + mj;
      }
      else
      {
        val = col[i] + mj;
      }
      col[i] = val;
      const double w = wptr[i];
      new_mean += w * val;
      if constexpr (Scale)
      {
        sum_sq += w * val * val;
      }
    }
    // pass 2: re-center (+re-scale)
    if constexpr (Scale)
    {
      const double var = sum_sq - new_mean * new_mean;
      const double new_et = std::sqrt(std::max(0.0, var));
      const double inv_et = 1.0 / std::max(new_et, PCA_TOL);
      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] = (col[i] - new_mean) * inv_et;
      }
      et(j) = new_et;
    }
    else
    {
      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] -= new_mean;
      }
    }
    mean_p(j) = new_mean;
  }
}

// ---------------------------------------------------------------------------
// eliminate the scale branch for the post processing, which takes a suprising
// amount of time when the data is large and there's many missing values
// ---------------------------------------------------------------------------
template <bool Scale>
static inline void finalize_impute(
    const arma::mat &X, const arma::mat &Xhat, const arma::mat &fittedX,
    const arma::umat &miss, const arma::rowvec &et, const arma::rowvec &mean_p,
    arma::vec &imputed_vals, double &sse_out)
{
  const arma::uword nrX = X.n_rows;
  const arma::uword ncX = X.n_cols;
  const double *xp = X.memptr();
  const double *xhp = Xhat.memptr();
  const double *fp = fittedX.memptr();
  const arma::uword *mp = miss.memptr();
  const double *mnp = mean_p.memptr();
  const double *etp = Scale ? et.memptr() : nullptr;
  double *ivp = imputed_vals.memptr();
  arma::uword iv_idx = 0;
  double sse = 0.0;

  for (arma::uword j = 0; j < ncX; ++j)
  {
    const double mj = mnp[j];
    const arma::uword off = j * nrX;
    if constexpr (Scale)
    {
      const double ej = etp[j];
      for (arma::uword i = 0; i < nrX; ++i)
      {
        const arma::uword k = off + i;
        if (mp[k])
        {
          ivp[iv_idx++] = xhp[k] * ej + mj;
        }
        else
        {
          const double d = fp[k] * ej + mj - xp[k];
          sse += d * d;
        }
      }
    }
    else
    {
      for (arma::uword i = 0; i < nrX; ++i)
      {
        const arma::uword k = off + i;
        if (mp[k])
        {
          ivp[iv_idx++] = xhp[k] + mj;
        }
        else
        {
          const double d = fp[k] + mj - xp[k];
          sse += d * d;
        }
      }
    }
  }
  sse_out = sse;
}

// [[Rcpp::export]]
Rcpp::List pca_imp_internal_cpp(
    const arma::mat &X,
    const arma::umat &miss,
    const arma::uword ncp,
    const bool scale,
    const bool regularized,
    const double threshold,
    const arma::uword init,
    const arma::uword maxiter,
    const arma::uword miniter,
    const arma::rowvec &row_w,
    const double coeff_ridge)
{
  const arma::uword nrX = X.n_rows;
  const arma::uword ncX = X.n_cols;
  const double nrX_d = static_cast<double>(nrX);
  const double ncX_d = static_cast<double>(ncX);
  double old = arma::datum::inf;
  double objective = 0.0;
  const arma::uvec missing = arma::find(miss);

  const arma::mat not_miss = arma::conv_to<arma::mat>::from(1 - miss);
  const arma::rowvec denom = row_w * not_miss;
  const arma::rowvec weighted_sum = row_w * X;
  const arma::rowvec sum_sq = row_w * (X % X);
  arma::rowvec mean_p = weighted_sum / denom;
  arma::rowvec et = arma::sqrt(sum_sq / denom - mean_p % mean_p);

  arma::mat Xhat = X.each_row() - mean_p;
  if (scale)
  {
    Xhat.each_row() /= et;
  }
  if (init == 0)
  {
    Xhat.elem(missing).zeros();
  }
  else
  {
    Rcpp::NumericVector rand_vals = Rcpp::rnorm(missing.n_elem);
    Xhat.elem(missing) = Rcpp::as<arma::vec>(rand_vals);
  }

  arma::mat fittedX = Xhat;
  double trace_val;
  double sigma2 = 0.0;
  // compute once, reuse every iteration
  const arma::colvec row_w_col = row_w.t();
  const arma::colvec sqrt_row_w = arma::sqrt(row_w_col);
  const arma::colvec inv_sqrt_row_w = 1.0 / sqrt_row_w;
  const double *sptr = sqrt_row_w.memptr();
  const double *isptr = inv_sqrt_row_w.memptr();
  const double *wptr = row_w.memptr();
  // precompute combined weight*mask for the hot objective loop:
  //   w_mask(i,j) = row_w(i) * not_miss(i,j)
  // used as: objective += w_mask(i,j) * (Xhat(i,j) - fittedX(i,j))^2
  arma::mat w_mask = not_miss;
  w_mask.each_col() %= row_w_col;
  const double *wmp = w_mask.memptr();
  // constants for SVD path
  const bool tall = (nrX >= ncX);
  const arma::uword k = ncp + 1; // clamped in R
  // loop-carried buffers. Note: no need to size upfront since armadillo sizes
  // them during the first iteration and subsequence iterations are no op
  arma::mat U, V, X_work, AA_NxN, Vn, eigvecs, diff;
  arma::vec vs_top, d_inv, eigvals, lambda_shrinked;
  // size of AA_NxN passed to the eig solver
  GramWorkspace gram_ws;
  if (!gram_ws.init(nrX, ncX, tall, k))
  {
    Rcpp::stop("workspace query failed");
  }
  AA_NxN.set_size(gram_ws.n, gram_ws.n);
  // pre-compute constants for sigma2
  const double ncp_d = static_cast<double>(ncp);
  const double min_dim = std::min(ncX_d, nrX_d - 1.0);
  const double denom_sigma = (ncX_d - ncp_d) * (nrX_d - 1.0 - ncp_d);
  const double scale_factor = (nrX_d * ncX_d) / min_dim;
  // Rcpp::Timer timer("pca_imp_gram");
  // Rcpp::Timer::ScopedTimer scpdtmr(timer, "pca_imp_internal_cpp_total"); // measure total time
  for (arma::uword nb_iter = 1;; ++nb_iter)
  {
    Xhat.elem(missing) = fittedX.elem(missing);

    // timer.tic("restandardize");
    if (scale)
    {
      restandardize<true>(Xhat, mean_p, et, wptr, nrX, ncX);
    }
    else
    {
      restandardize<false>(Xhat, mean_p, et, wptr, nrX, ncX);
    }
    // timer.toc("restandardize");
    // timer.tic("svd");
    SVD_triplet(Xhat, sptr, isptr, ncp, tall,
                vs_top, trace_val, U, V, X_work, AA_NxN,
                d_inv, Vn, eigvals, eigvecs, gram_ws);
    // timer.toc("svd");
    // timer.tic("post_svd");
    sigma2 = 0.0;
    if (regularized)
    {
      const double sum_top_sq = arma::dot(vs_top.head(ncp), vs_top.head(ncp));
      const double sum_tail_sq = trace_val - sum_top_sq;
      sigma2 = scale_factor * sum_tail_sq / denom_sigma;
      sigma2 = std::min(sigma2 * coeff_ridge, vs_top(ncp) * vs_top(ncp));
    }
    const auto vs_h = vs_top.head(ncp);
    lambda_shrinked = (vs_h % vs_h - sigma2) /
                      arma::clamp(vs_h, PCA_TOL, arma::datum::inf);
    // timer.toc("post_svd");
    // timer.tic("reconstruct");
    scale_cols_inplace(U.memptr(), lambda_shrinked.memptr(), U.n_rows, ncp);
    fittedX = U * V.t();
    // timer.toc("reconstruct");

    // timer.tic("objective");
    {
      const double *xp = Xhat.memptr();
      const double *fp = fittedX.memptr();
      const arma::uword ntot = nrX * ncX;
      double acc = 0.0;
      for (arma::uword idx = 0; idx < ntot; ++idx)
      {
        const double d = xp[idx] - fp[idx];
        acc += wmp[idx] * d * d;
      }
      objective = acc;
    }
    // timer.toc("objective");
    // convergence
    const double criterion = (old > 0.0 && !std::isinf(old))
                                 ? std::abs(1.0 - objective / old)
                                 : 1.0;
    old = objective;

    if (criterion < threshold && nb_iter > miniter)
    {
      break;
    }
    if (nb_iter > maxiter)
    {
      Rcpp::warning("Stopped after " + std::to_string(maxiter) + " iterations");
      break;
    }
  }
  // timer.tic("postprocessing");
  // final unstandardize
  arma::vec imputed_vals(missing.n_elem);
  double sse = 0.0;
  if (scale)
  {
    finalize_impute<true>(X, Xhat, fittedX, miss, et, mean_p, imputed_vals, sse);
  }
  else
  {
    finalize_impute<false>(X, Xhat, fittedX, miss, et, mean_p, imputed_vals, sse);
  }
  const double mse = sse / static_cast<double>(X.n_elem - missing.n_elem);
  // timer.toc("postprocessing");

  return Rcpp::List::create(
      Rcpp::Named("imputed_vals") = imputed_vals,
      Rcpp::Named("mse") = mse);
}
