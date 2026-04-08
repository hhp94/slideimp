#include "eig_sym_sel.h"
#include "loc_timer.h"
#if defined(_OPENMP)
#include <omp.h>
#endif

// ---------------------------------------------------------------------------
// eliminate the scale branch so inner loops auto-vectorize in the while hot
// loop.
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
                          const double *obs_sumsq,
                          const int cores)
{
  double *xptr = Xhat.memptr();
  const double *fptr = fittedX.memptr();
  const arma::uword n_mc = miss_cols.n_elem;
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
  for (arma::uword cidx = 0; cidx < n_mc; ++cidx)
  {
    const arma::uword j = miss_cols[cidx];
    const arma::uword off = j * nrX;
    double *col = xptr + off;
    const double *fcol = fptr + off;

    const double old_mu = mean_p(j);
    const double old_et = Scale ? et(j) : 1.0;
    // gather: write fitted values into Xhat at missing rows and accumulate their
    // unstandardized contributions.
    const arma::uword beg = miss_rows_offsets[cidx];
    const arma::uword end = miss_rows_offsets[cidx + 1];
    double imp_sum = 0.0;
    double imp_sumsq = 0.0;
    for (arma::uword p = beg; p < end; ++p)
    {
      const arma::uword i = miss_rows_flat[p];
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
    }

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
    arma::mat X,
    const arma::umat &miss,
    const arma::uword ncp,
    const bool scale,
    const bool regularized,
    const double threshold,
    const arma::uword init,
    const arma::uword maxiter,
    const arma::uword miniter,
    arma::rowvec row_w,
    const double coeff_ridge,
    const int cores = 1)
{
  LOC_TIMER_OBJ(pca_imp_gram);
  LOC_TIC(pca_imp_gram, "pca_imp_internal_cpp_total");

  const arma::uword nrX = X.n_rows;
  const arma::uword ncX = X.n_cols;
  const double nrX_d = static_cast<double>(nrX);
  const double ncX_d = static_cast<double>(ncX);
  double old = arma::datum::inf;
  const arma::uvec missing = arma::find(miss);
  X.elem(missing).zeros();
  row_w /= arma::accu(row_w);

  const arma::mat not_miss = arma::conv_to<arma::mat>::from(1 - miss);
  const arma::rowvec denom = row_w * not_miss;
  const arma::rowvec weighted_sum = row_w * X;
  const arma::rowvec sum_sq = row_w * (X % X);
  const double *wsptr = weighted_sum.memptr(); // neded for restandardization
  const double *ssptr = sum_sq.memptr();
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
  const double *xhp = Xhat.memptr();
  const double *fxp = fittedX.memptr();
  // constants for SVD path
  const bool tall = (nrX >= ncX);
  const arma::uword k = ncp + 1; // clamped in R
  // loop-carried buffers. Note: no need to size upfront since armadillo sizes
  // them during the first iteration and subsequence iterations are no op
  arma::mat U, V, X_work, AA_NxN, Vn, eigvecs;
  arma::vec vs_top, d_inv, eigvals, lambda_shrinked;
  // size of AA_NxN passed to the eig solver
  GramWorkspace gram_ws;
  if (!gram_ws.init(nrX, ncX, tall, k))
  {
    Rcpp::stop("workspace query failed");
  }
  AA_NxN.set_size(gram_ws.n, gram_ws.n);
  // set sizes of other buffers
  U.set_size(nrX, ncp);
  V.set_size(ncX, ncp);
  Vn.set_size(tall ? ncX : nrX, ncp);
  fittedX.set_size(nrX, ncX);
  if (tall)
  {
    X_work.set_size(nrX, ncX);
  }
  eigvecs.set_size(tall ? ncX : nrX, k);
  eigvals.set_size(k);
  vs_top.set_size(k);
  d_inv.set_size(ncp);
  lambda_shrinked.set_size(ncp);
  // pre-compute constants for sigma2
  const double ncp_d = static_cast<double>(ncp);
  const double min_dim = std::min(ncX_d, nrX_d - 1.0);
  const double denom_sigma = (ncX_d - ncp_d) * (nrX_d - 1.0 - ncp_d);
  const double scale_factor = (nrX_d * ncX_d) / min_dim;
  // pre-compute per-column missing structure and observed weighted sums. This
  // doesn't change so we can just calculate it once. Columns with zero missing
  // entries are skipped entirely in restandardize
  const arma::urowvec col_nmiss_row = arma::sum(miss, 0);
  arma::uvec miss_cols_idx = arma::find(col_nmiss_row > 0);
  const arma::uword n_mc = miss_cols_idx.n_elem;
  arma::uvec miss_rows_offsets(n_mc + 1, arma::fill::zeros);
  miss_rows_offsets.tail(n_mc) = arma::cumsum(col_nmiss_row.elem(miss_cols_idx));
  arma::uvec miss_rows_flat(miss_rows_offsets[n_mc]);
  const arma::uword *mp = miss.memptr();
  {
    for (arma::uword cidx = 0; cidx < n_mc; ++cidx)
    {
      const arma::uword j = miss_cols_idx[cidx];
      const arma::uword off = j * nrX;
      arma::uword p = miss_rows_offsets[cidx];
      for (arma::uword i = 0; i < nrX; ++i)
      {
        if (mp[off + i])
          miss_rows_flat[p++] = i;
      }
    }
  }
  arma::uvec col_iv_offsets(ncX + 1, arma::fill::zeros);
  for (arma::uword j = 0; j < ncX; ++j)
  {
    col_iv_offsets[j + 1] = col_iv_offsets[j] + col_nmiss_row[j];
  }
  // objective loop
  double objective = 0.0;
  arma::vec col_obj(ncX);
  double *cop = col_obj.memptr();
  // hot while loop
  for (arma::uword nb_iter = 1;; ++nb_iter)
  {
    LOC_TIC(pca_imp_gram, "restandardize");
    if (scale)
    {
      impute_restandardize<true>(Xhat, fittedX, mean_p, et, wptr, nrX,
                                 miss_cols_idx, miss_rows_flat, miss_rows_offsets,
                                 wsptr, ssptr, cores);
    }
    else
    {
      impute_restandardize<false>(Xhat, fittedX, mean_p, et, wptr, nrX,
                                  miss_cols_idx, miss_rows_flat, miss_rows_offsets,
                                  wsptr, ssptr, cores);
    }
    LOC_TOC(pca_imp_gram, "restandardize");

    LOC_TIC(pca_imp_gram, "svd");
    SVD_triplet(Xhat, sptr, isptr, ncp, tall,
                vs_top, trace_val, U, V, X_work, AA_NxN,
                d_inv, Vn, eigvals, eigvecs, gram_ws LOC_TIMER_ARG(pca_imp_gram));
    LOC_TOC(pca_imp_gram, "svd");

    LOC_TIC(pca_imp_gram, "post_svd");
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
    LOC_TOC(pca_imp_gram, "post_svd");

    LOC_TIC(pca_imp_gram, "reconstruct");
    scale_cols_inplace(U.memptr(), lambda_shrinked.memptr(), U.n_rows, ncp);
    fittedX = U * V.t();
    LOC_TOC(pca_imp_gram, "reconstruct");
    LOC_TIC(pca_imp_gram, "objective");
    {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(static)
#endif
      for (arma::uword j = 0; j < ncX; ++j)
      {
        const arma::uword off = j * nrX;
        const double *xh = xhp + off;
        const double *fx = fxp + off;
        const double *wm = wmp + off;
        double s = 0.0;
        for (arma::uword i = 0; i < nrX; ++i)
        {
          const double d = xh[i] - fx[i];
          s += wm[i] * d * d;
        }
        cop[j] = s;
      }
      objective = arma::accu(col_obj);
    }
    LOC_TOC(pca_imp_gram, "objective");

    // convergence
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
  LOC_TIC(pca_imp_gram, "postprocessing");
  arma::vec imputed_vals(missing.n_elem);
  double *ivp = imputed_vals.memptr();
  const double *xp = X.memptr();
  double sse = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sse) num_threads(cores) schedule(static)
#endif
  for (arma::uword j = 0; j < ncX; ++j)
  {
    const double mj = mean_p[j];
    const double ej = scale ? et[j] : 1.0;
    const arma::uword off = j * nrX;
    arma::uword iv_idx = col_iv_offsets[j];
    for (arma::uword i = 0; i < nrX; ++i)
    {
      const arma::uword kk = off + i;
      if (mp[kk])
      {
        ivp[iv_idx++] = xhp[kk] * ej + mj;
      }
      else
      {
        const double d = fxp[kk] * ej + mj - xp[kk];
        sse += d * d;
      }
    }
  }
  const double mse = sse / static_cast<double>(X.n_elem - missing.n_elem);
  LOC_TOC(pca_imp_gram, "postprocessing");

  LOC_TOC(pca_imp_gram, "pca_imp_internal_cpp_total");

  return Rcpp::List::create(
      Rcpp::Named("imputed_vals") = imputed_vals,
      Rcpp::Named("mse") = mse);
}
