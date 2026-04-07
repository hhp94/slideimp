#include "eig_sym_sel.h"
#include "loc_timer.h"

// ---------------------------------------------------------------------------
// eliminate the scale branch so inner loops auto-vectorize in the while hot loop
// ---------------------------------------------------------------------------
template <bool Scale>
void restandardize(arma::mat &Xhat,
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
  const arma::uword n_mc = miss_cols.n_elem;
  for (arma::uword cidx = 0; cidx < n_mc; ++cidx)
  {
    const arma::uword j = miss_cols[cidx];
    double *col = xptr + j * nrX;
    const double mj = mean_p(j);
    const double ej = Scale ? et(j) : 1.0;
    // pass 1: unstandardize entire column
    if constexpr (Scale)
    {
      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] = col[i] * ej + mj;
      }
    }
    else
    {
      for (arma::uword i = 0; i < nrX; ++i)
      {
        col[i] += mj;
      }
    }
    // gather: accumulate only over missing rows
    const arma::uword beg = miss_rows_offsets[cidx];
    const arma::uword end = miss_rows_offsets[cidx + 1];
    double imp_sum = 0.0;
    double imp_sumsq = 0.0;
    for (arma::uword p = beg; p < end; ++p)
    {
      const arma::uword i = miss_rows_flat[p];
      const double v = col[i];
      const double w = wptr[i];
      imp_sum += w * v;
      if constexpr (Scale)
      {
        imp_sumsq += w * v * v;
      }
    }
    const double new_mean = obs_sum[j] + imp_sum;
    // pass 2: re-standardize
    if constexpr (Scale)
    {
      const double new_var = obs_sumsq[j] + imp_sumsq - new_mean * new_mean;
      const double new_et = std::sqrt(std::max(0.0, new_var));
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
        col[i] -= new_mean;
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
  LOC_TIMER_OBJ(pca_imp_gram);
  LOC_TIC(pca_imp_gram, "pca_imp_internal_cpp_total");

  const arma::uword nrX = X.n_rows;
  const arma::uword ncX = X.n_cols;
  const arma::uword ntot = nrX * ncX;
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
  arma::uvec miss_cols_idx;
  {
    std::vector<arma::uword> tmp;
    tmp.reserve(ncX);
    for (arma::uword j = 0; j < ncX; ++j)
    {
      if (col_nmiss_row[j] > 0)
      {
        tmp.push_back(j);
      }
    }
    miss_cols_idx = arma::uvec(tmp);
  }
  const arma::uword n_mc = miss_cols_idx.n_elem;
  arma::uvec miss_rows_offsets(n_mc + 1);
  miss_rows_offsets[0] = 0;
  for (arma::uword c = 0; c < n_mc; ++c) {
    miss_rows_offsets[c + 1] = miss_rows_offsets[c] + col_nmiss_row[miss_cols_idx[c]];
  }
  arma::uvec miss_rows_flat(miss_rows_offsets[n_mc]);
  arma::vec obs_sum_vec(ncX, arma::fill::zeros);
  arma::vec obs_sumsq_vec(ncX, arma::fill::zeros);
  {
    const double *xp_init = X.memptr();
    const arma::uword *mp_init = miss.memptr();
    for (arma::uword cidx = 0; cidx < n_mc; ++cidx)
    {
      const arma::uword j = miss_cols_idx[cidx];
      const arma::uword off = j * nrX;
      arma::uword p = miss_rows_offsets[cidx];
      double os = 0.0, oss = 0.0;
      for (arma::uword i = 0; i < nrX; ++i)
      {
        if (mp_init[off + i])
        {
          miss_rows_flat[p++] = i;
        }
        else
        {
          const double w = wptr[i];
          const double x = xp_init[off + i];
          os += w * x;
          oss += w * x * x;
        }
      }
      obs_sum_vec[j] = os;
      obs_sumsq_vec[j] = oss;
    }
  }
  const double *obs_sum_ptr = obs_sum_vec.memptr();
  const double *obs_sumsq_ptr = obs_sumsq_vec.memptr();
  // hot while loop
  for (arma::uword nb_iter = 1;; ++nb_iter)
  {
    LOC_TIC(pca_imp_gram, "restandardize");
    Xhat.elem(missing) = fittedX.elem(missing);
    if (scale)
    {
      restandardize<true>(Xhat, mean_p, et, wptr, nrX,
                          miss_cols_idx, miss_rows_flat, miss_rows_offsets,
                          obs_sum_ptr, obs_sumsq_ptr);
    }
    else
    {
      restandardize<false>(Xhat, mean_p, et, wptr, nrX,
                           miss_cols_idx, miss_rows_flat, miss_rows_offsets,
                           obs_sum_ptr, obs_sumsq_ptr);
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

    LOC_TIC(pca_imp_gram, "reconstruct+objective");
    scale_cols_inplace(U.memptr(), lambda_shrinked.memptr(), U.n_rows, ncp);
    fittedX = U * V.t();
    {
      double acc = 0.0;
      for (arma::uword idx = 0; idx < ntot; ++idx)
      {
        const double d = xhp[idx] - fxp[idx];
        acc += wmp[idx] * d * d;
      }
      objective = acc;
    }
    LOC_TOC(pca_imp_gram, "reconstruct+objective");

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
  LOC_TIC(pca_imp_gram, "postprocessing");
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
  LOC_TOC(pca_imp_gram, "postprocessing");

  LOC_TOC(pca_imp_gram, "pca_imp_internal_cpp_total");

  return Rcpp::List::create(
      Rcpp::Named("imputed_vals") = imputed_vals,
      Rcpp::Named("mse") = mse);
}
