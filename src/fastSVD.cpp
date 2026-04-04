#include "partial_eig.h"

// symmetric eigendecomposition on the cross-product matrix
void fastSVD_triplet(const arma::mat &Xhat,
                     const arma::colvec &sqrt_row_w,
                     const arma::colvec &inv_sqrt_row_w,
                     arma::uword ncp,
                     bool tall,
                     arma::vec &vs_top,
                     double &trace_val,
                     arma::mat &U,
                     arma::mat &V,
                     arma::mat &X_work,
                     arma::mat &AA_NxN,
                     bool partial,
                     double tol = 1e-15)
{
  // weight rows: X_work = diag(sqrt_row_w) * Xhat
  X_work = Xhat;
  X_work.each_col() %= sqrt_row_w;

  // form the smaller cross-product
  if (tall)
  {
    AA_NxN = X_work.t() * X_work; // cols x cols
  }
  else
  {
    AA_NxN = X_work * X_work.t(); // rows x rows
  }

  arma::uword dim = AA_NxN.n_rows;
  arma::uword k = std::min(ncp + 1, dim);
  arma::vec eigvals;
  arma::mat eigvecs;
  bool success;

  if (partial)
  {
    success = partial_eig_sym(eigvals, eigvecs, AA_NxN, k);
    if (!success)
    {
      Rcpp::stop("`partial_eig_sym` failed to converge");
    }
  }
  else
  {
    success = arma::eig_sym(eigvals, eigvecs, AA_NxN, "dc");
    if (!success)
    {
      Rcpp::stop("`arma::eig_sym` failed to converge");
    }
    // take largest k (already ascending -> tail + reverse)
    eigvals = eigvals.tail(k);
    eigvecs = eigvecs.tail_cols(k);
    eigvals = arma::reverse(eigvals);
    eigvecs = arma::fliplr(eigvecs);
  }

  // trace of the cross-product (needed for sigma2 estimate)
  trace_val = arma::trace(AA_NxN);

  // singular values from eigenvalues (clamp negatives from numerical noise)
  vs_top = arma::sqrt(arma::clamp(eigvals, 0.0, arma::datum::inf));

  // extract top ncp eigenvectors and corresponding singular values
  arma::mat Vn = eigvecs.head_cols(ncp);
  arma::vec d = vs_top.head(ncp);

  // safe reciprocal of singular values
  arma::vec d_inv(ncp, arma::fill::zeros);
  for (arma::uword i = 0; i < ncp; ++i)
  {
    if (d(i) > tol)
    {
      d_inv(i) = 1.0 / d(i);
    }
  }

  // recover the other factor: U = X * V * D^{-1}  or  V = X' * U * D^{-1}
  if (tall)
  {
    // AA_NxN was X'X (cols x cols), eigvecs are right singular vectors
    U = X_work * (Vn.each_row() % d_inv.t());
    V = Vn;
  }
  else
  {
    // AA_NxN was XX' (rows x rows), eigvecs are left singular vectors
    V = X_work.t() * (Vn.each_row() % d_inv.t());
    U = Vn;
  }

  // undo row weighting on U: U = diag(1/sqrt_row_w) * U
  U.each_col() %= inv_sqrt_row_w;
}

// ---------------------------------------------------------------------------
// template to eliminate the scale branch so inner loops auto-vectorize in the
// while hot loop
// ---------------------------------------------------------------------------
template <bool Scale>
void restandardize(arma::mat &Xhat,
                   arma::rowvec &mean_p,
                   arma::rowvec &et,
                   const double *wptr,
                   arma::uword nrX,
                   arma::uword ncX)
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
      double w = wptr[i];
      new_mean += w * val;
      if constexpr (Scale)
      {
        sum_sq += w * val * val;
      }
    }
    // pass 2: re-center (+ re-scale)
    if constexpr (Scale)
    {
      double var = sum_sq - new_mean * new_mean;
      double new_et = std::sqrt(std::max(0.0, var));
      new_et = std::max(new_et, 0.0);
      const double inv_et = 1.0 / new_et;
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

// [[Rcpp::export]]
Rcpp::List pca_imp_internal_cpp(
    const arma::mat &X,
    const arma::umat &miss,
    const arma::uword ncp,
    bool scale,
    bool regularized,
    double threshold,
    arma::uword init,
    arma::uword maxiter,
    arma::uword miniter,
    const arma::rowvec &row_w,
    double coeff_ridge,
    bool partial)
{
  arma::uword nrX = X.n_rows;
  arma::uword ncX = X.n_cols;
  double nrX_d = static_cast<double>(nrX);
  double ncX_d = static_cast<double>(ncX);
  arma::uword nb_iter = 1;
  double old = arma::datum::inf;
  double objective = 0.0;
  arma::uvec missing = arma::find(miss);

  arma::mat not_miss = arma::conv_to<arma::mat>::from(1 - miss);
  arma::rowvec denom = row_w * not_miss;
  arma::rowvec weighted_sum = row_w * X;
  arma::rowvec sum_sq = row_w * (X % X);
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
  arma::mat U(nrX, ncp, arma::fill::zeros);
  arma::mat V(ncX, ncp, arma::fill::zeros);
  arma::vec vs_top;
  double trace_val;
  double min_dim = std::min(ncX_d, nrX_d - 1.0);
  double sigma2 = 0.0;
  arma::vec lambda_shrinked(ncp);

  // pre-allocate loop temporary
  arma::mat diff(nrX, ncX);

  // compute once, reuse every iteration
  arma::colvec sqrt_row_w = arma::sqrt(row_w.t());
  arma::colvec inv_sqrt_row_w = 1.0 / sqrt_row_w;
  arma::colvec row_w_col = row_w.t();
  const double *wptr = row_w_col.memptr();

  // buffers for the SVD routine
  arma::mat X_work(nrX, ncX);
  bool tall = (nrX >= ncX);
  arma::uword small_dim = tall ? ncX : nrX;
  arma::mat AA_NxN(small_dim, small_dim);

  // pre-compute constants for sigma2
  double denom_sigma = (nrX_d - 1.0) * ncX_d - (nrX_d - 1.0) * ncp - ncX_d * ncp + static_cast<double>(ncp * ncp);
  double scale_factor = (nrX_d * ncX_d) / min_dim;

  while (nb_iter > 0)
  {
    // update Xhat
    Xhat.elem(missing) = fittedX.elem(missing);
    if (scale)
    {
      restandardize<true>(Xhat, mean_p, et, wptr, nrX, ncX);
    }
    else
    {
      restandardize<false>(Xhat, mean_p, et, wptr, nrX, ncX);
    }
    // SVD
    fastSVD_triplet(Xhat, sqrt_row_w, inv_sqrt_row_w, ncp, tall,
                    vs_top, trace_val, U, V, X_work, AA_NxN, partial);
    // sigma2 via trace identity
    sigma2 = 0.0;
    if (regularized)
    {
      double sum_top_sq = arma::dot(vs_top.head(ncp), vs_top.head(ncp));
      double sum_tail_sq = trace_val - sum_top_sq;
      sigma2 = scale_factor * sum_tail_sq / denom_sigma;
      sigma2 = std::min(sigma2 * coeff_ridge, vs_top(ncp) * vs_top(ncp));
    }

    // shrinkage: lambda = (d^2 - sigma2)
    arma::vec vs_head = vs_top.head(ncp);
    lambda_shrinked = (vs_head % vs_head - sigma2) /
                      arma::clamp(vs_head, 1e-15, arma::datum::inf);

    // fitted reconstruction (U is overwritten next iteration)
    U.each_row() %= lambda_shrinked.t();
    fittedX = U * V.t();

    // objective on observed entries (contiguous mask, no scatter)
    diff = Xhat - fittedX;
    diff %= not_miss;
    objective = arma::dot(row_w_col, arma::sum(diff % diff, 1));

    // convergence
    double criterion = (old > 0.0 && !std::isinf(old)) ? std::abs(1.0 - objective / old) : 1.0;
    old = objective;
    ++nb_iter;

    if (criterion < threshold && nb_iter > miniter)
    {
      nb_iter = 0;
    }
    if (nb_iter > maxiter)
    {
      nb_iter = 0;
      Rcpp::warning("Stopped after " + std::to_string(maxiter) + " iterations");
    }
  }

  // final unstandardize
  if (scale)
  {
    Xhat.each_row() %= et;
  }
  Xhat.each_row() += mean_p;
  if (scale)
  {
    fittedX.each_row() %= et;
  }
  fittedX.each_row() += mean_p;

  // MSE on observed values
  diff = fittedX - X;
  diff %= not_miss;
  double n_observed = static_cast<double>(X.n_elem - missing.n_elem);
  double mse = arma::accu(diff % diff) / n_observed;

  arma::vec imputed_vals = Xhat.elem(missing);

  return Rcpp::List::create(
      Rcpp::Named("imputed_vals") = imputed_vals,
      Rcpp::Named("mse") = mse);
}
