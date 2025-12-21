#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Codes copied from bootSVD and FactoMiner::svd.triplet. All pre-conditioning done in R
arma::mat genQ_cpp(arma::uword n)
{
  Rcpp::NumericVector rv = Rcpp::rnorm(n * n);
  arma::mat normal_mat(rv.begin(), n, n, false, true);
  arma::mat Q, R;
  arma::qr_econ(Q, R, normal_mat);
  return Q;
}

void qrSVD_cpp(const arma::mat &A, arma::mat &U, arma::vec &d, arma::mat &V, arma::uword lim_attempts = 50)
{
  arma::uword n = A.n_rows;
  arma::uword p = A.n_cols;
  bool success = arma::svd_econ(U, d, V, A);
  if (success)
  {
    return;
  }
  arma::uword attempt = 0;
  while (!success && attempt < lim_attempts)
  {
    attempt++;
    arma::mat Q_n = genQ_cpp(n);
    arma::mat Q_p = (n == p) ? Q_n : genQ_cpp(p);
    arma::mat precond = Q_n.t() * A * Q_p;
    success = arma::svd_econ(U, d, V, precond);
    if (success)
    {
      U = Q_n * U;
      V = Q_p * V;
    }
  }
  if (!success)
  {
    Rcpp::warning("SVD failed to converge after %d preconditioning attempts", lim_attempts);
  }
}

void fastSVD_triplet(const arma::mat &A, const arma::rowvec &row_w, arma::uword ncp,
                     arma::vec &vs, arma::mat &U, arma::mat &V, double tol = 1e-15)
{
  // normalize row weights and scale X by multiplying each row by corresponding sqrt(row_w) weight
  arma::rowvec sqrt_row_w = arma::sqrt(row_w);
  arma::mat X = A.each_col() % sqrt_row_w.t();

  arma::uword rows = X.n_rows;
  arma::uword cols = X.n_cols;
  bool tall = (rows >= cols);

  arma::mat AA_NxN;
  if (!tall)
  {
    AA_NxN = X * X.t();
  }
  else
  {
    AA_NxN = X.t() * X;
  }

  arma::mat Vn_full;
  arma::vec d_squared;
  arma::mat V_unused;
  qrSVD_cpp(AA_NxN, Vn_full, d_squared, V_unused);

  arma::vec d_full = arma::sqrt(arma::clamp(d_squared, 0.0, arma::datum::inf));

  arma::uword vs_len = std::min(cols, rows - 1);
  vs = d_full.subvec(0, vs_len - 1);

  arma::mat Vn = Vn_full.cols(0, ncp - 1);
  arma::vec d = d_full.subvec(0, ncp - 1);

  arma::vec d_inv = arma::zeros<arma::vec>(ncp);
  for (arma::uword i = 0; i < ncp; i++)
  {
    if (d(i) > tol)
    {
      d_inv(i) = 1.0 / d(i);
    }
  }

  arma::mat Vp;
  if (!tall)
  {
    Vp = X.t() * (Vn * arma::diagmat(d_inv));
  }
  else
  {
    Vp = X * (Vn * arma::diagmat(d_inv));
  }

  if (!tall)
  {
    U = Vn;
    V = Vp;
  }
  else
  {
    U = Vp;
    V = Vn;
  }

  // sign flip based on colSums(V)
  arma::rowvec col_sums = arma::sum(V, 0);
  for (arma::uword i = 0; i < ncp; i++)
  {
    if (col_sums(i) < 0)
    {
      U.col(i) *= -1.0;
      V.col(i) *= -1.0;
    }
  }

  U.each_col() /= sqrt_row_w.t();
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
    double coeff_ridge)
{
  double nrX = X.n_rows;
  double ncX = X.n_cols;
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
  arma::mat U;
  arma::mat V;
  arma::vec vs;
  double min_dim = std::min(ncX, nrX - 1);
  double sigma2 = 0.0;
  arma::vec lambda_shrinked;
  // double total_w = arma::sum(row_w); Guaranteed to be 1 in the R wrapper
  while (nb_iter > 0)
  {
    // update and unstandardize Xhat
    Xhat.elem(missing) = fittedX.elem(missing);
    if (scale)
    {
      Xhat.each_row() %= et;
    }
    Xhat.each_row() += mean_p;
    // re standardize Xhat
    // mean_p = (row_w * Xhat) / total_w;
    mean_p = (row_w * Xhat);
    Xhat.each_row() -= mean_p;
    // et = arma::sqrt((row_w * (Xhat % Xhat)) / total_w);
    et = arma::sqrt(row_w * (Xhat % Xhat));
    if (scale)
    {
      Xhat.each_row() /= et;
    }

    fastSVD_triplet(Xhat, row_w, ncp, vs, U, V);
    // sigma2 calculation
    sigma2 = 0.0;
    if (regularized)
    {
      arma::vec vs_tail = vs.subvec(ncp, vs.n_elem - 1);
      arma::vec vs_tail_sq = vs_tail % vs_tail;
      double denom_sigma = (nrX - 1) * ncX - (nrX - 1) * ncp - ncX * ncp + ncp * ncp;
      double scale_factor = (nrX * ncX) / min_dim;
      sigma2 = scale_factor * arma::sum(vs_tail_sq) / denom_sigma;
      sigma2 = std::min(sigma2 * coeff_ridge, vs(ncp) * vs(ncp));
    }
    // lambda_shrinked calculation
    arma::vec vs_head = vs.subvec(0, ncp - 1);
    lambda_shrinked = (vs_head % vs_head - sigma2) / vs_head;

    // fittedX calculation
    // fittedX <- tcrossprod(t(t(U[,1:ncp] * row.w) * lambda.shrinked), V[,1:ncp])
    // U[,1:ncp] * row.w: multiply each row by corresponding weight
    arma::mat U_ncp = U.head_cols(ncp);
    arma::mat U_weighted = U_ncp.each_col() % row_w.t();
    // t(t(...) * lambda.shrinked) sweep * results by corresponding lambda
    arma::mat U_lambda = U_weighted.each_row() % lambda_shrinked.t();
    // tcrossprod with V[,1:ncp]
    arma::mat V_ncp = V.head_cols(ncp);
    fittedX = U_lambda * V_ncp.t();
    // divide by row.w
    fittedX.each_col() /= row_w.t();

    // now we calculate objective
    arma::mat diff = Xhat - fittedX;
    diff.elem(missing).zeros();
    // sum(diff^2 * row.w) - each row weighted by row_w
    objective = arma::dot(row_w, arma::sum(diff % diff, 1));

    // Check convergence
    // if old is Inf, then set criterion to 1 which will not trigger convergence, then just
    // calculate criterion as usual
    double criterion = (old > 0.0 && !std::isinf(old)) ? std::abs(1.0 - objective / old) : 1.0;
    old = objective;
    nb_iter++;
    if ((criterion < threshold) && (nb_iter > miniter))
    {
      nb_iter = 0;
    }
    if ((objective < threshold) && (nb_iter > miniter))
    {
      nb_iter = 0;
    }
    if (nb_iter > maxiter)
    {
      nb_iter = 0;
      Rcpp::warning("Stopped after " + std::to_string(maxiter) + " iterations");
    }
  }
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

  // calculate MSE on observed values which is just the MSE between fittedX and X
  arma::mat diff_obs = fittedX - X;
  diff_obs.elem(missing).zeros(); // ignore missing positions
  double n_observed = static_cast<double>(X.n_elem - missing.n_elem);
  double mse = arma::accu(diff_obs % diff_obs) / n_observed;

  // only extract imputed values instead of returning the double full matrix
  arma::vec imputed_vals = Xhat.elem(missing);

  return Rcpp::List::create(
      Rcpp::Named("imputed_vals") = imputed_vals,
      Rcpp::Named("mse") = mse);
}

