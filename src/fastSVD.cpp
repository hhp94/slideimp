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

Rcpp::List qrSVD_cpp(const arma::mat &A, arma::uword lim_attempts = 50)
{
  arma::mat U;
  arma::vec d;
  arma::mat V;
  arma::uword n = A.n_rows;
  arma::uword p = A.n_cols;

  bool success = arma::svd_econ(U, d, V, A);
  if (success)
  {
    return Rcpp::List::create(
        Rcpp::Named("u") = U,
        Rcpp::Named("d") = d,
        Rcpp::Named("v") = V);
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
    Rcpp::stop("SVD failed to converge after %d preconditioning attempts", lim_attempts);
  }

  return Rcpp::List::create(
      Rcpp::Named("u") = U,
      Rcpp::Named("d") = d,
      Rcpp::Named("v") = V);
}

// [[Rcpp::export]]
Rcpp::List fastSVD_triplet_cpp(const arma::mat &A, const arma::vec &row_w, arma::uword ncp, double tol = 1e-15)
{
  // mormalize row weights and scale X by multiply each row by corresponding sqrt(row.w) weight
  arma::vec row_w_norm = row_w / arma::sum(row_w);
  arma::vec sqrt_row_w = arma::sqrt(row_w_norm);

  arma::mat X = A.each_col() % sqrt_row_w;

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
  Rcpp::List svdAA = qrSVD_cpp(AA_NxN);
  arma::mat Vn_full = Rcpp::as<arma::mat>(svdAA["u"]);
  arma::vec d_squared = Rcpp::as<arma::vec>(svdAA["d"]);
  arma::vec d_full = arma::sqrt(arma::clamp(d_squared, 0.0, arma::datum::inf));
  arma::uword vs_len = std::min(cols, rows - 1);
  arma::vec vs = d_full.subvec(0, vs_len - 1);
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
  arma::mat U, V;
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
  U.each_col() /= sqrt_row_w;

  return Rcpp::List::create(
      Rcpp::Named("vs") = vs,
      Rcpp::Named("U") = U,
      Rcpp::Named("V") = V);
}
