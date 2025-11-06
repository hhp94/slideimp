// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Col Mins/Max. min if min == 0, max otherwise
// [[Rcpp::export]]
arma::rowvec col_min_max(const arma::mat& mat, const int min = 0)
{
  arma::rowvec result(mat.n_cols);
  if(min == 0) {
    result = arma::min(mat, 0);
  } else {
    result = arma::max(mat, 0);
  }
  return result;
}

// [[Rcpp::export]]
arma::rowvec col_vars(const arma::mat& mat) {
  arma::rowvec result(mat.n_cols);
  for (arma::uword i = 0; i < mat.n_cols; ++i) {
    arma::vec col = mat.col(i);
    arma::vec valid = col.elem(arma::find_finite(col));
    if (valid.n_elem < 1) {
      result(i) = arma::datum::nan;
    } else {
      result(i) = arma::var(valid);
    }
  }
  return result;
}
