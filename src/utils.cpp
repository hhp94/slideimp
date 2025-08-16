// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Col Mins/Max. min if min == 0, max otherwise
// [[Rcpp::export]]
arma::rowvec colMMs(const arma::mat& mat, const int min = 0)
{
  arma::rowvec result;
  if(min == 0) {
    result = arma::min(mat, 0);
  } else {
    result = arma::max(mat, 0);
  }
  return result;
}
