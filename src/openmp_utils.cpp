#include <Rcpp.h>
using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
bool has_openmp() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}
