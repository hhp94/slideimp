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

// [[Rcpp::export]]
int get_max_threads() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}
