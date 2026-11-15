// mipca_crps.cpp -- crps_block(), the fair ensemble CRPS behind crps_all()
// (Ferro 2014, fair scores).
//
// BITWISE against this R form:
//   ds <- t(apply(d, 1, sort)); w <- 2 * seq_len(M) - M - 1
//   rowMeans(abs(d - truth)) - as.vector(ds %*% w) / (M * (M - 1))
// Each piece keeps R's mechanism: a sort (a permutation), a long-double mean
// in ascending column order as rowMeans does, and R's own dgemv.

#include <Rcpp.h>
#include <R_ext/BLAS.h>

#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector crps_block(Rcpp::NumericMatrix d,
                                  Rcpp::NumericVector truth) {
  const int N = d.nrow();
  const int M = d.ncol();
  if (M < 2) {
    // The fair CRPS spread term is undefined at one draw.
    Rcpp::stop("`d` must have at least two columns of draws.");
  }
  if (truth.size() != N) {
    Rcpp::stop("`truth` must have one value per row of `d`.");
  }

  const double* dp = d.begin();
  const double* yp = truth.begin();

  for (int i = 0; i < N; ++i) {
    if (!R_finite(yp[i])) Rcpp::stop("`truth` contains non-finite values.");
  }

  // One gathered pass per row: validate, mean, sort. The accumulation must stay
  // in ascending j, rowMeans' order. +0.0 and -0.0 tie but differ in bytes, the
  // one input where a sort may not match R's bitwise.
  std::vector<double> ds(static_cast<size_t>(N) * M);
  std::vector<double> buf(M);
  std::vector<double> mad(N);
  for (int i = 0; i < N; ++i) {
    long double acc = 0.0L;
    const double yi = yp[i];
    for (int j = 0; j < M; ++j) {
      const double v = dp[static_cast<R_xlen_t>(j) * N + i];
      // Before the sort: std::sort over a NaN is undefined behaviour.
      if (!R_finite(v)) Rcpp::stop("`d` contains non-finite values.");
      buf[j] = v;
      acc += std::fabs(v - yi);
    }
    mad[i] = static_cast<double>(acc / static_cast<long double>(M));
    std::sort(buf.begin(), buf.end());
    for (int j = 0; j < M; ++j) ds[static_cast<R_xlen_t>(j) * N + i] = buf[j];
  }

  // Integer arithmetic until the - 1, as R forms w; every value is exact.
  std::vector<double> w(M);
  for (int j = 0; j < M; ++j) {
    w[j] = static_cast<double>(2 * (j + 1) - M) - 1.0;
  }

  std::vector<double> g(N);
  {
    const char trans = 'N';
    const double one = 1.0, zero = 0.0;
    const int ione = 1;
    F77_CALL(dgemv)
    (&trans, &N, &M, &one, ds.data(), &N, w.data(), &ione, &zero, g.data(),
     &ione FCONE);
  }

  // A product of doubles, not ints: M * (M - 1) overflows int past M = 46341.
  // Exact either way, and equal to the R form's.
  const double mm1 =
      static_cast<double>(M) * static_cast<double>(M - 1);
  Rcpp::NumericVector out(N);
  for (int i = 0; i < N; ++i) out[i] = mad[i] - g[i] / mm1;
  return out;
}
