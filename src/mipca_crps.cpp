#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
Rcpp::NumericVector crps_block(Rcpp::NumericMatrix d,
                               Rcpp::NumericVector truth) {
  const int N = d.nrow();
  const int M = d.ncol();
  if (M < 2) {
    Rcpp::stop("'d' must have at least two columns of draws");
  }
  if (truth.size() != N) {
    Rcpp::stop("'truth' must have one value per row of 'd'");
  }
  if (N == 0) {
    return Rcpp::NumericVector(0);
  }

  const double *dp = d.begin();
  const double *yp = truth.begin();

  for (int i = 0; i < N; ++i) {
    if (!R_finite(yp[i])) Rcpp::stop("'truth' contains non-finite values");
  }

  std::vector<double> ds(static_cast<size_t>(N) * M);
  std::vector<double> buf(M);
  std::vector<double> mad(N);
  for (int i = 0; i < N; ++i) {
    double acc = 0.0;
    const double yi = yp[i];
    for (int j = 0; j < M; ++j) {
      const double v = dp[static_cast<R_xlen_t>(j) * N + i];
      if (!R_finite(v)) Rcpp::stop("'d' contains non-finite values");
      buf[j] = v;
      acc += std::fabs(v - yi);
    }
    mad[i] = acc / static_cast<double>(M);
    std::sort(buf.begin(), buf.end());
    for (int j = 0; j < M; ++j) ds[static_cast<R_xlen_t>(j) * N + i] = buf[j];
  }

  std::vector<double> w(M);
  for (int j = 0; j < M; ++j) {
    w[j] = 2.0 * (j + 1) - M - 1.0;
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

  const double mm1 = static_cast<double>(M) * static_cast<double>(M - 1);
  Rcpp::NumericVector out(N);
  for (int i = 0; i < N; ++i) out[i] = mad[i] - g[i] / mm1;
  return out;
}
