// Review of the wrapper layer (dev/review_eigen_wrappers.md): direct calls of
// src/si_lapack.h that dev/oracle_pairs.cpp does not make. Plain Rcpp, no
// Armadillo. dev/review_eigen_wrappers_probe.R stages this file in one
// folder with src/si_lapack.h and src/si_lapack.cpp and compiles it with
// Rcpp::sourceCpp(), which builds si_lapack.cpp as a unit of its own.
// Nothing here decides anything; the R script prints what comes back.
//
// The plugin is registered by that script and puts R's LAPACK and BLAS on
// the link line, as src/Makevars does; plain Rcpp does not.
// [[Rcpp::plugins(review_lapack)]]

#include <Rcpp.h>
#include "si_lapack.h"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace
{

  // A value no eigenvalue of the probe's inputs can take.
  constexpr double sentinel = -777.25;

  // Deterministic values in (-1, 1), so the fuzz needs no R RNG.
  struct Lcg
  {
    unsigned long long s = 88172645463325252ULL;
    double next()
    {
      s ^= s << 13;
      s ^= s >> 7;
      s ^= s << 17;
      return static_cast<double>(s % 2000001ULL) / 1000000.0 - 1.0;
    }
    void fill(std::vector<double> &v)
    {
      for (double &x : v)
      {
        x = next();
      }
    }
  };

  std::size_t at_least_one(const std::size_t n)
  {
    return (n > 0) ? n : 1;
  }

} // namespace

// si_la::syevr_range with w given n cells, of which the header promises only
// the first iu - il + 1 are written. The rest start as the sentinel; how
// many of them the call changed is returned.
// [[Rcpp::export]]
Rcpp::List probe_syevr(std::string uplo, Rcpp::NumericMatrix A, int il, int iu)
{
  const std::size_t n = static_cast<std::size_t>(A.nrow());
  const std::size_t count = static_cast<std::size_t>(iu - il + 1);
  std::vector<double> a(A.begin(), A.end());
  std::vector<double> w(n, sentinel);
  Rcpp::NumericMatrix Z(static_cast<int>(n), static_cast<int>(count));
  const int status =
      si_la::syevr_range(uplo.at(0), n, static_cast<std::size_t>(il),
                         static_cast<std::size_t>(iu), a.data(), w.data(),
                         Z.begin());
  int changed = 0;
  for (std::size_t i = count; i < n; ++i)
  {
    if (!(w[i] == sentinel))
    {
      ++changed;
    }
  }
  return Rcpp::List::create(
      Rcpp::Named("status") = status,
      Rcpp::Named("values") =
          Rcpp::NumericVector(w.begin(), w.begin() + count),
      Rcpp::Named("past_count") =
          Rcpp::NumericVector(w.begin() + count, w.end()),
      Rcpp::Named("changed_past_count") = changed,
      Rcpp::Named("vectors") = Z);
}

// [[Rcpp::export]]
Rcpp::List probe_syevd(std::string uplo, Rcpp::NumericMatrix A)
{
  const std::size_t n = static_cast<std::size_t>(A.nrow());
  Rcpp::NumericMatrix V = Rcpp::clone(A);
  Rcpp::NumericVector w(static_cast<int>(n));
  const int status = si_la::syevd(uplo.at(0), n, V.begin(), w.begin());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("values") = w,
                            Rcpp::Named("vectors") = V);
}

// [[Rcpp::export]]
Rcpp::List probe_solve_tri(std::string uplo, Rcpp::NumericMatrix A,
                           Rcpp::NumericMatrix B)
{
  Rcpp::NumericMatrix X = Rcpp::clone(B);
  double rcond = -1.0;
  bool approx = false;
  const int status =
      si_la::solve_tri(uplo.at(0), static_cast<std::size_t>(A.nrow()),
                       static_cast<std::size_t>(B.ncol()), A.begin(),
                       X.begin(), &rcond, &approx);
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("value") = X,
                            Rcpp::Named("rcond") = rcond,
                            Rcpp::Named("approx") = approx);
}

// [[Rcpp::export]]
Rcpp::List probe_svd(Rcpp::NumericMatrix A)
{
  Rcpp::NumericMatrix W = Rcpp::clone(A);
  Rcpp::NumericVector s(std::min(A.nrow(), A.ncol()));
  const int status =
      si_la::svd_values(static_cast<std::size_t>(A.nrow()),
                        static_cast<std::size_t>(A.ncol()), W.begin(),
                        s.begin());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("value") = s);
}

// A 1 x 1 product of a row a (1 x k) and a column b (k x 1), by the three
// routes a product can take: dgemv 'T' on b with the row as the vector
// (what si_la::matmul and Armadillo choose), dgemv 'N' on a with b as the
// vector (the test order swapped), and dgemm.
// [[Rcpp::export]]
Rcpp::NumericVector probe_dot_routes(Rcpp::NumericVector a,
                                     Rcpp::NumericVector b)
{
  const std::size_t k = static_cast<std::size_t>(a.size());
  double by_row = 0.0;
  double by_col = 0.0;
  double by_gemm = 0.0;
  si_la::gemv('T', k, 1, 1.0, b.begin(), k, a.begin(), 0.0, &by_row);
  si_la::gemv('N', 1, k, 1.0, a.begin(), 1, b.begin(), 0.0, &by_col);
  si_la::gemm('N', 'N', 1, 1, k, 1.0, a.begin(), 1, b.begin(), k, 0.0,
              &by_gemm, 1);
  return Rcpp::NumericVector::create(Rcpp::Named("row_first") = by_row,
                                     Rcpp::Named("column_first") = by_col,
                                     Rcpp::Named("dgemm") = by_gemm);
}

// Refusals dev/oracle_pairs.cpp does not try: the leading-dimension checks
// for a transposed operand, the allocation failure, and two boundaries. Each
// returns before a buffer is used, or uses only what the buffers hold.
// [[Rcpp::export]]
Rcpp::IntegerVector probe_refusals()
{
  const std::size_t int_max =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  double a[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double b[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  double c[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double w[3] = {0.0, 0.0, 0.0};
  double rc = 0.0;
  bool approx = false;
  return Rcpp::IntegerVector::create(
      // op(A) = A' with A 3 x 2: lda must be at least 3.
      Rcpp::Named("gemm_T_lda") =
          si_la::gemm('T', 'N', 2, 2, 3, 1.0, a, 2, b, 3, 0.0, c, 2),
      // op(B) = B' with B 3 x 2: ldb must be at least 3.
      Rcpp::Named("gemm_T_ldb") =
          si_la::gemm('N', 'T', 2, 3, 2, 1.0, a, 2, b, 2, 0.0, c, 2),
      Rcpp::Named("gemv_T_lda") =
          si_la::gemv('T', 3, 2, 1.0, a, 2, b, 0.0, c),
      // A' A with A 3 x 2: lda must be at least 3.
      Rcpp::Named("syrk_T_lda") =
          si_la::syrk('U', 'T', 2, 3, 1.0, a, 2, 0.0, c, 2),
      // An empty inner dimension with a leading dimension of 0.
      Rcpp::Named("gemm_k0_ld0") =
          si_la::gemm('N', 'N', 2, 2, 0, 1.0, a, 2, b, 0, 0.5, c, 2),
      Rcpp::Named("solve_alloc") =
          si_la::solve_tri('U', int_max, int_max, a, b, &rc, &approx),
      Rcpp::Named("syevd_n_32767") = si_la::syevd('U', 32767, a, w),
      Rcpp::Named("syevr_n_82595525") =
          si_la::syevr_range('U', 82595525, 1, 1, a, w, c),
      Rcpp::Named("err_arg") = si_la::err_arg,
      Rcpp::Named("err_size") = si_la::err_size,
      Rcpp::Named("err_alloc") = si_la::err_alloc,
      Rcpp::Named("err_nonfinite") = si_la::err_nonfinite);
}

// Every wrapper over every combination of small dimensions, 0 included, with
// leading dimensions passed exactly as the operand's row count (so 0 for an
// operand with no rows). Returns how many calls were made and how many did
// not return 0. An argument that reaches xerbla shows as an R error, or as
// a line the BLAS prints, in the output of the process that runs this.
// [[Rcpp::export]]
Rcpp::IntegerVector probe_fuzz(int max_dim)
{
  const std::size_t top = static_cast<std::size_t>(max_dim);
  const char nt[2] = {'N', 'T'};
  const char ul[2] = {'U', 'L'};
  const double betas[3] = {0.0, 1.0, 0.5};
  Lcg rng;
  int calls = 0;
  int nonzero = 0;
  auto tally = [&](const int status) {
    ++calls;
    if (status != 0)
    {
      ++nonzero;
    }
  };

  for (const char ta : nt)
    for (const char tb : nt)
      for (std::size_t m = 0; m <= top; ++m)
        for (std::size_t n = 0; n <= top; ++n)
          for (std::size_t k = 0; k <= top; ++k)
            for (const double beta : betas)
            {
              const std::size_t ar = (ta == 'N') ? m : k;
              const std::size_t ac = (ta == 'N') ? k : m;
              const std::size_t br = (tb == 'N') ? k : n;
              const std::size_t bc = (tb == 'N') ? n : k;
              std::vector<double> a(at_least_one(ar * ac));
              std::vector<double> b(at_least_one(br * bc));
              std::vector<double> c(at_least_one(m * n));
              rng.fill(a);
              rng.fill(b);
              rng.fill(c);
              tally(si_la::gemm(ta, tb, m, n, k, 1.0, a.data(), ar, b.data(),
                                br, beta, c.data(), m));
            }

  for (const char t : nt)
    for (std::size_t m = 0; m <= top; ++m)
      for (std::size_t n = 0; n <= top; ++n)
        for (const double beta : betas)
        {
          std::vector<double> a(at_least_one(m * n));
          std::vector<double> x(at_least_one((t == 'N') ? n : m));
          std::vector<double> y(at_least_one((t == 'N') ? m : n));
          rng.fill(a);
          rng.fill(x);
          rng.fill(y);
          tally(si_la::gemv(t, m, n, 1.0, a.data(), m, x.data(), beta,
                            y.data()));
        }

  for (const char u : ul)
    for (const char t : nt)
      for (std::size_t n = 0; n <= top; ++n)
        for (std::size_t k = 0; k <= top; ++k)
          for (const double beta : betas)
          {
            const std::size_t ar = (t == 'N') ? n : k;
            const std::size_t ac = (t == 'N') ? k : n;
            std::vector<double> a(at_least_one(ar * ac));
            std::vector<double> c(at_least_one(n * n));
            rng.fill(a);
            rng.fill(c);
            tally(si_la::syrk(u, t, n, k, 1.0, a.data(), ar, beta, c.data(),
                              n));
          }

  for (const bool tr : {false, true})
    for (std::size_t ar = 0; ar <= top; ++ar)
      for (std::size_t ac = 0; ac <= top; ++ac)
        for (std::size_t bc = 0; bc <= top; ++bc)
          for (const double beta : betas)
          {
            const std::size_t br = tr ? ar : ac;
            const std::size_t cr = tr ? ac : ar;
            std::vector<double> a(at_least_one(ar * ac));
            std::vector<double> b(at_least_one(br * bc));
            std::vector<double> c(at_least_one(cr * bc));
            rng.fill(a);
            rng.fill(b);
            rng.fill(c);
            tally(si_la::matmul(tr, -1.0, a.data(), ar, ac, b.data(), br, bc,
                                beta, c.data()));
          }

  for (std::size_t m = 0; m <= top; ++m)
    for (std::size_t n = 0; n <= top; ++n)
    {
      const std::size_t p = std::min(m, n);
      std::vector<double> x(at_least_one(m * n));
      std::vector<double> q(at_least_one(m * p));
      std::vector<double> r(at_least_one(p * n));
      rng.fill(x);
      tally(si_la::qr_econ(m, n, x.data(), q.data(), r.data()));

      std::vector<double> s(at_least_one(p));
      rng.fill(x);
      tally(si_la::svd_values(m, n, x.data(), s.data()));
    }

  for (const char u : ul)
    for (std::size_t n = 0; n <= top; ++n)
    {
      std::vector<double> a(at_least_one(n * n));
      std::vector<double> w(at_least_one(n));
      rng.fill(a);
      tally(si_la::syevd(u, n, a.data(), w.data()));

      for (std::size_t il = 1; il <= n; ++il)
        for (std::size_t iu = il; iu <= n; ++iu)
        {
          // w and z at n cells and n columns, whatever the range.
          std::vector<double> z(at_least_one(n * n));
          rng.fill(a);
          tally(si_la::syevr_range(u, n, il, iu, a.data(), w.data(),
                                   z.data()));
        }

      for (std::size_t nrhs = 0; nrhs <= top; ++nrhs)
        for (const int kind : {0, 1, 2})
        {
          // 0: random triangle; 1: a zero on the diagonal; 2: all zero.
          std::vector<double> b(at_least_one(n * nrhs));
          rng.fill(a);
          rng.fill(b);
          for (std::size_t i = 0; i < n; ++i)
          {
            a[i + i * n] += 3.0;
          }
          if (kind == 1 && n > 0)
          {
            a[(n - 1) + (n - 1) * n] = 0.0;
          }
          if (kind == 2)
          {
            std::fill(a.begin(), a.end(), 0.0);
          }
          double rc = 0.0;
          bool approx = false;
          tally(si_la::solve_tri(u, n, nrhs, a.data(), b.data(), &rc,
                                 &approx));
        }
    }

  return Rcpp::IntegerVector::create(Rcpp::Named("calls") = calls,
                                     Rcpp::Named("nonzero_status") = nonzero);
}
