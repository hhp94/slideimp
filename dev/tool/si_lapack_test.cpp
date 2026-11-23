// The test of the wrapper layer's own (dev/eigen_conversion.md at 4f2b669,
// "Rulings on the last items"): each function of src/si_lapack.h exported to R
// as it is, with guard cells after every buffer it writes. dev/tool/si_lapack_test.R
// stages this file beside src/si_lapack.h and src/si_lapack.cpp, compiles it
// with Rcpp::sourceCpp() and decides; nothing here decides whether a result is
// right. sourceCpp() builds si_lapack.cpp as a translation unit of its own,
// because the header of the same name is included here.
//
// Every wrapper export takes its output buffers from R, appends guard cells
// holding a sentinel to each, calls the wrapper on the copy, and returns
// `status`, each buffer as the wrapper left it without its guard cells, and
// one `guard_<buffer>` flag per buffer: whether every guard cell still holds
// the sentinel. A NaN written there fails the flag, since NaN equals
// nothing. Sizes come from R as doubles, so a size past int can be asked
// for with a buffer of one cell: the wrapper must refuse it before it reads
// or writes one.
//
// probe_dsyevr_past() is the one export that calls LAPACK itself: dsyevr
// with range 'I', as si_la::syevr_range calls it, but with every cell of W
// it is given filled with the sentinel, so it can report how many cells past
// iu - il + 1 dsyevr wrote. It shows that the cases the guard cells after
// syevr_range's eigenvalues are for do make dsyevr write there.

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif

#include "si_lapack.h"

#include <cstddef>
#include <string>
#include <vector>

namespace
{

  constexpr double sentinel = -777.25;

  // Guard cells after every output buffer.
  constexpr std::size_t n_guard = 16;

  // Past this order the extra guard cells syevr_range's buffers get are not
  // added: those are the size cases, which hand it a buffer of one cell.
  constexpr std::size_t big_order = 100000;

  std::size_t as_size(const double d)
  {
    if (!(d >= 0.0))
    {
      Rcpp::stop("a size must be a nonnegative number");
    }
    return static_cast<std::size_t>(d);
  }

  char as_char(const std::string &s)
  {
    if (s.size() != 1)
    {
      Rcpp::stop("an option must be one character");
    }
    return s[0];
  }

  // The cells of x followed by `guard` cells holding the sentinel.
  std::vector<double> guarded(const Rcpp::NumericVector &x,
                              const std::size_t guard)
  {
    std::vector<double> v(x.begin(), x.end());
    v.resize(v.size() + guard, sentinel);
    return v;
  }

  bool guard_intact(const std::vector<double> &v, const std::size_t len)
  {
    for (std::size_t i = len; i < v.size(); ++i)
    {
      if (!(v[i] == sentinel))
      {
        return false;
      }
    }
    return true;
  }

  Rcpp::NumericVector body(const std::vector<double> &v,
                           const std::size_t len)
  {
    return Rcpp::NumericVector(v.begin(), v.begin() + len);
  }

  // A buffer of what R passed, with its guard cells, and the length that was
  // passed.
  struct Out
  {
    std::size_t len;
    std::vector<double> v;
    Out(const Rcpp::NumericVector &x, const std::size_t guard)
        : len(static_cast<std::size_t>(x.size())), v(guarded(x, guard)) {}
    double *data() { return v.data(); }
    Rcpp::NumericVector cells() const { return body(v, len); }
    bool intact() const { return guard_intact(v, len); }
  };

} // namespace

// [[Rcpp::export]]
Rcpp::List status_codes()
{
  return Rcpp::List::create(
      Rcpp::Named("ok") = si_la::ok,
      Rcpp::Named("err_size") = si_la::err_size,
      Rcpp::Named("err_arg") = si_la::err_arg,
      Rcpp::Named("err_nonfinite") = si_la::err_nonfinite,
      Rcpp::Named("err_count") = si_la::err_count,
      Rcpp::Named("err_alloc") = si_la::err_alloc);
}

// PCA_IMP_DIAGNOSTICS as this unit sees it. The wrapper unit is compiled
// with the same flags; if it were not, svd_values would not link.
// [[Rcpp::export]]
int diagnostics_built()
{
  return PCA_IMP_DIAGNOSTICS;
}

// [[Rcpp::export]]
Rcpp::List t_gemm(const std::string &ta, const std::string &tb,
                  const double m, const double n, const double k,
                  const double alpha,
                  const Rcpp::NumericVector &a, const double lda,
                  const Rcpp::NumericVector &b, const double ldb,
                  const double beta,
                  const Rcpp::NumericVector &c, const double ldc)
{
  Out c_out(c, n_guard);
  const int status = si_la::gemm(
      as_char(ta), as_char(tb), as_size(m), as_size(n), as_size(k), alpha,
      a.begin(), as_size(lda), b.begin(), as_size(ldb), beta,
      c_out.data(), as_size(ldc));
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("c") = c_out.cells(),
                            Rcpp::Named("guard_c") = c_out.intact());
}

// [[Rcpp::export]]
Rcpp::List t_gemv(const std::string &trans, const double m, const double n,
                  const double alpha,
                  const Rcpp::NumericVector &a, const double lda,
                  const Rcpp::NumericVector &x, const double beta,
                  const Rcpp::NumericVector &y)
{
  Out y_out(y, n_guard);
  const int status = si_la::gemv(
      as_char(trans), as_size(m), as_size(n), alpha, a.begin(), as_size(lda),
      x.begin(), beta, y_out.data());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("y") = y_out.cells(),
                            Rcpp::Named("guard_y") = y_out.intact());
}

// [[Rcpp::export]]
Rcpp::List t_syrk(const std::string &uplo, const std::string &trans,
                  const double n, const double k, const double alpha,
                  const Rcpp::NumericVector &a, const double lda,
                  const double beta,
                  const Rcpp::NumericVector &c, const double ldc)
{
  Out c_out(c, n_guard);
  const int status = si_la::syrk(
      as_char(uplo), as_char(trans), as_size(n), as_size(k), alpha,
      a.begin(), as_size(lda), beta, c_out.data(), as_size(ldc));
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("c") = c_out.cells(),
                            Rcpp::Named("guard_c") = c_out.intact());
}

// w and z get n guard cells more than the others: dsyevr may write every
// eigenvalue of its search interval, up to n of them, before it drops the
// surplus, and a wrapper that handed it the caller's w would write them
// into the guard cells rather than past the buffer.
// [[Rcpp::export]]
Rcpp::List t_syevr_range(const std::string &uplo, const double n,
                         const double il, const double iu,
                         const Rcpp::NumericVector &a,
                         const Rcpp::NumericVector &w,
                         const Rcpp::NumericVector &z)
{
  const std::size_t n_s = as_size(n);
  const std::size_t extra = (n_s <= big_order) ? n_s : 0;
  Out a_out(a, n_guard);
  Out w_out(w, n_guard + extra);
  Out z_out(z, n_guard + extra);
  const int status = si_la::syevr_range(
      as_char(uplo), n_s, as_size(il), as_size(iu), a_out.data(),
      w_out.data(), z_out.data());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("w") = w_out.cells(),
                            Rcpp::Named("z") = z_out.cells(),
                            Rcpp::Named("guard_a") = a_out.intact(),
                            Rcpp::Named("guard_w") = w_out.intact(),
                            Rcpp::Named("guard_z") = z_out.intact());
}

// [[Rcpp::export]]
Rcpp::List t_syevd(const std::string &uplo, const double n,
                   const Rcpp::NumericVector &a,
                   const Rcpp::NumericVector &w)
{
  Out a_out(a, n_guard);
  Out w_out(w, n_guard);
  const int status = si_la::syevd(as_char(uplo), as_size(n), a_out.data(),
                                  w_out.data());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("a") = a_out.cells(),
                            Rcpp::Named("w") = w_out.cells(),
                            Rcpp::Named("guard_a") = a_out.intact(),
                            Rcpp::Named("guard_w") = w_out.intact());
}

// x is copied first, so the R object is never handed to the wrapper, and
// `x` comes back as the wrapper left the copy: it takes x as const.
// [[Rcpp::export]]
Rcpp::List t_qr_econ(const double m, const double n,
                     const Rcpp::NumericVector &x,
                     const Rcpp::NumericVector &q,
                     const Rcpp::NumericVector &r)
{
  Out x_in(x, n_guard);
  Out q_out(q, n_guard);
  Out r_out(r, n_guard);
  const int status = si_la::qr_econ(as_size(m), as_size(n), x_in.data(),
                                    q_out.data(), r_out.data());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("x") = x_in.cells(),
                            Rcpp::Named("q") = q_out.cells(),
                            Rcpp::Named("r") = r_out.cells(),
                            Rcpp::Named("guard_q") = q_out.intact(),
                            Rcpp::Named("guard_r") = r_out.intact());
}

// [[Rcpp::export]]
Rcpp::List t_solve_tri_right(const std::string &uplo, const double m,
                             const double n,
                             const Rcpp::NumericVector &a,
                             const Rcpp::NumericVector &b)
{
  Out b_out(b, n_guard);
  const int status = si_la::solve_tri_right(as_char(uplo), as_size(m),
                                            as_size(n), a.begin(),
                                            b_out.data());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("b") = b_out.cells(),
                            Rcpp::Named("guard_b") = b_out.intact());
}

// [[Rcpp::export]]
Rcpp::List t_svd_values(const double m, const double n,
                        const Rcpp::NumericVector &a,
                        const Rcpp::NumericVector &s)
{
  Out a_out(a, n_guard);
  Out s_out(s, n_guard);
  const int status = si_la::svd_values(as_size(m), as_size(n), a_out.data(),
                                       s_out.data());
  return Rcpp::List::create(Rcpp::Named("status") = status,
                            Rcpp::Named("s") = s_out.cells(),
                            Rcpp::Named("guard_a") = a_out.intact(),
                            Rcpp::Named("guard_s") = s_out.intact());
}

// dsyevr as si_la::syevr_range calls it (jobz 'V', range 'I', abstol 0, the
// queried workspace raised to 26 n and 10 n), with W n cells of the
// sentinel. Returns info and the number of cells past iu - il + 1 that no
// longer hold the sentinel. a is copied, not destroyed.
// [[Rcpp::export]]
Rcpp::List probe_dsyevr_past(const std::string &uplo, const int n,
                             const int il, const int iu,
                             const Rcpp::NumericVector &a)
{
  const char ul = as_char(uplo);
  std::vector<double> a_c(a.begin(), a.end());
  std::vector<double> w(static_cast<std::size_t>(n), sentinel);
  const int count = iu - il + 1;
  std::vector<double> z(static_cast<std::size_t>(n) * count);
  std::vector<int> isuppz(2 * static_cast<std::size_t>(count));
  const char jobz = 'V';
  const char range = 'I';
  const double vl = 0.0;
  const double vu = 0.0;
  const double abstol = 0.0;
  int m_out = 0;
  int info = 0;
  double work_q = 0.0;
  int iwork_q = 0;
  int lwork = -1;
  int liwork = -1;
  F77_CALL(dsyevr)(&jobz, &range, &ul, &n, a_c.data(), &n, &vl, &vu, &il,
                   &iu, &abstol, &m_out, w.data(), z.data(), &n,
                   isuppz.data(), &work_q, &lwork, &iwork_q, &liwork,
                   &info FCONE FCONE FCONE);
  if (info != 0)
  {
    return Rcpp::List::create(Rcpp::Named("info") = info,
                              Rcpp::Named("past") = NA_INTEGER);
  }
  lwork = static_cast<int>(work_q);
  if (lwork < 26 * n)
  {
    lwork = 26 * n;
  }
  liwork = iwork_q;
  if (liwork < 10 * n)
  {
    liwork = 10 * n;
  }
  std::vector<double> work(static_cast<std::size_t>(lwork));
  std::vector<int> iwork(static_cast<std::size_t>(liwork));
  F77_CALL(dsyevr)(&jobz, &range, &ul, &n, a_c.data(), &n, &vl, &vu, &il,
                   &iu, &abstol, &m_out, w.data(), z.data(), &n,
                   isuppz.data(), work.data(), &lwork, iwork.data(), &liwork,
                   &info FCONE FCONE FCONE);
  int past = 0;
  for (int i = count; i < n; ++i)
  {
    if (!(w[static_cast<std::size_t>(i)] == sentinel))
    {
      ++past;
    }
  }
  return Rcpp::List::create(Rcpp::Named("info") = info,
                            Rcpp::Named("past") = past);
}
