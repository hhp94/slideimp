// R's BLAS and LAPACK for the PCA chain; src/si_lapack.h is the interface.
//
// Writing R Extensions binds this file in four places:
//
// - "Fortran character strings": USE_FC_LEN_T is defined before any R header
//   is included, so R_ext/BLAS.h and R_ext/Lapack.h declare the hidden length
//   of every Fortran character argument, and each such argument is passed its
//   length with FCONE. Every character argument here is a single character,
//   the only case FCONE covers; ilaenv, whose routine name is a longer
//   string, is not called.
// - "Calling C from Fortran and vice versa": every routine is reached through
//   F77_CALL, never by a hand-written trailing underscore.
// - "Numerical analysis subroutines": only routines R's own headers declare
//   are called, and none is declared here.
// - "Portable C and C++ code": a standard header is included for everything
//   used from it, every standard name carries std::, the package's header is
//   included in quotes, and no R header sits inside an extern "C" block.
//
// This file includes nothing of Rcpp, Armadillo or Eigen, and R headers only
// for those declarations: no function here calls the R API.

#ifndef USE_FC_LEN_T
#define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif

#include "si_lapack.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <vector>

namespace
{

  constexpr std::size_t int_max =
      static_cast<std::size_t>(std::numeric_limits<int>::max());

  constexpr double eps = std::numeric_limits<double>::epsilon();

  bool to_int(const std::size_t v, int &out)
  {
    if (v > int_max)
    {
      return false;
    }
    out = static_cast<int>(v);
    return true;
  }

  // Whether a length computed in double fits int. A routine's documented
  // minimum is tested with this before the routine is asked for its optimum:
  // past int, LAPACK's own arithmetic for the answer has already wrapped.
  bool fits_int(const double len)
  {
    return len <= static_cast<double>(int_max);
  }

  // The length of a workspace: what the routine's query reported, raised to
  // the documented minimum. false when that does not fit int. A query that
  // reported nothing usable (0, a negative value, NaN) leaves the minimum.
  bool work_len(const double queried, const double minimum, int &out)
  {
    double len = (minimum > 1.0) ? minimum : 1.0;
    if (queried > len)
    {
      len = queried;
    }
    if (!fits_int(len))
    {
      return false;
    }
    out = static_cast<int>(len);
    return true;
  }

  template <typename T>
  bool alloc(std::vector<T> &v, const std::size_t len)
  {
    try
    {
      v.resize(len);
    }
    catch (const std::exception &)
    {
      return false;
    }
    return true;
  }

  bool all_finite(const double *x, const std::size_t len)
  {
    for (std::size_t i = 0; i < len; ++i)
    {
      if (!std::isfinite(x[i]))
      {
        return false;
      }
    }
    return true;
  }

  // Whether the uplo triangle of the n x n matrix a, diagonal included, holds
  // only finite values. The other triangle is not looked at: LAPACK does not
  // read it, and a caller may leave it unset.
  bool tri_finite(const char uplo, const std::size_t n, const double *a)
  {
    for (std::size_t j = 0; j < n; ++j)
    {
      const std::size_t first = (uplo == 'U') ? 0 : j;
      const std::size_t len = (uplo == 'U') ? j + 1 : n - j;
      if (!all_finite(a + j * n + first, len))
      {
        return false;
      }
    }
    return true;
  }

  bool is_nt(const char c)
  {
    return c == 'N' || c == 'T';
  }

  bool is_ul(const char c)
  {
    return c == 'U' || c == 'L';
  }

  // x = beta * x the way BLAS scales a result: beta = 0 assigns 0, so a NaN
  // or Inf already in x does not survive, and beta = 1 leaves x alone.
  void scale_range(double *x, const std::size_t len, const double beta)
  {
    if (beta == 1.0)
    {
      return;
    }
    if (beta == 0.0)
    {
      std::fill(x, x + len, 0.0);
      return;
    }
    for (std::size_t i = 0; i < len; ++i)
    {
      x[i] *= beta;
    }
  }

  // An upper bound on the levels of dgelsd's divide and conquer for an n x n
  // problem. LAPACK derives the count from SMLSIZ, which it reads through
  // ilaenv; this is the count for SMLSIZ = 1, and a larger SMLSIZ needs no
  // more, so a workspace sized with it is long enough without asking.
  double gelsd_levels(const std::size_t n)
  {
    const double l = std::floor(std::log2(static_cast<double>(n) / 2.0)) + 1.0;
    return (l > 0.0) ? l : 0.0;
  }

} // namespace

namespace si_la
{

  int gemm(const char trans_a, const char trans_b,
           const std::size_t m, const std::size_t n, const std::size_t k,
           const double alpha,
           const double *a, const std::size_t lda,
           const double *b, const std::size_t ldb,
           const double beta,
           double *c, const std::size_t ldc)
  {
    if (!is_nt(trans_a) || !is_nt(trans_b))
    {
      return err_arg;
    }
    if (m == 0 || n == 0)
    {
      return ok;
    }
    if (ldc < m)
    {
      return err_arg;
    }
    if (k == 0)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        scale_range(c + j * ldc, m, beta);
      }
      return ok;
    }
    const std::size_t a_rows = (trans_a == 'N') ? m : k;
    const std::size_t b_rows = (trans_b == 'N') ? k : n;
    if (lda < a_rows || ldb < b_rows)
    {
      return err_arg;
    }
    int m_i, n_i, k_i, lda_i, ldb_i, ldc_i;
    if (!to_int(m, m_i) || !to_int(n, n_i) || !to_int(k, k_i) ||
        !to_int(lda, lda_i) || !to_int(ldb, ldb_i) || !to_int(ldc, ldc_i))
    {
      return err_size;
    }
    F77_CALL(dgemm)(&trans_a, &trans_b, &m_i, &n_i, &k_i, &alpha,
                    a, &lda_i, b, &ldb_i, &beta, c, &ldc_i FCONE FCONE);
    return ok;
  }

  int gemv(const char trans,
           const std::size_t m, const std::size_t n,
           const double alpha,
           const double *a, const std::size_t lda,
           const double *x,
           const double beta,
           double *y)
  {
    if (!is_nt(trans))
    {
      return err_arg;
    }
    const std::size_t len_x = (trans == 'N') ? n : m;
    const std::size_t len_y = (trans == 'N') ? m : n;
    if (len_y == 0)
    {
      return ok;
    }
    if (len_x == 0)
    {
      scale_range(y, len_y, beta);
      return ok;
    }
    if (lda < m)
    {
      return err_arg;
    }
    int m_i, n_i, lda_i;
    if (!to_int(m, m_i) || !to_int(n, n_i) || !to_int(lda, lda_i))
    {
      return err_size;
    }
    const int inc = 1;
    F77_CALL(dgemv)(&trans, &m_i, &n_i, &alpha, a, &lda_i, x, &inc,
                    &beta, y, &inc FCONE);
    return ok;
  }

  int syrk(const char uplo, const char trans,
           const std::size_t n, const std::size_t k,
           const double alpha,
           const double *a, const std::size_t lda,
           const double beta,
           double *c, const std::size_t ldc)
  {
    if (!is_ul(uplo) || !is_nt(trans))
    {
      return err_arg;
    }
    if (n == 0)
    {
      return ok;
    }
    if (ldc < n)
    {
      return err_arg;
    }
    if (k == 0)
    {
      for (std::size_t j = 0; j < n; ++j)
      {
        if (uplo == 'U')
        {
          scale_range(c + j * ldc, j + 1, beta);
        }
        else
        {
          scale_range(c + j * ldc + j, n - j, beta);
        }
      }
      return ok;
    }
    const std::size_t a_rows = (trans == 'N') ? n : k;
    if (lda < a_rows)
    {
      return err_arg;
    }
    int n_i, k_i, lda_i, ldc_i;
    if (!to_int(n, n_i) || !to_int(k, k_i) ||
        !to_int(lda, lda_i) || !to_int(ldc, ldc_i))
    {
      return err_size;
    }
    F77_CALL(dsyrk)(&uplo, &trans, &n_i, &k_i, &alpha, a, &lda_i,
                    &beta, c, &ldc_i FCONE FCONE);
    return ok;
  }

  // The order of the tests is Armadillo's (glue_times::apply): a product
  // with one row is taken before one with one column.
  int matmul(const bool trans_a,
             const double alpha,
             const double *a, const std::size_t a_rows, const std::size_t a_cols,
             const double *b, const std::size_t b_rows, const std::size_t b_cols,
             const double beta,
             double *c)
  {
    const std::size_t c_rows = trans_a ? a_cols : a_rows;
    const std::size_t inner = trans_a ? a_rows : a_cols;
    if (inner != b_rows)
    {
      return err_arg;
    }
    if (c_rows == 0 || b_cols == 0)
    {
      return ok;
    }
    if (inner == 0)
    {
      scale_range(c, c_rows * b_cols, beta);
      return ok;
    }
    if (c_rows == 1)
    {
      // op(A) is one row, held as the vector a: the product is (B' a)'.
      return gemv('T', b_rows, b_cols, alpha, b, b_rows, a, beta, c);
    }
    if (b_cols == 1)
    {
      return gemv(trans_a ? 'T' : 'N', a_rows, a_cols, alpha, a, a_rows, b,
                  beta, c);
    }
    return gemm(trans_a ? 'T' : 'N', 'N', c_rows, b_cols, inner, alpha,
                a, a_rows, b, b_rows, beta, c, c_rows);
  }

  int syevr_range(const char uplo, const std::size_t n,
                  const std::size_t il, const std::size_t iu,
                  double *a, double *w, double *z)
  {
    if (!is_ul(uplo))
    {
      return err_arg;
    }
    if (il < 1 || il > iu || iu > n)
    {
      return err_arg;
    }
    int n_i, il_i, iu_i;
    if (!to_int(n, n_i) || !to_int(il, il_i) || !to_int(iu, iu_i))
    {
      return err_size;
    }
    const std::size_t count = iu - il + 1;

    // dsyevr documents 26 n and 10 n as the least it accepts.
    const double lwork_min = 26.0 * static_cast<double>(n);
    const double liwork_min = 10.0 * static_cast<double>(n);
    if (!fits_int(lwork_min))
    {
      return err_size;
    }

    // dsyevr does not refuse NaN or Inf: asked for every eigenpair it may not
    // return, and otherwise it may report success with NaN in the result.
    if (!tri_finite(uplo, n, a))
    {
      return err_nonfinite;
    }

    // dsyevr documents W as n cells whatever the range. For a proper subset
    // it first stores every eigenvalue its search interval holds and then
    // drops those outside il..iu, so a range that cuts a group of equal
    // eigenvalues writes more than count of them. The caller's w has count
    // cells; the n cells are held here.
    std::vector<int> isuppz;
    std::vector<double> w_all;
    if (!alloc(isuppz, 2 * count) || !alloc(w_all, n))
    {
      return err_alloc;
    }

    const char jobz = 'V';
    const char range = 'I';
    const double vl = 0.0;
    const double vu = 0.0;
    const double abstol = 0.0;
    int m_out = 0;
    int info = 0;

    double work_q = 0.0;
    int iwork_q = 0;
    const int lwork_query = -1;
    const int liwork_query = -1;
    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n_i, a, &n_i, &vl, &vu,
                     &il_i, &iu_i, &abstol, &m_out, w_all.data(), z, &n_i,
                     isuppz.data(), &work_q, &lwork_query,
                     &iwork_q, &liwork_query, &info FCONE FCONE FCONE);
    if (info != 0)
    {
      return info;
    }

    int lwork, liwork;
    if (!work_len(work_q, lwork_min, lwork) ||
        !work_len(static_cast<double>(iwork_q), liwork_min, liwork))
    {
      return err_size;
    }
    std::vector<double> work;
    std::vector<int> iwork;
    if (!alloc(work, static_cast<std::size_t>(lwork)) ||
        !alloc(iwork, static_cast<std::size_t>(liwork)))
    {
      return err_alloc;
    }

    F77_CALL(dsyevr)(&jobz, &range, &uplo, &n_i, a, &n_i, &vl, &vu,
                     &il_i, &iu_i, &abstol, &m_out, w_all.data(), z, &n_i,
                     isuppz.data(), work.data(), &lwork,
                     iwork.data(), &liwork, &info FCONE FCONE FCONE);
    if (info != 0)
    {
      return info;
    }
    if (m_out < 0 || static_cast<std::size_t>(m_out) != count)
    {
      return err_count;
    }
    std::copy(w_all.begin(), w_all.begin() + count, w);
    return ok;
  }

  int syevd(const char uplo, const std::size_t n, double *a, double *w)
  {
    if (!is_ul(uplo))
    {
      return err_arg;
    }
    if (n == 0)
    {
      return ok;
    }
    int n_i;
    if (!to_int(n, n_i))
    {
      return err_size;
    }

    // dsyevd documents these for jobz 'V'. In double, since 2 n^2 leaves int
    // long before n does.
    const double nd = static_cast<double>(n);
    const double lwork_min = 1.0 + 6.0 * nd + 2.0 * nd * nd;
    const double liwork_min = 3.0 + 5.0 * nd;
    if (!fits_int(lwork_min))
    {
      return err_size;
    }

    // dsyevd does not refuse NaN or Inf: it may report success with NaN in
    // the result.
    if (!tri_finite(uplo, n, a))
    {
      return err_nonfinite;
    }

    const char jobz = 'V';
    int info = 0;

    double work_q = 0.0;
    int iwork_q = 0;
    const int lwork_query = -1;
    const int liwork_query = -1;
    F77_CALL(dsyevd)(&jobz, &uplo, &n_i, a, &n_i, w, &work_q, &lwork_query,
                     &iwork_q, &liwork_query, &info FCONE FCONE);
    if (info != 0)
    {
      return info;
    }

    int lwork, liwork;
    if (!work_len(work_q, lwork_min, lwork) ||
        !work_len(static_cast<double>(iwork_q), liwork_min, liwork))
    {
      return err_size;
    }
    std::vector<double> work;
    std::vector<int> iwork;
    if (!alloc(work, static_cast<std::size_t>(lwork)) ||
        !alloc(iwork, static_cast<std::size_t>(liwork)))
    {
      return err_alloc;
    }

    F77_CALL(dsyevd)(&jobz, &uplo, &n_i, a, &n_i, w, work.data(), &lwork,
                     iwork.data(), &liwork, &info FCONE FCONE);
    return info;
  }

  int qr_econ(const std::size_t m, const std::size_t n, const double *x,
              double *q, double *r)
  {
    if (m < n)
    {
      return err_arg;
    }
    if (n == 0)
    {
      return ok;
    }
    int m_i, n_i;
    if (!to_int(m, m_i) || !to_int(n, n_i))
    {
      return err_size;
    }

    // dgeqrf and dorgqr both work in q, which has x's shape.
    std::copy(x, x + m * n, q);

    std::vector<double> tau;
    if (!alloc(tau, n))
    {
      return err_alloc;
    }

    // Each routine is asked for its own workspace. Both document max(1, n)
    // as the least they accept.
    const double lwork_min = static_cast<double>(n);
    const int lwork_query = -1;
    int info = 0;
    double work_q = 0.0;
    int lwork;
    std::vector<double> work;

    F77_CALL(dgeqrf)(&m_i, &n_i, q, &m_i, tau.data(), &work_q,
                     &lwork_query, &info);
    if (info != 0)
    {
      return info;
    }
    if (!work_len(work_q, lwork_min, lwork))
    {
      return err_size;
    }
    if (!alloc(work, static_cast<std::size_t>(lwork)))
    {
      return err_alloc;
    }
    F77_CALL(dgeqrf)(&m_i, &n_i, q, &m_i, tau.data(), work.data(), &lwork,
                     &info);
    if (info != 0)
    {
      return info;
    }

    // r is n x n: the upper triangle dgeqrf left in q, zero below it.
    for (std::size_t j = 0; j < n; ++j)
    {
      std::copy(q + j * m, q + j * m + j + 1, r + j * n);
      std::fill(r + j * n + j + 1, r + j * n + n, 0.0);
    }

    work_q = 0.0;
    F77_CALL(dorgqr)(&m_i, &n_i, &n_i, q, &m_i, tau.data(), &work_q,
                     &lwork_query, &info);
    if (info != 0)
    {
      return info;
    }
    if (!work_len(work_q, lwork_min, lwork))
    {
      return err_size;
    }
    if (!alloc(work, static_cast<std::size_t>(lwork)))
    {
      return err_alloc;
    }
    F77_CALL(dorgqr)(&m_i, &n_i, &n_i, q, &m_i, tau.data(), work.data(),
                     &lwork, &info);
    return info;
  }

  int solve_tri(const char uplo, const std::size_t n, const std::size_t nrhs,
                const double *a, double *b,
                double *rcond, bool *approx)
  {
    *rcond = 0.0;
    *approx = false;
    if (!is_ul(uplo))
    {
      return err_arg;
    }
    if (n == 0 || nrhs == 0)
    {
      return ok;
    }
    int n_i, nrhs_i;
    if (!to_int(n, n_i) || !to_int(nrhs, nrhs_i))
    {
      return err_size;
    }

    // The least-squares fallback starts from the right-hand side as given,
    // not from what dtrtrs left in b.
    std::vector<double> rhs;
    if (!alloc(rhs, n * nrhs))
    {
      return err_alloc;
    }
    std::copy(b, b + n * nrhs, rhs.begin());

    const char trans = 'N';
    const char diag = 'N';
    int info = 0;
    F77_CALL(dtrtrs)(&uplo, &trans, &diag, &n_i, &nrhs_i, a, &n_i, b, &n_i,
                     &info FCONE FCONE FCONE);

    // info > 0 is a zero on the diagonal: nothing was solved, and the
    // estimate stays 0.
    bool solved = (info == 0);
    if (solved)
    {
      std::vector<double> work;
      std::vector<int> iwork;
      if (!alloc(work, 3 * n) || !alloc(iwork, n))
      {
        return err_alloc;
      }
      const char norm = '1';
      double estimate = 0.0;
      F77_CALL(dtrcon)(&norm, &uplo, &diag, &n_i, a, &n_i, &estimate,
                       work.data(), iwork.data(), &info FCONE FCONE FCONE);
      *rcond = (info == 0) ? estimate : 0.0;
      solved = !(*rcond < eps) && !std::isnan(*rcond);
    }
    if (solved)
    {
      return ok;
    }

    // The triangle on its own, the other one zero.
    std::vector<double> tri;
    if (!alloc(tri, n * n))
    {
      return err_alloc;
    }
    for (std::size_t j = 0; j < n; ++j)
    {
      const std::size_t first = (uplo == 'U') ? 0 : j;
      const std::size_t last = (uplo == 'U') ? j + 1 : n;
      std::copy(a + j * n + first, a + j * n + last, tri.data() + j * n + first);
    }
    if (!all_finite(tri.data(), n * n) || !all_finite(rhs.data(), n * nrhs))
    {
      return err_nonfinite;
    }

    // Floors under what dgelsd reports, from its documented minima: SMLSIZ =
    // 25 in the terms that grow with it, and the level bound above in the
    // terms that shrink with it. A library that does not report an integer
    // workspace length leaves iwork_q at 0, and the floor stands.
    const double nd = static_cast<double>(n);
    const double levels = gelsd_levels(n);
    const double lwork_min = 12.0 * nd + 2.0 * nd * 25.0 + 8.0 * nd * levels +
                             nd * static_cast<double>(nrhs) + 26.0 * 26.0;
    const double liwork_min = 3.0 * nd * levels + 11.0 * nd;
    if (!fits_int(lwork_min) || !fits_int(liwork_min))
    {
      return err_size;
    }

    std::vector<double> sv;
    if (!alloc(sv, n))
    {
      return err_alloc;
    }
    int m_i = n_i;
    int lda = n_i;
    int ldb = n_i;
    int rank = 0;
    double rcond_sv = nd * eps;

    double work_q = 0.0;
    int iwork_q = 0;
    int lwork_query = -1;
    F77_CALL(dgelsd)(&m_i, &n_i, &nrhs_i, tri.data(), &lda, rhs.data(), &ldb,
                     sv.data(), &rcond_sv, &rank, &work_q, &lwork_query,
                     &iwork_q, &info);
    if (info != 0)
    {
      return info;
    }

    int lwork, liwork;
    if (!work_len(work_q, lwork_min, lwork) ||
        !work_len(static_cast<double>(iwork_q), liwork_min, liwork))
    {
      return err_size;
    }
    std::vector<double> work;
    std::vector<int> iwork;
    if (!alloc(work, static_cast<std::size_t>(lwork)) ||
        !alloc(iwork, static_cast<std::size_t>(liwork)))
    {
      return err_alloc;
    }

    F77_CALL(dgelsd)(&m_i, &n_i, &nrhs_i, tri.data(), &lda, rhs.data(), &ldb,
                     sv.data(), &rcond_sv, &rank, work.data(), &lwork,
                     iwork.data(), &info);
    if (info != 0)
    {
      return info;
    }
    std::copy(rhs.begin(), rhs.end(), b);
    *approx = true;
    return ok;
  }

#if PCA_IMP_DIAGNOSTICS
  int svd_values(const std::size_t m, const std::size_t n, double *a, double *s)
  {
    const std::size_t p = std::min(m, n);
    if (p == 0)
    {
      return ok;
    }
    int m_i, n_i;
    if (!to_int(m, m_i) || !to_int(n, n_i))
    {
      return err_size;
    }
    if (!all_finite(a, m * n))
    {
      return err_nonfinite;
    }

    // dgesdd documents these for jobz 'N'.
    const double pd = static_cast<double>(p);
    const double lwork_min =
        3.0 * pd + std::max(static_cast<double>(std::max(m, n)), 7.0 * pd);
    if (!fits_int(lwork_min))
    {
      return err_size;
    }
    std::vector<int> iwork;
    if (!alloc(iwork, 8 * p))
    {
      return err_alloc;
    }

    // jobz 'N' reads neither u nor vt, but each needs a leading dimension of
    // at least 1.
    const char jobz = 'N';
    double u = 0.0;
    double vt = 0.0;
    const int ld_one = 1;
    int info = 0;

    double work_q = 0.0;
    const int lwork_query = -1;
    F77_CALL(dgesdd)(&jobz, &m_i, &n_i, a, &m_i, s, &u, &ld_one, &vt,
                     &ld_one, &work_q, &lwork_query, iwork.data(),
                     &info FCONE);
    if (info != 0)
    {
      return info;
    }
    int lwork;
    if (!work_len(work_q, lwork_min, lwork))
    {
      return err_size;
    }
    std::vector<double> work;
    if (!alloc(work, static_cast<std::size_t>(lwork)))
    {
      return err_alloc;
    }

    F77_CALL(dgesdd)(&jobz, &m_i, &n_i, a, &m_i, s, &u, &ld_one, &vt,
                     &ld_one, work.data(), &lwork, iwork.data(), &info FCONE);
    return info;
  }
#endif // PCA_IMP_DIAGNOSTICS

} // namespace si_la
