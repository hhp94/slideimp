// The Armadillo oracle (dev/eigen_conversion.md, "The Armadillo oracle"):
// each wrapper of src/si_lapack.h beside the Armadillo call it stands in
// for, exported to R as one function per pair. dev/oracle_run.R compiles
// this file and runs every pair; nothing here decides whether a pair agrees.
//
// dev/oracle_run.R stages this file in one folder with src/si_lapack.h and
// src/si_lapack.cpp and compiles it with Rcpp::sourceCpp(). sourceCpp()
// builds si_lapack.cpp as a translation unit of its own, because the header
// of the same name is included here, so the wrappers are compiled the way
// the package compiles them: R's BLAS and LAPACK headers and no Armadillo.
// This unit holds Armadillo and never includes R_ext/Lapack.h.
//
// Every export returns list(wrapper = , arma = ). The wrapper side holds
// `status` and the values; the Armadillo side holds `failed` and the values.
//
// [[Rcpp::depends(RcppArmadillo)]]

// As the package is built (src/Makevars).
#define ARMA_64BIT_WORD 1

// Armadillo writes its warnings to ARMA_CERR_STREAM. Pointing that at a
// stream of our own lets a pair report whether Armadillo warned, which is
// how it says it fell back to an approximate solution.
#include <sstream>
#include <string>
namespace oracle
{
  inline std::ostringstream arma_log;
}
#define ARMA_CERR_STREAM oracle::arma_log

#include <RcppArmadillo.h>
#include "si_lapack.h"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <limits>
#include <vector>

namespace
{

  std::string take_arma_log()
  {
    std::string s = oracle::arma_log.str();
    oracle::arma_log.str(std::string());
    oracle::arma_log.clear();
    return s;
  }

  // BLAS wants a leading dimension of at least 1 even for an operand with no
  // rows; Armadillo's own calls never pass such an operand.
  arma::blas_int ld(const arma::uword n_rows)
  {
    return static_cast<arma::blas_int>(std::max<arma::uword>(1, n_rows));
  }

} // namespace

// [[Rcpp::export]]
Rcpp::List oracle_info()
{
  return Rcpp::List::create(
      Rcpp::Named("armadillo") = arma::arma_version::as_string(),
      Rcpp::Named("warn_level") = static_cast<int>(ARMA_WARN_LEVEL),
      Rcpp::Named("uword_bytes") = static_cast<int>(sizeof(arma::uword)),
      Rcpp::Named("blas_int_bytes") = static_cast<int>(sizeof(arma::blas_int)));
}

// si_la::gemm against arma::dgemm_, the declaration the PCA chain calls
// today. C0 is the content of C before the call.
// [[Rcpp::export]]
Rcpp::List pair_gemm(std::string trans_a, std::string trans_b,
                     arma::mat A, arma::mat B,
                     double alpha, double beta, arma::mat C0)
{
  const char ta = trans_a.at(0);
  const char tb = trans_b.at(0);
  const arma::uword m = (ta == 'N') ? A.n_rows : A.n_cols;
  const arma::uword k = (ta == 'N') ? A.n_cols : A.n_rows;
  const arma::uword n = (tb == 'N') ? B.n_cols : B.n_rows;
  if (C0.n_rows != m || C0.n_cols != n)
  {
    Rcpp::stop("pair_gemm: C0 is not op(A) rows x op(B) cols");
  }

  arma::mat Cw = C0;
  const int status = si_la::gemm(ta, tb, m, n, k, alpha,
                                 A.memptr(), A.n_rows, B.memptr(), B.n_rows,
                                 beta, Cw.memptr(), Cw.n_rows);

  arma::mat Ca = C0;
  const arma::blas_int mi = static_cast<arma::blas_int>(m);
  const arma::blas_int ni = static_cast<arma::blas_int>(n);
  const arma::blas_int ki = static_cast<arma::blas_int>(k);
  const arma::blas_int lda = ld(A.n_rows);
  const arma::blas_int ldb = ld(B.n_rows);
  const arma::blas_int ldc = ld(Ca.n_rows);
  arma::dgemm_(&ta, &tb, &mi, &ni, &ki, &alpha, A.memptr(), &lda,
               B.memptr(), &ldb, &beta, Ca.memptr(), &ldc, 1, 1);

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = Cw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = false, Rcpp::Named("value") = Ca));
}

// si_la::gemm on a block inside a larger buffer, against the Armadillo
// expression for the same block. The two forms are the ones the chain uses:
//   "tn_rows": C.rows(j0, j0 + nj - 1) = A.cols(j0, j0 + nj - 1).t() * B
//   "nt_cols": C.cols(j0, j0 + nj - 1) = A * B.rows(j0, j0 + nj - 1).t()
// [[Rcpp::export]]
Rcpp::List pair_gemm_block(std::string kind, arma::mat A, arma::mat B,
                           int j0, int nj, arma::mat C0)
{
  const arma::uword first = static_cast<arma::uword>(j0);
  const arma::uword count = static_cast<arma::uword>(nj);
  const arma::uword last = first + count - 1;
  arma::mat Cw = C0;
  arma::mat Ca = C0;
  int status = 0;
  if (kind == "tn_rows")
  {
    status = si_la::gemm('T', 'N', count, B.n_cols, A.n_rows, 1.0,
                         A.colptr(first), A.n_rows, B.memptr(), B.n_rows,
                         0.0, Cw.memptr() + first, Cw.n_rows);
    Ca.rows(first, last) = A.cols(first, last).t() * B;
  }
  else if (kind == "nt_cols")
  {
    status = si_la::gemm('N', 'T', A.n_rows, count, A.n_cols, 1.0,
                         A.memptr(), A.n_rows, B.memptr() + first, B.n_rows,
                         0.0, Cw.colptr(first), Cw.n_rows);
    Ca.cols(first, last) = A * B.rows(first, last).t();
  }
  else
  {
    Rcpp::stop("pair_gemm_block: unknown kind");
  }
  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = Cw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = false, Rcpp::Named("value") = Ca));
}

// si_la::gemv against arma::dgemv_. With an empty x the wrapper returns
// beta * y, which dgemv does not do, so that case is compared with the
// Armadillo product instead and needs alpha = 1, beta = 0.
// [[Rcpp::export]]
Rcpp::List pair_gemv(std::string trans, arma::mat A, arma::vec x,
                     double alpha, double beta, arma::vec y0)
{
  const char t = trans.at(0);
  const arma::uword len_x = (t == 'N') ? A.n_cols : A.n_rows;
  const arma::uword len_y = (t == 'N') ? A.n_rows : A.n_cols;
  if (x.n_elem != len_x || y0.n_elem != len_y)
  {
    Rcpp::stop("pair_gemv: x or y0 has the wrong length");
  }

  arma::vec yw = y0;
  const int status = si_la::gemv(t, A.n_rows, A.n_cols, alpha, A.memptr(),
                                 A.n_rows, x.memptr(), beta, yw.memptr());

  arma::vec ya = y0;
  if (len_x == 0)
  {
    if (alpha != 1.0 || beta != 0.0)
    {
      Rcpp::stop("pair_gemv: an empty x is compared at alpha = 1, beta = 0");
    }
    if (t == 'N')
    {
      ya = A * x;
    }
    else
    {
      ya = A.t() * x;
    }
  }
  else
  {
    const arma::blas_int mi = static_cast<arma::blas_int>(A.n_rows);
    const arma::blas_int ni = static_cast<arma::blas_int>(A.n_cols);
    const arma::blas_int lda = ld(A.n_rows);
    const arma::blas_int inc = 1;
    arma::dgemv_(&t, &mi, &ni, &alpha, A.memptr(), &lda, x.memptr(), &inc,
                 &beta, ya.memptr(), &inc, 1);
  }

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = yw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = false, Rcpp::Named("value") = ya));
}

// si_la::syrk against arma::dsyrk_. The whole of C is returned, so the
// triangle the call must leave alone is compared too.
// [[Rcpp::export]]
Rcpp::List pair_syrk(std::string uplo, std::string trans, arma::mat A,
                     double alpha, double beta, arma::mat C0)
{
  const char ul = uplo.at(0);
  const char t = trans.at(0);
  const arma::uword n = (t == 'N') ? A.n_rows : A.n_cols;
  const arma::uword k = (t == 'N') ? A.n_cols : A.n_rows;
  if (C0.n_rows != n || C0.n_cols != n)
  {
    Rcpp::stop("pair_syrk: C0 is not n x n");
  }

  arma::mat Cw = C0;
  const int status = si_la::syrk(ul, t, n, k, alpha, A.memptr(), A.n_rows,
                                 beta, Cw.memptr(), Cw.n_rows);

  arma::mat Ca = C0;
  const arma::blas_int ni = static_cast<arma::blas_int>(n);
  const arma::blas_int ki = static_cast<arma::blas_int>(k);
  const arma::blas_int lda = ld(A.n_rows);
  const arma::blas_int ldc = ld(Ca.n_rows);
  arma::dsyrk_(&ul, &t, &ni, &ki, &alpha, A.memptr(), &lda, &beta,
               Ca.memptr(), &ldc, 1, 1);

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = Cw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = false, Rcpp::Named("value") = Ca));
}

// si_la::syrk writing a diagonal block of a larger Gram, as the cached tall
// Gram does: the upper triangle of G(j0.., j0..) = X.cols(j0..)' X.cols(j0..).
// Against the Armadillo product of the same block; only the upper triangle
// of the block is returned from either side.
// [[Rcpp::export]]
Rcpp::List pair_syrk_block(arma::mat X, int j0)
{
  const arma::uword first = static_cast<arma::uword>(j0);
  const arma::uword nc = X.n_cols - first;
  arma::mat G(X.n_cols, X.n_cols, arma::fill::zeros);
  const int status = si_la::syrk('U', 'T', nc, X.n_rows, 1.0,
                                 X.colptr(first), X.n_rows, 0.0,
                                 G.colptr(first) + first, G.n_rows);
  const arma::mat Gw =
      arma::trimatu(G.submat(first, first, X.n_cols - 1, X.n_cols - 1));

  const arma::mat Xc = X.cols(first, X.n_cols - 1);
  const arma::mat Ga = arma::trimatu(Xc.t() * Xc);

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = Gw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = false, Rcpp::Named("value") = Ga));
}

// si_la::matmul against Armadillo's operator*, in the three forms the LOBPCG
// code writes: "set" is C = op(A) * B, "minus" is C -= op(A) * B, "plus" is
// C += op(A) * B.
// [[Rcpp::export]]
Rcpp::List pair_matmul(bool trans_a, std::string mode, arma::mat A,
                       arma::mat B, arma::mat C0)
{
  double alpha = 1.0;
  double beta = 0.0;
  if (mode == "minus")
  {
    alpha = -1.0;
    beta = 1.0;
  }
  else if (mode == "plus")
  {
    beta = 1.0;
  }
  else if (mode != "set")
  {
    Rcpp::stop("pair_matmul: unknown mode");
  }

  arma::mat Cw = C0;
  const int status = si_la::matmul(trans_a, alpha, A.memptr(), A.n_rows,
                                   A.n_cols, B.memptr(), B.n_rows, B.n_cols,
                                   beta, Cw.memptr());

  arma::mat Ca = C0;
  bool failed = false;
  take_arma_log();
  try
  {
    if (mode == "set")
    {
      if (trans_a)
      {
        Ca = A.t() * B;
      }
      else
      {
        Ca = A * B;
      }
    }
    else if (mode == "minus")
    {
      if (trans_a)
      {
        Ca -= A.t() * B;
      }
      else
      {
        Ca -= A * B;
      }
    }
    else
    {
      if (trans_a)
      {
        Ca += A.t() * B;
      }
      else
      {
        Ca += A * B;
      }
    }
  }
  catch (const std::exception &)
  {
    failed = true;
  }
  const std::string log = take_arma_log();

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = Cw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = failed, Rcpp::Named("value") = Ca,
          Rcpp::Named("log") = log));
}

// si_la::syevr_range against arma::eig_sym(). Armadillo has no dsyevr call:
// eig_sym() is dsyevd, a different algorithm, so this pair agrees to a
// tolerance and not to the bit. Only the uplo triangle of A is meaningful;
// the reference is built from that triangle alone.
// [[Rcpp::export]]
Rcpp::List pair_syevr(std::string uplo, arma::mat A, int il, int iu)
{
  const char ul = uplo.at(0);
  const arma::uword n = A.n_rows;
  const arma::uword count = static_cast<arma::uword>(iu - il + 1);

  arma::mat Aw = A;
  arma::vec w(count);
  arma::mat Z(n, count);
  const int status =
      si_la::syevr_range(ul, n, static_cast<std::size_t>(il),
                         static_cast<std::size_t>(iu), Aw.memptr(),
                         w.memptr(), Z.memptr());

  const arma::mat S = (ul == 'U') ? arma::mat(arma::symmatu(A))
                                  : arma::mat(arma::symmatl(A));
  arma::vec all_w;
  arma::mat all_Z;
  take_arma_log();
  const bool ok = arma::eig_sym(all_w, all_Z, S);
  const std::string log = take_arma_log();
  arma::vec wa;
  arma::mat Za;
  if (ok)
  {
    wa = all_w.subvec(static_cast<arma::uword>(il - 1),
                      static_cast<arma::uword>(iu - 1));
    Za = all_Z.cols(static_cast<arma::uword>(il - 1),
                    static_cast<arma::uword>(iu - 1));
  }

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("values") = w,
          Rcpp::Named("vectors") = Z),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = !ok, Rcpp::Named("values") = wa,
          Rcpp::Named("vectors") = Za, Rcpp::Named("sym") = S,
          Rcpp::Named("log") = log));
}

// si_la::syevd against arma::dsyevd_, the declaration the LOBPCG
// Rayleigh-Ritz calls today: a workspace query, then the solve. `extra`
// lengthens both workspaces past what the query asked for, which the chain
// does when a buffer was sized by an earlier, larger problem.
// [[Rcpp::export]]
Rcpp::List pair_syevd(std::string uplo, arma::mat A, int extra)
{
  const char ul = uplo.at(0);
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);

  arma::mat Aw = A;
  arma::vec ww(A.n_rows);
  const int status = si_la::syevd(ul, A.n_rows, Aw.memptr(), ww.memptr());

  arma::mat Aa = A;
  arma::vec wa(A.n_rows);
  arma::blas_int info = 0;
  if (n > 0)
  {
    const char jobz = 'V';
    arma::blas_int lwork = -1;
    arma::blas_int liwork = -1;
    double work_q = 0.0;
    arma::blas_int iwork_q = 0;
    arma::dsyevd_(&jobz, &ul, &n, Aa.memptr(), &n, wa.memptr(), &work_q,
                  &lwork, &iwork_q, &liwork, &info, 1, 1);
    if (info == 0)
    {
      lwork = static_cast<arma::blas_int>(work_q) + extra;
      liwork = iwork_q + extra;
      std::vector<double> work(static_cast<std::size_t>(lwork));
      std::vector<arma::blas_int> iwork(static_cast<std::size_t>(liwork));
      arma::dsyevd_(&jobz, &ul, &n, Aa.memptr(), &n, wa.memptr(),
                    work.data(), &lwork, iwork.data(), &liwork, &info, 1, 1);
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("values") = ww,
          Rcpp::Named("vectors") = Aw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = (info != 0), Rcpp::Named("values") = wa,
          Rcpp::Named("vectors") = Aa));
}

// si_la::qr_econ against arma::qr_econ().
// [[Rcpp::export]]
Rcpp::List pair_qr_econ(arma::mat X)
{
  const arma::uword p = std::min(X.n_rows, X.n_cols);
  arma::mat Qw(X.n_rows, p);
  arma::mat Rw(p, X.n_cols);
  const int status = si_la::qr_econ(X.n_rows, X.n_cols, X.memptr(),
                                    Qw.memptr(), Rw.memptr());

  arma::mat Qa;
  arma::mat Ra;
  take_arma_log();
  const bool ok = arma::qr_econ(Qa, Ra, X);
  const std::string log = take_arma_log();

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("q") = Qw,
          Rcpp::Named("r") = Rw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = !ok, Rcpp::Named("q") = Qa,
          Rcpp::Named("r") = Ra, Rcpp::Named("log") = log));
}

// si_la::solve_tri against arma::solve() on trimatu() or trimatl(), default
// options. The Armadillo side also reports rcond() of the same triangle and
// whether solve() wrote its "approx solution" warning.
// [[Rcpp::export]]
Rcpp::List pair_solve_tri(std::string uplo, arma::mat A, arma::mat B)
{
  const char ul = uplo.at(0);

  arma::mat Xw = B;
  double rcond_w = 0.0;
  bool approx = false;
  const int status = si_la::solve_tri(ul, A.n_rows, B.n_cols, A.memptr(),
                                      Xw.memptr(), &rcond_w, &approx);

  arma::mat Xa;
  bool failed = false;
  take_arma_log();
  try
  {
    if (ul == 'U')
    {
      Xa = arma::solve(arma::trimatu(A), B);
    }
    else
    {
      Xa = arma::solve(arma::trimatl(A), B);
    }
  }
  catch (const std::exception &)
  {
    failed = true;
  }
  const std::string log = take_arma_log();
  double rcond_a = NA_REAL;
  if (A.n_elem > 0)
  {
    rcond_a = (ul == 'U') ? arma::rcond(arma::trimatu(A))
                          : arma::rcond(arma::trimatl(A));
    take_arma_log();
  }

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = Xw,
          Rcpp::Named("rcond") = rcond_w, Rcpp::Named("approx") = approx),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = failed, Rcpp::Named("value") = Xa,
          Rcpp::Named("rcond") = rcond_a,
          Rcpp::Named("approx") =
              (log.find("approx solution") != std::string::npos),
          Rcpp::Named("log") = log));
}

// The triangular solve as src/lobpcg_warm.h writes it,
//   AQ = solve(trimatl(R.t()), AQ.t()).t()
// against si_la::solve_tri given the same two transposes, written out.
// [[Rcpp::export]]
Rcpp::List pair_solve_chain(arma::mat R, arma::mat AQ)
{
  const arma::mat Rt = R.t();
  arma::mat Bt = AQ.t();
  double rcond_w = 0.0;
  bool approx = false;
  const int status = si_la::solve_tri('L', Rt.n_rows, Bt.n_cols, Rt.memptr(),
                                      Bt.memptr(), &rcond_w, &approx);
  const arma::mat Xw = Bt.t();

  arma::mat Xa;
  bool failed = false;
  take_arma_log();
  try
  {
    Xa = arma::solve(arma::trimatl(R.t()), AQ.t()).t();
  }
  catch (const std::exception &)
  {
    failed = true;
  }
  const std::string log = take_arma_log();

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = Xw,
          Rcpp::Named("rcond") = rcond_w, Rcpp::Named("approx") = approx),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = failed, Rcpp::Named("value") = Xa,
          Rcpp::Named("approx") =
              (log.find("approx solution") != std::string::npos),
          Rcpp::Named("log") = log));
}

// si_la::svd_values against arma::svd(), the form that returns the singular
// values alone.
// [[Rcpp::export]]
Rcpp::List pair_svd(arma::mat X)
{
  arma::mat Xw = X;
  arma::vec sw(std::min(X.n_rows, X.n_cols));
  const int status =
      si_la::svd_values(X.n_rows, X.n_cols, Xw.memptr(), sw.memptr());

  arma::vec sa;
  bool failed = false;
  take_arma_log();
  try
  {
    sa = arma::svd(X);
  }
  catch (const std::exception &)
  {
    failed = true;
  }
  const std::string log = take_arma_log();

  return Rcpp::List::create(
      Rcpp::Named("wrapper") = Rcpp::List::create(
          Rcpp::Named("status") = status, Rcpp::Named("value") = sw),
      Rcpp::Named("arma") = Rcpp::List::create(
          Rcpp::Named("failed") = failed, Rcpp::Named("value") = sa,
          Rcpp::Named("log") = log));
}

// What each wrapper returns for a call it must refuse before any Fortran
// routine is entered. No Armadillo side: Armadillo would raise, or let BLAS
// reach xerbla. The oversized calls pass dimensions past int with one-cell
// buffers, which is safe only because the wrapper returns before using them.
// [[Rcpp::export]]
Rcpp::IntegerVector wrapper_refusals()
{
  const std::size_t huge =
      static_cast<std::size_t>(std::numeric_limits<int>::max()) + 1;
  double a[4] = {1.0, 0.0, 0.0, 1.0};
  double b[4] = {1.0, 0.0, 0.0, 1.0};
  double c[4] = {0.0, 0.0, 0.0, 0.0};
  double w[2] = {0.0, 0.0};
  double rc = 0.0;
  bool approx = false;

  Rcpp::IntegerVector out = Rcpp::IntegerVector::create(
      Rcpp::Named("gemm_trans") =
          si_la::gemm('X', 'N', 2, 2, 2, 1.0, a, 2, b, 2, 0.0, c, 2),
      Rcpp::Named("gemm_conj") =
          si_la::gemm('N', 'C', 2, 2, 2, 1.0, a, 2, b, 2, 0.0, c, 2),
      Rcpp::Named("gemm_lda") =
          si_la::gemm('N', 'N', 2, 2, 2, 1.0, a, 1, b, 2, 0.0, c, 2),
      Rcpp::Named("gemm_ldb") =
          si_la::gemm('N', 'N', 2, 2, 2, 1.0, a, 2, b, 1, 0.0, c, 2),
      Rcpp::Named("gemm_ldc") =
          si_la::gemm('N', 'N', 2, 2, 2, 1.0, a, 2, b, 2, 0.0, c, 1),
      Rcpp::Named("gemm_size") =
          si_la::gemm('N', 'N', huge, 1, 1, 1.0, a, huge, b, 1, 0.0, c, huge),
      Rcpp::Named("gemv_trans") =
          si_la::gemv('X', 2, 2, 1.0, a, 2, b, 0.0, c),
      Rcpp::Named("gemv_lda") =
          si_la::gemv('N', 2, 2, 1.0, a, 1, b, 0.0, c),
      Rcpp::Named("gemv_size") =
          si_la::gemv('N', huge, 1, 1.0, a, huge, b, 0.0, c),
      Rcpp::Named("syrk_uplo") =
          si_la::syrk('X', 'N', 2, 2, 1.0, a, 2, 0.0, c, 2),
      Rcpp::Named("syrk_trans") =
          si_la::syrk('U', 'C', 2, 2, 1.0, a, 2, 0.0, c, 2),
      Rcpp::Named("syrk_lda") =
          si_la::syrk('U', 'N', 2, 2, 1.0, a, 1, 0.0, c, 2),
      Rcpp::Named("syrk_ldc") =
          si_la::syrk('U', 'N', 2, 2, 1.0, a, 2, 0.0, c, 1),
      Rcpp::Named("syrk_size") =
          si_la::syrk('U', 'N', huge, 1, 1.0, a, huge, 0.0, c, huge),
      Rcpp::Named("matmul_inner") =
          si_la::matmul(false, 1.0, a, 2, 2, b, 1, 2, 0.0, c),
      Rcpp::Named("matmul_size") =
          si_la::matmul(false, 1.0, a, huge, 2, b, 2, 2, 0.0, c),
      Rcpp::Named("syevr_uplo") = si_la::syevr_range('X', 2, 1, 2, a, w, c),
      Rcpp::Named("syevr_il_zero") = si_la::syevr_range('U', 2, 0, 1, a, w, c),
      Rcpp::Named("syevr_il_past_iu") =
          si_la::syevr_range('U', 2, 2, 1, a, w, c),
      Rcpp::Named("syevr_iu_past_n") =
          si_la::syevr_range('U', 2, 1, 3, a, w, c));
  // An Rcpp vector is created from at most 20 named values.
  Rcpp::IntegerVector more = Rcpp::IntegerVector::create(
      Rcpp::Named("syevr_size") = si_la::syevr_range('U', huge, 1, 1, a, w, c),
      Rcpp::Named("syevr_work_size") =
          si_la::syevr_range('U', 100000000, 1, 1, a, w, c),
      Rcpp::Named("syevd_uplo") = si_la::syevd('X', 2, a, w),
      Rcpp::Named("syevd_size") = si_la::syevd('U', huge, a, w),
      Rcpp::Named("syevd_work_size") = si_la::syevd('U', 40000, a, w),
      Rcpp::Named("qr_size") = si_la::qr_econ(huge, 1, a, b, c),
      Rcpp::Named("solve_uplo") =
          si_la::solve_tri('X', 2, 1, a, b, &rc, &approx),
      Rcpp::Named("solve_size") =
          si_la::solve_tri('U', huge, 1, a, b, &rc, &approx),
      Rcpp::Named("svd_size") = si_la::svd_values(huge, 1, a, w));

  Rcpp::IntegerVector all(out.size() + more.size());
  Rcpp::CharacterVector names(out.size() + more.size());
  const Rcpp::CharacterVector out_names = out.names();
  const Rcpp::CharacterVector more_names = more.names();
  for (R_xlen_t i = 0; i < out.size(); ++i)
  {
    all[i] = out[i];
    names[i] = out_names[i];
  }
  for (R_xlen_t i = 0; i < more.size(); ++i)
  {
    all[out.size() + i] = more[i];
    names[out.size() + i] = more_names[i];
  }
  all.names() = names;
  return all;
}

// The status codes of src/si_lapack.h, so the runner names them and does not
// repeat their numbers.
// [[Rcpp::export]]
Rcpp::IntegerVector wrapper_codes()
{
  return Rcpp::IntegerVector::create(
      Rcpp::Named("ok") = si_la::ok,
      Rcpp::Named("err_size") = si_la::err_size,
      Rcpp::Named("err_arg") = si_la::err_arg,
      Rcpp::Named("err_nonfinite") = si_la::err_nonfinite,
      Rcpp::Named("err_count") = si_la::err_count,
      Rcpp::Named("err_alloc") = si_la::err_alloc);
}
