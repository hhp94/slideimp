#ifndef SLIDEIMP_SI_LAPACK_H
#define SLIDEIMP_SI_LAPACK_H

// R's BLAS and LAPACK behind plain functions: double pointers, dimensions and
// an int status. Every BLAS and LAPACK call the PCA chain makes goes through
// this header; src/si_lapack.cpp is the one file that includes R_ext/BLAS.h
// and R_ext/Lapack.h for it. This header declares no BLAS or LAPACK symbol
// and uses no Armadillo, Eigen or R type, so any translation unit can
// include it.
//
// Conventions, for every function here:
//
// - Matrices are column-major. A dimension is a std::size_t; one that does
//   not fit the int BLAS and LAPACK take, or a workspace length that does
//   not, is refused with err_size, never narrowed.
// - The return value is a status. 0 is success; anything else is a failure,
//   and what the output buffers then hold is unspecified. A value LAPACK
//   returned in info is passed on as it is (positive: the routine's own
//   failure). The wrapper's own refusals are the err_* codes below, far from
//   anything LAPACK returns.
// - Nothing here calls the R API, prints, or throws: raising the R error is
//   the caller's. A failed allocation is err_alloc.
// - Arguments are checked before the Fortran routine is entered, so an
//   illegal one is err_arg here. It never reaches xerbla, which under R
//   raises an R error from inside the library and unwinds no C++ frame.
// - Input and output buffers must not overlap unless a function says so.
//
// dev/pca_imp_structure.md, "What the wrappers reproduce", records the
// Armadillo call each of these stands in for, as read in armadillo_bits.
// dev/oracle_run.R checks each against that call.

#include <cstddef>

namespace si_la
{

  inline constexpr int ok = 0;
  // A size or a workspace length does not fit int.
  inline constexpr int err_size = -1001;
  // An option character, a leading dimension or a pair of dimensions the
  // routine would refuse.
  inline constexpr int err_arg = -1002;
  // The input holds NaN or Inf, where the Armadillo call refused it.
  inline constexpr int err_nonfinite = -1003;
  // LAPACK returned without error but with fewer eigenpairs than asked for.
  inline constexpr int err_count = -1004;
  // Workspace could not be allocated.
  inline constexpr int err_alloc = -1005;

  // dgemm: C = alpha * op(A) * op(B) + beta * C, where op(X) is X for 'N' and
  // its transpose for 'T'. op(A) is m x k, op(B) is k x n, C is m x n; lda,
  // ldb and ldc are the leading dimensions of the buffers. When m or n is 0
  // nothing is written. When k is 0 the product is empty and C becomes
  // beta * C, as dgemm defines it.
  int gemm(char trans_a, char trans_b,
           std::size_t m, std::size_t n, std::size_t k,
           double alpha,
           const double *a, std::size_t lda,
           const double *b, std::size_t ldb,
           double beta,
           double *c, std::size_t ldc);

  // dgemv: y = alpha * op(A) * x + beta * y, with A m x n and unit strides.
  // For 'N', x has n elements and y has m; for 'T', x has m and y has n.
  // When x is empty y becomes beta * y; dgemv itself would leave y alone.
  int gemv(char trans,
           std::size_t m, std::size_t n,
           double alpha,
           const double *a, std::size_t lda,
           const double *x,
           double beta,
           double *y);

  // dsyrk: C = alpha * A * A' + beta * C for trans 'N' (A is n x k), or
  // C = alpha * A' * A + beta * C for trans 'T' (A is k x n). C is n x n and
  // only its uplo triangle ('U' or 'L') is read or written. When k is 0 that
  // triangle becomes beta * C.
  int syrk(char uplo, char trans,
           std::size_t n, std::size_t k,
           double alpha,
           const double *a, std::size_t lda,
           double beta,
           double *c, std::size_t ldc);

  // C = alpha * op(A) * B + beta * C for two whole matrices, sent to the BLAS
  // routine Armadillo's operator* sends it to: dgemv when the product is one
  // row or one column, dgemm otherwise. op(A) is A' when trans_a, else A. C
  // is (rows of op(A)) x b_cols, stored without padding. err_arg when the
  // inner dimensions differ.
  //
  // Not as Armadillo: for square operands up to 4 x 4 Armadillo runs its own
  // loops, and for A' * A on one object it calls dsyrk; here both go to dgemm
  // or dgemv.
  int matmul(bool trans_a,
             double alpha,
             const double *a, std::size_t a_rows, std::size_t a_cols,
             const double *b, std::size_t b_rows, std::size_t b_cols,
             double beta,
             double *c);

  // dsyevr, jobz 'V', range 'I', abstol 0: eigenpairs il..iu of the symmetric
  // n x n matrix whose uplo triangle is in a, counted from 1 in ascending
  // order, so il = n - k + 1 and iu = n are the k largest. a is destroyed.
  // w receives iu - il + 1 eigenvalues, ascending; z is n x (iu - il + 1) and
  // receives the eigenvectors in the same order. err_arg unless
  // 1 <= il <= iu <= n.
  int syevr_range(char uplo, std::size_t n, std::size_t il, std::size_t iu,
                  double *a, double *w, double *z);

  // dsyevd, jobz 'V': all eigenpairs of the symmetric n x n matrix whose uplo
  // triangle is in a. On success a holds the eigenvectors and w the n
  // eigenvalues, ascending. n = 0 succeeds and writes nothing.
  int syevd(char uplo, std::size_t n, double *a, double *w);

  // dgeqrf then dorgqr: the economical QR of the m x n matrix x, as
  // arma::qr_econ() returns it. With p = min(m, n), q is m x p with
  // orthonormal columns and r is p x n, upper triangular or trapezoidal, and
  // x = q * r. p = 0 succeeds and writes nothing. The sign of each column is
  // LAPACK's: the diagonal of r is not made positive.
  int qr_econ(std::size_t m, std::size_t n, const double *x,
              double *q, double *r);

  // The solve of arma::solve(arma::trimatu(A), B) and of trimatl, with
  // default options: X with tri(A) * X = B, where tri(A) is the uplo triangle
  // ('U' or 'L') of the n x n matrix a and the other triangle is not read. b
  // is n x nrhs and is overwritten with X.
  //
  // dtrtrs solves and dtrcon estimates the reciprocal condition number in
  // the 1-norm, returned in *rcond. When dtrtrs finds a zero on the diagonal,
  // or *rcond is below DBL_EPSILON or NaN, the triangular answer is discarded
  // and b receives the minimum-norm least-squares solution from dgelsd, with
  // singular values under n * DBL_EPSILON times the largest treated as zero;
  // *approx is then true and the status is still 0. Armadillo prints a
  // warning to its error stream at that point and goes on with that
  // solution; whether to accept it is the caller's decision. err_nonfinite
  // when that fallback is entered with NaN or Inf in tri(A) or B.
  int solve_tri(char uplo, std::size_t n, std::size_t nrhs,
                const double *a, double *b,
                double *rcond, bool *approx);

  // dgesdd, jobz 'N': the min(m, n) singular values of the m x n matrix a,
  // descending, into s. a is destroyed. An empty matrix succeeds and writes
  // nothing. err_nonfinite when a holds NaN or Inf.
  int svd_values(std::size_t m, std::size_t n, double *a, double *s);

} // namespace si_la

#endif // SLIDEIMP_SI_LAPACK_H
