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
// dev/pca_imp_structure.md, "What the wrappers do", has what each of these
// does around its routine, and why. dev/si_lapack_test.R checks each against
// R's own functions.

#include <cstddef>

#ifndef PCA_IMP_DIAGNOSTICS
# define PCA_IMP_DIAGNOSTICS 0
#endif

namespace si_la {

inline constexpr int ok = 0;
// A size or a workspace length does not fit int.
inline constexpr int err_size = -1001;
// An option character, a leading dimension or a pair of dimensions the
// routine would refuse.
inline constexpr int err_arg = -1002;
// The input holds NaN or Inf, in a function that says it refuses them.
inline constexpr int err_nonfinite = -1003;
// LAPACK returned without error but with another number of eigenpairs than
// asked for.
inline constexpr int err_count = -1004;
// Workspace could not be allocated.
inline constexpr int err_alloc = -1005;

// dgemm: C = alpha * op(A) * op(B) + beta * C, where op(X) is X for 'N' and
// its transpose for 'T'. op(A) is m x k, op(B) is k x n, C is m x n; lda,
// ldb and ldc are the leading dimensions of the buffers. When m or n is 0
// nothing is written. When k is 0 the product is empty and C becomes
// beta * C, as dgemm defines it.
int gemm(char trans_a, char trans_b, std::size_t m, std::size_t n,
         std::size_t k, double alpha, const double *a, std::size_t lda,
         const double *b, std::size_t ldb, double beta, double *c,
         std::size_t ldc);

// dgemv: y = alpha * op(A) * x + beta * y, with A m x n and unit strides.
// For 'N', x has n elements and y has m; for 'T', x has m and y has n.
// When x is empty y becomes beta * y; dgemv itself would leave y alone.
int gemv(char trans, std::size_t m, std::size_t n, double alpha,
         const double *a, std::size_t lda, const double *x, double beta,
         double *y);

// dsyrk: C = alpha * A * A' + beta * C for trans 'N' (A is n x k), or
// C = alpha * A' * A + beta * C for trans 'T' (A is k x n). C is n x n and
// only its uplo triangle ('U' or 'L') is read or written. When k is 0 that
// triangle becomes beta * C.
int syrk(char uplo, char trans, std::size_t n, std::size_t k, double alpha,
         const double *a, std::size_t lda, double beta, double *c,
         std::size_t ldc);

// dsyevr, jobz 'V', range 'I', abstol 0: eigenpairs il..iu of the symmetric
// n x n matrix whose uplo triangle is in a, counted from 1 in ascending
// order, so il = n - k + 1 and iu = n are the k largest. a is destroyed.
// w receives iu - il + 1 eigenvalues, ascending; z is n x (iu - il + 1) and
// receives the eigenvectors in the same order. Neither needs more cells
// than that: the n cells dsyevr documents for its eigenvalues are held
// inside. err_arg unless 1 <= il <= iu <= n. err_nonfinite when the uplo
// triangle of a holds NaN or Inf, tested before dsyevr is entered; the
// other triangle is not read and may hold anything.
int syevr_range(char uplo, std::size_t n, std::size_t il, std::size_t iu,
                double *a, double *w, double *z);

// dsyevd, jobz 'V': all eigenpairs of the symmetric n x n matrix whose uplo
// triangle is in a. On success a holds the eigenvectors and w the n
// eigenvalues, ascending. n = 0 succeeds and writes nothing. err_nonfinite
// when the uplo triangle of a holds NaN or Inf, tested before dsyevd is
// entered; the other triangle is not read and may hold anything.
int syevd(char uplo, std::size_t n, double *a, double *w);

// dgeqrf then dorgqr: the economical QR of the m x n matrix x, m >= n. q
// is m x n with orthonormal columns, r is n x n and upper triangular, and
// x = q * r. err_arg when m < n. n = 0 succeeds and writes nothing. The
// sign of each column is LAPACK's: the diagonal of r is not made positive.
int qr_econ(std::size_t m, std::size_t n, const double *x, double *q,
            double *r);

// dtrtrs, trans 'N', diag 'N': X with tri(A) * X = B, where tri(A) is the
// uplo triangle ('U' or 'L') of the n x n matrix a and the other triangle
// is not read. b is n x nrhs and is overwritten with X. A positive status
// i is dtrtrs's: diagonal element i of tri(A) is exactly zero, and b is as
// it was given. Nothing estimates the condition of tri(A), and NaN or Inf
// in either operand is not refused: it goes through to X.
int solve_tri(char uplo, std::size_t n, std::size_t nrhs, const double *a,
              double *b);

// dgesdd, jobz 'N': the min(m, n) singular values of the m x n matrix a,
// descending, into s. a is destroyed. An empty matrix succeeds and writes
// nothing. err_nonfinite when a holds NaN or Inf. Only the diagnostics
// build calls it, so only that build compiles it, and the package as
// shipped does not reference dgesdd.
#if PCA_IMP_DIAGNOSTICS
int svd_values(std::size_t m, std::size_t n, double *a, double *s);
#endif

} // namespace si_la

#endif // SLIDEIMP_SI_LAPACK_H
