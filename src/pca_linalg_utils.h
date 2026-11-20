#ifndef SLIDEIMP_PCA_LINALG_UTILS_H
#define SLIDEIMP_PCA_LINALG_UTILS_H

// Small dense helpers of the PCA chain, over column-major buffers with
// std::size_t dimensions. The two products go to dgemm through
// src/si_lapack.h and return its status; nothing here uses Armadillo or the
// R API.

#include "si_lapack.h"
#include <cmath>
#include <cstddef>

inline constexpr double PCA_TOL = 1e-15;

namespace pca_detail {
// sklearn svd_flip convention: for each of the k columns of the n x k
// matrix V, force the largest-magnitude entry to be non-negative.
inline void canonicalize_signs(double *V, const std::size_t n,
                               const std::size_t k) {
  if (n == 0) return;
  for (std::size_t j = 0; j < k; ++j) {
    double *col = V + j * n;
    std::size_t imax = 0;
    double amax = std::abs(col[0]);
    for (std::size_t i = 1; i < n; ++i) {
      const double a = std::abs(col[i]);
      if (a > amax) {
        amax = a;
        imax = i;
      }
    }
    if (col[imax] < 0.0) {
      for (std::size_t i = 0; i < n; ++i) col[i] = -col[i];
    }
  }
}

// col[i] = (col[i] * rw[i]) * cs[j] -- row scaling then column scaling,
// fused into one pass.
inline void scale_rows_cols_inplace(double *M, const double *rw,
                                    const double *cs, const std::size_t nr,
                                    const std::size_t nc) {
  for (std::size_t j = 0; j < nc; ++j) {
    const double cj = cs[j];
    double *col = M + j * nr;
    for (std::size_t i = 0; i < nr; ++i) {
      col[i] = (col[i] * rw[i]) * cj;
    }
  }
}

inline void copy_scale_rows(double *dst, const double *src, const double *w,
                            const std::size_t nr, const std::size_t nc) {
  for (std::size_t j = 0; j < nc; ++j) {
    double *dcol = dst + j * nr;
    const double *scol = src + j * nr;
    for (std::size_t i = 0; i < nr; ++i) {
      dcol[i] = scol[i] * w[i];
    }
  }
}

// Copy the strict upper triangle of the n x n matrix M onto its lower one.
inline void mirror_upper_to_lower(double *M, const std::size_t n) {
  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t i = 0; i < j; ++i) {
      M[j + i * n] = M[i + j * n];
    }
  }
}

// C = A * B, A a_rows x a_cols, B b_rows x b_cols, C a_rows x b_cols, none
// padded. C must already hold a_rows * b_cols cells. Always dgemm, never
// dgemv, whatever the shape. When a_cols is 0 the product is empty and C
// is zero. Returns the si_la status: err_arg when the inner dimensions
// differ.
[[nodiscard]] inline int gemm_nn(const double *A, const std::size_t a_rows,
                                 const std::size_t a_cols, const double *B,
                                 const std::size_t b_rows,
                                 const std::size_t b_cols, double *C) {
  if (a_cols != b_rows) return si_la::err_arg;
  return si_la::gemm('N', 'N', a_rows, b_cols, a_cols, 1.0, A, a_rows, B,
                     b_rows, 0.0, C, a_rows);
}

// Rows j0 .. j0+nj-1 of C = (columns j0 .. j0+nj-1 of A)' * B, written in
// place: A is a_rows x a_cols, B is a_rows x b_cols, and C is a_cols x
// b_cols, already allocated; no other row of C is touched. A column block
// of a column-major matrix is contiguous with leading dimension a_rows, and
// a row block of C is an offset with leading dimension a_cols, so the
// product is one dgemm whatever the block. When a_rows is 0 the block of C
// is zero. Returns the si_la status: err_arg when B does not have a_rows
// rows or the block runs past A.
[[nodiscard]] inline int
gemm_tn_block(const double *A, const std::size_t a_rows,
              const std::size_t a_cols, const std::size_t j0,
              const std::size_t nj, const double *B, const std::size_t b_rows,
              const std::size_t b_cols, double *C) {
  if (b_rows != a_rows || nj > a_cols || j0 > a_cols - nj)
    return si_la::err_arg;
  return si_la::gemm('T', 'N', nj, b_cols, a_rows, 1.0, A + j0 * a_rows, a_rows,
                     B, b_rows, 0.0, C + j0, a_cols);
}

} // namespace pca_detail

#endif // SLIDEIMP_PCA_LINALG_UTILS_H
