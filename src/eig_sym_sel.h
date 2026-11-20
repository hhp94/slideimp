#ifndef SLIDEIMP_EIG_SYM_SEL_H
#define SLIDEIMP_EIG_SYM_SEL_H

#include <algorithm>
#include <cstddef>
#include <utility>
#include "pca_linalg_utils.h"
#include "si_lapack.h"

// the triangle of the Gram dsyevr reads: the upper one, which
// form_weighted_gram() writes (src/gram_ops.h passes 'U' to each dsyrk).
inline constexpr char PEIG_UPLO = 'U';

// ---------------------------------------------------------------------------
// the eigenpairs asked of dsyevr: the top topk of an n x n symmetric matrix,
// which dsyevr counts from 1 in ascending order as il..n.
// ---------------------------------------------------------------------------
struct EigSymRange {
  std::size_t n = 0;
  std::size_t topk = 0;
  std::size_t il = 0;

  [[nodiscard]] bool init(const std::size_t dim_n,
                          const std::size_t topk_request) {
    if (topk_request == 0 || topk_request > dim_n) {
      return false;
    }
    n = dim_n;
    topk = topk_request;
    il = n - topk + 1;
    return true;
  }
};

// ---------------------------------------------------------------------------
// top-k symmetric eigendecomposition via LAPACK dsyevr, through
// si_la::syevr_range.
//
// A is n x n, and n must be rng.n; only its upper triangle is read.
// eigvals receives rng.topk eigenvalues, descending, and eigvecs, which is
// rng.n x rng.topk, the eigenvectors in the same order, each with its sign
// set by canonicalize_signs(). The status is si_la's: err_arg when n is not
// rng.n, otherwise what si_la::syevr_range returned.
//
// destroys A: dsyevr overwrites the upper triangle it reads, so callers must
// pass a matrix they can discard (the Gram is rebuilt every outer iteration).
// ---------------------------------------------------------------------------
[[nodiscard]] inline int eig_sym_sel(double *eigvals, double *eigvecs,
                                     double *A, const std::size_t n,
                                     const EigSymRange &rng) {
  if (n != rng.n) {
    return si_la::err_arg;
  }

  const int st =
      si_la::syevr_range(PEIG_UPLO, rng.n, rng.il, rng.n, A, eigvals, eigvecs);
  if (st != si_la::ok) {
    return st;
  }

  // flip ascending -> descending in place.
  for (std::size_t i = 0, j = rng.topk - 1; i < j; ++i, --j) {
    std::swap(eigvals[i], eigvals[j]);
    std::swap_ranges(eigvecs + i * rng.n, eigvecs + (i + 1) * rng.n,
                     eigvecs + j * rng.n);
  }

  pca_detail::canonicalize_signs(eigvecs, rng.n, rng.topk);
  return si_la::ok;
}

#endif // SLIDEIMP_EIG_SYM_SEL_H
