#ifndef SPECTRA_TOPK_EIG_H
#define SPECTRA_TOPK_EIG_H
//
// Top-k symmetric eigensolver backed by Spectra (implicitly restarted
// Lanczos), used by the experimental pca_imp_spectra() backend.
//
// The Gram matrix built by form_weighted_gram() only has a valid upper
// triangle (dsyrk with uplo = 'U'); DenseSymMatProd<double, Eigen::Upper>
// reads exactly that triangle, so no mirroring pass is needed and A is left
// untouched for a dsyevr fallback.
//

#include <RcppArmadillo.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <algorithm>
#include <random>
#include "pca_linalg_utils.h"

struct SpectraOptions
{
  double tol = 1e-10;
  int maxiter = 1000;
  // Krylov subspace dimension. 0 = auto: min(n, max(2*nev + 1, 20)),
  // matching RSpectra's default.
  int ncv = 0;
};

struct SpectraState
{
  // dominant eigenvector from the previous outer iteration, used as the
  // starting residual for the next Lanczos factorization.
  arma::vec v0;

  bool seeded_for(arma::uword n) const
  {
    return v0.n_elem == n;
  }
};

struct SpectraResult
{
  bool converged = false;
  int nconv = 0;
  int niter = 0;
  int nops = 0;
  double max_rel_res = std::numeric_limits<double>::quiet_NaN();
};

inline bool spectra_topk_eig(arma::vec &eigvals,
                             arma::mat &eigvecs,
                             const arma::mat &A,
                             const arma::uword k,
                             SpectraState &state,
                             const SpectraOptions &opt,
                             SpectraResult &res)
{
  const arma::uword n = A.n_rows;

  // Spectra requires 1 <= nev <= n - 1 and nev < ncv <= n.
  if (k < 1 || k >= n)
  {
    return false;
  }

  const int nev = static_cast<int>(k);
  int ncv = (opt.ncv > 0)
                ? opt.ncv
                : std::min(static_cast<int>(n), std::max(2 * nev + 1, 20));
  ncv = std::max(ncv, nev + 1);
  ncv = std::min(ncv, static_cast<int>(n));

  Eigen::Map<const Eigen::MatrixXd> M(A.memptr(), n, n);
  Spectra::DenseSymMatProd<double, Eigen::Upper> op(M);
  Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
                         Spectra::DenseSymMatProd<double, Eigen::Upper>>
      eigs(&op, nev, ncv);

  if (state.seeded_for(n))
  {
    // Warm-starting Lanczos with the *exact* previous dominant eigenvector is
    // degenerate when the matrix is (nearly) unchanged: the first residual
    // beta collapses toward 0 and the factorization can silently produce
    // spurious Ritz values that still report as converged. A small
    // deterministic perturbation keeps the start rich in the dominant
    // direction while guaranteeing a well-conditioned first Lanczos step.
    arma::vec v0 = state.v0;
    std::mt19937 gen(12345u);
    std::normal_distribution<double> nd(0.0, 1.0);
    const double eps = 1e-3 / std::sqrt(static_cast<double>(n));
    for (arma::uword i = 0; i < n; ++i)
    {
      v0[i] += eps * nd(gen);
    }
    eigs.init(v0.memptr());
  }
  else
  {
    eigs.init();
  }

  const int nconv = static_cast<int>(
      eigs.compute(opt.maxiter, opt.tol, Spectra::LARGEST_ALGE));

  res.nconv = nconv;
  res.niter = static_cast<int>(eigs.num_iterations());
  res.nops = static_cast<int>(eigs.num_operations());
  res.converged = (eigs.info() == Spectra::SUCCESSFUL && nconv >= nev);

  if (!res.converged)
  {
    return false;
  }

  eigvals.set_size(k);
  eigvecs.set_size(n, k);
  Eigen::Map<Eigen::VectorXd>(eigvals.memptr(), nev) = eigs.eigenvalues();
  Eigen::Map<Eigen::MatrixXd>(eigvecs.memptr(), n, nev) = eigs.eigenvectors();

  pca_detail::canonicalize_signs(eigvecs);

  // warm-start the next solve with the current dominant eigenvector.
  state.v0 = eigvecs.col(0);

  return true;
}

// explicit residual check: max_j ||A x_j - lambda_j x_j|| / max(|lambda_1|,
// tiny). Guards against the silent-wrong-convergence mode above. Intended to
// run OUTSIDE any timed region (costs k matvecs). Only A's upper triangle
// need be valid.
inline double spectra_max_rel_res(const arma::mat &A,
                                  const arma::vec &eigvals,
                                  const arma::mat &eigvecs)
{
  const arma::uword k = eigvals.n_elem;
  if (k == 0 || eigvecs.n_cols != k)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const arma::mat As = arma::symmatu(A);
  const double denom = std::max(std::abs(eigvals[0]), PCA_TOL);
  double worst = 0.0;

  for (arma::uword j = 0; j < k; ++j)
  {
    const arma::vec xj = eigvecs.col(j);
    const double r = arma::norm(As * xj - eigvals[j] * xj, 2) / denom;
    worst = std::max(worst, r);
  }

  return worst;
}

#endif // SPECTRA_TOPK_EIG_H
