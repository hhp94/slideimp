#ifndef EIG_SYM_SEL_H
#define EIG_SYM_SEL_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cstring>
#include <vector>
#include "pca_linalg_utils.h"

extern "C"
{
  void dsyevr_(const char *jobz, const char *range, const char *uplo,
               const arma::blas_int *n, double *a, const arma::blas_int *lda,
               const double *vl, const double *vu,
               const arma::blas_int *il, const arma::blas_int *iu,
               const double *abstol, arma::blas_int *m,
               double *w, double *z, const arma::blas_int *ldz,
               arma::blas_int *isuppz,
               double *work, arma::blas_int *lwork,
               arma::blas_int *iwork, arma::blas_int *liwork,
               arma::blas_int *info);
}

// ---------------------------------------------------------------------------
// constants
// ---------------------------------------------------------------------------
inline constexpr char PEIG_JOBZ = 'V';
inline constexpr char PEIG_RANGE = 'I';
inline constexpr char PEIG_UPLO = 'U';
inline constexpr double PEIG_ABSTOL = 0.0;
inline constexpr double PEIG_VL = 0.0;
inline constexpr double PEIG_VU = 0.0;

// ---------------------------------------------------------------------------
// dsyevr workspace
// ---------------------------------------------------------------------------
struct EigSymWorkspace
{
  arma::blas_int n = 0;
  arma::blas_int ldc = 0;
  arma::blas_int topk = 0;
  arma::blas_int il = 0;
  arma::blas_int iu = 0;
  arma::blas_int lwork = 0;
  arma::blas_int liwork = 0;
  std::vector<double> work;
  std::vector<arma::blas_int> iwork;
  std::vector<arma::blas_int> isuppz;

  bool init(arma::uword dim_n, arma::uword topk_request)
  {
    n = static_cast<arma::blas_int>(dim_n);
    ldc = n;
    topk = static_cast<arma::blas_int>(topk_request);
    if (topk <= 0 || topk > n)
    {
      return false;
    }
    il = n - topk + 1;
    iu = n;
    isuppz.assign(2 * static_cast<size_t>(topk), 0);

    arma::blas_int m_out = 0;
    arma::blas_int info = 0;
    arma::blas_int lw_flag = -1;
    arma::blas_int liw_flag = -1;
    double work_opt = 0.0;
    arma::blas_int iwork_opt = 0;
    double dummy_a = 0.0;
    double dummy_w = 0.0;
    double dummy_z = 0.0;

    dsyevr_(&PEIG_JOBZ, &PEIG_RANGE, &PEIG_UPLO, &n,
            &dummy_a, &ldc,
            &PEIG_VL, &PEIG_VU, &il, &iu,
            &PEIG_ABSTOL, &m_out,
            &dummy_w, &dummy_z, &ldc,
            isuppz.data(),
            &work_opt, &lw_flag,
            &iwork_opt, &liw_flag,
            &info);

    if (info != 0)
    {
      return false;
    }

    arma::blas_int lwork_q = static_cast<arma::blas_int>(work_opt);
    arma::blas_int liwork_q = iwork_opt;

    const arma::blas_int lwork_min = 26 * n;
    const arma::blas_int liwork_min = 10 * n;
    if (lwork_q < lwork_min)
    {
      lwork_q = lwork_min;
    }
    if (liwork_q < liwork_min)
    {
      liwork_q = liwork_min;
    }
    lwork = lwork_q;
    liwork = liwork_q;
    work.resize(static_cast<size_t>(lwork));
    iwork.resize(static_cast<size_t>(liwork));
    return true;
  }
};

// ---------------------------------------------------------------------------
// top-k symmetric eigendecomposition via LAPACK dsyevr.
// ---------------------------------------------------------------------------
inline bool eig_sym_sel(arma::vec &eigvals,
                        arma::mat &eigvecs,
                        arma::mat &A,
                        EigSymWorkspace &ws)
{
  if (static_cast<arma::uword>(ws.n) != A.n_rows || A.n_rows != A.n_cols)
  {
    return false;
  }

  arma::blas_int m_out = 0;
  arma::blas_int info = 0;

  eigvals.set_size(ws.topk);
  eigvecs.set_size(ws.n, ws.topk);

  dsyevr_(&PEIG_JOBZ, &PEIG_RANGE, &PEIG_UPLO, &ws.n,
          A.memptr(), &ws.ldc,
          &PEIG_VL, &PEIG_VU, &ws.il, &ws.iu,
          &PEIG_ABSTOL, &m_out,
          eigvals.memptr(), eigvecs.memptr(), &ws.ldc,
          ws.isuppz.data(),
          ws.work.data(), &ws.lwork,
          ws.iwork.data(), &ws.liwork,
          &info);

  if (info != 0 || m_out != ws.topk)
  {
    return false;
  }

  // flip ascending -> descending in place.
  {
    double *ev = eigvals.memptr();
    double *Z = eigvecs.memptr();
    const std::ptrdiff_t n_stride = ws.n;
    for (arma::blas_int i = 0, j = ws.topk - 1; i < j; ++i, --j)
    {
      std::swap(ev[i], ev[j]);
      std::swap_ranges(Z + static_cast<std::ptrdiff_t>(i) * n_stride,
                       Z + static_cast<std::ptrdiff_t>(i + 1) * n_stride,
                       Z + static_cast<std::ptrdiff_t>(j) * n_stride);
    }
  }

  pca_detail::canonicalize_signs(eigvecs);
  return true;
}

#endif // EIG_SYM_SEL_H
