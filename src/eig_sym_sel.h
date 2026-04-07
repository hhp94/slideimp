#ifndef EIG_SYM_SEL
#define EIG_SYM_SEL

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include "loc_timer.h"

extern "C"
{
  void dsyevr_(const char *jobz, const char *range, const char *uplo,
               const int *n, double *a, const int *lda,
               const double *vl, const double *vu,
               const int *il, const int *iu,
               const double *abstol, int *m,
               double *w, double *z, const int *ldz,
               int *isuppz,
               double *work, int *lwork,
               int *iwork, int *liwork,
               int *info);
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

inline constexpr double SYRK_ALPHA = 1.0;
inline constexpr double SYRK_BETA = 0.0;

// ---------------------------------------------------------------------------
// shared workspace, called only once so it's not hot
// ---------------------------------------------------------------------------
struct GramWorkspace
{
  // dsyrk:
  char trans = 'N';
  arma::blas_int nr = 0;  // rows of input A
  arma::blas_int nc = 0;  // cols of input A
  arma::blas_int n = 0;   // == ldc, order of C
  arma::blas_int k = 0;   // contraction dim
  arma::blas_int lda = 0; // leading dim of A (its row count)
  arma::blas_int ldc = 0; // == n
  // dsyevr: top-k eigenpairs of C
  arma::blas_int topk = 0;
  arma::blas_int il = 0;
  arma::blas_int iu = 0;
  arma::blas_int lwork = 0;
  arma::blas_int liwork = 0;
  std::vector<double> work;
  std::vector<int> iwork;
  std::vector<int> isuppz;

  bool init(arma::uword nrows_A, arma::uword ncols_A, bool tall,
            arma::uword topk_request)
  {
    nr = static_cast<arma::blas_int>(nrows_A);
    nc = static_cast<arma::blas_int>(ncols_A);

    trans = tall ? 'T' : 'N';
    n = tall ? nc : nr;
    k = tall ? nr : nc;
    lda = nr;
    ldc = n;

    topk = static_cast<arma::blas_int>(topk_request);
    il = n - topk + 1;
    iu = n;

    isuppz.assign(2 * static_cast<size_t>(topk), 0);

    // dsyevr workspace query
    int m_out = 0;
    int info = 0;
    int lw_flag = -1;
    int liw_flag = -1;
    double work_opt = 0.0;
    int iwork_opt = 0;
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

    int lwork_q = static_cast<int>(work_opt);
    int liwork_q = iwork_opt;

    // floor at documented minima to guard under-reports.
    const int lwork_min = 26 * n;
    const int liwork_min = 10 * n;
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
// LAPACK syrk_upper wrapper, dispatched by ws.trans
// ---------------------------------------------------------------------------
inline void syrk_upper(const arma::mat &A, arma::mat &C,
                       const GramWorkspace &ws)
{
  arma::dsyrk_(&PEIG_UPLO, &ws.trans, &ws.n, &ws.k,
               &SYRK_ALPHA, A.memptr(), &ws.lda,
               &SYRK_BETA, C.memptr(), &ws.ldc,
               1, 1);
}

// ---------------------------------------------------------------------------
// top-k symmetric eigendecomposition via LAPACK dsyevr
// ---------------------------------------------------------------------------
inline bool partial_eig_sym(arma::vec &eigvals,
                            arma::mat &eigvecs,
                            arma::mat &A,
                            GramWorkspace &ws)
{
  if (static_cast<arma::uword>(ws.n) != A.n_rows || A.n_rows != A.n_cols)
  {
    return false;
  }

  int m_out = 0;
  int info = 0;

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

  // dsyevr returns ascending order; flip to descending in place.
  {
    double *ev = eigvals.memptr();
    for (int i = 0, j = ws.topk - 1; i < j; ++i, --j)
    {
      std::swap(ev[i], ev[j]);
    }
  }
  {
    double *Z = eigvecs.memptr();
    const std::ptrdiff_t n_stride = ws.n;
    for (int i = 0, j = ws.topk - 1; i < j; ++i, --j)
    {
      std::swap_ranges(Z + static_cast<std::ptrdiff_t>(i) * n_stride,
                       Z + static_cast<std::ptrdiff_t>(i + 1) * n_stride,
                       Z + static_cast<std::ptrdiff_t>(j) * n_stride);
    }
  }

  return true;
}

constexpr double PCA_TOL = 1e-15;

// helpers to address the heavy overhead in svd copy creating functions like %=
static inline void scale_cols_inplace(double *M, const double *s,
                                      arma::uword nr, arma::uword nc)
{
  for (arma::uword j = 0; j < nc; ++j)
  {
    const double sj = s[j];
    double *col = M + j * nr;
    for (arma::uword i = 0; i < nr; ++i)
    {
      col[i] *= sj;
    }
  }
}

static inline void scale_rows_inplace(double *M, const double *w,
                                      arma::uword nr, arma::uword nc)
{
  for (arma::uword j = 0; j < nc; ++j)
  {
    double *col = M + j * nr;
    for (arma::uword i = 0; i < nr; ++i)
    {
      col[i] *= w[i];
    }
  }
}

static inline void copy_scale_rows(double *dst, const double *src,
                                   const double *w,
                                   arma::uword nr, arma::uword nc)
{
  for (arma::uword j = 0; j < nc; ++j)
  {
    double *dcol = dst + j * nr;
    const double *scol = src + j * nr;
    for (arma::uword i = 0; i < nr; ++i)
    {
      dcol[i] = scol[i] * w[i];
    }
  }
}

// ---------------------------------------------------------------------------
// SVD with symmetric eigs_sym on gram matrix
// ---------------------------------------------------------------------------
inline void SVD_triplet(const arma::mat &Xhat,
                        const double *sw,
                        const double *isw,
                        const arma::uword ncp,
                        const bool tall,
                        arma::vec &vs_top,
                        double &trace_val,
                        arma::mat &U,
                        arma::mat &V,
                        arma::mat &X_work,
                        arma::mat &AA_NxN,
                        arma::vec &d_inv,
                        arma::mat &Vn,
                        arma::vec &eigvals,
                        arma::mat &eigvecs,
                        GramWorkspace &gram_ws
                            LOC_TIMER_PARAM(timer))
{
  const arma::uword nr = static_cast<arma::uword>(gram_ws.nr);
  const arma::uword nc = static_cast<arma::uword>(gram_ws.nc);

    LOC_TIC(timer, "form_gram");
  if (tall)
  {
    X_work.set_size(nr, nc);
    copy_scale_rows(X_work.memptr(), Xhat.memptr(), sw, nr, nc);
    syrk_upper(X_work, AA_NxN, gram_ws);
  }
  else
  {
    syrk_upper(Xhat, AA_NxN, gram_ws);
    double *M = AA_NxN.memptr();
    const arma::uword ng = AA_NxN.n_rows;
    for (arma::uword j = 0; j < ng; ++j)
    {
      const double sj = sw[j];
      double *col = M + j * ng;
      for (arma::uword i = 0; i <= j; ++i)
      {
        col[i] *= sw[i] * sj;
      }
    }
  }
  LOC_TOC(timer, "form_gram");

  trace_val = arma::trace(AA_NxN);

  LOC_TIC(timer, "eig");
  if (!partial_eig_sym(eigvals, eigvecs, AA_NxN, gram_ws))
  {
    Rcpp::stop("`partial_eig_sym` failed to converge");
  }
  LOC_TOC(timer, "eig");

  const arma::uword ke = eigvals.n_elem;
  vs_top.set_size(ke);
  d_inv.set_size(ncp);
  {
    const double *ep = eigvals.memptr();
    double *vp = vs_top.memptr();
    double *dp = d_inv.memptr();
    for (arma::uword i = 0; i < ke; ++i)
    {
      const double e = ep[i];
      vp[i] = (e > 0.0) ? std::sqrt(e) : 0.0;
    }
    for (arma::uword i = 0; i < ncp; ++i)
    {
      dp[i] = (vp[i] > PCA_TOL) ? (1.0 / vp[i]) : 0.0;
    }
  }

  LOC_TIC(timer, "recover_factors");
  if (tall)
  {
    // V = eigvecs[:,0:ncp]; Vn(i,j) = V(i,j) * d_inv(j).
    const arma::uword n = eigvecs.n_rows;
    V.set_size(n, ncp);
    Vn.set_size(n, ncp);
    const double *EV = eigvecs.memptr();
    double *Vp = V.memptr();
    double *Vnp = Vn.memptr();
    const double *dp = d_inv.memptr();
    for (arma::uword j = 0; j < ncp; ++j)
    {
      const double dj = dp[j];
      const double *ecol = EV + j * n;
      double *vcol = Vp + j * n;
      double *ncol = Vnp + j * n;
      for (arma::uword i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        vcol[i] = e;
        ncol[i] = e * dj;
      }
    }
    U = X_work * Vn;
  }
  else
  {
    const arma::uword n = eigvecs.n_rows;
    U.set_size(n, ncp);
    Vn.set_size(n, ncp);
    const double *EV = eigvecs.memptr();
    double *Up = U.memptr();
    double *Vnp = Vn.memptr();
    const double *dp = d_inv.memptr();
    for (arma::uword j = 0; j < ncp; ++j)
    {
      const double dj = dp[j];
      const double *ecol = EV + j * n;
      double *ucol = Up + j * n;
      double *ncol = Vnp + j * n;
      for (arma::uword i = 0; i < n; ++i)
      {
        const double e = ecol[i];
        ucol[i] = e;
        ncol[i] = e * dj * sw[i];
      }
    }
    V = Xhat.t() * Vn;
  }
  // Undo row weighting on U in place.
  scale_rows_inplace(U.memptr(), isw, U.n_rows, U.n_cols);
  LOC_TOC(timer, "recover_factors");
}

#endif // EIG_SYM_SEL
