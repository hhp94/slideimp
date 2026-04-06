#ifndef PARTIAL_EIG_H
#define PARTIAL_EIG_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

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

#endif // PARTIAL_EIG_H
