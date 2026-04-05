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
    // shared: the n x n Gram matrix C, upper triangle only
    int n = 0;
    int ldc = 0; // == n
    int nr = 0;  // rows of A (== lda_A)
    int nc = 0;  // cols of A

    // dsyrk: C := A A^T (trans='N') or A^T A (trans='T'), upper triangle
    char trans = 'N';
    int contract = 0; // contraction dim
    int lda_A = 0;    // leading dim of A (always its row count)

    // dsyevr: top-k eigenpairs of C
    int topk = 0;
    int il = 0;
    int iu = 0;
    int lwork = 0;
    int liwork = 0;
    std::vector<double> work;
    std::vector<int> iwork;
    std::vector<int> isuppz;

    // nrows_A, ncols_A: dims of the input matrix A (e.g. Xhat).
    // tall=false -> C = A A^T (n = nrows_A)
    // tall=true -> C = A^T A (n = ncols_A)
    // topk_request: number of top eigenpairs needed (e.g. ncp + 1)
    bool init(arma::uword nrows_A, arma::uword ncols_A, bool tall,
              arma::uword topk_request)
    {
        nr = static_cast<int>(nrows_A);
        nc = static_cast<int>(ncols_A);

        trans = tall ? 'T' : 'N';
        n = tall ? nc : nr;
        contract = tall ? nr : nc;
        lda_A = nr;
        ldc = n;

        topk = static_cast<int>(topk_request);
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
    arma::blas_int n_i = ws.n;
    arma::blas_int k_i = ws.contract;
    arma::blas_int lda = ws.lda_A;
    arma::blas_int ldc = ws.ldc;
    arma::dsyrk_(&PEIG_UPLO, &ws.trans, &n_i, &k_i,
                 &SYRK_ALPHA, A.memptr(), &lda,
                 &SYRK_BETA, C.memptr(), &ldc,
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
    if (static_cast<arma::uword>(ws.n) != A.n_rows)
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
