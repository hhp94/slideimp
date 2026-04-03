#ifndef PARTIAL_EIG_H
#define PARTIAL_EIG_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
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

// -----------------------------------------------------------------------------
// partial_eig_sym
//  Computes the largest k eigenvalues/eigenvectors of a dense symmetric
//  matrix using dsyevr (range='I'). Returns them in descending order.
//
//  Preconditions: A is square, k_request > 0, A.n_rows > 0.
//  WARNING: A is overwritten (upper triangle destroyed).
// -----------------------------------------------------------------------------
inline bool partial_eig_sym(arma::vec &eigvals,
                            arma::mat &eigvecs,
                            arma::mat &A,
                            arma::uword k_request)
{
        const int n = static_cast<int>(A.n_rows);
        const int k = static_cast<int>(std::min(k_request, static_cast<arma::uword>(n)));

        int il = n - k + 1; // 1-based lowest index (ascending order)
        int iu = n;         // 1-based highest index
        int lda = n;
        int ldz = n;

        char jobz = 'V';
        char range = 'I';
        char uplo = 'U';
        double abstol = 0.0;
        double vl = 0.0; // unused but some LAPACK builds dereference
        double vu = 0.0;

        int m_out = 0;
        int info = 0;

        eigvals.set_size(k);
        eigvecs.set_size(n, k);
        std::vector<int> isuppz(2 * k);

        // --- workspace query -------------------------------------------------
        int lwork_q = -1;
        int liwork_q = -1;
        double work_query;
        int iwork_query;

        dsyevr_(&jobz, &range, &uplo, &n,
                A.memptr(), &lda,
                &vl, &vu, &il, &iu,
                &abstol, &m_out,
                eigvals.memptr(), eigvecs.memptr(), &ldz,
                isuppz.data(),
                &work_query, &lwork_q,
                &iwork_query, &liwork_q,
                &info);

        if (info != 0)
        {
                return false;
        }

        int lwork = static_cast<int>(std::ceil(work_query));
        int liwork = iwork_query;
        std::vector<double> work(lwork);
        std::vector<int> iwork(liwork);

        // --- actual computation ----------------------------------------------
        dsyevr_(&jobz, &range, &uplo, &n,
                A.memptr(), &lda,
                &vl, &vu, &il, &iu,
                &abstol, &m_out,
                eigvals.memptr(), eigvecs.memptr(), &ldz,
                isuppz.data(),
                work.data(), &lwork,
                iwork.data(), &liwork,
                &info);

        if (info != 0)
        {
                return false;
        }

        // If dsyevr didn't converge on all requested eigenvalues the trailing
        // entries of eigvals/eigvecs are uninitialised — fail rather than
        // silently return garbage after the reverse below.
        if (m_out != k)
        {
                return false;
        }

        // dsyevr returns ascending order — flip to descending
        eigvals = arma::reverse(eigvals);
        eigvecs = arma::fliplr(eigvecs);

        return true;
}

#endif // PARTIAL_EIG_H
