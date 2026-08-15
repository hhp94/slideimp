//
// Standalone Spectra backend for pca_imp-style EM PCA imputation (dev bench).
//
// Runs a plain (non-incremental) version of the pca_imp_internal_cpp EM loop
// with the eigensolver swapped to Spectra (implicitly restarted Lanczos, via
// RSpectra's bundled headers): impute -> restandardize -> weighted Gram ->
// Spectra top-k eig -> reconstruct -> objective, using the *same* loss
// formula and convergence criterion as armaSVD.cpp. No fixed-column Gram
// cache, no column permutation, no solver fallback — the gradual, simple
// version. Convergence behavior (n_iter, final objective) should roughly
// match pca_imp; block timings come from rcpptimer (load_all1(timer = TRUE)),
// written to `pca_imp_spectra` next to pca_imp's `pca_imp_gram`.
//
// Imputed values are not returned; this is a diagnostics/timing harness.
//

#include "gram_ops.h"
#include "spectra_topk_eig.h"
#include "svd_triplet.h" // split_eigvec_block, pca_detail gemm helpers
#include "matrix_checks.h"
#include "loc_timer.h"
#include <limits>

// C = A * B.t() into a preallocated C (same dgemm shape as armaSVD.cpp's
// reconstruct_cols over all columns).
static inline void gemm_nt_full(const arma::mat &A,
                                const arma::mat &B,
                                arma::mat &C)
{
  const char trN = 'N';
  const char trT = 'T';
  const double alpha = 1.0;
  const double beta = 0.0;
  const arma::blas_int m = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int n = static_cast<arma::blas_int>(B.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int lda = m;
  const arma::blas_int ldb = static_cast<arma::blas_int>(B.n_rows);
  const arma::blas_int ldc = static_cast<arma::blas_int>(C.n_rows);
  arma::dgemm_(&trN, &trT, &m, &n, &k,
               &alpha,
               A.memptr(), &lda,
               B.memptr(), &ldb,
               &beta,
               C.memptr(), &ldc,
               1, 1);
}

// [[Rcpp::export]]
Rcpp::List pca_imp_spectra_cpp(
    const arma::mat &obj,
    const arma::uvec &eligible_idx,
    const arma::uword ncp,
    const bool scale,
    const bool regularized,
    const double threshold,
    const arma::uword maxiter,
    const arma::uword miniter,
    arma::rowvec row_w,
    const double coeff_ridge,
    const double spectra_tol,
    const arma::uword spectra_maxiter,
    const arma::uword spectra_ncv)
{
  LOC_TIMER_OBJ(pca_imp_spectra);
  LOC_TIC(pca_imp_spectra, "pca_imp_spectra_cpp_total");
  stop_on_inf(obj);

  const arma::uword nrX = obj.n_rows;
  const arma::uword n_elig = eligible_idx.n_elem;

  if (ncp == 0)
  {
    Rcpp::stop("ncp must be >= 1");
  }
  if (maxiter == 0)
  {
    Rcpp::stop("maxiter must be >= 1");
  }
  if (miniter > maxiter)
  {
    Rcpp::stop("miniter must be <= maxiter");
  }
  if (n_elig == 0)
  {
    Rcpp::stop("eligible_idx is empty");
  }
  if (row_w.n_elem != nrX)
  {
    Rcpp::stop("row_w length must equal n_rows(obj)");
  }
  if (eligible_idx.max() >= obj.n_cols)
  {
    Rcpp::stop("eligible_idx contains out-of-bounds column indices");
  }
  if (spectra_maxiter == 0)
  {
    Rcpp::stop("spectra_maxiter must be >= 1");
  }
  if (!row_w.is_finite() || arma::any(row_w < 0.0))
  {
    Rcpp::stop("row_w must be finite and non-negative");
  }
  const double row_w_sum = arma::accu(row_w);
  if (!(row_w_sum > 0.0))
  {
    Rcpp::stop("row_w must have positive sum");
  }
  row_w /= row_w_sum;
  const double *wptr = row_w.memptr();

  const bool tall = (nrX >= n_elig);
  const arma::uword k = regularized ? ncp + 1 : ncp;
  const arma::uword n_gram = tall ? n_elig : nrX;

  // Spectra needs 1 <= nev <= n_gram - 1.
  if (k >= n_gram)
  {
    Rcpp::stop("ncp too large for the Spectra solver: need ncp + regularized < min(n_rows, n_eligible_cols)");
  }

  const double nrX_d = static_cast<double>(nrX);
  const double n_elig_d = static_cast<double>(n_elig);
  const double ncp_d = static_cast<double>(ncp);

  const double denom_sigma = (n_elig_d - ncp_d) * (nrX_d - 1.0 - ncp_d);
  if (regularized && !(denom_sigma > 0.0))
  {
    Rcpp::stop("regularized=true requires n_elig > ncp and nrX - 1 > ncp");
  }
  const double min_dim = std::min(n_elig_d, nrX_d - 1.0);
  const double scale_factor = (nrX_d * n_elig_d) / min_dim;

  // -------------------------------------------------------------------------
  // pass 1: per-column missing counts and flat missing-row index.
  // -------------------------------------------------------------------------
  LOC_TIC(pca_imp_spectra, "pass1_count_missing");
  arma::uvec col_nmiss(n_elig, arma::fill::zeros);
  for (arma::uword j = 0; j < n_elig; ++j)
  {
    const double *p = obj.colptr(eligible_idx[j]);
    arma::uword c = 0;
    for (arma::uword i = 0; i < nrX; ++i)
    {
      c += static_cast<arma::uword>(std::isnan(p[i]));
    }
    col_nmiss[j] = c;
  }

  arma::uvec miss_offsets(n_elig + 1, arma::fill::zeros);
  for (arma::uword j = 0; j < n_elig; ++j)
  {
    miss_offsets[j + 1] = miss_offsets[j] + col_nmiss[j];
  }
  const arma::uword n_missing = miss_offsets[n_elig];
  arma::uvec miss_rows(n_missing);
  LOC_TOC(pca_imp_spectra, "pass1_count_missing");

  // -------------------------------------------------------------------------
  // pass 2: copy eligible columns, record missing rows, weighted moments
  // over observed entries, initial standardization with missing -> 0
  // (mean-impute start, matching init = 0).
  // -------------------------------------------------------------------------
  LOC_TIC(pca_imp_spectra, "pass2_scan");
  arma::mat Xhat(nrX, n_elig);

  for (arma::uword j = 0; j < n_elig; ++j)
  {
    const double *src = obj.colptr(eligible_idx[j]);
    double *dst = Xhat.colptr(j);
    arma::uword p = miss_offsets[j];

    double d_acc = 0.0, ws_acc = 0.0;
    for (arma::uword i = 0; i < nrX; ++i)
    {
      const double v = src[i];
      if (std::isnan(v))
      {
        miss_rows[p++] = i;
      }
      else
      {
        const double w = wptr[i];
        d_acc += w;
        ws_acc += w * v;
      }
    }
    if (!(d_acc > 0.0))
    {
      Rcpp::stop("eligible column has zero observed weight (column %u)",
                 static_cast<unsigned>(eligible_idx[j] + 1));
    }
    const double mu = ws_acc / d_acc;

    double s_eff = 1.0;
    if (scale)
    {
      double ssq = 0.0;
      for (arma::uword i = 0; i < nrX; ++i)
      {
        const double v = src[i];
        if (!std::isnan(v))
        {
          const double vc = v - mu;
          ssq += wptr[i] * vc * vc;
        }
      }
      const double sd = std::sqrt(std::max(0.0, ssq / d_acc));
      s_eff = (sd < PCA_TOL) ? PCA_TOL : sd;
    }

    const double inv_s = 1.0 / s_eff;
    for (arma::uword i = 0; i < nrX; ++i)
    {
      const double v = src[i];
      dst[i] = std::isnan(v) ? 0.0 : (v - mu) * inv_s;
    }
  }
  LOC_TOC(pca_imp_spectra, "pass2_scan");

  // -------------------------------------------------------------------------
  // workspaces.
  // -------------------------------------------------------------------------
  const arma::colvec sqrt_row_w = arma::sqrt(row_w.t());
  arma::colvec inv_sqrt_row_w(nrX);
  for (arma::uword i = 0; i < nrX; ++i)
  {
    const double s = sqrt_row_w[i];
    inv_sqrt_row_w[i] = (s > PCA_TOL) ? (1.0 / s) : (1.0 / PCA_TOL);
  }
  const double *sptr = sqrt_row_w.memptr();
  const double *isptr = inv_sqrt_row_w.memptr();

  GramWorkspace gram_ws;
  gram_ws.init(nrX, n_elig, tall);
  GramCache gram_cache; // inactive: full syrk each iteration

  arma::mat AA_NxN(n_gram, n_gram, arma::fill::zeros);
  arma::mat X_work;
  if (tall)
  {
    X_work.set_size(nrX, n_elig);
  }
  arma::mat fittedX(nrX, n_elig, arma::fill::zeros);
  arma::mat U, V, Vn, eigvecs;
  arma::vec eigvals, vs_top, d_inv;
  arma::vec lambda_shrinked(ncp);

  SpectraOptions opt;
  opt.tol = spectra_tol;
  opt.maxiter = static_cast<int>(spectra_maxiter);
  opt.ncv = static_cast<int>(spectra_ncv);
  SpectraState state;

  std::vector<int> eig_niter, eig_nops;
  std::vector<double> eig_rel_res;
  eig_niter.reserve(64);
  eig_nops.reserve(64);
  eig_rel_res.reserve(64);

  double old = arma::datum::inf;
  double objective = 0.0;
  arma::uword n_iter_final = 0;
  double criterion_final = arma::datum::nan;
  bool converged = false;

  // -------------------------------------------------------------------------
  // EM loop.
  // -------------------------------------------------------------------------
  for (arma::uword nb_iter = 1;; ++nb_iter)
  {
    if (nb_iter % 5 == 0)
    {
      Rcpp::checkUserInterrupt();
    }

    // impute missing entries from the previous fit (zeros on the first
    // iteration), then re-standardize each column with the imputed values
    // included — the plain equivalent of impute_restandardize(): all moments
    // are computed in the current standardized frame, so the update composes
    // with the frame chain without tracking absolute means/sds.
    LOC_TIC(pca_imp_spectra, "restandardize");
    for (arma::uword j = 0; j < n_elig; ++j)
    {
      double *col = Xhat.colptr(j);
      const double *fit = fittedX.colptr(j);

      const arma::uword beg = miss_offsets[j];
      const arma::uword end = miss_offsets[j + 1];
      for (arma::uword p = beg; p < end; ++p)
      {
        const arma::uword i = miss_rows[p];
        col[i] = fit[i];
      }

      double m1 = 0.0, m2 = 0.0;
      for (arma::uword i = 0; i < nrX; ++i)
      {
        const double v = col[i];
        const double w = wptr[i];
        m1 += w * v;
        m2 += w * v * v;
      }

      if (scale)
      {
        const double var = std::max(0.0, m2 - m1 * m1);
        const double sd = std::sqrt(var);
        const double s_eff = (sd < PCA_TOL) ? PCA_TOL : sd;
        const double inv_s = 1.0 / s_eff;
        for (arma::uword i = 0; i < nrX; ++i)
        {
          col[i] = (col[i] - m1) * inv_s;
        }
      }
      else
      {
        for (arma::uword i = 0; i < nrX; ++i)
        {
          col[i] -= m1;
        }
      }
    }
    LOC_TOC(pca_imp_spectra, "restandardize");

    LOC_TIC(pca_imp_spectra, "form_gram");
    if (tall)
    {
      pca_detail::copy_scale_rows(X_work.memptr(), Xhat.memptr(), sptr,
                                  nrX, n_elig);
    }
    form_weighted_gram(Xhat, sptr, tall, X_work, AA_NxN, gram_ws, gram_cache);
    LOC_TOC(pca_imp_spectra, "form_gram");

    const double trace_val = arma::trace(AA_NxN);

    LOC_TIC(pca_imp_spectra, "eig");
    SpectraResult sres;
    const bool ok = spectra_topk_eig(eigvals, eigvecs, AA_NxN, k,
                                     state, opt, sres);
    LOC_TOC(pca_imp_spectra, "eig");

    // residual verification outside the timed block: Spectra can report
    // convergence with spurious Ritz values in degenerate warm starts.
    double rel_res = std::numeric_limits<double>::quiet_NaN();
    if (ok)
    {
      rel_res = spectra_max_rel_res(AA_NxN, eigvals, eigvecs);
    }
    eig_niter.push_back(sres.niter);
    eig_nops.push_back(sres.nops);
    eig_rel_res.push_back(rel_res);

    if (!ok || !(rel_res < 1e-6))
    {
      Rcpp::stop("Spectra eigensolver failed at EM iteration %u (nconv = %d, max_rel_res = %g)",
                 static_cast<unsigned>(nb_iter), sres.nconv, rel_res);
    }

    // ---- post_svd: same shrinkage and loss pieces as armaSVD.cpp ----------
    LOC_TIC(pca_imp_spectra, "post_svd");
    const arma::uword ke = eigvals.n_elem;
    vs_top.set_size(ke);
    d_inv.set_size(ncp);
    for (arma::uword i = 0; i < ke; ++i)
    {
      const double e = eigvals[i];
      vs_top[i] = (e > 0.0) ? std::sqrt(e) : 0.0;
    }
    for (arma::uword i = 0; i < ncp; ++i)
    {
      d_inv[i] = (vs_top[i] > PCA_TOL) ? (1.0 / vs_top[i]) : 0.0;
    }

    const double sum_top_sq = arma::dot(vs_top.head(ncp), vs_top.head(ncp));
    const double tail = std::max(0.0, trace_val - sum_top_sq);
    double sigma2 = 0.0;
    if (regularized)
    {
      sigma2 = scale_factor * tail / denom_sigma;
      sigma2 = std::min(sigma2 * coeff_ridge, vs_top(ncp) * vs_top(ncp));
    }

    double top_corr = 0.0;
    for (arma::uword i = 0; i < ncp; ++i)
    {
      const double v = vs_top[i];
      const double lam = (v > PCA_TOL) ? (v - sigma2 / v) : 0.0;
      lambda_shrinked[i] = lam;
      const double d = v - lam;
      top_corr += d * d;
    }
    const double objective_full = tail + top_corr;
    LOC_TOC(pca_imp_spectra, "post_svd");

    // ---- recover factors in weighted space (as SVD_triplet) ---------------
    LOC_TIC(pca_imp_spectra, "recover_factors");
    if (tall)
    {
      split_eigvec_block(eigvecs, V, Vn, d_inv, nullptr, ncp);
      pca_detail::gemm_nn(X_work, Vn, U);
    }
    else
    {
      split_eigvec_block(eigvecs, U, Vn, d_inv, sptr, ncp);
      pca_detail::gemm_tn(Xhat, Vn, V);
    }
    LOC_TOC(pca_imp_spectra, "recover_factors");

    // ---- reconstruct full fitted matrix -----------------------------------
    LOC_TIC(pca_imp_spectra, "reconstruct");
    pca_detail::scale_rows_cols_inplace(U.memptr(), isptr,
                                        lambda_shrinked.memptr(),
                                        U.n_rows, ncp);
    gemm_nt_full(U, V, fittedX);
    LOC_TOC(pca_imp_spectra, "reconstruct");

    // ---- objective: same formula as armaSVD.cpp ---------------------------
    LOC_TIC(pca_imp_spectra, "objective");
    double miss_contrib = 0.0;
    for (arma::uword j = 0; j < n_elig; ++j)
    {
      const double *xh = Xhat.colptr(j);
      const double *fx = fittedX.colptr(j);
      const arma::uword beg = miss_offsets[j];
      const arma::uword end = miss_offsets[j + 1];
      double s = 0.0;
      for (arma::uword p = beg; p < end; ++p)
      {
        const arma::uword i = miss_rows[p];
        const double d = xh[i] - fx[i];
        s += wptr[i] * d * d;
      }
      miss_contrib += s;
    }
    objective = std::max(0.0, objective_full - miss_contrib);
    LOC_TOC(pca_imp_spectra, "objective");

    const double criterion = std::abs(1.0 - objective / old);
    old = objective;

    bool stop_now = false;
    if (nb_iter >= miniter)
    {
      if (criterion < threshold || objective < threshold)
      {
        stop_now = true;
        converged = true;
      }
    }
    if (!stop_now && nb_iter >= maxiter)
    {
      stop_now = true;
      converged = false;
      Rcpp::warning("Stopped after " + std::to_string(maxiter) + " iterations");
    }

    if (stop_now)
    {
      n_iter_final = nb_iter;
      criterion_final = criterion;
      break;
    }
  }

  LOC_TOC(pca_imp_spectra, "pca_imp_spectra_cpp_total");

  return Rcpp::List::create(
      Rcpp::Named("n_iter") = static_cast<int>(n_iter_final),
      Rcpp::Named("objective") = objective,
      Rcpp::Named("criterion_final") = criterion_final,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("eigvals") = eigvals,
      Rcpp::Named("n_gram") = n_gram,
      Rcpp::Named("tall") = tall,
      Rcpp::Named("k") = k,
      Rcpp::Named("eig_niter") = eig_niter,
      Rcpp::Named("eig_nops") = eig_nops,
      Rcpp::Named("eig_max_rel_res") = eig_rel_res);
}
