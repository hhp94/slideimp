#ifndef PCA_IMP_DIAGNOSTICS
#define PCA_IMP_DIAGNOSTICS 0
#endif

#include <Rcpp.h>
#include "svd_triplet.h"
#include "matrix_checks.h"
#include "loc_timer.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

// Types, for the whole file. A row index and a column index of obj come from
// R and are int; a count of rows, columns, components or iterations is held
// in std::size_t once it has been checked, so that an offset into
// column-major storage, such as j * nrX, is computed in std::size_t. Matrices
// are std::vector<double>, column-major and unpadded, zero-filled when sized.

template <typename F>
static inline void for_each_missing_in_col(std::size_t cidx,
                                           const std::vector<std::size_t> &offsets,
                                           const std::vector<int> &rows_flat,
                                           F &&f)
{
  const std::size_t beg = offsets[cidx];
  const std::size_t end = offsets[cidx + 1];
  for (std::size_t p = beg; p < end; ++p)
  {
    f(static_cast<std::size_t>(rows_flat[p]));
  }
}

// ---------------------------------------------------------------------------
// impute_restandardize
//
// changing columns occupy positions [n_fixed, n_fixed + n_mc) of the
// permuted Xhat, so column indices are derived rather than passed. et is
// read and written only when Scale.
// ---------------------------------------------------------------------------
template <bool Scale>
void impute_restandardize(double *xptr,
                          const double *fxp,
                          double *mean_p,
                          [[maybe_unused]] double *et,
                          const double *wptr,
                          const double *swptr,
                          double *X_chg_scaled,
                          const std::size_t nrX,
                          const std::size_t n_fixed,
                          const std::size_t n_mc,
                          const std::vector<int> &miss_rows_flat,
                          const std::vector<std::size_t> &miss_rows_offsets,
                          const double *mu0,
                          const double *obs_sum_c,
                          const double *obs_sumsq_c)
{
  for (std::size_t cidx = 0; cidx < n_mc; ++cidx)
  {
    const std::size_t j = n_fixed + cidx;
    const std::size_t off = j * nrX;

    double *col = xptr + off;
    const double *fit_col = fxp + off;

    // Compact changing-column weighted buffer.
    //
    // Tall:
    //   points into column n_fixed of X_work
    //
    // Wide cached:
    //   points into gram_cache.X_chg_scaled
    //
    // Wide non-cached:
    //   nullptr
    double *scaled_col = X_chg_scaled ? (X_chg_scaled + cidx * nrX) : nullptr;

    const double old_mu = mean_p[j];

    // old_mu stays within O(sd) of mu0[j], so this subtraction is exact
    // (Sterbenz) and dmu_old is an O(sd) quantity.
    const double dmu_old = old_mu - mu0[j];

    double old_et;
    if constexpr (Scale)
    {
      old_et = et[j];
    }
    else
    {
      old_et = 1.0;
    }

    double imp_sum_c = 0.0;
    [[maybe_unused]] double imp_sumsq_c = 0.0; // accumulated only when Scale

    for_each_missing_in_col(cidx, miss_rows_offsets, miss_rows_flat,
                            [&](std::size_t i)
                            {
                              const double v_std = fit_col[i];
                              col[i] = v_std;

                              const double w = wptr[i];

                              if constexpr (Scale)
                              {
                                // centered imputed value v - mu0[j], formed
                                // directly at O(sd); never materialize the
                                // O(mean) value v itself.
                                const double v_c = v_std * old_et + dmu_old;
                                imp_sum_c += w * v_c;
                                imp_sumsq_c += w * v_c * v_c;
                              }
                              else
                              {
                                const double v_c = v_std + dmu_old;
                                imp_sum_c += w * v_c;
                              }
                            });

    // total weight over observed + imputed rows is exactly 1 (row_w was
    // normalized), so the centered sums add up to the full mean offset.
    const double dmu_new = obs_sum_c[j] + imp_sum_c;
    const double new_mu = mu0[j] + dmu_new;

    if constexpr (Scale)
    {
      const double new_var = obs_sumsq_c[j] + imp_sumsq_c - dmu_new * dmu_new;
      const double new_et = std::sqrt(std::max(0.0, new_var));
      const double new_et_eff = std::max(new_et, PCA_TOL);

      const double inv_new = 1.0 / new_et_eff;
      const double a = old_et * inv_new;
      // (old_mu - new_mu) * inv_new. The mu0 offset cancels.
      const double b = (dmu_old - dmu_new) * inv_new;

      if (scaled_col)
      {
        for (std::size_t i = 0; i < nrX; ++i)
        {
          const double z = col[i] * a + b;
          col[i] = z;
          scaled_col[i] = z * swptr[i];
        }
      }
      else
      {
        for (std::size_t i = 0; i < nrX; ++i)
        {
          col[i] = col[i] * a + b;
        }
      }

      et[j] = new_et_eff;
    }
    else
    {
      // old_mu - new_mu. The mu0 offset cancels.
      const double b = dmu_old - dmu_new;

      if (scaled_col)
      {
        for (std::size_t i = 0; i < nrX; ++i)
        {
          const double z = col[i] + b;
          col[i] = z;
          scaled_col[i] = z * swptr[i];
        }
      }
      else
      {
        for (std::size_t i = 0; i < nrX; ++i)
        {
          col[i] += b;
        }
      }
    }

    mean_p[j] = new_mu;
  }
}

// Starting values at the missing positions of Xhat and of fittedX, which the
// caller hands over zero-filled: before the first dgemm only those positions
// of fittedX are read.
static inline void initialize_missing_start(double *Xhat,
                                            double *fittedX,
                                            const int init,
                                            const std::size_t nrX,
                                            const std::size_t n_fixed,
                                            const std::size_t n_mc,
                                            const std::vector<int> &miss_rows_flat,
                                            const std::vector<std::size_t> &miss_rows_offsets)
{
  if (n_mc == 0)
  {
    return;
  }

  if (init == 0)
  {
    for (std::size_t ci = 0; ci < n_mc; ++ci)
    {
      const std::size_t j = n_fixed + ci;

      double *col = col_ptr(Xhat, j, nrX);
      double *fit_col = col_ptr(fittedX, j, nrX);

      for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                              [&](std::size_t i)
                              {
                                col[i] = 0.0;
                                fit_col[i] = 0.0;
                              });
    }

    return;
  }

  // the caller has refused more missing cells than an int counts.
  Rcpp::NumericVector rand_vals =
      Rcpp::rnorm(static_cast<int>(miss_rows_flat.size()));

  for (std::size_t ci = 0; ci < n_mc; ++ci)
  {
    const std::size_t j = n_fixed + ci;

    double *col = col_ptr(Xhat, j, nrX);
    double *fit_col = col_ptr(fittedX, j, nrX);

    std::size_t p = miss_rows_offsets[ci];

    for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                            [&](std::size_t i)
                            {
                              const double z = rand_vals[p++];
                              col[i] = z;
                              fit_col[i] = z;
                            });
  }
}

static inline double missing_weighted_residual_ss(
    const double *xhp,
    const double *fxp,
    const double *wptr,
    const std::size_t nrX,
    const std::size_t n_fixed,
    const std::size_t n_mc,
    const std::vector<int> &miss_rows_flat,
    const std::vector<std::size_t> &miss_rows_offsets)
{
  double out = 0.0;

  for (std::size_t ci = 0; ci < n_mc; ++ci)
  {
    const std::size_t j = n_fixed + ci;
    const std::size_t off = j * nrX;

    const double *xh = xhp + off;
    const double *fx = fxp + off;

    double s = 0.0;

    for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                            [&](std::size_t i)
                            {
                              const double d = xh[i] - fx[i];
                              s += wptr[i] * d * d;
                            });

    out += s;
  }

  return out;
}

// dst[:, j0 .. j0+nj) = U * V[j0 .. j0+nj, :]', by dgemm on the strided row
// block of V. U is u_rows x u_cols, V has v_rows rows and u_cols columns,
// and dst has u_rows rows.
static inline void reconstruct_cols(const double *U,
                                    const std::size_t u_rows,
                                    const std::size_t u_cols,
                                    const double *V,
                                    const std::size_t v_rows,
                                    double *dst,
                                    const std::size_t j0,
                                    const std::size_t nj)
{
  if (nj == 0)
  {
    return;
  }
  const int st = si_la::gemm('N', 'T', u_rows, nj, u_cols,
                             1.0,
                             U, u_rows,
                             V + j0, v_rows,
                             0.0,
                             col_ptr(dst, j0, u_rows), u_rows);
  if (st != si_la::ok)
  {
    Rcpp::stop("internal error: dgemm refused its arguments (status %d)", st);
  }
}

// Internal. `eligible_idx` holds 0-based columns of `obj`. Returns the
// imputed cells as `loc`, an integer matrix of each one's 1-based row and
// column in `obj`, and `value`, in the column-major order of the eligible
// columns with the complete ones dropped.
// [[Rcpp::export]]
Rcpp::List pca_imp_internal_cpp(
    const Rcpp::NumericMatrix &obj,
    const Rcpp::IntegerVector &eligible_idx,
    const int ncp,
    const bool scale,
    const bool regularized,
    const double threshold,
    const int init,
    const int maxiter,
    const int miniter,
    const Rcpp::NumericVector &row_w,
    const double coeff_ridge,
    const int solver,
    const int warmup_iters,
    const double lobpcg_tol,
    const int lobpcg_maxiter,
    const int trace_iter = 0)
{
  LOC_TIMER_OBJ(pca_imp_gram);
  LOC_TIC(pca_imp_gram, "pca_imp_internal_cpp_total");
  // No Inf scan here. pca_imp() calls this once per `nb.init` on an unchanged
  // `obj`, so a guard at this level would re-scan the whole matrix on every
  // initialization. check_finite() in pca_imp() is the gate, and it runs once.

  // -------------------------------------------------------------------------
  // basic argument validation before we set up workspaces. NA_integer_ is
  // INT_MIN, so every lower bound below refuses it too.
  // -------------------------------------------------------------------------
  const std::size_t nrX = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_elig = static_cast<std::size_t>(eligible_idx.size());

  if (maxiter < 1)
  {
    Rcpp::stop("maxiter must be >= 1");
  }
  if (miniter < 0)
  {
    Rcpp::stop("miniter must be >= 0");
  }
  if (miniter > maxiter)
  {
    Rcpp::stop("miniter must be <= maxiter");
  }
  if (ncp < 1)
  {
    Rcpp::stop("ncp must be >= 1");
  }
  if (init < 0)
  {
    Rcpp::stop("init must be >= 0");
  }
  if (trace_iter < 0)
  {
    Rcpp::stop("trace_iter must be >= 0");
  }
  if (solver < 0 || solver > 2)
  {
    Rcpp::stop("solver must be 0 (exact), 1 (lobpcg), or 2 (auto)");
  }

  const bool lobpcg_requested = (solver != 0);

  if (lobpcg_requested && lobpcg_maxiter < 1)
  {
    Rcpp::stop("lobpcg_maxiter must be >= 1 when solver is 1 (lobpcg) or 2 (auto)");
  }
  if (lobpcg_requested && warmup_iters < 0)
  {
    Rcpp::stop("warmup_iters must be >= 0 when solver is 1 (lobpcg) or 2 (auto)");
  }

  if (n_elig == 0)
  {
    Rcpp::stop("eligible_idx is empty");
  }
  if (static_cast<std::size_t>(row_w.size()) != nrX)
  {
    Rcpp::stop("row_w length must equal n_rows(obj)");
  }
  stop_on_bad_col_index(eligible_idx, obj.ncol(), "eligible_idx");

  const std::size_t n_cp = static_cast<std::size_t>(ncp);
  const std::size_t max_iter = static_cast<std::size_t>(maxiter);
  const std::size_t min_iter = static_cast<std::size_t>(miniter);
  const std::size_t trace_every = static_cast<std::size_t>(trace_iter);

  const std::size_t min_problem_dim = std::min(nrX, n_elig);
  if (n_cp >= min_problem_dim)
  {
    Rcpp::stop("ncp + 1 must be <= min(n_rows, n_eligible_cols)");
  }

  // prevent nrX * n_elig overflow later.
  if (nrX > (std::numeric_limits<std::size_t>::max)() / n_elig)
  {
    Rcpp::stop("problem dimensions are too large");
  }

  const bool tall = (nrX >= n_elig);
  const std::size_t k = regularized ? n_cp + 1 : n_cp;
  const std::size_t n_gram = tall ? n_elig : nrX;

  // LAPACK workspace guard. The dimensions themselves fit an int: nrX is an
  // R dimension, and n_cp and n_gram are below one.
  {
    const std::size_t int_max =
        static_cast<std::size_t>((std::numeric_limits<int>::max)());

    // si_la::syevr_range gives dsyevr at least 26*n and 10*n workspace.
    if (n_gram > int_max / 26)
    {
      Rcpp::stop("problem dimension too large for LAPACK dsyevr workspace");
    }

    // LOBPCG internally solves small Rayleigh-Ritz problems up to ~3*k.
    // only apply this guard when the requested solver can actually use
    // LOBPCG (lobpcg_maxiter >= 1 was already enforced above).
    if (lobpcg_requested)
    {
      const long double rr_dim = 3.0L * static_cast<long double>(k);
      const long double dsyevd_lwork = 1.0L + 6.0L * rr_dim + 2.0L * rr_dim * rr_dim;
      const long double dsyevd_liwork = 3.0L + 5.0L * rr_dim;

      if (rr_dim > static_cast<long double>(int_max) ||
          dsyevd_lwork > static_cast<long double>(int_max) ||
          dsyevd_liwork > static_cast<long double>(int_max))
      {
        Rcpp::stop("ncp too large for LOBPCG internal LAPACK workspace");
      }
    }
  }

  // the kernel's own copy of the weights, normalized below.
  std::vector<double> w(row_w.begin(), row_w.end());

  for (std::size_t i = 0; i < nrX; ++i)
  {
    if (!std::isfinite(w[i]))
    {
      Rcpp::stop("row_w must be finite");
    }
  }
  double row_w_sum = 0.0;
  for (std::size_t i = 0; i < nrX; ++i)
  {
    if (w[i] < 0.0)
    {
      Rcpp::stop("row_w must be non-negative");
    }
    row_w_sum += w[i];
  }
  if (!(row_w_sum > 0.0))
  {
    Rcpp::stop("row_w must have positive sum");
  }

  const double nrX_d = static_cast<double>(nrX);
  const double n_elig_d = static_cast<double>(n_elig);
  const double ncp_d = static_cast<double>(n_cp);

  const double denom_sigma = (n_elig_d - ncp_d) * (nrX_d - 1.0 - ncp_d);
  if (regularized && !(denom_sigma > 0.0))
  {
    Rcpp::stop("regularized=true requires n_elig > ncp and nrX - 1 > ncp");
  }

  double old = std::numeric_limits<double>::infinity();

  for (std::size_t i = 0; i < nrX; ++i)
  {
    w[i] /= row_w_sum;
  }
  const double *wptr = w.data();

  const double *objp = obj.begin();
  const int *elig = eligible_idx.begin();

  LOC_TIC(pca_imp_gram, "pass1_count_missing");
  std::vector<std::size_t> col_nmiss(n_elig, 0);
  for (std::size_t j_local = 0; j_local < n_elig; ++j_local)
  {
    const double *p =
        col_ptr(objp, static_cast<std::size_t>(elig[j_local]), nrX);
    std::size_t c = 0;
    for (std::size_t i = 0; i < nrX; ++i)
    {
      c += static_cast<std::size_t>(std::isnan(p[i]));
    }
    col_nmiss[j_local] = c;
  }
  LOC_TOC(pca_imp_gram, "pass1_count_missing");

  // perm lists the complete columns first, then the ones with a hole, each
  // group in the order eligible_idx gives.
  std::vector<std::size_t> perm;
  perm.reserve(n_elig);
  for (std::size_t j = 0; j < n_elig; ++j)
  {
    if (col_nmiss[j] == 0)
    {
      perm.push_back(j);
    }
  }
  const std::size_t n_fixed = perm.size();
  for (std::size_t j = 0; j < n_elig; ++j)
  {
    if (col_nmiss[j] > 0)
    {
      perm.push_back(j);
    }
  }
  const std::size_t n_mc = n_elig - n_fixed;

  std::vector<int> eligible_idx_perm(n_elig);
  for (std::size_t jp = 0; jp < n_elig; ++jp)
  {
    eligible_idx_perm[jp] = elig[perm[jp]];
  }

  // after the permutation, changing columns are exactly positions
  // [n_fixed, n_elig); loops derive local column indices as n_fixed + ci.
  std::vector<std::size_t> miss_rows_offsets(n_mc + 1, 0);
  for (std::size_t ci = 0; ci < n_mc; ++ci)
  {
    miss_rows_offsets[ci + 1] =
        miss_rows_offsets[ci] + col_nmiss[perm[n_fixed + ci]];
  }
  const std::size_t n_missing = miss_rows_offsets[n_mc];

  // `loc` is an R matrix, whose dims are int.
  if (n_missing > static_cast<std::size_t>((std::numeric_limits<int>::max)()))
  {
    Rcpp::stop(
        std::string("The eligible columns hold ") +
        std::to_string(n_missing) +
        " missing values; the result has a row for each, and an R "
        "matrix has at most " +
        std::to_string((std::numeric_limits<int>::max)()) +
        " rows.");
  }
  std::vector<int> miss_rows_flat(n_missing);

  std::vector<double> Xhat(nrX * n_elig);
  std::vector<double> fittedX(nrX * n_elig);

  std::vector<double> denom(n_elig);
  std::vector<double> mu0(n_elig);
  std::vector<double> obs_sum_c(n_elig);
  std::vector<double> obs_sumsq_c(n_elig);

  LOC_TIC(pca_imp_gram, "pass2_scan");
  {
    for (std::size_t jp = 0; jp < n_elig; ++jp)
    {
      const double *src =
          col_ptr(objp, static_cast<std::size_t>(eligible_idx_perm[jp]), nrX);
      double *dst = col_ptr(Xhat.data(), jp, nrX);

      const bool has_miss = (jp >= n_fixed);
      const std::size_t cidx = has_miss ? (jp - n_fixed) : 0;
      std::size_t p = has_miss ? miss_rows_offsets[cidx] : 0;

      // pass A: observed weight, raw weighted sum, missing-row indices.
      double d_acc = 0.0, ws_acc = 0.0;
      for (std::size_t i = 0; i < nrX; ++i)
      {
        const double v = src[i];
        if (std::isnan(v))
        {
          if (has_miss)
          {
            miss_rows_flat[p++] = static_cast<int>(i);
          }
        }
        else
        {
          const double wi = wptr[i];
          d_acc += wi;
          ws_acc += wi * v;
        }
      }
      denom[jp] = d_acc;

      // per-column reference mean for centered moment tracking. The
      // d_acc == 0 case is rejected right after this loop; 0.0 just keeps
      // the pass-B arithmetic finite until then.
      const double mu0_j = (d_acc > 0.0) ? (ws_acc / d_acc) : 0.0;
      mu0[jp] = mu0_j;

      // pass B: centered moments at O(sd) magnitude (avoids the
      // E[x^2] - mu^2 catastrophic cancellation), plus the Xhat write.
      double osc_acc = 0.0, ossq_acc = 0.0;
      for (std::size_t i = 0; i < nrX; ++i)
      {
        const double v = src[i];
        if (std::isnan(v))
        {
          dst[i] = 0.0;
        }
        else
        {
          dst[i] = v;
          const double vc = v - mu0_j;
          const double wi = wptr[i];
          osc_acc += wi * vc;
          ossq_acc += wi * vc * vc;
        }
      }
      obs_sum_c[jp] = osc_acc; // ~0 by construction; kept for exactness
      obs_sumsq_c[jp] = ossq_acc;
    }
  }
  LOC_TOC(pca_imp_gram, "pass2_scan");

  // every eligible column must have at least one observed row with positive
  // weight, else mean_p (and the downstream Gram) is NaN.
  for (std::size_t j = 0; j < n_elig; ++j)
  {
    if (!(denom[j] > 0.0))
    {
      Rcpp::stop("eligible column has zero observed weight (column %d)",
                 eligible_idx_perm[j] + 1);
    }
  }

  const double *mu0ptr = mu0.data();
  const double *oscptr = obs_sum_c.data();
  const double *ossqptr = obs_sumsq_c.data();

  // the column means, and when scaling the column scales floored at PCA_TOL.
  // et stays empty otherwise, and nothing reads it.
  std::vector<double> mean_p(n_elig);
  std::vector<double> et;
  if (scale)
  {
    et.resize(n_elig);
  }
  for (std::size_t j = 0; j < n_elig; ++j)
  {
    const double mean_resid = obs_sum_c[j] / denom[j];
    mean_p[j] = mu0[j] + mean_resid;
    if (scale)
    {
      const double var = obs_sumsq_c[j] / denom[j] - mean_resid * mean_resid;
      const double e = std::sqrt(var < 0.0 ? 0.0 : var);
      et[j] = e < PCA_TOL ? PCA_TOL : e;
    }
  }

  if (scale)
  {
    for (std::size_t j = 0; j < n_elig; ++j)
    {
      double *col = col_ptr(Xhat.data(), j, nrX);
      const double mu = mean_p[j];
      const double inv_et = 1.0 / et[j];

      for (std::size_t i = 0; i < nrX; ++i)
      {
        col[i] = (col[i] - mu) * inv_et;
      }
    }
  }
  else
  {
    for (std::size_t j = 0; j < n_elig; ++j)
    {
      double *col = col_ptr(Xhat.data(), j, nrX);
      const double mu = mean_p[j];

      for (std::size_t i = 0; i < nrX; ++i)
      {
        col[i] -= mu;
      }
    }
  }

  initialize_missing_start(Xhat.data(), fittedX.data(), init, nrX,
                           n_fixed, n_mc, miss_rows_flat, miss_rows_offsets);

  // ---------------------------------------------------------------------------
  // SVD / solve buffers.
  // ---------------------------------------------------------------------------
  double trace_val;
  double sigma2 = 0.0;
  std::vector<double> sqrt_row_w(nrX);
  std::vector<double> inv_sqrt_row_w(nrX);
  for (std::size_t i = 0; i < nrX; ++i)
  {
    const double s = std::sqrt(w[i]);
    sqrt_row_w[i] = s;
    inv_sqrt_row_w[i] = (s > PCA_TOL) ? (1.0 / s) : (1.0 / PCA_TOL);
  }
  const double *sptr = sqrt_row_w.data();
  const double *isptr = inv_sqrt_row_w.data();

  GramWorkspace gram_ws;
  gram_ws.init(nrX, n_elig, tall);
  EigSymRange eig_range;
  if (!eig_range.init(n_gram, k))
  {
    Rcpp::stop("internal error: eigenpair count out of range");
  }

  // every buffer SVD_triplet() reads or writes is sized here, once, and
  // written through raw pointers; the sizes are the ones its header lists.
  // In the wide case V's fixed rows stay unwritten until postprocessing, and
  // X_work is empty.
  std::vector<double> AA_NxN(n_gram * n_gram);
  std::vector<double> X_work(tall ? nrX * n_elig : 0);
  std::vector<double> U(nrX * n_cp);
  std::vector<double> V(n_elig * n_cp);
  std::vector<double> Vn(n_gram * n_cp);
  std::vector<double> eigvecs(n_gram * k);
  std::vector<double> eigvals(k);
  std::vector<double> vs_top(k);
  std::vector<double> d_inv(n_cp);
  std::vector<double> lambda_shrinked(n_cp);

  LOC_TIC(pca_imp_gram, "gram_cache_init");
  GramCache gram_cache;
  const int gram_init_st =
      gram_cache_init(gram_cache, Xhat.data(), sptr, nrX, n_fixed, n_mc,
                      tall, X_work.data());
  LOC_TOC(pca_imp_gram, "gram_cache_init");
  if (gram_init_st != si_la::ok)
  {
    Rcpp::stop("internal error: dsyrk refused its arguments (status %d)",
               gram_init_st);
  }
  double *weighted_chg_ptr = gram_changing_block(
      gram_cache, X_work.data(), nrX, n_fixed, n_mc, tall);

  HybridEigContext hyb_ctx;
  hyb_ctx.eig_range = &eig_range;
  hyb_ctx.solver_mode = static_cast<HybridSolverMode>(solver);
  hyb_ctx.warmup_iters = lobpcg_requested ? warmup_iters : 0;
  hyb_ctx.lobpcg_opt.tol = lobpcg_tol;
  hyb_ctx.lobpcg_opt.maxiter = lobpcg_requested ? lobpcg_maxiter : 0;
#if PCA_IMP_DIAGNOSTICS
  hyb_ctx.init_logs(max_iter);
#endif

  const double min_dim = std::min(n_elig_d, nrX_d - 1.0);
  const double scale_factor = (nrX_d * n_elig_d) / min_dim;

  // ---------------------------------------------------------------------------
  // hot loop prep
  // ---------------------------------------------------------------------------
  double objective = 0.0;
  const double *xhp = Xhat.data();
  double *etptr = scale ? et.data() : nullptr;

  std::size_t n_iter_final = 0;
  double criterion_final = std::numeric_limits<double>::quiet_NaN();
  bool converged = false;

  struct PrecGuard
  {
    std::streamsize old;

    PrecGuard() : old(Rcpp::Rcout.precision()) {}

    ~PrecGuard()
    {
      Rcpp::Rcout.precision(old);
    }
  };

  PrecGuard prec_guard;

  const int iter_w = 4;
  const int obj_w = 8;
  const int ratio_w = 12;

  if (trace_every > 0)
  {
    Rcpp::Rcout.precision(6);
    Rcpp::Rcout
        << "iter : "
        << std::right << std::setw(obj_w) << "obj"
        << " "
        << std::right << std::setw(ratio_w) << "ratio"
        << std::endl;
  }

  auto print_iter = [&](std::size_t it, double obj_val, double crit)
  {
    Rcpp::Rcout
        << std::setw(iter_w) << it << " : "
        << std::setw(obj_w) << obj_val
        << " "
        << std::setw(ratio_w) << crit
        << std::endl;
  };

#if PCA_IMP_DIAGNOSTICS
  // eigval_hist is n_cp x max_iter. V_prev holds the leading n_cp
  // eigenvectors of the iteration before, and is empty during the first.
  std::vector<double> eigval_hist(n_cp * max_iter);
  std::vector<double> obj_hist(max_iter);
  std::vector<double> subspace_cos_hist(max_iter);
  std::vector<double> V_prev;
  std::vector<double> V_cross(n_cp * n_cp);
  std::vector<double> V_cross_sv(n_cp);
#endif

  for (std::size_t nb_iter = 1;; ++nb_iter)
  {
    if (nb_iter % 5 == 0)
    {
      Rcpp::checkUserInterrupt();
    }
    LOC_TIC(pca_imp_gram, "restandardize");
    if (scale)
    {
      impute_restandardize<true>(Xhat.data(), fittedX.data(),
                                 mean_p.data(), etptr,
                                 wptr, sptr, weighted_chg_ptr, nrX,
                                 n_fixed, n_mc,
                                 miss_rows_flat, miss_rows_offsets,
                                 mu0ptr, oscptr, ossqptr);
    }
    else
    {
      impute_restandardize<false>(Xhat.data(), fittedX.data(),
                                  mean_p.data(), etptr,
                                  wptr, sptr, weighted_chg_ptr, nrX,
                                  n_fixed, n_mc,
                                  miss_rows_flat, miss_rows_offsets,
                                  mu0ptr, oscptr, ossqptr);
    }
    LOC_TOC(pca_imp_gram, "restandardize");

    LOC_TIC(pca_imp_gram, "svd");
    SVD_triplet(Xhat.data(), sptr, n_cp, tall,
                vs_top.data(), trace_val, U.data(), V.data(),
                X_work.data(), AA_NxN.data(),
                d_inv.data(), Vn.data(), eigvals.data(), eigvecs.data(),
                gram_ws, gram_cache,
                hyb_ctx, static_cast<int>(nb_iter - 1) LOC_TIMER_ARG(pca_imp_gram));
    LOC_TOC(pca_imp_gram, "svd");

    LOC_TIC(pca_imp_gram, "post_svd");
    sigma2 = 0.0;
    double sum_top_sq = 0.0;
    for (std::size_t cp = 0; cp < n_cp; ++cp)
    {
      sum_top_sq += vs_top[cp] * vs_top[cp];
    }
    // clamp tail energy once and use consistently in both sigma2 (was unclamped)
    // and the objective.
    const double tail = std::max(0.0, trace_val - sum_top_sq);
    if (regularized)
    {
      sigma2 = scale_factor * tail / denom_sigma;
      sigma2 = std::min(sigma2 * coeff_ridge, vs_top[n_cp] * vs_top[n_cp]);
    }

    double top_corr = 0.0;
    {
      const double *vp = vs_top.data();
      double *lp = lambda_shrinked.data();
      for (std::size_t i = 0; i < n_cp; ++i)
      {
        const double v = vp[i];
        const double lam = (v > PCA_TOL) ? (v - sigma2 / v) : 0.0;
        lp[i] = lam;
        const double d = v - lam;
        top_corr += d * d;
      }
    }
    const double objective_full = tail + top_corr;
    LOC_TOC(pca_imp_gram, "post_svd");

    LOC_TIC(pca_imp_gram, "reconstruct");

    // after this, U contains diag(isw) %*% U %*% diag(lambda_shrinked):
    // the row de-weighting of SVD_triplet's weighted-space U fused with the
    // shrunk column scaling.
    pca_detail::scale_rows_cols_inplace(U.data(),
                                        isptr,
                                        lambda_shrinked.data(),
                                        nrX,
                                        n_cp);
    // reconstruct only the changing columns of fittedX; the fixed block is
    // never read inside the loop and is filled once after it exits.
    reconstruct_cols(U.data(), nrX, n_cp, V.data(), n_elig, fittedX.data(),
                     n_fixed, n_mc);

    LOC_TOC(pca_imp_gram, "reconstruct");

    LOC_TIC(pca_imp_gram, "objective");
    {
      const double miss_contrib =
          missing_weighted_residual_ss(Xhat.data(),
                                       fittedX.data(),
                                       wptr,
                                       nrX,
                                       n_fixed,
                                       n_mc,
                                       miss_rows_flat,
                                       miss_rows_offsets);

      objective = std::max(0.0, objective_full - miss_contrib);
    }
    LOC_TOC(pca_imp_gram, "objective");

#if PCA_IMP_DIAGNOSTICS
    std::copy(vs_top.begin(), vs_top.begin() + n_cp,
              eigval_hist.begin() + (nb_iter - 1) * n_cp);
    obj_hist[nb_iter - 1] = objective;
    if (!V_prev.empty())
    {
      // the smallest singular value of V_prev' * V_new: the cosine of the
      // largest angle between the two subspaces.
      int st = si_la::gemm('T', 'N', n_cp, n_cp, n_gram, 1.0,
                           V_prev.data(), n_gram, eigvecs.data(), n_gram,
                           0.0, V_cross.data(), n_cp);
      if (st != si_la::ok)
      {
        Rcpp::stop("internal error: dgemm refused its arguments (status %d)",
                   st);
      }
      st = si_la::svd_values(n_cp, n_cp, V_cross.data(), V_cross_sv.data());
      if (st > 0 || st == si_la::err_nonfinite)
      {
        // dgesdd did not converge, or the product holds NaN or Inf.
        Rcpp::stop("pca_imp diagnostics: the SVD of the subspace product failed");
      }
      if (st == si_la::err_alloc)
      {
        Rcpp::stop("pca_imp diagnostics: could not allocate the workspace of dgesdd");
      }
      if (st != si_la::ok)
      {
        Rcpp::stop("internal error: dgesdd refused its arguments (status %d)", st);
      }
      subspace_cos_hist[nb_iter - 1] =
          *std::min_element(V_cross_sv.begin(), V_cross_sv.end());
    }
    else
    {
      subspace_cos_hist[nb_iter - 1] = std::numeric_limits<double>::quiet_NaN();
    }
    V_prev.assign(eigvecs.begin(), eigvecs.begin() + n_gram * n_cp);
#endif

    // missMDA convergence test:
    //
    // criterion <- abs(1 - objective / old)
    // old <- objective
    //
    // missMDA increments nb.iter before testing `nb.iter > 5`.
    // Since `nb_iter` here is the number of completed iterations,
    // `nb_iter >= miniter` with default miniter = 5 reproduces that behavior.
    const double criterion = std::abs(1.0 - objective / old);
    old = objective;

    bool stop_now = false;
    if (nb_iter >= min_iter)
    {
      if (criterion < threshold || objective < threshold)
      {
        stop_now = true;
        converged = true;
      }
    }
    if (!stop_now && nb_iter >= max_iter)
    {
      stop_now = true;
      converged = false;
      // deliberately silent. The std::string overload of Rcpp::warning is a
      // bare Rf_warning, so under options(warn = 2) - or any calling handler
      // that invokes a restart - R longjmps straight out of this frame.
      // END_RCPP catches C++ exceptions, not longjmps, so Xhat, fittedX,
      // X_work, AA_NxN, the eigen workspaces, PrecGuard and the wrapper's
      // RNGScope would all be skipped. `converged` is already in the returned
      // list and pca_imp() warns on it once, for the winning restart only.
    }

    if (trace_every > 0)
    {
      const bool periodic = (nb_iter == 1) || (nb_iter % trace_every == 0);
      if (periodic || stop_now)
        print_iter(nb_iter, objective, criterion);
    }

    if (stop_now)
    {
      n_iter_final = nb_iter;
      // A fit stopped by its objective alone has no criterion to report:
      // there the ratio is of two values near zero, NaN, Inf or anything
      // between by rounding.
      criterion_final =
          (converged && !(criterion < threshold)) ? NA_REAL : criterion;
      break;
    }
  }

  // ---------------------------------------------------------------------------
  // postprocessing.
  // ---------------------------------------------------------------------------
  LOC_TIC(pca_imp_gram, "postprocessing");
  // the changing columns of fittedX hold the final reconstruction from the
  // last hot-loop DGEMM (U carries the isw/lambda_shrinked scaling); the
  // fixed block was skipped in the loop and is completed here for the SSE.
  // In the wide case SVD_triplet() also left the fixed rows of V unwritten:
  // they come from the last Vn and the fixed columns of Xhat, which nothing
  // has written since pass2_scan. The tall case returns V complete.
  if (!tall)
  {
    const int st = pca_detail::gemm_tn_block(
        Xhat.data(), nrX, n_elig, 0, n_fixed, Vn.data(),
        n_gram, n_cp, V.data());
    if (st != si_la::ok)
    {
      Rcpp::stop("internal error: dgemm refused its arguments (status %d)", st);
    }
  }
  reconstruct_cols(U.data(), nrX, n_cp, V.data(), n_elig, fittedX.data(),
                   0, n_fixed);
  const double *fxp = fittedX.data();

  double sse = 0.0;
  for (std::size_t j = 0; j < n_elig; ++j)
  {
    const double ej = scale ? et[j] : 1.0;
    const double ej2 = ej * ej;
    const std::size_t off = j * nrX;
    const double *xh = xhp + off;
    const double *fx = fxp + off;
    double s = 0.0;
    for (std::size_t i = 0; i < nrX; ++i)
    {
      const double d = xh[i] - fx[i];
      s += d * d;
    }
    sse += s * ej2;
  }
  for (std::size_t ci = 0; ci < n_mc; ++ci)
  {
    const std::size_t j_local = n_fixed + ci;
    const double ej = scale ? et[j_local] : 1.0;
    const double ej2 = ej * ej;
    const std::size_t off = j_local * nrX;
    const double *xh = xhp + off;
    const double *fx = fxp + off;
    double s = 0.0;
    for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                            [&](std::size_t i)
                            {
                              const double d = xh[i] - fx[i];
                              s += d * d;
                            });
    sse -= s * ej2;
  }
  const std::size_t n_obs = nrX * n_elig - n_missing;
  const double mse = (n_obs > 0)
                         ? sse / static_cast<double>(n_obs)
                         : std::numeric_limits<double>::quiet_NaN();

  // n_missing fits an int: refused above otherwise.
  Rcpp::IntegerMatrix loc(static_cast<int>(n_missing), 2);
  Rcpp::colnames(loc) = Rcpp::CharacterVector::create("row", "col");
  Rcpp::NumericVector value(static_cast<R_xlen_t>(n_missing));
  {
    int *loc_row = loc.begin();
    int *loc_col = loc_row + n_missing;
    double *val = value.begin();
    std::size_t iv = 0;
    for (std::size_t ci = 0; ci < n_mc; ++ci)
    {
      const std::size_t j_local = n_fixed + ci;
      const double ej = scale ? et[j_local] : 1.0;
      const double mj = mean_p[j_local];
      const int orig_col_1based = eligible_idx_perm[j_local] + 1;

      // missMDA::imputePCA() returns completeObs using the final Xhat at
      // missing positions, not the final fittedX. fittedX is still used above
      // for the nb.init selection MSE.
      const double *xh = xhp + j_local * nrX;

      for_each_missing_in_col(ci, miss_rows_offsets, miss_rows_flat,
                              [&](std::size_t i)
                              {
                                loc_row[iv] = static_cast<int>(i) + 1;
                                loc_col[iv] = orig_col_1based;
                                val[iv] = xh[i] * ej + mj;
                                ++iv;
                              });
    }
  }

  LOC_TOC(pca_imp_gram, "postprocessing");
  LOC_TOC(pca_imp_gram, "pca_imp_internal_cpp_total");
  const int n_iter_out = static_cast<int>(n_iter_final);
#if PCA_IMP_DIAGNOSTICS
  Rcpp::NumericMatrix eigval_hist_out(ncp, n_iter_out);
  std::copy(eigval_hist.begin(), eigval_hist.begin() + n_cp * n_iter_final,
            eigval_hist_out.begin());
  obj_hist.resize(n_iter_final);
  subspace_cos_hist.resize(n_iter_final);

  hyb_ctx.truncate_logs(n_iter_final);

  return Rcpp::List::create(
      Rcpp::Named("loc") = loc,
      Rcpp::Named("value") = value,
      Rcpp::Named("mse") = mse,
      Rcpp::Named("n_iter") = n_iter_out,
      Rcpp::Named("criterion_final") = criterion_final,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("solver_chosen") = hyb_ctx.solver_chosen_code(),
      Rcpp::Named("eigval_hist") = eigval_hist_out,
      Rcpp::Named("obj_hist") = obj_hist,
      Rcpp::Named("subspace_cos_hist") = subspace_cos_hist,
      Rcpp::Named("n_exact") = hyb_ctx.n_exact,
      Rcpp::Named("n_lobpcg_ok") = hyb_ctx.n_lobpcg_ok,
      Rcpp::Named("n_lobpcg_bad") = hyb_ctx.n_lobpcg_bad,
      Rcpp::Named("auto_n_probe_iter") = hyb_ctx.auto_n_probe_iter,
      Rcpp::Named("auto_min_exact_iter") = hyb_ctx.auto_min_exact_iter,
      Rcpp::Named("auto_exact_n") = hyb_ctx.auto_exact_n,
      Rcpp::Named("auto_lobpcg_n") = hyb_ctx.auto_lobpcg_n,
      Rcpp::Named("auto_exact_mean_sec") = hyb_ctx.auto_exact_mean_sec(),
      Rcpp::Named("auto_lobpcg_mean_sec") = hyb_ctx.auto_lobpcg_mean_sec(),
      Rcpp::Named("eig_path") = hyb_ctx.path_log,
      Rcpp::Named("eig_reason") = hyb_ctx.reason_log,
      Rcpp::Named("eig_lobpcg_iter") = hyb_ctx.lobpcg_iter_log,
      Rcpp::Named("eig_lobpcg_fail_iter") = hyb_ctx.lobpcg_fail_iter_log,
      Rcpp::Named("eig_lobpcg_max_rel_res") = hyb_ctx.lobpcg_max_rel_res_log);
#else
  return Rcpp::List::create(
      Rcpp::Named("loc") = loc,
      Rcpp::Named("value") = value,
      Rcpp::Named("mse") = mse,
      Rcpp::Named("n_iter") = n_iter_out,
      Rcpp::Named("criterion_final") = criterion_final,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("solver_chosen") = hyb_ctx.solver_chosen_code(),
      Rcpp::Named("n_exact") = hyb_ctx.n_exact,
      Rcpp::Named("n_lobpcg_ok") = hyb_ctx.n_lobpcg_ok,
      Rcpp::Named("n_lobpcg_bad") = hyb_ctx.n_lobpcg_bad);
#endif
}
