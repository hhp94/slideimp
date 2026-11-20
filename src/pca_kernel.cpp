#include <Rcpp.h>
#include "svd_triplet.h"
#include "matrix_checks.h"
#include "loc_timer.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <limits>
#include <ostream>
#include <vector>

// Types, for the whole file. A row index and a column index of obj come from
// R and are int; a count of rows, columns, components or iterations is held
// in std::size_t once it has been checked, so that an offset into
// column-major storage, such as j * nrX, is computed in std::size_t. Matrices
// are std::vector<double>, column-major and unpadded, zero-filled when sized.

namespace {

// the EM loop asks for a user interrupt once in this many iterations
constexpr std::size_t INTERRUPT_EVERY = 5;

// The eligible columns as the kernel holds them: permuted so that the
// n_fixed complete ones come first and the n_mc with a hole after them,
// each group in the order eligible_idx gives, every column nrX rows long.
// Column j holds column eligible_idx_perm[j] of obj. Changing column ci is
// column n_fixed + ci, and its holes are at the 0-based rows
// miss_rows_flat[miss_rows_offsets[ci] .. miss_rows_offsets[ci + 1]), in
// ascending order. Built once, by lay_out_columns(), and completed by
// scan_columns(), which fills miss_rows_flat; read-only from there on.
struct ColumnLayout {
  std::size_t nrX = 0;
  std::size_t n_fixed = 0;
  std::size_t n_mc = 0;
  std::vector<int> eligible_idx_perm;         // n_fixed + n_mc
  std::vector<std::size_t> miss_rows_offsets; // n_mc + 1
  std::vector<int> miss_rows_flat;            // miss_rows_offsets[n_mc]
};

// Per column of the layout: the observed rows' weight denom, a reference
// mean mu0, and the weighted sums obs_sum_c and obs_sumsq_c of the observed
// values' deviations from mu0 and of their squares, which pass2_scan takes
// and nothing changes after; then the column mean mean_p and, when scaling,
// the column scale et floored at PCA_TOL, which the standardization sets
// and every EM iteration updates. et stays empty otherwise, and nothing
// reads it. scale_of(j) is the scale of column j: et[j], or 1 when not
// scaling.
struct ColumnMoments {
  std::vector<double> denom;
  std::vector<double> mu0;
  std::vector<double> obs_sum_c;
  std::vector<double> obs_sumsq_c;
  std::vector<double> mean_p;
  std::vector<double> et;

  ColumnMoments(const std::size_t n_elig, const bool scale)
      : denom(n_elig), mu0(n_elig), obs_sum_c(n_elig), obs_sumsq_c(n_elig),
        mean_p(n_elig), et(scale ? n_elig : 0) {}

  double scale_of(const std::size_t j) const {
    return et.empty() ? 1.0 : et[j];
  }
};

template <typename F>
void for_each_missing_in_col(const ColumnLayout &layout, std::size_t cidx,
                             F &&f) {
  const std::size_t beg = layout.miss_rows_offsets[cidx];
  const std::size_t end = layout.miss_rows_offsets[cidx + 1];
  for (std::size_t p = beg; p < end; ++p) {
    f(static_cast<std::size_t>(layout.miss_rows_flat[p]));
  }
}

// ---------------------------------------------------------------------------
// basic argument validation before we set up workspaces. NA_integer_ is
// INT_MIN, so every lower bound below refuses it too. Of row_w only the
// length is checked here; normalized_row_weights() checks its values.
// ---------------------------------------------------------------------------
void check_arguments(const Rcpp::NumericMatrix &obj,
                     const Rcpp::IntegerVector &eligible_idx, const int ncp,
                     const bool regularized, const int init, const int maxiter,
                     const int miniter, const Rcpp::NumericVector &row_w,
                     const int solver, const int warmup_iters,
                     const int lobpcg_maxiter, const int trace_iter) {
  const std::size_t nrX = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_elig = static_cast<std::size_t>(eligible_idx.size());

  if (maxiter < 1) {
    Rcpp::stop("'maxiter' must be >= 1");
  }
  if (miniter < 0) {
    Rcpp::stop("'miniter' must be >= 0");
  }
  if (miniter > maxiter) {
    Rcpp::stop("'miniter' must be <= 'maxiter'");
  }
  if (ncp < 1) {
    Rcpp::stop("'ncp' must be >= 1");
  }
  if (init < 0) {
    Rcpp::stop("'init' must be >= 0");
  }
  if (trace_iter < 0) {
    Rcpp::stop("'trace_iter' must be >= 0");
  }
  if (solver < 0 || solver > 2) {
    Rcpp::stop("'solver' must be 0 (exact), 1 (lobpcg) or 2 (auto)");
  }

  const bool lobpcg_requested = (solver != 0);

  if (lobpcg_requested && lobpcg_maxiter < 1) {
    Rcpp::stop("'lobpcg_maxiter' must be >= 1 when 'solver' is 1 (lobpcg) or 2 "
               "(auto)");
  }
  if (lobpcg_requested && warmup_iters < 0) {
    Rcpp::stop(
        "'warmup_iters' must be >= 0 when 'solver' is 1 (lobpcg) or 2 (auto)");
  }

  if (n_elig == 0) {
    Rcpp::stop("'eligible_idx' is empty");
  }
  if (static_cast<std::size_t>(row_w.size()) != nrX) {
    Rcpp::stop("'row_w' length must equal nrow('obj')");
  }
  stop_on_bad_col_index(eligible_idx, obj.ncol(), "eligible_idx");

  const std::size_t n_cp = static_cast<std::size_t>(ncp);

  const std::size_t min_problem_dim = std::min(nrX, n_elig);
  if (n_cp >= min_problem_dim) {
    Rcpp::stop("'ncp' + 1 must be <= min(nrow('obj'), length('eligible_idx'))");
  }

  // prevent nrX * n_elig overflow later.
  if (nrX > std::numeric_limits<std::size_t>::max() / n_elig) {
    Rcpp::stop("problem dimensions are too large");
  }

  const bool tall = (nrX >= n_elig);
  const std::size_t k = regularized ? n_cp + 1 : n_cp;
  const std::size_t n_gram = tall ? n_elig : nrX;

  // LAPACK workspace guard: the least workspace of each eigensolver the call
  // can reach, as src/si_lapack.h defines it, must fit an int. The
  // dimensions themselves fit one: nrX is an R dimension, and n_cp and
  // n_gram are below one.
  const double int_max = static_cast<double>(std::numeric_limits<int>::max());

  // the exact solver, dsyevr, on the Gram.
  if (si_la::syevr_lwork_min(n_gram) > int_max) {
    Rcpp::stop("problem dimension too large for LAPACK dsyevr workspace");
  }

  // LOBPCG's Rayleigh-Ritz problems, by dsyevd, of order at most 3 k: only
  // when the requested solver can use LOBPCG (lobpcg_maxiter >= 1 was
  // enforced above).
  if (lobpcg_requested) {
    const std::size_t rr_dim = 3 * k;
    if (static_cast<double>(rr_dim) > int_max ||
        si_la::syevd_lwork_min(rr_dim) > int_max ||
        si_la::syevd_liwork_min(rr_dim) > int_max) {
      Rcpp::stop("'ncp' is too large for the LAPACK workspace of LOBPCG");
    }
  }
}

// row_w is checked where it lies, and the kernel normalizes its own copy,
// which this returns. check_arguments() has checked its length.
std::vector<double> normalized_row_weights(const Rcpp::NumericVector &row_w) {
  const std::size_t nrX = static_cast<std::size_t>(row_w.size());
  const double *row_w_ptr = row_w.begin();
  for (std::size_t i = 0; i < nrX; ++i) {
    if (!std::isfinite(row_w_ptr[i])) {
      Rcpp::stop("'row_w' must be finite");
    }
  }
  double row_w_sum = 0.0;
  for (std::size_t i = 0; i < nrX; ++i) {
    if (row_w_ptr[i] < 0.0) {
      Rcpp::stop("'row_w' must be non-negative");
    }
    row_w_sum += row_w_ptr[i];
  }
  if (!(row_w_sum > 0.0)) {
    Rcpp::stop("'row_w' must have a positive sum");
  }
  std::vector<double> w(row_w_ptr, row_w_ptr + nrX);
  for (std::size_t i = 0; i < nrX; ++i) {
    w[i] /= row_w_sum;
  }
  return w;
}

// pass 1: the missing cells of each eligible column of obj, in the order
// eligible_idx gives.
std::vector<std::size_t> count_missing(const double *objp, const int *elig,
                                       const std::size_t nrX,
                                       const std::size_t n_elig) {
  std::vector<std::size_t> col_nmiss(n_elig, 0);
  for (std::size_t j_local = 0; j_local < n_elig; ++j_local) {
    const double *p =
        col_ptr(objp, static_cast<std::size_t>(elig[j_local]), nrX);
    std::size_t c = 0;
    for (std::size_t i = 0; i < nrX; ++i) {
      c += static_cast<std::size_t>(std::isnan(p[i]));
    }
    col_nmiss[j_local] = c;
  }
  return col_nmiss;
}

// The layout of the eligible columns, from their counts of missing cells;
// elig is eligible_idx. miss_rows_flat comes back sized and zero.
ColumnLayout lay_out_columns(const std::vector<std::size_t> &col_nmiss,
                             const int *elig, const std::size_t nrX) {
  const std::size_t n_elig = col_nmiss.size();

  // perm lists the complete columns first, then the ones with a hole, each
  // group in the order eligible_idx gives.
  std::vector<std::size_t> perm;
  perm.reserve(n_elig);
  for (std::size_t j = 0; j < n_elig; ++j) {
    if (col_nmiss[j] == 0) {
      perm.push_back(j);
    }
  }
  const std::size_t n_fixed = perm.size();
  for (std::size_t j = 0; j < n_elig; ++j) {
    if (col_nmiss[j] > 0) {
      perm.push_back(j);
    }
  }
  const std::size_t n_mc = n_elig - n_fixed;

  ColumnLayout layout{nrX,
                      n_fixed,
                      n_mc,
                      std::vector<int>(n_elig),
                      std::vector<std::size_t>(n_mc + 1, 0),
                      {}};
  for (std::size_t jp = 0; jp < n_elig; ++jp) {
    layout.eligible_idx_perm[jp] = elig[perm[jp]];
  }

  // after the permutation, changing columns are exactly positions
  // [n_fixed, n_elig); loops derive local column indices as n_fixed + ci.
  for (std::size_t ci = 0; ci < n_mc; ++ci) {
    layout.miss_rows_offsets[ci + 1] =
        layout.miss_rows_offsets[ci] + col_nmiss[perm[n_fixed + ci]];
  }
  const std::size_t n_missing = layout.miss_rows_offsets[n_mc];

  // `loc` is an R matrix, whose dims are int.
  if (n_missing > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    Rcpp::stop("the columns to impute hold %d missing values, but the result "
               "has a row for each and an R matrix holds at most %d rows. "
               "Impute the columns in groups with group_imp()",
               n_missing, std::numeric_limits<int>::max());
  }
  // scan_columns() fills it.
  layout.miss_rows_flat.assign(n_missing, 0);
  return layout;
}

// pass 2, over the columns of obj in the order of the layout: each is copied
// into its column of Xhat, 0 at a hole, the rows of its holes go to
// layout.miss_rows_flat, and its observed weight and moments to mom.
void scan_columns(const double *objp, const double *wptr, ColumnLayout &layout,
                  double *Xhat, ColumnMoments &mom) {
  const std::size_t nrX = layout.nrX;
  const std::size_t n_fixed = layout.n_fixed;
  const std::size_t n_elig = layout.eligible_idx_perm.size();
  for (std::size_t jp = 0; jp < n_elig; ++jp) {
    const double *src = col_ptr(
        objp, static_cast<std::size_t>(layout.eligible_idx_perm[jp]), nrX);
    double *dst = col_ptr(Xhat, jp, nrX);

    // pass1_count_missing read this column: a complete one (jp < n_fixed)
    // has no NaN, and a changing one has as many as its slice of
    // miss_rows_flat holds.
    std::size_t p =
        (jp >= n_fixed) ? layout.miss_rows_offsets[jp - n_fixed] : 0;

    // pass A: observed weight, raw weighted sum, missing-row indices.
    double d_acc = 0.0, ws_acc = 0.0;
    for (std::size_t i = 0; i < nrX; ++i) {
      const double v = src[i];
      if (std::isnan(v)) {
        layout.miss_rows_flat[p++] = static_cast<int>(i);
      } else {
        const double wi = wptr[i];
        d_acc += wi;
        ws_acc += wi * v;
      }
    }
    mom.denom[jp] = d_acc;

    // per-column reference mean for centered moment tracking. The
    // d_acc == 0 case is rejected right after this loop; 0.0 just keeps
    // the pass-B arithmetic finite until then.
    const double mu0_j = (d_acc > 0.0) ? (ws_acc / d_acc) : 0.0;
    mom.mu0[jp] = mu0_j;

    // pass B: centered moments at O(sd) magnitude (avoids the
    // E[x^2] - mu^2 catastrophic cancellation), plus the Xhat write.
    double osc_acc = 0.0, ossq_acc = 0.0;
    for (std::size_t i = 0; i < nrX; ++i) {
      const double v = src[i];
      if (std::isnan(v)) {
        dst[i] = 0.0;
      } else {
        dst[i] = v;
        const double vc = v - mu0_j;
        const double wi = wptr[i];
        osc_acc += wi * vc;
        ossq_acc += wi * vc * vc;
      }
    }
    mom.obs_sum_c[jp] = osc_acc; // ~0 by construction; kept for exactness
    mom.obs_sumsq_c[jp] = ossq_acc;
  }
}

// The column means and, when scaling, the column scales, from the moments
// scan_columns() took; then Xhat is centered, and scaled, in place.
void standardize_columns(double *Xhat, const ColumnLayout &layout,
                         const bool scale, ColumnMoments &mom) {
  const std::size_t nrX = layout.nrX;
  const std::size_t n_elig = layout.eligible_idx_perm.size();

  // every eligible column must have at least one observed row with positive
  // weight, else mean_p (and the downstream Gram) is NaN.
  for (std::size_t j = 0; j < n_elig; ++j) {
    if (!(mom.denom[j] > 0.0)) {
      Rcpp::stop("an eligible column has zero observed weight (column %d)",
                 layout.eligible_idx_perm[j] + 1);
    }
  }

  for (std::size_t j = 0; j < n_elig; ++j) {
    const double mean_resid = mom.obs_sum_c[j] / mom.denom[j];
    mom.mean_p[j] = mom.mu0[j] + mean_resid;
    if (scale) {
      const double var =
          mom.obs_sumsq_c[j] / mom.denom[j] - mean_resid * mean_resid;
      const double e = std::sqrt(var < 0.0 ? 0.0 : var);
      mom.et[j] = e < PCA_TOL ? PCA_TOL : e;
    }
  }

  if (scale) {
    for (std::size_t j = 0; j < n_elig; ++j) {
      double *col = col_ptr(Xhat, j, nrX);
      const double mu = mom.mean_p[j];
      const double inv_et = 1.0 / mom.et[j];

      for (std::size_t i = 0; i < nrX; ++i) {
        col[i] = (col[i] - mu) * inv_et;
      }
    }
  } else {
    for (std::size_t j = 0; j < n_elig; ++j) {
      double *col = col_ptr(Xhat, j, nrX);
      const double mu = mom.mean_p[j];

      for (std::size_t i = 0; i < nrX; ++i) {
        col[i] -= mu;
      }
    }
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
void impute_restandardize(double *xptr, const double *fxp, ColumnMoments &mom,
                          const double *wptr, const double *swptr,
                          double *X_chg_scaled, const ColumnLayout &layout) {
  const std::size_t nrX = layout.nrX;
  const std::size_t n_fixed = layout.n_fixed;
  const std::size_t n_mc = layout.n_mc;
  double *mean_p = mom.mean_p.data();
  [[maybe_unused]] double *et = mom.et.data();
  const double *mu0 = mom.mu0.data();
  const double *obs_sum_c = mom.obs_sum_c.data();
  [[maybe_unused]] const double *obs_sumsq_c = mom.obs_sumsq_c.data();
  for (std::size_t cidx = 0; cidx < n_mc; ++cidx) {
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
    if constexpr (Scale) {
      old_et = et[j];
    } else {
      old_et = 1.0;
    }

    double imp_sum_c = 0.0;
    [[maybe_unused]] double imp_sumsq_c = 0.0; // accumulated only when Scale

    for_each_missing_in_col(layout, cidx, [&](std::size_t i) {
      const double v_std = fit_col[i];
      col[i] = v_std;

      const double w = wptr[i];

      if constexpr (Scale) {
        // centered imputed value v - mu0[j], formed
        // directly at O(sd); never materialize the
        // O(mean) value v itself.
        const double v_c = v_std * old_et + dmu_old;
        imp_sum_c += w * v_c;
        imp_sumsq_c += w * v_c * v_c;
      } else {
        const double v_c = v_std + dmu_old;
        imp_sum_c += w * v_c;
      }
    });

    // total weight over observed + imputed rows is exactly 1 (row_w was
    // normalized), so the centered sums add up to the full mean offset.
    const double dmu_new = obs_sum_c[j] + imp_sum_c;
    const double new_mu = mu0[j] + dmu_new;

    if constexpr (Scale) {
      const double new_var = obs_sumsq_c[j] + imp_sumsq_c - dmu_new * dmu_new;
      const double new_et = std::sqrt(std::max(0.0, new_var));
      const double new_et_eff = std::max(new_et, PCA_TOL);

      const double inv_new = 1.0 / new_et_eff;
      const double a = old_et * inv_new;
      // (old_mu - new_mu) * inv_new. The mu0 offset cancels.
      const double b = (dmu_old - dmu_new) * inv_new;

      if (scaled_col) {
        for (std::size_t i = 0; i < nrX; ++i) {
          const double z = col[i] * a + b;
          col[i] = z;
          scaled_col[i] = z * swptr[i];
        }
      } else {
        for (std::size_t i = 0; i < nrX; ++i) {
          col[i] = col[i] * a + b;
        }
      }

      et[j] = new_et_eff;
    } else {
      // old_mu - new_mu. The mu0 offset cancels.
      const double b = dmu_old - dmu_new;

      if (scaled_col) {
        for (std::size_t i = 0; i < nrX; ++i) {
          const double z = col[i] + b;
          col[i] = z;
          scaled_col[i] = z * swptr[i];
        }
      } else {
        for (std::size_t i = 0; i < nrX; ++i) {
          col[i] += b;
        }
      }
    }

    mean_p[j] = new_mu;
  }
}

// Starting values at the missing positions of Xhat and of fittedX, which the
// caller hands over zero-filled: before the first dgemm only those positions
// of fittedX are read. A start of zero, init == 0, is written to Xhat only,
// since fittedX already holds it.
void initialize_missing_start(double *Xhat, double *fittedX, const int init,
                              const ColumnLayout &layout) {
  const std::size_t nrX = layout.nrX;
  const std::size_t n_fixed = layout.n_fixed;
  const std::size_t n_mc = layout.n_mc;
  if (n_mc == 0) {
    return;
  }

  if (init == 0) {
    for (std::size_t ci = 0; ci < n_mc; ++ci) {
      const std::size_t j = n_fixed + ci;

      double *col = col_ptr(Xhat, j, nrX);

      for_each_missing_in_col(layout, ci, [&](std::size_t i) { col[i] = 0.0; });
    }

    return;
  }

  // the caller has refused more missing cells than an int counts.
  Rcpp::NumericVector rand_vals =
      Rcpp::rnorm(static_cast<int>(layout.miss_rows_flat.size()));

  for (std::size_t ci = 0; ci < n_mc; ++ci) {
    const std::size_t j = n_fixed + ci;

    double *col = col_ptr(Xhat, j, nrX);
    double *fit_col = col_ptr(fittedX, j, nrX);

    std::size_t p = layout.miss_rows_offsets[ci];

    for_each_missing_in_col(layout, ci, [&](std::size_t i) {
      const double z = rand_vals[p++];
      col[i] = z;
      fit_col[i] = z;
    });
  }
}

double missing_weighted_residual_ss(const double *xhp, const double *fxp,
                                    const double *wptr,
                                    const ColumnLayout &layout) {
  const std::size_t nrX = layout.nrX;
  const std::size_t n_fixed = layout.n_fixed;
  const std::size_t n_mc = layout.n_mc;
  double out = 0.0;

  for (std::size_t ci = 0; ci < n_mc; ++ci) {
    const std::size_t j = n_fixed + ci;
    const std::size_t off = j * nrX;

    const double *xh = xhp + off;
    const double *fx = fxp + off;

    double s = 0.0;

    for_each_missing_in_col(layout, ci, [&](std::size_t i) {
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
void reconstruct_cols(const double *U, const std::size_t u_rows,
                      const std::size_t u_cols, const double *V,
                      const std::size_t v_rows, double *dst,
                      const std::size_t j0, const std::size_t nj) {
  if (nj == 0) {
    return;
  }
  const int st =
      si_la::gemm('N', 'T', u_rows, nj, u_cols, 1.0, U, u_rows, V + j0, v_rows,
                  0.0, col_ptr(dst, j0, u_rows), u_rows);
  stop_if_gemm_refused(st);
}

// Restores, when it goes out of scope, the precision Rcpp::Rcout had when it
// was made.
struct PrecGuard {
  std::streamsize old;

  PrecGuard() : old(Rcpp::Rcout.precision()) {}

  ~PrecGuard() { Rcpp::Rcout.precision(old); }
};

// The trace of the EM loop on Rcpp::Rcout, in columns of these widths.
// start_trace() sets the precision, for a PrecGuard made before it to
// restore, and prints the header; print_iter() prints one iteration.
const int iter_w = 4;
const int obj_w = 8;
const int ratio_w = 12;

void start_trace() {
  Rcpp::Rcout.precision(6);
  Rcpp::Rcout << "iter : " << std::right << std::setw(obj_w) << "obj"
              << " " << std::right << std::setw(ratio_w) << "ratio"
              << std::endl;
}

void print_iter(std::size_t it, double obj_val, double crit) {
  Rcpp::Rcout << std::setw(iter_w) << it << " : " << std::setw(obj_w) << obj_val
              << " " << std::setw(ratio_w) << crit << std::endl;
}

#if PCA_IMP_DIAGNOSTICS
// The diagnostics build only. V_prev and V_new are n_gram x n_cp, and
// V_cross (n_cp x n_cp) and V_cross_sv (n_cp) are scratch the caller sized
// once.
double subspace_cosine(const std::vector<double> &V_prev, const double *V_new,
                       const std::size_t n_gram, const std::size_t n_cp,
                       std::vector<double> &V_cross,
                       std::vector<double> &V_cross_sv) {
  // the smallest singular value of V_prev' * V_new: the cosine of the
  // largest angle between the two subspaces.
  int st = si_la::gemm('T', 'N', n_cp, n_cp, n_gram, 1.0, V_prev.data(), n_gram,
                       V_new, n_gram, 0.0, V_cross.data(), n_cp);
  stop_if_gemm_refused(st);
  st = si_la::svd_values(n_cp, n_cp, V_cross.data(), V_cross_sv.data());
  if (st > 0 || st == si_la::err_nonfinite) {
    // dgesdd did not converge, or the product holds NaN or Inf.
    Rcpp::stop("pca_imp diagnostics: the SVD of the subspace product failed");
  }
  if (st == si_la::err_alloc) {
    Rcpp::stop(
        "pca_imp diagnostics: could not allocate the workspace of dgesdd");
  }
  if (st != si_la::ok) {
    Rcpp::stop("internal error: dgesdd refused its arguments (status %d)", st);
  }
  return *std::min_element(V_cross_sv.begin(), V_cross_sv.end());
}
#endif

// The fixed block of fittedX, which the EM loop skips, completed for the
// SSE; the changing columns of fittedX hold the final reconstruction from
// the last hot-loop DGEMM (U carries the isw/lambda_shrinked scaling).
// In the wide case SVD_triplet() also left the fixed rows of V unwritten:
// they come from the last Vn and the fixed columns of Xhat, which nothing
// has written since pass2_scan. The tall case returns V complete.
void complete_fixed_block(const double *Xhat, double *fittedX,
                          SolveBuffers &buf, const GramDims &gram_dims,
                          const std::size_t n_fixed, const std::size_t n_cp) {
  const std::size_t nrX = gram_dims.nr;
  const std::size_t n_elig = gram_dims.nc;
  const std::size_t n_gram = gram_dims.n;
  if (!gram_dims.tall) {
    const int st =
        pca_detail::gemm_tn_block(Xhat, nrX, n_elig, 0, n_fixed, buf.Vn.data(),
                                  n_gram, n_cp, buf.V.data());
    stop_if_gemm_refused(st);
  }
  reconstruct_cols(buf.U.data(), nrX, n_cp, buf.V.data(), n_elig, fittedX, 0,
                   n_fixed);
}

// The mean squared residual over the observed cells, the nb.init selection
// statistic: the residual over every cell in the units of obj (times et^2
// per column), less that over the holes.
double observed_mse(const double *xhp, const double *fxp,
                    const ColumnMoments &mom, const ColumnLayout &layout) {
  const std::size_t nrX = layout.nrX;
  const std::size_t n_fixed = layout.n_fixed;
  const std::size_t n_mc = layout.n_mc;
  const std::size_t n_elig = layout.eligible_idx_perm.size();
  const std::size_t n_missing = layout.miss_rows_flat.size();

  double sse = 0.0;
  for (std::size_t j = 0; j < n_elig; ++j) {
    const double ej = mom.scale_of(j);
    const double ej2 = ej * ej;
    const std::size_t off = j * nrX;
    const double *xh = xhp + off;
    const double *fx = fxp + off;
    double s = 0.0;
    for (std::size_t i = 0; i < nrX; ++i) {
      const double d = xh[i] - fx[i];
      s += d * d;
    }
    sse += s * ej2;
  }
  for (std::size_t ci = 0; ci < n_mc; ++ci) {
    const std::size_t j_local = n_fixed + ci;
    const double ej = mom.scale_of(j_local);
    const double ej2 = ej * ej;
    const std::size_t off = j_local * nrX;
    const double *xh = xhp + off;
    const double *fx = fxp + off;
    double s = 0.0;
    for_each_missing_in_col(layout, ci, [&](std::size_t i) {
      const double d = xh[i] - fx[i];
      s += d * d;
    });
    sse -= s * ej2;
  }
  // n_obs >= n_elig >= 1: standardize_columns() refused a column with no
  // observed weight, so every column has an observed cell.
  const std::size_t n_obs = nrX * n_elig - n_missing;
  return sse / static_cast<double>(n_obs);
}

// The imputed cells as they go back to R: `loc`, an integer matrix of each
// one's 1-based row and column in obj, and `value`, in the order of the
// layout's holes.
struct ImputedCells {
  Rcpp::IntegerMatrix loc;
  Rcpp::NumericVector value;
};

ImputedCells imputed_cells(const double *xhp, const ColumnMoments &mom,
                           const ColumnLayout &layout) {
  const std::size_t nrX = layout.nrX;
  const std::size_t n_fixed = layout.n_fixed;
  const std::size_t n_mc = layout.n_mc;
  const std::size_t n_missing = layout.miss_rows_flat.size();

  // n_missing fits an int: lay_out_columns() refused it otherwise.
  Rcpp::IntegerMatrix loc(static_cast<int>(n_missing), 2);
  Rcpp::colnames(loc) = Rcpp::CharacterVector::create("row", "col");
  Rcpp::NumericVector value(static_cast<R_xlen_t>(n_missing));
  int *loc_row = loc.begin();
  int *loc_col = loc_row + n_missing;
  double *val = value.begin();
  std::size_t iv = 0;
  for (std::size_t ci = 0; ci < n_mc; ++ci) {
    const std::size_t j_local = n_fixed + ci;
    const double ej = mom.scale_of(j_local);
    const double mj = mom.mean_p[j_local];
    const int orig_col_1based = layout.eligible_idx_perm[j_local] + 1;

    // missMDA::imputePCA() returns completeObs using the final Xhat at
    // missing positions, not the final fittedX. fittedX is still used
    // for the nb.init selection MSE.
    const double *xh = xhp + j_local * nrX;

    for_each_missing_in_col(layout, ci, [&](std::size_t i) {
      loc_row[iv] = static_cast<int>(i) + 1;
      loc_col[iv] = orig_col_1based;
      val[iv] = xh[i] * ej + mj;
      ++iv;
    });
  }
  return ImputedCells{loc, value};
}

} // namespace

// Internal. `eligible_idx` holds 0-based columns of `obj`. Returns the
// imputed cells as `loc`, an integer matrix of each one's 1-based row and
// column in `obj`, and `value`, in the column-major order of the eligible
// columns with the complete ones dropped.
// [[Rcpp::export]]
Rcpp::List
pca_imp_internal_cpp(const Rcpp::NumericMatrix &obj,
                     const Rcpp::IntegerVector &eligible_idx, const int ncp,
                     const bool scale, const bool regularized,
                     const double threshold, const int init, const int maxiter,
                     const int miniter, const Rcpp::NumericVector &row_w,
                     const double coeff_ridge, const int solver,
                     const int warmup_iters, const double lobpcg_tol,
                     const int lobpcg_maxiter, const int trace_iter = 0) {
  LOC_TIMER_OBJ(pca_kernel);
  LOC_TIC(pca_kernel, "pca_imp_internal_cpp_total");
  // No Inf scan here. pca_imp() calls this once per `nb.init` on an unchanged
  // `obj`, so a guard at this level would re-scan the whole matrix on every
  // initialization. check_finite() in pca_imp() is the gate, and it runs once.

  check_arguments(obj, eligible_idx, ncp, regularized, init, maxiter, miniter,
                  row_w, solver, warmup_iters, lobpcg_maxiter, trace_iter);
  const std::vector<double> w = normalized_row_weights(row_w);

  const std::size_t nrX = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_elig = static_cast<std::size_t>(eligible_idx.size());
  const std::size_t n_cp = static_cast<std::size_t>(ncp);
  const std::size_t max_iter = static_cast<std::size_t>(maxiter);
  const std::size_t min_iter = static_cast<std::size_t>(miniter);
  const std::size_t trace_every = static_cast<std::size_t>(trace_iter);
  const bool lobpcg_requested = (solver != 0);
  const bool tall = (nrX >= n_elig);
  const std::size_t k = regularized ? n_cp + 1 : n_cp;
  const std::size_t n_gram = tall ? n_elig : nrX;

  const double nrX_d = static_cast<double>(nrX);
  const double n_elig_d = static_cast<double>(n_elig);
  const double ncp_d = static_cast<double>(n_cp);

  const double denom_sigma = (n_elig_d - ncp_d) * (nrX_d - 1.0 - ncp_d);
  if (regularized && !(denom_sigma > 0.0)) {
    Rcpp::stop("'regularized' = TRUE requires length('eligible_idx') > 'ncp' "
               "and nrow('obj') - 1 > 'ncp'");
  }

  const double *objp = obj.begin();
  const int *elig = eligible_idx.begin();

  LOC_TIC(pca_kernel, "pass1_count_missing");
  const std::vector<std::size_t> col_nmiss =
      count_missing(objp, elig, nrX, n_elig);
  LOC_TOC(pca_kernel, "pass1_count_missing");

  ColumnLayout layout = lay_out_columns(col_nmiss, elig, nrX);

  std::vector<double> Xhat(nrX * n_elig);
  std::vector<double> fittedX(nrX * n_elig);

  ColumnMoments mom(n_elig, scale);

  LOC_TIC(pca_kernel, "pass2_scan");
  scan_columns(objp, w.data(), layout, Xhat.data(), mom);
  LOC_TOC(pca_kernel, "pass2_scan");

  standardize_columns(Xhat.data(), layout, scale, mom);

  initialize_missing_start(Xhat.data(), fittedX.data(), init, layout);

  // ---------------------------------------------------------------------------
  // SVD / solve buffers.
  // ---------------------------------------------------------------------------
  std::vector<double> sqrt_row_w(nrX);
  std::vector<double> inv_sqrt_row_w(nrX);
  for (std::size_t i = 0; i < nrX; ++i) {
    const double s = std::sqrt(w[i]);
    sqrt_row_w[i] = s;
    inv_sqrt_row_w[i] = (s > PCA_TOL) ? (1.0 / s) : (1.0 / PCA_TOL);
  }

  GramDims gram_dims;
  gram_dims.init(nrX, n_elig, tall);
  EigSymRange eig_range;
  if (!eig_range.init(n_gram, k)) {
    Rcpp::stop("internal error: eigenpair count out of range");
  }

  // every buffer SVD_triplet() reads or writes is sized here, once, by
  // SolveBuffers, with the sizes its header lists, and the hot code reads
  // them through raw pointers. In the wide case V's fixed rows stay
  // unwritten until postprocessing, and X_work is empty.
  LOC_TIC(pca_kernel, "alloc_solve_buffers");
  SolveBuffers buf(gram_dims, n_cp, k);
  std::vector<double> lambda_shrinked(n_cp);
  LOC_TOC(pca_kernel, "alloc_solve_buffers");

  LOC_TIC(pca_kernel, "gram_cache_init");
  GramCache gram_cache;
  const int gram_init_st =
      gram_cache_init(gram_cache, Xhat.data(), sqrt_row_w.data(), gram_dims,
                      layout.n_fixed, buf.X_work.data());
  LOC_TOC(pca_kernel, "gram_cache_init");
  if (gram_init_st != si_la::ok) {
    Rcpp::stop("internal error: dsyrk refused its arguments (status %d)",
               gram_init_st);
  }
  double *weighted_chg_ptr =
      gram_changing_block(gram_cache, buf.X_work.data(), gram_dims);

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
  std::size_t n_iter_final = 0;
  double criterion_final = std::numeric_limits<double>::quiet_NaN();
  bool converged = false;

  PrecGuard prec_guard;

  if (trace_every > 0) {
    start_trace();
  }

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

  double old = std::numeric_limits<double>::infinity();
  for (std::size_t nb_iter = 1;; ++nb_iter) {
    if (nb_iter % INTERRUPT_EVERY == 0) {
      Rcpp::checkUserInterrupt();
    }
    LOC_TIC(pca_kernel, "restandardize");
    if (scale) {
      impute_restandardize<true>(Xhat.data(), fittedX.data(), mom, w.data(),
                                 sqrt_row_w.data(), weighted_chg_ptr, layout);
    } else {
      impute_restandardize<false>(Xhat.data(), fittedX.data(), mom, w.data(),
                                  sqrt_row_w.data(), weighted_chg_ptr, layout);
    }
    LOC_TOC(pca_kernel, "restandardize");

    LOC_TIC(pca_kernel, "svd");
    const double trace_val = SVD_triplet(
        Xhat.data(), sqrt_row_w.data(), n_cp, buf, gram_dims, gram_cache,
        hyb_ctx, static_cast<int>(nb_iter - 1) LOC_TIMER_ARG(pca_kernel));
    LOC_TOC(pca_kernel, "svd");

    LOC_TIC(pca_kernel, "post_svd");
    double sigma2 = 0.0;
    double sum_top_sq = 0.0;
    for (std::size_t cp = 0; cp < n_cp; ++cp) {
      sum_top_sq += buf.vs_top[cp] * buf.vs_top[cp];
    }
    // the tail energy is clamped at zero once, and both sigma2 and the
    // objective use the clamped value.
    const double tail = std::max(0.0, trace_val - sum_top_sq);
    if (regularized) {
      sigma2 = scale_factor * tail / denom_sigma;
      sigma2 =
          std::min(sigma2 * coeff_ridge, buf.vs_top[n_cp] * buf.vs_top[n_cp]);
    }

    double top_corr = 0.0;
    for (std::size_t i = 0; i < n_cp; ++i) {
      const double v = buf.vs_top[i];
      const double lam = (v > PCA_TOL) ? (v - sigma2 / v) : 0.0;
      lambda_shrinked[i] = lam;
      const double d = v - lam;
      top_corr += d * d;
    }
    const double objective_full = tail + top_corr;
    LOC_TOC(pca_kernel, "post_svd");

    LOC_TIC(pca_kernel, "reconstruct");

    // after this, U contains diag(isw) %*% U %*% diag(lambda_shrinked):
    // the row de-weighting of SVD_triplet's weighted-space U fused with the
    // shrunk column scaling.
    pca_detail::scale_rows_cols_inplace(buf.U.data(), inv_sqrt_row_w.data(),
                                        lambda_shrinked.data(), nrX, n_cp);
    // reconstruct only the changing columns of fittedX; the fixed block is
    // never read inside the loop and is filled once after it exits.
    reconstruct_cols(buf.U.data(), nrX, n_cp, buf.V.data(), n_elig,
                     fittedX.data(), layout.n_fixed, layout.n_mc);

    LOC_TOC(pca_kernel, "reconstruct");

    LOC_TIC(pca_kernel, "objective");
    const double miss_contrib = missing_weighted_residual_ss(
        Xhat.data(), fittedX.data(), w.data(), layout);
    const double objective = std::max(0.0, objective_full - miss_contrib);
    LOC_TOC(pca_kernel, "objective");

#if PCA_IMP_DIAGNOSTICS
    std::copy(buf.vs_top.begin(), buf.vs_top.begin() + n_cp,
              eigval_hist.begin() + (nb_iter - 1) * n_cp);
    obj_hist[nb_iter - 1] = objective;
    if (!V_prev.empty()) {
      subspace_cos_hist[nb_iter - 1] = subspace_cosine(
          V_prev, buf.eigvecs.data(), n_gram, n_cp, V_cross, V_cross_sv);
    } else {
      subspace_cos_hist[nb_iter - 1] = std::numeric_limits<double>::quiet_NaN();
    }
    V_prev.assign(buf.eigvecs.begin(), buf.eigvecs.begin() + n_gram * n_cp);
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
    if (nb_iter >= min_iter) {
      if (criterion < threshold || objective < threshold) {
        stop_now = true;
        converged = true;
      }
    }
    if (!stop_now && nb_iter >= max_iter) {
      stop_now = true;
      // converged stays false: the test above sets it only with stop_now.
      // deliberately silent. The std::string overload of Rcpp::warning is a
      // bare Rf_warning, so under options(warn = 2) - or any calling handler
      // that invokes a restart - R longjmps straight out of this frame.
      // END_RCPP catches C++ exceptions, not longjmps, so Xhat, fittedX,
      // X_work, AA_NxN, the eigen workspaces, PrecGuard and the wrapper's
      // RNGScope would all be skipped. `converged` is already in the returned
      // list and pca_imp() warns on it once, for the winning restart only.
    }

    if (trace_every > 0) {
      const bool periodic = (nb_iter == 1) || (nb_iter % trace_every == 0);
      if (periodic || stop_now) print_iter(nb_iter, objective, criterion);
    }

    if (stop_now) {
      n_iter_final = nb_iter;
      // On a converged fit criterion_final is the criterion that stopped it.
      // When the objective alone stopped the fit, the criterion did not, so
      // there is none to report and the value is NA. That can withhold a
      // well-defined ratio; it also covers one that is rounding noise, with
      // no tolerance needed to tell the two apart.
      criterion_final =
          (converged && !(criterion < threshold)) ? NA_REAL : criterion;
      break;
    }
  }

  // ---------------------------------------------------------------------------
  // postprocessing.
  // ---------------------------------------------------------------------------
  LOC_TIC(pca_kernel, "postprocessing");
  complete_fixed_block(Xhat.data(), fittedX.data(), buf, gram_dims,
                       layout.n_fixed, n_cp);
  const double mse = observed_mse(Xhat.data(), fittedX.data(), mom, layout);
  const ImputedCells imputed = imputed_cells(Xhat.data(), mom, layout);
  LOC_TOC(pca_kernel, "postprocessing");
  LOC_TOC(pca_kernel, "pca_imp_internal_cpp_total");
  const int n_iter_out = static_cast<int>(n_iter_final);
#if PCA_IMP_DIAGNOSTICS
  Rcpp::NumericMatrix eigval_hist_out(ncp, n_iter_out);
  std::copy(eigval_hist.begin(), eigval_hist.begin() + n_cp * n_iter_final,
            eigval_hist_out.begin());
  obj_hist.resize(n_iter_final);
  subspace_cos_hist.resize(n_iter_final);

  hyb_ctx.truncate_logs(n_iter_final);

  return Rcpp::List::create(
      Rcpp::Named("loc") = imputed.loc, Rcpp::Named("value") = imputed.value,
      Rcpp::Named("mse") = mse, Rcpp::Named("n_iter") = n_iter_out,
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
      Rcpp::Named("loc") = imputed.loc, Rcpp::Named("value") = imputed.value,
      Rcpp::Named("mse") = mse, Rcpp::Named("n_iter") = n_iter_out,
      Rcpp::Named("criterion_final") = criterion_final,
      Rcpp::Named("converged") = converged,
      Rcpp::Named("solver_chosen") = hyb_ctx.solver_chosen_code(),
      Rcpp::Named("n_exact") = hyb_ctx.n_exact,
      Rcpp::Named("n_lobpcg_ok") = hyb_ctx.n_lobpcg_ok,
      Rcpp::Named("n_lobpcg_bad") = hyb_ctx.n_lobpcg_bad);
#endif
}
