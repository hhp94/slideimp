// [[Rcpp::depends(RcppEigen, dqrng, BH)]]

#ifndef USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <RcppEigen.h>
#include "Spectra/SymEigsSolver.h"
#include <dqrng.h>
#include <dqrng_distribution.h>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>
#include <memory>
#include <chrono>
#include "loc_timer.h"

constexpr int SPECTRA_NCV_FLOOR = 20;

static int spectra_ncv_choose(int n, int k, int ncv_req) {
  if (ncv_req > 0) {
    if (ncv_req <= k || ncv_req > n)
      Rcpp::stop("'ncv' in spectra_opts() must be more than %d and at most "
                 "%d, the smaller dimension of the matrix. It is %d",
                 k, n, ncv_req);
    return ncv_req;
  }
  int ncv = 2 * k + 1;
  if (ncv < SPECTRA_NCV_FLOOR) ncv = SPECTRA_NCV_FLOOR;
  if (ncv > n) ncv = n;
  return ncv;
}

template <typename EigsT>
static int spectra_solve(EigsT &eigs, int n, int k, double tol, int maxitr,
                         double *values, double *vectors, double *start,
                         bool warm LOC_TIMER_PARAM(tm)) {
  LOC_TIC(tm, "s_init");
  if (warm)
    eigs.init(start);
  else
    eigs.init();
  LOC_TOC(tm, "s_init");
  LOC_TIC(tm, "s_compute");
  const int nconv = eigs.compute(Spectra::SortRule::LargestAlge, maxitr, tol,
                                 Spectra::SortRule::LargestAlge);
  LOC_TOC(tm, "s_compute");
  if (nconv < k) return nconv;

  LOC_TIC(tm, "s_vectors");
  const Eigen::VectorXd ev = eigs.eigenvalues();
  const Eigen::MatrixXd evec = eigs.eigenvectors();
  LOC_TOC(tm, "s_vectors");
  LOC_TIC(tm, "s_copy");
  std::copy(ev.data(), ev.data() + k, values);
  std::copy(evec.data(), evec.data() + static_cast<std::ptrdiff_t>(n) * k,
            vectors);
  std::fill(start, start + n, 0.0);
  for (int s = 0; s < k; ++s) {
    const double *col = vectors + static_cast<std::ptrdiff_t>(s) * n;
    for (int i = 0; i < n; ++i) start[i] += col[i];
  }
  LOC_TOC(tm, "s_copy");
  return nconv;
}

#ifdef LOC_TIMER
class TimedSymMatProd {
public:
  using Scalar = double;
  TimedSymMatProd(const Spectra::DenseSymMatProd<double> &op, Rcpp::Timer &t)
      : op_(op), t_(t) {}
  Eigen::Index rows() const { return op_.rows(); }
  Eigen::Index cols() const { return op_.cols(); }
  void perform_op(const double *x_in, double *y_out) const {
    t_.tic("s_matvec");
    op_.perform_op(x_in, y_out);
    t_.toc("s_matvec");
  }

private:
  const Spectra::DenseSymMatProd<double> &op_;
  Rcpp::Timer &t_;
};
#endif

static void sweep_cols(double *x, int n, int nd, double *etd) {
  const double inv_n = 1.0 / static_cast<double>(n);
  for (int c = 0; c < nd; ++c) {
    double *col = x + static_cast<std::ptrdiff_t>(c) * n;
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += col[i];
    const double mu = s * inv_n;
    for (int i = 0; i < n; ++i) col[i] -= mu;
    etd[c] += mu;
  }
}

static void moy_dot(const double *u, const double *v, const double *w,
                    const int *row, const int *col, int n, int nd, int S,
                    R_xlen_t m, double *ut, double *vwt, double *out) {
  for (int s = 0; s < S; ++s) {
    const double *ucol = u + static_cast<std::ptrdiff_t>(s) * n;
    for (int i = 0; i < n; ++i)
      ut[static_cast<std::ptrdiff_t>(i) * S + s] = ucol[i];
    const double ws = w[s];
    const double *vcol = v + static_cast<std::ptrdiff_t>(s) * nd;
    for (int j = 0; j < nd; ++j)
      vwt[static_cast<std::ptrdiff_t>(j) * S + s] = vcol[j] * ws;
  }
  for (R_xlen_t k = 0; k < m; ++k) {
    const double *a = ut + static_cast<std::ptrdiff_t>(row[k] - 1) * S;
    const double *b = vwt + static_cast<std::ptrdiff_t>(col[k] - 1) * S;
    double acc = 0.0;
    for (int s = 0; s < S; ++s) acc += a[s] * b[s];
    out[k] = acc;
  }
}

// [[Rcpp::export(.build_block)]]
void build_block(SEXP src, int src_off, int j1, int nb, int transformed,
                 Rcpp::NumericVector scratch, Rcpp::NumericMatrix Xd,
                 Rcpp::NumericVector g_fixed, Rcpp::NumericVector et,
                 Rcpp::NumericVector col_ss, Rcpp::IntegerVector mi_row,
                 Rcpp::IntegerVector cnt, Rcpp::IntegerVector col_ptr,
                 Rcpp::IntegerVector dirty_pos, int resume) {
  const int n = Xd.nrow();
  const int nd = Xd.ncol();
  const int p = cnt.size();
  if (nb < 1 || j1 < 1 || j1 + nb - 1 > p)
    Rcpp::stop("internal error: the block columns [%d, %d] fall outside 1..%d",
               j1, j1 + nb - 1, p);
  if (TYPEOF(src) != REALSXP)
    Rcpp::stop("internal error: 'src' must be a double matrix "
               "(storage.mode \"double\")");
  if (n < 1)
    Rcpp::stop("internal error: 'Xd' has no rows; the build needs at least "
               "one");
  if (src_off < 0 || Rf_xlength(src) % n != 0 ||
      Rf_xlength(src) < static_cast<R_xlen_t>(n) * (src_off + nb))
    Rcpp::stop("internal error: 'src' must have %d rows and at least %d "
               "columns; it has %.0f cells",
               n, src_off + nb, static_cast<double>(Rf_xlength(src)));
  if (et.size() != p || col_ss.size() != p || dirty_pos.size() != p ||
      col_ptr.size() != p + 1)
    Rcpp::stop("internal error: 'et', 'col_ss', 'dirty_pos' must have length "
               "p and 'col_ptr' p + 1");
  const bool has_g = g_fixed.size() > 0;
  if (resume && has_g)
    Rcpp::stop("internal error: a resume build must not accumulate a Gram; "
               "'g_fixed' comes off the saved state");
  if (has_g && g_fixed.size() != static_cast<R_xlen_t>(n) * n)
    Rcpp::stop("internal error: 'g_fixed' must be empty or n x n");
  if (has_g && scratch.size() < static_cast<R_xlen_t>(n) * nb)
    Rcpp::stop("internal error: 'scratch' must have at least %d x %d cells; "
               "it has %.0f",
               n, nb, static_cast<double>(scratch.size()));
  if (mi_row.size() != static_cast<R_xlen_t>(col_ptr[p]) - 1)
    Rcpp::stop("internal error: 'mi_row' must have col_ptr[p + 1] - 1 "
               "entries");

  const double *b = REAL(src) + static_cast<R_xlen_t>(src_off) * n;
  double *sc = has_g ? REAL(scratch) : nullptr;
  double *xd = REAL(Xd);
  int *rows_all = INTEGER(mi_row);
  int kc = 0;
  for (int jj = 0; jj < nb; ++jj) {
    const int j = j1 - 1 + jj;
    if (resume && dirty_pos[j] == 0) {
      if (cnt[j] != 0)
        Rcpp::stop("internal error: a resume build has no slot in 'Xd' for "
                   "column %d, which has holes",
                   j + 1);
      continue;
    }
    const double *col = b + static_cast<R_xlen_t>(jj) * n;

    double s = 0.0;
    R_xlen_t nobs = 0;
    int holes = 0;
    for (int i = 0; i < n; ++i) {
      const double x = col[i];
      if (ISNAN(x))
        ++holes;
      else {
        if (!R_FINITE(x)) {
          if (transformed)
            Rcpp::stop("the transform maps column %d, row %d to a non-finite "
                       "value. A transform must give a finite value for every "
                       "observed cell",
                       j + 1, i + 1);
          Rcpp::stop("'obj' contains non-finite values");
        }
        s += x;
        ++nobs;
      }
    }
    if (holes != cnt[j]) {
      if (transformed && holes > cnt[j])
        Rcpp::stop("the transform maps observed values in column %d to NA or "
                   "NaN (%d of them). A transform must give a finite value "
                   "for every observed cell",
                   j + 1, holes - cnt[j]);
      if (transformed)
        Rcpp::stop("the transform maps missing cells in column %d to a value "
                   "(%d of them). A transform must map NA to NA",
                   j + 1, cnt[j] - holes);
      Rcpp::stop("internal error: the scan and the block disagree on the "
                 "holes in column %d: %d and %d",
                 j + 1, cnt[j], holes);
    }
    const double mu = s / static_cast<double>(nobs);
    et[j] = mu;

    const int dp = dirty_pos[j];
    if (dp < 0 || dp > nd)
      Rcpp::stop("internal error: the slot of column %d is %d, outside 1..%d",
                 j + 1, dp, nd);
    if (holes > 0 && dp == 0)
      Rcpp::stop("internal error: a build has no slot in 'Xd' for column %d, "
                 "which has holes",
                 j + 1);
    if (dp == 0 && !has_g)
      Rcpp::stop("internal error: no slot in 'Xd' and no 'g_fixed' to fold "
                 "column %d into",
                 j + 1);
    double *dst = (dp > 0) ? xd + static_cast<R_xlen_t>(dp - 1) * n
                           : sc + static_cast<R_xlen_t>(kc) * n;
    int *rows = rows_all + (col_ptr[j] - 1);
    int k = 0;
    double ss = 0.0;
    for (int i = 0; i < n; ++i) {
      const double x = col[i];
      double v;
      if (ISNAN(x)) {
        v = 0.0;
        rows[k++] = i + 1;
      } else
        v = x - mu;
      dst[i] = v;
      ss += v * v;
    }
    col_ss[j] = ss;
    if (dp == 0) ++kc;
  }
  if (kc > 0) {
    const char uplo = 'L', trans = 'N';
    const double one = 1.0;
    F77_CALL(dsyrk)
    (&uplo, &trans, &n, &kc, &one, sc, &n, &one, REAL(g_fixed), &n FCONE FCONE);
  }
}

// [[Rcpp::export(.sym_from_lower)]]
void sym_from_lower(Rcpp::NumericMatrix g) {
  const int n = g.nrow();
  if (g.ncol() != n) Rcpp::stop("internal error: 'g' must be square");
  double *x = REAL(g);
  for (int c = 0; c < n; ++c)
    for (int r = c + 1; r < n; ++r)
      x[static_cast<R_xlen_t>(r) * n + c] = x[static_cast<R_xlen_t>(c) * n + r];
}

struct da_ctx {
  int n, nd, ng;
  bool tall, has_fixed;
  R_xlen_t nmiss, nstore;
  double *xp;
  double *dp;
  const double *gfp;
  const int *midx, *mrow, *mcol, *sidx, *scol;
  double sigma;
  const double *etd_in;
  bool resuming;
  const double *xtild_in;
  const double *eig_start_in;
  const double *disp;
  double *etd_out;
  double *xtild_out;
  double *eig_start_out;
  int S, K;
  int p;
  double dfP;
  double nu;
  int warmup, ndraws, thin;
  double gram_max_cond, spectra_tol;
  int ncv, spectra_maxitr;
  int chain_id;
  int refresh;
};

static double da_chain_run(const da_ctx &cx,
                           dqrng::random_64bit_generator &rng_eng);

static da_ctx da_ctx_build(int n, int nd, Rcpp::NumericVector g_fixed,
                           Rcpp::IntegerVector miss_idx_l,
                           Rcpp::IntegerVector mi_row,
                           Rcpp::IntegerVector mi_col_l,
                           Rcpp::IntegerVector store_idx,
                           Rcpp::IntegerVector mi_col_store, int S, int K,
                           int p, double dfP, double nu, int warmup, int ndraws,
                           int thin, double gram_max_cond, double spectra_tol,
                           int spectra_ncv, int spectra_maxitr, int refresh) {
  const R_xlen_t nmiss = miss_idx_l.size();
  const R_xlen_t nstore = store_idx.size();
  const bool has_fixed = g_fixed.size() > 0;

  if (nd < 1)
    Rcpp::stop("internal error: 'Xd' has no columns, so there is nothing to "
               "impute");

  const bool tall = !has_fixed && n > nd;
  const int ng = tall ? nd : n;

  const int k_solve = S + 1;
  if (S < 1 || k_solve >= ng)
    Rcpp::stop("internal error: the core solves 'S' + 1 eigenpairs on the "
               "Gram, which needs 'S' + 2 <= %d (the smaller dimension); got "
               "'S' = %d",
               ng, S);
  if (mi_row.size() != nmiss || mi_col_l.size() != nmiss)
    Rcpp::stop("internal error: 'mi_row', 'mi_col_l' and 'miss_idx_l' must "
               "agree");
  if (mi_col_store.size() != nstore)
    Rcpp::stop("internal error: 'store_idx' and 'mi_col_store' must agree");
  if (has_fixed &&
      g_fixed.size() != static_cast<R_xlen_t>(n) * static_cast<R_xlen_t>(n))
    Rcpp::stop("internal error: 'g_fixed' must be n x n or empty");
  if (dfP <= 0.0) Rcpp::stop("internal error: 'dfP' must be positive");
  if (!(nu > 2.0))
    Rcpp::stop("internal error: 'nu' must be > 2 (Inf for Gaussian noise)");
  if (K < 1) Rcpp::stop("internal error: 'K' must be positive");
  if (p < nd)
    Rcpp::stop("internal error: 'p' must be at least the dirty width 'nd'");
  if (warmup < 0) Rcpp::stop("internal error: 'warmup' must be >= 0");
  if (ndraws < 1) Rcpp::stop("internal error: 'ndraws' must be >= 1");
  if (thin < 1) Rcpp::stop("internal error: 'thin' must be >= 1");
  if (ndraws > (std::numeric_limits<int>::max() - warmup) / thin)
    Rcpp::stop("'iter_warmup' + 'iter_sampling' x 'thin' must be at most "
               "%d. Lower one of them",
               std::numeric_limits<int>::max());
  if (refresh < 0)
    Rcpp::stop("internal error: 'refresh' must be >= 0 (0 = no progress "
               "lines)");
  if (!(spectra_tol > 0.0))
    Rcpp::stop("'tol' in spectra_opts() must be a positive number. It is %g",
               spectra_tol);
  if (spectra_maxitr < 1)
    Rcpp::stop("'maxitr' in spectra_opts() must be at least 1. It is %d",
               spectra_maxitr);

  {
    const R_xlen_t ncell = static_cast<R_xlen_t>(n) * nd;
    for (R_xlen_t k = 0; k < nmiss; ++k) {
      if (miss_idx_l[k] < 1 || miss_idx_l[k] > ncell)
        Rcpp::stop("internal error: 'miss_idx_l' out of range (expects "
                   "1-based)");
      if (mi_row[k] < 1 || mi_row[k] > n)
        Rcpp::stop("internal error: 'mi_row' out of range (expects 1-based)");
      if (mi_col_l[k] < 1 || mi_col_l[k] > nd)
        Rcpp::stop("internal error: 'mi_col_l' out of range (expects "
                   "1-based)");
    }
    for (R_xlen_t j = 0; j < nstore; ++j) {
      if (store_idx[j] < 1 || store_idx[j] > nmiss)
        Rcpp::stop("internal error: 'store_idx' out of range (expects "
                   "1-based)");
      if (mi_col_store[j] < 1 || mi_col_store[j] > nd)
        Rcpp::stop("internal error: 'mi_col_store' out of range (expects "
                   "1-based)");
    }
  }

  da_ctx cx = {};
  cx.n = n;
  cx.nd = nd;
  cx.ng = ng;
  cx.tall = tall;
  cx.has_fixed = has_fixed;
  cx.nmiss = nmiss;
  cx.nstore = nstore;
  cx.gfp = has_fixed ? REAL(g_fixed) : NULL;
  cx.midx = INTEGER(miss_idx_l);
  cx.mrow = INTEGER(mi_row);
  cx.mcol = INTEGER(mi_col_l);
  cx.sidx = INTEGER(store_idx);
  cx.scol = INTEGER(mi_col_store);
  cx.S = S;
  cx.K = K;
  cx.p = p;
  cx.dfP = dfP;
  cx.nu = nu;
  cx.warmup = warmup;
  cx.ndraws = ndraws;
  cx.thin = thin;
  cx.gram_max_cond = gram_max_cond;
  cx.spectra_tol = spectra_tol;
  cx.ncv = spectra_ncv_choose(ng, k_solve, spectra_ncv);
  cx.spectra_maxitr = spectra_maxitr;
  cx.chain_id = 1;
  cx.refresh = refresh;
  return cx;
}

static void need_double(SEXP x, const char *what, int chain) {
  if (TYPEOF(x) != REALSXP)
    Rcpp::stop("internal error: a chain input must be of type double: '%s' of "
               "chain %d",
               what, chain);
}

static void need_list(SEXP x, const char *what) {
  if (TYPEOF(x) != VECSXP)
    Rcpp::stop("internal error: '%s' must be a list", what);
}

// The per-chain lists arrive as SEXP and are read with R's own API
// (dev/CPP.md, "No Rcpp::List built from a SEXP").
// [[Rcpp::export(.mipca_da_loop_multi)]]
Rcpp::List mipca_da_loop_multi(
    SEXP xd_list, Rcpp::NumericVector g_fixed, SEXP etd_list, SEXP xtild_list,
    SEXP eig_start_list, Rcpp::NumericVector sigma_in, Rcpp::NumericVector disp,
    bool resuming, Rcpp::CharacterVector rng_state_in,
    Rcpp::IntegerVector miss_idx_l, Rcpp::IntegerVector mi_row,
    Rcpp::IntegerVector mi_col_l, Rcpp::IntegerVector store_idx,
    Rcpp::IntegerVector mi_col_store, SEXP draws_list,
    Rcpp::IntegerVector streams, int S, int K, int p, double dfP, double nu,
    int warmup, int ndraws, int thin, double gram_max_cond, double spectra_tol,
    int spectra_ncv, int spectra_maxitr, int refresh) {
  need_list(xd_list, "xd_list");
  need_list(etd_list, "etd_list");
  need_list(xtild_list, "xtild_list");
  need_list(eig_start_list, "eig_start_list");
  need_list(draws_list, "draws_list");
  const int nc = static_cast<int>(Rf_xlength(xd_list));
  if (nc < 1)
    Rcpp::stop("internal error: the sampler needs at least one chain");
  if (Rf_xlength(draws_list) != nc || streams.size() != nc)
    Rcpp::stop("internal error: 'xd_list', 'draws_list' and 'streams' must "
               "have one entry per chain");
  if (Rf_xlength(etd_list) != nc || Rf_xlength(xtild_list) != nc ||
      Rf_xlength(eig_start_list) != nc || sigma_in.size() != nc)
    Rcpp::stop("internal error: 'etd_list', 'xtild_list', 'eig_start_list' "
               "and 'sigma_in' must have one entry per chain");
  if (resuming) {
    if (rng_state_in.size() != nc)
      Rcpp::stop("internal error: a resume must ship one serialized rng state "
                 "per chain; got %d states for %d chains",
                 static_cast<int>(rng_state_in.size()), nc);
  } else {
    if (rng_state_in.size() > 0)
      Rcpp::stop("internal error: a fresh run must ship no serialized rng "
                 "state; got %d",
                 static_cast<int>(rng_state_in.size()));
  }

  Rcpp::NumericMatrix xd0(VECTOR_ELT(xd_list, 0));
  const int n = xd0.nrow();
  const int nd = xd0.ncol();

  da_ctx shared = da_ctx_build(n, nd, g_fixed, miss_idx_l, mi_row, mi_col_l,
                               store_idx, mi_col_store, S, K, p, dfP, nu,
                               warmup, ndraws, thin, gram_max_cond, spectra_tol,
                               spectra_ncv, spectra_maxitr, refresh);

  if (resuming) {
    if (disp.size() > 0)
      Rcpp::stop("internal error: a dispersed init cannot be applied to a "
                 "resume; the saved state already holds where each chain "
                 "stopped");
  } else {
    if (disp.size() != nd)
      Rcpp::stop("internal error: a fresh run needs one dispersion scale per "
                 "dirty column (%d); got %d",
                 nd, static_cast<int>(disp.size()));
    for (int j = 0; j < nd; ++j)
      if (!R_finite(disp[j]) || disp[j] < 0.0)
        Rcpp::stop("internal error: 'disp' must be finite and non-negative; "
                   "entry %d is not",
                   j + 1);
    shared.disp = REAL(disp);
  }
  shared.resuming = resuming;

  if (!resuming) {
    for (int a = 0; a < nc; ++a) {
      if (streams[a] < 0)
        Rcpp::stop("internal error: 'streams' must be non-negative");
      for (int b = a + 1; b < nc; ++b)
        if (streams[a] == streams[b])
          Rcpp::stop("internal error: 'streams' must be unique; chains %d and "
                     "%d both got stream %d",
                     a + 1, b + 1, streams[a]);
    }
  }

  std::vector<da_ctx> ctxs(nc, shared);
  Rcpp::List etd_out_l(nc), xtild_out_l(nc), eig_start_out_l(nc);
  for (int c = 0; c < nc; ++c) {
    SEXP xd_c = VECTOR_ELT(xd_list, c);
    SEXP draws_c = VECTOR_ELT(draws_list, c);
    SEXP etd_c = VECTOR_ELT(etd_list, c);
    SEXP xtild_c = VECTOR_ELT(xtild_list, c);
    SEXP eig_start_c = VECTOR_ELT(eig_start_list, c);
    need_double(xd_c, "Xd", c + 1);
    need_double(draws_c, "draws", c + 1);
    need_double(etd_c, "etd", c + 1);
    need_double(xtild_c, "xtild_miss", c + 1);
    need_double(eig_start_c, "eig_start", c + 1);
    Rcpp::NumericMatrix xdc(xd_c);
    Rcpp::NumericMatrix drc(draws_c);
    Rcpp::NumericVector ec(etd_c);
    Rcpp::NumericVector xc(xtild_c);
    Rcpp::NumericVector sc(eig_start_c);
    if (xdc.nrow() != n || xdc.ncol() != nd)
      Rcpp::stop("internal error: 'Xd' of chain %d must be %d x %d like chain "
                 "1's",
                 c + 1, n, nd);
    if (drc.nrow() != shared.nstore || drc.ncol() != ndraws)
      Rcpp::stop("internal error: 'draws' of chain %d must be nstore x "
                 "ndraws",
                 c + 1);
    if (ec.size() != nd)
      Rcpp::stop("internal error: 'etd' of chain %d must have one entry per "
                 "column of 'Xd'",
                 c + 1);
    if (xc.size() != (resuming ? shared.nmiss : 0))
      Rcpp::stop(resuming
                     ? "internal error: 'xtild_miss' of chain %d must have one "
                       "entry per hole"
                     : "internal error: 'xtild_miss' of chain %d must be empty "
                       "on a fresh run; the chain draws its own start from "
                       "'disp'",
                 c + 1);
    if (sc.size() != (resuming ? shared.ng : 0))
      Rcpp::stop(resuming
                     ? "internal error: 'eig_start' of chain %d must have one "
                       "entry per row of the Gram"
                     : "internal error: 'eig_start' of chain %d must be empty "
                       "on a fresh run; the chain's first eigensolve starts "
                       "cold",
                 c + 1);
    for (R_xlen_t i = 0; i < sc.size(); ++i)
      if (!R_finite(sc[i]))
        Rcpp::stop("internal error: 'eig_start' of chain %d must be finite",
                   c + 1);
    ctxs[c].chain_id = c + 1;
    ctxs[c].xp = REAL(xdc);
    ctxs[c].dp = REAL(drc);
    ctxs[c].etd_in = REAL(ec);
    ctxs[c].xtild_in = resuming ? REAL(xc) : NULL;
    ctxs[c].eig_start_in = resuming ? REAL(sc) : NULL;
    ctxs[c].sigma = sigma_in[c];
    Rcpp::NumericVector eo(nd);
    Rcpp::NumericVector xo(shared.nmiss);
    Rcpp::NumericVector so(shared.ng);
    etd_out_l[c] = eo;
    xtild_out_l[c] = xo;
    eig_start_out_l[c] = so;
    ctxs[c].etd_out = REAL(eo);
    ctxs[c].xtild_out = REAL(xo);
    ctxs[c].eig_start_out = REAL(so);
  }

  typedef dqrng::random_64bit_generator::result_type rng_result_t;
  std::vector<std::unique_ptr<dqrng::random_64bit_generator>> gens;
  gens.reserve(nc);
  {
    dqrng::random_64bit_accessor acc;
    for (int c = 0; c < nc; ++c) {
      if (!resuming) {
        gens.push_back(acc.clone(static_cast<rng_result_t>(streams[c])));
      } else {
        std::unique_ptr<dqrng::random_64bit_generator> g = acc.clone(0);
        std::istringstream is(Rcpp::as<std::string>(rng_state_in[c]));
        is >> *g;
        if (is.fail())
          Rcpp::stop("update() cannot extend this fit: the saved random "
                     "number state of chain %d cannot be read. Fit again with "
                     "mipca() and 'keep_state = TRUE'",
                     c + 1);
        gens.push_back(std::move(g));
      }
    }
  }

#ifdef LOC_TIMER
  if (nc > 1)
    Rcpp::stop("an instrumented (LOC_TIMER) build refuses more than one "
               "chain; profiling is single-chain only");
#endif
  Rcpp::NumericVector sigma_out(nc);
  for (int c = 0; c < nc; ++c) sigma_out[c] = da_chain_run(ctxs[c], *gens[c]);

  Rcpp::CharacterVector rng_state(nc);
  for (int c = 0; c < nc; ++c) {
    std::ostringstream os;
    os << *gens[c];
    rng_state[c] = os.str();
  }

  return Rcpp::List::create(Rcpp::Named("sigma") = sigma_out,
                            Rcpp::Named("etd") = etd_out_l,
                            Rcpp::Named("xtild_miss") = xtild_out_l,
                            Rcpp::Named("eig_start") = eig_start_out_l,
                            Rcpp::Named("rng_state") = rng_state);
}

static inline double mt_chisq(dqrng::random_64bit_generator &g,
                              dqrng::normal_distribution &z, double d,
                              double c) {
  for (;;) {
    double x, v;
    do {
      x = z(g);
      v = 1.0 + c * x;
    } while (v <= 0.0);
    v = v * v * v;
    const double u = g.uniform01();
    const double x2 = x * x;
    if (u < 1.0 - 0.0331 * x2 * x2) return 2.0 * d * v;
    if (std::log(u) < 0.5 * x2 + d * (1.0 - v + std::log(v)))
      return 2.0 * d * v;
  }
}

static double da_chain_run(const da_ctx &cx,
                           dqrng::random_64bit_generator &rng_eng) {
  const int n = cx.n, nd = cx.nd, ng = cx.ng, S = cx.S, K = cx.K, p = cx.p;
  const int Sk = S + 1;
  const int warmup = cx.warmup, ndraws = cx.ndraws, thin = cx.thin;
  const int ncv = cx.ncv, spectra_maxitr = cx.spectra_maxitr;
  const bool tall = cx.tall, has_fixed = cx.has_fixed;
  const R_xlen_t nmiss = cx.nmiss, nstore = cx.nstore;
  const double dfP = cx.dfP, gram_max_cond = cx.gram_max_cond,
               spectra_tol = cx.spectra_tol;
  double *xp = cx.xp;
  double *dp = cx.dp;
  const double *gfp = cx.gfp;
  const int *midx = cx.midx, *mrow = cx.mrow, *mcol = cx.mcol, *sidx = cx.sidx,
            *scol = cx.scol;
  double sigma = cx.sigma;
  const int nopp = tall ? n : nd;

  const R_xlen_t nn = static_cast<R_xlen_t>(ng) * ng;
  std::vector<double> gram(nn);
  std::vector<double> evals(Sk), wmat(static_cast<R_xlen_t>(ng) * Sk);
  std::vector<double> proj(static_cast<R_xlen_t>(nopp) * S);
  std::vector<double> phi(S);
  std::vector<double> ut(static_cast<R_xlen_t>(n) * S);
  std::vector<double> vwt(static_cast<R_xlen_t>(nd) * S);
  std::vector<double> xtild(nmiss);
  std::vector<double> etd(cx.etd_in, cx.etd_in + nd);

  dqrng::normal_distribution z(0.0, 1.0);
  if (!cx.resuming) {
    for (R_xlen_t k = 0; k < nmiss; ++k)
      xtild[k] = cx.disp[mcol[k] - 1] * z(rng_eng);
  } else {
    std::copy(cx.xtild_in, cx.xtild_in + nmiss, xtild.begin());
  }

  std::vector<double> eig_start(ng);
  bool warm = cx.resuming;
  if (warm) std::copy(cx.eig_start_in, cx.eig_start_in + ng, eig_start.begin());

  const char uplo = 'L', trans_n = 'N', trans_t = 'T';
  const double one = 1.0, zero = 0.0;
  const int niter = warmup + ndraws * thin;

  const bool t_noise = R_FINITE(cx.nu);
  const double t_scale = t_noise ? std::sqrt((cx.nu - 2.0) / cx.nu) : 0.0;
  const double floor_scale =
      static_cast<double>(n) * static_cast<double>(p) / static_cast<double>(K);
  const double mt_nu = t_noise ? cx.nu : 3.0;
  const double mt_d = mt_nu / 2.0 - 1.0 / 3.0;
  const double mt_c = 1.0 / std::sqrt(9.0 * mt_d);

  const char *gram_trans = tall ? &trans_t : &trans_n;
  const char *proj_trans = tall ? &trans_n : &trans_t;
  const double *gram_beta = has_fixed ? &one : &zero;
  const double *moy_u = tall ? proj.data() : wmat.data();
  const double *moy_v = tall ? wmat.data() : proj.data();

  Eigen::Map<const Eigen::MatrixXd> gmap(gram.data(), ng, ng);
  Spectra::DenseSymMatProd<double> op(gmap);
#ifdef LOC_TIMER
  LOC_TIMER_OBJ(v2_da_phases);
  TimedSymMatProd top(op, v2_da_phases);
  Spectra::SymEigsSolver<TimedSymMatProd> eigs(top, Sk, ncv);
#else
  Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, Sk, ncv);
#endif

  const int refresh = cx.refresh, chain_id = cx.chain_id;
  const std::chrono::steady_clock::time_point t_start =
      std::chrono::steady_clock::now();
  int niter_w = 1;
  for (int t = niter; t >= 10; t /= 10) ++niter_w;

  for (int l = 1; l <= niter; ++l) {
    LOC_TIC(v2_da_phases, "iter");

    LOC_TIC(v2_da_phases, "istep");
    if (t_noise) {
      const double s = sigma * t_scale;
      for (R_xlen_t k = 0; k < nmiss; ++k) {
        const double zt = z(rng_eng);
        const double g2 = mt_chisq(rng_eng, z, mt_d, mt_c);
        xp[midx[k] - 1] = xtild[k] + s * zt * std::sqrt(mt_nu / g2);
      }
    } else {
      for (R_xlen_t k = 0; k < nmiss; ++k)
        xp[midx[k] - 1] = xtild[k] + sigma * z(rng_eng);
    }
    LOC_TOC(v2_da_phases, "istep");

    LOC_TIC(v2_da_phases, "store");
    if (l > warmup && (l - warmup) % thin == 0) {
      const R_xlen_t b = (l - warmup) / thin - 1;
      double *dcol = dp + b * nstore;
      for (R_xlen_t j = 0; j < nstore; ++j)
        dcol[j] = xp[midx[sidx[j] - 1] - 1] + etd[scol[j] - 1];
    }
    LOC_TOC(v2_da_phases, "store");

    LOC_TIC(v2_da_phases, "sweep");
    sweep_cols(xp, n, nd, etd.data());
    LOC_TOC(v2_da_phases, "sweep");

    LOC_TIC(v2_da_phases, "gram");

    if (has_fixed) {
      for (int c = 0; c < ng; ++c) {
        const R_xlen_t off = static_cast<R_xlen_t>(c) * ng;
        std::copy(gfp + off + c, gfp + off + ng, gram.data() + off + c);
      }
    }
    F77_CALL(dsyrk)
    (&uplo, gram_trans, &ng, &nopp, &one, xp, &n, gram_beta, gram.data(),
     &ng FCONE FCONE);

    double ss = 0.0;
    for (int i = 0; i < ng; ++i) ss += gram[static_cast<R_xlen_t>(i) * ng + i];
    LOC_TOC(v2_da_phases, "gram");

    LOC_TIC(v2_da_phases, "solve");
    const int nconv = spectra_solve(eigs, ng, Sk, spectra_tol, spectra_maxitr,
                                    evals.data(), wmat.data(), eig_start.data(),
                                    warm LOC_TIMER_ARG(v2_da_phases));
    LOC_TOC(v2_da_phases, "solve");
    if (nconv < Sk)
      Rcpp::stop("the PCA step did not converge at iteration %d. Raise "
                 "'maxitr' or 'ncv' in spectra_opts()",
                 l);
    warm = true;
    if (evals[S - 1] <= 0.0 || evals[0] / evals[S - 1] > gram_max_cond)
      Rcpp::stop("at iteration %d, component %d has almost no variance "
                 "(eigenvalue ratio %g). Lower 'rank'",
                 l, S, evals[0] / evals[S - 1]);

    LOC_TIC(v2_da_phases, "backout");

    double sd2 = 0.0;
    for (int s = 0; s < S; ++s) sd2 += evals[s];
    const double sigma2 = std::max(ss - sd2, 0.0) / dfP;
    sigma = std::sqrt(sigma2);

    const double floor_raw = floor_scale * sigma2;
    const double floor = std::min(floor_raw, evals[S]);
    double sphi = 0.0;
    for (int s = 0; s < S; ++s) {
      double ph = 1.0 - floor / evals[s];
      if (ph < 0.0) ph = 0.0;
      phi[s] = ph;
      sphi += ph;
    }

    F77_CALL(dgemm)
    (proj_trans, &trans_n, &nopp, &S, &ng, &one, xp, &n, wmat.data(), &ng,
     &zero, proj.data(), &nopp FCONE FCONE);
    LOC_TOC(v2_da_phases, "backout");

    LOC_TIC(v2_da_phases, "moy");
    moy_dot(moy_u, moy_v, phi.data(), mrow, mcol, n, nd, S, nmiss, ut.data(),
            vwt.data(), xtild.data());
    LOC_TOC(v2_da_phases, "moy");

    LOC_TIC(v2_da_phases, "pstep_noise");
    const double sd_t = std::sqrt(sigma2 * sphi / static_cast<double>(K));
    for (R_xlen_t k = 0; k < nmiss; ++k) xtild[k] += sd_t * z(rng_eng);
    LOC_TOC(v2_da_phases, "pstep_noise");

    LOC_TOC(v2_da_phases, "iter");

    if (refresh > 0 && (l == 1 || l % refresh == 0 || l == niter))
      Rcpp::message(Rcpp::wrap(
          tfm::format("Chain %d Iteration: %*d / %d [%3d%%] (%s)", chain_id,
                      niter_w, l, niter, static_cast<int>((100LL * l) / niter),
                      l <= warmup ? "Warmup" : "Sampling")));

    Rcpp::checkUserInterrupt();
  }

  if (refresh > 0) {
    const double secs = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - t_start)
                            .count();
    Rcpp::message(Rcpp::wrap(
        tfm::format("Chain %d finished in %.1f seconds.", chain_id, secs)));
  }

  std::copy(etd.begin(), etd.end(), cx.etd_out);
  std::copy(xtild.begin(), xtild.end(), cx.xtild_out);
  std::copy(eig_start.begin(), eig_start.end(), cx.eig_start_out);
  return sigma;
}
