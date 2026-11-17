// mipca.cpp -- the DA core: Spectra top-k, the streaming build, the chain loop
// and the multi-chain runner, in ONE translation unit.
//
// Never include RcppArmadillo.h in this translation unit: it clashes with
// R_ext/BLAS.h, which the Gram goes through. Spectra is the copy vendored in
// src/Spectra/ (upstream commit and file hashes in src/Spectra/MANIFEST), the
// 1.x API (selection rule an argument to compute(), op by reference,
// compute(selection, maxit, tol, sort)).
//
// [[Rcpp::depends(RcppEigen, dqrng, BH, sitmo)]]
#include <RcppEigen.h>
#include "Spectra/SymEigsSolver.h"
#include <R_ext/BLAS.h>
#include <dqrng.h>
#include <dqrng_distribution.h>
// The I-step's t draw: any boost::random distribution draws through dqrng's
// generator& seam (dqrng_distribution.h). BH is in LinkingTo for dqrng.
#include <boost/random/student_t_distribution.hpp>
#include <vector>
#include <stdexcept>
#include <cstdio>
#include <cmath>
#include <limits>
#include <sstream>
#include <iomanip>
#include <memory>
#include <chrono>
#include "loc_timer.h"

// ---- named constants --------------------------------------------------------
// Only numbers someone chose; the sampler's knobs cross the wire from R, which
// owns their defaults.

// The floor in ncv = min(n, max(2k + 1, 20)).
constexpr int SPECTRA_NCV_FLOOR = 20;

// ---- phase timing -----------------------------------------------------------
// LOC_TIC/LOC_TOC compile to nothing without -DLOC_TIMER (src/loc_timer.h).
// Phases do not partition an iteration: "iter" wraps the body, so its excess
// over the phase sum is time between phases.

// ---- top-k symmetric eigensolve ---------------------------------------------

// Picks ncv. The rule lives here and the default in R (ncv_req <= 0 means
// "apply the rule"); do not move either half.
static int spectra_ncv_choose(int n, int k, int ncv_req)
{
  if (ncv_req > 0)
  {
    if (ncv_req <= k || ncv_req > n)
      // k is the number of eigenpairs solved for, rank + 1 in the sampler,
      // and n the Gram side, the smaller dimension.
      Rcpp::stop("`ncv` in spectra_opts() must be more than %d and at most "
                 "%d, the smaller dimension of the matrix. It is %d.",
                 k, n, ncv_req);
    return ncv_req;
  }
  int ncv = 2 * k + 1;
  if (ncv < SPECTRA_NCV_FLOOR)
    ncv = SPECTRA_NCV_FLOOR;
  if (ncv > n)
    ncv = n;
  return ncv;
}

// Per-call preconditions, checked once at the loop's entry; not called from
// spectra_solve(), which runs every iteration.
static void spectra_check(double tol, int maxitr)
{
  if (!(tol > 0.0))
    Rcpp::stop("spectra tol must be positive; got %g.", tol);
  if (maxitr < 1)
    Rcpp::stop("spectra maxitr must be at least 1; got %d.", maxitr);
}

// Top-k eigenpairs, largest algebraic; returns the converged count. G is read
// only, so the loop hands over its live Gram uncopied. Either init fully
// resets the solver (src/Spectra/HermEigsBase.h), so a hoisted solver is
// bitwise a fresh one and consumes no R or dqrng stream. `warm` false starts
// from Spectra's SimpleRandom(0) vector; true starts from `start` (n
// doubles), which init copies before compute() runs. On convergence `start`
// is overwritten with the row sums of the k eigenvectors: the next solve's
// start, and the bytes the resume seam carries.
template <typename EigsT>
static int spectra_solve(EigsT &eigs, int n, int k, double tol, int maxitr,
                            double *values, double *vectors, double *start,
                            bool warm LOC_TIMER_PARAM(tm))
{
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
  if (nconv < k)
    return nconv;

  LOC_TIC(tm, "s_vectors");
  const Eigen::VectorXd ev = eigs.eigenvalues();
  const Eigen::MatrixXd evec = eigs.eigenvectors();
  LOC_TOC(tm, "s_vectors");
  LOC_TIC(tm, "s_copy");
  std::copy(ev.data(), ev.data() + k, values);
  std::copy(evec.data(), evec.data() + static_cast<std::ptrdiff_t>(n) * k,
            vectors);
  std::fill(start, start + n, 0.0);
  for (int s = 0; s < k; ++s)
  {
    const double *col = vectors + static_cast<std::ptrdiff_t>(s) * n;
    for (int i = 0; i < n; ++i)
      start[i] += col[i];
  }
  LOC_TOC(tm, "s_copy");
  return nconv;
}

#ifdef LOC_TIMER
// Timer builds only: forwards every matvec to the real operator and times it
// as "s_matvec". Two of them per solve run inside init() (Arnoldi::init),
// under "s_init"; the rest are compute()'s, so "s_compute" less "s_matvec"
// understates compute()'s non-matvec cost (Lanczos orthogonalization,
// restarts, the tridiagonal solves) by those two. Same arithmetic.
class TimedSymMatProd
{
public:
  using Scalar = double;
  TimedSymMatProd(const Spectra::DenseSymMatProd<double> &op, Rcpp::Timer &t)
      : op_(op), t_(t) {}
  Eigen::Index rows() const { return op_.rows(); }
  Eigen::Index cols() const { return op_.cols(); }
  void perform_op(const double *x_in, double *y_out) const
  {
    t_.tic("s_matvec");
    op_.perform_op(x_in, y_out);
    t_.toc("s_matvec");
  }

private:
  const Spectra::DenseSymMatProd<double> &op_;
  Rcpp::Timer &t_;
};
#endif

// ---- the sweep and moy kernels ----------------------------------------------
// Simple loops, vectorized or not at the compiler's call. Deliberate: the
// sweep multiplies by 1/n; the reductions feeding sigma^2 (trace, sum(evals),
// sum(phi)) stay long double, matching R's sum(), because sigma^2 is a
// difference of near-equal sums.

// MUTATES x in place (the dirty block); each column mean is added into etd.
static void sweep_cols(double *x, int n, int nd, double *etd)
{
  const double inv_n = 1.0 / static_cast<double>(n);
  for (int c = 0; c < nd; ++c)
  {
    double *col = x + static_cast<std::ptrdiff_t>(c) * n;
    double s = 0.0;
    for (int i = 0; i < n; ++i)
      s += col[i];
    const double mu = s * inv_n;
    for (int i = 0; i < n; ++i)
      col[i] -= mu;
    etd[c] += mu;
  }
}

// out[k] = sum_s u[row_k, s] * v[col_k, s] * w[s] over S-major copies, weights
// folded into v; one fused, serial pass over the holes.
static void moy_dot(const double *u, const double *v, const double *w,
                       const int *row, const int *col, int n, int nd, int S,
                       R_xlen_t m, double *ut, double *vwt, double *out)
{
  for (int s = 0; s < S; ++s)
  {
    const double *ucol = u + static_cast<std::ptrdiff_t>(s) * n;
    for (int i = 0; i < n; ++i)
      ut[static_cast<std::ptrdiff_t>(i) * S + s] = ucol[i];
    const double ws = w[s];
    const double *vcol = v + static_cast<std::ptrdiff_t>(s) * nd;
    for (int j = 0; j < nd; ++j)
      vwt[static_cast<std::ptrdiff_t>(j) * S + s] = vcol[j] * ws;
  }
  for (R_xlen_t k = 0; k < m; ++k)
  {
    const double *a = ut + static_cast<std::ptrdiff_t>(row[k] - 1) * S;
    const double *b = vwt + static_cast<std::ptrdiff_t>(col[k] - 1) * S;
    double acc = 0.0;
    for (int s = 0; s < S; ++s)
      acc += a[s] * b[s];
    out[k] = acc;
  }
}

// ---- the streaming build -----------------------------------------------------
// One column block per call, straight into the loop's inputs. A resume skips
// clean columns: nothing a resume reads comes from them. Holes are ISNAN, the
// same set as is.na(), mat_miss() and col_vars(); keep them in sync.

// MUTATES IN PLACE Xd, g_fixed (lower triangle), et, col_ss, mi_row, scratch.
// `src` is READ ONLY -- the caller's own matrix on the identity path -- and is
// a SEXP so an integer matrix is never coerced per block. BITWISE: et must
// equal .colMeans(na.rm = TRUE) and col_ss sum v * v, both long double in row
// order; build a reference before touching either.
// [[Rcpp::export(.build_block)]]
void build_block(SEXP src, int src_off, int j1, int nb, int transformed,
                    Rcpp::NumericVector scratch,
                    Rcpp::NumericMatrix Xd, Rcpp::NumericVector g_fixed,
                    Rcpp::NumericVector et, Rcpp::NumericVector col_ss,
                    Rcpp::IntegerVector mi_row, Rcpp::IntegerVector cnt,
                    Rcpp::IntegerVector col_ptr, Rcpp::IntegerVector dirty_pos,
                    int resume)
{
  const int n = Xd.nrow();
  const int nd = Xd.ncol();
  const int p = cnt.size();
  if (nb < 1 || j1 < 1 || j1 + nb - 1 > p)
    Rcpp::stop("block columns [%d, %d] fall outside 1..%d.", j1, j1 + nb - 1, p);
  if (TYPEOF(src) != REALSXP)
    Rcpp::stop("src must be a double matrix (storage.mode \"double\").");
  if (src_off < 0 || Rf_xlength(src) % n != 0 ||
      Rf_xlength(src) < static_cast<R_xlen_t>(n) * (src_off + nb))
    Rcpp::stop("src has %.0f cells; needs %d rows and at least %d columns.",
               static_cast<double>(Rf_xlength(src)), n, src_off + nb);
  if (et.size() != p || col_ss.size() != p || dirty_pos.size() != p ||
      col_ptr.size() != p + 1)
    Rcpp::stop("et, col_ss, dirty_pos must have length p and col_ptr p + 1.");
  const bool has_g = g_fixed.size() > 0;
  if (resume && has_g)
    Rcpp::stop("internal error: a resume build must not accumulate a Gram; "
               "g_fixed comes off the seam.");
  if (has_g && g_fixed.size() != static_cast<R_xlen_t>(n) * n)
    Rcpp::stop("g_fixed must be empty or n x n.");
  if (has_g && scratch.size() < static_cast<R_xlen_t>(n) * nb)
    Rcpp::stop("scratch has %.0f cells; needs at least %d x %d.",
               static_cast<double>(scratch.size()), n, nb);
  if (mi_row.size() != static_cast<R_xlen_t>(col_ptr[p]) - 1)
    Rcpp::stop("mi_row must have col_ptr[p + 1] - 1 entries.");

  const double *b = REAL(src) + static_cast<R_xlen_t>(src_off) * n;
  double *sc = has_g ? REAL(scratch) : nullptr;
  double *xd = REAL(Xd);
  int *rows_all = INTEGER(mi_row);
  int kc = 0; // clean columns compacted to scratch's front, for one dsyrk
  for (int jj = 0; jj < nb; ++jj)
  {
    const int j = j1 - 1 + jj; // 0-based global column
    if (resume && dirty_pos[j] == 0)
    {
      if (cnt[j] != 0)
        Rcpp::stop("internal error: column %d has holes and no slot in Xd.",
                   j + 1);
      continue;
    }
    const double *col = b + static_cast<R_xlen_t>(jj) * n;

    // Observed-cell mean, exactly as .colMeans(na.rm = TRUE) forms it.
    long double s = 0.0L;
    R_xlen_t nobs = 0;
    int holes = 0;
    for (int i = 0; i < n; ++i)
    {
      const double x = col[i];
      if (ISNAN(x))
        ++holes;
      else
      {
        if (!R_FINITE(x))
        {
          if (transformed)
            Rcpp::stop("The transform maps column %d, row %d to a non-finite "
                       "value.", j + 1, i + 1);
          Rcpp::stop("Z contains non-finite values.");
        }
        s += x;
        ++nobs;
      }
    }
    if (holes != cnt[j])
    {
      if (transformed && holes > cnt[j])
        Rcpp::stop("The transform maps %d observed value(s) in column %d to "
                   "NaN. A transform must give a finite value for every "
                   "observed cell.",
                   holes - cnt[j], j + 1);
      if (transformed)
        Rcpp::stop("The transform maps %d missing cell(s) in column %d to a "
                   "value. A transform must map NA to NA.",
                   cnt[j] - holes, j + 1);
      Rcpp::stop("internal error: column %d: the census counted %d holes, "
                 "the block holds %d.", j + 1, cnt[j], holes);
    }
    s /= nobs; // NaN when nobs == 0, as .colMeans; unreachable past admissibility
    const double mu = static_cast<double>(s);
    et[j] = mu;

    const int dp = dirty_pos[j];
    if (dp < 0 || dp > nd)
      Rcpp::stop("internal error: column %d: slot %d outside 1..%d.", j + 1, dp, nd);
    if (holes > 0 && dp == 0)
      Rcpp::stop("internal error: column %d has holes and no slot in Xd.", j + 1);
    if (dp == 0 && !has_g)
      Rcpp::stop("internal error: column %d has no slot in Xd and there is no "
                 "g_fixed to fold it into.", j + 1);
    double *dst = (dp > 0) ? xd + static_cast<R_xlen_t>(dp - 1) * n
                           : sc + static_cast<R_xlen_t>(kc) * n;
    int *rows = rows_all + (col_ptr[j] - 1);
    int k = 0;
    long double ss = 0.0L;
    for (int i = 0; i < n; ++i)
    {
      const double x = col[i];
      double v;
      if (ISNAN(x))
      {
        v = 0.0;
        rows[k++] = i + 1;
      }
      else
        v = x - mu;
      dst[i] = v;
      ss += static_cast<long double>(v) * v;
    }
    col_ss[j] = static_cast<double>(ss);
    if (dp == 0)
      ++kc;
  }
  if (kc > 0)
  {
    const char uplo = 'L', trans = 'N';
    const double one = 1.0;
    F77_CALL(dsyrk)
    (&uplo, &trans, &n, &kc, &one, sc, &n, &one, REAL(g_fixed), &n FCONE FCONE);
  }
}

// MUTATES g in place: mirrors a dsyrk('L') Gram's lower triangle onto the
// upper, as tcrossprod() would return it.
// [[Rcpp::export(.sym_from_lower)]]
void sym_from_lower(Rcpp::NumericMatrix g)
{
  const int n = g.nrow();
  if (g.ncol() != n)
    Rcpp::stop("g must be square.");
  double *x = REAL(g);
  for (int c = 0; c < n; ++c)
    for (int r = c + 1; r < n; ++r)
      x[static_cast<R_xlen_t>(r) * n + c] = x[static_cast<R_xlen_t>(c) * n + r];
}

// ---- the DA loop ------------------------------------------------------------
// da_chain_run() runs one chain on the calling thread; its only R API is the
// progress print and the interrupt check. Every normal draw is a scale times
// one unit normal, written x + s * z. A compiler may contract that into an
// FMA, so the draws are bitwise within one build and equal to rounding across
// builds.

// One chain's inputs as raw pointers. xp and dp are R-owned and MUTATED IN
// PLACE, so the caller must own them; the index vectors are read-only and
// shared across chains.
struct da_ctx
{
  int n, nd, ng;
  bool tall, has_fixed;
  R_xlen_t nmiss, nstore;
  double *xp;             // Xd, n x nd, mutated in place
  double *dp;             // draws, nstore x ndraws, written
  const double *gfp;      // g_fixed, NULL when !has_fixed
  const int *midx, *mrow, *mcol, *sidx, *scol;
  double sigma;           // initial sigma
  const double *etd_in;   // nd
  // The one bit the chain branches on for its start state; the three
  // pointers below follow it and are never tested in its place.
  bool resuming;
  const double *xtild_in; // nmiss; NULL on a fresh run
  // ng, the eigensolve's start vector; NULL on a fresh run, whose first
  // solve starts cold.
  const double *eig_start_in;
  // nd per-column scales already times `c`; NULL on a resume.
  const double *disp;
  double *etd_out;        // nd, written once at the end
  double *xtild_out;      // nmiss, written once at the end
  double *eig_start_out;  // ng, written once at the end
  int S, K;
  // FULL column count (clean + dirty), read only by the shrinkage floor; nd
  // is the dirty width. K = min(n - 1, p).
  int p;
  double dfP;
  double nu; // +Inf is Gaussian; finite > 2 is the variance-matched t
  int warmup, ndraws, thin;
  double gram_max_cond, spectra_tol;
  int ncv, spectra_maxitr;
  int chain_id; // 1-based, for progress lines only
  int refresh;  // iterations between progress lines; 0 = silent
};

static double da_chain_run(const da_ctx &cx,
                              dqrng::random_64bit_generator &rng_eng);

// Validation plus every chain-shared field. Per-chain pointers come back NULL,
// so a forgotten one faults loudly instead of aliasing chain 0.
static da_ctx da_ctx_build(int n, int nd,
                                 Rcpp::NumericVector g_fixed,
                                 Rcpp::IntegerVector miss_idx_l,
                                 Rcpp::IntegerVector mi_row,
                                 Rcpp::IntegerVector mi_col_l,
                                 Rcpp::IntegerVector store_idx,
                                 Rcpp::IntegerVector mi_col_store,
                                 int S, int K, int p, double dfP, double nu,
                                 int warmup, int ndraws, int thin,
                                 double gram_max_cond, double spectra_tol,
                                 int spectra_ncv, int spectra_maxitr,
                                 int refresh)
{
  const R_xlen_t nmiss = miss_idx_l.size();
  const R_xlen_t nstore = store_idx.size();
  const bool has_fixed = g_fixed.size() > 0;

  // The Gram goes on the smaller side; the cached clean block exists only on
  // the wide split, so has_fixed forces wide.
  const bool tall = !has_fixed && n > nd;
  const int ng = tall ? nd : n;

  // Mirrors .chk_rank() in R. The Gram-route bound is Spectra's own:
  // rank + 2 <= min(n, p). k_solve is da_chain_run()'s Sk.
  const int k_solve = S + 1;
  if (S < 1 || k_solve >= ng)
    Rcpp::stop("the core solves rank + 1 eigenpairs on the Gram, which needs "
               "rank + 2 <= %d (the smaller dimension); got rank %d.",
               ng, S);
  if (mi_row.size() != nmiss || mi_col_l.size() != nmiss)
    Rcpp::stop("mi_row, mi_col_l and miss_idx_l must agree.");
  if (mi_col_store.size() != nstore)
    Rcpp::stop("store_idx and mi_col_store must agree.");
  if (has_fixed &&
      g_fixed.size() != static_cast<R_xlen_t>(n) * static_cast<R_xlen_t>(n))
    Rcpp::stop("g_fixed must be n x n or empty.");
  if (dfP <= 0.0)
    Rcpp::stop("dfP must be positive.");
  if (!(nu > 2.0))
    Rcpp::stop("nu must be > 2 (Inf for Gaussian noise).");
  if (K < 1)
    Rcpp::stop("K must be positive.");
  if (p < nd)
    Rcpp::stop("p must be at least the dirty width nd.");
  // Mirrors mipca_core()'s floors, since this entry is directly callable.
  if (warmup < 0)
    Rcpp::stop("warmup must be >= 0.");
  if (ndraws < 1)
    Rcpp::stop("ndraws must be >= 1.");
  if (thin < 1)
    Rcpp::stop("thin must be >= 1.");
  // A division, so the check itself cannot overflow; needs thin >= 1 above.
  if (ndraws > (std::numeric_limits<int>::max() - warmup) / thin)
    // The loop length is an int.
    Rcpp::stop("`iter_warmup` + `iter_sampling` x `thin` must be at most "
               "%d. Lower one of them.",
               std::numeric_limits<int>::max());
  if (refresh < 0)
    Rcpp::stop("refresh must be >= 0 (0 = no progress lines).");
  spectra_check(spectra_tol, spectra_maxitr);

  // Every index checked here, so the hot loops can run on raw pointers.
  {
    const R_xlen_t ncell = static_cast<R_xlen_t>(n) * nd;
    for (R_xlen_t k = 0; k < nmiss; ++k)
    {
      if (miss_idx_l[k] < 1 || miss_idx_l[k] > ncell)
        Rcpp::stop("miss_idx_l out of range (expects 1-based).");
      if (mi_row[k] < 1 || mi_row[k] > n)
        Rcpp::stop("mi_row out of range (expects 1-based).");
      if (mi_col_l[k] < 1 || mi_col_l[k] > nd)
        Rcpp::stop("mi_col_l out of range (expects 1-based).");
    }
    for (R_xlen_t j = 0; j < nstore; ++j)
    {
      if (store_idx[j] < 1 || store_idx[j] > nmiss)
        Rcpp::stop("store_idx out of range (expects 1-based).");
      if (mi_col_store[j] < 1 || mi_col_store[j] > nd)
        Rcpp::stop("mi_col_store out of range (expects 1-based).");
    }
  }

  da_ctx cx;
  cx.n = n;
  cx.nd = nd;
  cx.ng = ng;
  cx.tall = tall;
  cx.has_fixed = has_fixed;
  cx.nmiss = nmiss;
  cx.nstore = nstore;
  cx.xp = NULL;
  cx.dp = NULL;
  cx.gfp = has_fixed ? REAL(g_fixed) : NULL;
  cx.midx = INTEGER(miss_idx_l);
  cx.mrow = INTEGER(mi_row);
  cx.mcol = INTEGER(mi_col_l);
  cx.sidx = INTEGER(store_idx);
  cx.scol = INTEGER(mi_col_store);
  cx.sigma = 0.0; // per-chain start state, set by the caller
  cx.etd_in = NULL;
  cx.resuming = false; // set by the entry below, with disp
  cx.xtild_in = NULL;
  cx.eig_start_in = NULL;
  cx.disp = NULL; // chain-shared, but set by the entry below once validated
  cx.etd_out = NULL;
  cx.xtild_out = NULL;
  cx.eig_start_out = NULL;
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

// A per-chain list element the ctx will point into must already be double.
static void need_double(SEXP x, const char *what, int chain)
{
  if (TYPEOF(x) != REALSXP)
    Rcpp::stop("chain %d: %s must be of type double.", chain, what);
}

// ---- the DA runner ----------------------------------------------------------
// The chains run in turn, each through da_chain_run() on its own stream.
// `resuming` is passed by the caller, never inferred from which inputs are
// empty. rng_state is raw operator<< text; R versions the seam.
// [[Rcpp::export(.mipca_da_loop_multi)]]
Rcpp::List mipca_da_loop_multi(Rcpp::List xd_list,
                                  Rcpp::NumericVector g_fixed,
                                  Rcpp::List etd_list,
                                  Rcpp::List xtild_list,
                                  Rcpp::List eig_start_list,
                                  Rcpp::NumericVector sigma_in,
                                  Rcpp::NumericVector disp,
                                  bool resuming,
                                  Rcpp::CharacterVector rng_state_in,
                                  Rcpp::IntegerVector miss_idx_l,
                                  Rcpp::IntegerVector mi_row,
                                  Rcpp::IntegerVector mi_col_l,
                                  Rcpp::IntegerVector store_idx,
                                  Rcpp::IntegerVector mi_col_store,
                                  Rcpp::List draws_list,
                                  Rcpp::IntegerVector streams,
                                  int S, int K, int p, double dfP, double nu,
                                  int warmup, int ndraws, int thin,
                                  double gram_max_cond, double spectra_tol,
                                  int spectra_ncv, int spectra_maxitr,
                                  int refresh)
{
  const int nc = static_cast<int>(xd_list.size());
  if (nc < 1)
    Rcpp::stop("need at least one chain.");
  if (draws_list.size() != nc || streams.size() != nc)
    Rcpp::stop("xd_list, draws_list and streams must have one entry per "
               "chain.");
  if (etd_list.size() != nc || xtild_list.size() != nc ||
      eig_start_list.size() != nc || sigma_in.size() != nc)
    Rcpp::stop("etd_list, xtild_list, eig_start_list and sigma_in must have "
               "one entry per chain.");
  if (resuming)
  {
    if (rng_state_in.size() != nc)
      Rcpp::stop("a resume must ship one serialized rng state per chain; "
                 "got %d for %d chain(s).",
                 static_cast<int>(rng_state_in.size()), nc);
  }
  else
  {
    if (rng_state_in.size() > 0)
      Rcpp::stop("a fresh run must ship no serialized rng state; got %d.",
                 static_cast<int>(rng_state_in.size()));
  }

  Rcpp::NumericMatrix xd0 = xd_list[0];
  const int n = xd0.nrow();
  const int nd = xd0.ncol();

  da_ctx shared =
      da_ctx_build(n, nd, g_fixed, miss_idx_l,
                      mi_row, mi_col_l, store_idx, mi_col_store, S, K, p,
                      dfP, nu, warmup, ndraws, thin, gram_max_cond, spectra_tol,
                      spectra_ncv, spectra_maxitr, refresh);

  if (resuming)
  {
    if (disp.size() > 0)
      Rcpp::stop("a dispersed init cannot be applied to a resume; the seam "
                 "already carries each chain's state.");
  }
  else
  {
    if (disp.size() != nd)
      Rcpp::stop("a fresh run needs one dispersion scale per dirty column "
                 "(%d); got %d.",
                 nd, static_cast<int>(disp.size()));
    for (int j = 0; j < nd; ++j)
      if (!R_finite(disp[j]) || disp[j] < 0.0)
        Rcpp::stop("disp must be finite and non-negative; entry %d is not.",
                   j + 1);
    shared.disp = REAL(disp);
  }
  shared.resuming = resuming;

  if (!resuming)
  {
    for (int a = 0; a < nc; ++a)
    {
      if (streams[a] < 0)
        Rcpp::stop("streams must be non-negative.");
      for (int b = a + 1; b < nc; ++b)
        if (streams[a] == streams[b])
          Rcpp::stop("streams must be unique; chains %d and %d both got "
                     "stream %d.",
                     a + 1, b + 1, streams[a]);
    }
  }

  std::vector<da_ctx> ctxs(nc, shared);
  Rcpp::List etd_out_l(nc), xtild_out_l(nc), eig_start_out_l(nc);
  for (int c = 0; c < nc; ++c)
  {
    // Doubles only, refused BEFORE the conversions below: an Rcpp vector
    // built from any other type owns a coerced copy that dies with this loop
    // body, and the ctx would keep a pointer into it.
    need_double(xd_list[c], "Xd", c + 1);
    need_double(draws_list[c], "draws", c + 1);
    need_double(etd_list[c], "etd", c + 1);
    need_double(xtild_list[c], "xtild_miss", c + 1);
    need_double(eig_start_list[c], "eig_start", c + 1);
    Rcpp::NumericMatrix xdc = xd_list[c];
    Rcpp::NumericMatrix drc = draws_list[c];
    Rcpp::NumericVector ec = etd_list[c];
    Rcpp::NumericVector xc = xtild_list[c];
    Rcpp::NumericVector sc = eig_start_list[c];
    if (xdc.nrow() != n || xdc.ncol() != nd)
      Rcpp::stop("chain %d: Xd must be %d x %d like chain 1's.", c + 1, n,
                 nd);
    if (drc.nrow() != shared.nstore || drc.ncol() != ndraws)
      Rcpp::stop("chain %d: draws must be nstore x ndraws.", c + 1);
    if (ec.size() != nd)
      Rcpp::stop("chain %d: etd must have one entry per column of Xd.",
                 c + 1);
    if (xc.size() != (resuming ? shared.nmiss : 0))
      Rcpp::stop(resuming
                     ? "chain %d: xtild_miss must have one entry per hole."
                     : "chain %d: xtild_miss must be empty on a fresh run; "
                       "the chain draws its own start from disp.",
                 c + 1);
    if (sc.size() != (resuming ? shared.ng : 0))
      Rcpp::stop(resuming
                     ? "chain %d: eig_start must have one entry per row of "
                       "the Gram."
                     : "chain %d: eig_start must be empty on a fresh run; "
                       "the chain's first eigensolve starts cold.",
                 c + 1);
    // Zero is left to Spectra, which refuses a start vector of norm zero.
    for (R_xlen_t i = 0; i < sc.size(); ++i)
      if (!R_finite(sc[i]))
        Rcpp::stop("chain %d: eig_start must be finite.", c + 1);
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

  // Every chain's generator is cloned here, before any chain runs: the
  // accessor reads the global XPtr. On a resume clone(0) only supplies the
  // type; operator>> overwrites its state.
  typedef dqrng::random_64bit_generator::result_type rng_result_t;
  std::vector<std::unique_ptr<dqrng::random_64bit_generator> > gens;
  gens.reserve(nc);
  {
    dqrng::random_64bit_accessor acc;
    for (int c = 0; c < nc; ++c)
    {
      if (!resuming)
      {
        gens.push_back(acc.clone(static_cast<rng_result_t>(streams[c])));
      }
      else
      {
        std::unique_ptr<dqrng::random_64bit_generator> g = acc.clone(0);
        std::istringstream is(Rcpp::as<std::string>(rng_state_in[c]));
        is >> *g;
        if (is.fail())
          Rcpp::stop("chain %d: serialized rng state failed to parse.",
                     c + 1);
        gens.push_back(std::move(g));
      }
    }
  }

  // An instrumented build stays single-chain: every chain would make its
  // own timer.
#ifdef LOC_TIMER
  if (nc > 1)
    Rcpp::stop("an instrumented (LOC_TIMER) build refuses more than one "
               "chain; profiling is single-chain only.");
#endif
  std::vector<double> sigma_out(nc);
  for (int c = 0; c < nc; ++c)
    sigma_out[c] = da_chain_run(ctxs[c], *gens[c]);

  Rcpp::CharacterVector rng_state(nc);
  for (int c = 0; c < nc; ++c)
  {
    std::ostringstream os;
    os << *gens[c];
    rng_state[c] = os.str();
  }

  return Rcpp::List::create(
      Rcpp::Named("sigma") = Rcpp::NumericVector(sigma_out.begin(),
                                                 sigma_out.end()),
      Rcpp::Named("etd") = etd_out_l,
      Rcpp::Named("xtild_miss") = xtild_out_l,
      Rcpp::Named("eig_start") = eig_start_out_l,
      Rcpp::Named("rng_state") = rng_state);
}

// One chain; returns the final sigma.
static double da_chain_run(const da_ctx &cx,
                              dqrng::random_64bit_generator &rng_eng)
{
  const int n = cx.n, nd = cx.nd, ng = cx.ng, S = cx.S, K = cx.K, p = cx.p;
  // S + 1 pairs. The extra pair's value is the floor's clamp; its vector,
  // wmat's last column, enters only the next solve's start (the row sums in
  // spectra_solve()). The back-out and moy read the first S columns.
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
  const int *midx = cx.midx, *mrow = cx.mrow, *mcol = cx.mcol,
            *sidx = cx.sidx, *scol = cx.scol;
  double sigma = cx.sigma;
  // The side of Xd the Gram is not on: ng is n when wide and nd when tall.
  const int nopp = tall ? n : nd;

  const R_xlen_t nn = static_cast<R_xlen_t>(ng) * ng;
  std::vector<double> gram(nn);
  std::vector<double> evals(Sk), wmat(static_cast<R_xlen_t>(ng) * Sk);
  // The factor the eigensolve did not produce: Xd' W when wide, Xd W when tall.
  std::vector<double> proj(static_cast<R_xlen_t>(nopp) * S);
  std::vector<double> phi(S);
  std::vector<double> ut(static_cast<R_xlen_t>(n) * S);
  std::vector<double> vwt(static_cast<R_xlen_t>(nd) * S);
  std::vector<double> xtild(nmiss);
  std::vector<double> etd(cx.etd_in, cx.etd_in + nd);

  // ---- start state: this chain's own dispersed draw, or the seam's ----------
  // Drawn here, on the chain's own clone: every chain starts itself, at a
  // dispersion multiple of each column's scale. A fresh chain consumes nmiss
  // normals and a resume none, so an extension's RNG stream continues
  // exactly where the run stopped. (The draws match an unbroken run to
  // rounding, not bitwise: the rebuilt panel drifts by rounding.)
  // The chain's one normal: every normal draw below is a scale times z.
  dqrng::normal_distribution z(0.0, 1.0);
  if (!cx.resuming)
  {
    for (R_xlen_t k = 0; k < nmiss; ++k)
      xtild[k] = cx.disp[mcol[k] - 1] * z(rng_eng);
  }
  else
  {
    std::copy(cx.xtild_in, cx.xtild_in + nmiss, xtild.begin());
  }

  // The eigensolve's start vector. A fresh chain's first solve is cold; every
  // later solve starts from the row sums of the eigenvectors before it, and a
  // resume from the seam's copy of them, so its first solve starts on the
  // bytes the unbroken run's would.
  std::vector<double> eig_start(ng);
  bool warm = cx.resuming;
  if (warm)
    std::copy(cx.eig_start_in, cx.eig_start_in + ng, eig_start.begin());

  const char uplo = 'L', trans_n = 'N', trans_t = 'T';
  const double one = 1.0, zero = 0.0;
  const int niter = warmup + ndraws * thin;

  // Loop-invariant factors, each the same expression the loop used to form.
  // t_scale is the variance match, sd = sigma * sqrt((nu - 2) / nu), read
  // only when the noise is t. floor_scale is the n p / K before sigma^2 in
  // the shrinkage floor.
  const bool t_noise = std::isfinite(cx.nu);
  const double t_scale = t_noise ? std::sqrt((cx.nu - 2.0) / cx.nu) : 0.0;
  const double floor_scale = static_cast<double>(n) * static_cast<double>(p) /
                             static_cast<double>(K);
  // Built once: the distribution holds its parameter and nothing else, so
  // one object draws what a fresh one per iteration would. Never drawn from
  // when the noise is Gaussian; the 1.0 only keeps the parameter legal.
  boost::random::student_t_distribution<double> tdist(t_noise ? cx.nu : 1.0);

  // The Gram and the back-out, by orientation. Wide: G = Xd Xd' ('N') and
  // proj = Xd' W ('T'). Tall: G = Xd' Xd ('T') and proj = Xd W ('N'). Either
  // way the rank-k update and the product run over nopp.
  const char *gram_trans = tall ? &trans_t : &trans_n;
  const char *proj_trans = tall ? &trans_n : &trans_t;
  const double *gram_beta = has_fixed ? &one : &zero;
  // moy's two factors: the n-side one first, the nd-side one second.
  const double *moy_u = tall ? proj.data() : wmat.data();
  const double *moy_v = tall ? wmat.data() : proj.data();

  // HOISTED solver: op holds a reference onto gram's buffer, so gram must
  // never be resized while it lives.
  Eigen::Map<const Eigen::MatrixXd> gmap(gram.data(), ng, ng);
  Spectra::DenseSymMatProd<double> op(gmap);
#ifdef LOC_TIMER
  LOC_TIMER_OBJ(v2_da_phases);
  TimedSymMatProd top(op, v2_da_phases);
  Spectra::SymEigsSolver<TimedSymMatProd> eigs(top, Sk, ncv);
#else
  Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double> > eigs(op, Sk, ncv);
#endif

  // Progress lines: each is built whole, written in one << and flushed.
  const int refresh = cx.refresh, chain_id = cx.chain_id;
  const std::chrono::steady_clock::time_point t_start =
      std::chrono::steady_clock::now();
  int niter_w = 1;
  for (int t = niter; t >= 10; t /= 10)
    ++niter_w;

  for (int l = 1; l <= niter; ++l)
  {
    LOC_TIC(v2_da_phases, "iter");

    // ---- I-step: impute every hole ------------------------------------------
    // Noise and scatter fused; the store below reads these doubles back from
    // xp. Finite nu: a variance-matched t at fixed nu, with the P-step
    // unchanged and no augmentation, so the chain targets no single model's
    // posterior; that is why a nu below 5 warns in R.
    LOC_TIC(v2_da_phases, "istep");
    if (t_noise)
    {
      const double s = sigma * t_scale;
      for (R_xlen_t k = 0; k < nmiss; ++k)
        xp[midx[k] - 1] = xtild[k] + s * tdist(rng_eng);
    }
    else
    {
      for (R_xlen_t k = 0; k < nmiss; ++k)
        xp[midx[k] - 1] = xtild[k] + sigma * z(rng_eng);
    }
    LOC_TOC(v2_da_phases, "istep");

    // ---- storage: the stored holes, de-centered with the CURRENT etd -------
    // Must run before the sweep below updates etd.
    LOC_TIC(v2_da_phases, "store");
    if (l > warmup && (l - warmup) % thin == 0)
    {
      const R_xlen_t b = (l - warmup) / thin - 1;
      double *dcol = dp + b * nstore;
      for (R_xlen_t j = 0; j < nstore; ++j)
        dcol[j] = xp[midx[sidx[j] - 1] - 1] + etd[scol[j] - 1];
    }
    LOC_TOC(v2_da_phases, "store");

    // ---- P-step -------------------------------------------------------------
    LOC_TIC(v2_da_phases, "sweep");
    sweep_cols(xp, n, nd, etd.data());
    LOC_TOC(v2_da_phases, "sweep");

    LOC_TIC(v2_da_phases, "gram");

    // G = g_fixed + Xd Xd', LOWER triangle only: must match Spectra's
    // DenseSymMatProd, which reads selfadjointView<Eigen::Lower>.
    if (has_fixed)
    {
      for (int c = 0; c < ng; ++c)
      {
        const R_xlen_t off = static_cast<R_xlen_t>(c) * ng;
        std::copy(gfp + off + c, gfp + off + ng, gram.data() + off + c);
      }
    }
    // lda is Xd's n in both orientations.
    F77_CALL(dsyrk)
    (&uplo, gram_trans, &ng, &nopp, &one, xp, &n, gram_beta, gram.data(),
     &ng FCONE FCONE);

    long double tr = 0.0L;
    for (int i = 0; i < ng; ++i)
      tr += gram[static_cast<R_xlen_t>(i) * ng + i];
    const double ss = static_cast<double>(tr);
    LOC_TOC(v2_da_phases, "gram");

    LOC_TIC(v2_da_phases, "solve");
    const int nconv = spectra_solve(eigs, ng, Sk, spectra_tol,
                                       spectra_maxitr, evals.data(),
                                       wmat.data(), eig_start.data(), warm
                                       LOC_TIMER_ARG(v2_da_phases));
    LOC_TOC(v2_da_phases, "solve");
    if (nconv < Sk)
    {
      char msg[160];
      // Spectra converged on nconv of the Sk eigenpairs of the Gram matrix.
      std::snprintf(msg, sizeof msg,
                    "The PCA step did not converge at iteration %d. Raise "
                    "`maxitr` or `ncv` in spectra_opts().",
                    l);
      throw std::runtime_error(msg);
    }
    warm = true; // converged, so eig_start now holds this solve's row sums
    if (evals[S - 1] <= 0.0 || evals[0] / evals[S - 1] > gram_max_cond)
    {
      char msg[160];
      // The Gram matrix is too ill-conditioned: the rank-th eigenvalue is
      // not positive, or too small beside the first.
      std::snprintf(msg, sizeof msg,
                    "At iteration %d, component %d has almost no variance "
                    "(eigenvalue ratio %g). Lower `rank`.",
                    l, S, evals[0] / evals[S - 1]);
      throw std::runtime_error(msg);
    }

    LOC_TIC(v2_da_phases, "backout");

    // sum(evals) is sum(d^2); d itself is never formed.
    long double sd2 = 0.0L;
    for (int s = 0; s < S; ++s)
      sd2 += static_cast<long double>(evals[s]);
    const double sigma2 = std::max(ss - static_cast<double>(sd2), 0.0) / dfP;
    sigma = std::sqrt(sigma2);

    // Shrinkage floor, Husson's form: missMDA::imputePCA, and Audigier, Husson
    // and Josse (2016, arXiv 1401.5747) eq. 4. The clamp makes phi >= 0, so
    // the guard cannot fire.
    const double floor_raw = floor_scale * sigma2;
    const double floor = std::min(floor_raw, evals[S]);
    long double sphi = 0.0L;
    for (int s = 0; s < S; ++s)
    {
      double ph = 1.0 - floor / evals[s];
      if (ph < 0.0)
        ph = 0.0; // clamped: no component sign flips
      phi[s] = ph;
      sphi += ph;
    }

    // UNSCALED on purpose: the textbook 1/d on this side cancels the d in its
    // d * phi weight, so phi alone weights moy = W diag(phi) (Xd' W)' (wide;
    // transposed when tall).
    F77_CALL(dgemm)
    (proj_trans, &trans_n, &nopp, &S, &ng, &one, xp, &n, wmat.data(), &ng,
     &zero, proj.data(), &nopp FCONE FCONE);
    LOC_TOC(v2_da_phases, "backout");

    // Overwrites xtild: legal only because the I-step above already read it.
    LOC_TIC(v2_da_phases, "moy");
    moy_dot(moy_u, moy_v, phi.data(), mrow, mcol, n, nd, S, nmiss, ut.data(),
            vwt.data(), xtild.data());
    LOC_TOC(v2_da_phases, "moy");

    LOC_TIC(v2_da_phases, "pstep_noise");
    const double sd_t =
        std::sqrt(sigma2 * static_cast<double>(sphi) / static_cast<double>(K));
    for (R_xlen_t k = 0; k < nmiss; ++k)
      xtild[k] += sd_t * z(rng_eng);
    LOC_TOC(v2_da_phases, "pstep_noise");

    LOC_TOC(v2_da_phases, "iter");

    if (refresh > 0 && (l == 1 || l % refresh == 0 || l == niter))
    {
      // Straight into Rcpp::Rcout, not snprintf into a fixed buffer: with a
      // run-time field width g++ 15 cannot bound the output and warns
      // (-Wformat-truncation), which R CMD check reports. Same bytes as
      // "Chain %d Iteration: %*d / %d [%3d%%] (%s)\n" while Rcout keeps
      // the classic locale, which nothing in this package changes.
      Rcpp::Rcout << "Chain " << chain_id << " Iteration: "
                  << std::setw(niter_w) << l << " / " << niter << " ["
                  << std::setw(3) << static_cast<int>((100LL * l) / niter)
                  << "%] (" << (l <= warmup ? "Warmup" : "Sampling")
                  << ")\n"
                  << std::flush;
    }

    Rcpp::checkUserInterrupt();
  }

  if (refresh > 0)
  {
    const double secs =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      t_start)
            .count();
    char line[96];
    std::snprintf(line, sizeof line, "Chain %d finished in %.1f seconds.\n",
                  chain_id, secs);
    Rcpp::Rcout << line << std::flush;
  }

  std::copy(etd.begin(), etd.end(), cx.etd_out);
  std::copy(xtild.begin(), xtild.end(), cx.xtild_out);
  std::copy(eig_start.begin(), eig_start.end(), cx.eig_start_out);
  return sigma;
}
