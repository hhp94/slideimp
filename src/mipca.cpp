// [[Rcpp::depends(RcppEigen, RSpectra, dqrng, BH, sitmo, RcppThread)]]
#include <RcppEigen.h>
#include <SymEigs.h>
#include <R_ext/BLAS.h>
#include <dqrng.h>
#include <dqrng_distribution.h>

#include <boost/random/student_t_distribution.hpp>
#include <RcppThread.h>
#include <vector>
#include <stdexcept>
#include <cstdio>
#include <cmath>
#include <limits>
#include <sstream>
#include <memory>
#include <chrono>
#include "loc_timer.h"

constexpr int SPECTRA_NCV_FLOOR = 20;

static int spectra_ncv_choose(int n, int k, int ncv_req)
{
  if (ncv_req > 0)
  {

    if (ncv_req <= k || ncv_req > n)
      Rcpp::stop("ncv must satisfy k < ncv <= %d (the Gram side), where k = "
                 "%d is the number of eigenpairs solved for (rank + 1 in the "
                 "sampler); got %d.",
                 n, k, ncv_req);
    return ncv_req;
  }
  int ncv = 2 * k + 1;
  if (ncv < SPECTRA_NCV_FLOOR)
    ncv = SPECTRA_NCV_FLOOR;
  if (ncv > n)
    ncv = n;
  return ncv;
}

static void spectra_check(double tol, int maxitr)
{
  if (!(tol > 0.0))
    Rcpp::stop("spectra tol must be positive; got %g.", tol);
  if (maxitr < 1)
    Rcpp::stop("spectra maxitr must be at least 1; got %d.", maxitr);
}

template <typename EigsT>
static int spectra_solve(EigsT &eigs, int n, int k, double tol, int maxitr,
                            double *values, double *vectors)
{
  eigs.init();
  const int nconv = eigs.compute(maxitr, tol, Spectra::LARGEST_ALGE);
  if (nconv < k)
    return nconv;

  const Eigen::VectorXd ev = eigs.eigenvalues();
  const Eigen::MatrixXd evec = eigs.eigenvectors();
  std::copy(ev.data(), ev.data() + k, values);
  std::copy(evec.data(), evec.data() + static_cast<std::ptrdiff_t>(n) * k,
            vectors);
  return nconv;
}

static int spectra_topk_v2_raw(const double *g, int n, int k, int ncv,
                               double tol, int maxitr, double *values,
                               double *vectors)
{
  Eigen::Map<const Eigen::MatrixXd> gmap(g, n, n);
  Spectra::DenseSymMatProd<double> op(gmap);
  Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
                         Spectra::DenseSymMatProd<double> >
      eigs(&op, k, ncv);
  return spectra_solve(eigs, n, k, tol, maxitr, values, vectors);
}

// [[Rcpp::export(.spectra_topk)]]
Rcpp::List spectra_topk(Rcpp::NumericMatrix g, int k, double tol, int ncv_req,
                           int maxitr)
{
  const int n = g.nrow();
  if (g.ncol() != n)
    Rcpp::stop("G must be square.");
  if (k < 1 || k >= n)
    Rcpp::stop("k must lie in [1, n - 1]; got %d with n = %d.", k, n);
  spectra_check(tol, maxitr);

  const int ncv = spectra_ncv_choose(n, k, ncv_req);
  Rcpp::NumericVector values(k);
  Rcpp::NumericMatrix vectors(n, k);
  const int nconv =
      spectra_topk_v2_raw(g.begin(), n, k, ncv, tol, maxitr, values.begin(),
                          vectors.begin());
  if (nconv < k)
    Rcpp::stop("spectra_topk: only %d of %d eigenpairs converged.", nconv,
               k);

  return Rcpp::List::create(Rcpp::Named("values") = values,
                            Rcpp::Named("vectors") = vectors,
                            Rcpp::Named("ncv") = ncv);
}

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
  int kc = 0;
  for (int jj = 0; jj < nb; ++jj)
  {
    const int j = j1 - 1 + jj;
    if (resume && dirty_pos[j] == 0)
    {

      if (cnt[j] != 0)
        Rcpp::stop("internal error: column %d has holes and no slot in Xd.",
                   j + 1);
      continue;
    }
    const double *col = b + static_cast<R_xlen_t>(jj) * n;

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
            Rcpp::stop("the transform maps column %d, row %d to a non-finite "
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
        Rcpp::stop("the transform maps %d observed value(s) in column %d to "
                   "NaN; a transform must be finite on every observed cell.",
                   holes - cnt[j], j + 1);
      if (transformed)
        Rcpp::stop("the transform maps %d missing cell(s) in column %d to a "
                   "value; a transform must map NA to NA.",
                   cnt[j] - holes, j + 1);
      Rcpp::stop("internal error: column %d: the census counted %d holes, "
                 "the block holds %d.", j + 1, cnt[j], holes);
    }
    s /= nobs;
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

struct da_ctx
{
  int n, nd, ng;
  bool tall, has_fixed;
  R_xlen_t nmiss, nstore;
  double *xp;
  double *dp;
  const double *gfp;
  const int *midx, *mrow, *mcol, *sidx, *scol;
  double sigma;
  const double *etd_in;
  const double *xtild_in;

  const double *disp;
  double *etd_out;
  double *xtild_out;
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

  const bool tall = !has_fixed && n > nd;
  const int ng = tall ? nd : n;

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

  if (warmup < 0)
    Rcpp::stop("warmup must be >= 0.");
  if (ndraws < 1)
    Rcpp::stop("ndraws must be >= 1.");
  if (thin < 1)
    Rcpp::stop("thin must be >= 1.");

  if (ndraws > (std::numeric_limits<int>::max() - warmup) / thin)
    Rcpp::stop("warmup + ndraws * thin exceeds the largest int; "
               "the loop length is an int.");
  if (refresh < 0)
    Rcpp::stop("refresh must be >= 0 (0 = no progress lines).");
  spectra_check(spectra_tol, spectra_maxitr);

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
  cx.sigma = 0.0;
  cx.etd_in = NULL;
  cx.xtild_in = NULL;
  cx.disp = NULL;
  cx.etd_out = NULL;
  cx.xtild_out = NULL;
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

// [[Rcpp::export(.mipca_da_loop_multi)]]
Rcpp::List mipca_da_loop_multi(Rcpp::List xd_list,
                                  Rcpp::NumericVector g_fixed,
                                  Rcpp::List etd_list,
                                  Rcpp::List xtild_list,
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
                                  int nthreads,
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
      sigma_in.size() != nc)
    Rcpp::stop("etd_list, xtild_list and sigma_in must have one entry per "
               "chain.");

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
  if (nthreads < 1)
    Rcpp::stop("nthreads must be at least 1.");

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
  Rcpp::List etd_out_l(nc), xtild_out_l(nc);
  for (int c = 0; c < nc; ++c)
  {
    Rcpp::NumericMatrix xdc = xd_list[c];
    Rcpp::NumericMatrix drc = draws_list[c];
    Rcpp::NumericVector ec = etd_list[c];
    Rcpp::NumericVector xc = xtild_list[c];
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
    ctxs[c].chain_id = c + 1;
    ctxs[c].xp = REAL(xdc);
    ctxs[c].dp = REAL(drc);
    ctxs[c].etd_in = REAL(ec);
    ctxs[c].xtild_in = REAL(xc);
    ctxs[c].sigma = sigma_in[c];
    Rcpp::NumericVector eo(nd);
    Rcpp::NumericVector xo(static_cast<R_xlen_t>(shared.nmiss));
    etd_out_l[c] = eo;
    xtild_out_l[c] = xo;
    ctxs[c].etd_out = REAL(eo);
    ctxs[c].xtild_out = REAL(xo);
  }

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

  std::vector<double> sigma_out(nc);
  if (nc == 1)
  {

    sigma_out[0] = da_chain_run(ctxs[0], *gens[0]);
  }
  else
  {
#ifdef LOC_TIMER

    Rcpp::stop("the multi-chain path refuses an instrumented (LOC_TIMER) "
               "build; profiling is single-chain only.");
#endif
    const int nt = nthreads < nc ? nthreads : nc;
    RcppThread::ThreadPool pool(nt);
    for (int c = 0; c < nc; ++c)
      pool.push([&ctxs, &gens, &sigma_out, c]
                { sigma_out[c] = da_chain_run(ctxs[c], *gens[c]); });

    pool.wait();
  }

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
      Rcpp::Named("rng_state") = rng_state);
}

static double da_chain_run(const da_ctx &cx,
                              dqrng::random_64bit_generator &rng_eng)
{
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
  const int *midx = cx.midx, *mrow = cx.mrow, *mcol = cx.mcol,
            *sidx = cx.sidx, *scol = cx.scol;
  double sigma = cx.sigma;

  const R_xlen_t nn = static_cast<R_xlen_t>(ng) * ng;
  std::vector<double> gram(nn);
  std::vector<double> evals(Sk), wmat(static_cast<R_xlen_t>(ng) * Sk);

  std::vector<double> proj(static_cast<R_xlen_t>(tall ? n : nd) * S);
  std::vector<double> phi(S);
  std::vector<double> ut(static_cast<R_xlen_t>(n) * S);
  std::vector<double> vwt(static_cast<R_xlen_t>(nd) * S);
  std::vector<double> xtild(nmiss);
  std::vector<double> etd(nd);

  std::copy(cx.etd_in, cx.etd_in + nd, etd.begin());

  if (cx.disp != NULL)
  {
    dqrng::normal_distribution z(0.0, 1.0);
    for (R_xlen_t k = 0; k < nmiss; ++k)
      xtild[k] = cx.disp[mcol[k] - 1] * z(rng_eng);
  }
  else
  {
    std::copy(cx.xtild_in, cx.xtild_in + nmiss, xtild.begin());
  }

  const char uplo = 'L', trans_n = 'N', trans_t = 'T';
  const double one = 1.0, zero = 0.0;
  const int niter = warmup + ndraws * thin;

  Eigen::Map<const Eigen::MatrixXd> gmap(gram.data(), ng, ng);
  Spectra::DenseSymMatProd<double> op(gmap);
  Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
                         Spectra::DenseSymMatProd<double> >
      eigs(&op, Sk, ncv);

  const int refresh = cx.refresh, chain_id = cx.chain_id;
  const std::chrono::steady_clock::time_point t_start =
      std::chrono::steady_clock::now();

  int niter_w = 1;
  for (int t = niter; t >= 10; t /= 10)
    ++niter_w;

  LOC_TIMER_OBJ(v2_da_phases);

  for (int l = 1; l <= niter; ++l)
  {
    LOC_TIC(v2_da_phases, "iter");

    LOC_TIC(v2_da_phases, "istep");
    if (std::isfinite(cx.nu))
    {
      const double s = sigma * std::sqrt((cx.nu - 2.0) / cx.nu);
      boost::random::student_t_distribution<double> tdist(cx.nu);
      for (R_xlen_t k = 0; k < nmiss; ++k)
        xp[midx[k] - 1] = xtild[k] + s * tdist(rng_eng);
    }
    else
    {
      dqrng::normal_distribution dist(0.0, sigma);
      for (R_xlen_t k = 0; k < nmiss; ++k)
        xp[midx[k] - 1] = xtild[k] + dist(rng_eng);
    }
    LOC_TOC(v2_da_phases, "istep");

    LOC_TIC(v2_da_phases, "store");
    if (l > warmup && (l - warmup) % thin == 0)
    {
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

    if (has_fixed)
    {
      for (int c = 0; c < ng; ++c)
      {
        const R_xlen_t off = static_cast<R_xlen_t>(c) * ng;
        std::copy(gfp + off + c, gfp + off + ng, gram.data() + off + c);
      }
    }

    F77_CALL(dsyrk)
    (&uplo, tall ? &trans_t : &trans_n, &ng, tall ? &n : &nd, &one, xp, &n,
     has_fixed ? &one : &zero, gram.data(), &ng FCONE FCONE);

    long double tr = 0.0L;
    for (int i = 0; i < ng; ++i)
      tr += gram[static_cast<R_xlen_t>(i) * ng + i];
    const double ss = static_cast<double>(tr);
    LOC_TOC(v2_da_phases, "gram");

    LOC_TIC(v2_da_phases, "solve");
    const int nconv = spectra_solve(eigs, ng, Sk, spectra_tol,
                                       spectra_maxitr, evals.data(),
                                       wmat.data());
    LOC_TOC(v2_da_phases, "solve");
    if (nconv < Sk)
    {
      char msg[160];
      std::snprintf(msg, sizeof msg,
                    "Spectra failed to converge on the Gram matrix at "
                    "iteration %d (%d of %d pairs).",
                    l, nconv, Sk);
      throw std::runtime_error(msg);
    }
    if (evals[S - 1] <= 0.0 || evals[0] / evals[S - 1] > gram_max_cond)
    {
      char msg[160];
      std::snprintf(msg, sizeof msg,
                    "Gram matrix is too ill-conditioned at iteration %d "
                    "(eigenvalue ratio %g).",
                    l, evals[0] / evals[S - 1]);
      throw std::runtime_error(msg);
    }

    LOC_TIC(v2_da_phases, "backout");

    long double sd2 = 0.0L;
    for (int s = 0; s < S; ++s)
      sd2 += static_cast<long double>(evals[s]);
    const double sigma2 = std::max(ss - static_cast<double>(sd2), 0.0) / dfP;
    sigma = std::sqrt(sigma2);

    const double floor_raw =
        (static_cast<double>(n) * static_cast<double>(p) /
         static_cast<double>(K)) * sigma2;
    const double floor = std::min(floor_raw, evals[S]);
    long double sphi = 0.0L;
    for (int s = 0; s < S; ++s)
    {
      double ph = 1.0 - floor / evals[s];
      if (ph < 0.0)
        ph = 0.0;
      phi[s] = ph;
      sphi += ph;
    }

    if (tall)
    {
      F77_CALL(dgemm)
      (&trans_n, &trans_n, &n, &S, &nd, &one, xp, &n, wmat.data(), &nd, &zero,
       proj.data(), &n FCONE FCONE);
    }
    else
    {
      F77_CALL(dgemm)
      (&trans_t, &trans_n, &nd, &S, &n, &one, xp, &n, wmat.data(), &n, &zero,
       proj.data(), &nd FCONE FCONE);
    }
    LOC_TOC(v2_da_phases, "backout");

    LOC_TIC(v2_da_phases, "moy");
    moy_dot(tall ? proj.data() : wmat.data(),
               tall ? wmat.data() : proj.data(), phi.data(), mrow, mcol, n,
               nd, S, nmiss, ut.data(), vwt.data(), xtild.data());
    LOC_TOC(v2_da_phases, "moy");

    LOC_TIC(v2_da_phases, "pstep_noise");
    const double sd_t =
        std::sqrt(sigma2 * static_cast<double>(sphi) / static_cast<double>(K));
    {
      dqrng::normal_distribution dist(0.0, sd_t);
      for (R_xlen_t k = 0; k < nmiss; ++k)
        xtild[k] += dist(rng_eng);
    }
    LOC_TOC(v2_da_phases, "pstep_noise");

    LOC_TOC(v2_da_phases, "iter");

    if (refresh > 0 && (l == 1 || l % refresh == 0 || l == niter))
    {
      char line[96];
      std::snprintf(line, sizeof line,
                    "Chain %d Iteration: %*d / %d [%3d%%] (%s)\n", chain_id,
                    niter_w, l, niter,
                    static_cast<int>((100LL * l) / niter),
                    l <= warmup ? "Warmup" : "Sampling");
      RcppThread::Rcout << line;
    }

    RcppThread::checkUserInterrupt();
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
    RcppThread::Rcout << line;
  }

  std::copy(etd.begin(), etd.end(), cx.etd_out);
  std::copy(xtild.begin(), xtild.end(), cx.xtild_out);
  return sigma;
}
