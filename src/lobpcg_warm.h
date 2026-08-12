#ifndef LOBPCG_WARM_H
#define LOBPCG_WARM_H
//
// Warm-start-only LOBPCG for the top-k eigenpairs of a real symmetric matrix.
//
// Logic:
//  - Caller provides a seeded LOBPCGState containing (X, P). X is n x k,
//  orthonormal, approximate top-k eigenvectors of the current A
//  (sign-canonicalized). P is n x k, intended to be a useful momentum direction;
//  we do NOT rely on X orth P or P^T P = I at entry, see "P invariant" below.
//  - lobpcg_solve computes AX = A*X and AP = A*P fresh at entry, since A
//  changes between outer calls. P is also re-projected onto X's orthogonal
//  complement and re-orthonormalized at entry, so the post-RR P from the
//  previous call (which is not orthonormal) is fine to commit.
//  - On converged == true, state is committed with the new (X, P). X is
//  sign-canonicalized; P is left as the post-RR search direction (not
//  orthonormal, not orth X), entry-time cleanup handles it.
//  - On converged == false, state is left untouched; caller must fall back
//  to the exact solver (dsyevr_) and re-seed.
//
// P invariant (committed state):
//  X is orthonormal and sign-canonical. P has the right shape (n x k);
//  it is NOT orthonormal nor strictly orth X. The entry-time projection
//  + QR re-orthonormalization rebuilds the clean basis we need for the
//  Rayleigh-Ritz.
//
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "pca_linalg_utils.h"

struct LOBPCGOptions
{
  double tol = 1e-8;
  int maxiter = 15;
  int stall_window = 5;
  double qr_rel_tol = 1e-12;
};

enum LOBPCGStatus
{
  LOBPCG_NOT_RUN = -1,
  LOBPCG_CONVERGED = 0,
  LOBPCG_INITIAL_RR_FAILED = 1,
  LOBPCG_STALLED = 2,
  LOBPCG_RESIDUAL_COLLAPSED = 3,
  LOBPCG_P_COLLAPSED = 4,
  LOBPCG_RR_FAILED = 5,
  LOBPCG_MAXITER = 6,
  LOBPCG_BAD_SEED = 7
};

struct LOBPCGResult
{
  arma::vec eigvals;
  int iterations = 0;
  bool converged = false;
  double max_rel_res = 0.0;
  int status = LOBPCG_NOT_RUN;
  int fail_iter = -1;
};

struct LOBPCGState
{
  arma::mat X;
  arma::mat P;
  std::vector<double> work;
  std::vector<arma::blas_int> iwork;

  bool seeded() const
  {
    return X.n_elem != 0 && P.n_elem != 0 &&
           X.n_rows == P.n_rows && X.n_cols == P.n_cols;
  }
};

namespace lobpcg_detail
{

  inline arma::uword qr_orthonormalize(arma::mat &V, double rel_tol,
                                       arma::mat *R_out = nullptr)
  {
    if (V.n_cols == 0)
    {
      return 0;
    }
    arma::mat Q, R;
    if (!arma::qr_econ(Q, R, V))
    {
      V.reset();
      return 0;
    }
    const arma::uword m = R.n_cols;
    arma::uword r = m;
    const double r00 = std::abs(R(0, 0));
    if (r00 == 0.0)
    {
      V.reset();
      return 0;
    }
    for (arma::uword i = 0; i < m; ++i)
    {
      if (std::abs(R(i, i)) < rel_tol * r00)
      {
        r = i;
        break;
      }
    }
    V = (r < m) ? Q.cols(0, r - 1) : Q;
    if (R_out)
    {
      *R_out = (r < m) ? R.submat(0, 0, r - 1, r - 1)
                       : R;
    }
    return r;
  }

  // QR-orthonormalize V, project it out of span(X), and QR again to tighten
  // orthogonality. required_cols > 0 additionally demands V keep that many
  // columns through both passes. On failure V may hold a partial result;
  // callers reset it as needed.
  inline bool qr_project_out_qr(arma::mat &V, const arma::mat &X,
                                double rel_tol, arma::uword required_cols)
  {
    if (qr_orthonormalize(V, rel_tol) == 0 ||
        (required_cols > 0 && V.n_cols != required_cols))
    {
      return false;
    }

    V -= X * (X.t() * V);

    if (qr_orthonormalize(V, rel_tol) == 0 ||
        (required_cols > 0 && V.n_cols != required_cols))
    {
      return false;
    }

    return true;
  }

  inline bool build_momentum_P(arma::mat &P,
                               const arma::mat &X_curr,
                               const arma::mat &X_prev,
                               const LOBPCGOptions &opt,
                               double min_projected_fro_norm = -1.0)
  {
    if (X_prev.n_rows != X_curr.n_rows ||
        X_prev.n_cols != X_curr.n_cols)
    {
      P.reset();
      return false;
    }

    P = X_prev - X_curr * (X_curr.t() * X_prev);

    if (min_projected_fro_norm >= 0.0 &&
        !(arma::norm(P, "fro") > min_projected_fro_norm))
    {
      P.reset();
      return false;
    }

    if (!qr_project_out_qr(P, X_curr, opt.qr_rel_tol, X_curr.n_cols))
    {
      P.reset();
      return false;
    }

    return true;
  }

  inline bool qr_orthonormalize_update_matvec(arma::mat &Q,
                                              arma::mat &AQ,
                                              const arma::mat &A,
                                              double rel_tol)
  {
    const arma::uword old_cols = Q.n_cols;

    arma::mat R;
    const arma::uword r = qr_orthonormalize(Q, rel_tol, &R);

    if (r == 0)
      return false;

    if (r == old_cols)
    {
      // Q_old = Q_new * R, so A*Q_new = A*Q_old * inv(R).
      AQ = arma::solve(arma::trimatl(R.t()), AQ.t()).t();
    }
    else
    {
      // QR dropped columns, so the old AQ cannot be transformed dimensionally.
      AQ = A * Q;
    }

    return true;
  }

  inline bool small_eig_desc(arma::mat &M, arma::vec &w, arma::mat &Z,
                             std::vector<double> &work_buf,
                             std::vector<arma::blas_int> &iwork_buf)
  {
    const arma::blas_int n = static_cast<arma::blas_int>(M.n_rows);
    if (n == 0)
    {
      w.reset();
      Z.reset();
      return true;
    }
    M = 0.5 * (M + M.t());
    w.set_size(n);
    arma::blas_int info = 0, lwork = -1, liwork = -1;
    const char jobz = 'V', uplo = 'U';
    double work_q = 0.0;
    arma::blas_int iwork_q = 0;

    arma::dsyevd_(&jobz, &uplo, &n, M.memptr(), &n, w.memptr(),
                  &work_q, &lwork, &iwork_q, &liwork, &info, 1, 1);
    if (info != 0)
    {
      return false;
    }

    if (work_buf.size() < static_cast<size_t>(work_q))
    {
      work_buf.resize(static_cast<size_t>(work_q));
    }
    if (iwork_buf.size() < static_cast<size_t>(iwork_q))
    {
      iwork_buf.resize(static_cast<size_t>(iwork_q));
    }
    lwork = static_cast<arma::blas_int>(work_buf.size());
    liwork = static_cast<arma::blas_int>(iwork_buf.size());

    arma::dsyevd_(&jobz, &uplo, &n, M.memptr(), &n, w.memptr(),
                  work_buf.data(), &lwork, iwork_buf.data(), &liwork, &info,
                  1, 1);
    if (info != 0)
    {
      return false;
    }

    Z = arma::fliplr(M);
    w = arma::reverse(w);
    return true;
  }

  inline void residual(const arma::mat &AX, const arma::mat &X,
                       const arma::vec &lam, arma::mat &R)
  {
    R.set_size(AX.n_rows, AX.n_cols);
    const arma::uword n = AX.n_rows, k = AX.n_cols;
    const double *ax = AX.memptr();
    const double *x = X.memptr();
    const double *l = lam.memptr();
    double *r = R.memptr();
    for (arma::uword j = 0; j < k; ++j)
    {
      const double lj = l[j];
      for (arma::uword i = 0; i < n; ++i)
      {
        r[i + j * n] = ax[i + j * n] - lj * x[i + j * n];
      }
    }
  }

  // Rayleigh-Ritz over the joint basis [X, B]. small_eig_desc symmetrizes
  // the assembled Gram, so the diagonal blocks need no pre-symmetrization.
  // On success X/AX are updated in place and (pp, app) hold the B-part of
  // the update (B*VB, AB*VB) for callers that commit it as the next P.
  inline bool rr_two_block(arma::mat &X, arma::mat &AX,
                           const arma::mat &B, const arma::mat &AB,
                           const arma::uword k, arma::vec &lam,
                           arma::mat &pp, arma::mat &app,
                           std::vector<double> &work_buf,
                           std::vector<arma::blas_int> &iwork_buf)
  {
    arma::mat XAX = X.t() * AX;
    arma::mat XAB = X.t() * AB;
    arma::mat BAB = B.t() * AB;

    const arma::uword sX = X.n_cols, sB = B.n_cols;
    const arma::uword sT = sX + sB;
    arma::mat gramA(sT, sT);
    gramA.submat(0, 0, sX - 1, sX - 1) = XAX;
    gramA.submat(0, sX, sX - 1, sT - 1) = XAB;
    gramA.submat(sX, 0, sT - 1, sX - 1) = XAB.t();
    gramA.submat(sX, sX, sT - 1, sT - 1) = BAB;

    arma::vec lam_s;
    arma::mat V_s;
    if (!small_eig_desc(gramA, lam_s, V_s, work_buf, iwork_buf))
    {
      return false;
    }
    lam = lam_s.head(k);
    const arma::mat Vtop = V_s.head_cols(k);
    const arma::mat VX = Vtop.rows(0, sX - 1);
    const arma::mat VB = Vtop.rows(sX, sT - 1);

    pp = B * VB;
    app = AB * VB;
    X = X * VX + pp;
    AX = AX * VX + app;
    return true;
  }

} // namespace lobpcg_detail

// ---------------------------------------------------------------------------
// seed LOBPCG state from two consecutive sign-canonicalized eigenvector
// blocks.
// ---------------------------------------------------------------------------
inline void seed_lobpcg_state(LOBPCGState &state,
                              const arma::mat &X_curr,
                              const arma::mat &X_prev,
                              const LOBPCGOptions &opt = LOBPCGOptions())
{
  state.X = X_curr;

  arma::mat P;

  const double no_min_temporal_norm = -1.0;

  if (!lobpcg_detail::build_momentum_P(
          P, X_curr, X_prev, opt, no_min_temporal_norm))
  {
    state.P.reset();
    return;
  }

  state.P = std::move(P);
}

// ---------------------------------------------------------------------------
// after a successful solve: rebuild P as the temporal momentum direction
// between two consecutive sign-canonicalized eigenvector blocks. When the
// temporal direction is degenerate (norm below the floor), the committed
// post-RR P is kept instead.
// ---------------------------------------------------------------------------
inline void refresh_momentum_P(LOBPCGState &state,
                               const arma::mat &X_curr,
                               const arma::mat &X_prev,
                               const LOBPCGOptions &opt = LOBPCGOptions())
{
  const double min_temporal_norm =
      1e-10 * std::sqrt(static_cast<double>(X_curr.n_cols));

  arma::mat P;
  if (lobpcg_detail::build_momentum_P(P, X_curr, X_prev, opt,
                                      min_temporal_norm))
  {
    state.P = std::move(P);
  }
}

// ---------------------------------------------------------------------------
// lobpcg_solve: warm start only, top-k, standard symmetric eigenproblem.
//
// requires a seeded state: an unseeded or mis-shaped (X, P) hard-errors via
// Rcpp::stop instead of falling back, so callers must check seeded() first.
// ---------------------------------------------------------------------------
inline LOBPCGResult lobpcg_solve(const arma::mat &A,
                                 LOBPCGState &state,
                                 const LOBPCGOptions &opt = LOBPCGOptions())
{
  using namespace lobpcg_detail;
  LOBPCGResult out;
  out.status = LOBPCG_MAXITER;

  if (!state.seeded())
    Rcpp::stop("lobpcg_solve: state not seeded with (X, P)");
  const arma::uword n = A.n_rows;
  const arma::uword k = state.X.n_cols;
  if (A.n_cols != n || state.X.n_rows != n)
    Rcpp::stop("lobpcg_solve: dimension mismatch");

  // work on local copies; state is committed only on success.
  arma::mat X = state.X;
  arma::mat P = state.P;

  // re-clean P at entry
  // the committed P from the previous solve is NOT orthonormal nor orth X
  // (it is the post-RR search direction pp = Ract*VR + Pact*VP). Project
  // out X, then orthonormalize.
  P -= X * (X.t() * P);
  bool have_P = qr_project_out_qr(P, X, opt.qr_rel_tol, k);
  if (!have_P)
  {
    P.reset();
  }
  // fresh matvecs: A has changed since the last call.
  arma::mat AX = A * X;
  arma::mat AP;
  if (have_P)
  {
    AP = A * P;
  }

  const double normA = std::max(1.0, arma::norm(A, "inf"));
  const double rel_tol_abs = opt.tol * normA;

  arma::vec lam;
  arma::mat rr_pp, rr_app; // B-part outputs of rr_two_block
  // initial RR over the JOINT basis [X, P] (not span(X) alone)
  if (have_P)
  {
    if (!rr_two_block(X, AX, P, AP, k, lam, rr_pp, rr_app,
                      state.work, state.iwork))
    {
      out.converged = false;
      out.status = LOBPCG_INITIAL_RR_FAILED;
      out.fail_iter = 0;
      return out;
    }
    // P, AP intentionally not rotated; loop's [X,R]-vs-P projection handles it.
  }
  else
  {
    // no usable P this call; fall back to X-only initial RR.
    arma::mat XAX = X.t() * AX;
    arma::mat V;
    if (!small_eig_desc(XAX, lam, V, state.work, state.iwork))
    {
      out.converged = false;
      out.status = LOBPCG_INITIAL_RR_FAILED;
      out.fail_iter = 0;
      return out;
    }
    X = X * V;
    AX = AX * V;
  }

  // periodically recompute AX/AP from scratch to shed accumulated
  // linear-update rounding error.
  constexpr int MATVEC_REFRESH = 20;

  arma::mat R;
  double smallest_res = std::numeric_limits<double>::infinity();
  int stall = 0;
  int it = 0;

  for (; it < opt.maxiter; ++it)
  {

    residual(AX, X, lam, R);
    const arma::vec rn = arma::sqrt(arma::sum(arma::square(R), 0)).t();
    const double max_res = rn.max();
    out.max_rel_res = max_res / normA;

    if (max_res <= rel_tol_abs)
    {
      out.converged = true;
      ++it;
      break;
    }

    // divergence guard: residual exploded well beyond what's possible
    // for a normalized eigenvector candidate. ||A*x - l*x|| / ||A|| <= 2
    // for any unit-norm x, so anything >> 2 means the iterate has lost
    // orthonormality and is amplifying. 10 leave headroom for
    // pathological mid-iteration excursions while still catching divergence.
    if (out.max_rel_res > 10.0)
    {
      out.status = LOBPCG_BAD_SEED;
      out.fail_iter = it;
      break;
    }

    if (max_res < smallest_res)
    {
      smallest_res = max_res;
      stall = 0;
    }
    else if (++stall >= opt.stall_window)
    {
      out.status = LOBPCG_STALLED;
      out.fail_iter = it;
      break;
    }

    const arma::uvec act_idx = arma::find(rn > rel_tol_abs);
    if (act_idx.is_empty())
    {
      out.converged = true;
      ++it;
      break;
    }
    arma::mat Ract = R.cols(act_idx);

    Ract -= X * (X.t() * Ract);
    if (!qr_project_out_qr(Ract, X, opt.qr_rel_tol, 0))
    {
      out.status = LOBPCG_RESIDUAL_COLLAPSED;
      out.fail_iter = it;
      break;
    }

    arma::mat ARact = A * Ract;
    // if we have no P this iteration, do RR on [X, Ract] only.
    if (!have_P)
    {
      if (!rr_two_block(X, AX, Ract, ARact, k, lam, rr_pp, rr_app,
                        state.work, state.iwork))
      {
        out.status = LOBPCG_RR_FAILED;
        out.fail_iter = it;
        break;
      }
      P = std::move(rr_pp);
      AP = std::move(rr_app);
      have_P = (P.n_cols > 0);
      if (((it + 1) % MATVEC_REFRESH) == 0)
      {
        AX = A * X;
        AP = A * P;
      }
      continue;
    }

    // standard 3-block path: project [X, Ract] out of P (linearity gives APact).
    // core of the algorithm
    arma::mat Pact, APact;
    {
      const arma::mat XtP = X.t() * P;
      const arma::mat RtP = Ract.t() * P;

      Pact = P - X * XtP - Ract * RtP;
      APact = AP - AX * XtP - ARact * RtP; // = A * Pact

      if (!qr_orthonormalize_update_matvec(Pact, APact, A, opt.qr_rel_tol))
      {
        out.status = LOBPCG_P_COLLAPSED;
        out.fail_iter = it;
        break;
      }

      const arma::mat XtP2 = X.t() * Pact;
      Pact -= X * XtP2;
      APact -= AX * XtP2;

      const arma::mat RtP2 = Ract.t() * Pact;
      Pact -= Ract * RtP2;
      APact -= ARact * RtP2;

      if (!qr_orthonormalize_update_matvec(Pact, APact, A, opt.qr_rel_tol))
      {
        out.status = LOBPCG_P_COLLAPSED;
        out.fail_iter = it;
        break;
      }
    }

    // diagonal blocks need no pre-symmetrization; small_eig_desc
    // symmetrizes the assembled Gram.
    arma::mat XAX_s = X.t() * AX;
    arma::mat XAR = X.t() * ARact;
    arma::mat XAP = X.t() * APact;
    arma::mat RAR = Ract.t() * ARact;
    arma::mat RAP = Ract.t() * APact;
    arma::mat PAP = Pact.t() * APact;

    const arma::uword sX = X.n_cols, sR = Ract.n_cols, sP = Pact.n_cols;
    const arma::uword sT = sX + sR + sP;

    arma::mat gramA(sT, sT);
    gramA.submat(0, 0, sX - 1, sX - 1) = XAX_s;
    gramA.submat(0, sX, sX - 1, sX + sR - 1) = XAR;
    gramA.submat(0, sX + sR, sX - 1, sT - 1) = XAP;
    gramA.submat(sX, 0, sX + sR - 1, sX - 1) = XAR.t();
    gramA.submat(sX, sX, sX + sR - 1, sX + sR - 1) = RAR;
    gramA.submat(sX, sX + sR, sX + sR - 1, sT - 1) = RAP;
    gramA.submat(sX + sR, 0, sT - 1, sX - 1) = XAP.t();
    gramA.submat(sX + sR, sX, sT - 1, sX + sR - 1) = RAP.t();
    gramA.submat(sX + sR, sX + sR, sT - 1, sT - 1) = PAP;

    arma::vec lam_s;
    arma::mat V_s;
    if (!small_eig_desc(gramA, lam_s, V_s, state.work, state.iwork))
    {
      out.status = LOBPCG_RR_FAILED;
      out.fail_iter = it;
      break;
    }

    lam = lam_s.head(k);
    const arma::mat Vtop = V_s.head_cols(k);
    const arma::mat VX = Vtop.rows(0, sX - 1);
    const arma::mat VR = Vtop.rows(sX, sX + sR - 1);
    const arma::mat VP = Vtop.rows(sX + sR, sT - 1);

    const arma::mat pp = Ract * VR + Pact * VP;
    const arma::mat app = ARact * VR + APact * VP;

    X = X * VX + pp;
    AX = AX * VX + app;
    P = pp;
    AP = app;
    if (((it + 1) % MATVEC_REFRESH) == 0)
    {
      AX = A * X;
      AP = A * P;
    }
  }

  out.iterations = it;
  out.eigvals = lam;

  if (out.converged)
  {
    out.status = LOBPCG_CONVERGED;
    pca_detail::canonicalize_signs(X);
    state.X = std::move(X);
    state.P = std::move(P);
  }

  return out;
}

#endif // LOBPCG_WARM_H
