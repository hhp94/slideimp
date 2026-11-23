#ifndef SLIDEIMP_LOBPCG_WARM_H
#define SLIDEIMP_LOBPCG_WARM_H
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
//  orthonormal, not orth X), entry-time cleanup handles it. P is empty in
//  the one case the P invariant below names.
//  - On converged == false, state is left untouched; caller must fall back
//  to the exact solver (dsyevr) and re-seed.
//
// P invariant (committed state):
//  X is orthonormal and sign-canonical. P is n x k; it is NOT orthonormal
//  nor strictly orth X. The entry-time projection + QR re-orthonormalization
//  rebuilds the clean basis we need for the Rayleigh-Ritz.
//  The exception: when P does not survive that cleanup, the first
//  Rayleigh-Ritz is over X alone and there is no P until the first
//  iteration builds one, so a call that converges at its first residual
//  test commits an empty P. The state is then not seeded() unless the
//  caller's refresh_momentum_P() builds a P from the temporal direction;
//  if it does not, the next call of hybrid_topk_eig_impl() takes the exact
//  path (HYB_REASON_NO_SEED under solver = lobpcg) and seeds the state
//  again from that solve, which can fail as after any exact solve.
//
// Storage: the iteration holds its matrices in lobpcg_detail::Block and
// enters BLAS and LAPACK only through src/si_lapack.h. The Gram comes in as
// a raw n x n buffer and is read in place, never copied; LOBPCGState holds
// its two blocks as Blocks, and LOBPCGResult's eigenvalues are a
// std::vector<double>.
//
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "pca_linalg_utils.h"
#include "si_lapack.h"

// tol and maxiter are R's defaults, those of lobpcg_control() in
// R/pca_imp.R; the kernel sets both from the call on every run.
struct LOBPCGOptions {
  double tol = 1e-9;
  int maxiter = 20;
  int stall_window = 5;
  double qr_rel_tol = 1e-12;
};

enum LOBPCGStatus {
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

namespace lobpcg_detail {

// A column-major matrix of doubles, nrow x ncol, unpadded, zero-filled when
// sized; a vector is one column. data() is the first cell, as for a
// std::vector, and the leading dimension is nrow.
struct Block {
  std::size_t nrow = 0;
  std::size_t ncol = 0;
  std::vector<double> cells;

  Block() = default;
  Block(const std::size_t rows, const std::size_t cols)
      : nrow(rows), ncol(cols), cells(rows * cols, 0.0) {}

  double *data() { return cells.data(); }
  const double *data() const { return cells.data(); }
  double &operator()(const std::size_t i, const std::size_t j) {
    return cells[i + j * nrow];
  }
  double operator()(const std::size_t i, const std::size_t j) const {
    return cells[i + j * nrow];
  }
  bool empty() const { return cells.empty(); }
  void clear() {
    nrow = 0;
    ncol = 0;
    cells.clear();
  }
};

// A copy of the rows x cols buffer at src, as a new block.
inline Block copy_of(const double *src, const std::size_t rows,
                     const std::size_t cols) {
  Block B(rows, cols);
  std::copy(src, src + rows * cols, B.data());
  return B;
}

} // namespace lobpcg_detail

struct LOBPCGResult {
  std::vector<double> eigvals;
  int iterations = 0;
  bool converged = false;
  double max_rel_res = 0.0;
  int status = LOBPCG_NOT_RUN;
  int fail_iter = -1;
};

struct LOBPCGState {
  lobpcg_detail::Block X;
  lobpcg_detail::Block P;

  bool seeded() const {
    return !X.empty() && !P.empty() && X.nrow == P.nrow && X.ncol == P.ncol;
  }
};

namespace lobpcg_detail {

// Whether a si_la status is a numerical failure of the factorization; false
// is success. A positive status is LAPACK's own info, or from
// si_la::solve_tri_right a zero on the diagonal, and err_nonfinite is the
// wrapper refusing a matrix that holds NaN or Inf: each call site
// treats both as the factorization failing. The rest are not the matrix's
// doing, and stop: err_alloc is memory, and err_arg and err_size cannot
// come from a call here. The three wrappers refuse with err_arg only a
// triangle other than 'U' or 'L' and a qr_econ of more columns than rows:
// each call passes 'U', and qr_orthonormalize() stops on such a block
// before qr_econ sees it. err_size is a dimension past int, which kernel
// setup refuses.
inline bool la_failed(const int status, const char *routine) {
  if (status == si_la::ok) {
    return false;
  }
  if (status > 0 || status == si_la::err_nonfinite) {
    return true;
  }
  if (status == si_la::err_alloc) {
    Rcpp::stop("lobpcg_solve: could not allocate the workspace of %s", routine);
  }
  Rcpp::stop("internal error: %s refused its arguments", routine);
}

// C = alpha * op(A) * B + beta * C, where op(A) is A' when trans_a and A
// is the a_rows x a_cols buffer at a. A block B of one column is a
// vector, and the product is dgemv; any other is dgemm. Operands that do
// not conform are an error here; no call passes them.
inline void product_into(const bool trans_a, const double alpha,
                         const double *a, const std::size_t a_rows,
                         const std::size_t a_cols, const Block &B,
                         const double beta, Block &C) {
  const std::size_t c_rows = trans_a ? a_cols : a_rows;
  const std::size_t inner = trans_a ? a_rows : a_cols;
  if (B.nrow != inner || C.nrow != c_rows || C.ncol != B.ncol) {
    Rcpp::stop("internal error: a product of blocks that do not conform");
  }
  const char trans = trans_a ? 'T' : 'N';
  const int st = (B.ncol == 1) ? si_la::gemv(trans, a_rows, a_cols, alpha, a,
                                             a_rows, B.data(), beta, C.data())
                               : si_la::gemm(trans, 'N', c_rows, B.ncol, inner,
                                             alpha, a, a_rows, B.data(), B.nrow,
                                             beta, C.data(), c_rows);
  if (st != si_la::ok) {
    Rcpp::stop("internal error: dgemm or dgemv refused its arguments");
  }
}

// A * B into a new block, A the a_rows x a_cols buffer at a: the Gram,
// read in place.
inline Block mul(const double *a, const std::size_t a_rows,
                 const std::size_t a_cols, const Block &B) {
  Block C(a_rows, B.ncol);
  product_into(false, 1.0, a, a_rows, a_cols, B, 0.0, C);
  return C;
}

inline Block mul(const Block &A, const Block &B) {
  return mul(A.data(), A.nrow, A.ncol, B);
}

// A' * B into a new block.
inline Block tmul(const Block &A, const Block &B) {
  Block C(A.ncol, B.ncol);
  product_into(true, 1.0, A.data(), A.nrow, A.ncol, B, 0.0, C);
  return C;
}

// C -= A * B as one BLAS call with alpha = -1 and beta = 1, for a C that is
// not inside the product. When it is, as in V -= X * (X' V), the product is
// formed first and subtracted cell by cell: project_out().
inline void sub_mul(Block &C, const Block &A, const Block &B) {
  product_into(false, -1.0, A.data(), A.nrow, A.ncol, B, 1.0, C);
}

// C -= D and C += D, cell by cell: how a product inside a sum or difference
// is added, once it is formed into a block of its own.
inline void sub_cells(Block &C, const Block &D) {
  if (C.nrow != D.nrow || C.ncol != D.ncol) {
    Rcpp::stop("internal error: LOBPCG blocks differ in size");
  }
  for (std::size_t i = 0; i < C.cells.size(); ++i) {
    C.cells[i] -= D.cells[i];
  }
}

inline void add_cells(Block &C, const Block &D) {
  if (C.nrow != D.nrow || C.ncol != D.ncol) {
    Rcpp::stop("internal error: LOBPCG blocks differ in size");
  }
  for (std::size_t i = 0; i < C.cells.size(); ++i) {
    C.cells[i] += D.cells[i];
  }
}

// A * B + C, the product formed first.
inline Block mul_plus(const Block &A, const Block &B, const Block &C) {
  Block D = mul(A, B);
  add_cells(D, C);
  return D;
}

// V -= X * (X' V): V with its part in span(X) taken out, X having
// orthonormal columns. V is inside the product, so the product is formed
// first, X' V innermost, and subtracted cell by cell.
inline void project_out(Block &V, const Block &X) {
  sub_cells(V, mul(X, tmul(X, V)));
}

// Rows r0 .. r0 + nr - 1 of the first nc columns of A, as a new block.
inline Block block_rows(const Block &A, const std::size_t r0,
                        const std::size_t nr, const std::size_t nc) {
  if (r0 > A.nrow || nr > A.nrow - r0 || nc > A.ncol) {
    Rcpp::stop("internal error: LOBPCG block out of range");
  }
  Block B(nr, nc);
  for (std::size_t j = 0; j < nc; ++j) {
    for (std::size_t i = 0; i < nr; ++i) {
      B(i, j) = A(r0 + i, j);
    }
  }
  return B;
}

// Write A, or with put_block_t its transpose, into C with its first cell
// at (r0, c0).
inline void put_block(Block &C, const std::size_t r0, const std::size_t c0,
                      const Block &A) {
  if (r0 > C.nrow || A.nrow > C.nrow - r0 || c0 > C.ncol ||
      A.ncol > C.ncol - c0) {
    Rcpp::stop("internal error: LOBPCG block out of range");
  }
  for (std::size_t j = 0; j < A.ncol; ++j) {
    for (std::size_t i = 0; i < A.nrow; ++i) {
      C(r0 + i, c0 + j) = A(i, j);
    }
  }
}

inline void put_block_t(Block &C, const std::size_t r0, const std::size_t c0,
                        const Block &A) {
  if (r0 > C.nrow || A.ncol > C.nrow - r0 || c0 > C.ncol ||
      A.nrow > C.ncol - c0) {
    Rcpp::stop("internal error: LOBPCG block out of range");
  }
  for (std::size_t j = 0; j < A.ncol; ++j) {
    for (std::size_t i = 0; i < A.nrow; ++i) {
      C(r0 + j, c0 + i) = A(i, j);
    }
  }
}

// The Frobenius norm as a plain sum of squares, unscaled where dnrm2 would
// scale as it sums: the value only feeds a threshold test, and the block is
// a projection of orthonormal columns, so the sum cannot overflow.
inline double frobenius_norm(const Block &A) {
  double ss = 0.0;
  for (const double v : A.cells) {
    ss += v * v;
  }
  return std::sqrt(ss);
}

inline bool all_finite(const std::vector<double> &x) {
  for (const double v : x) {
    if (!std::isfinite(v)) {
      return false;
    }
  }
  return true;
}

// The largest of n values, taken from -Inf with a > that NaN never passes,
// so NaN is skipped and n values of NaN give -Inf.
inline double max_skip_nan(const double *x, const std::size_t n) {
  double best = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < n; ++i) {
    if (x[i] > best) {
      best = x[i];
    }
  }
  return best;
}

// ||A||_inf of the a_rows x a_cols buffer at a: the largest row sum of
// absolute values, each row sum added from zero in column order, and the
// largest taken by max_skip_nan(). No BLAS.
inline double norm_inf(const double *a, const std::size_t a_rows,
                       const std::size_t a_cols) {
  std::vector<double> row_sum(a_rows, 0.0);
  for (std::size_t j = 0; j < a_cols; ++j) {
    const double *col = a + j * a_rows;
    for (std::size_t i = 0; i < a_rows; ++i) {
      row_sum[i] += std::abs(col[i]);
    }
  }
  return max_skip_nan(row_sum.data(), a_rows);
}

// The 2-norm of each column of R, as a plain sum of squares in row order,
// unscaled as for frobenius_norm(); the norms only feed the tests of the
// iteration and max_rel_res.
inline std::vector<double> column_norms(const Block &R) {
  std::vector<double> rn(R.ncol);
  for (std::size_t j = 0; j < R.ncol; ++j) {
    const double *col = R.data() + j * R.nrow;
    double ss = 0.0;
    for (std::size_t i = 0; i < R.nrow; ++i) {
      ss += col[i] * col[i];
    }
    rn[j] = std::sqrt(ss);
  }
  return rn;
}

// V orthonormalized in place by QR, keeping its columns up to the first
// whose pivot is below rel_tol times the first pivot; returns how many it
// kept. 0, with V cleared, when the QR fails or the first pivot is zero.
// When every column is kept, R_out receives the triangle R, V before the
// call being V after it times R; that is the one case its caller reads R
// in, and otherwise R_out is left as it was.
inline std::size_t qr_orthonormalize(Block &V, double rel_tol,
                                     Block *R_out = nullptr) {
  if (V.ncol == 0) {
    return 0;
  }
  // Every block here has at most k columns and n rows, k <= n. With more
  // columns than rows R would have fewer diagonal cells than the loop
  // below reads.
  if (V.ncol > V.nrow) {
    Rcpp::stop("internal error: a LOBPCG block has more columns than rows");
  }
  // si_la::qr_econ takes Q as v_rows x v_cols and R as v_cols x v_cols, and
  // writes every cell of both.
  const std::size_t v_rows = V.nrow;
  const std::size_t v_cols = V.ncol;
  Block Q(v_rows, v_cols);
  Block R(v_cols, v_cols);
  if (la_failed(si_la::qr_econ(v_rows, v_cols, V.data(), Q.data(), R.data()),
                "dgeqrf or dorgqr")) {
    V.clear();
    return 0;
  }
  std::size_t r = v_cols;
  const double r00 = std::abs(R(0, 0));
  if (r00 == 0.0) {
    V.clear();
    return 0;
  }
  for (std::size_t i = 0; i < v_cols; ++i) {
    if (std::abs(R(i, i)) < rel_tol * r00) {
      r = i;
      break;
    }
  }
  // keep the first r columns of Q, which lead its buffer.
  if (R_out && r == v_cols) {
    *R_out = std::move(R);
  }
  Q.cells.resize(v_rows * r);
  Q.ncol = r;
  V = std::move(Q);
  return r;
}

// QR-orthonormalize V, project it out of span(X), and QR again to tighten
// orthogonality. required_cols > 0 additionally demands V keep that many
// columns through both passes. On failure V may hold a partial result;
// callers reset it as needed.
inline bool qr_project_out_qr(Block &V, const Block &X, double rel_tol,
                              std::size_t required_cols) {
  if (qr_orthonormalize(V, rel_tol) == 0 ||
      (required_cols > 0 && V.ncol != required_cols)) {
    return false;
  }

  project_out(V, X);

  if (qr_orthonormalize(V, rel_tol) == 0 ||
      (required_cols > 0 && V.ncol != required_cols)) {
    return false;
  }

  return true;
}

inline bool build_momentum_P(Block &P, const Block &X_curr, const Block &X_prev,
                             const LOBPCGOptions &opt,
                             double min_projected_fro_norm) {
  if (X_prev.nrow != X_curr.nrow || X_prev.ncol != X_curr.ncol) {
    P.clear();
    return false;
  }

  // P = X_prev - X_curr * (X_curr' X_prev)
  P = X_prev;
  project_out(P, X_curr);

  if (min_projected_fro_norm >= 0.0 &&
      !(frobenius_norm(P) > min_projected_fro_norm)) {
    P.clear();
    return false;
  }

  if (!qr_project_out_qr(P, X_curr, opt.qr_rel_tol, X_curr.ncol)) {
    P.clear();
    return false;
  }

  return true;
}

// A is the n x n Gram, read in place.
inline bool qr_orthonormalize_update_matvec(Block &Q, Block &AQ,
                                            const double *A,
                                            const std::size_t n,
                                            double rel_tol) {
  const std::size_t old_cols = Q.ncol;

  Block R;
  const std::size_t r = qr_orthonormalize(Q, rel_tol, &R);

  if (r == 0) {
    return false;
  }

  if (r == old_cols) {
    // Q_old = Q_new * R, so A*Q_new = A*Q_old * inv(R): solve X * R = AQ
    // for X = A*Q_new, in place.
    if (AQ.ncol != old_cols) {
      Rcpp::stop("internal error: A*Q and Q differ in columns");
    }
    // A zero on the diagonal of R is the block collapsing, as a QR that
    // drops every column is.
    if (la_failed(
            si_la::solve_tri_right('U', AQ.nrow, AQ.ncol, R.data(), AQ.data()),
            "dtrsm")) {
      return false;
    }
  } else {
    // QR dropped columns, so the old AQ cannot be transformed dimensionally.
    AQ = mul(A, n, n, Q);
  }

  return true;
}

// All eigenpairs of the small symmetric M, descending, by dsyevd through
// si_la::syevd. M is symmetrized, and on success its cells become Z's and
// M comes back empty. false when dsyevd fails or M holds NaN or Inf, which
// the wrapper refuses before LAPACK; the caller then reports a failed
// Rayleigh-Ritz.
inline bool small_eig_desc(Block &M, Block &w, Block &Z) {
  const std::size_t n = M.nrow;
  if (n == 0) {
    w.clear();
    Z.clear();
    return true;
  }
  if (M.ncol != n) {
    Rcpp::stop("internal error: a Rayleigh-Ritz Gram is not square");
  }
  // M = 0.5 * (M + M'), every cell from M as it was.
  Block S(n, n);
  for (std::size_t j = 0; j < n; ++j) {
    for (std::size_t i = 0; i < n; ++i) {
      S(i, j) = (M(i, j) + M(j, i)) * 0.5;
    }
  }
  M = std::move(S);
  w = Block(n, 1);
  if (la_failed(si_la::syevd('U', n, M.data(), w.data()), "dsyevd")) {
    return false;
  }

  // dsyevd's order is ascending: reverse the columns and the values in
  // place, and hand the columns over as Z.
  for (std::size_t i = 0, j = n - 1; i < j; ++i, --j) {
    std::swap_ranges(M.data() + i * n, M.data() + (i + 1) * n,
                     M.data() + j * n);
  }
  std::reverse(w.cells.begin(), w.cells.end());
  Z = std::exchange(M, Block());
  return true;
}

inline void residual(const Block &AX, const Block &X, const Block &lam,
                     Block &R) {
  if (R.nrow != AX.nrow || R.ncol != AX.ncol) {
    R = Block(AX.nrow, AX.ncol);
  }
  const std::size_t n = AX.nrow, k = AX.ncol;
  const double *ax = AX.data();
  const double *x = X.data();
  const double *l = lam.data();
  double *r = R.data();
  for (std::size_t j = 0; j < k; ++j) {
    const double lj = l[j];
    for (std::size_t i = 0; i < n; ++i) {
      r[i + j * n] = ax[i + j * n] - lj * x[i + j * n];
    }
  }
}

// Rayleigh-Ritz over the joint basis [X, B], keeping the top sX Ritz pairs,
// as many as X has columns. small_eig_desc symmetrizes the assembled Gram,
// so the diagonal blocks need no pre-symmetrization. On success X/AX are
// updated in place and (pp, app) hold the B-part of the update (B*VB,
// AB*VB) for callers that commit it as the next P.
inline bool rr_two_block(Block &X, Block &AX, const Block &B, const Block &AB,
                         Block &lam, Block &pp, Block &app) {
  const Block XAX = tmul(X, AX);
  const Block XAB = tmul(X, AB);
  const Block BAB = tmul(B, AB);

  const std::size_t sX = X.ncol, sB = B.ncol;
  const std::size_t sT = sX + sB;
  Block gramA(sT, sT);
  put_block(gramA, 0, 0, XAX);
  put_block(gramA, 0, sX, XAB);
  put_block_t(gramA, sX, 0, XAB);
  put_block(gramA, sX, sX, BAB);

  Block lam_s;
  Block V_s;
  if (!small_eig_desc(gramA, lam_s, V_s)) {
    return false;
  }
  lam = block_rows(lam_s, 0, sX, 1);
  const Block VX = block_rows(V_s, 0, sX, sX);
  const Block VB = block_rows(V_s, sX, sB, sX);

  pp = mul(B, VB);
  app = mul(AB, VB);
  X = mul_plus(X, VX, pp);
  AX = mul_plus(AX, VX, app);
  return true;
}

} // namespace lobpcg_detail

// ---------------------------------------------------------------------------
// seed LOBPCG state from two consecutive sign-canonicalized eigenvector
// blocks.
// ---------------------------------------------------------------------------
inline void seed_lobpcg_state(LOBPCGState &state,
                              const lobpcg_detail::Block &X_curr,
                              const lobpcg_detail::Block &X_prev,
                              const LOBPCGOptions &opt) {
  state.X = X_curr;

  const double no_min_temporal_norm = -1.0;

  // build_momentum_P() leaves state.P empty when it fails.
  lobpcg_detail::build_momentum_P(state.P, X_curr, X_prev, opt,
                                  no_min_temporal_norm);
}

// ---------------------------------------------------------------------------
// after a successful solve: rebuild P as the temporal momentum direction
// between two consecutive sign-canonicalized eigenvector blocks. When the
// temporal direction is degenerate (norm below the floor), the committed
// post-RR P is kept instead.
// ---------------------------------------------------------------------------
inline void refresh_momentum_P(LOBPCGState &state,
                               const lobpcg_detail::Block &X_curr,
                               const lobpcg_detail::Block &X_prev,
                               const LOBPCGOptions &opt) {
  // the temporal direction is degenerate below a Frobenius norm of
  // TEMPORAL_NORM_FLOOR times the square root of its column count.
  constexpr double TEMPORAL_NORM_FLOOR = 1e-10;
  const double min_temporal_norm =
      TEMPORAL_NORM_FLOOR * std::sqrt(static_cast<double>(X_curr.ncol));

  lobpcg_detail::Block P;
  if (lobpcg_detail::build_momentum_P(P, X_curr, X_prev, opt,
                                      min_temporal_norm)) {
    state.P = std::move(P);
  }
}

// ---------------------------------------------------------------------------
// lobpcg_solve: warm start only, top-k, standard symmetric eigenproblem.
//
// requires a seeded state: an unseeded or mis-shaped (X, P) hard-errors via
// Rcpp::stop instead of falling back, so callers must check seeded() first.
// a is the n x n Gram, column-major, both triangles valid; it is read in
// place, never copied or written.
// ---------------------------------------------------------------------------
inline LOBPCGResult lobpcg_solve(const double *a, const std::size_t n,
                                 LOBPCGState &state, const LOBPCGOptions &opt) {
  using namespace lobpcg_detail;
  LOBPCGResult out;
  out.status = LOBPCG_MAXITER;

  if (!state.seeded()) {
    Rcpp::stop("internal error: lobpcg_solve: state not seeded with (X, P)");
  }
  const std::size_t k = state.X.ncol;
  if (state.X.nrow != n) {
    Rcpp::stop("internal error: lobpcg_solve: dimension mismatch");
  }

  // work on local copies; state is committed only on success.
  Block X = state.X;
  Block P = state.P;

  // re-clean P at entry
  // the committed P from the previous solve is NOT orthonormal nor orth X
  // (it is the post-RR search direction pp = Ract*VR + Pact*VP). Project
  // out X, then orthonormalize.
  project_out(P, X);
  bool have_P = qr_project_out_qr(P, X, opt.qr_rel_tol, k);
  if (!have_P) {
    P.clear();
  }
  // fresh matvecs: A has changed since the last call.
  Block AX = mul(a, n, n, X);
  Block AP;
  if (have_P) {
    AP = mul(a, n, n, P);
  }

  // ||A||_inf scales both the convergence test below and the divergence
  // guard's max_rel_res, so it has to be A's actual norm. Substitute 1.0 only
  // for the degenerate Gram whose norm is zero or non-finite, where the ratio
  // is undefined; a zero Gram has an exactly zero residual, so any positive
  // tolerance converges on the first pass either way. Flooring at 1.0
  // unconditionally would make the test ABSOLUTE for every Gram with a norm
  // under 1 - the ordinary case under scale = FALSE on bounded data, where
  // A's diagonal is a weighted column variance - and would leave the
  // divergence guard unreachable, since max_rel_res would be under-reported
  // by the same factor.
  const double normA_inf = norm_inf(a, n, n);
  const double normA =
      (normA_inf > 0.0 && std::isfinite(normA_inf)) ? normA_inf : 1.0;
  const double rel_tol_abs = opt.tol * normA;

  Block lam;
  Block rr_pp, rr_app; // B-part outputs of rr_two_block
  // initial RR over the JOINT basis [X, P] (not span(X) alone)
  if (have_P) {
    if (!rr_two_block(X, AX, P, AP, lam, rr_pp, rr_app)) {
      out.status = LOBPCG_INITIAL_RR_FAILED;
      out.fail_iter = 0;
      return out;
    }
    // P, AP intentionally not rotated; loop's [X,R]-vs-P projection handles it.
  } else {
    // no usable P this call; fall back to X-only initial RR.
    Block XAX = tmul(X, AX);
    Block V;
    if (!small_eig_desc(XAX, lam, V)) {
      out.status = LOBPCG_INITIAL_RR_FAILED;
      out.fail_iter = 0;
      return out;
    }
    X = mul(X, V);
    AX = mul(AX, V);
  }

  // periodically recompute AX/AP from scratch to shed accumulated
  // linear-update rounding error.
  constexpr int MATVEC_REFRESH = 20;

  // the divergence guard's bound on max_rel_res, in the loop below.
  constexpr double DIVERGED_REL_RES = 10.0;

  Block R;
  double smallest_res = std::numeric_limits<double>::infinity();
  int stall = 0;
  int it = 0;

  for (; it < opt.maxiter; ++it) {

    residual(AX, X, lam, R);
    const std::vector<double> rn = column_norms(R);
    const double max_res = max_skip_nan(rn.data(), rn.size());
    out.max_rel_res = max_res / normA;

    // A residual norm that is NaN or Inf is the iterate diverging, whatever
    // the largest finite norm says: max_skip_nan() passes over NaN, so
    // without this a block of NaN would pass the convergence test below.
    if (!all_finite(rn)) {
      out.status = LOBPCG_BAD_SEED;
      out.fail_iter = it;
      break;
    }

    if (max_res <= rel_tol_abs) {
      out.converged = true;
      ++it;
      break;
    }

    // divergence guard: residual exploded well beyond what's possible
    // for a normalized eigenvector candidate. ||A*x - l*x|| / ||A|| <= 2
    // for any unit-norm x, so anything >> 2 means the iterate has lost
    // orthonormality and is amplifying. 10 leaves headroom for
    // pathological mid-iteration excursions while still catching divergence.
    if (out.max_rel_res > DIVERGED_REL_RES) {
      out.status = LOBPCG_BAD_SEED;
      out.fail_iter = it;
      break;
    }

    if (max_res < smallest_res) {
      smallest_res = max_res;
      stall = 0;
    } else if (++stall >= opt.stall_window) {
      out.status = LOBPCG_STALLED;
      out.fail_iter = it;
      break;
    }

    // the active residuals: the columns whose norm is above the tolerance,
    // in order.
    std::vector<std::size_t> act_idx;
    for (std::size_t j = 0; j < rn.size(); ++j) {
      if (rn[j] > rel_tol_abs) {
        act_idx.push_back(j);
      }
    }
    if (act_idx.empty()) {
      out.converged = true;
      ++it;
      break;
    }
    Block Ract(R.nrow, act_idx.size());
    for (std::size_t c = 0; c < act_idx.size(); ++c) {
      std::copy(R.data() + act_idx[c] * R.nrow,
                R.data() + (act_idx[c] + 1) * R.nrow, Ract.data() + c * R.nrow);
    }

    project_out(Ract, X);
    if (!qr_project_out_qr(Ract, X, opt.qr_rel_tol, 0)) {
      out.status = LOBPCG_RESIDUAL_COLLAPSED;
      out.fail_iter = it;
      break;
    }

    Block ARact = mul(a, n, n, Ract);
    // if we have no P this iteration, do RR on [X, Ract] only.
    if (!have_P) {
      if (!rr_two_block(X, AX, Ract, ARact, lam, rr_pp, rr_app)) {
        out.status = LOBPCG_RR_FAILED;
        out.fail_iter = it;
        break;
      }
      // rr_pp has as many columns as X, k >= 1.
      P = std::move(rr_pp);
      AP = std::move(rr_app);
      have_P = true;
      if (((it + 1) % MATVEC_REFRESH) == 0) {
        AX = mul(a, n, n, X);
        AP = mul(a, n, n, P);
      }
      continue;
    }

    // standard 3-block path: project [X, Ract] out of P (linearity gives APact).
    Block Pact, APact;
    {
      const Block XtP = tmul(X, P);
      const Block RtP = tmul(Ract, P);

      // Pact = P - X * XtP - Ract * RtP, and APact = A * Pact by linearity,
      // each product formed first. Pact and APact take over the cells of P
      // and AP, which are not read again: below they become pp and app, and
      // every exit before that leaves without committing P.
      Pact = std::exchange(P, Block());
      sub_cells(Pact, mul(X, XtP));
      sub_cells(Pact, mul(Ract, RtP));
      APact = std::exchange(AP, Block());
      sub_cells(APact, mul(AX, XtP));
      sub_cells(APact, mul(ARact, RtP));

      if (!qr_orthonormalize_update_matvec(Pact, APact, a, n, opt.qr_rel_tol)) {
        out.status = LOBPCG_P_COLLAPSED;
        out.fail_iter = it;
        break;
      }

      const Block XtP2 = tmul(X, Pact);
      sub_mul(Pact, X, XtP2);
      sub_mul(APact, AX, XtP2);

      const Block RtP2 = tmul(Ract, Pact);
      sub_mul(Pact, Ract, RtP2);
      sub_mul(APact, ARact, RtP2);

      if (!qr_orthonormalize_update_matvec(Pact, APact, a, n, opt.qr_rel_tol)) {
        out.status = LOBPCG_P_COLLAPSED;
        out.fail_iter = it;
        break;
      }
    }

    // diagonal blocks need no pre-symmetrization; small_eig_desc
    // symmetrizes the assembled Gram.
    const Block XAX = tmul(X, AX);
    const Block XAR = tmul(X, ARact);
    const Block XAP = tmul(X, APact);
    const Block RAR = tmul(Ract, ARact);
    const Block RAP = tmul(Ract, APact);
    const Block PAP = tmul(Pact, APact);

    const std::size_t sX = X.ncol, sR = Ract.ncol, sP = Pact.ncol;
    const std::size_t sT = sX + sR + sP;

    Block gramA(sT, sT);
    put_block(gramA, 0, 0, XAX);
    put_block(gramA, 0, sX, XAR);
    put_block(gramA, 0, sX + sR, XAP);
    put_block_t(gramA, sX, 0, XAR);
    put_block(gramA, sX, sX, RAR);
    put_block(gramA, sX, sX + sR, RAP);
    put_block_t(gramA, sX + sR, 0, XAP);
    put_block_t(gramA, sX + sR, sX, RAP);
    put_block(gramA, sX + sR, sX + sR, PAP);

    Block lam_s;
    Block V_s;
    if (!small_eig_desc(gramA, lam_s, V_s)) {
      out.status = LOBPCG_RR_FAILED;
      out.fail_iter = it;
      break;
    }

    lam = block_rows(lam_s, 0, k, 1);
    const Block VX = block_rows(V_s, 0, sX, k);
    const Block VR = block_rows(V_s, sX, sR, k);
    const Block VP = block_rows(V_s, sX + sR, sP, k);

    // pp = Ract * VR + Pact * VP, and app likewise, each product formed
    // first.
    Block pp = mul(Ract, VR);
    add_cells(pp, mul(Pact, VP));
    Block app = mul(ARact, VR);
    add_cells(app, mul(APact, VP));

    X = mul_plus(X, VX, pp);
    AX = mul_plus(AX, VX, app);
    P = std::move(pp);
    AP = std::move(app);
    if (((it + 1) % MATVEC_REFRESH) == 0) {
      AX = mul(a, n, n, X);
      AP = mul(a, n, n, P);
    }
  }

  out.iterations = it;
  out.eigvals = std::move(lam.cells);

  if (out.converged) {
    out.status = LOBPCG_CONVERGED;
    pca_detail::canonicalize_signs(X.data(), X.nrow, X.ncol);
    state.X = std::move(X);
    state.P = std::move(P);
  }

  return out;
}

#endif // SLIDEIMP_LOBPCG_WARM_H
