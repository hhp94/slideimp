// mipca_diag.cpp -- the per-cell kernel: classic split-Rhat with the Gelman
// fold, and one Geyer ESS.
//
// MUST STAY a transcription of posterior 1.7.0's rhat_basic / ess_basic,
// diffable line for line; equal to rounding, not bitwise (posterior's
// autocovariance is an FFT, this one a direct sum).
// Output must stay BITWISE identical across tile sizes: no cross-tile
// reduction, no state carried between cells.

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

// ---- named constants --------------------------------------------------------
// The bare literals in ess_one() and diag_na() are posterior's own (.ess(),
// is_constant()), left as written for the line-for-line diff.

// Lags of autocovariance per refill; sets how far ahead lags are computed,
// never the order a lag is summed in.
constexpr int ACOV_LAG_BLOCK = 32;

// Workspace, allocated once per call and reused across cells; it must carry
// no state from one cell to the next.
struct Ws {
  std::vector<double> fold;  // S x C, |cell - median(cell)|
  std::vector<double> xs;    // ni x nc, split chains of cell
  std::vector<double> xf;    // ni x nc, split chains of fold
  std::vector<double> med;   // S x C, median scratch (nth_element mutates)
  std::vector<double> cm;    // nc, per-column means
  std::vector<double> cen;   // ni x nc, centered columns
  std::vector<double> am;    // ni, autocovariance averaged over columns
  std::vector<double> rho;   // ni, Geyer rho ladder
};

// The NA policy: any non-finite value, or no variation at all.
bool diag_na(const double* x, int n) {
  for (int i = 0; i < n; ++i) {
    if (!R_finite(x[i])) return true;
  }
  double lo = x[0], hi = x[0];
  for (int i = 1; i < n; ++i) {
    lo = std::min(lo, x[i]);
    hi = std::max(hi, x[i]);
  }
  return (hi - lo) < DBL_EPSILON;
}

// Halve every chain; an odd S drops the middle draw. Column order must match
// R's cbind(first halves, second halves).
void split_chains(const double* x, int S, int C, double* out) {
  if (S == 1) {
    std::copy(x, x + C, out);
    return;
  }
  const int ni = S / 2;
  const int start2 = (S + 1) / 2;  // == ceiling(S/2); odd S skips the middle
  for (int c = 0; c < C; ++c) {
    const double* col = x + static_cast<size_t>(c) * S;
    std::copy(col, col + ni, out + static_cast<size_t>(c) * ni);
    std::copy(col + start2, col + S, out + static_cast<size_t>(C + c) * ni);
  }
}

// As stats::median(). Finite input only: nth_element over a NaN is undefined.
double median_all(const double* x, int n, std::vector<double>& scratch) {
  std::copy(x, x + n, scratch.begin());
  const int h = n / 2;
  std::nth_element(scratch.begin(), scratch.begin() + h, scratch.begin() + n);
  const double m1 = scratch[h];
  if (n % 2 == 1) return m1;
  const double m2 = *std::max_element(scratch.begin(), scratch.begin() + h);
  return (m1 + m2) / 2.0;
}

// Split-Rhat for one already-split ni x nc matrix.
double rhat_one(const double* x, int ni, int nc, Ws& w) {
  if (nc < 2 || diag_na(x, ni * nc)) return NA_REAL;
  double wsum = 0.0;
  for (int j = 0; j < nc; ++j) {
    const double* col = x + static_cast<size_t>(j) * ni;
    double m = 0.0;
    for (int t = 0; t < ni; ++t) m += col[t];
    m /= ni;
    w.cm[j] = m;
    double acc = 0.0;
    for (int t = 0; t < ni; ++t) acc += (col[t] - m) * (col[t] - m);
    wsum += acc / (ni - 1.0);  // ni == 1 gives 0/0 = NaN
  }
  const double W = wsum / nc;
  double mm = 0.0;
  for (int j = 0; j < nc; ++j) mm += w.cm[j];
  mm /= nc;
  double B = 0.0;
  for (int j = 0; j < nc; ++j) B += (w.cm[j] - mm) * (w.cm[j] - mm);
  B /= (nc - 1.0);
  return std::sqrt((ni * B / W + ni - 1.0) / ni);
}

// Geyer initial-positive plus initial-monotone ESS for one already-split
// ni x nc matrix. Reads stay below ni through the `t < ni - 5` bound.
double ess_one(const double* x, int ni, int nc, Ws& w) {
  if (ni < 3 || diag_na(x, ni * nc)) return NA_REAL;
  for (int j = 0; j < nc; ++j) {
    const double* col = x + static_cast<size_t>(j) * ni;
    double m = 0.0;
    for (int t = 0; t < ni; ++t) m += col[t];
    m /= ni;
    w.cm[j] = m;
    double* c = w.cen.data() + static_cast<size_t>(j) * ni;
    for (int t = 0; t < ni; ++t) c[t] = col[t] - m;
  }
  // am[k] = mean over columns of the biased autocovariance at lag k.
  double* am = w.am.data();
  int have = 0;
  auto ensure = [&](int k) {
    while (have <= k) {
      const int hi = std::min(have + ACOV_LAG_BLOCK, ni);
      for (int lag = have; lag < hi; ++lag) am[lag] = 0.0;
      for (int j = 0; j < nc; ++j) {
        const double* c = w.cen.data() + static_cast<size_t>(j) * ni;
        for (int lag = have; lag < hi; ++lag) {
          double acc = 0.0;
          for (int u = 0; u + lag < ni; ++u) acc += c[u] * c[u + lag];
          am[lag] += acc;
        }
      }
      for (int lag = have; lag < hi; ++lag) am[lag] /= nc * (double)ni;
      have = hi;
    }
  };
  ensure(1);
  const double mean_var = am[0] * ni / (ni - 1.0);
  double var_plus = mean_var * (ni - 1.0) / ni;
  if (nc > 1) {
    double mm = 0.0;
    for (int j = 0; j < nc; ++j) mm += w.cm[j];
    mm /= nc;
    double acc = 0.0;
    for (int j = 0; j < nc; ++j) acc += (w.cm[j] - mm) * (w.cm[j] - mm);
    var_plus += acc / (nc - 1.0);
  }
  std::fill(w.rho.begin(), w.rho.begin() + ni, 0.0);
  double* rho = w.rho.data();
  int t = 0;
  double rho_even = 1.0;
  rho[0] = 1.0;
  double rho_odd = 1.0 - (mean_var - am[1]) / var_plus;
  rho[1] = rho_odd;
  while (t < ni - 5 && !std::isnan(rho_even + rho_odd) &&
         rho_even + rho_odd > 0) {
    t += 2;
    ensure(t + 1);
    rho_even = 1.0 - (mean_var - am[t]) / var_plus;
    rho_odd = 1.0 - (mean_var - am[t + 1]) / var_plus;
    if (rho_even + rho_odd >= 0) {
      rho[t] = rho_even;
      rho[t + 1] = rho_odd;
    }
  }
  const int max_t = t;
  if (rho_even > 0) rho[max_t] = rho_even;
  t = 0;
  while (t <= max_t - 4) {
    t += 2;
    if (rho[t] + rho[t + 1] > rho[t - 2] + rho[t - 1]) {
      rho[t] = (rho[t - 2] + rho[t - 1]) / 2.0;
      rho[t + 1] = rho[t];
    }
  }
  const double ess = nc * (double)ni;
  // posterior's `rho[1:max_t]` quirk at max_t == 0, kept on purpose; R
  // reports ESS as NA where the loop above cannot run (.cell_diag_read()).
  double s = 0.0;
  if (max_t == 0) {
    s = rho[0];
  } else {
    for (int i = 0; i < max_t; ++i) s += rho[i];
  }
  const double tau = -1.0 + 2.0 * s + rho[max_t];
  return ess / std::max(tau, 1.0 / std::log10(ess));
}

// Sized once per call; every tile reuses it.
void ws_alloc(Ws& w, int S, int C, int ni, int nc) {
  w.fold.resize(static_cast<size_t>(S) * C);
  w.xs.resize(static_cast<size_t>(ni) * nc);
  w.xf.resize(static_cast<size_t>(ni) * nc);
  w.med.resize(static_cast<size_t>(S) * C);
  w.cm.resize(nc);
  w.cen.resize(static_cast<size_t>(ni) * nc);
  w.am.resize(ni);
  w.rho.resize(ni);
}

// One tile: gather, then per-cell compute, writing only its own slice of outp
// (rhat at [i], ess at [ncell + i]).
void diag_tile(int i0, int nt, const int* rows0,
               const std::vector<const double*>& base, int nmiss, int S,
               int C, int ni, int nc, int ncell, double* outp, Ws& w,
               double* tb, int* trow) {
  const int n_full = S * C;
  for (int k = 0; k < nt; ++k) trow[k] = rows0[i0 + k];

  // The gather: per (chain, draw) slab, the tile's rows in one ascending burst.
  for (int c = 0; c < C; ++c) {
    const size_t coff = static_cast<size_t>(c) * S;
    for (int t = 0; t < S; ++t) {
      const double* p = base[c] + static_cast<size_t>(t) * nmiss;
      for (int k = 0; k < nt; ++k) {
        tb[static_cast<size_t>(k) * n_full + coff + t] = p[trow[k]];
      }
    }
  }

  for (int k = 0; k < nt; ++k) {
    const double* cell = tb + static_cast<size_t>(k) * n_full;
    split_chains(cell, S, C, w.xs.data());

    bool bad = false;
    for (int q = 0; q < n_full; ++q) {
      if (!R_finite(cell[q])) {
        bad = true;
        break;
      }
    }
    double rhat;
    if (bad) {
      // As in R: a non-finite cell poisons the fold arm, so rhat is NA. ESS
      // judges the SPLIT matrix alone, so an odd S's middle draw can differ.
      rhat = NA_REAL;
    } else {
      const double med = median_all(cell, n_full, w.med);
      for (int q = 0; q < n_full; ++q) w.fold[q] = std::fabs(cell[q] - med);
      split_chains(w.fold.data(), S, C, w.xf.data());
      const double a = rhat_one(w.xs.data(), ni, nc, w);
      const double b = rhat_one(w.xf.data(), ni, nc, w);
      rhat = (ISNAN(a) || ISNAN(b)) ? NA_REAL : std::max(a, b);
    }
    outp[i0 + k] = rhat;
    outp[static_cast<size_t>(ncell) + i0 + k] = ess_one(w.xs.data(), ni, nc, w);
  }
}

}  // namespace

// Draws and row indices only, never the map or a holdout: R picks the cells
// (.diag_sample(), .tune_diag_draw()), so this cannot know a held-out cell.
// [[Rcpp::export]]
Rcpp::NumericMatrix mipca_cell_diag(Rcpp::List draws,
                                       Rcpp::IntegerVector rows,
                                       int tile_cells) {
  const int C = draws.size();
  if (C < 1) Rcpp::stop("draws must hold at least one chain.");
  if (tile_cells < 1) Rcpp::stop("tile_cells must be a positive count.");
  std::vector<Rcpp::NumericMatrix> ms;
  ms.reserve(C);
  for (int c = 0; c < C; ++c) {
    ms.emplace_back(Rcpp::as<Rcpp::NumericMatrix>(draws[c]));
  }
  const int nmiss = ms[0].nrow();
  const int S = ms[0].ncol();
  for (int c = 1; c < C; ++c) {
    if (ms[c].nrow() != nmiss || ms[c].ncol() != S) {
      Rcpp::stop("chains disagree on dimensions.");
    }
  }
  // rows0 is the 0-based copy the tiles read.
  const int ncell = rows.size();
  std::vector<int> rows0(ncell);
  for (int i = 0; i < ncell; ++i) {
    if (rows[i] == NA_INTEGER || rows[i] < 1 || rows[i] > nmiss) {
      Rcpp::stop("rows must be 1-based indices into the hole map.");
    }
    if (i > 0 && rows[i] <= rows[i - 1]) {
      Rcpp::stop("rows must be strictly ascending.");
    }
    rows0[i] = rows[i] - 1;
  }

  const int ni = (S == 1) ? 1 : S / 2;
  const int nc = (S == 1) ? C : 2 * C;
  const int n_full = S * C;

  std::vector<const double*> base(C);
  for (int c = 0; c < C; ++c) base[c] = ms[c].begin();

  const int tile = std::min(tile_cells, std::max(ncell, 1));
  const int ntile = (ncell + tile - 1) / tile;  // 0 when ncell == 0

  Rcpp::NumericMatrix out(ncell, 2);
  double* outp = out.begin();

  // The tiles run in turn on one workspace, an interrupt check before each.
  Ws w;
  ws_alloc(w, S, C, ni, nc);
  std::vector<double> tb(static_cast<size_t>(tile) * n_full);
  std::vector<int> trow(tile);
  for (int ti = 0; ti < ntile; ++ti) {
    Rcpp::checkUserInterrupt();
    const int i0 = ti * tile;
    diag_tile(i0, std::min(tile, ncell - i0), rows0.data(), base, nmiss, S, C,
              ni, nc, ncell, outp, w, tb.data(), trow.data());
  }

  out.attr("dimnames") = Rcpp::List::create(
      R_NilValue, Rcpp::CharacterVector::create("rhat", "ess"));
  return out;
}
