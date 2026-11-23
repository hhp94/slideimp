#include <Rcpp.h>

#include "matrix_checks.h"

#include <cmath>
#include <cstddef>

// The column scan: everything the entry points read off their input, in one pass
// over the columns `cols` (0-based, in the order the results come back) of a
// double matrix, read in place. Indices and offsets follow
// src/matrix_checks.h.
//
// Per column read: the missing count (NA and NaN), whether it holds Inf or
// -Inf, and, with `want_var`, the sample variance by Welford's update. Inf
// is not missing: it enters the arithmetic, so a column holding one has
// variance NaN, as in base R, and a column with fewer than two observed
// values has variance NA. With `want_rows`, the missing count of each row
// over the columns read. With `tails`, a pair (below, above), the count of
// each column's cells below `below` and above `above`; NA and NaN are in
// neither, as with colSums(obj < below, na.rm = TRUE).
//
// The tail counts take a loop of their own over each column, right after the
// column's scan loop, so the column is read again from cache, and the
// scan loops keep no branch for them.
//
// It counts and never refuses: each caller applies its own policy in R
// (R/col_scan.R). A matrix that is not double is refused, not converted, since
// converting it would copy it.
// [[Rcpp::export]]
Rcpp::List col_scan_cpp(SEXP obj, const Rcpp::IntegerVector &cols,
                        const bool want_var, const bool want_rows,
                        const Rcpp::NumericVector &tails) {
  if (TYPEOF(obj) != REALSXP || !Rf_isMatrix(obj)) {
    Rcpp::stop("internal error: col_scan_cpp() takes a double matrix");
  }
  const R_xlen_t n_tails = tails.size();
  if (n_tails != 0 && n_tails != 2) {
    Rcpp::stop("internal error: col_scan_cpp() takes no tails or two");
  }
  const bool want_tails = n_tails == 2;
  const double below = want_tails ? tails[0] : 0.0;
  const double above = want_tails ? tails[1] : 0.0;
  const Rcpp::NumericMatrix x(obj);
  stop_on_bad_col_index(cols, x.ncol(), "cols");
  const std::size_t n_rows = static_cast<std::size_t>(x.nrow());
  const R_xlen_t n_read = cols.size();

  Rcpp::IntegerVector n_miss(Rcpp::no_init(n_read));
  Rcpp::LogicalVector has_inf(Rcpp::no_init(n_read));
  Rcpp::NumericVector var(Rcpp::no_init(want_var ? n_read : 0));
  Rcpp::IntegerVector row_miss(want_rows ? x.nrow() : 0); // zero-filled
  Rcpp::IntegerVector n_below(Rcpp::no_init(want_tails ? n_read : 0));
  Rcpp::IntegerVector n_above(Rcpp::no_init(want_tails ? n_read : 0));

  const double *base = x.begin();
  int *nm = n_miss.begin();
  int *hi = has_inf.begin();
  double *vo = var.begin();
  int *rm = row_miss.begin();
  int *n_lo = n_below.begin();
  int *n_up = n_above.begin();

  for (R_xlen_t k = 0; k < n_read; ++k) {
    const double *col =
        col_ptr(base, static_cast<std::size_t>(cols[k]), n_rows);

    bool inf = false;

    if (want_var) {
      std::size_t n = 0;
      double mean = 0.0;
      double mean2 = 0.0;
      for (std::size_t i = 0; i < n_rows; ++i) {
        const double v = col[i];
        if (std::isnan(v)) {
          if (want_rows) {
            ++rm[i];
          }
          continue;
        }
        ++n;
        inf = inf || std::isinf(v);
        const double delta = v - mean;
        mean += delta / static_cast<double>(n);
        const double delta2 = v - mean;
        mean2 += delta * delta2;
      }
      nm[k] = static_cast<int>(n_rows - n);
      vo[k] = (n < 2) ? NA_REAL : mean2 / static_cast<double>(n - 1);
    } else {
      // Counting alone takes no branch on a hole. With one, a matrix whose
      // holes are scattered counted slower than this loop does.
      int miss = 0;
      for (std::size_t i = 0; i < n_rows; ++i) {
        const double v = col[i];
        const int hole = std::isnan(v);
        miss += hole;
        inf = inf | std::isinf(v);
        if (want_rows) {
          rm[i] += hole;
        }
      }
      nm[k] = miss;
    }
    hi[k] = inf ? 1 : 0;

    if (want_tails) {
      // a comparison with NaN is false, so a hole counts in neither tail
      int lo = 0;
      int up = 0;
      for (std::size_t i = 0; i < n_rows; ++i) {
        lo += col[i] < below;
        up += col[i] > above;
      }
      n_lo[k] = lo;
      n_up[k] = up;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("n_miss") = n_miss, Rcpp::Named("has_inf") = has_inf,
      Rcpp::Named("var") = want_var ? static_cast<SEXP>(var) : R_NilValue,
      Rcpp::Named("row_miss") =
          want_rows ? static_cast<SEXP>(row_miss) : R_NilValue,
      Rcpp::Named("n_below") =
          want_tails ? static_cast<SEXP>(n_below) : R_NilValue,
      Rcpp::Named("n_above") =
          want_tails ? static_cast<SEXP>(n_above) : R_NilValue);
}
