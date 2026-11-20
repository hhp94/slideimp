#ifndef SLIDEIMP_MATRIX_CHECKS_H
#define SLIDEIMP_MATRIX_CHECKS_H

#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <string>

// Helpers shared by the kernels that read an R matrix in place:
// src/mat_stats.cpp, src/sample_each_rep_cpp.cpp and the K-NN kernel
// (src/imputed_value.h, src/impute_knn_brute.cpp).
//
// Indices that come from R, and counts bounded by nrow() or ncol(), are int.
// An offset into column-major storage is taken in std::size_t, through
// col_ptr(), so a matrix past 2^31 cells is read correctly whatever integer
// type the caller holds its column and row count in.

// The first cell of column `col` in column-major storage of `n_rows` rows.
template <typename T>
inline T *col_ptr(T *base, const std::size_t col, const std::size_t n_rows) {
  return base + col * n_rows;
}

// The value at position i of an integer vector, for a refusal: "NA" for
// NA_integer_, the number otherwise.
static inline std::string int_label(const int v) {
  return v == NA_INTEGER ? std::string("NA") : std::to_string(v);
}

// Refuses a 0-based column index outside 0 .. n_cols - 1. NA_integer_ is
// INT_MIN, so it is refused too.
static inline void stop_on_bad_col_index(const Rcpp::IntegerVector &v,
                                         const int n_cols, const char *name) {
  for (R_xlen_t i = 0; i < v.size(); ++i) {
    if (v[i] < 0 || v[i] >= n_cols) {
      Rcpp::stop("'%s' holds %s at position %d, outside the columns of 'obj' "
                 "(0 .. %d)",
                 name, int_label(v[i]), i + 1, n_cols - 1);
    }
  }
}

// Refuses a negative count. NA_integer_ is INT_MIN, so it is refused too.
static inline void stop_on_negative(const Rcpp::IntegerVector &v,
                                    const char *name) {
  for (R_xlen_t i = 0; i < v.size(); ++i) {
    if (v[i] < 0) {
      Rcpp::stop("'%s' holds %s at position %d. It must be 0 or more", name,
                 int_label(v[i]), i + 1);
    }
  }
}

static inline void stop_on_inf(const Rcpp::NumericMatrix &obj) {
  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());
  const double *x = obj.begin();
  for (std::size_t c = 0; c < n_cols; ++c) {
    const double *col = col_ptr(x, c, n_rows);
    for (std::size_t r = 0; r < n_rows; ++r) {
      if (std::isinf(col[r])) {
        Rcpp::stop("'obj' holds an infinite value at row %d, column %d. "
                   "Infinite values are not supported",
                   r + 1, c + 1);
      }
    }
  }
}

#endif // SLIDEIMP_MATRIX_CHECKS_H
