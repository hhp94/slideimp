// [[Rcpp::depends(BH, bigmemory)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>
#if defined(_OPENMP)
#include <omp.h>
#endif

// This files deals with the bigmem matrices by modifying them in place to avoid
// realizing them in R.

// [[Rcpp::export]]
void bigmem_impute_colmeans(SEXP pBigMat, const std::vector<size_t> &col_indices, int cores = 1)
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  Rcpp::XPtr<BigMatrix> pMat(pBigMat);
  MatrixAccessor<double> mat(*pMat);

  size_t n_cols = col_indices.size();
  size_t n_rows = pMat->nrow();
  size_t total_cols = pMat->ncol();

  if (n_cols == 0 || n_rows == 0)
  {
    return;
  }
  // Convert to 0-based column indices with bounds check
  std::vector<size_t> cpp_cols(n_cols);
  for (size_t j = 0; j < n_cols; ++j)
  {
    if (col_indices[j] == 0 || col_indices[j] > total_cols)
    {
      throw std::runtime_error("Invalid column index");
    }
    cpp_cols[j] = col_indices[j] - 1;
  }
  // Calculate column means for all requested columns
  std::vector<double> col_means(n_cols);
  std::vector<bool> has_na(n_cols);
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t j = 0; j < n_cols; ++j)
  {
    const size_t col = cpp_cols[j];
    double sum = 0.0;
    size_t n_valid = 0;

    for (size_t i = 0; i < n_rows; ++i)
    {
      double val = mat[col][i];
      if (!std::isnan(val))
      {
        sum += val;
        n_valid++;
      }
    }
    col_means[j] = (n_valid > 0) ? sum / n_valid : std::numeric_limits<double>::quiet_NaN();
    has_na[j] = (n_valid < n_rows);
  }
  // Then loop through and impute missing values
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t j = 0; j < n_cols; ++j)
  {
    if (!has_na[j] || std::isnan(col_means[j]))
    {
      continue;
    }
    const size_t col = cpp_cols[j];
    const double mean_val = col_means[j];
    for (size_t i = 0; i < n_rows; ++i)
    {
      if (std::isnan(mat[col][i]))
      {
        mat[col][i] = mean_val;
      }
    }
  }
}

// Divide columns at certain regions by counts_vec
// [[Rcpp::export]]
void bigmem_avg(
    SEXP pBigMat,
    const std::vector<size_t> &start,
    const std::vector<size_t> &end,
    std::vector<int> &counts_vec,
    int cores = 1)
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif
  Rcpp::XPtr<BigMatrix> pMat(pBigMat);
  MatrixAccessor<double> mat(*pMat);
  size_t n_rows = pMat->nrow();
  size_t n_cols = pMat->ncol();

  // Check that counts_vec has the correct length
  if (counts_vec.size() != n_cols)
  {
    throw std::runtime_error("counts_vec must have the same length as number of columns");
  }

  // Check that start and end vectors have the same length
  if (start.size() != end.size())
  {
    throw std::runtime_error("start and end vectors must have the same length");
  }

  // If start is length 0 vec then immediate return
  if (start.size() == 0)
  {
    return;
  }

  // Validation: check all start <= end and find min/max bounds
  size_t min_start = std::numeric_limits<size_t>::max();
  size_t max_end = 0;

  for (size_t i = 0; i < start.size(); ++i)
  {
    // Check for 1-based indexing validity
    if (start[i] == 0 || end[i] == 0)
    {
      throw std::runtime_error("Column indices must be >= 1 (1-based indexing)");
    }

    // Check start <= end
    if (start[i] > end[i])
    {
      throw std::runtime_error("start index must be <= end index");
    }

    // Track min and max for bounds checking
    min_start = std::min(min_start, start[i]);
    max_end = std::max(max_end, end[i]);
  }

  // Bounds check against matrix dimensions
  if (min_start > n_cols || max_end > n_cols)
  {
    throw std::runtime_error("Column indices exceed matrix dimensions");
  }

  // Process each range. Parallel over each range.
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < start.size(); ++i)
  {
    // Convert to 0-based indexing
    const size_t start_col = start[i] - 1;
    const size_t end_col = end[i] - 1;

    for (size_t col = start_col; col <= end_col; ++col)
    {
      const double divisor = static_cast<double>(counts_vec[col]);
      // Divide all values in this column by the denominator
#ifdef _OPENMP
#pragma omp simd
#endif
      for (size_t row = 0; row < n_rows; ++row)
      {
        mat[col][row] /= divisor;
      }
    }
  }
}
