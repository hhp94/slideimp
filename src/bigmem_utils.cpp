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

void bigmem_avg(
    SEXP pBigMat, const std::vector<size_t> &start,
    const std::vector<size_t> &end, std::vector<int> &denom_vec,
    int cores = 1)
{
#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

  Rcpp::XPtr<BigMatrix> pMat(pBigMat);
  MatrixAccessor<double> mat(*pMat);
  size_t n_rows = pMat->nrow();
  size_t n_cols = pMat->ncol();
  // Check that denom_vec has the correct length
  if (denom_vec.size() != n_cols)
  {
    throw std::runtime_error("denom_vec must have the same length as number of columns");
  }
  // Check that start and end vectors have the same length
  if (start.size() != end.size())
  {
    throw std::runtime_error("start and end vectors must have the same length");
  }
  // Process each range. Parallel over range instead of columns.
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (size_t i = 0; i < start.size(); ++i)
  {
    // Convert to 0-based indexing and validate bounds
    if (start[i] == 0 || start[i] > n_cols || end[i] == 0 || end[i] > n_cols)
    {
      throw std::runtime_error("Invalid column index");
    }
    const size_t start_col = start[i] - 1;
    const size_t end_col = end[i] - 1;

    if (start_col > end_col)
    {
      throw std::runtime_error("start index must be <= end index");
    }
    for (size_t col = start_col; col <= end_col; ++col)
    {
      // No need to check for division by zero
      double divisor = static_cast<double>(denom_vec[col]);
      // Divide all values in this column by the denominator
      for (size_t row = 0; row < n_rows; ++row)
      {
        mat[col][row] /= divisor;
      }
    }
  }
}
