#ifndef MATRIX_CHECKS_H
#define MATRIX_CHECKS_H

#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <string>

// Shared by check_inf() (src/mat_stats.cpp) and the K-NN kernel
// (src/impute_knn_brute.cpp). Column offsets are size_t, so a matrix past
// 2^31 cells is read correctly.
static inline void stop_on_inf(const Rcpp::NumericMatrix &obj)
{
    const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
    const std::size_t n_cols = static_cast<std::size_t>(obj.ncol());
    const double *x = obj.begin();
    for (std::size_t c = 0; c < n_cols; ++c)
    {
        const double *col = x + c * n_rows;
        for (std::size_t r = 0; r < n_rows; ++r)
        {
            if (std::isinf(col[r]))
            {
                Rcpp::stop(
                    std::string("Infinite value found at row ") +
                    std::to_string(r + 1) +
                    ", column " +
                    std::to_string(c + 1) +
                    ". Infinite values are not supported.");
            }
        }
    }
}

#endif // MATRIX_CHECKS_H
