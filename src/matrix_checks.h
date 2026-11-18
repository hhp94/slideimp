#ifndef MATRIX_CHECKS_H
#define MATRIX_CHECKS_H

// Rcpp.h, not RcppArmadillo.h: the helpers here take a pointer to
// column-major storage and its dims, so Armadillo callers (memptr(), n_rows,
// n_cols) and Rcpp ones (begin(), nrow(), ncol()) share them. A translation
// unit that also uses Armadillo includes RcppArmadillo.h before this header,
// since RcppArmadillo.h refuses to come after Rcpp.h.
#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <string>

static inline void stop_on_inf(const double *x,
                               const std::size_t n_rows,
                               const std::size_t n_cols)
{
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
