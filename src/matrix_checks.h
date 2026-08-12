#ifndef MATRIX_CHECKS_H
#define MATRIX_CHECKS_H

#include <RcppArmadillo.h>
#include <cmath>
#include <string>

static inline void stop_on_inf(const arma::mat &obj)
{
    const arma::uword n_rows = obj.n_rows;
    const arma::uword n_cols = obj.n_cols;

    for (arma::uword c = 0; c < n_cols; ++c)
    {
        const double *col = obj.colptr(c);
        for (arma::uword r = 0; r < n_rows; ++r)
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
