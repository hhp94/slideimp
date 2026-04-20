#include "imputed_value.h"
#include <mlpack.h>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <cmath>
#include <cstring>
#include <stdexcept>
#if defined(_OPENMP)
#include <omp.h>
#endif

// [[Rcpp::export]]
arma::mat impute_knn_mlpack(
    const arma::mat &obj,
    const arma::uword k,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const int method,
    const double dist_pow,
    const int cores = 1)
{
    GroupLayout layout{grp_impute.n_elem, grp_miss_no_imp.n_elem, grp_complete.n_elem};
    const arma::uword n_rows = obj.n_rows;

    arma::mat obj_reordered(n_rows, layout.n_working());
    arma::mat nmiss_masked(n_rows, layout.n_imp + layout.n_mni);

    // mean-fill variant of copy_with_mask
    auto copy_mean_fill = [&](arma::uword local_pos, arma::uword orig_pos)
    {
        const double *src = obj.colptr(orig_pos);
        double *dst = obj_reordered.colptr(local_pos);
        double *mask_dst = nmiss_masked.colptr(local_pos);

        // first pass: accumulate mean over finite entries, write mask
        double sum = 0.0;
        arma::uword n_finite = 0;
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            double v = src[r];
            bool finite = std::isfinite(v);
            mask_dst[r] = finite ? 1.0 : 0.0;
            if (finite)
            {
                sum += v;
                ++n_finite;
            }
        }
        // second pass: write mean-filled values to dst
        double mean_val = (n_finite > 0) ? sum / static_cast<double>(n_finite) : 0.0;
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            double v = src[r];
            dst[r] = std::isfinite(v) ? v : mean_val;
        }
    };

    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        copy_mean_fill(i, grp_impute(i));
    }
    arma::uvec col_offsets;
    std::vector<arma::uvec> rows_to_impute_vec;
    arma::mat result = initialize_result_matrix(
        nmiss_masked, grp_impute, layout, col_offsets, rows_to_impute_vec);
    if (result.n_rows == 0)
    {
        return result;
    }
    for (arma::uword p = 0; p < layout.n_mni; ++p)
    {
        copy_mean_fill(layout.mni_start() + p, grp_miss_no_imp(p));
    }
    // group 3: NaN-free by construction, straight memcpy
    for (arma::uword p = 0; p < layout.n_complete; ++p)
    {
        std::memcpy(
            obj_reordered.colptr(layout.complete_start() + p),
            obj.colptr(grp_complete(p)),
            n_rows * sizeof(double));
    }
    // query matrix = group 1 columns (the ones we're imputing)
    arma::mat query_mat = obj_reordered.cols(0, layout.n_imp - 1);

    arma::umat resultingNeighbors;
    arma::mat resultingDistances;

    if (method == 0)
    {
        using KNNType = mlpack::NeighborSearch<
            mlpack::NearestNeighborSort,
            mlpack::EuclideanDistance,
            arma::mat,
            mlpack::BallTree>;
        KNNType nn(obj_reordered);
        nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
    }
    else if (method == 1)
    {
        using KNNType = mlpack::NeighborSearch<
            mlpack::NearestNeighborSort,
            mlpack::ManhattanDistance,
            arma::mat,
            mlpack::BallTree>;
        KNNType nn(obj_reordered);
        nn.Search(query_mat, k + 1, resultingNeighbors, resultingDistances);
    }
    else
    {
        throw std::invalid_argument("Invalid method: 0=Euclid, 1=Manhattan");
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores) schedule(dynamic)
#endif
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        arma::uword actual_neighbors = std::min(k, resultingNeighbors.n_rows - 1);
        if (actual_neighbors == 0)
            continue;

        arma::uvec nn_columns = resultingNeighbors(arma::span(1, actual_neighbors), i);
        arma::vec nn_dists = resultingDistances(arma::span(1, actual_neighbors), i);
        arma::vec weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);

        impute_column_values(
            result, obj_reordered, nmiss_masked, layout,
            col_offsets(i),
            nn_columns, weights,
            rows_to_impute_vec[i]);
    }

    return result;
}
