#include "imputed_value.h"

// parallelism
#include <RcppThread.h> // RcppThread::parallelFor

// mlpack
#include <mlpack.h>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>

// standard library
#include <cmath>   // std::isnan
#include <cstring> // std::memcpy
#include <stdexcept>

// [[Rcpp::export]]
arma::mat impute_knn_mlpack(
    const arma::mat &obj,
    const arma::uword k,
    const arma::uvec &grp_impute,
    const arma::uvec &grp_miss_no_imp,
    const arma::uvec &grp_complete,
    const int method,
    const double dist_pow,
    int cores = 1)
{
    validate_knn_inputs(
        obj, k, grp_impute, grp_miss_no_imp, grp_complete, method, dist_pow);
    stop_on_inf(obj);
    GroupLayout layout{grp_impute.n_elem, grp_miss_no_imp.n_elem, grp_complete.n_elem};
    const arma::uword n_rows = obj.n_rows;

    // obj_reordered holds ALL three groups (mlpack needs one contiguous matrix
    // for the ball tree). `nmiss_masked` covers only groups 1+2 — group 3 has no
    // missing entries by construction.
    arma::mat obj_reordered(n_rows, layout.n_working());
    MaskMat nmiss_masked(n_rows, layout.n_masked());
    arma::uvec n_col_valid(layout.n_masked(), arma::fill::zeros);

    // mean-fill variant of copy_with_mask: mlpack's distance metrics are
    // mask-unaware, so we need a sensible numeric value in place of NaN.
    // Tracks valid count in the same pass so initialize_result_matrix can
    // derive missing counts
    auto copy_mean_fill = [&](arma::uword local_pos, arma::uword orig_pos)
    {
        const double *src = obj.colptr(orig_pos);
        double *dst = obj_reordered.colptr(local_pos);
        mask_t *mask_dst = nmiss_masked.colptr(local_pos);

        double sum = 0.0;
        arma::uword n_finite = 0;
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            double v = src[r];
            mask_t finite = !std::isnan(v);
            mask_dst[r] = finite;
            n_finite += finite;
            sum += finite ? v : 0.0;
        }
        double mean_val = (n_finite > 0) ? sum / static_cast<double>(n_finite) : 0.0;
        for (arma::uword r = 0; r < n_rows; ++r)
        {
            dst[r] = mask_dst[r] ? src[r] : mean_val;
        }
        n_col_valid(local_pos) = n_finite;
    };

    // group 1
    for (arma::uword i = 0; i < layout.n_imp; ++i)
    {
        copy_mean_fill(i, grp_impute(i));
    }
    arma::uvec col_offsets;
    std::vector<arma::uvec> rows_to_impute_vec;
    arma::mat result = initialize_result_matrix(
        nmiss_masked, grp_impute, layout, n_col_valid, col_offsets, rows_to_impute_vec);

    if (result.n_rows == 0)
    {
        return result;
    }

    // group 2
    for (arma::uword p = 0; p < layout.n_mni; ++p)
    {
        copy_mean_fill(layout.mni_start() + p, grp_miss_no_imp(p));
    }

    // group 3: NaN-free by construction, straight memcpy.
    // Order matters: impute_column_values de-tags group 3 neighbors as
    // grp_complete(local - complete_start), so physical layout here MUST match
    // grp_complete's order.
    for (arma::uword p = 0; p < layout.n_complete; ++p)
        std::memcpy(
            obj_reordered.colptr(layout.complete_start() + p),
            obj.colptr(grp_complete(p)),
            n_rows * sizeof(double));

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

    cores = std::max(1, cores);
    const size_t n_threads = static_cast<size_t>(cores);
    const size_t n_batches = static_cast<size_t>(cores);

    RcppThread::parallelFor(
        0,
        layout.n_imp,
        [&](arma::uword i)
        {
            if (resultingNeighbors.n_rows <= 1)
            {
                return;
            }
            const arma::uword actual_neighbors = std::min<arma::uword>(k, resultingNeighbors.n_rows - 1);
            if (actual_neighbors == 0)
            {
                return;
            }

            arma::uvec nn_columns = resultingNeighbors(arma::span(1, actual_neighbors), i);
            arma::vec nn_dists = resultingDistances(arma::span(1, actual_neighbors), i);
            arma::vec weights = 1.0 / arma::pow(nn_dists + epsilon, dist_pow);

            impute_column_values(
                result, obj_reordered, nmiss_masked, layout,
                col_offsets(i),
                nn_columns, weights,
                rows_to_impute_vec[i],
                obj, grp_complete);
        },
        n_threads, n_batches);

    return result;
}
