#' K-Nearest Neighbor Imputation for Numeric Matrices
#'
#' Impute missing values in a numeric matrix using full K-nearest neighbors (K-NN).
#'
#' @param obj A numeric matrix with samples in rows and features in columns.
#' @param k Integer. Number of nearest neighbors to use for K-NN imputation.
#' @param colmax Numeric scalar between `0` and `1`. Columns with a missing-data
#'  proportion greater than `colmax` are excluded from the main imputation
#'  method. Excluded columns are left unchanged unless `post_imp = TRUE`, in
#'  which case remaining missing values are replaced by column means when
#'  possible.
#' @param method Character. K-NN imputation distance method: either `"euclidean"`
#'  or `"manhattan"`.
#' @param cores Integer. Number of cores to use for K-NN imputation. Defaults
#'   to `1`.
#' @param post_imp Logical. If `TRUE`, replace missing values remaining after
#'  the main imputation method with column means when possible.
#' @param subset Optional character or integer vector specifying columns to
#'  target for imputation. If `NULL`, all eligible columns are targeted.
#' @param dist_pow Numeric. Power used to penalize more distant neighbors in
#'  the weighted average, applied to the distance defined by `method`, which
#'  is a mean squared difference for `"euclidean"` and a mean absolute
#'  difference for `"manhattan"`. `dist_pow = 0` gives an unweighted average
#'  of the nearest neighbors. The full weighting rule is given under Details
#'  in [knn_imp()].
#' @param na_check Logical. If `TRUE`, check whether the returned matrix still
#'  contains missing values.
#' @param .progress Logical. If `TRUE`, show imputation progress.
#'
#' @returns A numeric matrix of the same dimensions as `obj`, with missing
#' values imputed. The returned object has class `slideimp_results`.
#'
#' @details
#' `knn_imp()` performs imputation column-wise, treating rows as observations
#' and columns as features.
#'
#' Nearest neighbors are found using brute-force K-NN. The distance between
#' the target column and a candidate column is computed over the rows where
#' both are observed, then divided by the number of those rows, so that
#' candidates overlapping the target on different numbers of rows stay
#' comparable. `method = "euclidean"` averages squared differences and does
#' not take a square root; `method = "manhattan"` averages absolute
#' differences. A candidate that shares no observed row with the target is
#' not a neighbor and is never used, so a column may be imputed from fewer
#' than `k` neighbors, or left missing when nothing overlaps it at all.
#'
#' When `dist_pow > 0`, imputed values are distance-weighted averages of the
#' selected neighbors: a neighbor at distance `d` is weighted
#' `(d_min / d)^dist_pow`, where `d_min` is the smallest positive distance
#' among the neighbors chosen for that column. `d_min` is common to every
#' weight and cancels out of the average; it serves only to keep the largest
#' weight at `1`.
#'
#' Because the euclidean distance above is a mean squared difference,
#' `dist_pow` acts on that squared quantity. Weighting by
#' `1 / d^dist_pow` there is equivalent to weighting by
#' `1 / d^(2 * dist_pow)` on the square-rooted euclidean distance that is
#' more commonly written. A given `dist_pow` is therefore a harsher penalty
#' under `"euclidean"` than under `"manhattan"`, and a value tuned for one
#' `method` does not carry the same meaning under the other.
#'
#' Neighbors are selected once per column, but a neighbor that is itself
#' missing at the row being imputed contributes nothing to that row: it drops
#' out of the weighted sum and out of the total weight, so the average is
#' taken over whichever of the `k` are observed there. The contributing set
#' can therefore differ from row to row within one column, and a row where
#' none of the selected neighbors is observed is left missing.
#'
#' A neighbor at distance zero, one that agrees with the target on every
#' co-observed row, would otherwise take an infinite weight. Such neighbors
#' are given weight `1` while the remaining neighbors are weighted against a
#' small positive floor, so exact matches dominate the average without the
#' other neighbors dropping out of it entirely. When every selected neighbor
#' is at distance zero, all weights are `1`.
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @section K-NN performance optimization:
#' - Use `subset` when only specific columns need imputation.
#' - Use grouped or sliding-window imputation for very large matrices.
#'
#' @references
#' Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
#' Botstein D, Altman RB (2001). Missing value estimation methods for DNA
#' microarrays. *Bioinformatics*, 17(6), 520-525.
#' \doi{10.1093/bioinformatics/17.6.520}
#'
#' @examples
#' set.seed(123)
#' obj <- sim_mat(20, 20, perc_col_na = 1)$input
#' sum(is.na(obj))
#'
#' # Select `k` with `tune_imp()`.
#' result <- knn_imp(obj, k = 10, .progress = FALSE)
#' result
#'
#' @export
knn_imp <- function(
  obj,
  k,
  colmax = 0.9,
  method = c("euclidean", "manhattan"),
  cores = 1,
  post_imp = TRUE,
  subset = NULL,
  dist_pow = 0,
  na_check = TRUE,
  .progress = FALSE
) {
  # Pre-conditioning
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    min.rows = 1,
    min.cols = 2,
    null.ok = FALSE,
    .var.name = "obj"
  )
  check_finite(obj)
  method <- match.arg(method)
  checkmate::assert_int(k, lower = 1, upper = ncol(obj) - 1, .var.name = "k")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  checkmate::assert_number(
    dist_pow,
    lower = 0,
    finite = TRUE,
    .var.name = "dist_pow"
  )
  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_flag(na_check, .var.name = "na_check")

  subset <- resolve_subset(subset, obj)
  if (is.null(subset)) {
    # resolve_subset() has already reported the empty subset. The input is
    # returned unchanged, but still as a slideimp_results as @returns promises.
    return(
      new_slideimp_results(
        obj,
        "knn",
        fallback = FALSE,
        post_imp = post_imp,
        na_check = na_check
      )
    )
  }

  # compute per-column missingness
  cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)

  # early exit if no missingness in subset columns
  if (!any(cmiss[subset] > 0)) {
    cli::cli_inform(
      "No missing values in subset columns. Returning input unchanged."
    )
    return(
      new_slideimp_results(
        obj,
        "knn",
        fallback = FALSE,
        post_imp = post_imp,
        na_check = na_check
      )
    )
  }

  miss_rate <- cmiss / nrow(obj)

  # partitioning: determine eligible columns under colmax threshold.
  # Inclusive, matching pca_imp() and the roxygen ("greater than colmax").
  eligible <- miss_rate <= min(colmax, 1)
  n_elig <- sum(eligible)

  if (k > n_elig - 1L) {
    cli::cli_abort(
      c(
        "{.arg k} ({k}) exceeds usable columns ({n_elig}).",
        "i" = "Reduce {.arg k} or relax {.arg colmax} to admit more columns."
      ),
      class = "slideimp_infeasible"
    )
  }

  eligible_idx <- which(eligible)
  has_miss_idx <- which(cmiss > 0L)
  eligible_has_miss <- intersect(eligible_idx, has_miss_idx)
  eligible_complete <- setdiff(eligible_idx, has_miss_idx)

  grp_impute <- sort(intersect(eligible_has_miss, subset))
  grp_miss_no_imp <- sort(setdiff(eligible_has_miss, subset))
  grp_complete <- sort(eligible_complete)

  if (length(grp_impute) == 0L) {
    cli::cli_abort(
      c(
        "All subset columns with missing values exceed {.arg colmax} ({colmax}).",
        "i" = "Relax {.arg colmax} to admit columns with more missingness."
      ),
      class = "slideimp_infeasible"
    )
  }

  method <- switch(method, "euclidean" = 0L, "manhattan" = 1L)

  # impute using brute-force K-NN
  imputed_values <- impute_knn_brute(
    obj = obj,
    k = k,
    grp_impute = as.integer(grp_impute - 1L),
    grp_miss_no_imp = as.integer(grp_miss_no_imp - 1L),
    grp_complete = as.integer(grp_complete - 1L),
    method = method,
    dist_pow = dist_pow,
    cores = cores,
    pb = .progress
  )

  # convert NaN values back to NA due to C++ handling
  imputed_values[is.nan(imputed_values)] <- NA_real_

  imp_indices <- cbind(imputed_values[, 1], imputed_values[, 2])
  obj[imp_indices] <- imputed_values[, 3]

  # post-imputation: fill any remaining NAs with column means.
  # Only two kinds of subset column can still hold one: a column K-NN was not
  # allowed to touch (missing but over colmax), and a column K-NN could not
  # finish (no candidate shared an observed row, leaving NA in the third
  # result column). Both are known here without another pass over obj, and
  # naming them keeps mean_imp_col() from copying the whole matrix to fill
  # nothing, which is the usual case.
  if (post_imp) {
    na_cols <- sort(unique(c(
      setdiff(intersect(subset, has_miss_idx), grp_impute),
      imputed_values[is.na(imputed_values[, 3]), 2]
    )))
    if (length(na_cols) > 0L) {
      obj <- mean_imp_col(obj, subset = na_cols, cores = cores)
    }
  }

  return(
    new_slideimp_results(
      obj,
      "knn",
      fallback = FALSE,
      post_imp = post_imp,
      na_check = na_check
    )
  )
}
