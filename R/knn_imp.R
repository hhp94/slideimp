#' K-Nearest Neighbors Imputation for Numeric Matrices
#'
#' Impute missing values in a numeric matrix using k-nearest neighbors
#' (K-NN).
#'
#' @details
#' `knn_imp()` performs imputation column-wise, treating rows as observations
#' and columns as features.
#'
#' When `dist_pow > 0`, imputed values are computed as distance-weighted
#' averages. Weights are inverse distances raised to the power of `dist_pow`.
#'
#' If `tree = TRUE`, nearest neighbors are found with a ball tree via the
#' `mlpack` package. This can be faster for some large, low-missingness data
#' sets, but it requires initially filling missing values with column means,
#' which can introduce bias when missingness is high.
#'
#' @section Performance optimization:
#' - `tree = FALSE` uses brute-force K-NN. This is always safe and is often
#'   faster for small to moderate data sets or high-dimensional data.
#' - `tree = TRUE` uses ball-tree K-NN. Consider this only when run time is
#'   prohibitive and missingness is low, for example less than 5%.
#' - Use `subset` when only specific columns need imputation.
#'
#' @param obj A numeric matrix with samples in rows and features in columns.
#' @param k Integer. Number of nearest neighbors to use for imputation.
#' @param colmax Numeric scalar between `0` and `1`. Columns with a missing-data
#'   proportion greater than `colmax` are not imputed.
#' @param method Character. Distance metric for nearest-neighbor calculation:
#'   either `"euclidean"` or `"manhattan"`.
#' @param cores Integer. Number of cores to use for parallel computation.
#'   Defaults to `1`.
#' @param post_imp Logical. If `TRUE`, replace any remaining missing values
#'   with column means after K-NN imputation.
#' @param subset Optional character or integer vector specifying columns to
#'   impute. If `NULL`, all eligible columns are imputed.
#' @param dist_pow Numeric. Power used to penalize more distant neighbors in
#'   the weighted average. `dist_pow = 0` gives an unweighted average of the
#'   nearest neighbors.
#' @param tree Logical. If `FALSE`, use brute-force K-NN. If `TRUE`, use
#'   ball-tree K-NN via `mlpack`.
#' @param na_check Logical. If `TRUE`, check whether the result still contains
#'   missing values.
#' @param .progress Logical. If `TRUE`, show imputation progress.
#'
#' @returns A numeric matrix of the same dimensions as `obj`, with missing
#' values imputed. The returned object has class `slideimp_results`.
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
  tree = FALSE,
  na_check = TRUE,
  .progress = TRUE
) {
  # Pre-conditioning
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, null.ok = FALSE, .var.name = "obj")
  check_finite(obj)
  method <- match.arg(method)
  checkmate::assert_int(k, lower = 1, upper = ncol(obj) - 1, .var.name = "k")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  stopifnot(length(dist_pow) == 1, dist_pow >= 0, !is.infinite(dist_pow))
  checkmate::assert_flag(tree, .var.name = "tree")
  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_flag(na_check, .var.name = "na_check")

  subset <- resolve_subset(subset, obj)
  if (is.null(subset)) {
    cli::cli_inform("No columns to impute. Returning input unchanged.")
    return(obj)
  }

  # compute per-column missingness
  cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)

  # early exit if no missingness in subset columns
  if (!any(cmiss[subset] > 0)) {
    cli::cli_inform("No missing values in subset columns. Returning input unchanged.")
    return(obj)
  }

  miss_rate <- cmiss / nrow(obj)

  # partitioning: determine eligible columns (under colmax threshold)
  eligible <- miss_rate < min(colmax, 1)
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

  # Column groups as ORIGINAL (0-based) indices into obj.
  # The C++ side will use cmiss to know which columns have missingness and
  # skip/handle accordingly; it also uses `eligible` to restrict the working set.
  eligible_idx <- which(eligible) # 1-based orig indices of eligible cols
  has_miss_idx <- which(cmiss > 0L) # 1-based orig indices with any NA
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

  method <- switch(method,
    "euclidean" = 0L,
    "manhattan" = 1L
  )

  # Impute
  if (!tree) {
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
  } else {
    imputed_values <- impute_knn_mlpack(
      obj = obj,
      k = k,
      grp_impute = as.integer(grp_impute - 1L),
      grp_miss_no_imp = as.integer(grp_miss_no_imp - 1L),
      grp_complete = as.integer(grp_complete - 1L),
      method = method,
      dist_pow = dist_pow,
      cores = cores
    )
  }

  # convert NaN values back to NA due to C++ handling
  imputed_values[is.nan(imputed_values)] <- NA_real_

  # Column indices returned by C++ are already original-matrix indices now,
  # so no remapping through orig_indices is needed.
  imp_indices <- cbind(imputed_values[, 1], imputed_values[, 2])
  obj[imp_indices] <- imputed_values[, 3]

  # post-imputation: fill any remaining NAs with column means
  if (post_imp) {
    obj <- mean_imp_col(obj, subset = subset, cores = cores)
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
