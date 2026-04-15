check_cache_memory <- function(n_miss_cols, max_cache) {
  if (n_miss_cols <= 1L) {
    return(invisible(NULL))
  }
  n <- as.numeric(n_miss_cols)
  cache_gb <- (n * (n - 1) / 2) * 8 / 1024^3
  if (cache_gb > max_cache) {
    stop(sprintf(
      paste0(
        "Cache would require %.1f GB for %d missing columns, ",
        "which exceeds `max_cache` (%.1f GB). ",
        "Increase `max_cache` or set `max_cache = 0` to disable the cache."
      ),
      cache_gb, n_miss_cols, max_cache
    ))
  }
  invisible(NULL)
}

#' K-Nearest Neighbor Imputation for Numeric Matrices
#'
#' Imputes missing values in a numeric matrix using k-nearest neighbors (KNN).
#'
#' @details
#' This function performs imputation **column-wise** (using rows as observations).
#'
#' When `dist_pow > 0`, imputed values are computed as distance-weighted averages
#' where weights are inverse distances raised to the power of `dist_pow`.
#'
#' The `tree` parameter (when `TRUE`) uses a BallTree for faster neighbor search
#' via `{mlpack}` but **requires pre-filling** missing values with column means.
#' This can introduce a small bias when missingness is high.
#'
#' @section Performance Optimization:
#' - **`tree = FALSE`** (default, brute-force KNN): Always safe and usually faster for
#'   small-to-moderate data or high-dimensional cases.
#' - **`tree = TRUE`** (BallTree KNN): Only use when imputation runtime becomes prohibitive
#'   and missingness is low (<5% missing). Tree construction has overhead.
#' - **Subset imputation**: Use the `subset` parameter for efficiency when only
#'   specific columns need imputation (e.g., epigenetic clocks CpGs).
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param k Number of nearest neighbors for imputation. 10 is a good starting point.
#' @param colmax A number from 0 to 1. Threshold of column-wise missing data rate above which imputation is skipped.
#' @param method Either "euclidean" (default) or "manhattan". Distance metric for nearest neighbor calculation.
#' @param cores Number of cores for KNN parallelization (OpenMP). On macOS, OpenMP may need additional compiler configuration.
#' @param post_imp Whether to impute remaining missing values (those that failed imputation)
#' using column means.
#' @param subset Character vector of column names or integer vector of column
#' indices specifying which columns to impute.
#' @param dist_pow The amount of penalization for further away nearest neighbors in the weighted average.
#' `dist_pow = 0` (default) is the simple average of the nearest neighbors.
#' @param tree Logical. `FALSE` (default) = brute-force K-NN. `TRUE` = use `{mlpack}` BallTree.
#' @param max_cache Maximum allowed cache size in GB (default `4`). When
#' greater than `0`, pairwise distances between columns with missing values
#' are pre-computed and cached, which is faster for moderate-sized data but
#' uses O(m^2) memory where m is the number of columns with missing values.
#' Set to `0` to disable caching and trade speed for lower memory usage.
#' @param na_check Boolean. Check if there are `NA` leftover in the results or not.
#'
#' @returns A numeric matrix of the same dimensions as `obj` with missing values imputed.
#'
#' @references
#' Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
#' Botstein D, Altman RB (2001).
#' Missing value estimation methods for DNA microarrays.
#' Bioinformatics 17(6): 520-525.
#'
#' @examples
#' # Basic K-NN imputation
#' obj <- sim_mat(20, 20, perc_col_na = 1)$input
#' sum(is.na(obj))
#' result <- knn_imp(obj, k = 10)
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
  max_cache = 4,
  na_check = TRUE
) {
  # Pre-conditioning
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, null.ok = FALSE, .var.name = "obj")
  if (!anyNA(obj)) {
    return(obj)
  }
  method <- match.arg(method)
  checkmate::assert_int(k, lower = 1, upper = ncol(obj) - 1, .var.name = "k")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  stopifnot(length(dist_pow) == 1, dist_pow >= 0, !is.infinite(dist_pow))
  checkmate::assert_flag(tree, .var.name = "tree")
  checkmate::assert_flag(na_check, .var.name = "na_check")
  checkmate::assert_number(max_cache, lower = 0, finite = TRUE, null.ok = FALSE, .var.name = "max_cache")

  subset <- resolve_subset(subset, obj)
  if (is.null(subset)) {
    return(obj)
  }

  # partitioning
  miss <- is.na(obj)
  cmiss <- colSums(miss)
  miss_rate <- cmiss / nrow(obj)

  eligible <- miss_rate < min(colmax, 1)
  pre_imp_cols <- obj[, eligible, drop = FALSE]
  pre_imp_miss <- miss[, eligible, drop = FALSE]

  n_elig <- ncol(pre_imp_cols)
  if (k > n_elig - 1L) {
    if (post_imp) {
      cli::cli_inform(
        "{.arg k} ({k}) exceeds usable columns ({n_elig}). Falling back to mean imputation."
      )
      obj <- mean_imp_col(obj, subset = subset, cores = cores)
    }
    # post_imp ran mean_imp_col but colmax-excluded cols may still have NAs -> check
    # no post_imp -> definitely still has NAs -> skip check
    return(
      as_slideimp_results(
        obj,
        "knn",
        fallback = TRUE,
        post_imp = post_imp,
        na_check = na_check && post_imp
      )
    )
  }

  # column groups (1-based local indices into pre_imp_cols)
  orig_indices <- which(eligible)
  local_cmiss <- cmiss[eligible]
  local_has_miss <- which(local_cmiss > 0L)
  local_in_subset <- which(orig_indices %in% subset)
  grp_impute <- sort(intersect(local_has_miss, local_in_subset))
  grp_miss_no_imp <- sort(setdiff(local_has_miss, local_in_subset))
  grp_complete <- which(local_cmiss == 0L)

  cache <- max_cache > 0
  if (cache) {
    check_cache_memory(
      n_miss_cols = length(grp_impute),
      max_cache = max_cache
    )
  }

  # Impute
  if (!tree) {
    # Important: pre-fill with zero (required for auto-vectorization in brute-force path)
    pre_imp_cols[pre_imp_miss] <- 0.0
    imputed_values <- impute_knn_brute(
      obj = pre_imp_cols,
      nmiss = !pre_imp_miss,
      k = k,
      grp_impute = as.integer(grp_impute - 1L),
      grp_miss_no_imp = as.integer(grp_miss_no_imp - 1L),
      grp_complete = as.integer(grp_complete - 1L),
      method = switch(method,
        "euclidean" = 0L,
        "manhattan" = 1L
      ),
      dist_pow = dist_pow,
      cores = cores,
      cache = cache
    )
  } else {
    imputed_values <- impute_knn_mlpack(
      obj = mean_imp_col(pre_imp_cols),
      nmiss = !pre_imp_miss,
      k = k,
      grp_impute = as.integer(grp_impute - 1L),
      method = switch(method,
        "euclidean" = 0L,
        "manhattan" = 1L
      ),
      dist_pow = dist_pow,
      cores = cores
    )
  }

  # convert NaN values back to NA due to C++ handling
  imputed_values[is.nan(imputed_values)] <- NA_real_

  # map column indices from pre_imp_cols back to original matrix columns
  imp_indices <- cbind(imputed_values[, 1], orig_indices[imputed_values[, 2]])
  obj[imp_indices] <- imputed_values[, 3]

  # post-imputation: fill any remaining NAs with column means
  if (post_imp) {
    obj <- mean_imp_col(obj, subset = subset, cores = cores)
  }

  return(
    as_slideimp_results(
      obj,
      "knn",
      fallback = FALSE,
      post_imp = post_imp,
      na_check = na_check
    )
  )
}
