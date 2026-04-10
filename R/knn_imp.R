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
#' @description
#' Imputes missing values in numeric matrices using full k-nearest neighbor imputation.
#'
#' @details
#' This function performs **column-wise** nearest neighbor imputation.
#'
#' When `dist_pow > 0`, imputed values are computed as distance-weighted averages
#' where weights are inverse distances raised to the power of `dist_pow`.
#'
#' The `tree` parameter (when `TRUE`) uses a BallTree for faster neighbor search
#' via `{mlpack}` but **requires pre-filling** missing values with column means.
#' This can introduce a small bias when missingness is high.
#'
#' @section Performance Optimization:
#' - **`tree = TRUE`** (BallTree): Only use when imputation runtime becomes prohibitive
#'   and missingness is low (<5% missing). Tree construction has overhead.
#' - **`tree = FALSE`** (default, brute-force): Always safe and usually faster for
#'   small-to-moderate data or high-dimensional cases.
#' - **Subset imputation**: Use the `subset` parameter for efficiency when only
#'   specific columns need imputation (e.g., epigenetic clocks CpGs).
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param k Number of nearest neighbors for imputation. 10 is a good starting point.
#' @param colmax A number from 0 to 1. Threshold of missing data above which K-NN imputation is skipped.
#' @param method Either "euclidean" (default) or "manhattan". Distance metric for nearest neighbor calculation.
#' @param cores Number of cores to parallelize over.
#' @param post_imp Whether to impute remaining missing values (those that failed K-NN imputation)
#' using column means (default = `TRUE`).
#' @param subset Character vector of column names or integer vector of column
#' indices specifying which columns to impute.
#' @param dist_pow The amount of penalization for further away nearest neighbors in the weighted average.
#' `dist_pow = 0` (default) is the simple average of the nearest neighbors.
#' @param tree Logical. `FALSE` (default) = brute-force K-NN. `TRUE` = use mlpack BallTree.
#' @param max_cache Maximum allowed cache size in GB (default `4`). When
#' greater than `0`, pairwise distances between columns with missing values
#' are pre-computed and cached, which is faster for moderate-sized data but
#' uses O(m^2) memory where m is the number of columns with missing values.
#' Set to `0` to disable caching and trade speed for lower memory usage on
#' very wide data.
#'
#' @return A numeric matrix of the same dimensions as `obj` with missing values imputed.
#'
#' @references
#' Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu (2002).
#' Diagnosis of multiple cancer types by shrunken centroids of gene expression
#' PNAS 99: 6567-6572. Available at www.pnas.org
#'
#' @examples
#' data(khanmiss1)
#' sum(is.na(khanmiss1))
#'
#' # Basic K-NN imputation (khanmiss1 has genes in rows, so transpose)
#' t_khanmiss1 <- t(khanmiss1)
#' result <- knn_imp(t_khanmiss1, k = 5)
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
  max_cache = 4
) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, null.ok = FALSE, .var.name = "obj")
  checkmate::assert_int(k, lower = 1, upper = ncol(obj) - 1, .var.name = "k")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  checkmate::assert(
    checkmate::check_character(
      subset,
      min.len = 0, max.len = ncol(obj), any.missing = FALSE, unique = TRUE, null.ok = TRUE
    ),
    checkmate::check_integerish(
      subset,
      lower = 1, upper = ncol(obj), min.len = 0, max.len = ncol(obj),
      any.missing = FALSE, null.ok = TRUE, unique = TRUE
    ),
    combine = "or",
    .var.name = "subset"
  )
  stopifnot(length(dist_pow) == 1, dist_pow >= 0, !is.infinite(dist_pow))
  checkmate::assert_flag(tree, .var.name = "tree")
  checkmate::assert_number(max_cache, lower = 0, finite = TRUE, null.ok = FALSE, .var.name = "max_cache")

  # Subset logic
  if (is.null(subset)) {
    subset <- seq_len(ncol(obj))
  } else if (length(subset) == 0) {
    subset <- integer(0)
  } else if (is.character(subset)) {
    if (is.null(colnames(obj))) {
      stop("`subset` contains characters but `obj` doesn't have column names")
    }
    matched_indices <- match(subset, colnames(obj), nomatch = NA)
    if (anyNA(matched_indices)) {
      message("Feature(s) in `subset` not found in `colnames(obj)` and is dropped")
    }
    subset <- matched_indices[!is.na(matched_indices)]
  }
  if (length(subset) == 0) {
    message("No features in subset detected. No imputation was performed.")
    return(obj)
  }
  miss <- is.na(obj)
  cmiss <- colSums(miss)
  miss_rate <- cmiss / nrow(obj)
  if (any(miss_rate == 1)) {
    stop("Col(s) with all missing detected. Remove before proceed")
  }

  # columns with more missing than colmax are not included at all
  eligible <- miss_rate < colmax
  pre_imp_cols <- obj[, eligible, drop = FALSE]
  pre_imp_miss <- miss[, eligible, drop = FALSE]

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
  if (post_imp && anyNA(obj[, subset, drop = FALSE])) {
    na_indices <- which(is.na(obj[, subset, drop = FALSE]), arr.ind = TRUE)
    sub_means <- colMeans(obj[, subset, drop = FALSE], na.rm = TRUE)
    i_vec <- na_indices[, 1]
    jj_vec <- na_indices[, 2]
    obj[cbind(i_vec, subset[jj_vec])] <- sub_means[jj_vec]
  }

  class(obj) <- c("SlideImpImputedMatrix", class(obj))
  attr(obj, "imp_method") <- "knn"
  return(obj)
}
