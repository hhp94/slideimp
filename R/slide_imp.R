# Given a start and end vectors. Give the counts_vec that counts the number
# of times that position is covered and also the regions where counts are > 1
find_overlap_regions <- function(start, end) {
  max_pos <- max(end)
  delta <- rep(0, max_pos + 1)
  delta[start] <- delta[start] + 1
  valid_ends <- end + 1 <= max_pos + 1
  delta[end[valid_ends] + 1] <- delta[end[valid_ends] + 1] - 1
  counts_vec <- cumsum(delta)[seq_len(max_pos)]
  # extract regions where counts > 1
  overlaps <- counts_vec > 1
  rle_over <- rle(overlaps)
  # compute start and end positions of the TRUE runs
  ends <- cumsum(rle_over$lengths)
  starts <- ends - rle_over$lengths + 1
  # return the region and counts_vec
  return(
    list(
      region = cbind(
        start = starts[rle_over$values],
        end = ends[rle_over$values]
      ),
      counts_vec = counts_vec
    )
  )
}

#' Sliding Window K-NN or PCA Imputation
#'
#' @description
#' Performs sliding window K-NN or PCA imputation of large numeric matrices column-wise. This
#' method assumes that columns are meaningfully sorted by `location`.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#' @param location A sorted numeric vector of length `ncol(obj)` giving the
#'   position of each column (e.g., genomic coordinates). Used to define
#'   sliding windows.
#' @param window_size Window width in the same units as `location`.
#' @param overlap_size Overlap between consecutive windows in the same units
#'   as `location`. Must be less than `window_size`. Default is `0`.
#' @param min_window_n Minimum number of columns a window must contain to be
#'   imputed. Windows smaller than this are dropped. `k` and `ncp` must also
#'   be smaller than `min_window_n`.
#' @param knn_method Either "euclidean" (default) or "manhattan". Distance metric for nearest neighbor calculation.
#' @param pca_method "regularized" by default or "EM".
#' @param .progress Show progress bar (default = `TRUE`).
#'
#' @details
#' The sliding window approach divides the input matrix into smaller, overlapping
#' segments based on `location` values and applies imputation to each window
#' independently. Values in overlapping areas are averaged across windows to
#' produce the final imputed result.
#'
#' Windows are determined by `find_windows()`, which walks along the sorted
#' `location` vector and groups columns whose positions fall within
#' `window_size` of each other.
#'
#' Specify `k` and related arguments to use K-NN, `ncp` and related arguments for PCA.
#'
#' @inherit knn_imp return
#'
#' @examples
#' # Generate sample data with missing values with 20 samples and 100 columns
#' # where the column order is sorted (i.e., by genomic position)
#' set.seed(1234)
#' beta_matrix <- t(sim_mat(100, 20)$input)
#' location <- 1:100
#'
#' # Sliding Window K-NN imputation by specifying `k`
#' imputed_knn <- slide_imp(
#'   beta_matrix,
#'   location = location,
#'   k = 5,
#'   window_size = 50,
#'   overlap_size = 10,
#'   scale = FALSE # This argument belongs to PCA imputation and will be ignored
#' )
#' imputed_knn
#'
#' # Sliding Window PCA imputation by specifying `ncp`
#' pca_knn <- slide_imp(
#'   beta_matrix,
#'   location = location,
#'   ncp = 2,
#'   window_size = 50,
#'   overlap_size = 10
#' )
#' pca_knn
#'
#' @export
slide_imp <- function(
  obj,
  location,
  window_size,
  overlap_size = 0,
  min_window_n,
  # KNN-specific parameters
  k = NULL,
  colmax = 0.9,
  knn_method = c("euclidean", "manhattan"),
  cores = 1,
  post_imp = TRUE,
  dist_pow = 0,
  subset = NULL,
  # PCA-specific parameters
  ncp = NULL,
  scale = TRUE,
  pca_method = c("regularized", "EM"),
  coeff.ridge = 1,
  seed = NULL,
  row.w = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  # Others
  .progress = TRUE
) {
  # minimal pre-conditioning to avoid code fragility
  if (sum(c(is.null(k), is.null(ncp))) != 1L) {
    stop("Specify either 'k' for K-NN imputation or 'ncp' for PCA imputation. Not both nor neither.")
  }
  imp_method <- if (!is.null(k)) {
    "knn"
  } else {
    "pca"
  }
  # Pre-conditioning ----
  checkmate::assert_matrix(obj, mode = "numeric", col.names = "unique", null.ok = FALSE, .var.name = "obj")
  cn <- colnames(obj)
  checkmate::assert_numeric(location,
    len = ncol(obj), any.missing = FALSE, sorted = TRUE,
    finite = TRUE, null.ok = FALSE, .var.name = "location"
  )
  checkmate::assert_number(window_size,
    lower = .Machine$double.eps, finite = TRUE,
    null.ok = FALSE, .var.name = "window_size"
  )
  checkmate::assert_number(overlap_size,
    lower = 0, finite = TRUE,
    null.ok = FALSE, .var.name = "overlap_size"
  )
  if (overlap_size >= window_size) {
    stop("`overlap_size` must be strictly less than `window_size`.")
  }
  checkmate::assert_int(min_window_n, lower = 2L, null.ok = FALSE, .var.name = "min_window_n")
  if (imp_method == "knn") {
    knn_method <- match.arg(knn_method)
    checkmate::assert_int(k, lower = 1, upper = min_window_n - 1L, null.ok = FALSE, .var.name = "k")
  } else if (imp_method == "pca") {
    pca_method <- match.arg(pca_method)
    checkmate::assert_int(ncp,
      lower = 1, upper = min(min_window_n - 1L, min(nrow(obj), ncol(obj)) - 1),
      .var.name = "ncp"
    )
    obj_vars <- col_vars(obj)
    if (any(obj_vars < .Machine$double.eps | is.na(obj_vars))) {
      stop("Features with zero variance after na.rm not permitted for PCA Imputation. Try 'col_vars(obj)'")
    }
    rm(obj_vars)
  }
  checkmate::assert_flag(.progress, .var.name = ".progress", null.ok = FALSE)

  # Windowing Logic ----
  windows <- find_windows(location, window_size, overlap_size)
  start <- windows$start
  end <- windows$end

  # drop windows that are too small. Temporary solution
  window_n <- end - start + 1L
  keep <- window_n >= min_window_n
  if (!any(keep)) {
    stop(
      "All windows have fewer than `min_window_n` (", min_window_n,
      ") columns. Consider increasing `window_size` or decreasing `min_window_n`."
    )
  }
  if (any(!keep) && .progress) {
    message(sprintf("Dropping %d window(s) with fewer than %d columns.", sum(!keep), min_window_n))
  }
  start <- start[keep]
  end <- end[keep]

  # PCA will always impute all the values
  if (imp_method == "knn") {
    if (!is.null(subset)) {
      if (is.character(subset)) {
        stopifnot("`subset` are characters but `obj` doesn't have colnames" = !is.null(cn))
        matched <- match(subset, cn, nomatch = NA)
        subset <- matched[!is.na(matched)]
        subset <- sort(subset)
      } else {
        subset <- as.integer(subset)
        if (any(!subset %in% seq_len(ncol(obj)))) {
          stop("Invalid indices in `subset`: must be between 1 and ncol(obj)")
        }
        subset <- sort(unique(subset))
      }
      if (length(subset) == 0) {
        stop("`subset` is not found in `colnames(obj)`")
      }
    } else {
      subset <- seq_len(ncol(obj))
    }
    # Calculate where the subset indices in each window
    subset_list <- lapply(seq_along(start), function(i) {
      first <- findInterval(start[i] - 1, subset) + 1
      last <- findInterval(end[i], subset)
      if (first <= last) {
        subset[first:last] - start[i] + 1
      } else {
        integer(0)
      }
    })
  }
  # Overlap regions to average over
  overlap <- find_overlap_regions(start, end)

  # Sliding Imputation ----
  result <- matrix(0, nrow = nrow(obj), ncol = ncol(obj), dimnames = list(rownames(obj), colnames(obj)))
  if (.progress) {
    message("Step 1/2: Imputing")
    n_windows <- length(start)
    n_steps <- max(1, round(n_windows / 20))
  }

  for (i in seq_along(start)) {
    if (.progress && (i %% n_steps == 0 || i == n_windows || i == 1)) {
      message(sprintf(" Processing window %d of %d", i, n_windows))
    }
    window_cols <- start[i]:end[i]
    if (imp_method == "knn") {
      imputed_window <- knn_imp(
        obj = obj[, window_cols, drop = FALSE],
        k = k,
        colmax = colmax,
        cores = cores,
        method = knn_method,
        post_imp = post_imp,
        dist_pow = dist_pow,
        subset = subset_list[[i]]
      )
    } else if (imp_method == "pca") {
      imputed_window <- pca_imp(
        obj = obj[, window_cols, drop = FALSE],
        ncp = ncp,
        scale = scale,
        method = pca_method,
        coeff.ridge = coeff.ridge,
        seed = seed,
        nb.init = nb.init,
        maxiter = maxiter,
        miniter = miniter,
        row.w = row.w
      )
    }

    result[, window_cols] <- result[, window_cols] + imputed_window
  }

  if (.progress) {
    message("Step 2/2: Averaging overlapping regions")
  }
  counts_vec <- overlap$counts_vec
  # Pad to ncol(obj) in case trailing columns are uncovered
  if (length(counts_vec) < ncol(obj)) {
    counts_vec <- c(counts_vec, rep(0L, ncol(obj) - length(counts_vec)))
  }
  # Average only where multiple windows overlap
  gt_1_idx <- which(counts_vec > 1)
  if (length(gt_1_idx) > 0) {
    result[, gt_1_idx] <- sweep(result[, gt_1_idx, drop = FALSE], 2, counts_vec[gt_1_idx], "/")
  }
  # Restore original values for columns not covered by any window
  uncovered <- which(counts_vec == 0)
  if (length(uncovered) > 0) {
    result[, uncovered] <- obj[, uncovered]
    if (.progress) {
      message(sprintf("Note: %d column(s) not covered by any window; original values retained.", length(uncovered)))
    }
  }

  # Post-imputation ----
  if (imp_method == "knn" && post_imp) {
    if (.progress) {
      message("Post-imputation: filling remaining NAs with column means")
    }
    if (anyNA(result[, subset, drop = FALSE])) {
      na_indices <- which(is.na(result[, subset, drop = FALSE]), arr.ind = TRUE)
      sub_means <- colMeans(result[, subset, drop = FALSE], na.rm = TRUE)
      i_vec <- na_indices[, 1]
      jj_vec <- na_indices[, 2]
      j_vec <- subset[jj_vec]
      result[cbind(i_vec, j_vec)] <- sub_means[jj_vec]
    }
  }

  class(result) <- c("ImputedMatrix", class(result))
  attr(result, "imp_method") <- imp_method
  return(result)
}
