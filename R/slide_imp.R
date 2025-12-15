# Given a start and end vectors. Give the counts_vec that counts the number
# of times that position is covered and also the regions where counts are > 1
find_overlap_regions <- function(start, end) {
  max_pos <- max(end)
  delta <- rep(0, max_pos + 1)
  delta[start] <- delta[start] + 1
  valid_ends <- end + 1 <= max_pos + 1
  delta[end[valid_ends] + 1] <- delta[end[valid_ends] + 1] - 1
  counts_vec <- cumsum(delta)[seq_len(max_pos)]
  storage.mode(counts_vec) <- "integer"
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
#' Performs K-NN or PCA imputation large numeric matrices using a sliding
#' window approach column-wise. This method assumes that columns are meaningfully sorted.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#' @param n_feat Number of features in a window.
#' @param n_overlap Number of overlapping features between two windows.
#' @param knn_method Either "euclidean" (default) or "manhattan". Distance metric for nearest neighbor calculation.
#' @param pca_method "Regularized" by default or "EM".
#' @param .progress Show progress bar (default = `TRUE`).
#'
#' @details
#' The sliding window approach divides the input matrix into smaller, overlapping
#' segments and applies imputation to each window independently. Values in overlapping
#' areas are averaged across windows to produce the final imputed result.
#' This approach assumes that features (columns) are sorted meaningfully (e.g.,
#' by genomic position, time, etc.).
#'
#' Specify `k` and related arguments to use K-NN, `ncp` and related arguments for PCA.
#'
#' @examples
#' # Generate sample data with missing values with 20 samples and 100 columns
#' # where the column order is sorted (e.g., by genomic position or time)
#'
#' set.seed(1234)
#' beta_matrix <- t(sim_mat(100, 20)$input)
#'
#' # Sliding Window K-NN imputation by specifying `k`
#' imputed_knn <- slide_imp(
#'   beta_matrix,
#'   k = 5,
#'   n_feat = 50,
#'   n_overlap = 10
#' )
#' imputed_knn
#'
#' # Sliding Window PCA imputation by specifying `ncp`
#' pca_knn <- slide_imp(
#'   beta_matrix,
#'   ncp = 2,
#'   n_feat = 50,
#'   n_overlap = 10
#' )
#' pca_knn
#'
#' @export
slide_imp <- function(
  obj,
  n_feat,
  n_overlap,
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
  pca_method = c("regularized", "em"),
  coeff.ridge = 1,
  seed = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  # Others
  .progress = TRUE
) {
  if (sum(c(is.null(k), is.null(ncp))) != 1L) {
    stop("Specify either 'k' for K-NN imputation or 'ncp' for PCA imputation. Not both nor neither.")
  }
  imp_method <- if (!is.null(k)) {
    "knn"
  } else {
    "pca"
  }
  # Pre-conditioning ----
  checkmate::assert_matrix(obj, mode = "numeric", row.names = "named", col.names = "unique", null.ok = FALSE, .var.name = "obj")
  cn <- colnames(obj)
  checkmate::assert_int(n_feat, lower = 2, upper = ncol(obj), null.ok = FALSE, .var.name = "n_feat")
  checkmate::assert_int(n_overlap, lower = 0, upper = n_feat - 1, null.ok = FALSE, .var.name = "n_overlap")
  if (imp_method == "knn") {
    knn_method <- match.arg(knn_method)
    checkmate::assert_int(k, lower = 1, upper = n_feat - 1, null.ok = FALSE, .var.name = "k")
    checkmate::assert_number(colmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "colmax")
    checkmate::assert_int(cores, lower = 1, null.ok = FALSE, .var.name = "cores")
    checkmate::assert_flag(post_imp, .var.name = "post_imp", null.ok = FALSE)
    checkmate::assert_number(dist_pow, lower = 0, null.ok = FALSE, .var.name = "dist_pow")
    checkmate::assert(
      checkmate::check_character(subset, min.len = 1, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
      checkmate::check_integerish(subset, lower = 1, upper = ncol(obj), min.len = 1, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
      combine = "or",
      .var.name = "subset"
    )
  } else if (imp_method == "pca") {
    pca_method <- match.arg(pca_method)
    checkmate::assert_flag(scale, .var.name = "scale")
    checkmate::assert_int(ncp, lower = 1, upper = min(n_feat, nrow(obj)), .var.name = "ncp")
    checkmate::assert_number(coeff.ridge, .var.name = "coeff.ridge")
    checkmate::assert_number(seed, null.ok = TRUE, .var.name = "seed")
    # checkmate::assert_numeric(row.w, lower = 0, upper = 1, any.missing = FALSE, len = nrow(obj), null.ok = TRUE, .var.name = "row.w")
    # checkmate::assert_integerish(ind.sup, lower = 1, upper = nrow(obj), any.missing = FALSE, len = nrow(obj), null.ok = TRUE, .var.name = "ind.sup")
    checkmate::assert_int(nb.init, lower = 1, .var.name = "nb.init")
    checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
    checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")
    obj_vars <- col_vars(obj)
    if (
      any(abs(obj_vars) < .Machine$double.eps^0.5 | is.nan(obj_vars) | is.na(obj_vars))
    ) {
      stop("Features with zero variance after na.rm not permitted for PCA Imputation. Try `col_vars(obj)`")
    }
    rm(obj_vars)
  }
  checkmate::assert_flag(.progress, .var.name = ".progress", null.ok = FALSE)

  # Windowing Logic ----
  idx <- 1
  max_step <- ceiling((ncol(obj) - idx) / (n_feat - n_overlap))
  step <- 0:max_step
  start <- idx + (step * n_feat) - (step * n_overlap)
  end <- start + n_feat - 1
  # Overshoot
  n_overshoot <- sum(end > ncol(obj))
  corrected_length <- length(end) - n_overshoot
  start <- start[1:corrected_length]
  end <- end[1:corrected_length]
  end[corrected_length] <- ncol(obj)
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
  result <- matrix(0, nrow = nrow(obj), ncol = ncol(obj), dimnames = list(row.names(obj), colnames(obj)))
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
        miniter = miniter
        # row.w = row.w,
        # ind.sup = ind.sup
      )
    }

    result[, window_cols] <- result[, window_cols] + imputed_window
  }

  if (.progress) {
    message("Step 2/2: Averaging overlapping regions")
  }
  result <- sweep(result, 2, overlap$counts_vec, "/")

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
