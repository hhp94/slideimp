# given a start and end vectors. Give the counts_vec that counts the number
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

#' Sliding-Window K-NN or PCA Imputation
#'
#' Perform sliding-window K-NN or PCA imputation on a numeric matrix whose
#' columns are meaningfully ordered.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#'
#' @param location A sorted numeric vector of length `ncol(obj)` giving the
#'   position of each column, such as genomic coordinates.
#' @param window_size Numeric. Window width in the same units as `location`.
#' @param overlap_size Numeric. Overlap between consecutive windows in the
#'   same units as `location`. Must be less than `window_size`. Ignored when
#'   `flank = TRUE`.
#' @param flank Logical. If `FALSE`, imputation uses sliding windows across the
#'   full matrix. If `TRUE`, one window of width `window_size` is created for
#'   each feature listed in `subset`; `overlap_size` is ignored. `subset` must
#'   be supplied when `flank = TRUE`.
#' @param min_window_n Integer. Minimum number of columns a window must contain
#'   to be considered for imputation. For non-dry runs, the selected `k` or
#'   `ncp` value must be feasible for the usable columns in each retained
#'   window after applying `colmax`.
#' @param subset Optional character or integer vector specifying columns to
#'   impute. If `NULL`, all eligible columns are imputed. Required when
#'   `flank = TRUE`.
#' @param dry_run Logical. If `TRUE`, skip imputation and return a
#'   `slideimp_tbl` describing the windows that would be used after all
#'   filtering rules are applied. In this mode, `k` and `ncp` are not required.
#' @param method Character or `NULL`. For K-NN imputation, one of
#'   `"euclidean"` or `"manhattan"`. For PCA imputation, one of
#'   `"regularized"` or `"EM"`. If `NULL`, the corresponding backend default
#'   is used.
#' @param cores Integer. Number of cores to use for K-NN imputation.
#' @param on_infeasible Character. One of `"skip"`, `"error"`, or `"mean"`.
#'   Controls behavior when a window is infeasible for imputation, for example
#'   when `k` or `ncp` exceeds the number of usable columns after applying
#'   `colmax`.
#'
#' @details
#' The sliding-window approach divides the input matrix into smaller segments
#' based on `location` values and applies imputation to each window
#' independently. Values in overlapping regions are averaged across windows to
#' produce the final imputed result.
#'
#' Two windowing modes are supported:
#'
#' * `flank = FALSE`: greedily partition `location` into windows of width
#'   `window_size` with the requested `overlap_size` between consecutive
#'   windows.
#' * `flank = TRUE`: create one window per feature in `subset`, centered on
#'   that feature using the supplied `window_size`.
#'
#' Specify `k` and related arguments to use [knn_imp()], or `ncp` and related
#' arguments to use [pca_imp()].
#'
#' @inheritSection pca_imp Performance tips
#'
#' @returns If `dry_run = FALSE`, a numeric matrix of the same dimensions as
#'   `obj`, with missing values imputed. The returned object has class
#'   `slideimp_results`.
#'
#'   If `dry_run = TRUE`, a data frame of class `slideimp_tbl` with columns
#'   `start`, `end`, and `window_n`, plus `subset_local` and, when
#'   `flank = TRUE`, `target`.
#'
#' @examples
#' set.seed(1234)
#'
#' # Example data with 20 samples and 100 ordered columns
#' beta_matrix <- sim_mat(20, 100)$input
#' location <- 1:100
#'
#' # First perform a dry run to inspect the calculated windows
#' window_statistics <- slide_imp(
#'   beta_matrix,
#'   location = location,
#'   window_size = 50,
#'   overlap_size = 10,
#'   min_window_n = 10,
#'   dry_run = TRUE,
#'   .progress = FALSE
#' )
#' window_statistics
#'
#' # Sliding-window K-NN imputation
#' imputed_knn <- slide_imp(
#'   beta_matrix,
#'   location = location,
#'   k = 5,
#'   window_size = 50,
#'   overlap_size = 10,
#'   min_window_n = 10,
#'   .progress = FALSE
#' )
#' imputed_knn
#'
#' # Sliding-window PCA imputation
#' imputed_pca <- slide_imp(
#'   beta_matrix,
#'   location = location,
#'   ncp = 2,
#'   window_size = 50,
#'   overlap_size = 10,
#'   min_window_n = 10,
#'   .progress = FALSE
#' )
#' imputed_pca
#'
#' # K-NN imputation with flanking windows
#' imputed_flank <- slide_imp(
#'   beta_matrix,
#'   location = location,
#'   k = 2,
#'   window_size = 30,
#'   flank = TRUE,
#'   subset = c(10, 30, 70),
#'   min_window_n = 5,
#'   .progress = FALSE
#' )
#' imputed_flank
#'
#' @export
slide_imp <- function(
  obj,
  location,
  window_size,
  overlap_size = 0,
  flank = FALSE,
  min_window_n,
  subset = NULL,
  dry_run = FALSE,
  # K-NN-specific parameters
  k = NULL,
  cores = 1,
  dist_pow = 0,
  # PCA-specific parameters
  ncp = NULL,
  scale = TRUE,
  coeff.ridge = 1,
  threshold = 1e-6,
  seed = NULL,
  row.w = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  solver = c("auto", "exact", "lobpcg"),
  lobpcg_control = NULL,
  clamp = NULL,
  # Shared
  method = NULL,
  .progress = TRUE,
  colmax = 0.9,
  post_imp = TRUE,
  na_check = TRUE,
  on_infeasible = c("skip", "error", "mean")
) {
  checkmate::assert_flag(dry_run, .var.name = "dry_run", null.ok = FALSE)
  # minimal pre-conditioning to avoid code fragility
  if (!dry_run) {
    if (sum(c(is.null(k), is.null(ncp))) != 1L) {
      stop("Specify either 'k' for K-NN imputation or 'ncp' for PCA imputation. Not both nor neither.")
    }
    imp_method <- if (!is.null(k)) "knn" else "pca"
  } else {
    # k/ncp irrelevant when only computing window statistics
    imp_method <- NA_character_
  }
  if (!dry_run) {
    on_infeasible <- match.arg(on_infeasible)
  }
  # Pre-conditioning ----
  checkmate::assert_matrix(obj, mode = "numeric", null.ok = FALSE, .var.name = "obj")
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
  checkmate::assert_flag(flank, .var.name = "flank", null.ok = FALSE)
  if (flank) {
    stopifnot("`subset` must be provided when `flank = TRUE`." = !is.null(subset))
    overlap_size <- 0
  }
  checkmate::assert_int(min_window_n, lower = 2L, null.ok = FALSE, .var.name = "min_window_n")

  if (!dry_run && imp_method == "knn") {
    method <- if (is.null(method)) "euclidean" else match.arg(method, c("euclidean", "manhattan"))
    checkmate::assert_int(k, lower = 1L, upper = min_window_n - 1L, null.ok = FALSE, .var.name = "k")
  } else if (!dry_run && imp_method == "pca") {
    method <- if (is.null(method)) "regularized" else match.arg(method, c("regularized", "EM"))
    solver <- match.arg(solver)
    checkmate::assert_int(ncp,
      lower = 1, upper = min(min_window_n - 1L, min(nrow(obj), ncol(obj)) - 1L),
      .var.name = "ncp"
    )
    # other PCA arguments, including `lobpcg_control`, are checked in pca_imp().
    # do not resolve `lobpcg_control` here. The first successful PCA window needs
    # the actual window size so pca_imp can apply solver = "auto". Its chosen
    # concrete solver is then reused for later windows.
  }

  checkmate::assert_flag(.progress, .var.name = ".progress", null.ok = FALSE)
  checkmate::assert_flag(na_check, .var.name = "na_check")

  # resolve subset to sorted integer indices (needed before windowing when flank = TRUE)
  subset <- resolve_subset(subset, obj, sort = TRUE)
  if (is.null(subset)) {
    return(obj)
  }

  # windowing Logic ----
  if (flank) {
    windows <- find_windows_flank(location, subset, window_size)
    start <- windows$start
    end <- windows$end
    # each window has exactly one target column; subset_local gives its local index
    subset_list <- as.list(windows$subset_local)
  } else {
    windows <- find_windows(location, window_size, overlap_size)
    start <- windows$start
    end <- windows$end
  }

  # drop windows that are too small
  window_n <- end - start + 1L
  keep <- window_n >= min_window_n
  if (!any(keep)) {
    stop(
      "All windows have fewer than `min_window_n` (", min_window_n,
      ") columns. Consider increasing `window_size` or decreasing `min_window_n`."
    )
  }
  if (any(!keep) && .progress) {
    message(
      sprintf("Dropping %d window(s) with fewer than %d columns.", sum(!keep), min_window_n)
    )
  }
  start <- start[keep]
  end <- end[keep]

  if (flank) {
    subset_list <- subset_list[keep]
    # track which original target columns (global indices) survived filtering
    target_cols <- subset[keep]
  }

  if (!flank) {
    # build per-window local subset indices. Used by knn_imp; for pca it only
    # drives the "drop windows that cover no subset columns" filter below so
    # that non-flank behavior matches flank mode.
    subset_list <- lapply(seq_along(start), function(i) {
      first <- findInterval(start[i] - 1, subset) + 1
      last <- findInterval(end[i], subset)
      if (first <= last) {
        subset[first:last] - start[i] + 1
      } else {
        integer(0)
      }
    })

    # drop windows that cover no subset columns (parity with flank mode).
    keep_sub <- lengths(subset_list) > 0L
    if (!any(keep_sub)) {
      stop("No windows cover any column in `subset`. Consider increasing `window_size`.")
    }
    if (any(!keep_sub) && .progress) {
      message(
        sprintf("Dropping %d window(s) covering no `subset` columns.", sum(!keep_sub))
      )
    }
    start <- start[keep_sub]
    end <- end[keep_sub]
    subset_list <- subset_list[keep_sub]

    # overlap regions to average over (computed on surviving windows)
    overlap <- find_overlap_regions(start, end)
  }

  # check for all NA/inf column
  for (i in seq_along(start)) {
    window_cols <- start[i]:end[i]
    check_finite(obj[, window_cols, drop = FALSE])
  }

  # early return: window statistics only ----
  if (dry_run) {
    out <- data.frame(
      start = start,
      end = end,
      window_n = end - start + 1L
    )
    if (flank) {
      out$target <- target_cols
    }
    out$subset_local <- subset_list
    class(out) <- c("slideimp_tbl", "data.frame")
    return(out)
  }

  # Sliding Imputation ----
  result <- obj

  copy_max_bytes <- 500 * 1024^2
  copy_max_cols_per_chunk <- as.integer(copy_max_bytes %/% (max(1L, nrow(result)) * 8))
  copy_col_chunk <- min(10000L, max(1L, copy_max_cols_per_chunk))

  if (length(subset) > 0L) {
    for (s in seq(1L, length(subset), by = copy_col_chunk)) {
      e <- min(s + copy_col_chunk - 1L, length(subset))
      result[, subset[s:e]] <- 0
    }
  }

  if (.progress) {
    message("Step 1/2: Imputing")
    n_windows <- length(start)
    n_steps <- max(1, round(n_windows / 20))
  }

  # meta data
  fallback_flags <- logical(length(start))
  skipped_flags <- logical(length(start))

  # for PCA with solver = "auto", let the first successful PCA window choose
  # the concrete solver, then reuse that solver for the remaining windows.
  pca_solver_current <- if (imp_method == "pca") solver else NULL
  pca_solver_lock_window <- NA
  for (i in seq_along(start)) {
    if (.progress && (i %% n_steps == 0 || i == n_windows || i == 1)) {
      message(sprintf(" Processing window %d of %d", i, n_windows))
    }
    window_cols <- start[i]:end[i]
    sub_mat <- obj[, window_cols, drop = FALSE]

    imputed_window <- tryCatch(
      suppressMessages(
        if (imp_method == "knn") {
          knn_imp(
            obj = sub_mat, k = k, colmax = colmax, cores = cores,
            method = method, post_imp = post_imp, dist_pow = dist_pow,
            na_check = FALSE, .progress = FALSE,
            subset = subset_list[[i]]
          )
        } else {
          pca_imp(
            obj = sub_mat, ncp = ncp, scale = scale, method = method,
            coeff.ridge = coeff.ridge, threshold = threshold, seed = seed,
            nb.init = nb.init, maxiter = maxiter, miniter = miniter,
            row.w = row.w, lobpcg_control = lobpcg_control, solver = pca_solver_current,
            na_check = FALSE, colmax = colmax, post_imp = post_imp,
            clamp = clamp
          )
        }
      ),
      slideimp_infeasible = function(e) {
        switch(on_infeasible,
          error = stop(e),
          skip = structure(sub_mat, fallback = TRUE, skipped = TRUE),
          mean = structure(
            mean_imp_col(sub_mat, subset = subset_list[[i]], cores = cores),
            fallback = TRUE
          )
        )
      }
    )

    fallback_flags[i] <- isTRUE(attr(imputed_window, "fallback"))
    skipped_flags[i] <- isTRUE(attr(imputed_window, "skipped"))

    if (
      imp_method == "pca" && solver == "auto" && pca_solver_current == "auto" &&
        !fallback_flags[i]
    ) {
      chosen_solver <- attr(imputed_window, "solver_chosen", exact = TRUE)
      if (
        is.character(chosen_solver) && length(chosen_solver) == 1L && chosen_solver %in% c("exact", "lobpcg")
      ) {
        pca_solver_current <- chosen_solver
        pca_solver_lock_window <- i
      }
    }

    if (skipped_flags[i]) {
      next
    }

    local_idx <- subset_list[[i]]

    if (length(local_idx) > 0L) {
      global_idx <- window_cols[local_idx]
      result[, global_idx] <- result[, global_idx, drop = FALSE] + imputed_window[, local_idx, drop = FALSE]
    }
  }

  if (.progress) {
    message("Step 2/2: Averaging overlapping regions")
  }

  if (flank) {
    # in flank mode each target column is imputed exactly once, no averaging needed.
    # Columns not targeted by any surviving window, or skipped, get original values.
    skipped_targets <- target_cols[skipped_flags]
    uncovered <- union(setdiff(subset, target_cols), skipped_targets)

    if (length(uncovered) > 0) {
      for (s in seq(1L, length(uncovered), by = copy_col_chunk)) {
        e <- min(s + copy_col_chunk - 1L, length(uncovered))
        cols <- uncovered[s:e]
        result[, cols] <- obj[, cols, drop = FALSE]
      }

      if (.progress) {
        message(sprintf("Note: %d column(s) not covered by any window; original values retained.", length(uncovered)))
      }
    }
  } else {
    counts_vec <- overlap$counts_vec
    if (length(counts_vec) < ncol(obj)) {
      counts_vec <- c(counts_vec, rep(0L, ncol(obj) - length(counts_vec)))
    }
    # remove skipped-window contributions
    if (any(skipped_flags)) {
      for (i in which(skipped_flags)) {
        cols_i <- start[i]:end[i]
        counts_vec[cols_i] <- counts_vec[cols_i] - 1L
      }
    }
    # average only target columns where multiple non-skipped windows overlap.
    target_counts <- counts_vec[subset]

    gt_1_idx <- subset[target_counts > 1L]
    if (length(gt_1_idx) > 0) {
      result[, gt_1_idx] <- sweep(result[, gt_1_idx, drop = FALSE], 2, counts_vec[gt_1_idx], "/")
    }

    # restore original values for target columns not covered by any non-skipped window.
    uncovered <- subset[target_counts == 0L]
    if (length(uncovered) > 0) {
      for (s in seq(1L, length(uncovered), by = copy_col_chunk)) {
        e <- min(s + copy_col_chunk - 1L, length(uncovered))
        cols <- uncovered[s:e]
        result[, cols] <- obj[, cols, drop = FALSE]
      }
      if (.progress) {
        message(sprintf("Note: %d column(s) not covered by any window; original values retained.", length(uncovered)))
      }
    }
  }

  if (na_check) {
    if (flank) {
      # all requested targets, not just non-skipped ones
      imputed_cols <- subset
    } else {
      imputed_cols <- subset
    }
    has_remaining_na <- FALSE
    if (length(imputed_cols) > 0L) {
      for (s in seq(1L, length(imputed_cols), by = copy_col_chunk)) {
        e <- min(s + copy_col_chunk - 1L, length(imputed_cols))
        if (anyNA(result[, imputed_cols[s:e], drop = FALSE])) {
          has_remaining_na <- TRUE
          break
        }
      }
    }
  } else {
    has_remaining_na <- NULL
  }
  fallback_windows <- which(fallback_flags)

  class(result) <- c("slideimp_results", class(result))
  attr(result, "imp_method") <- imp_method
  attr(result, "metacaller") <- "slide_imp"
  attr(result, "fallback") <- fallback_windows
  attr(result, "fallback_action") <- on_infeasible # "skip" or "mean"
  attr(result, "has_remaining_na") <- has_remaining_na
  attr(result, "flank") <- flank
  attr(result, "post_imp") <- post_imp
  if (imp_method == "pca") {
    attr(result, "solver_requested") <- solver
    attr(result, "solver_chosen") <- if (solver == "auto") {
      if (is.na(pca_solver_lock_window)) NA else pca_solver_current
    } else {
      solver
    }
    attr(result, "solver_lock_window") <- pca_solver_lock_window
  }

  result
}
