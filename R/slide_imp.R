# Number of windows that cover each of the `n` columns. `start` and `end` are
# 1-based inclusive column indices, one pair per window.
#
# tabulate() counts every occurrence of an index. An indexed update of the form
# `delta[end + 1] <- delta[end + 1] - 1` applies a repeated index once, so two
# windows that end on the same column are counted as one and every later column
# carries a count one too high. Nothing here may assume the ends are distinct.
count_window_cover <- function(start, end, n) {
  bins <- n + 1L
  cumsum(tabulate(start, bins) - tabulate(end + 1L, bins))[seq_len(n)]
}

# The `dry_run = TRUE` return value. One row per window; `target` is present in
# flank mode only. Shared by the empty-subset exit and the main path so that
# both hand back the same shape.
slide_window_tbl <- function(start, end, subset_local, target = NULL) {
  out <- data.frame(
    start = start,
    end = end,
    window_n = end - start + 1L
  )
  if (!is.null(target)) {
    out$target <- target
  }
  out$subset_local <- subset_local
  class(out) <- c("slideimp_tbl", "data.frame")
  out
}

# The `dry_run = FALSE` return value: the shared constructor plus the
# attributes only slide_imp() carries. Shared by the empty-subset exit and the
# main path so that a caller cannot tell from the result which one ran.
#
# Passing the matrix through a function does not copy it: measured at 0 extra
# blocks on a 76 MB matrix, against 1 for a forced copy.
slide_imp_results <- function(
  result,
  imp_method,
  fallback,
  on_infeasible,
  flank,
  post_imp,
  solver = NULL,
  solver_chosen = NA_character_,
  solver_lock_window = NA_integer_,
  mean_filled_cols = integer(0),
  mean_filled_cells = mean_filled_cells_matrix()
) {
  out <- new_slideimp_results(
    result,
    imp_method,
    fallback = fallback,
    post_imp = post_imp,
    mean_filled_cols = mean_filled_cols,
    mean_filled_cells = mean_filled_cells
  )
  attr(out, "metacaller") <- "slide_imp"
  attr(out, "fallback_action") <- on_infeasible
  attr(out, "flank") <- flank
  if (imp_method == "pca") {
    attr(out, "solver_requested") <- solver
    attr(out, "solver_chosen") <- solver_chosen
    attr(out, "solver_lock_window") <- solver_lock_window
  }
  out
}

#' Sliding-Window K-NN or PCA Imputation
#'
#' Perform sliding-window K-NN or PCA imputation on a numeric matrix whose
#' columns are meaningfully ordered. Not intended for Illumina DNA
#' methylation microarrays, use [group_imp()] instead.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#'
#' @param location A sorted numeric vector of length `ncol(obj)` giving the
#'   position of each column, such as genomic coordinates.
#' @param window_size Numeric. Window width in the same units as `location`.
#'   When `flank = TRUE`, the distance covered on each side of a feature.
#' @param overlap_size Numeric. Overlap between consecutive windows in the
#'   same units as `location`. Must be less than `window_size`. Ignored when
#'   `flank = TRUE`.
#' @param flank Logical. If `FALSE`, imputation uses sliding windows across the
#'   full matrix. If `TRUE`, one window is created for each feature listed in
#'   `subset`, and the window extends `window_size` on each side of that
#'   feature. `overlap_size` is ignored. `subset` must be supplied when
#'   `flank = TRUE`.
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
#'   windows. A window that would add no column to the window before it is not
#'   created.
#' * `flank = TRUE`: create one window per feature in `subset`. Each window
#'   covers the `location` values within `window_size` of that feature, on
#'   both sides, so it spans twice `window_size`.
#'
#' Specify `k` and related arguments to use [knn_imp()], or `ncp` and related
#' arguments to use [pca_imp()].
#'
#' `on_infeasible = "skip"` does not impute an infeasible window, `"mean"`
#' fills it with column means, and `"error"` stops at the first such window.
#' With `"skip"` and `"mean"`, `slide_imp()` raises one warning of class
#' `slideimp_fallback` that lists the affected windows.
#'
#' An infeasible window adds nothing to the `mean_filled_cols` and
#' `mean_filled_cells` attributes of the result. Where windows overlap,
#' `mean_filled_cells` lists a cell that any window filled with the column
#' mean.
#'
#' @inheritSection slideimp-package Column mean fill
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @inheritSection pca_imp PCA Performance tips
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
  on_infeasible = c("skip", "error", "mean")
) {
  # conditions raised from inside a handler report this frame, not the handler
  slide_env <- environment()
  checkmate::assert_flag(dry_run, .var.name = "dry_run", null.ok = FALSE)
  # minimal pre-conditioning to avoid code fragility
  if (!dry_run) {
    if (sum(c(is.null(k), is.null(ncp))) != 1L) {
      cli::cli_abort(
        "Specify either {.arg k} for K-NN imputation or {.arg ncp} for PCA imputation. Not both nor neither."
      )
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
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    null.ok = FALSE,
    .var.name = "obj"
  )
  checkmate::assert_numeric(
    location,
    len = ncol(obj),
    any.missing = FALSE,
    sorted = TRUE,
    finite = TRUE,
    null.ok = FALSE,
    .var.name = "location"
  )
  checkmate::assert_number(
    window_size,
    lower = .Machine$double.eps,
    finite = TRUE,
    null.ok = FALSE,
    .var.name = "window_size"
  )
  checkmate::assert_number(
    overlap_size,
    lower = 0,
    finite = TRUE,
    null.ok = FALSE,
    .var.name = "overlap_size"
  )
  checkmate::assert_flag(flank, .var.name = "flank", null.ok = FALSE)
  if (flank) {
    if (is.null(subset)) {
      cli::cli_abort(
        "{.arg subset} must be provided when {.code flank = TRUE}."
      )
    }
    # subset is validated later on. `overlap_size` is documented as ignored in
    # this mode, so it is not compared against `window_size` either.
    overlap_size <- 0
  } else if (overlap_size >= window_size) {
    cli::cli_abort(
      "{.arg overlap_size} must be strictly less than {.arg window_size}."
    )
  }
  checkmate::assert_int(
    min_window_n,
    lower = 2L,
    null.ok = FALSE,
    .var.name = "min_window_n"
  )

  if (!dry_run && imp_method == "knn") {
    method <- if (is.null(method)) {
      "euclidean"
    } else {
      match.arg(method, c("euclidean", "manhattan"))
    }
    checkmate::assert_int(
      k,
      lower = 1L,
      upper = min_window_n - 1L,
      null.ok = FALSE,
      .var.name = "k"
    )
  } else if (!dry_run && imp_method == "pca") {
    method <- if (is.null(method)) {
      "regularized"
    } else {
      match.arg(method, c("regularized", "EM"))
    }
    solver <- match.arg(solver)
    # every window is handed to pca_imp(), which caps `ncp` at
    # `min(nrow(obj) - 2, n_elig - 1)` and raises `slideimp_infeasible` above
    # it. The row half of that cap is knowable here and has to be enforced
    # here: `nrow(obj) - 1L` passes this check and then makes EVERY window
    # infeasible, which `on_infeasible = "skip"` turns into the input returned
    # unimputed. The column half depends on `colmax` and per-window variance,
    # so `min_window_n - 1L` is the tightest static bound available for it.
    max_ncp <- min(min_window_n - 1L, ncol(obj) - 1L, nrow(obj) - 2L)
    if (max_ncp < 1L) {
      cli::cli_abort(
        c(
          "PCA imputation is infeasible with the current data and settings.",
          "i" = "Number of rows: {nrow(obj)}. PCA needs at least 3.",
          "i" = "Number of columns: {ncol(obj)}.",
          "i" = "{.arg min_window_n}: {min_window_n}."
        ),
        class = "slideimp_infeasible"
      )
    }
    checkmate::assert_int(
      ncp,
      lower = 1,
      upper = max_ncp,
      .var.name = "ncp"
    )
    # other PCA arguments, including `lobpcg_control`, are checked in pca_imp().
    # do not resolve `lobpcg_control` here. The first successful PCA window needs
    # the actual window size so pca_imp can apply solver = "auto". Its chosen
    # concrete solver is then reused for later windows.
  }

  checkmate::assert_flag(.progress, .var.name = ".progress", null.ok = FALSE)

  # resolve subset to sorted integer indices (needed before windowing when flank = TRUE)
  subset <- resolve_subset(subset, obj, sort = TRUE)
  if (is.null(subset)) {
    # resolve_subset() has already reported the empty subset. Nothing is
    # imputed, but the return value keeps the shape `@returns` promises: no
    # windows for a dry run, the input as a `slideimp_results` otherwise.
    if (dry_run) {
      return(slide_window_tbl(
        start = integer(0),
        end = integer(0),
        subset_local = list(),
        target = if (flank) integer(0) else NULL
      ))
    }
    return(slide_imp_results(
      obj,
      imp_method,
      fallback = integer(0),
      on_infeasible = on_infeasible,
      flank = flank,
      post_imp = post_imp,
      solver = if (imp_method == "pca") solver else NULL
    ))
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
    cli::cli_abort(c(
      "All windows have fewer than {.arg min_window_n} ({min_window_n}) columns.",
      "i" = "Increase {.arg window_size} or decrease {.arg min_window_n}."
    ))
  }
  if (any(!keep) && .progress) {
    n_drop <- sum(!keep)
    cli::cli_inform(
      "Dropping {n_drop} window{?s} with fewer than {min_window_n} columns."
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
      cli::cli_abort(c(
        "No windows cover any column in {.arg subset}.",
        "i" = "Increase {.arg window_size}."
      ))
    }
    if (any(!keep_sub) && .progress) {
      n_drop_sub <- sum(!keep_sub)
      cli::cli_inform(
        "Dropping {n_drop_sub} window{?s} covering no {.arg subset} columns."
      )
    }
    start <- start[keep_sub]
    end <- end[keep_sub]
    subset_list <- subset_list[keep_sub]

    # how many surviving windows cover each column, to average over
    window_counts <- count_window_cover(start, end, ncol(obj))
  }

  # every bulk column pass below copies at most this many columns at a time
  copy_max_bytes <- 500 * 1024^2
  copy_max_cols_per_chunk <- as.integer(
    copy_max_bytes %/% (max(1L, nrow(obj)) * 8)
  )
  copy_col_chunk <- min(10000L, max(1L, copy_max_cols_per_chunk))

  # validate covered columns once (not a submatrix per window)
  covered_cols <- collapse::funique(
    sequence(end - start + 1L, from = start),
    sort = TRUE
  )

  if (length(covered_cols) == ncol(obj)) {
    check_finite(obj)
  } else {
    # Inf: whole matrix so the reported row/col refer to obj
    check_inf(obj)
    # all-NA: only columns a window will impute. check_inf() has cleared the
    # whole matrix, so the one refusal left is the all-NA column, whose message
    # carries no position. Chunking therefore cannot change what is reported.
    for (s in seq(1L, length(covered_cols), by = copy_col_chunk)) {
      e <- min(s + copy_col_chunk - 1L, length(covered_cols))
      check_finite(obj[, covered_cols[s:e], drop = FALSE])
    }
  }

  # early return: window statistics only ----
  if (dry_run) {
    return(slide_window_tbl(
      start = start,
      end = end,
      subset_local = subset_list,
      target = if (flank) target_cols else NULL
    ))
  }

  # Sliding Imputation ----
  result <- obj

  if (length(subset) > 0L) {
    for (s in seq(1L, length(subset), by = copy_col_chunk)) {
      e <- min(s + copy_col_chunk - 1L, length(subset))
      result[, subset[s:e]] <- 0
    }
  }

  n_windows <- length(start)
  if (.progress) {
    cli::cli_inform("Step 1/2: Imputing")
    n_steps <- max(1, round(n_windows / 20))
  }

  # meta data
  fallback_flags <- logical(length(start))
  skipped_flags <- logical(length(start))
  # what `post_imp` set to a column mean inside each window, in `obj`'s own
  # column indices. NULL for a window that recorded nothing, which is most.
  mean_filled <- vector("list", length(start))

  # for PCA with solver = "auto", let the first successful PCA window choose
  # the concrete solver, then reuse that solver for the remaining windows.
  pca_solver_current <- if (imp_method == "pca") solver else NULL
  pca_solver_lock_window <- NA_integer_
  for (i in seq_along(start)) {
    if (.progress && (i %% n_steps == 0 || i == n_windows || i == 1)) {
      cli::cli_inform("Processing window {i} of {n_windows}")
    }
    window_cols <- start[i]:end[i]
    sub_mat <- obj[, window_cols, drop = FALSE]
    local_idx <- subset_list[[i]]

    # pca_imp() has no `subset`. It imputes the whole window, and it raises
    # `slideimp_infeasible` for a window whose only missing values sit in
    # ineligible columns. When the requested columns of this window are already
    # complete there is nothing to take from it, so it is not called.
    # knn_imp() makes the same exit itself through its `subset` argument.
    targets_complete <- imp_method == "pca" &&
      !anyNA(
        if (length(local_idx) == ncol(sub_mat)) {
          sub_mat
        } else {
          sub_mat[, local_idx, drop = FALSE]
        }
      )

    imputed_window <- if (targets_complete) {
      sub_mat
    } else {
      tryCatch(
        suppressMessages(
          if (imp_method == "knn") {
            knn_imp(
              obj = sub_mat,
              k = k,
              colmax = colmax,
              cores = cores,
              method = method,
              post_imp = post_imp,
              dist_pow = dist_pow,
              .progress = FALSE,
              subset = local_idx
            )
          } else {
            pca_imp(
              obj = sub_mat,
              ncp = ncp,
              scale = scale,
              method = method,
              coeff.ridge = coeff.ridge,
              threshold = threshold,
              seed = seed,
              nb.init = nb.init,
              maxiter = maxiter,
              miniter = miniter,
              row.w = row.w,
              lobpcg_control = lobpcg_control,
              solver = pca_solver_current,
              colmax = colmax,
              post_imp = post_imp,
              clamp = clamp,
              .progress = FALSE
            )
          }
        ),
        slideimp_infeasible = function(e) {
          switch(
            on_infeasible,
            error = cli::cli_abort(
              c(
                if (flank) {
                  "Window {i} (target column {target_cols[i]}, columns {start[i]} to {end[i]}) is infeasible for imputation."
                } else {
                  "Window {i} (columns {start[i]} to {end[i]}) is infeasible for imputation."
                },
                "i" = "Use {.code on_infeasible = \"skip\"} or {.code on_infeasible = \"mean\"} to handle infeasible windows."
              ),
              parent = e,
              class = "slideimp_infeasible",
              call = slide_env
            ),
            skip = structure(sub_mat, fallback = TRUE, skipped = TRUE),
            mean = structure(
              mean_imp_col(sub_mat, subset = local_idx, cores = cores),
              fallback = TRUE
            )
          )
        }
      )
    }

    fallback_flags[i] <- isTRUE(attr(imputed_window, "fallback"))
    skipped_flags[i] <- isTRUE(attr(imputed_window, "skipped"))

    if (
      imp_method == "pca" &&
        solver == "auto" &&
        pca_solver_current == "auto" &&
        !fallback_flags[i]
    ) {
      # only a COMPLETED auto probe is a verdict worth reusing. pca_imp() also
      # reports solver_chosen = "exact" when its own size heuristic demoted the
      # window before the kernel ran (n_gram < 250, or k_eig / n_gram > 0.10),
      # and when the EM loop ended before the probe finished. Pinning on either
      # hands every later window a decision that was never made - and a pinned
      # solver bypasses the per-window size demotion and the warmup_iters raise
      # that solver = "auto" applies. Leaving `pca_solver_current` at "auto"
      # lets the next window decide for itself, which costs nothing measurable:
      # the probe is cheap relative to the window it is deciding about.
      probed <- isTRUE(attr(imputed_window, "solver_probed", exact = TRUE))
      chosen_solver <- attr(imputed_window, "solver_chosen", exact = TRUE)
      if (
        probed &&
          is.character(chosen_solver) &&
          length(chosen_solver) == 1L &&
          chosen_solver %in% c("exact", "lobpcg")
      ) {
        pca_solver_current <- chosen_solver
        pca_solver_lock_window <- i
      }
    }

    if (skipped_flags[i]) {
      next
    }

    if (length(local_idx) > 0L) {
      global_idx <- window_cols[local_idx]
      result[, global_idx] <- result[, global_idx, drop = FALSE] +
        imputed_window[, local_idx, drop = FALSE]
    }

    # a window the backend did not run - complete targets, or a mean fallback,
    # which `fallback` already reports - carries no record and yields nothing
    filled <- local_mean_filled(imputed_window, local_idx, window_cols)
    if (length(filled$cols) > 0L || length(filled$row) > 0L) {
      mean_filled[[i]] <- filled
    }
  }

  if (.progress) {
    cli::cli_inform("Step 2/2: Averaging overlapping regions")
  }

  if (flank) {
    # in flank mode each target column is imputed exactly once, no averaging needed.
    # Columns not targeted by any surviving window, or skipped, get original values.
    no_window <- setdiff(subset, target_cols)
    uncovered <- union(no_window, target_cols[skipped_flags])
    n_no_window <- length(no_window)
  } else {
    n_no_window <- sum(window_counts[subset] == 0L)
    # remove skipped-window contributions
    counts_vec <- window_counts
    for (i in which(skipped_flags)) {
      cols_i <- start[i]:end[i]
      counts_vec[cols_i] <- counts_vec[cols_i] - 1L
    }
    # average only target columns where multiple non-skipped windows overlap.
    target_counts <- counts_vec[subset]

    gt_1_idx <- subset[target_counts > 1L]
    if (length(gt_1_idx) > 0L) {
      for (s in seq(1L, length(gt_1_idx), by = copy_col_chunk)) {
        e <- min(s + copy_col_chunk - 1L, length(gt_1_idx))
        cols <- gt_1_idx[s:e]
        result[, cols] <- sweep(
          result[, cols, drop = FALSE],
          2,
          counts_vec[cols],
          "/"
        )
      }
    }

    # target columns not covered by any non-skipped window.
    uncovered <- subset[target_counts == 0L]
  }

  # restore original values where no window contributed.
  if (length(uncovered) > 0L) {
    for (s in seq(1L, length(uncovered), by = copy_col_chunk)) {
      e <- min(s + copy_col_chunk - 1L, length(uncovered))
      cols <- uncovered[s:e]
      result[, cols] <- obj[, cols, drop = FALSE]
    }
  }
  # skipped windows are reported by the warning below, not here.
  if (.progress && n_no_window > 0L) {
    cli::cli_inform(c(
      "No window covers {n_no_window} requested column{?s}.",
      "i" = "{cli::qty(n_no_window)}{?This column keeps/These columns keep} the original values."
    ))
  }

  fallback_windows <- which(fallback_flags)
  n_fallback <- length(fallback_windows)
  if (n_fallback > 0L) {
    # raised whatever `.progress` says: the default `on_infeasible = "skip"`
    # otherwise hands back unimputed columns with no condition at all.
    cli::cli_warn(
      c(
        "{n_fallback} of {n_windows} window{?s} {cli::qty(n_fallback)}{?is/are} infeasible for imputation.",
        "x" = if (on_infeasible == "skip") {
          "{cli::qty(n_fallback)}Skipped window{?s}: {fmt_trunc(fallback_windows)}. Missing values in a skipped window can remain {.code NA}."
        } else {
          "{cli::qty(n_fallback)}Window{?s} filled with column means: {fmt_trunc(fallback_windows)}."
        },
        "i" = "Use {.code on_infeasible = \"error\"} to see why a window is infeasible."
      ),
      class = "slideimp_fallback"
    )
  }

  mean_filled <- combine_mean_filled(mean_filled, nrow(obj))

  slide_imp_results(
    result,
    imp_method,
    fallback = fallback_windows,
    on_infeasible = on_infeasible, # "skip" or "mean"
    flank = flank,
    post_imp = post_imp,
    solver = if (imp_method == "pca") solver else NULL,
    solver_chosen = if (imp_method != "pca") {
      NA_character_
    } else if (solver != "auto") {
      solver
    } else if (is.na(pca_solver_lock_window)) {
      NA_character_
    } else {
      pca_solver_current
    },
    solver_lock_window = pca_solver_lock_window,
    # a column's missing rate does not depend on the window, so a column held
    # out by `colmax` is held out in every window that ran. A cell the method
    # left can differ between windows: it is recorded when any window filled
    # it with the column mean, even where another window's value is averaged in.
    mean_filled_cols = mean_filled$cols,
    mean_filled_cells = mean_filled$cells
  )
}
