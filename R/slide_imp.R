# Number of windows that cover each of the `n` columns. `start` and `end` are
# 1-based inclusive column indices, one pair per window.
#
# Windows can end on the same column, so the ends are counted with tabulate():
# an indexed update such as `delta[end + 1] <- delta[end + 1] - 1` applies a
# repeated index once.
count_window_cover <- function(start, end, n) {
  bins <- n + 1L
  cumsum(tabulate(start, bins) - tabulate(end + 1L, bins))[seq_len(n)]
}

# A window's value at each hole of its targets, column by column and by row
# within a column, NA where the window left the hole. `res` is what
# knn_imp_parts() or pca_imp_parts() returned for the window with its
# targets as `subset`, or the mean fallback in the same form, and `holes`
# is.na() of the window's targets. fill_imp_parts() writes the parts into a
# copy of the targets, validating them in the window, so the copy is the
# size of the targets.
window_hole_values <- function(res, holes) {
  fill_imp_parts(res)[holes]
}

# The row of each position `at` in a column-major matrix of `n_row` rows, as
# integer, which a part's map must be. which() returns double positions for
# a matrix past 2^31 - 1 cells, and the row taken from them would stay double.
hole_row_index <- function(at, n_row) {
  as.integer((at - 1) %% n_row + 1)
}

# The parts of a slide_imp() result, from the windows' values summed at the
# holes of the `subset` columns. `acc` holds the sums, one per hole, column by
# column in `subset` order and by row within a column; `hole_n` the holes of
# each `subset` column, and `n_ran` the number of windows that ran with it as
# a target. `seen` marks the `subset` columns some window held as a target,
# and `rows` the rows of their holes, in the order of `acc`: the windows
# read them, so no second pass over the input is needed. A hole's value is
# its sum divided by its column's `n_ran`, which is exact for a column one
# window ran over. A column no window ran over has no part, so it comes back
# as it is. A sum holds NA where some window left the hole, and the hole
# stays NA: such columns are one `left_na` part, whose map names the holes
# every window filled, and the other columns one part that fills every
# hole. `kind` is the parts' kind. The first check is on the arguments:
# slide_imp() counts `n_ran` only at targets it has already marked in
# `seen`, so from there it holds by construction.
slide_parts <- function(subset, hole_n, acc, n_ran, seen, rows, kind) {
  if (any(n_ran > 0L & !seen)) {
    stop(
      "slideimp internal, slide_imp: a window ran over a column no window held",
      call. = FALSE
    )
  }
  ran <- n_ran > 0L & hole_n > 0L
  cols <- subset[ran]
  if (length(cols) == 0L) {
    return(list())
  }
  ran_hole <- rep.int(n_ran > 0L, hole_n)
  values <- if (all(ran_hole)) acc else acc[ran_hole]
  values <- values / rep.int(n_ran[ran], hole_n[ran])
  ran_seen <- rep.int(n_ran[seen] > 0L, hole_n[seen])
  map <- cbind(
    row = if (all(ran_seen)) rows else rows[ran_seen],
    col = rep.int(cols, hole_n[ran])
  )
  rm(ran_hole, ran_seen)
  left <- is.na(values)
  if (!any(left)) {
    return(list(new_part(cols = cols, map = map, values = values, kind = kind)))
  }
  left_cols <- unique(map[left, "col"])
  in_left <- map[, "col"] %in% left_cols
  list(
    new_part(
      cols = setdiff(cols, left_cols),
      map = map[!in_left, , drop = FALSE],
      values = values[!in_left],
      kind = kind
    ),
    new_part(
      cols = left_cols,
      map = map[in_left & !left, , drop = FALSE],
      values = values[in_left & !left],
      kind = kind,
      left_na = TRUE
    )
  )
}

# The `dry_run = TRUE` return value. One row per window; `target` is present in
# flank mode only.
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
# attributes only slide_imp() carries.
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
#' Impute missing values in a numeric matrix whose columns are ordered by
#' position, using K-NN or PCA imputation within sliding windows. For Illumina
#' DNA methylation microarrays, use [group_imp()] instead.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#'
#' @param location Numeric. The position of each column of `obj`, such as a
#'   genomic coordinate, sorted in ascending order.
#' @param window_size Numeric. The width of each window, in the units of
#'   `location`. With `flank = TRUE`, the distance each window reaches on each
#'   side of its feature.
#' @param overlap_size Numeric. The overlap between consecutive windows, in the
#'   units of `location`, less than `window_size`. Ignored when `flank = TRUE`.
#' @param flank Logical. If `FALSE`, slide windows across all columns. If
#'   `TRUE`, build one window around each column in `subset`, as described
#'   under Details.
#' @param min_window_n Integer. The smallest number of columns, at least `2`,
#'   that a window needs to be kept. `k` and `ncp` must be less than
#'   `min_window_n`.
#' @param k Integer or `NULL`. Number of nearest neighbors for K-NN
#'   imputation.
#' @param ncp Integer or `NULL`. Number of components for PCA imputation.
#' @param subset Optional character or integer vector. The columns to impute
#'   and return, by name or by number, or `NULL` for every column. Must be
#'   supplied when `flank = TRUE`. Defaults to `NULL`.
#' @param dry_run Logical. If `TRUE`, return the windows that `min_window_n`
#'   and `subset` leave, as a `slideimp_tbl`, without imputing. `k` and `ncp`
#'   are not needed then.
#' @param method Character or `NULL`. One of `"euclidean"` or `"manhattan"`
#'   for K-NN imputation, or `"regularized"` or `"EM"` for PCA imputation.
#'   `NULL` uses `"euclidean"` for K-NN and `"regularized"` for PCA.
#' @param on_infeasible Character. One of `"skip"`, `"error"`, or `"mean"`.
#'   The action for a window that is infeasible for imputation, for example
#'   because `k` or `ncp` is too large for the columns that `colmax` leaves.
#'   Defaults to `"skip"`.
#'
#' @details
#' `slide_imp()` splits the columns of `obj` into windows by `location` and
#' imputes each window on its own. There are two windowing modes:
#'
#' * `flank = FALSE`: greedily partition `location` into windows of width
#'   `window_size` with the requested `overlap_size` between consecutive
#'   windows. A window that would add no column to the window before it is not
#'   created. A column in more than one window takes the mean of their imputed
#'   values.
#' * `flank = TRUE`: create one window per feature in `subset`. Each window
#'   covers the `location` values within `window_size` of that feature, on
#'   both sides, so it spans twice `window_size`.
#'
#' Supply exactly one of `k` and `ncp`. `k` imputes each window with
#' [knn_imp()], and `ncp` imputes it with [pca_imp()]. The arguments of the
#' other method are ignored.
#'
#' The result holds the `subset` columns in the order they have in `obj`.
#' In each window, every column the method can use takes part in the
#' imputation, in `subset` or not: for K-NN every column within `colmax`, and
#' for PCA every column within `colmax` whose variance is not zero. A name
#' that is not in `colnames(obj)`, or a `subset` that selects no column, is an
#' error.
#'
#' `on_infeasible = "skip"` does not impute an infeasible window, `"mean"`
#' fills it with column means, and `"error"` stops at the first such window.
#' With `"skip"` and `"mean"`, `slide_imp()` raises one warning of class
#' `slideimp_fallback` that lists the affected windows.
#'
#' An infeasible window adds nothing to the `mean_filled_cols` and
#' `mean_filled_cells` attributes. Where windows overlap, `mean_filled_cells`
#' lists each cell that at least one window filled with its column mean.
#'
#' With `solver = "auto"`, the first window that completes the solver probe
#' picks the solver, and every later window uses it.
#'
#' `slide_imp()` rejects a column with no observed values only when a window
#' covers it. A column that no window covers is returned unchanged.
#'
#' @inheritSection slideimp-package Column mean fill
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @inheritSection pca_imp PCA Performance tips
#'
#' @returns If `dry_run = FALSE`, a numeric matrix of class `slideimp_results`
#'   with the rows of `obj` and the columns that `subset` selects, with
#'   missing values imputed.
#'
#'   If `dry_run = TRUE`, a data frame of class `slideimp_tbl` with one row per
#'   window and the columns `start` and `end` (the first and last column of
#'   `obj` in the window), `window_n`, `subset_local` and, when
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
  # k/ncp are irrelevant to a dry run.
  imp_method <- NA_character_
  if (!dry_run) {
    if (sum(c(is.null(k), is.null(ncp))) != 1L) {
      cli::cli_abort(
        "Supply exactly one of {.arg k} (K-NN imputation) or {.arg ncp} (PCA imputation)."
      )
    }
    imp_method <- if (!is.null(k)) "knn" else "pca"
    on_infeasible <- match.arg(on_infeasible)
  }
  # Pre-conditioning ----
  checkmate::assert_matrix(
    obj,
    mode = "double",
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
        "Supply {.arg subset} when {.code flank = TRUE}."
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
    # each window's pca_imp() caps `ncp` at `min(nrow(obj) - 2, n_elig - 1)`
    # and raises `slideimp_infeasible` above it. The row half is enforced
    # here, or every window is infeasible and "skip" returns the input
    # unimputed. The column half depends on `colmax` and per-window variance,
    # so `min_window_n - 1L` is the static bound for it.
    max_ncp <- min(min_window_n - 1L, ncol(obj) - 1L, nrow(obj) - 2L)
    if (max_ncp < 1L) {
      cli::cli_abort(
        c(
          "PCA imputation is infeasible: {.arg obj} needs at least 3 rows and 2 columns.",
          "x" = "{.arg obj} has {nrow(obj)} row{?s} and {ncol(obj)} column{?s}."
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
    # other PCA arguments are checked in pca_imp(). `lobpcg_control` is not
    # resolved here: solver = "auto" needs each window's size.
  }

  checkmate::assert_flag(.progress, .var.name = ".progress", null.ok = FALSE)

  # the columns the result holds, ascending: `result_cols` is NULL for every
  # column, and `subset` then lists every column, for the windows. Needed
  # before windowing when flank = TRUE. resolve_narrow_subset() refuses an
  # unknown name and an empty `subset`, in this frame's name.
  result_cols <- resolve_narrow_subset(subset, obj)
  subset <- if (is.null(result_cols)) seq_len(ncol(obj)) else result_cols

  # windowing Logic ----
  windows <- if (flank) {
    find_windows_flank(location, subset, window_size)
  } else {
    find_windows(location, window_size, overlap_size)
  }

  # drop windows that are too small
  keep <- windows$end - windows$start + 1L >= min_window_n
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
  start <- windows$start[keep]
  end <- windows$end[keep]

  # `target_j[i]` is the position in `subset` of window i's first target: a
  # window's targets are consecutive in `subset`, so its holes are one run of
  # `acc` below. findInterval() runs once over all windows, not per window: it
  # costs a pass over `subset` per call.
  if (flank) {
    # each window has exactly one target column; subset_local gives its local index
    subset_list <- as.list(windows$subset_local[keep])
    # track which original target columns (global indices) survived filtering
    target_cols <- subset[keep]
    target_j <- which(keep)
  } else {
    # build per-window local subset indices, ascending: the `subset` each
    # window passes to the backend. They also drive the "drop windows that
    # cover no subset columns" filter below, so that non-flank behavior
    # matches flank mode.
    first <- findInterval(start - 1, subset) + 1
    last <- findInterval(end, subset)
    subset_list <- lapply(seq_along(start), function(i) {
      if (first[i] <= last[i]) {
        subset[first[i]:last[i]] - start[i] + 1
      } else {
        integer(0)
      }
    })

    # drop windows that cover no subset columns (parity with flank mode).
    keep_sub <- lengths(subset_list) > 0L
    if (!any(keep_sub)) {
      cli::cli_abort(c(
        "No windows cover any column in {.arg subset}.",
        "i" = "Increase {.arg window_size} or decrease {.arg min_window_n}."
      ))
    }
    if (any(!keep_sub) && .progress) {
      n_drop_sub <- sum(!keep_sub)
      cli::cli_inform(
        "Dropping {n_drop_sub} window{?s} with no column in {.arg subset}."
      )
    }
    start <- start[keep_sub]
    end <- end[keep_sub]
    subset_list <- subset_list[keep_sub]
    target_j <- as.integer(first[keep_sub])
    rm(first, last)
  }

  # One scan of the whole matrix, before any window runs, so a dry run
  # refuses what a full run does. Each window's scan is this one indexed
  # (scan_cols()), and the variances serve only PCA windows. An Inf is
  # refused wherever it sits, and a column with no observed value only where a
  # window reads it.
  scan <- col_scan(obj, var = identical(imp_method, "pca"))
  covered_cols <- which(count_window_cover(start, end, ncol(obj)) > 0L)
  refuse_nonfinite(scan, empty_j = covered_cols)

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
  n_windows <- length(start)
  if (.progress) {
    cli::cli_inform("Step 1/2: Imputing")
    n_steps <- max(1, round(n_windows / 20))
  }

  # The windows hand back parts, and no matrix the size of `obj` is built
  # until the one write at the end. `acc` holds one double per hole of the
  # `subset` columns, column by column in `subset` order and by row within a
  # column, and each window that runs adds its value at every hole of its
  # targets, in window order, from 0. A hole the window left adds NA, so it
  # stays NA. An observed cell is never written. `hole_off[j]`
  # counts the holes before `subset[j]`, in double, since on a large input
  # they can pass 2^31 - 1.
  hole_n <- scan$n_miss[subset]
  hole_off <- c(0, cumsum(as.double(hole_n)))
  acc <- numeric(hole_off[[length(hole_off)]])
  # The rows of those holes, kept from the first window that holds a column
  # as a target, so the parts' maps need no second pass over `obj`. The
  # windows run in order of `location`, and each one's targets are
  # consecutive in `subset`, so a window's new targets are those after
  # `n_seen`, the last `subset` column a window held. `hole_rows` holds one
  # piece per window, in `subset` order, and `seen` marks the columns read.
  hole_rows <- vector("list", length(start))
  seen <- logical(length(subset))
  n_seen <- 0L
  # the number of windows that ran with each `subset` column as a target,
  # which its holes' sums are divided by. A skipped window does not count; a
  # mean-filled one does.
  n_ran <- integer(length(subset))

  # meta data
  fallback_flags <- logical(length(start))
  # what `post_imp` set to a column mean inside each window, numbered as the
  # result's columns, which are `subset`. NULL for a window that recorded
  # nothing, which is most.
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
    sub_scan <- scan_cols(scan, window_cols)
    local_idx <- as.integer(subset_list[[i]])
    # The window's targets are consecutive in `subset`, from its `j`-th
    # column: a window takes every `subset` column it covers, and in flank
    # mode its one target. Their holes are therefore one run of `acc`.
    j <- target_j[[i]]
    tj <- j - 1L + seq_along(local_idx)
    if (any(subset[tj] != window_cols[local_idx])) {
      stop(
        "slideimp internal, slide_imp: the targets of window ",
        i,
        " are not consecutive in 'subset'",
        call. = FALSE
      )
    }
    holes <- is.na(sub_mat[, local_idx, drop = FALSE])
    new_t <- tj > n_seen
    if (any(new_t)) {
      # the rows of the holes of the new targets, column by column, from
      # their positions in `holes`
      at_new <- which(holes[, new_t, drop = FALSE])
      hole_rows[[i]] <- hole_row_index(at_new, nrow(holes))
      seen[tj[new_t]] <- TRUE
      n_seen <- tj[[length(tj)]]
    }

    # Both internals take the window's targets as `subset` and return parts
    # in the window's numbering, which together cover the targets. A window
    # whose targets are complete returns no part, through the internal's own
    # early exit. An infeasible window under "skip" or "mean" is NULL here.
    res <- tryCatch(
      suppressMessages(
        if (imp_method == "knn") {
          knn_imp_parts(
            obj = sub_mat,
            scan = sub_scan,
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
          pca_imp_parts(
            obj = sub_mat,
            scan = sub_scan,
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
            subset = local_idx,
            .progress = FALSE
          )
        }
      ),
      slideimp_infeasible = function(e) {
        if (on_infeasible == "error") {
          cli::cli_abort(
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
          )
        }
        NULL
      }
    )

    if (is.null(res)) {
      fallback_flags[i] <- TRUE
      if (on_infeasible == "skip") {
        next
      }
      # on_infeasible = "mean": each hole of the targets set to its column
      # mean, as the internals' own result, so the same door reads it
      res <- list(
        obj = sub_mat,
        parts = list(fallback_part(
          sub_mat,
          local_idx[hole_n[tj] > 0L],
          "mean"
        )),
        keep = local_idx
      )
    }

    n_ran[tj] <- n_ran[tj] + 1L
    n_holes <- hole_off[[j + length(local_idx)]] - hole_off[[j]]
    if (n_holes > 0) {
      at <- hole_off[[j]] + seq_len(n_holes)
      acc[at] <- acc[at] + window_hole_values(res, holes)
    }

    # a mean fallback chooses no solver and carries no record: `fallback`
    # already reports it
    if (fallback_flags[i]) {
      next
    }

    # `pca_solver_current` starts at `solver` and leaves "auto" only here, so
    # it is "auto" only while `solver = "auto"` has not locked
    if (
      imp_method == "pca" &&
        pca_solver_current == "auto" &&
        isTRUE(res$info$solver_probed)
    ) {
      # only a completed probe locks the solver. pca_imp() also reports
      # solver_chosen = "exact" when its size heuristic demoted the window or
      # the EM loop ended before the probe finished; locking on that would
      # skip the per-window demotion and warmup_iters raise of "auto".
      pca_solver_current <- res$info$solver_chosen
      pca_solver_lock_window <- i
    }

    # The internals' result columns are the window's targets, and the
    # window's `t`-th target is the result's column `tj[t]`.
    info <- res$info
    if (
      length(info$mean_filled_cols) > 0L || nrow(info$mean_filled_cells) > 0L
    ) {
      mean_filled[[i]] <- mean_filled_record(info, tj)
    }
  }
  rm(res, sub_mat, sub_scan, holes)
  hole_rows <- unlist(hole_rows, use.names = FALSE)
  if (is.null(hole_rows)) {
    hole_rows <- integer(0)
  }

  if (.progress) {
    cli::cli_inform("Step 2/2: Averaging overlapping windows")
  }

  # every window marks its targets in `seen` before it runs, and a window
  # holds as a target every `subset` column it covers, so the columns no
  # window covers are the ones never seen
  n_no_window <- sum(!seen)

  # one write: the averaged holes into one copy of the columns the result
  # holds. A column no window ran over comes back as it is.
  parts <- slide_parts(subset, hole_n, acc, n_ran, seen, hole_rows, imp_method)
  rm(acc, hole_rows)
  result <- fill_parts(obj, parts, keep = result_cols, n_miss = scan$n_miss)
  rm(parts)

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
    # a pinned `solver`, or the one "auto" locked on; NA while "auto" never
    # locked
    solver_chosen = if (imp_method == "pca" && pca_solver_current != "auto") {
      pca_solver_current
    } else {
      NA_character_
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
