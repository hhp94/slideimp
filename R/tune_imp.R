# for each rep, we are sampling such that for each col we don't:
# * Break row room budget
# * Break col room budget
# * Introduce zcv columns
# only then, after collecting n_cols, do we return a rep. The row room budget
# is more tricky, since we are iterating column wise, we have to update
# row room budget after each selected column. Convert to arma is straightforawrd
sample_each_rep <- function(
  obj, # arma::mat &obj
  pool_idx,
  na_per_col,
  row_room,
  col_room,
  max_attempts
) {
  # arma seed is hooked to R's seed. So we won't have to seed in arma.
  n_cols <- length(na_per_col) # total cols needed
  total_na <- sum(na_per_col) # total na values needed. Needed to track 2D index

  for (iter in seq_len(max_attempts)) {
    # resetting state by reshuffling column order.
    shuffled <- pool_idx[sample.int(length(pool_idx))]
    attempt_row_room <- row_room
    row_out <- integer(total_na)
    col_out <- integer(total_na)
    written <- 0L
    slot <- 0L # each column occupy a "slot"
    for (col_idx in shuffled) {
      if (slot >= n_cols) {
        # if we filled n_cols, then no longer need to scan the available pool
        break
      }
      needed <- na_per_col[slot + 1L] # how many row NA is needed for this col

      # column budget (colmax). This column don't have enough rows left for us given budget
      if (col_room[col_idx] < needed) {
        next
      }

      candidate <- obj[, col_idx] # slice column
      obs_row <- which(!is.na(candidate)) # arma::find_finite(candidate)
      # protect 2 distinct observed values so column stays non-ZV
      uniq_vals <- unique(candidate[obs_row])
      if (length(uniq_vals) < 2L) {
        next
      }
      kept_vals <- sample(uniq_vals, size = 2L)
      # match keep the first two sampled missing values. This brings `kept_global` back to
      # to global row space. But we probably can find_unique() and sample the index in arma.
      kept_global <- obs_row[match(kept_vals, candidate[obs_row])]
      left_over <- setdiff(obs_row, kept_global)
      # row budget (rowmax)
      left_over <- left_over[attempt_row_room[left_over] > 0L]
      if (length(left_over) < needed) {
        next
      }
      # avoid sample() scalar trap
      selected <- left_over[sample.int(length(left_over), needed)]
      # commit
      slot <- slot + 1L
      idx <- (written + 1L):(written + needed)
      row_out[idx] <- selected
      col_out[idx] <- col_idx
      written <- written + needed
      attempt_row_room[selected] <- attempt_row_room[selected] - 1L
    }

    if (slot == n_cols) {
      return(cbind(row = row_out, col = col_out))
    }
  }

  cli::cli_abort("Failed to sample NA locations after {max_attempts} attempts.")
}

#' @export
sample_na_loc_stratified <- function(obj,
                                     n_cols = NULL,
                                     n_rows = 1L,
                                     num_na = NULL,
                                     n_reps = 1L,
                                     rowmax = 0.9,
                                     colmax = 0.9,
                                     subset = NULL,
                                     max_attempts = 10L) {
  checkmate::assert_matrix(obj, min.rows = 1, min.cols = 1, .var.name = "obj")
  checkmate::assert_count(n_cols, positive = TRUE, null.ok = TRUE, .var.name = "n_cols")
  checkmate::assert_count(n_rows, positive = TRUE, .var.name = "n_rows")
  checkmate::assert_count(num_na, positive = TRUE, null.ok = TRUE, .var.name = "num_na")
  checkmate::assert_count(n_reps, positive = TRUE, .var.name = "n_reps")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_count(max_attempts, positive = TRUE, .var.name = "max_attempts")

  if (is.null(n_cols) && is.null(num_na)) {
    cli::cli_abort("Either {.arg n_cols} or {.arg num_na} must be supplied.")
  }
  if (!is.null(num_na) && num_na < n_rows) {
    cli::cli_abort("{.arg num_na} ({num_na}) must be >= {.arg n_rows} ({n_rows}).")
  }

  # resolve (n_cols, na_per_col, max_row_miss)
  if (!is.null(num_na)) {
    n_cols <- as.integer(num_na %/% n_rows)
    na_per_col <- rep.int(as.integer(n_rows), n_cols)
    remainder <- num_na %% n_rows
    # last `remainder` columns get a +1 bump; rev() so bumped buckets
    # end up at the tail and are consumed last (smallest-first greedy)
    na_per_col[seq_len(remainder)] <- na_per_col[seq_len(remainder)] + 1L
    na_per_col <- rev(na_per_col)
  } else {
    na_per_col <- rep.int(as.integer(n_rows), n_cols)
    num_na <- n_rows * n_cols
  }
  max_row_miss <- max(na_per_col)
  # for each column, we keep 2 values, so max_row_miss can't exceed this hard bound
  checkmate::assert_true(max_row_miss <= nrow(obj) - 2L, .var.name = "max_row_miss <= nrow(obj)")

  # resolve `subset` into integer pool_idx
  pool_idx <- if (is.null(subset)) {
    seq_len(ncol(obj))
  } else if (is.character(subset)) {
    checkmate::assert_character(
      subset,
      any.missing = FALSE, min.len = 1, unique = TRUE, .var.name = "subset"
    )
    if (is.null(colnames(obj))) {
      cli::cli_abort("{.arg subset} is character but {.arg obj} has no colnames.")
    }
    missing_cols <- setdiff(subset, colnames(obj))
    if (length(missing_cols)) {
      cli::cli_abort(
        "{.arg subset} contains colnames not in {.arg obj}: {fmt_trunc(missing_cols, 6)}"
      )
    }
    match(subset, colnames(obj))
  } else if (is.numeric(subset)) {
    checkmate::assert_integerish(
      subset,
      lower = 1, upper = ncol(obj),
      any.missing = FALSE, min.len = 1, unique = TRUE, .var.name = "subset"
    )
    as.integer(subset)
  } else {
    cli::cli_abort("{.arg subset} must be numeric, character, or NULL.")
  }

  if (n_cols > length(pool_idx)) {
    cli::cli_abort(
      "Cannot stratify across {n_cols} columns; {.arg subset} has only {length(pool_idx)}."
    )
  }

  # pre-injection state and global feasibility (checked across the full obj,
  # since imputation uses all columns — untouched cols must also be healthy)
  not_na_mat <- !is.na(obj)
  max_allowed_col_miss <- floor(nrow(obj) * colmax)
  max_allowed_row_miss <- floor(ncol(obj) * rowmax)
  current_col_miss <- nrow(obj) - colSums(not_na_mat)
  current_row_miss <- ncol(obj) - rowSums(not_na_mat)
  current_col_vars <- col_vars(obj)

  if (any(current_col_miss > max_allowed_col_miss)) {
    bad <- which(current_col_miss > max_allowed_col_miss)
    cli::cli_abort(
      "Columns already exceed {.arg colmax} before injection: {fmt_trunc(bad, 6)}"
    )
  }
  if (any(current_row_miss > max_allowed_row_miss)) {
    bad <- which(current_row_miss > max_allowed_row_miss)
    cli::cli_abort(
      "Rows already exceed {.arg rowmax} before injection: {fmt_trunc(bad, 6)}"
    )
  }
  if (any(is.na(current_col_vars) | current_col_vars <= 0)) {
    bad <- which(is.na(current_col_vars) | current_col_vars <= 0)
    cli::cli_abort(
      "Columns already have zero variance before injection: {fmt_trunc(bad, 6)}"
    )
  }

  row_room <- max_allowed_row_miss - current_row_miss
  col_room <- max_allowed_col_miss - current_col_miss

  replicate(
    n_reps,
    sample_each_rep(
      obj = obj,
      pool_idx = pool_idx,
      na_per_col = na_per_col,
      row_room = row_room,
      col_room = col_room,
      max_attempts = max_attempts
    ),
    simplify = FALSE
  )
}

#' Resolve NA locations for tune_imp
#'
#' If `na_loc` is NULL, generates `n_reps` random NA location matrices via
#' `sample_na_loc_stratified()`. Otherwise normalizes the user-supplied
#' positions into a list of 2-column (row, col) integer matrices and
#' bounds-checks them.
#'
#' Accepted `na_loc` shapes:
#' - 2-column integerish matrix of (row, col) pairs (treated as a single rep)
#' - integer vector of linear positions (treated as a single rep)
#' - list whose elements are either of the above
#'
#' @noRd
#' @keywords internal
resolve_na_loc <- function(
  obj,
  na_loc,
  n_reps,
  num_na,
  n_cols,
  n_rows,
  rowmax,
  colmax,
  subset,
  max_attempts
) {
  if (is.null(na_loc)) {
    return(sample_na_loc_stratified(
      obj = obj,
      n_cols = n_cols,
      n_rows = n_rows,
      num_na = num_na,
      n_reps = n_reps,
      rowmax = rowmax,
      colmax = colmax,
      subset = subset,
      max_attempts = max_attempts
    ))
  }

  if (!is.list(na_loc)) {
    if (is.matrix(na_loc) || is.numeric(na_loc)) {
      na_loc <- list(na_loc)
    } else {
      stop("`na_loc` must be a 2-col matrix, an integer vector, or a list of these.",
        call. = FALSE
      )
    }
  }

  checkmate::assert_list(na_loc, min.len = 1, .var.name = "na_loc")

  nr <- nrow(obj)
  nc <- ncol(obj)
  n_total <- length(obj)

  for (i in seq_along(na_loc)) {
    elem <- na_loc[[i]]
    nm <- sprintf("na_loc[[%d]]", i)
    if (is.matrix(elem)) {
      checkmate::assert_matrix(
        elem,
        mode = "integerish", ncols = 2, min.rows = 1, .var.name = nm
      )
      checkmate::assert_true(
        all(elem[, 1] >= 1 & elem[, 1] <= nr),
        .var.name = paste0(nm, " row indices")
      )
      checkmate::assert_true(
        all(elem[, 2] >= 1 & elem[, 2] <= nc),
        .var.name = paste0(nm, " col indices")
      )
    } else {
      checkmate::assert_integerish(
        elem,
        lower = 1, upper = n_total, any.missing = FALSE,
        min.len = 1, unique = TRUE, .var.name = nm
      )
    }
  }

  elem_sizes <- vapply(
    na_loc,
    function(x) if (is.matrix(x)) nrow(x) else length(x),
    numeric(1)
  )
  if (length(unique(elem_sizes)) != 1L) {
    stop("All elements in `na_loc` must specify the same number of NA positions.",
      call. = FALSE
    )
  }

  na_loc
}

#' Tune Parameters for Imputation Methods
#'
#' Tunes hyperparameters for imputation methods such as [slide_imp()], [knn_imp()],
#' [pca_imp()], or user-supplied custom functions by repeated cross-validation.
#'
#' @details
#' The function supports tuning for built-in imputation methods ("slide_imp",
#' "knn_imp", "pca_imp") or custom functions provided via `.f`.
#'
#' When `.f` is a character string, the columns in `parameters` are validated
#' against the chosen method's requirements:
#' - `"knn_imp"`: requires `k` in `parameters`
#' - `"pca_imp"`: requires `ncp` in `parameters`
#' - `"slide_imp"`: requires `window_size`, `overlap_size`, and `min_window_n`,
#' plus exactly one of `k` or `ncp`
#'
#' When `.f` is a custom function, the columns in `parameters` must correspond
#' to the arguments of `.f` (excluding the `obj` argument). The custom function
#' must accept `obj` (a numeric matrix) as its first argument and return a
#' numeric matrix of identical dimensions.
#'
#' Tuning results can be evaluated using the \pkg{yardstick} package or
#' [compute_metrics()].
#'
#' @inheritParams slide_imp
#' @inheritParams group_imp
#'
#' @param parameters A data.frame specifying parameter combinations to tune,
#'   where each column represents a parameter accepted by `.f` (excluding `obj`).
#'   List columns are supported for complex parameters. Duplicate rows are
#'   automatically removed. `NULL` is treated as tuning the function with its
#'   default parameters.
#' @param n_reps Integer. Number of repetitions for random NA injection
#'   (default `1`).
#' @param rowmax,colmax Numbers between 0 and 1. NA injection cannot create
#'   rows/columns with a higher proportion of missing values than these
#'   thresholds.
#' @param n_cols Integer. The number of columns to receive injected NAs per
#'   repetition. Ignored when `num_na` is supplied (in which case `n_cols` is
#'   derived as `as.integer(num_na %/% n_rows)`). Must be provided if
#'   `num_na` is `NULL`. Ignored when `na_loc` is supplied.
#' @param n_rows Integer. The target number of NAs per column (default `1L`).
#'   - When `num_na` is supplied: used as the base size. Most columns receive
#'     exactly `n_rows` NAs; `num_na %% n_rows` columns receive `n_rows + 1`.
#'   - When `num_na` is `NULL`: every selected column receives exactly
#'     `n_rows` NAs.
#'   Ignored when `na_loc` is supplied.
#' @param num_na Integer. Total number of missing values to inject per
#'   repetition. If supplied, `n_cols` is computed automatically and the NAs
#'   are distributed as evenly as possible using `n_rows` as the base
#'   (`num_na` must be `>= n_rows`). If omitted, exactly `n_cols * n_rows`
#'   NAs are injected. At least one of `n_cols` or `num_na` must be supplied
#'   when `na_loc` is `NULL`. Ignored when `na_loc` is supplied.
#' @param subset Optional integer or character vector restricting which columns
#'   of `obj` are eligible for NA injection.
#'   - If `NULL` (default): all columns are eligible.
#'   - If character: values must exist in `colnames(obj)`.
#'   - If integer/numeric: values must be valid 1-based column indices.
#'   The vector must be unique and must contain at least `n_cols` columns
#'   (or the number derived from `num_na`).
#' @param max_attempts Integer. Maximum number of resampling attempts per
#'   repetition before giving up due to row-budget exhaustion (default `10`).
#' @param na_loc Optional. Pre-defined missing value locations to bypass random
#'   NA injection. Accepted formats include:
#'   - A two-column integer matrix (row, column indices).
#'   - A numeric vector of linear locations.
#'   - A list where each element is one of the above formats (one per repetition).
#' @param .progress Logical. Show a progress bar during tuning (default `TRUE`).
#' @param cores Controls the number of cores to parallelize over for K-NN and
#'   sliding-window K-NN imputation with OpenMP. For other methods, use
#'   `mirai::daemons()` instead.
#' @param location Required only for `slide_imp`. Numeric vector of column
#'   locations.
#' @param pin_blas Logical. Pin BLAS threads to 1 during parallel tuning
#'   (default `FALSE`).
#' @param .f The imputation method to tune. Either a character string
#'   (`"knn_imp"`, `"pca_imp"`, or `"slide_imp"`) or a custom function.
#'
#' @return A data.frame of class `c("slideimp_tune", "slideimp_tbl", "data.frame")`
#'   containing:
#'   - All columns from `parameters`
#'   - `param_set`: Unique parameter set ID
#'   - `rep`: Repetition index (from `n_reps`)
#'   - `result`: A nested list-column of data.frames with `truth` and `estimate`
#'     columns
#'   - `error`: Character column with failure reason (or `NA` if successful)
#'
#' @examples
#' data(khanmiss1)
#' obj <- t(khanmiss1)[1:20, sample.int(nrow(khanmiss1), size = 200)]
#'
#' # Tune full K-NN imputation
#' parameters <- data.frame(k = c(5, 10))
#'
#' # With random NA injection
#' results <- tune_imp(obj, parameters, .f = "knn_imp", n_reps = 1, num_na = 20)
#'
#' # Compute metrics on results
#' compute_metrics(results)
#'
#' # Tune with fixed NA positions (2 repetitions)
#' na_positions <- list(
#'   matrix(c(1, 2, 3, 1, 1, 1), ncol = 2), # Rows 1-3 in column 1
#'   matrix(c(2, 3, 4, 2, 2, 2), ncol = 2) # Rows 2-4 in column 2
#' )
#' results_fixed <- tune_imp(
#'   obj,
#'   data.frame(k = 10),
#'   .f = "knn_imp",
#'   na_loc = na_positions
#' )
#'
#' compute_metrics(results_fixed)
#'
#' # Custom imputation function example
#' custom_imp <- function(obj, mean = 0, sd = 1) {
#'   na_pos <- is.na(obj)
#'   obj[na_pos] <- stats::rnorm(sum(na_pos), mean = mean, sd = sd)
#'   obj
#' }
#'
#' # Setup 2 cores for parallelization (mirai)
#' mirai::daemons(2)
#' parameters_custom <- data.frame(mean = c(0, 0, 1), sd = c(1, 2, 1))
#' results_custom <- tune_imp(
#'   obj,
#'   parameters_custom,
#'   .f = custom_imp,
#'   n_reps = 2,
#'   num_na = 20
#' )
#' mirai::daemons(0)
#' compute_metrics(results_custom)
#'
#' @export
tune_imp <- function(
  obj,
  parameters = NULL,
  .f,
  n_reps = 1,
  na_loc = NULL,
  n_cols = NULL,
  n_rows = 1,
  num_na = NULL,
  rowmax = 0.9,
  colmax = 0.9,
  subset = NULL,
  max_attempts = 10L,
  .progress = TRUE,
  cores = 1,
  location = NULL,
  pin_blas = FALSE
) {
  # pre-conditioning
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    min.rows = 1,
    min.cols = 2,
    null.ok = FALSE,
    .var.name = "obj"
  )
  nr <- nrow(obj)
  nc <- ncol(obj)

  checkmate::assert_true(sum(is.infinite(obj)) == 0, .var.name = "obj")
  if (is.null(parameters)) {
    parameters_is_null <- TRUE
    parameters <- data.frame(.placeholder = TRUE)
  } else {
    parameters_is_null <- FALSE
    checkmate::assert_data_frame(
      parameters,
      any.missing = FALSE,
      all.missing = FALSE,
      min.rows = 1,
      col.names = "unique",
      .var.name = "parameters",
      null.ok = FALSE
    )
    parameters <- unique(parameters)
  }
  checkmate::assert_flag(pin_blas, null.ok = FALSE, .var.name = "pin_blas")

  # validate .f
  if (is.character(.f)) {
    stopifnot(
      "`.f` must be 'slide_imp', 'knn_imp', or 'pca_imp'." =
        .f %in% c("slide_imp", "knn_imp", "pca_imp") && length(.f) == 1
    )

    if (.f == "slide_imp") {
      if (is.null(location)) {
        stop("`slide_imp` requires the `location` argument to be provided.")
      }
      checkmate::assert_numeric(
        location,
        len = nc,
        any.missing = FALSE,
        sorted = TRUE,
        finite = TRUE,
        null.ok = FALSE,
        .var.name = "location"
      )
      if (!all(c("window_size", "overlap_size", "min_window_n") %in% names(parameters))) {
        stop("`slide_imp` requires `window_size`, `overlap_size`, and `min_window_n` in `parameters`.")
      }
      if ("flank" %in% names(parameters)) {
        stop("Tuning with `flank` passed as parameters is not yet supported.")
      }
      if ("k" %in% names(parameters) && "ncp" %in% names(parameters)) {
        stop("`slide_imp` requires exactly one of `k` or `ncp` in `parameters`, not both.")
      }
      if (!any(c("k", "ncp") %in% names(parameters))) {
        stop("`slide_imp` requires either `k` or `ncp` in `parameters`.")
      }
      if (!".progress" %in% names(parameters)) {
        parameters$.progress <- FALSE
      }
    } else if (.f == "knn_imp") {
      if (!"k" %in% names(parameters)) {
        stop("`knn_imp` requires `k` in `parameters`.")
      }
      if ("ncp" %in% names(parameters)) {
        stop("`knn_imp` does not accept `ncp`. Did you mean `.f = 'pca_imp'`?")
      }
    } else if (.f == "pca_imp") {
      if (!"ncp" %in% names(parameters)) {
        stop("`pca_imp` requires `ncp` in `parameters`.")
      }
      if ("k" %in% names(parameters)) {
        stop("`pca_imp` does not accept `k`. Did you mean `.f = 'knn_imp'`?")
      }
    }

    message(sprintf("Tuning %s", .f))
  } else if (is.function(.f)) {
    message("Tuning custom function")
  } else {
    stop("`.f` must be a function or one of 'slide_imp', 'knn_imp', 'pca_imp'.")
  }

  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")

  # parallelization mode flags
  is_knn_mode <- {
    (is.character(.f) && .f == "knn_imp") ||
      (is.character(.f) && .f == "slide_imp" && "k" %in% names(parameters))
  }

  message("Step 1/2: Resolving NA locations")

  na_loc <- resolve_na_loc(
    obj = obj,
    na_loc = na_loc,
    n_reps = n_reps,
    num_na = num_na,
    n_cols = n_cols,
    n_rows = n_rows,
    rowmax = rowmax,
    colmax = colmax,
    subset = subset,
    max_attempts = max_attempts
  )
  n_reps <- length(na_loc)

  .rowid <- seq_len(nrow(parameters))

  indices <- expand.grid(
    param_set = .rowid,
    rep_id = seq_len(n_reps)
  )

  parallelize <- tryCatch(mirai::require_daemons(), error = function(e) FALSE)
  if (cores > 1 && is_knn_mode) {
    if (!has_openmp()) {
      message("OpenMP not available. KNN will run single-threaded.")
      cores <- 1
    } else if (parallelize) {
      message(
        "Both `cores > 1` and `mirai::daemons()` detected. ",
        "Setting `cores = 1` to avoid nested parallelism. ",
        "Parallelization will be handled by `mirai`."
      )
      cores <- 1
    }
  }

  if (!parallelize && cores > 1 && !is_knn_mode) {
    message(
      sprintf(
        "cores = %d but is ignored for non-KNN imputation. Call `mirai::daemons(%d)` to set up parallelization.",
        cores,
        cores
      )
    )
  }

  if (parallelize || (is_knn_mode && cores > 1)) {
    message("Running Mode: parallel...")
  } else {
    message("Running Mode: sequential...")
  }

  is_slide <- is.character(.f) && .f == "slide_imp"
  fixed_args <- list()
  if (is_slide) fixed_args$location <- location
  if (is_knn_mode) fixed_args$cores <- cores

  parameters_list <- lapply(split(parameters, f = as.factor(.rowid)), function(row) {
    row_list <- as.list(row)
    if (parameters_is_null) {
      row_list$.placeholder <- NULL
    }
    row_list <- lapply(row_list, function(x) {
      if (is.list(x) && length(x) == 1) {
        x[[1]]
      } else {
        x
      }
    })
    row_list
  })

  if (is.character(.f)) {
    impute_f <- switch(.f,
      slide_imp = slide_imp,
      knn_imp   = knn_imp,
      pca_imp   = pca_imp
    )
  } else {
    impute_f <- .f
    pre <- obj
    na_positions <- na_loc[[1]]
    pre[na_positions] <- NA
    param_vec <- parameters_list[[1]]
    probe_result <- do.call(
      impute_f,
      c(list(obj = pre), fixed_args, param_vec)
    )
    tryCatch(
      {
        checkmate::assert_matrix(
          probe_result,
          mode = "numeric",
          nrows = nr, ncols = nc, null.ok = FALSE,
          .var.name = "imputed_result"
        )
        checkmate::assert_true(
          sum(is.infinite(probe_result)) == 0,
          .var.name = "imputed_result"
        )
      },
      error = function(e) {
        stop(
          "Custom `.f` must accept `obj` as its first argument and return a ",
          "numeric matrix with the same dimensions as `obj` and no Inf values.",
          call. = FALSE
        )
      }
    )
  }

  message("Step 2/2: Tuning\n")

  iter <- seq_len(nrow(indices))

  if (parallelize) {
    check_pin_blas(pin_blas)
    crated_fn <- carrier::crate(
      function(i) {
        if (pin_blas) {
          RhpcBLASctl::blas_set_num_threads(1)
          RhpcBLASctl::omp_set_num_threads(1)
        }
        tryCatch(
          {
            pre <- obj
            na_positions <- na_loc[[indices[i, "rep_id", drop = TRUE]]]
            pre[na_positions] <- NA
            truth_vec <- obj[na_positions]
            param_vec <- parameters_list[[indices[i, "param_set", drop = TRUE]]]
            imputed_result <- do.call(
              impute_f,
              args = c(list(obj = pre), fixed_args, param_vec)
            )
            estimate_vec <- imputed_result[na_positions]
            data.frame(truth = truth_vec, estimate = estimate_vec)
          },
          error = function(e) {
            out <- data.frame(truth = numeric(), estimate = numeric())
            attr(out, "error") <- conditionMessage(e)
            out
          }
        )
      },
      impute_f = impute_f,
      obj = obj,
      na_loc = na_loc,
      indices = indices,
      nr = nr,
      nc = nc,
      parameters_list = parameters_list,
      fixed_args = fixed_args,
      pin_blas = pin_blas
    )
    result_list <- mirai::mirai_map(iter, crated_fn)[.progress = .progress]
  } else {
    run_sequential <- function(i) {
      tryCatch(
        {
          pre <- obj
          na_positions <- na_loc[[indices[i, "rep_id", drop = TRUE]]]
          pre[na_positions] <- NA
          truth_vec <- obj[na_positions]
          param_vec <- parameters_list[[indices[i, "param_set", drop = TRUE]]]
          imputed_result <- do.call(
            impute_f,
            args = c(list(obj = pre), fixed_args, param_vec)
          )
          estimate_vec <- imputed_result[na_positions]
          data.frame(truth = truth_vec, estimate = estimate_vec)
        },
        error = function(e) {
          out <- data.frame(truth = numeric(), estimate = numeric())
          attr(out, "error") <- conditionMessage(e)
          out
        }
      )
    }
    if (.progress) pb <- cli::cli_progress_bar(total = length(iter))
    result_list <- lapply(iter, function(i) {
      out <- run_sequential(i)
      if (.progress) cli::cli_progress_update(id = pb)
      out
    })
    if (.progress) cli::cli_progress_done(id = pb)
  }

  error_vec <- vapply(result_list, function(x) {
    e <- attr(x, "error")
    if (is.null(e)) {
      NA_character_
    } else {
      e
    }
  }, character(1))
  result_df <- cbind(
    parameters[indices$param_set, , drop = FALSE],
    indices,
    data.frame(result = I(result_list), error = error_vec)
  )

  if (any(!is.na(result_df$error))) {
    warning("Some tuning iterations failed. Check the 'error' column in the output.",
      call. = FALSE
    )
  }

  if (parameters_is_null) {
    result_df$.placeholder <- NULL
  }

  class(result_df) <- c("slideimp_tune", "slideimp_tbl", "data.frame")
  return(result_df)
}

# compute_metrics ----
calc_mae <- function(truth, estimate) {
  mean(abs(truth - estimate), na.rm = TRUE)
}

calc_rmse <- function(truth, estimate) {
  sqrt(mean((truth - estimate)^2, na.rm = TRUE))
}

calc_rsq <- function(truth, estimate) {
  valid <- !is.na(truth) & !is.na(estimate)
  if (sum(valid) < 2) {
    return(NA_real_)
  }
  stats::cor(truth[valid], estimate[valid])^2
}

calc_rsq_trad <- function(truth, estimate) {
  valid <- !is.na(truth) & !is.na(estimate)
  if (sum(valid) < 2) {
    return(NA_real_)
  }
  truth_valid <- truth[valid]
  estimate_valid <- estimate[valid]
  ss_res <- sum((truth_valid - estimate_valid)^2)
  mean_truth <- mean(truth_valid)
  ss_tot <- sum((truth_valid - mean_truth)^2)
  if (abs(ss_tot) < .Machine$double.eps) {
    return(NA_real_)
  }
  1 - (ss_res / ss_tot)
}

calc_mape <- function(truth, estimate) {
  mean(abs((truth - estimate) / truth), na.rm = TRUE) * 100
}

calc_bias <- function(truth, estimate) {
  mean(estimate - truth, na.rm = TRUE)
}

calc_all_metrics <- function(x, metric_fns) {
  estimates <- vapply(
    metric_fns,
    function(fn) fn(x$truth, x$estimate),
    numeric(1)
  )

  data.frame(
    .metric = names(metric_fns),
    .estimator = "standard",
    .estimate = estimates
  )
}

#' Compute Prediction Accuracy Metrics
#'
#' Computes prediction accuracy metrics for results from [tune_imp()].
#'
#' For alternative or faster metrics, see the `{yardstick}` package.
#'
#' @param results A `slideimp_tune` data.frame from [tune_imp()] containing
#' a `result` list-column with data.frames that have `truth` and `estimate` columns.
#'
#' @param metrics A character vector of metric names to compute. Defaults
#' to `c("mae", "rmse")`. Also available: `"mape"`, `"bias"`, `"rsq"`, and `"rsq_trad"`.
#'
#' @return A data.frame with the original parameters and un-nested metrics
#' (`.metric`, `.estimator`, `.estimate`).
#'
#' @examples
#' data(khanmiss1)
#'
#' set.seed(1234)
#' results <- tune_imp(
#'   obj = t(khanmiss1),
#'   parameters = data.frame(k = 10),
#'   .f = "knn_imp",
#'   n_reps = 1,
#'   num_na = 20
#' )
#'
#' compute_metrics(results)
#'
#' @export
compute_metrics <- function(results, metrics = c("mae", "rmse")) {
  UseMethod("compute_metrics")
}

#' @rdname compute_metrics
#' @export
compute_metrics.data.frame <- function(results, metrics = c("mae", "rmse")) {
  if (!"result" %in% names(results)) {
    stop("`results` must contain a 'result' column.")
  }
  first_result <- results$result[[1]]
  if (!is.data.frame(first_result) ||
    !all(c("truth", "estimate") %in% names(first_result))) {
    stop("Each element of 'result' must be a data.frame with 'truth' and 'estimate' columns.")
  }
  compute_metrics.slideimp_tune(results, metrics = metrics)
}

#' @rdname compute_metrics
#' @export
compute_metrics.slideimp_tune <- function(results, metrics = c("mae", "rmse")) {
  checkmate::assert_character(metrics, unique = TRUE)
  .metrics_list <- list(
    mae = calc_mae,
    rmse = calc_rmse,
    rsq = calc_rsq,
    mape = calc_mape,
    bias = calc_bias,
    rsq_trad = calc_rsq_trad
  )
  invalid_metrics <- setdiff(metrics, names(.metrics_list))
  if (length(invalid_metrics) > 0) {
    stop(
      "Unknown metrics: ", paste(invalid_metrics, collapse = ", "), "\n",
      "Available metrics: ", paste(names(.metrics_list), collapse = ", ")
    )
  }
  metric_fns <- .metrics_list[metrics]
  # always compute n and n_miss per result element
  results$n <- vapply(results$result, nrow, integer(1))
  results$n_miss <- vapply(results$result, function(x) sum(is.na(x$estimate)), integer(1))
  metrics_list <- lapply(
    results$result,
    \(x) calc_all_metrics(x, metric_fns = metric_fns)
  )
  keep_cols <- setdiff(names(results), "result")
  n_metrics <- vapply(metrics_list, nrow, integer(1))
  row_data <- results[rep(seq_len(nrow(results)), n_metrics), keep_cols, drop = FALSE]
  metric_df <- collapse::rowbind(metrics_list)
  out <- cbind(row_data, metric_df)
  rownames(out) <- NULL
  out
}


# #' Convert 2D Positions to Linear Indices
# #'
# #' This function converts 2D positions (row and column indices, i.e., first and
# #' second columns) in a matrix to their corresponding linear (1D) positions.
# #'
# #' @param pos_2d A numeric matrix with exactly 2 columns (first for rows,
# #'   second for columns) and at least 1 row. Each entry must be a positive
# #'   integer within the matrix bounds.
# #' @param nrow The number of rows in the matrix.
# #' @param ncol The number of columns in the matrix.
# #' @return A numeric vector of linear positions corresponding to the input 2D
# #'   positions.
# #'
# #' @keywords internal
# #' @noRd
# grid_to_linear <- function(pos_2d, nrow, ncol) {
#   checkmate::assert_matrix(pos_2d, ncols = 2, min.rows = 1, mode = "integerish", .var.name = "pos_2d")
#   checkmate::assert_int(nrow, lower = 1, "nrow")
#   checkmate::assert_int(ncol, lower = 1, "ncol")
#
#   row <- pos_2d[, 1]
#   col <- pos_2d[, 2]
#
#   checkmate::assert_true(all(row >= 1 & row <= nrow), .var.name = "row indices")
#   checkmate::assert_true(all(col >= 1 & col <= ncol), .var.name = "column indices")
#
#   linear_pos <- (col - 1) * nrow + row
#   return(linear_pos)
# }
