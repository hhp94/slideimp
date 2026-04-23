#' Sample Missing Value Locations with Constraints
#'
#' Samples indices for `NA` injection into a matrix while maintaining row/column
#' missing value budgets and avoiding zero-variance columns.
#'
#' @details
#' The function uses a greedy stochastic search for valid `NA` locations. It
#' ensures that:
#'
#' - Total missingness per row and column does not exceed `rowmax` and `colmax`.
#' - At least two distinct observed values are preserved in every column to
#'   ensure the column maintains non-zero variance.
#'
#' @inheritParams tune_imp
#'
#' @returns A list of length `n_reps`. Each element is a two-column integer
#' matrix (`row`, `col`) representing the coordinates of the sampled `NA`
#' locations.
#'
#' @examples
#' mat <- matrix(runif(100), nrow = 10)
#'
#' # Sample 5 `NA` across 5 columns (1 per column)
#' locs <- sample_na_loc(mat, n_cols = 5, n_rows = 1)
#' locs
#'
#' # Inject the `NA` from the first repetition
#' mat[locs[[1]]] <- NA
#' mat
#'
#' @export
sample_na_loc <- function(
  obj,
  n_cols = NULL,
  n_rows = 2L,
  num_na = NULL,
  n_reps = 1L,
  rowmax = 0.9,
  colmax = 0.9,
  na_col_subset = NULL,
  max_attempts = 100
) {
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
    cli::cli_abort(
      "{.arg num_na} ({num_na}) must be >= {.arg n_rows} ({n_rows})
     (each selected column needs at least one NA)."
    )
  }

  # resolve (n_cols, na_per_col, max_row_miss)
  if (!is.null(num_na)) {
    n_cols <- as.integer(num_na %/% n_rows)
    remainder <- num_na %% n_rows
    # distribute remainder across the n_cols buckets: each gets an extra
    # `remainder %/% n_cols`, and the first `remainder %% n_cols` get +1 more.
    base_extra <- remainder %/% n_cols
    leftover <- remainder %% n_cols
    na_per_col <- rep.int(as.integer(n_rows + base_extra), n_cols)
    if (leftover > 0L) {
      na_per_col[seq_len(leftover)] <- na_per_col[seq_len(leftover)] + 1L
    }
    # sort() so bumped buckets end up at the tail and are consumed last
    # (smallest-first greedy)
    na_per_col <- sort(na_per_col)
  } else {
    na_per_col <- rep.int(as.integer(n_rows), n_cols)
    num_na <- n_rows * n_cols
  }
  max_row_miss <- max(na_per_col)
  # for each column, we keep 2 values, so max_row_miss can't exceed this hard bound
  max_row_miss <- max(na_per_col)
  # for each column, we keep 2 values, so max_row_miss can't exceed this hard bound
  if (max_row_miss > nrow(obj) - 2L) {
    cli::cli_abort(c(
      "Too many NAs requested per column.",
      "x" = "Maximum of {max_row_miss} NAs per column (from {.arg n_rows} = {n_rows}),
             but must keep at least 2 observed values per column to preserve variance.",
      "i" = "Matrix has only {nrow(obj)} rows."
    ))
  }

  # resolve `na_col_subset` into integer pool_idx
  pool_idx <- if (is.null(na_col_subset)) {
    seq_len(ncol(obj))
  } else if (is.character(na_col_subset)) {
    checkmate::assert_character(
      na_col_subset,
      any.missing = FALSE, min.len = 1, unique = TRUE, .var.name = "na_col_subset"
    )
    if (is.null(colnames(obj))) {
      cli::cli_abort("{.arg na_col_subset} is character but {.arg obj} has no colnames.")
    }
    missing_cols <- setdiff(na_col_subset, colnames(obj))
    if (length(missing_cols)) {
      cli::cli_abort(
        "{.arg na_col_subset} contains colnames not in {.arg obj}: {fmt_trunc(missing_cols, 6)}"
      )
    }
    match(na_col_subset, colnames(obj))
  } else if (is.numeric(na_col_subset)) {
    checkmate::assert_integerish(
      na_col_subset,
      lower = 1, upper = ncol(obj),
      any.missing = FALSE, min.len = 1, unique = TRUE, .var.name = "na_col_subset"
    )
    as.integer(na_col_subset)
  } else {
    cli::cli_abort("{.arg na_col_subset} must be numeric, character, or NULL.")
  }

  if (n_cols > length(pool_idx)) {
    if (!is.null(num_na)) {
      cli::cli_abort(c(
        "Cannot place {num_na} NAs (requires {n_cols} columns given {.arg n_rows} = {n_rows}).",
        "i" = "{.arg na_col_subset} provides only {length(pool_idx)} column{?s}."
      ))
    } else {
      cli::cli_abort(c(
        "Cannot stratify across {n_cols} columns.",
        "i" = "{.arg na_col_subset} provides only {length(pool_idx)} column{?s}."
      ))
    }
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
    cli::cli_abort(c(
      "Some columns already exceed {.arg colmax} before injection.",
      "x" = "Problematic columns: {fmt_trunc(bad, 6)}"
    ))
  }
  if (any(current_row_miss > max_allowed_row_miss)) {
    bad <- which(current_row_miss > max_allowed_row_miss)
    cli::cli_abort(c(
      "Some rows already exceed {.arg rowmax} before injection.",
      "x" = "Problematic rows: {fmt_trunc(bad, 6)}"
    ))
  }
  if (any(is.na(current_col_vars) | current_col_vars <= 0)) {
    bad <- which(is.na(current_col_vars) | current_col_vars <= 0)
    cli::cli_abort(c(
      "Some columns already have zero (or NA) variance before injection.",
      "x" = "Problematic columns: {fmt_trunc(bad, 6)}",
      "i" = "NA injection requires at least 2 distinct observed values per column."
    ))
  }

  row_room <- max_allowed_row_miss - current_row_miss
  col_room <- max_allowed_col_miss - current_col_miss

  pool_idx0 <- pool_idx - 1L
  lapply(seq_len(n_reps), function(x) {
    result <- sample_each_rep_cpp(
      obj = obj,
      pool_idx_in = pool_idx0,
      na_per_col = na_per_col,
      row_room = row_room,
      col_room = col_room,
      max_attempts = max_attempts
    )
    colnames(result) <- c("row", "col")
    result
  })
}

#' Resolve NA Locations for `tune_imp()`
#'
#' If `na_loc` is NULL, generates `n_reps` random NA location matrices via
#' `sample_na_loc()`. Otherwise normalizes the user-supplied
#' positions into a list of 2-column (row, col) integer matrices and
#' bounds-checks them.
#'
#' @inheritParams tune_imp
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
  na_col_subset,
  max_attempts
) {
  if (is.null(na_loc) && is.null(n_cols) && is.null(num_na)) {
    n_total <- nrow(obj) * ncol(obj)
    pct_based <- max(as.integer(n_rows), as.integer(round(0.05 * n_total)))
    num_na <- min(500L, pct_based)
    pct_actual <- round(100 * num_na / n_total, 1)
    cli::cli_inform(c(
      "i" = "Using default {.arg num_na} = {num_na} (~{pct_actual}% of cells).",
      " " = "Increase for more reliability or decrease if missing is dense."
    ))
  }

  if (is.null(na_loc)) {
    return(sample_na_loc(
      obj = obj,
      n_cols = n_cols,
      n_rows = n_rows,
      num_na = num_na,
      n_reps = n_reps,
      rowmax = rowmax,
      colmax = colmax,
      na_col_subset = na_col_subset,
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

  # canonical form is going to be a 2D index
  na_loc <- lapply(na_loc, function(elem) {
    if (is.matrix(elem)) {
      storage.mode(elem) <- "integer"
      colnames(elem) <- c("row", "col")
      elem
    } else {
      lp <- as.integer(elem)
      cbind(
        row = ((lp - 1L) %% nr) + 1L,
        col = ((lp - 1L) %/% nr) + 1L
      )
    }
  })

  return(na_loc)
}

#' Tune Parameters for Imputation Methods
#'
#' Tunes hyperparameters for imputation methods such as [slide_imp()],
#' [knn_imp()], [pca_imp()], or user-supplied custom functions by repeated
#' cross-validation. For [group_imp()], tune [knn_imp()] or [pca_imp()] on
#' a single group.
#'
#' @inheritParams slide_imp
#' @inheritParams group_imp
#'
#' @param parameters A data.frame specifying parameter combinations to tune,
#'   where each column represents a parameter accepted by `.f` (excluding
#'   `obj`). List columns are supported for complex parameters. Duplicate
#'   rows are automatically removed. `NULL` is treated as tuning the function
#'   with its default parameters.
#' @param n_reps Integer. Number of repetitions for random NA injection
#'   (default `1`).
#' @param rowmax,colmax Numbers between 0 and 1. NA injection cannot create
#'   rows/columns with a higher proportion of missing values than these
#'   thresholds.
#' @param n_cols Integer. The number of columns to receive injected `NA` per
#'   repetition. Ignored when `num_na` is supplied (in which case `n_cols`
#'   is derived as `num_na %/% n_rows`). Must be provided if `num_na` is
#'   `NULL`. Ignored in [tune_imp()] when `na_loc` is supplied.
#' @param n_rows Integer. The target number of `NA` values per column
#'   (default `2L`).
#'   - When `num_na` is supplied: used as the base size. Most columns
#'     receive exactly `n_rows` missing values; `num_na %% n_rows` columns
#'     receive one extra. If there's only one column, it receives all
#'     the remainder.
#'   - When `num_na` is `NULL`: every selected column receives exactly
#'     `n_rows` `NA`.
#'   Ignored in [tune_imp()] when `na_loc` is supplied.
#'
#' @param num_na Integer. Total number of missing values to inject per
#'   repetition. If supplied, `n_cols` is computed automatically and
#'   missing values are distributed as evenly as possible, using `n_rows`
#'   as the base size (`num_na` must be at least `n_rows`). If omitted
#'   but `n_cols` is supplied, the total injected is `n_cols * n_rows`.
#'   If `num_na`, `n_cols`, and `na_loc` are all `NULL`, [tune_imp()]
#'   defaults to roughly 5% of total cells, capped at 500.
#'   [sample_na_loc()] has no default. Ignored in [tune_imp()] when
#'   `na_loc` is supplied.
#'
#' @param na_col_subset Optional integer or character vector restricting which
#'   columns of `obj` are eligible for NA injection.
#'   - If `NULL` (default): all columns are eligible.
#'   - If character: values must exist in `colnames(obj)`.
#'   - If integer/numeric: values must be valid 1-based column indices.
#'   The vector must be unique and must contain at least `n_cols` columns
#'   (or the number derived from `num_na`).
#'   Ignored in [tune_imp()] when `na_loc` is supplied.
#'
#' @param max_attempts Integer. Maximum number of resampling attempts per
#'   repetition before giving up due to row-budget exhaustion (default `100`).
#' @param na_loc Optional. Pre-defined missing value locations to bypass
#'   random NA injection with [sample_na_loc()]. Accepted formats include:
#'   - A two-column integer matrix (row, column indices).
#'   - A numeric vector of linear locations.
#'   - A list where each element is one of the above formats (one per
#'     repetition).
#'
#' @param .progress Logical. Show a progress bar during tuning
#'   (default `TRUE`).
#'
#' @param cores Controls the number of cores to parallelize over for K-NN
#'   and sliding-window K-NN imputation with OpenMP. For other methods, use
#'   `mirai::daemons()` instead.
#'
#' @param location Required only for `slide_imp`. Numeric vector of column
#'   locations.
#'
#' @param pin_blas Logical. Pin BLAS threads to 1 during parallel tuning
#'   (default `FALSE`).
#'
#' @param .f Either `"knn_imp"`, `"pca_imp"`, `"slide_imp"`, or a custom
#'   function specifying the imputation method to tune.
#'
#' @details
#' The function supports tuning for built-in methods (`"slide_imp"`,
#' `"knn_imp"`, `"pca_imp"`) or custom functions provided via `.f`.
#'
#' When `.f` is a character string, the columns in `parameters` are validated
#' against the chosen method's requirements:
#' - `"knn_imp"`: requires `k` in `parameters`
#' - `"pca_imp"`: requires `ncp` in `parameters`
#' - `"slide_imp"`: requires `window_size`, `overlap_size`, and
#'   `min_window_n`, plus exactly one of `k` or `ncp`
#'
#' When `.f` is a custom function, the columns in `parameters` must
#' correspond to the arguments of `.f` (excluding the `obj` argument). The
#' custom function must accept `obj` (a numeric matrix) as its first argument
#' and return a numeric matrix of identical dimensions.
#'
#' Tuning results can be evaluated using the `yardstick` package or
#' [compute_metrics()].
#'
#' @inheritSection group_imp Parallelization
#'
#' @returns A `data.frame` of class `slideimp_tune` containing:
#'   - `...`: All columns originally provided in `parameters`.
#'   - `param_set`: An integer ID for the unique parameter combination.
#'   - `rep_id`: An integer indicating the repetition index.
#'   - `result`: A nested list-column where each element is a data.frame
#'     containing `truth` (original values) and `estimate` (imputed values).
#'   - `error`: A character column containing the error message if the
#'     iteration failed, otherwise `NA`.
#'
#' @examples
#' # Setup example data. Increase `num_na` (500) and `n_reps` (10-30) in real
#' # analyses
#' obj <- sim_mat(10, 50)$input
#'
#' # 1. Tune K-NN imputation with random NA injection
#' params_knn <- data.frame(k = c(2, 4))
#' results <- tune_imp(obj, params_knn, .f = "knn_imp", n_reps = 1, num_na = 10)
#' compute_metrics(results)
#'
#' # 2. Tune with fixed NA positions
#' na_positions <- list(
#'   matrix(c(1, 2, 3, 1, 1, 1), ncol = 2),
#'   matrix(c(2, 3, 4, 2, 2, 2), ncol = 2)
#' )
#'
#' results_fixed <- tune_imp(
#'   obj,
#'   data.frame(k = 2),
#'   .f = "knn_imp",
#'   na_loc = na_positions
#' )
#'
#' # 3. Custom imputation function
#' custom_fill <- function(obj, val = 0) {
#'   obj[is.na(obj)] <- val
#'   obj
#' }
#' tune_imp(obj, data.frame(val = c(0, 1)), .f = custom_fill, num_na = 10)
#' @examplesIf interactive() && requireNamespace("mirai", quietly = TRUE)
#' # 4. Parallel tuning (requires mirai package)
#' mirai::daemons(2)
#' parameters_custom <- data.frame(mean = c(0, 1), sd = c(1, 1))
#'
#' # Define a simple custom function for illustration
#' custom_imp <- function(obj, mean, sd) {
#'   na_pos <- is.na(obj)
#'   obj[na_pos] <- stats::rnorm(sum(na_pos), mean = mean, sd = sd)
#'   obj
#' }
#'
#' results_p <- tune_imp(
#'   obj, parameters_custom,
#'   .f = custom_imp, n_reps = 1, num_na = 10
#' )
#' mirai::daemons(0) # Close workers
#' @export
tune_imp <- function(
  obj,
  parameters = NULL,
  .f,
  na_loc = NULL,
  num_na = NULL,
  n_reps = 1,
  n_cols = NULL,
  n_rows = 2,
  rowmax = 0.9,
  colmax = 0.9,
  na_col_subset = NULL,
  max_attempts = 100,
  .progress = TRUE,
  cores = 1,
  location = NULL,
  pin_blas = FALSE
) {
  # this function sometime leaves behind big matrices. This can help clean up.
  on.exit(gc(verbose = FALSE), add = TRUE)
  # pre-conditioning
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    min.rows = 1,
    min.cols = 2,
    null.ok = FALSE,
    .var.name = "obj"
  )
  check_finite(obj)
  nr <- nrow(obj)
  nc <- ncol(obj)
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

  # parallelization + subset-injection mode flag.
  # TRUE when the underlying imputer is knn_imp(), in which case we:
  # (a) disallow user-supplied `subset` in `parameters`, and
  # (b) inject a per-rep `subset` derived from the NA locations so that
  #   only injected-NA columns are imputed (neighbors across all cols).
  is_knn_mode <- is.character(.f) && (
    .f == "knn_imp" ||
      (.f == "slide_imp" && "k" %in% names(parameters))
  )

  if (is_knn_mode && "subset" %in% names(parameters)) {
    stop(
      "`subset` cannot be supplied in `parameters` for K-NN-based tuning.",
      call. = FALSE
    )
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
    na_col_subset = na_col_subset,
    max_attempts = max_attempts
  )
  n_reps <- length(na_loc)

  # `na_loc` is guaranteed to be a list of 2-col (row, col) integer
  # matrices (see resolve_na_loc). For K-NN-based tuning, precompute the
  # unique injected columns per rep to pass as `subset` at call time.
  subset_per_rep <- if (is_knn_mode) {
    lapply(na_loc, function(m) {
      sort.int(unique(as.integer(m[, 2L])))
    })
  } else {
    NULL
  }

  .rowid <- seq_len(nrow(parameters))

  indices <- expand.grid(
    param_set = .rowid,
    rep_id = seq_len(n_reps)
  )

  # `truth_vec` depends only on `rep_id`, so precompute once per rep.
  rep_ids <- indices$rep_id
  param_sets <- indices$param_set
  truth_list <- lapply(na_loc, function(pos) obj[pos])

  parallelize <- tryCatch(mirai::require_daemons(), error = function(e) FALSE)
  if (is_knn_mode) {
    if (cores > 1 && parallelize) {
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
        "cores = %d but is ignored for non-K-NN imputation. Call `mirai::daemons(%d)` to set up parallelization.",
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
    # run the function once to validate the output
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
    # share obj across daemons via shared memory instead of crating. Each worker
    # attaches on demand and materializes a local `pre`.
    big_obj <- bigmemory::as.big.matrix(obj, shared = TRUE)
    big_obj_desc <- bigmemory::describe(big_obj)
    on.exit(
      {
        rm(big_obj)
        gc()
      },
      add = TRUE
    )
    crated_fn <- carrier::crate(
      function(i) {
        if (pin_blas) {
          RhpcBLASctl::blas_set_num_threads(1)
          RhpcBLASctl::omp_set_num_threads(1)
        }
        tryCatch(
          {
            rep_id <- rep_ids[i]
            ps <- param_sets[i]
            na_positions <- na_loc[[rep_id]]
            truth_vec <- truth_list[[rep_id]]
            src <- bigmemory::attach.big.matrix(big_obj_desc)
            pre <- src[, ]
            rm(src)
            pre[na_positions] <- NA
            call_args <- c(
              list(obj = pre),
              fixed_args,
              parameters_list[[ps]]
            )
            if (!is.null(subset_per_rep)) {
              call_args$subset <- subset_per_rep[[rep_id]]
            }
            imputed_values <- do.call(impute_f, call_args)[na_positions]
            data.frame(truth = truth_vec, estimate = imputed_values)
          },
          error = function(e) {
            out <- data.frame(truth = numeric(), estimate = numeric())
            attr(out, "error") <- conditionMessage(e)
            out
          }
        )
      },
      impute_f = impute_f,
      big_obj_desc = big_obj_desc,
      truth_list = truth_list,
      na_loc = na_loc,
      rep_ids = rep_ids,
      param_sets = param_sets,
      parameters_list = parameters_list,
      fixed_args = fixed_args,
      subset_per_rep = subset_per_rep,
      pin_blas = pin_blas
    )
    result_list <- mirai::mirai_map(iter, crated_fn)[.progress = .progress]
  } else {
    # sequential: `indices` is built by expand.grid with param_set varying
    # fastest, so consecutive iterations share rep_id. Build `pre` once per
    # rep and reuse it across all param_sets in that rep.
    run_one <- function(pre, ps, rep_id) {
      call_args <- c(
        list(obj = pre),
        fixed_args,
        parameters_list[[ps]]
      )
      if (!is.null(subset_per_rep)) {
        call_args$subset <- subset_per_rep[[rep_id]]
      }
      imputed_result <- do.call(impute_f, call_args)
      data.frame(
        truth = truth_list[[rep_id]],
        estimate = imputed_result[na_loc[[rep_id]]]
      )
    }
    result_list <- vector("list", length(iter))
    if (.progress) pb <- cli::cli_progress_bar(total = length(iter))
    current_rep <- 0L
    pre <- NULL
    for (i in iter) {
      rep_id <- rep_ids[i]
      ps <- param_sets[i]
      if (rep_id != current_rep) {
        pre <- NULL
        pre <- obj
        pre[na_loc[[rep_id]]] <- NA
        current_rep <- rep_id
      }
      result_list[[i]] <- tryCatch(
        run_one(pre, ps, rep_id),
        error = function(e) {
          out <- data.frame(truth = numeric(), estimate = numeric())
          attr(out, "error") <- conditionMessage(e)
          out
        }
      )
      if (.progress) cli::cli_progress_update(id = pb)
    }
    if (.progress) cli::cli_progress_done(id = pb)
    rm(pre)
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
#' @param results A `slideimp_tune` data.frame from [tune_imp()]. Must
#' contain a `result` list-column with data.frames that have `truth` and
#' `estimate` columns.
#'
#' @param metrics A character vector of metric names to compute. Defaults
#' to `c("mae", "rmse")`. Also available: `"mape"`, `"bias"`, `"rsq"`, and
#' `"rsq_trad"`.
#'
#' @returns A data.frame with the original parameters along with
#'   unnested metrics: `.metric`, `.estimator`, and `.estimate`.
#'
#' @examples
#' obj <- sim_mat(100, 100)$input
#'
#' set.seed(1234)
#' results <- tune_imp(
#'   obj = obj,
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
# #' @returns A numeric vector of linear positions corresponding to the input 2D
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

# # sample_each_rep (R) ----
# # for each rep, we are sampling such that for each col we don't:
# # * Break row room budget
# # * Break col room budget
# # * Introduce zcv columns
# # only then, after collecting n_cols, do we return a rep. The row room budget
# # is more tricky, since we are iterating column wise, we have to update
# # row room budget after each selected column. Convert to arma is straightforawrd
# sample_each_rep <- function(
#   obj, # arma::mat &obj
#   pool_idx, # arma::uvec &pool_idx (convert to zero)
#   na_per_col, # arma::uvec &na_per_col
#   row_room, # arma::uvec &row_room
#   col_room, # arma::uvec &col_room
#   max_attempts # arma::uword max_attempts
# ) {
#   # arma seed is hooked to R's seed. So we won't have to seed in arma.
#   n_cols <- length(na_per_col) # total cols needed
#   total_na <- sum(na_per_col) # total na values needed. Needed to track 2D index
#
#   for (iter in seq_len(max_attempts)) {
#     # resetting state by reshuffling column order.
#     shuffled <- pool_idx[sample.int(length(pool_idx))]
#     attempt_row_room <- row_room
#     row_out <- integer(total_na)
#     col_out <- integer(total_na)
#     written <- 0L
#     slot <- 0L # each column occupy a "slot"
#     for (col_idx in shuffled) {
#       if (slot >= n_cols) {
#         # if we filled n_cols, then no longer need to scan the available pool
#         break
#       }
#       needed <- na_per_col[slot + 1L] # how many row NA is needed for this col
#       # column budget (colmax). This column don't have enough rows left for us given budget
#       if (col_room[col_idx] < needed) {
#         next
#       }
#       candidate <- obj[, col_idx] # slice column
#       obs_row <- which(!is.na(candidate)) # arma::find_finite(candidate)
#       # protect 2 distinct observed values so column stays non-ZV
#       uniq_vals <- unique(candidate[obs_row])
#       if (length(uniq_vals) < 2L) {
#         next
#       }
#       kept_vals <- sample(uniq_vals, size = 2L)
#       # match keep the first two sampled missing values. This brings `kept_global` back to
#       # to global row space. But we probably can find_unique() and sample the index in arma.
#       kept_global <- obs_row[match(kept_vals, candidate[obs_row])]
#       left_over <- setdiff(obs_row, kept_global)
#       # row budget (rowmax)
#       left_over <- left_over[attempt_row_room[left_over] > 0L]
#       if (length(left_over) < needed) {
#         next
#       }
#       # avoid sample() scalar trap
#       selected <- left_over[sample.int(length(left_over), needed)]
#       # commit
#       slot <- slot + 1L
#       idx <- (written + 1L):(written + needed)
#       row_out[idx] <- selected
#       col_out[idx] <- col_idx
#       written <- written + needed
#       attempt_row_room[selected] <- attempt_row_room[selected] - 1L
#     }
#
#     if (slot == n_cols) {
#       return(cbind(row = row_out, col = col_out))
#     }
#   }
#
#   cli::cli_abort("Failed to sample NA locations after {max_attempts} attempts.")
# }

#' Legacy Function
#'
#' @keywords internal
#' @noRd
inject_na <- function(
  obj,
  num_na = NULL,
  rowmax = 0.9,
  colmax = 0.9,
  check_sd = FALSE,
  max_iter = 1000
) {
  # subset the matrix to the specified features and samples
  na_mat <- !is.na(obj)
  # check if existing NA pattern already exceeds thresholds
  max_col_miss <- floor(nrow(na_mat) * colmax)
  max_row_miss <- floor(ncol(na_mat) * rowmax)
  # calculate current missingness
  current_col_miss <- nrow(na_mat) - colSums(na_mat)
  current_row_miss <- ncol(na_mat) - rowSums(na_mat)
  # check if any columns already exceed the threshold
  bad_cols <- which(current_col_miss > max_col_miss)
  if (length(bad_cols) > 0) {
    stop("Some columns have missing > colmax before na injections")
  }
  # check if any rows already exceed the threshold
  bad_rows <- which(current_row_miss > max_row_miss)
  if (length(bad_rows) > 0) {
    stop("Some rows have missing > rowmax before na injections")
  }
  not_na <- which(na_mat)
  # ensure 'num_na' does not exceed the number of available non-NA elements
  if (num_na > length(not_na)) {
    stop(
      sprintf(
        "'num_na' (%d) exceeds the number of available non-NA elements (%d).
        Adjust 'num_na' or increase feature/sample group size.",
        num_na,
        length(not_na)
      )
    )
  }

  if (check_sd) {
    # initial check for zero variance
    cvars <- col_vars(obj)
    if (any(cvars < .Machine$double.eps | is.na(cvars))) {
      stop("Zero variance columns detected before na injections")
    }
  }

  # initialize variables for the while loop
  c_miss <- TRUE
  r_miss <- TRUE
  sd_ok <- !check_sd # TRUE if not checking, FALSE if we need to check

  na_loc <- NULL
  iter <- 0
  # inject `NA` while ensuring missingness thresholds and iter are not exceeded
  while (c_miss || r_miss || !sd_ok) {
    iter <- iter + 1
    if (iter > max_iter) {
      stop(
        "NA injection failed. Adjust 'num_na' or increase 'colmax' and 'rowmax'."
      )
    }
    na_mat_test <- na_mat
    na_loc <- sample(not_na, size = num_na)
    na_mat_test[na_loc] <- FALSE

    # calculate the counts of missing values in columns and rows
    col_miss_count <- nrow(na_mat_test) - colSums(na_mat_test)
    row_miss_count <- ncol(na_mat_test) - rowSums(na_mat_test)
    # check if any column or row exceeds the missingness thresholds
    c_miss <- any(col_miss_count > max_col_miss)
    r_miss <- any(row_miss_count > max_row_miss)

    if (check_sd) {
      # this is the real obj, not is.na(obj)
      obj_test <- obj
      obj_test[na_loc] <- NA
      cvars <- col_vars(obj_test)
      sd_ok <- !any(cvars < .Machine$double.eps | is.na(cvars))
    }
  }
  return(na_loc)
}
