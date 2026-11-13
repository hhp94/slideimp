#' Sample Missing-Value Locations with Constraints
#'
#' Sample matrix indices for `NA` injection while respecting row and column
#' missingness limits and avoiding zero-variance columns.
#'
#' @details
#' The function uses a greedy stochastic search for valid `NA` locations. It
#' ensures that:
#'
#' - Total missingness per row and column does not exceed `rowmax` and
#'   `colmax`.
#' - At least two distinct observed values are preserved in every affected
#'   column.
#' - A column with a variance below `.Machine$double.eps` counts as constant
#'   and is never selected. [pca_imp()] uses the same threshold to leave a
#'   column out of its fit.
#'
#' A column that cannot take `n_rows` missing values under these limits is
#' never selected. Such a column does not stop the sampling.
#'
#' When both `n_cols` and `num_na` are `NULL`, `num_na` defaults to 5% of the
#' cells of `obj`, with a maximum of 500. The default never asks for more
#' columns than the number that can take `n_rows` missing values.
#'
#' @param obj A numeric matrix.
#' @param n_cols Integer or `NULL`. Number of columns to receive injected
#'   missing values.
#' @param n_rows Integer. Target number of missing values to inject per selected
#'   column.
#' @param num_na Integer or `NULL`. Total number of missing values to inject per
#'   repetition. If supplied, `n_cols` is derived from `num_na` and `n_rows`,
#'   and missing values are distributed as evenly as possible across columns.
#' @param n_reps Integer. Number of independent repetitions.
#' @param rowmax Numeric scalar between `0` and `1`. Maximum allowed
#'   missing-data proportion per row after injection.
#' @param colmax Numeric scalar between `0` and `1`. Maximum allowed
#'   missing-data proportion per column after injection.
#' @param na_col_subset Optional integer or character vector restricting which
#'   columns are eligible for missing-value injection.
#' @param max_attempts Integer. Maximum number of resampling attempts per
#'   repetition before giving up.
#'
#' @returns A list of length `n_reps`. Each element is a two-column integer
#' matrix with row and column indices for sampled `NA` locations.
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(runif(100), nrow = 10)
#'
#' # Sample 5 missing values across 5 columns
#' locs <- sample_na_loc(mat, n_cols = 5, n_rows = 1)
#' locs
#'
#' # Inject the missing values from the first repetition
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
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    min.rows = 1,
    min.cols = 1,
    .var.name = "obj"
  )
  # Inf makes col_vars() NaN, which would read as a constant column
  check_inf(obj)
  checkmate::assert_count(
    n_cols,
    positive = TRUE,
    null.ok = TRUE,
    .var.name = "n_cols"
  )
  checkmate::assert_count(n_rows, positive = TRUE, .var.name = "n_rows")
  checkmate::assert_count(
    num_na,
    positive = TRUE,
    null.ok = TRUE,
    .var.name = "num_na"
  )
  checkmate::assert_count(n_reps, positive = TRUE, .var.name = "n_reps")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_count(
    max_attempts,
    positive = TRUE,
    .var.name = "max_attempts"
  )

  n <- nrow(obj)
  # the frame of this call, for errors raised inside a closure below
  fn_env <- environment()
  # with neither given, `num_na` gets a default further down, once the number
  # of columns that can take `n_rows` NAs is known
  auto_num_na <- is.null(n_cols) && is.null(num_na)
  num_na_supplied <- !is.null(num_na)

  if (num_na_supplied && num_na < n_rows) {
    cli::cli_abort(c(
      "{.arg num_na} ({num_na}) must be at least {.arg n_rows} ({n_rows}).",
      "i" = "Each selected column receives at least {.arg n_rows} NAs."
    ))
  }

  # resolve (n_cols, na_per_col)
  if (num_na_supplied) {
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
  } else if (!auto_num_na) {
    na_per_col <- rep.int(as.integer(n_rows), n_cols)
    num_na <- n_rows * n_cols
  }
  # a column keeps 2 observed values, so it cannot take more than `n - 2` NAs.
  # The default `num_na` fills whole columns of `n_rows`.
  max_per_col <- if (auto_num_na) n_rows else max(na_per_col)
  if (max_per_col > n - 2L) {
    cli::cli_abort(
      c(
        "Too many NAs requested per column.",
        "x" = "The request puts {max_per_col} NA{?s} in one column, but {.arg obj} has {n} row{?s}.",
        "i" = "Each column must keep at least 2 observed values. Decrease {.arg n_rows} or {.arg num_na}."
      ),
      class = "slideimp_infeasible"
    )
  }

  # resolve `na_col_subset` into integer pool_idx
  pool_idx <- if (is.null(na_col_subset)) {
    seq_len(ncol(obj))
  } else if (is.character(na_col_subset)) {
    checkmate::assert_character(
      na_col_subset,
      any.missing = FALSE,
      min.len = 1,
      unique = TRUE,
      .var.name = "na_col_subset"
    )
    if (is.null(colnames(obj))) {
      cli::cli_abort(
        "{.arg na_col_subset} is character but {.arg obj} has no colnames."
      )
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
      lower = 1,
      upper = ncol(obj),
      any.missing = FALSE,
      min.len = 1,
      unique = TRUE,
      .var.name = "na_col_subset"
    )
    # round(): an integerish value can sit just under its integer
    as.integer(round(na_col_subset))
  } else {
    cli::cli_abort("{.arg na_col_subset} must be numeric, character, or NULL.")
  }

  n_pool <- length(pool_idx)
  if (!auto_num_na && n_cols > n_pool) {
    lead <- if (num_na_supplied) {
      "Cannot place {num_na} NAs (requires {n_cols} columns given {.arg n_rows} = {n_rows})."
    } else {
      "Cannot place NAs in {n_cols} columns."
    }
    pool_size <- if (is.null(na_col_subset)) {
      "{.arg obj} has only {n_pool} column{?s}."
    } else {
      "{.arg na_col_subset} provides only {n_pool} column{?s}."
    }
    cli::cli_abort(c(lead, "i" = pool_size), class = "slideimp_infeasible")
  }

  # pre-injection state. mat_miss() counts without an n x p logical mask.
  col_miss <- unname(mat_miss(obj, col = TRUE))
  row_miss <- unname(mat_miss(obj, col = FALSE))
  max_col_miss <- max_miss_count(n, colmax)
  max_row_miss <- max_miss_count(ncol(obj), rowmax)

  if (any(row_miss > max_row_miss)) {
    bad <- which(row_miss > max_row_miss)
    cli::cli_abort(c(
      "Some rows already exceed {.arg rowmax} before injection.",
      "x" = "Problematic rows: {fmt_trunc(bad, 6)}"
    ))
  }
  # never negative after the check above. It crosses into C++ as unsigned.
  row_room <- max_row_miss - row_miss

  # The most NAs each pool column can take: its room under `colmax`, less the
  # 2 observed values it has to keep. A column over `colmax`, with fewer than
  # 2 observed values, or constant can take none. Such a column is left out
  # of the sampling. It does not stop the call: only the columns that receive
  # NAs have to qualify, and knn_imp() and pca_imp() apply their own `colmax`
  # to the rest.
  # "Constant" is pca_imp()'s eligibility test, variance below double.eps,
  # not variance zero. pca_imp() refuses a call in which every column with NAs
  # is ineligible, so with .f = "pca_imp" a tune_imp() repetition that blanked
  # only near-constant columns would fail.
  pool_miss <- col_miss[pool_idx]
  pool_var <- unname(col_vars(obj))[pool_idx]
  pool_cap <- pmin(max_col_miss - pool_miss, n - pool_miss - 2L)
  pool_cap[
    is.na(pool_var) | pool_var < .Machine$double.eps | pool_cap < 0L
  ] <- 0L

  if (auto_num_na) {
    num_na <- default_num_na(length(obj), n_rows, sum(pool_cap >= n_rows))
    pct <- signif(100 * num_na / length(obj), 2)
    cli::cli_inform(c(
      "Using default {.arg num_na} = {num_na} (~{pct}% of cells).",
      "i" = "Increase it for more reliable results, or decrease it if missing values are dense."
    ))
    n_cols <- num_na %/% n_rows
    na_per_col <- rep.int(as.integer(n_rows), n_cols)
  }

  # Feasibility, exact at the column level: the i-th largest demand needs i
  # columns that can take it. Only the interplay with the row budgets is left
  # to the stochastic search in C++.
  fits <- pool_cap >= na_per_col[1L] # `na_per_col` is ascending
  demand <- rev(na_per_col)
  cap_at_rank <- sort(pool_cap[fits], decreasing = TRUE)[seq_len(n_cols)]
  cap_at_rank[is.na(cap_at_rank)] <- 0L
  short <- which(demand > cap_at_rank)
  if (length(short) > 0L) {
    n_need <- short[1L]
    need <- demand[n_need]
    n_can <- sum(pool_cap >= need)
    cli::cli_abort(
      c(
        "The requested NAs do not fit in the available columns.",
        "x" = "{n_need} column{?s} must take at least {need} NA{?s} each, but only {n_can} column{?s} can.",
        "i" = "Decrease {.arg num_na}, {.arg n_cols} or {.arg n_rows}, or increase {.arg colmax}."
      ),
      class = "slideimp_infeasible"
    )
  }

  pool_fit0 <- pool_idx[fits] - 1L
  cap_fit <- pool_cap[fits]
  lapply(seq_len(n_reps), function(rep_id) {
    result <- sample_each_rep_cpp(
      obj = obj,
      pool_idx_in = pool_fit0,
      na_per_col = na_per_col,
      row_room = row_room,
      pool_cap = cap_fit,
      max_attempts = max_attempts
    )
    # zero rows: no attempt filled every column
    if (nrow(result) == 0L) {
      cli::cli_abort(
        c(
          "No valid NA locations were found for repetition {rep_id} in {max_attempts} attempt{?s}.",
          "i" = "Decrease {.arg num_na}, {.arg n_cols} or {.arg n_rows}, or increase {.arg rowmax}, {.arg colmax} or {.arg max_attempts}."
        ),
        class = "slideimp_infeasible",
        call = fn_env
      )
    }
    # arma::umat is a double matrix under 64-bit indices
    storage.mode(result) <- "integer"
    colnames(result) <- c("row", "col")
    result
  })
}

#' Largest Missing Count Allowed by a Proportion Cap
#'
#' The largest count `k` out of `n` with `k / n <= prop`. This is the
#' comparison `knn_imp()` and `pca_imp()` make against `colmax`, so a column
#' they admit is never over the cap here. `floor(n * prop)` alone can land one
#' short, because the product can fall just under an integer.
#'
#' @noRd
#' @keywords internal
max_miss_count <- function(n, prop) {
  k <- floor(n * prop)
  k <- k + ((k + 1) / n <= prop) - (k / n > prop)
  as.integer(k)
}

#' Default `num_na` for `sample_na_loc()`
#'
#' Used when neither `num_na` nor `n_cols` is given: 5% of the cells, at most
#' 500, in whole columns of `n_rows`, and never more columns than `n_fit`, the
#' number of columns that can take `n_rows` NAs. `n_total` is a double, since
#' the cell count of a large cohort passes `.Machine$integer.max`.
#'
#' @noRd
#' @keywords internal
default_num_na <- function(n_total, n_rows, n_fit) {
  target <- min(500, round(0.05 * n_total))
  n_cols <- min(max(1, target %/% n_rows), max(1, n_fit))
  as.integer(n_cols * n_rows)
}

#' Resolve NA Locations for `tune_imp()`
#'
#' If `na_loc` is NULL, generates `n_reps` random NA location matrices via
#' `sample_na_loc()`. Otherwise normalizes the user-supplied
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
  na_col_subset,
  max_attempts
) {
  if (is.null(na_loc)) {
    # with `num_na` and `n_cols` both NULL, sample_na_loc() picks the default
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
      cli::cli_abort(
        "{.arg na_loc} must be a 2-col matrix, an integer vector, or a list of these."
      )
    }
  }

  checkmate::assert_list(na_loc, min.len = 1, .var.name = "na_loc")

  nr <- nrow(obj)
  nc <- ncol(obj)
  n_total <- length(obj)

  # validate each element, then put it in the canonical form: a 2-column
  # integer matrix of (row, col). round() before the integer conversion: an
  # integerish value can sit just under its integer (`0.29 * 100`), and
  # as.integer() alone truncates it to the cell before.
  for (i in seq_along(na_loc)) {
    elem <- na_loc[[i]]
    nm <- sprintf("na_loc[[%d]]", i)
    if (is.matrix(elem)) {
      checkmate::assert_matrix(
        elem,
        mode = "integerish",
        any.missing = FALSE,
        ncols = 2,
        min.rows = 1,
        .var.name = nm
      )
      # columns are positional, except for a matrix that names them the
      # other way round
      if (identical(colnames(elem), c("col", "row"))) {
        elem <- elem[, c("row", "col"), drop = FALSE]
      }
      checkmate::assert_true(
        all(elem[, 1] >= 1 & elem[, 1] <= nr),
        .var.name = paste0(nm, " row indices")
      )
      checkmate::assert_true(
        all(elem[, 2] >= 1 & elem[, 2] <= nc),
        .var.name = paste0(nm, " col indices")
      )
      elem <- round(elem)
      storage.mode(elem) <- "integer"
      colnames(elem) <- c("row", "col")
    } else {
      checkmate::assert_integerish(
        elem,
        lower = 1,
        upper = n_total,
        any.missing = FALSE,
        min.len = 1,
        unique = TRUE,
        .var.name = nm
      )
      lp <- as.integer(round(elem))
      elem <- cbind(
        row = ((lp - 1L) %% nr) + 1L,
        col = ((lp - 1L) %/% nr) + 1L
      )
    }
    na_loc[[i]] <- elem
  }

  elem_sizes <- vapply(na_loc, nrow, integer(1))
  if (length(unique(elem_sizes)) != 1L) {
    cli::cli_abort(
      "All elements in {.arg na_loc} must specify the same number of NA positions."
    )
  }

  # a repeated cell would be scored once per repeat
  for (i in seq_along(na_loc)) {
    # as linear positions, in double: one hash pass, and no overflow
    elem <- na_loc[[i]]
    dup <- which(duplicated((elem[, 2L] - 1) * nr + elem[, 1L]))
    if (length(dup) > 0L) {
      cli::cli_abort(c(
        "{.arg na_loc} positions must be unique.",
        "x" = "{.code na_loc[[{i}]]} has {length(dup)} repeated position{?s}: {fmt_trunc(dup, 6)}."
      ))
    }
  }

  for (i in seq_along(na_loc)) {
    on_hole <- is.na(obj[na_loc[[i]]])
    if (any(on_hole)) {
      bad <- which(on_hole)
      cli::cli_abort(c(
        "{.arg na_loc} positions must point at observed cells.",
        "x" = "{.code na_loc[[{i}]]} has {length(bad)} position{?s} on a cell
               that is already NA or NaN: {fmt_trunc(bad, 6)}.",
        "i" = "Use observed cells only, or omit {.arg na_loc} to sample them."
      ))
    }
  }

  return(na_loc)
}

#' Tune Imputation Method Parameters
#'
#' Tune method specific hyperparameters by repeatedly masking observed values,
#' imputing them, and comparing the imputed values with the original values.
#'
#' @param parameters A `data.frame` specifying parameter combinations to tune.
#'   Each column should be a parameter accepted by `.f`, excluding `obj`.
#'   List-columns are supported for complex parameters. Duplicate rows are
#'   removed. `NULL` is treated as a single parameter set with no additional
#'   arguments, which is useful for functions whose required arguments all have
#'   defaults.
#' @param .f One of `"knn_imp"`, `"pca_imp"`, or `"slide_imp"`, or a custom
#'   imputation function.
#' @param na_loc Optional predefined missing-value locations. Accepted formats
#'   are a two-column integer matrix of row and column indices, a numeric vector
#'   of linear positions, or a list whose elements are either of those formats.
#' @param .progress Logical. If `TRUE`, show progress during tuning.
#' @param cores Integer. Number of cores to use for K-NN and sliding-window
#'   K-NN imputation. For other methods, use `mirai::daemons()`.
#' @param location Numeric vector of column locations. Required when
#'   `.f = "slide_imp"`.
#' @param pin_blas Logical. If `TRUE`, pin BLAS threads to 1 on the `mirai`
#'   daemons to reduce thread contention.
#' @inheritParams sample_na_loc
#'
#' @details
#' Built-in methods can be selected by passing `.f = "knn_imp"`,
#' `.f = "pca_imp"`, or `.f = "slide_imp"`. A custom function can also be
#' supplied. Custom functions must accept `obj` as their first argument and
#' return a numeric matrix with the same dimensions as `obj`. The return value
#' is checked at every iteration. An iteration that returns something else is
#' recorded in the `error` column.
#'
#' The columns of `parameters` for a custom function must not use the names
#' of the columns that `tune_imp()` and [compute_metrics()] add: `param_set`,
#' `rep_id`, `result`, `error`, `n`, `n_miss`, `.metric`, `.estimator`, and
#' `.estimate`.
#'
#' When `mirai` daemons are set, a custom function runs on the daemons. There
#' it finds the packages that are attached on the daemons. It does not find
#' variables from the session that defined it. Pass such values through
#' `parameters`.
#'
#' When `na_loc` is supplied, `num_na`, `n_reps`, `n_cols`, `n_rows`, `rowmax`,
#' `colmax`, `na_col_subset`, and `max_attempts` are ignored. The supplied
#' positions are used as given. They are not checked against `rowmax` or
#' `colmax`.
#'
#' When `.f` is a character string, columns in `parameters` are validated
#' against the selected method:
#'
#' - `"knn_imp"` requires `k`.
#' - `"pca_imp"` requires `ncp`.
#' - `"slide_imp"` requires `window_size`, `overlap_size`, and `min_window_n`,
#'   plus exactly one of `k` or `ncp`.
#'
#' To tune parameters for grouped imputation, tune [knn_imp()] or [pca_imp()]
#' on representative groups, then pass the selected parameters to [group_imp()].
#'
#' The top-level `rowmax` and `colmax` arguments control random missing-value
#' injection performed by [sample_na_loc()]. To tune or pass an imputation
#' method's own `colmax` argument, include a `colmax` column in `parameters`.
#'
#' Tuning results can be summarized with [compute_metrics()] or evaluated with
#' external packages such as `yardstick`.
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @inheritSection group_imp Parallelization
#'
#' @inheritSection pca_imp PCA Performance tips
#'
#' @returns A data frame of class `slideimp_tune` containing:
#'   - columns originally provided in `parameters`;
#'   - `param_set`, an integer ID for each unique parameter combination;
#'   - `rep_id`, an integer repetition index;
#'   - `result`, a list-column where each element is a data frame containing
#'     `truth` and `estimate` columns;
#'   - `error`, a character column containing the error message if the
#'     iteration failed, otherwise `NA`.
#'
#' @examples
#' set.seed(123)
#'
#' # Simulate some data
#' obj <- sim_mat(10, 50)$input
#'
#' # Tune K-NN imputation with random missing-value injection.
#' # Use larger `num_na` and `n_reps` values for real analyses.
#' params_knn <- data.frame(k = c(2, 4))
#' results <- tune_imp(
#'   obj,
#'   params_knn,
#'   .f = "knn_imp",
#'   n_reps = 1,
#'   num_na = 10,
#'   .progress = FALSE
#' )
#' compute_metrics(results)
#'
#' # Tune with fixed missing-value positions
#' na_positions <- list(
#'   matrix(c(1, 2, 3, 1, 1, 1), ncol = 2),
#'   matrix(c(2, 3, 4, 2, 2, 2), ncol = 2)
#' )
#'
#' results_fixed <- tune_imp(
#'   obj,
#'   data.frame(k = 2),
#'   .f = "knn_imp",
#'   na_loc = na_positions,
#'   .progress = FALSE
#' )
#'
#' # Custom imputation function
#' custom_fill <- function(obj, val = 0) {
#'   obj[is.na(obj)] <- val
#'   obj
#' }
#'
#' tune_imp(
#'   obj,
#'   data.frame(val = c(0, 1)),
#'   .f = custom_fill,
#'   num_na = 10,
#'   .progress = FALSE
#' )
#'
#' @examplesIf interactive() && requireNamespace("mirai", quietly = TRUE)
#' # Parallel tuning with mirai
#' mirai::daemons(2)
#'
#' parameters_custom <- data.frame(mean = c(0, 1), sd = c(1, 1))
#'
#' custom_imp <- function(obj, mean, sd) {
#'   na_pos <- is.na(obj)
#'   obj[na_pos] <- stats::rnorm(sum(na_pos), mean = mean, sd = sd)
#'   obj
#' }
#'
#' results_p <- tune_imp(
#'   obj,
#'   parameters_custom,
#'   .f = custom_imp,
#'   n_reps = 1,
#'   num_na = 10,
#'   .progress = FALSE
#' )
#'
#' mirai::daemons(0)
#'
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
  # only gc on longer running calls
  .start <- proc.time()[["elapsed"]]
  on.exit(
    {
      if (proc.time()[["elapsed"]] - .start > 0.5) {
        gc(verbose = FALSE)
      }
    },
    add = TRUE
  )

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
    if (!(length(.f) == 1 && .f %in% c("slide_imp", "knn_imp", "pca_imp"))) {
      cli::cli_abort(
        "{.arg .f} must be one of {.val slide_imp}, {.val knn_imp}, or {.val pca_imp}."
      )
    }
    fn_name <- .f
    # Args excluded from tuning across all .f:
    #   obj       - supplied by tune_imp
    #   subset    - injected by tune_imp from na_loc (for slide_imp/knn_imp)
    #   cores     - passed via tune_imp's `cores` arg, not tuneable
    #   .progress - tune_imp owns progress reporting
    # slide_imp-specific exclusions:
    #   location  - passed via tune_imp's `location` arg
    #   dry_run   - returns a data.frame, breaks matrix contract
    if (!parameters_is_null) {
      # the other method's tuning parameter is the likely mix-up, so it gets
      # its own message before the general check on parameter names
      if (.f == "knn_imp") {
        if ("ncp" %in% names(parameters)) {
          cli::cli_abort(
            "{.fn knn_imp} does not accept {.arg ncp}. Did you mean {.code .f = 'pca_imp'}?"
          )
        }
        check_unknown_params(
          param_names = names(parameters),
          mode = "knn",
          arg = "parameters"
        )
      } else if (.f == "pca_imp") {
        if ("k" %in% names(parameters)) {
          cli::cli_abort(
            "{.fn pca_imp} does not accept {.arg k}. Did you mean {.code .f = 'knn_imp'}?"
          )
        }
        check_unknown_params(
          param_names = names(parameters),
          mode = "pca",
          arg = "parameters"
        )
      } else if (.f == "slide_imp") {
        allowed_params <- slide_tune_params()
        unknown_params <- setdiff(names(parameters), allowed_params)

        if (length(unknown_params) > 0) {
          cli::cli_abort(c(
            "{cli::qty(length(unknown_params))}Unknown parameter{?s} in {.arg parameters} for {.fn {fn_name}}:",
            "x" = "{fmt_trunc(unknown_params, 10)}",
            "i" = "Allowed: {.arg {allowed_params}}."
          ))
        }
      }
    }

    if (.f == "slide_imp") {
      if (is.null(location)) {
        cli::cli_abort(
          "Tuning {.fn slide_imp} requires the {.arg location} argument provided to `tune_imp(location = ...)`."
        )
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
      if (
        !all(
          c("window_size", "overlap_size", "min_window_n") %in%
            names(parameters)
        )
      ) {
        cli::cli_abort(
          "{.fn slide_imp} requires {.arg window_size}, {.arg overlap_size}, and {.arg min_window_n} in {.arg parameters}."
        )
      }
      if ("k" %in% names(parameters) && "ncp" %in% names(parameters)) {
        cli::cli_abort(
          "{.fn slide_imp} requires exactly one of {.arg k} or {.arg ncp} in {.arg parameters}, not both."
        )
      }
      if (!any(c("k", "ncp") %in% names(parameters))) {
        cli::cli_abort(
          "{.fn slide_imp} requires either {.arg k} or {.arg ncp} in {.arg parameters}."
        )
      }
    } else if (.f == "knn_imp") {
      if (!"k" %in% names(parameters)) {
        cli::cli_abort("{.fn knn_imp} requires {.arg k} in {.arg parameters}.")
      }
    } else if (.f == "pca_imp") {
      if (!"ncp" %in% names(parameters)) {
        cli::cli_abort(
          "{.fn pca_imp} requires {.arg ncp} in {.arg parameters}."
        )
      }
    }

    # The methods draw their own progress bar by default and `tune_imp()` owns
    # progress here. `.progress = FALSE` is passed with the fixed arguments, so
    # a `.progress` column cannot stay: the call would get the argument twice.
    if (".progress" %in% names(parameters)) {
      cli::cli_alert_info(
        "{.arg .progress} in {.arg parameters} is suppressed."
      )
      parameters$.progress <- NULL
      parameters <- unique(parameters)
    }

    cli::cli_inform("Tuning {.fn {fn_name}}")
  } else if (is.function(.f)) {
    # a parameter with one of these names would be overwritten, or would
    # leave two columns of the same name in the output
    clash <- intersect(
      names(parameters),
      c(
        "param_set",
        "rep_id",
        "result",
        "error",
        "n",
        "n_miss",
        ".metric",
        ".estimator",
        ".estimate"
      )
    )
    if (length(clash) > 0) {
      cli::cli_abort(c(
        "{.arg parameters} uses column names that the output of {.fn tune_imp} or {.fn compute_metrics} needs.",
        "x" = "{cli::qty(length(clash))}Name{?s} in use: {.val {clash}}.",
        "i" = "{cli::qty(length(clash))}Rename the column{?s} and the matching argument{?s} of {.arg .f}."
      ))
    }
    cli::cli_inform("Tuning {.field custom function}")
  } else {
    cli::cli_abort(
      "{.arg .f} must be a function or one of {.val slide_imp}, {.val knn_imp}, {.val pca_imp}."
    )
  }

  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_integerish(
    cores,
    lower = 1,
    len = 1,
    any.missing = FALSE,
    null.ok = FALSE,
    .var.name = "cores"
  )

  # Two separate flags:
  # - is_subset_mode: slide_imp (both PCA & KNN) + knn_imp -> controls subset
  # injection by infering automatically from `na_loc`
  # - is_knn_mode: only paths that benefit from RcppThread cores
  is_subset_mode <- is.character(.f) && .f %in% c("slide_imp", "knn_imp")
  is_knn_mode <- is.character(.f) &&
    (.f == "knn_imp" ||
      (.f == "slide_imp" && "k" %in% names(parameters)))

  cli::cli_inform("Step 1/2: {.field Resolving NA locations}")

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

  # Pre-compute subset per rep only when needed (now includes all slide_imp)
  subset_per_rep <- if (is_subset_mode) {
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
  rep_ids <- indices$rep_id
  param_sets <- indices$param_set
  truth_list <- lapply(na_loc, function(pos) obj[pos])

  parallelize <- mirai_daemons_set()

  if (is_knn_mode) {
    if (cores > 1 && parallelize) {
      cli::cli_inform(c(
        "!" = "Both {.arg cores} > 1 and active {.pkg mirai} daemons detected.",
        "i" = "Setting {.arg cores} = 1 to avoid nested parallelism. Distribution will be handled by {.pkg mirai}."
      ))
      cores <- 1
    }
  }

  if (!parallelize && cores > 1 && !is_knn_mode) {
    cli::cli_inform(c(
      "!" = "{.arg cores} = {cores} is ignored for non-K-NN imputation.",
      "i" = "Call {.code mirai::daemons({cores})} to enable parallel execution."
    ))
  }

  if (parallelize) {
    cli::cli_inform("Running mode: {.strong mirai}")
  } else if (is_knn_mode && cores > 1) {
    cli::cli_inform("Running mode: {.strong threaded} ({cores} cores)")
  } else {
    cli::cli_inform("Running mode: {.strong sequential}")
  }

  is_slide <- is.character(.f) && .f == "slide_imp"
  is_custom <- is.function(.f)
  fixed_args <- list()
  if (!is_custom) {
    # only the masked cells are read back, so the methods' own progress bar
    # is turned off
    fixed_args$.progress <- FALSE
  }
  if (is_slide) {
    fixed_args$location <- location
  }
  if (is_knn_mode) {
    fixed_args$cores <- cores
  }

  parameters_list <- lapply(
    split(parameters, f = as.factor(.rowid)),
    function(row) {
      row_list <- as.list(row)
      if (parameters_is_null) {
        row_list$.placeholder <- NULL
      }
      row_list <- lapply(row_list, function(x) {
        # expand.grid() makes a character column a factor, and match.arg()
        # in the methods rejects a factor
        if (is.factor(x)) {
          return(as.character(x))
        }
        if (is.list(x) && length(x) == 1) {
          x[[1]]
        } else {
          x
        }
      })
      row_list
    }
  )

  if (is.character(.f)) {
    impute_f <- switch(
      .f,
      slide_imp = slide_imp,
      knn_imp = knn_imp,
      pca_imp = pca_imp
    )
  } else {
    # run the function once to validate the output
    impute_f <- .f
    pre <- obj
    na_positions <- na_loc[[1]]
    pre[na_positions] <- NA
    param_vec <- parameters_list[[1]]
    tryCatch(
      {
        probe_result <- do.call(impute_f, c(list(obj = pre), param_vec))
        checkmate::assert_matrix(
          probe_result,
          mode = "numeric",
          nrows = nr,
          ncols = nc,
          null.ok = FALSE,
          .var.name = "imputed_result"
        )
        checkmate::assert_true(
          sum(is.infinite(probe_result)) == 0,
          .var.name = "imputed_result"
        )
      },
      error = function(e) {
        cli::cli_abort(
          c(
            "Custom {.arg .f} must accept {.arg obj} as its first argument and return a numeric matrix with the same dimensions as {.arg obj} and no {.val Inf} values.",
            "i" = "Original error: {conditionMessage(e)}"
          )
        )
      }
    )
  }

  cli::cli_inform("Step 2/2: {.field Tuning}")

  iter <- seq_len(nrow(indices))

  # What is wrong with the return of a custom `.f`, or NULL when it can be
  # scored. The probe above sees parameter set 1 only, so every iteration
  # asks: a matrix of another shape still answers a `[row, col]` lookup, and
  # its values would be scored as if they were the imputed cells. Base R
  # only, and no free variable, so the copy sent to a daemon needs nothing
  # from this package.
  result_problem <- function(x, nr, nc, na_positions) {
    if (!is.matrix(x) || !is.numeric(x)) {
      return("`.f` did not return a numeric matrix.")
    }
    if (nrow(x) != nr || ncol(x) != nc) {
      return(sprintf(
        "`.f` returned a %d x %d matrix, but `obj` is %d x %d.",
        nrow(x),
        ncol(x),
        nr,
        nc
      ))
    }
    if (any(is.infinite(x[na_positions]))) {
      return("`.f` returned `Inf` for a masked cell.")
    }
    NULL
  }

  n_task_failed <- 0L
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
    crated_fn <- .crate(
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
            pre <- src[,]
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
            imputed_result <- do.call(impute_f, call_args)
            if (is_custom) {
              problem <- result_problem(imputed_result, nr, nc, na_positions)
              if (!is.null(problem)) {
                stop(problem, call. = FALSE)
              }
            }
            data.frame(
              truth = truth_vec,
              estimate = imputed_result[na_positions]
            )
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
      pin_blas = pin_blas,
      is_custom = is_custom,
      result_problem = result_problem,
      nr = nr,
      nc = nc,
      # The crate re-parents a custom `.f`, so that the frame that defined it
      # (which may hold the whole matrix) is not sent with every task. With
      # `baseenv()` as the parent it would then see base R and nothing else:
      # not `stats`, not a package attached on the daemons. The global
      # environment is sent by reference, so this parent resolves on the
      # daemon, through the daemon's own search path.
      .parent_env = if (is_custom) globalenv() else baseenv()
    )
    result_list <- collect_mirai_map(
      mirai::mirai_map(iter, crated_fn),
      .progress
    )
    # A task that fails outside its own `tryCatch()` comes back as a mirai
    # error value in place of the data frame: its daemon stopped, or the task
    # could not be read there. Record it like any other failed iteration.
    task_failed <- which(vapply(result_list, mirai::is_error_value, logical(1)))
    n_task_failed <- length(task_failed)
    for (i in task_failed) {
      out <- data.frame(truth = numeric(), estimate = numeric())
      attr(out, "error") <- paste0(
        "The mirai task failed: ",
        mirai_error_text(result_list[[i]])
      )
      result_list[[i]] <- out
    }
  } else {
    # sequential: `indices` is built by expand.grid with param_set varying
    # fastest, so consecutive iterations share rep_id. Build `pre` once per
    # rep and reuse it across all param_sets in that rep.
    run_sequential <- function(pre, ps, rep_id) {
      call_args <- c(
        list(obj = pre),
        fixed_args,
        parameters_list[[ps]]
      )
      if (!is.null(subset_per_rep)) {
        call_args$subset <- subset_per_rep[[rep_id]]
      }
      imputed_result <- do.call(impute_f, call_args)
      if (is_custom) {
        problem <- result_problem(imputed_result, nr, nc, na_loc[[rep_id]])
        if (!is.null(problem)) {
          stop(problem, call. = FALSE)
        }
      }
      data.frame(
        truth = truth_list[[rep_id]],
        estimate = imputed_result[na_loc[[rep_id]]]
      )
    }
    result_list <- vector("list", length(iter))
    if (.progress) {
      pb <- cli::cli_progress_bar(name = "Tuning", total = length(iter))
    }
    current_rep <- 0L
    pre <- NULL
    for (i in iter) {
      rep_id <- rep_ids[i]
      ps <- param_sets[i]
      if (rep_id != current_rep) {
        pre <- obj
        pre[na_loc[[rep_id]]] <- NA
        current_rep <- rep_id
      }
      result_list[[i]] <- tryCatch(
        run_sequential(pre, ps, rep_id),
        error = function(e) {
          out <- data.frame(truth = numeric(), estimate = numeric())
          attr(out, "error") <- conditionMessage(e)
          out
        }
      )
      if (.progress) cli::cli_progress_update(id = pb)
    }
    if (.progress) cli::cli_progress_done(id = pb)
  }

  error_vec <- vapply(
    result_list,
    function(x) {
      e <- attr(x, "error")
      if (is.null(e)) {
        NA_character_
      } else {
        e
      }
    },
    character(1)
  )

  result_df <- cbind(
    parameters[indices$param_set, , drop = FALSE],
    indices,
    data.frame(result = I(result_list), error = error_vec)
  )

  if (any(!is.na(result_df$error))) {
    cli::cli_warn(c(
      "Some tuning iterations may have failed. Check the {.field error} column in the output.",
      if (n_task_failed > 0L) {
        c(
          "x" = "{n_task_failed} iteration{?s} did not return from the {.pkg mirai} daemons. A daemon may have stopped."
        )
      }
    ))
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
  truth_valid <- truth[valid]
  estimate_valid <- estimate[valid]
  # a correlation with a constant is undefined, and cor() warns about it
  # once per row of the tuning output
  if (
    all(truth_valid == truth_valid[1L]) ||
      all(estimate_valid == estimate_valid[1L])
  ) {
    return(NA_real_)
  }
  stats::cor(truth_valid, estimate_valid)^2
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
  # Drop the cells with no estimate, and only those. `na.rm = TRUE` on the
  # ratio would also drop a 0 / 0 cell, and the mean would then be taken over
  # fewer cells than `n - n_miss` reports. A true value of 0 leaves the result
  # not finite: Inf, or NaN when the estimate is 0 too.
  valid <- !is.na(truth) & !is.na(estimate)
  mean(abs((truth[valid] - estimate[valid]) / truth[valid])) * 100
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
#' Compute prediction accuracy metrics for results from [tune_imp()].
#'
#' @param results A `slideimp_tune` data frame from [tune_imp()]. Must contain
#'   a `result` list-column whose elements are data frames with `truth` and
#'   `estimate` columns.
#' @param metrics Character vector of metric names to compute. Defaults to
#'   `c("mae", "rmse")`. Available metrics are `"mae"`, `"rmse"`, `"mape"`,
#'   `"bias"`, `"rsq"`, and `"rsq_trad"`.
#'
#' @details
#' `"mape"` divides by the true value. It is not finite (`Inf` or `NaN`) for a
#' repetition in which a masked value is exactly `0`. Use another metric for
#' data that contain zeros.
#'
#' `"rsq"` is `NA` when the true values or the estimates of a repetition are
#' constant.
#'
#' @returns A data frame containing the original parameter columns along with
#' unnested metric columns: `.metric`, `.estimator`, and `.estimate`.
#'
#' @examples
#' set.seed(1234)
#' obj <- sim_mat(20, 30)$input
#'
#' results <- tune_imp(
#'   obj = obj,
#'   parameters = data.frame(k = 5),
#'   .f = "knn_imp",
#'   n_reps = 1,
#'   num_na = 10,
#'   .progress = FALSE
#' )
#'
#' compute_metrics(results)
#'
#' @export
compute_metrics <- function(results, metrics = c("mae", "rmse")) {
  UseMethod("compute_metrics")
}

#' @rdname compute_metrics
#' @method compute_metrics data.frame
#' @export
compute_metrics.data.frame <- function(results, metrics = c("mae", "rmse")) {
  if (!"result" %in% names(results)) {
    cli::cli_abort("{.arg results} must contain a {.field result} column.")
  }
  first_result <- results$result[[1]]
  if (
    !is.data.frame(first_result) ||
      !all(c("truth", "estimate") %in% names(first_result))
  ) {
    cli::cli_abort(
      "Each element of the {.field result} column must be a data frame with {.field truth} and {.field estimate} columns."
    )
  }
  compute_metrics.slideimp_tune(results, metrics = metrics)
}

#' @rdname compute_metrics
#' @method compute_metrics slideimp_tune
#' @export
compute_metrics.slideimp_tune <- function(results, metrics = c("mae", "rmse")) {
  checkmate::assert_character(
    metrics,
    min.len = 1,
    any.missing = FALSE,
    unique = TRUE,
    .var.name = "metrics"
  )
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
    cli::cli_abort(c(
      "{cli::qty(length(invalid_metrics))}Unknown metric{?s} in {.arg metrics}: {fmt_trunc(invalid_metrics, 6)}.",
      "i" = "Available metrics: {.val {names(.metrics_list)}}."
    ))
  }
  metric_fns <- .metrics_list[metrics]
  # always compute n and n_miss per result element
  results$n <- vapply(results$result, nrow, integer(1))
  results$n_miss <- vapply(
    results$result,
    function(x) sum(is.na(x$estimate)),
    integer(1)
  )
  metrics_list <- lapply(
    results$result,
    \(x) calc_all_metrics(x, metric_fns = metric_fns)
  )
  keep_cols <- setdiff(names(results), "result")
  n_metrics <- vapply(metrics_list, nrow, integer(1))
  row_data <- results[
    rep(seq_len(nrow(results)), n_metrics),
    keep_cols,
    drop = FALSE
  ]
  metric_df <- collapse::rowbind(metrics_list)
  out <- cbind(row_data, metric_df)
  rownames(out) <- NULL
  out
}
