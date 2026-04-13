#' Stratified Sampling of NA Locations
#'
#' Randomly samples matrix positions for NA injection in a **stratified** manner
#' (balanced by column). For each replication, it randomly selects `n_cols`
#' columns (from all columns or the restricted `subset` if supplied) and then
#' samples exactly `n_rows` row positions **within each** of those columns.
#'
#' @inheritParams slide_imp
#' @inheritParams tune_imp
#'
#' @param n_cols Integer. Number of columns to sample per replication.
#' @param n_rows Integer. Number of rows to sample **per selected column**.
#' @param n_reps Integer. Number of independent replications to generate.
#' @param subset Integer or character vector (optional). Column indices (or
#' column names) to sample from. If supplied, `n_cols` columns are randomly
#' sampled from this pool for **each** replication (instead of from all
#' columns). If `NULL` (default), samples from all columns.
#'
#' @return A list of length `n_reps`. Each element is an integer matrix with
#' two columns named `"row"` and `"col"` containing the 1-based row and
#' column indices of the positions where NAs should be injected. Each matrix
#' has exactly `n_rows * n_cols` rows.
#'
#' @seealso [sample_na_loc()]
#'
#' @examples
#' set.seed(123)
#'
#' # Create a sample matrix (complete data)
#' mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#'
#' # Sample 2 columns with 15 NAs each, repeated 5 times
#' na_locs <- sample_na_loc_stratified(mat, n_cols = 2, n_rows = 15, n_reps = 5)
#'
#' # Check structure
#' length(na_locs) # 5 replications
#' nrow(na_locs[[1]]) # 30 positions per replication
#' table(na_locs[[1]][, "col"]) # exactly 15 per column
#'
#' # Example usage: inject NAs into a test matrix
#' mat_test <- mat
#' loc <- na_locs[[1]]
#' mat_test[cbind(loc[, "row"], loc[, "col"])] <- NA
#'
#' # Using a restricted pool of columns (still randomly samples n_cols from
#' # the pool for each replication)
#' na_locs_fixed <- sample_na_loc_stratified(
#'   mat,
#'   n_cols = 2,
#'   n_rows = 10,
#'   n_reps = 3,
#'   subset = c(1, 5, 8)
#' )
#' @export
sample_na_loc_stratified <- function(obj, n_cols, n_rows, n_reps, subset = NULL) {
  checkmate::assert_matrix(obj, min.rows = 1, min.cols = 1, .var.name = "obj")
  checkmate::assert_int(n_cols, lower = 1, .var.name = "n_cols")
  checkmate::assert_int(n_rows, lower = 1, .var.name = "n_rows")
  checkmate::assert_int(n_reps, lower = 1, .var.name = "n_reps")
  checkmate::assert_true(n_rows <= nrow(obj), .var.name = "n_rows <= nrow(obj)")

  # resolve `subset` into an integer pool of column indices to sample from.
  pool <- if (is.null(subset)) {
    seq_len(ncol(obj))
  } else if (is.character(subset)) {
    checkmate::assert_character(subset, any.missing = FALSE, min.len = 1, unique = TRUE, .var.name = "subset")
    if (is.null(colnames(obj))) {
      stop("`subset` is character but `obj` has no colnames.", call. = FALSE)
    }
    if (!all(subset %in% colnames(obj))) {
      missing <- setdiff(subset, colnames(obj))
      stop("`subset` contains colnames not in `obj`: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    match(subset, colnames(obj))
  } else {
    checkmate::assert_integerish(
      subset,
      lower = 1, upper = ncol(obj),
      any.missing = FALSE, min.len = 1, unique = TRUE,
      .var.name = "subset"
    )
    as.integer(subset)
  }

  checkmate::assert_true(n_cols <= length(pool), .var.name = "n_cols <= length(subset)")

  replicate(
    n = n_reps,
    {
      chosen_cols <- sample(pool, size = n_cols)
      values <- lapply(chosen_cols, function(x) {
        chosen_rows <- sample.int(nrow(obj), size = n_rows)
        matrix(
          c(chosen_rows, rep(x, n_rows)),
          ncol = 2,
          dimnames = list(NULL, c("row", "col"))
        )
      })
      do.call(rbind, values)
    },
    simplify = FALSE
  )
}

#' Sample Missing Locations
#'
#' This helper function randomly selects positions in a matrix to inject a
#' specified number of NA values, ensuring that the injection does not exceed
#' specified missingness thresholds for rows and columns. It attempts to find a
#' valid set of positions within a maximum number of iterations.
#'
#' @inheritParams slide_imp
#' @inheritParams tune_imp
#'
#' @param num_na The number of missing values to inject **per replication**.
#' @param n_reps Integer. Number of independent replications (default `1`).
#' @param rowmax Number between 0 and 1. NA injection cannot create rows with a
#'   higher proportion of missing values than this threshold.
#' @param colmax Number between 0 and 1. NA injection cannot create columns with
#'   a higher proportion of missing values than this threshold.
#' @param check_sd Logical. If `TRUE`, also ensure that no column becomes
#'   zero-variance after injection.
#' @param max_iter Maximum number of iterations to attempt finding valid NA
#'   positions (default `1000`).
#'
#' @return A list of length `n_reps`. Each element is an integer vector of
#' linear (1-based) indices indicating the positions in the matrix where NAs
#' should be injected.
#'
#' @details
#' The function first checks that the existing missingness in `obj` already
#' respects the `rowmax` and `colmax` thresholds. It then repeatedly samples
#' random positions from the currently non-NA elements and verifies that adding
#' the new NAs would not violate the row/column missingness limits (or create
#' zero-variance columns when `check_sd = TRUE`). Sampling continues until a
#' valid set of positions is found or `max_iter` is reached (in which case an
#' error is thrown).
#'
#' @examples
#' set.seed(123)
#'
#' # Create a sample matrix (no existing NAs)
#' mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
#'
#' # Basic usage - inject 10 NAs (defaults: rowmax/colmax = 0.9,
#' # check_sd = FALSE, n_reps = 1)
#' na_locs <- sample_na_loc(mat, num_na = 10)
#' na_indices <- na_locs[[1]]   # linear indices for the single replication
#'
#' # Multiple replications
#' na_locs_multi <- sample_na_loc(mat, num_na = 15, n_reps = 3)
#' @export
sample_na_loc <- function(
  obj,
  num_na,
  n_reps = 1,
  rowmax = 0.9,
  colmax = 0.9,
  check_sd = FALSE,
  max_iter = 1000
) {
  # TODO: improve efficiency. The current rejection-sampling loop rebuilds
  # na_mat_test and recomputes col/row sums on every iteration, which is
  # wasteful for large matrices or tight thresholds. We can use rejection
  # resampling for a much better function.
  checkmate::assert_matrix(obj, min.rows = 1, min.cols = 1, .var.name = "obj")
  checkmate::assert_int(num_na, lower = 1, .var.name = "num_na")
  checkmate::assert_int(n_reps, lower = 1, .var.name = "n_reps")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(check_sd, .var.name = "check_sd")
  checkmate::assert_int(max_iter, lower = 1, .var.name = "max_iter")

  na_mat <- !is.na(obj)
  max_col_miss <- floor(nrow(na_mat) * colmax)
  max_row_miss <- floor(ncol(na_mat) * rowmax)
  current_col_miss <- nrow(na_mat) - colSums(na_mat)
  current_row_miss <- ncol(na_mat) - rowSums(na_mat)
  bad_cols <- which(current_col_miss > max_col_miss)
  if (length(bad_cols) > 0) {
    stop("Some columns have missing > colmax before na injections")
  }
  bad_rows <- which(current_row_miss > max_row_miss)
  if (length(bad_rows) > 0) {
    stop("Some rows have missing > rowmax before na injections")
  }
  not_na <- which(na_mat)
  if (num_na >= length(not_na)) {
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
    cvars <- col_vars(obj)
    if (any(cvars < .Machine$double.eps | is.na(cvars))) {
      stop("Zero variance columns detected before na injections")
    }
  }

  replicate(
    n = n_reps,
    {
      c_miss <- TRUE
      r_miss <- TRUE
      sd_ok <- !check_sd

      na_loc <- NULL
      iter <- 0
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

        col_miss_count <- nrow(na_mat_test) - colSums(na_mat_test)
        row_miss_count <- ncol(na_mat_test) - rowSums(na_mat_test)
        c_miss <- any(col_miss_count > max_col_miss)
        r_miss <- any(row_miss_count > max_row_miss)

        if (check_sd) {
          obj_test <- obj
          obj_test[na_loc] <- NA
          cvars <- col_vars(obj_test)
          sd_ok <- !any(cvars < .Machine$double.eps | is.na(cvars))
        }
      }
      na_loc
    },
    simplify = FALSE
  )
}

#' Convert 2D Positions to Linear Indices
#'
#' This function converts 2D positions (row and column indices, i.e., first and
#' second columns) in a matrix to their corresponding linear (1D) positions.
#'
#' @param pos_2d A numeric matrix with exactly 2 columns (first for rows,
#'   second for columns) and at least 1 row. Each entry must be a positive
#'   integer within the matrix bounds.
#' @param nrow The number of rows in the matrix.
#' @param ncol The number of columns in the matrix.
#' @return A numeric vector of linear positions corresponding to the input 2D
#'   positions.
#'
#' @keywords internal
#' @noRd
grid_to_linear <- function(pos_2d, nrow, ncol) {
  checkmate::assert_matrix(pos_2d, ncols = 2, min.rows = 1, mode = "integerish", .var.name = "pos_2d")
  checkmate::assert_int(nrow, lower = 1, "nrow")
  checkmate::assert_int(ncol, lower = 1, "ncol")

  row <- pos_2d[, 1]
  col <- pos_2d[, 2]

  checkmate::assert_true(all(row >= 1 & row <= nrow), .var.name = "row indices")
  checkmate::assert_true(all(col >= 1 & col <= ncol), .var.name = "column indices")

  linear_pos <- (col - 1) * nrow + row
  return(linear_pos)
}

#' Resolve NA locations for tune_imp
#'
#' If `na_loc` is NULL, generates `n_reps` random NA location vectors via
#' `sample_na_loc()`. Otherwise normalizes the user-supplied positions into
#' a list of linear-index integer vectors and bounds-checks them.
#'
#' Accepted `na_loc` shapes:
#' - integer vector of linear positions (treated as a single rep)
#' - 2-column integerish matrix of (row, col) pairs (treated as a single rep)
#' - list whose elements are either of the above
#'
#' @keywords internal
#' @noRd
resolve_na_loc <- function(
  obj,
  na_loc,
  n_reps,
  num_na,
  rowmax,
  colmax,
  check_sd,
  max_iter
) {
  if (is.null(na_loc)) {
    return(sample_na_loc(
      obj = obj,
      num_na = num_na,
      n_reps = n_reps,
      rowmax = rowmax,
      colmax = colmax,
      check_sd = check_sd,
      max_iter = max_iter
    ))
  }

  # wrap a single vector/matrix into a length-1 list
  if (!is.list(na_loc)) {
    if (is.matrix(na_loc) || is.numeric(na_loc) || is.integer(na_loc)) {
      na_loc <- list(na_loc)
    } else {
      stop("`na_loc` must be a vector, a 2-col matrix, or a list of these.")
    }
  }

  checkmate::assert_list(
    na_loc,
    min.len = 1,
    types = c("matrix", "integerish"),
    .var.name = "na_loc"
  )

  nr <- nrow(obj)
  nc <- ncol(obj)
  n_total <- length(obj)

  na_loc <- lapply(seq_along(na_loc), function(i) {
    elem <- na_loc[[i]]
    if (is.matrix(elem) && ncol(elem) == 2) {
      grid_to_linear(elem, nrow = nr, ncol = nc)
    } else {
      elem
    }
  })

  elem_lengths <- vapply(na_loc, length, numeric(1))
  if (length(unique(elem_lengths)) != 1) {
    stop("All elements in `na_loc` must have the same length")
  }
  for (i in seq_along(na_loc)) {
    checkmate::assert_integerish(
      na_loc[[i]],
      lower = 1,
      upper = n_total,
      any.missing = FALSE,
      min.len = 1,
      unique = TRUE,
      null.ok = FALSE,
      .var.name = sprintf("na_loc[[%d]]", i)
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
#' @param na_loc Optional. Pre-defined missing value locations to bypass random
#'   NA injection. Accepted formats include:
#'   - A two-column integer matrix (row, column indices).
#'   - A numeric vector of linear locations.
#'   - A list where each element is one of the above formats (one per repetition).
#' @param num_na The number of missing values used to estimate prediction
#'   quality. Ignored if `na_loc` is provided (default `500`).
#' @param rowmax,colmax Numbers between 0 and 1. NA injection cannot create
#'   rows/columns with a higher proportion of missing values than these
#'   thresholds.
#' @param check_sd Logical. Check that no column becomes zero-variance after
#'   NA injection (default `TRUE`).
#' @param max_iter Maximum number of iterations to attempt finding valid NA
#'   positions (default `1000`).
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
#'   matrix(c(2, 3, 4, 2, 2, 2), ncol = 2)  # Rows 2-4 in column 2
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
  num_na = 500,
  rowmax = 0.9,
  colmax = 0.9,
  check_sd = TRUE,
  max_iter = 1000,
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

  checkmate::assert_count(n_reps, positive = TRUE, .var.name = "n_reps")
  if (is.null(na_loc)) {
    checkmate::assert_count(num_na, positive = TRUE, null.ok = FALSE, .var.name = "num_na")
    checkmate::assert_count(max_iter, positive = TRUE, null.ok = FALSE, .var.name = "max_iter")
  }
  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_flag(check_sd, .var.name = "check_sd")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")

  # parallelization mode flags (need check_sd finalized before resolving na_loc)
  is_knn_mode <- {
    (is.character(.f) && .f == "knn_imp") ||
      (is.character(.f) && .f == "slide_imp" && "k" %in% names(parameters))
  }
  if (is.character(.f)) {
    if ((.f == "slide_imp" && "ncp" %in% names(parameters)) || .f == "pca_imp") {
      check_sd <- TRUE
    }
  }

  message("Step 1/2: Resolving NA locations")

  na_loc <- resolve_na_loc(
    obj = obj,
    na_loc = na_loc,
    n_reps = n_reps,
    num_na = num_na,
    rowmax = rowmax,
    colmax = colmax,
    check_sd = check_sd,
    max_iter = max_iter
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
#' @param results A slideimp_tune object from [tune_imp()] containing a `result` column
#' with data.frames that have `truth` and `estimate` columns.
#' @param metrics A character vector of metric names to compute. Defaults
#' to `c("mae", "rmse")`. Also available: `"mape"`, `"bias"`, `"rsq"`, and `"rsq_trad"`.
#'
#' @return A data.frame with the original parameters and un-nested metrics
#' (`.metric`, `.estimator`, `.estimate`).
#'
#' @examples
#' data(khanmiss1)
#' set.seed(1234)
#' results <- tune_imp(
#'   obj = t(khanmiss1),
#'   parameters = data.frame(k = 10),
#'   .f = "knn_imp",
#'   rep = 1,
#'   num_na = 20
#' )
#'
#' compute_metrics(results)
#'
#' @export
compute_metrics <- function(results, metrics = c("mae", "rmse")) {
  UseMethod("compute_metrics")
}

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
