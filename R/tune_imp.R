#' Inject NA Values into a Matrix
#'
#' This helper function randomly selects positions in a matrix to inject a specified number of NA values,
#' ensuring that the injection does not exceed specified missingness thresholds for rows and columns.
#' It attempts to find a valid set of positions within a maximum number of iterations.
#'
#' @inheritParams slide_imp
#' @inheritParams tune_imp
#'
#' @param num_na The number of missing values to inject.
#' @param rowmax Number between 0 and 1. NA injection cannot create rows with a higher proportion of
#'   missing values than this threshold.
#' @param colmax Number between 0 and 1. NA injection cannot create columns with a higher proportion of
#'   missing values than this threshold.
#' @param check_sd Logical. If `TRUE`, also ensure that no column becomes zero-variance after injection.
#' @param max_iter Maximum number of iterations to attempt finding valid NA positions (default `1000`).
#'
#' @return A vector of integer indices indicating the positions in the matrix
#' where NAs should be injected.
#'
#' @details
#' The function first checks that the existing missingness in `obj` already respects the `rowmax`
#' and `colmax` thresholds. It then repeatedly samples random positions from the currently non-NA
#' elements and verifies that adding the new NAs would not violate the row/column missingness limits
#' (or create zero-variance columns when `check_sd = TRUE`). Sampling continues until a valid set of
#' positions is found or `max_iter` is reached (in which case an error is thrown).
#'
#' @examples
#' #' set.seed(123)
#' #'
#' #' # Create a sample matrix (no existing NAs)
#' #' mat <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' #'
#' #' # Basic usage - inject 10 NAs (defaults: rowmax/colmax = 0.9, check_sd = FALSE)
#' #' na_loc <- inject_na(mat, num_na = 10)
#'
#' @export
inject_na <- function(
  obj,
  num_na,
  rowmax = 0.9,
  colmax = 0.9,
  check_sd = FALSE,
  max_iter = 1000
) {
  checkmate::assert_matrix(obj, min.rows = 1, min.cols = 1, .var.name = "obj")
  checkmate::assert_int(num_na, lower = 1, .var.name = "num_na")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(check_sd, .var.name = "check_sd")
  checkmate::assert_int(max_iter, lower = 1, .var.name = "max_iter")

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
  # inject NAs while ensuring missingness thresholds and iter are not exceeded
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

#' Convert 2D Positions to Linear Indices
#'
#' This function converts 2D positions (row and column indices, i.e., first and second columns)
#' in a matrix to their corresponding linear (1D) positions.
#'
#' @param pos_2d A numeric matrix with exactly 2 columns (first for rows, second for columns)
#' and at least 1 row. Each entry must be a positive integer within the matrix bounds.
#' @param nrow The number of rows in the matrix.
#' @param ncol The number of columns in the matrix.
#' @return A numeric vector of linear positions corresponding to the input 2D positions.
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

#' Tune Parameters for Imputation Methods
#'
#' Tunes hyperparameters for imputation methods such as [slide_imp()], [knn_imp()], [pca_imp()],
#' or user-supplied custom functions by repeated cross-validation.
#'
#' @details
#' The function supports tuning for built-in imputation methods ("slide_imp", "knn_imp", "pca_imp")
#' or custom functions provided via `.f`.
#'
#' When `.f` is a character string, the columns in `parameters` are validated against the chosen
#' method's requirements:
#'   - `"knn_imp"`: requires `k` in `parameters`
#'   - `"pca_imp"`: requires `ncp` in `parameters`
#'   - `"slide_imp"`: requires `window_size` and `overlap_size`, plus exactly one of `k` or `ncp`
#'
#' When `.f` is a custom function, the columns in `parameters` must correspond to the arguments of `.f`
#' (excluding the `obj` argument). The custom function must accept `obj` (a numeric matrix) as its
#' first argument and return a numeric matrix of identical dimensions.
#'
#' Tuning results can be evaluated using the `{yardstick}` package or [compute_metrics()].
#'
#' @inheritParams slide_imp
#' @param num_na The number of missing values used to estimate prediction quality.
#' @param rowmax Number between 0 to 1. NA injection cannot create rows with more missing % than this number.
#' @param colmax Number between 0 to 1. NA injection cannot create cols with more missing % than this number.
#' @param check_sd Check if after NA injections zero variance columns are created or not.
#' @param max_iter Maximum number of iterations to attempt finding valid NA positions (default to 1000).
#' @param parameters A data.frame specifying parameter combinations to tune, where each column
#' represents a parameter accepted by `.f` (excluding `obj`). List columns are supported
#' for complex parameters. Duplicate rows are automatically removed. `NULL` would be equal to tune the
#' function with default parameters.
#' @param rep Either an integer specifying the number of repetitions for random NA injection, or
#' a list defining fixed NA positions for each repetition (in which case `num_na` is ignored).
#' The list elements can be one of the following formats:
#'   - A two-column integer matrix. The first column is the row index, the second column is the column index.
#'   Each row is an missing value.
#'   - A numeric vector specifying linear locations of NAs.
#' @param .f The imputation method to tune. Either a character string (`"knn_imp"`, `"pca_imp"`,
#' or `"slide_imp"`) specifying a built-in method, or a custom function. Custom functions must
#' accept `obj` as the first argument, accept the arguments in `parameters`, and return a matrix
#' with the same dimensions as `obj`.
#' @param cores Controls the number of cores to parallelize over for K-NN and sliding window K-NN imputation with OpenMP only.
#' To setup parallelization for K-NN without OpenMP, PCA, and sliding window PCA imputation, use `mirai::daemons()`.
#'
#' @return A data.frame with columns from `parameters`, plus `param_set` (unique parameter set ID),
#' `rep` (repetition index), and `result` (a nested data.frame containing `truth` and `estimate`
#' columns for true and imputed values, respectively).
#'
#' @examples
#' data(khanmiss1)
#' obj <- t(khanmiss1)[1:20, sample.int(nrow(khanmiss1), size = 200)]
#'
#' # Tune full K-NN imputation
#' parameters <- data.frame(k = c(5, 10))
#'
#' # With random NA injection
#' results <- tune_imp(obj, parameters, .f = "knn_imp", rep = 1, num_na = 20)
#'
#' # Compute metrics on results
#' compute_metrics(results)
#'
#' # Tune with fixed NA positions (2 repetitions)
#' # Positions must not be NA in the original `obj`
#' na_positions <- list(
#'   matrix(c(1, 2, 3, 1, 1, 1), ncol = 2), # Rows 1-3 in column 1
#'   matrix(c(2, 3, 4, 2, 2, 2), ncol = 2) # Rows 2-4 in column 2
#' )
#' results_fixed <- tune_imp(
#'   obj,
#'   data.frame(k = 10),
#'   .f = "knn_imp",
#'   rep = na_positions
#' )
#'
#' compute_metrics(results_fixed)
#'
#' # Custom imputation function example, with 2 cores parallelization with `mirai::daemons()`
#' # Note that the functions have to be in package::function form. See `?carrier::crate`.
#' custom_imp <- function(obj, mean = 0, sd = 1) {
#'   na_pos <- is.na(obj)
#'   obj[na_pos] <- stats::rnorm(sum(na_pos), mean = mean, sd = sd)
#'   obj
#' }
#'
#' mirai::daemons(2) # Setup 2 cores for parallelization
#' parameters_custom <- data.frame(mean = c(0, 0, 1), sd = c(1, 2, 1))
#' results_custom <- tune_imp(
#'   obj,
#'   parameters_custom,
#'   .f = custom_imp,
#'   rep = 2,
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
  rep = 1,
  num_na = 100,
  rowmax = 0.9,
  colmax = 0.9,
  check_sd = TRUE,
  max_iter = 1000,
  .progress = TRUE,
  cores = 1,
  location = NULL
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
  # if parameters is NULL, create a 1-row placeholder
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

  # validate .f
  if (is.character(.f)) {
    stopifnot(
      "`.f` must be 'slide_imp', 'knn_imp', or 'pca_imp'." =
        .f %in% c("slide_imp", "knn_imp", "pca_imp") && length(.f) == 1
    )

    if (.f == "slide_imp") {
      # validate location argument (required for slide_imp)
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
      # validate required tuning parameters
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

  if (is.numeric(rep)) {
    checkmate::assert_count(rep, positive = TRUE, .var.name = "rep")
    checkmate::assert_count(num_na, positive = TRUE, null.ok = FALSE, .var.name = "num_na")
    checkmate::assert_count(max_iter, positive = TRUE, null.ok = FALSE, .var.name = "max_iter")
    rep_is_list <- FALSE
    n_reps <- rep
  } else if (is.list(rep)) {
    checkmate::assert_list(rep, min.len = 1, types = c("matrix", "integerish"), .var.name = "rep")
    # convert any 2D positions to linear positions
    rep <- lapply(seq_along(rep), \(i) {
      elem <- rep[[i]]
      if (is.matrix(elem) && ncol(elem) == 2) {
        grid_to_linear(elem, nrow = nr, ncol = nc)
      } else {
        elem
      }
    })
    # validate as linear positions
    elem_lengths <- vapply(rep, length, numeric(1))
    if (length(unique(elem_lengths)) != 1) {
      stop("All elements in `rep` list must have the same length")
    }
    for (i in seq_along(rep)) {
      checkmate::assert_integerish(
        rep[[i]],
        lower = 1,
        upper = length(obj),
        any.missing = FALSE,
        min.len = 1,
        unique = TRUE,
        null.ok = FALSE,
        .var.name = sprintf("rep[[%d]]", i)
      )
    }
    rep_is_list <- TRUE
    n_reps <- length(rep)
  } else {
    stop("`rep` must be either a positive integer or a list of NA location vectors")
  }
  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_flag(check_sd, .var.name = "check_sd")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")

  .rowid <- seq_len(nrow(parameters))

  indices <- expand.grid(
    param_set = .rowid,
    rep = seq_len(n_reps)
  )

  message("Step 1/2: Injecting NA")

  # parallelization
  is_knn_mode <- {
    (is.character(.f) && .f == "knn_imp") ||
      (is.character(.f) && .f == "slide_imp" && "k" %in% names(parameters))
  }
  if (is.character(.f)) {
    if ((.f == "slide_imp" && "ncp" %in% names(parameters)) || .f == "pca_imp") {
      check_sd <- TRUE
    }
  }

  # Generate or use NA injection locations
  if (rep_is_list) {
    # Use the provided list of NA locations
    na_loc <- rep
  } else {
    # Generate NA injection locations for each repetition
    na_loc <- replicate(
      n = n_reps,
      inject_na(
        obj = obj,
        num_na = num_na,
        rowmax = rowmax,
        colmax = colmax,
        check_sd = check_sd,
        max_iter = max_iter
      ),
      simplify = FALSE
    )
  }

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

  # Build extra fixed args for slide_imp (location is not a tuning parameter)
  is_slide <- is.character(.f) && .f == "slide_imp"
  fixed_args <- list()
  if (is_slide) fixed_args$location <- location
  if (is_knn_mode) fixed_args$cores <- cores

  # Parameter list
  parameters_list <- lapply(split(parameters, f = as.factor(.rowid)), function(row) {
    row_list <- as.list(row)
    # remove placeholder column
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

  # validate the output of custom function once and convert character .f to function
  if (is.character(.f)) {
    impute_f <- switch(.f,
      slide_imp = slide_imp,
      knn_imp   = knn_imp,
      pca_imp   = pca_imp
    )
  } else {
    # custom function: probe once on the main process to validate output shape
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
    # crated function only needed for parallel branch
    crated_fn <- carrier::crate(
      function(i) {
        # pin BLAS/OMP threads inside each worker. Doesn't modify user's BLAS.
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
        tryCatch(
          {
            pre <- obj
            na_positions <- na_loc[[indices[i, "rep", drop = TRUE]]]
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
      fixed_args = fixed_args
    )
    result_list <- mirai::mirai_map(iter, crated_fn)[.progress = .progress]
  } else {
    # sequential: plain closure, no crate, no BLAS pinning
    run_one <- function(i) {
      tryCatch(
        {
          pre <- obj
          na_positions <- na_loc[[indices[i, "rep", drop = TRUE]]]
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
      out <- run_one(i)
      if (.progress) cli::cli_progress_update(id = pb)
      out
    })
    if (.progress) cli::cli_progress_done(id = pb)
  }

  # combine parameters with results
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
  results$metrics <- lapply(
    results$result,
    \(x) calc_all_metrics(x, metric_fns = metric_fns)
  )

  keep_cols <- setdiff(names(results), c("result", "metrics"))

  out <- do.call(rbind, lapply(seq_len(nrow(results)), function(i) {
    row_data <- results[i, keep_cols, drop = FALSE]
    metric_df <- results$metrics[[i]]
    row_data <- row_data[rep(1, nrow(metric_df)), , drop = FALSE]
    rownames(row_data) <- NULL
    cbind(row_data, metric_df)
  }))

  as.data.frame(out)
}
