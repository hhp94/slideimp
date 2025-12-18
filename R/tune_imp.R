#' Inject NA Values into a Matrix
#'
#' This helper function randomly selects positions in a matrix to inject a specified number of NA values,
#' ensuring that the injection does not exceed specified missingness thresholds for rows and columns.
#' It attempts to find a valid set of positions within a maximum number of iterations.
#'
#' @inheritParams slide_imp
#' @param num_na The number of missing values used to estimate prediction quality.
#' @param rowmax Number between 0 to 1. NA injection cannot create rows with more missing % than this number.
#' @param colmax Number between 0 to 1. NA injection cannot create cols with more missing % than this number.
#' @param check_sd Check if after NA injections zero variance columns are created or not.
#' @param max_iter Maximum number of iterations to attempt finding valid NA positions (default to 1000).
#'
#' @return A vector of integer indices indicating the positions in the matrix
#' where NAs should be injected.
#'
#' @details
#' The function uses the `num_na` parameter to determine the number of NAs to inject.
#' It then repeatedly samples random positions from existing non-NA elements and checks if injecting NAs
#' at those positions would exceed the missingness thresholds for any row or column (accounting for existing NAs).
#' If no valid set is found within `max_iter` attempts, an error is thrown.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(1:100, nrow = 10, ncol = 10)
#' # Inject 10 NAs
#' na_positions <- inject_na(mat, num_na = 10)
#' mat[na_positions] <- NA
#' }
#' @keywords internal
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
    if (any(abs(cvars) < .Machine$double.eps^0.5 | is.nan(cvars) | is.na(cvars))) {
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
      sd_ok <- !any(abs(cvars) < .Machine$double.eps^0.5 | is.nan(cvars) | is.na(cvars))
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
#' @description
#' Tunes hyperparameters for imputation methods such as [slide_imp()], [knn_imp()], [pca_imp()],
#' or user-supplied custom functions by repeated cross-validation.
#'
#' @details
#' The function supports tuning for built-in imputation methods ("slide_imp", "knn_imp", "pca_imp")
#' or custom functions provided via `.f`.
#'
#' When using a custom `.f`, the columns in `parameters` must correspond to the arguments of `.f`
#' (excluding the `obj` argument). The custom function must accept `obj` (a numeric matrix) as its
#' first argument and return a numeric matrix of identical dimensions.
#'
#' Tuning results can be evaluated using metrics from the `{yardstick}` package or via the
#' [compute_metrics()] helper function.
#'
#' @inheritParams inject_na
#' @inheritParams slide_imp
#' @param parameters A data.frame specifying parameter combinations to tune, where each column
#'   represents a parameter accepted by `.f` (excluding `obj`). List columns are supported
#'   for complex parameters. Duplicate rows are automatically removed. When `.f = NULL`, the
#'   imputation method is inferred from the column names:
#'   - `k` → K-NN imputation ([knn_imp()])
#'   - `ncp` → PCA imputation ([pca_imp()])
#'   - `k` or `ncp` with `n_feat` and `n_overlap` → sliding window imputation ([slide_imp()])
#' @param rep Either an integer specifying the number of repetitions for random NA injection, or
#'   a list defining fixed NA positions for each repetition (in which case `num_na` is ignored).
#'   The list elements can be one of the following formats:
#'   - A two-column integer matrix. The first column is the row index, the second column is the column index.
#'   Each row is an missing value.
#'   - A numeric vector specifying linear locations of NAs.
#' @param .f Custom function to tune. Must accept `obj` as the first argument, accept the arguments in `parameters`,
#' and return a matrix with the same dimension as `obj` (default = `NULL`).
#' @param cores Controls the number of cores to parallelize over for K-NN and sliding window K-NN only.
#' To setup parallelization for PCA and sliding window PCA, use [mirai::daemons()].
#'
#' @return A [tibble::tibble()] with columns from `parameters`, plus `param_set` (unique parameter set ID),
#'   `rep` (repetition index), and `result` (a nested tibble containing `truth` and `estimate`
#'   columns for true and imputed values, respectively).
#'
#' @examples
#' data(khanmiss1)
#' obj <- t(khanmiss1)[1:20, sample.int(nrow(khanmiss1), size = 200)]
#'
#' # Tune full K-NN imputation
#' parameters <- data.frame(k = c(5, 10))
#'
#' # With random NA injection
#' results <- tune_imp(obj, parameters, rep = 1, num_na = 20)
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
#'   rep = na_positions
#' )
#'
#' compute_metrics(results_fixed)
#'
#' # Custom imputation function example
#' custom_imp <- function(obj, mean = 0, sd = 1) {
#'   na_pos <- is.na(obj)
#'   obj[na_pos] <- rnorm(sum(na_pos), mean = mean, sd = sd)
#'   obj
#' }
#' parameters_custom <- data.frame(mean = c(0, 0, 1), sd = c(1, 2, 1))
#' results_custom <- tune_imp(
#'   obj,
#'   parameters_custom,
#'   .f = custom_imp,
#'   rep = 1,
#'   num_na = 20
#' )
#'
#' compute_metrics(results_custom)
#'
#' @export
tune_imp <- function(
  obj,
  parameters,
  .f = NULL,
  rep = 1,
  num_na = 100,
  rowmax = 0.9,
  colmax = 0.9,
  check_sd = FALSE,
  max_iter = 1000,
  .progress = TRUE,
  cores = 1
) {
  fun <- NULL
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

  # .f logics
  if (is.null(.f)) {
    has_k <- "k" %in% names(parameters)
    has_ncp <- "ncp" %in% names(parameters)
    has_slide_params <- any(c("n_feat", "n_overlap") %in% names(parameters))
    if (has_slide_params) {
      if (!all(c("n_feat", "n_overlap") %in% names(parameters))) {
        stop("`parameters` must have both `n_feat` and `n_overlap` for `slide_imp` tuning.")
      }
      if (has_k && has_ncp) {
        stop(
          "For sliding window imputation (slide_imp), cannot have both `k` and `ncp` in `parameters`.\n",
          "Provide only one:\n",
          "- `k` for slide_imp with K-NN\n",
          "- `ncp` for slide_imp with PCA"
        )
      } else if (!has_k && !has_ncp) {
        stop(
          "For sliding window imputation (slide_imp), must provide either `k` or `ncp` in `parameters` along with `n_feat` and `n_overlap`.\n",
          "Provide:\n",
          "- `k` for slide_imp with K-NN\n",
          "- `ncp` for slide_imp with PCA"
        )
      } else {
        .f <- "slide_imp"
      }
    } else if (has_k && has_ncp) {
      stop(
        "Both `k` and `ncp` found in `parameters`.\n",
        "Specify either:\n",
        "- `k` for K-NN imputation (knn_imp)\n",
        "- `ncp` for PCA imputation (pca_imp)\n",
        "- `n_feat` and `n_overlap` plus either `k` or `ncp` for sliding window imputation (slide_imp)"
      )
    } else if (has_k) {
      .f <- "knn_imp"
    } else if (has_ncp) {
      .f <- "pca_imp"
    } else {
      stop(
        "Cannot infer imputation method from `parameters` columns.\n",
        "Either specify `.f` directly, or include one of:\n",
        "- `k` for K-NN imputation (knn_imp)\n",
        "- `ncp` for PCA imputation (pca_imp)\n",
        "- `n_feat` and `n_overlap` plus either `k` or `ncp` for sliding window imputation (slide_imp)"
      )
    }
  }
  # validate .f
  stopifnot(
    "`.f` must be a function or 'slide_imp' or 'knn_imp' or 'pca_imp'." = (
      is.function(.f) || (is.character(.f) && (.f %in% c("slide_imp", "knn_imp", "pca_imp")) && length(.f) == 1)
    )
  )

  if (is.function(.f)) {
    message("Tuning custom function")
  } else {
    message(sprintf("Tuning %s", .f))
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
    purrr::walk(seq_along(rep), \(i) {
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
    })
    rep_is_list <- TRUE
    n_reps <- length(rep)
  } else {
    stop("`rep` must be either a positive integer or a list of NA location vectors")
  }
  checkmate::assert_count(max_iter, positive = TRUE, .var.name = "max_iter")
  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_flag(check_sd, .var.name = "check_sd")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")

  .rowid <- seq_len(nrow(parameters))
  parameters_list <- lapply(split(parameters, f = as.factor(.rowid)), function(row) {
    row_list <- as.list(row)
    row_list <- lapply(row_list, function(x) {
      if (is.list(x) && length(x) == 1) {
        x[[1]]
      } else {
        x
      }
    })
    row_list
  })

  indices <- tibble::as_tibble(expand.grid(
    param_set = .rowid,
    rep = seq_len(n_reps)
  ))

  if (.progress) {
    message("Step 1/2: Injecting NA")
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
  # Setup parallelization
  parallelize <- tryCatch(mirai::require_daemons(), error = function(e) FALSE)

  if (is.character(.f)) {
    if (.f == "slide_imp") {
      parameters$.progress <- FALSE
    }
    if ((.f == "slide_imp" && "k" %in% names(parameters)) || .f == "knn_imp") {
      # for `slide_imp` knn mode. Don't parallel through mirai.
      parameters$cores <- cores
      parallelize <- FALSE
    }
    if ((.f == "slide_imp" && "ncp" %in% names(parameters)) || .f == "pca_imp") {
      check_sd <- TRUE
    }
  }
  if (parallelize) {
    fn <- purrr::in_parallel
  } else {
    fn <- function(x, ...) {
      x
    }
  }
  if (.progress) {
    if (parallelize) {
      message("Running in parallel...")
    } else {
      if (cores > 1) {
        warning(
          sprintf(
            "cores = %d but running **sequential**. Call `mirai::daemons(%d)` to set up the parallelization",
            cores,
            cores
          )
        )
      }
    }
  }

  # Create the crated function based on the type of imputation
  # Determine the target function and validation requirements
  if (is.character(.f)) {
    if (.f == "slide_imp") {
      target_function <- slide_imp
      validate_output <- FALSE
    } else if (.f == "knn_imp") {
      target_function <- knn_imp
      validate_output <- FALSE
    } else if (.f == "pca_imp") {
      target_function <- pca_imp
      validate_output <- FALSE
    }
  } else {
    target_function <- .f
    validate_output <- TRUE
  }
  # Create a single unified crated function
  crated_fn <- fn(
    function(i) {
      tryCatch(
        {
          # Create matrix with injected NAs
          pre <- obj
          na_positions <- na_loc[[indices[i, "rep", drop = TRUE]]]
          pre[na_positions] <- NA
          # Get true values
          truth_vec <- obj[na_positions]
          # Get parameters for this iteration
          param_vec <- parameters_list[[indices[i, "param_set", drop = TRUE]]]
          # Run imputation function
          imputed_result <- do.call(
            target_function,
            args = c(list(obj = pre), param_vec)
          )
          # validate result if it's a custom function
          if (validate_output) {
            checkmate::assert_matrix(
              imputed_result,
              mode = "numeric",
              nrows = nr,
              ncols = nc,
              null.ok = FALSE,
              .var.name = "imputed_result"
            )
            checkmate::assert_true(sum(is.infinite(imputed_result)) == 0, .var.name = "imputed_result")
          }
          estimate_vec <- imputed_result[na_positions]
          tibble::tibble(truth = truth_vec, estimate = estimate_vec)
        },
        error = function(e) {
          message(e)
          tibble::tibble(truth = numeric(), estimate = numeric())
        }
      )
    },
    target_function = target_function,
    validate_output = validate_output,
    obj = obj,
    na_loc = na_loc,
    indices = indices,
    nr = nr,
    nc = nc,
    parameters_list = parameters_list
  )
  if (.progress) {
    message("Step 2/2: Tuning")
  }
  # Execute the mapping with the crated function
  result_list <- purrr::map(seq_len(nrow(indices)), crated_fn, .progress = .progress)
  # Combine parameters with results
  result_df <- tibble::as_tibble(cbind(
    parameters[indices$param_set, , drop = FALSE],
    indices,
    tibble::tibble(result = result_list)
  ))

  class(result_df) <- c("TuneImp", class(result_df))
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
  if (abs(ss_tot) < .Machine$double.eps^0.5) {
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

  tibble::tibble(
    .metric = names(metric_fns),
    .estimator = "standard",
    .estimate = estimates
  )
}

#' Compute Prediction Accuracy Metrics
#'
#' Calculates prediction accuracy metrics for each result in a [tune_imp()] output.
#' Use the `{yardstick}` package for other and faster metrics.
#'
#' @param results A tibble from [tune_imp()] containing a `result` column
#' with tibbles that have `truth` and `estimate` columns.
#' @param metrics A character vector of metric names to compute. Available options
#' are: `"mae"`, `"rmse"`, `"rsq"`, `"mape"`, `"bias"`, `"calc_rsq_trad"`. Defaults
#' to `c("mae", "rmse", "rsq")`.
#'
#' @return A tibble with the original parameters and unnested metrics
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
compute_metrics <- function(results, metrics = c("mae", "rmse", "rsq")) {
  UseMethod("compute_metrics")
}

#' @export
compute_metrics.TuneImp <- function(results, metrics = c("mae", "rmse", "rsq")) {
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
  results$metrics <- lapply(
    results$result,
    \(x) calc_all_metrics(x, metric_fns = metric_fns)
  )
  keep_cols <- setdiff(names(results), c("result", "metrics"))

  out <- do.call(rbind, lapply(seq_len(nrow(results)), function(i) {
    row_data <- results[i, keep_cols, drop = FALSE]
    row_data <- row_data[rep(1, nrow(results$metrics[[i]])), , drop = FALSE]
    row_data <- cbind(row_data, results$metrics[[i]])
    row_data
  }))

  return(tibble::as_tibble(out))
}
