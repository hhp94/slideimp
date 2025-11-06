#' Inject NA Values into a Matrix
#'
#' This helper function randomly selects positions in a matrix to inject a specified number of NA values,
#' ensuring that the injection does not exceed specified missingness thresholds for rows and columns.
#' It attempts to find a valid set of positions within a maximum number of iterations.
#'
#' @inheritParams slide_imp
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param num_na The number of missing values used to estimate prediction quality.
#' @param max_iter Maximum number of iterations to attempt finding valid NA positions (default: 1000).
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
#' @noRd
inject_na <- function(
    obj,
    num_na = NULL,
    rowmax = 0.9,
    colmax = 0.9,
    max_iter = 1000) {
  # subset the matrix to the specified features and samples
  na_mat <- !is.na(obj)
  # check if existing NA pattern already exceeds thresholds
  max_col_miss <- floor(nrow(na_mat) * colmax)
  max_row_miss <- floor(ncol(na_mat) * rowmax)
  # calculate current missingness
  current_col_miss <- nrow(na_mat) - colSums(na_mat)
  current_row_miss <- ncol(na_mat) - rowSums(na_mat)
  # Check if any columns already exceed the threshold
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
  # initialize variables for the while loop
  c_miss <- TRUE
  r_miss <- TRUE
  na_loc <- NULL
  iter <- 0
  # inject NAs while ensuring missingness thresholds and iter are not exceeded
  while (c_miss || r_miss) {
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
  }
  return(na_loc)
}

#' Tune Parameters for Imputation Methods
#'
#' @description
#' Tunes hyperparameters for imputation methods like [slide_imp()], [knn_imp()], [pca_imp()],
#' or custom functions by injecting missing values into the dataset and evaluating
#' performance across parameter combinations.
#'
#' @details
#' Supports tuning built-in methods ('slide_imp', 'knn_imp', 'pca_imp') or custom functions via `.f`.
#'
#' For custom `.f`, `parameters` columns must match `.f`'s arguments (excluding `obj`). The
#' function must take `obj` as an argument and return a numeric matrix of the same dimensions.
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param parameters Data frame of parameter combinations to tune where each column is
#' a parameter accepted by `.f`. Duplicates are removed. Supports list columns.
#' @param .f Imputation function: "knn_imp" (default), "slide_imp", "pca_imp", or custom function.
#' @param rep Positive integer for random NA injections per combination (default 1), or
#' list of integer vectors specifying linear NA positions of a matrix (all unique, same length).
#' @param num_na Number of NAs to inject if `rep` is integer (default 100). Ignored if `rep` is list.
#' @param max_iter Max iterations to find valid NA positions (default 1000).
#' @param colmax Max proportion of NAs per column (0-1, default 1).
#' @param rowmax Max proportion of NAs per row (0-1, default 1).
#' @param .progress Show progress bar (default FALSE).
#' @param cores Number of cores for parallelization (default 1).
#'
#' @return Tibble with parameter columns, `param_set` (ID), `rep` (repetition), and `result` (nested tibble of `truth` and `estimate`).
#'
#' @examples
#' data(khanmiss1)
#'
#' parameters <- data.frame(
#'   n_feat = c(100, 100, 100),
#'   k = c(5, 10, 10),
#'   n_overlap = c(10, 10, 10),
#'   knn_method = "euclidean",
#'   post_imp = FALSE
#' )
#'
#' set.seed(1234)
#' obj <- t(khanmiss1)[1:20, sample.int(nrow(khanmiss1), size = 200)]
#'
#' # Random NA injection
#' results <- tune_imp(obj, parameters, .f = "slide_imp", rep = 1, num_na = 20)
#'
#' # Specific NA positions
#' obj_complete <- obj
#' obj_complete[is.na(obj_complete)] <- 0
#' na_positions <- list(
#'   sample(1:length(obj_complete), 20, replace = FALSE),
#'   sample(1:length(obj_complete), 20, replace = FALSE),
#'   sample(1:length(obj_complete), 20, replace = FALSE)
#' )
#' results_fixed <- tune_imp(obj_complete, data.frame(k = 10), .f = "knn_imp", rep = na_positions)
#'
#' # Custom imputation
#' custom_imp <- function(obj, mean = 0, sd = 1) {
#'   na_pos <- is.na(obj)
#'   obj[na_pos] <- rnorm(sum(na_pos), mean = mean, sd = sd)
#'   obj
#' }
#' parameters_custom <- data.frame(mean = c(0, 0, 1), sd = c(1, 2, 1))
#' results_custom <- tune_imp(obj, parameters_custom, .f = custom_imp, rep = 1, num_na = 20)
#'
#' # Analyze the results with yardstick
#' @examplesIf requireNamespace("dplyr", quietly = TRUE) && requireNamespace("yardstick", quietly = TRUE) && requireNamespace("tidyr", quietly = TRUE)
#' met_set <- yardstick::metric_set(yardstick::mae, yardstick::rmse, yardstick::rsq)
#' results$metrics <- lapply(results$result, \(x) met_set(x, truth = truth, estimate = estimate))
#' tidyr::unnest(dplyr::select(results, -result), metrics)
#'
#' results_custom$metrics <- lapply(results_custom$result, \(x) met_set(x, truth = truth, estimate = estimate))
#' tidyr::unnest(dplyr::select(results_custom, -result), metrics)
#'
#' @export
tune_imp <- function(
    obj,
    parameters,
    .f = "knn_imp",
    rep = 1,
    num_na = 100,
    rowmax = 0.9,
    colmax = 0.9,
    max_iter = 1000,
    .progress = FALSE,
    cores = 1) {
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
  checkmate::assert_true(sum(is.infinite(obj)) == 0, .var.name = "obj")
  stopifnot(
    "`.f` must be a function or 'slide_imp' or 'knn_imp' or 'pca_imp'." = (
      is.function(.f) || (is.character(.f) && (.f %in% c("slide_imp", "knn_imp", "pca_imp")) && length(.f) == 1)
    )
  )
  if (is.numeric(rep)) {
    checkmate::assert_count(rep, positive = TRUE, .var.name = "rep")
    checkmate::assert_count(num_na, positive = TRUE, null.ok = FALSE, .var.name = "num_na")
    checkmate::assert_count(max_iter, positive = TRUE, null.ok = FALSE, .var.name = "max_iter")
    rep_is_list <- FALSE
    n_reps <- rep
  } else if (is.list(rep)) {
    checkmate::assert_list(rep, types = "integerish", unique = TRUE, min.len = 1, .var.name = "rep")
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
  checkmate::assert_number(rowmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")
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
  if (is.character(.f)) {
    if (.f == "slide_imp") {
      parameters$.progress <- FALSE
      if (cores > 1 && "k" %in% names(parameters)) {
        # for `slide_imp` knn mode. Don't parallel through mirai.
        parameters$cores <- cores
        cores <- 1
      }
    }
  }
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
        max_iter = max_iter
      ),
      simplify = FALSE
    )
  }
  # Setup parallelization
  if (cores > 1) {
    tryCatch(
      mirai::require_daemons(),
      error = function(e) {
        stop(sprintf(
          "%d cores requested, but no mirai daemon is setup. Call `mirai::daemons(%d)` to set up the parallelization",
          cores, cores
        ))
      }
    )
    fn <- purrr::in_parallel
  } else {
    fn <- function(x, ...) {
      x
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
              nrows = nrow(obj),
              ncols = ncol(obj),
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
    parameters_list = parameters_list
  )
  # Execute the mapping with the crated function
  result_list <- purrr::map(seq_len(nrow(indices)), crated_fn, .progress = .progress)
  # Combine parameters with results
  result_df <- tibble::as_tibble(cbind(
    parameters[indices$param_set, , drop = FALSE],
    indices,
    tibble::tibble(result = result_list)
  ))
  return(result_df)
}
