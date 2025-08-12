#' Inject NA Values into a Matrix
#'
#' This helper function randomly selects positions in a matrix to inject a specified number of NA values,
#' ensuring that the injection does not exceed specified missingness thresholds for rows and columns.
#' It attempts to find a valid set of positions within a maximum number of iterations.
#'
#' @inheritParams SlideKnn
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
inject_na <- function(
    obj,
    num_na = 100,
    rowmax = 0.9,
    colmax = 0.9,
    max_iter = 1000) {
  # subset the matrix to the specified features and samples
  na_mat <- !is.na(obj)
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

  # max allowed missing counts per column and row
  max_col_miss <- floor(nrow(na_mat) * colmax)
  max_row_miss <- floor(ncol(na_mat) * rowmax)

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

#' Tune Parameters for [SlideKnn()]/[knn_imp()]/Custom Imputation
#'
#' @description
#' This function tunes the parameters for the [SlideKnn()] or [knn_imp()] imputation methods by injecting
#' missing values into the dataset multiple times and evaluating the imputation performance for different
#' parameter combinations. Can also tune custom imputation functions. Can also accept list of NA locations.
#'
#' @details
#' This function allows tuning of hyperparameters for matrix imputation methods, including the built-in
#' 'SlideKnn' and 'knn_imp', or a custom function provided to `.f`.
#'
#' For a custom function in `.f`, the `parameters` data.frame must have columns whose names match the
#' argument names of `.f` (excluding `obj`). The custom function must take `obj` as its first input
#' argument and output a numeric matrix of the same dimensions as `obj`.
#'
#' For the built-in methods ('SlideKnn' or 'knn_imp'), certain parameters are required in `parameters`:
#' - For 'SlideKnn': `n_feat`, `k`, and `n_overlap`
#' - For 'knn_imp': `k`
#'
#' Default values are set for optional parameters if not provided (e.g., `method = "euclidean"`,
#' `post_imp = FALSE`).
#'
#' **Note:** The `nboot` parameter is always internally set to 1 for tuning purposes.
#'
#' @inheritParams SlideKnn
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#'   Note: keep `obj` small since this function doesn't support `bigmemory`.
#' @param parameters A data frame specifying the parameter combinations to tune. Duplicated rows are
#'   removed. The required columns depend on `.f`; see [knn_imp()] or [SlideKnn()] for details about
#'   the parameters. Any `nboot` values in this data frame will be ignored.
#' @param .f The imputation function to tune. Can be the string "SlideKnn" (default), "knn_imp", or
#'   a custom function. See details.
#' @param rep Either:
#'   - A positive integer specifying the number of repetitions for randomly injecting missing values
#'     to evaluate each parameter combination (default is 1).
#'   - A list of integer vectors, where each vector contains the positions (1-indexed) in the matrix
#'     where NAs should be injected. All vectors must have the same length, and all elements must be
#'     unique (no duplicate NA location sets). The length of the list determines the number of repetitions.
#' @param num_na The number of missing values to inject randomly when `rep` is an integer.
#'   Must be a positive integer when `rep` is an integer. This parameter is ignored (with a warning)
#'   when `rep` is a list. Default is NULL.
#'
#' @inheritParams inject_na
#'
#' @return A tibble containing:
#' - All parameter columns from the input `parameters` data frame
#' - `param_set`: Integer identifier for each unique parameter combination
#' - `rep`: The repetition number (1 to length of `rep` if list, or 1 to `rep` if integer)
#' - `result`: A nested tibble with columns `truth` (original values) and `estimate` (imputed values)
#'
#' @seealso [knn_imp()], [SlideKnn()], [inject_na()]
#'
#' @examples
#' data(khanmiss1)
#'
#' parameters <- data.frame(
#'   n_feat = c(100, 100, 100),
#'   k = c(5, 10, 10),
#'   n_overlap = c(10, 10, 10),
#'   method = "euclidean",
#'   # Set post_imp to FALSE to estimate just the k-NN imputation quality
#'   post_imp = FALSE
#' )
#'
#' set.seed(1234)
#' # Tune SlideKnn function on a subset of khanmiss1
#' obj <- t(khanmiss1)[1:20, sample.int(nrow(khanmiss1), size = 200)]
#' anyNA(obj)
#'
#' # Method 1: Random NA injection with integer rep
#' results <- tune_imp(obj, parameters, .f = "SlideKnn", rep = 1, num_na = 20)
#'
#' # Method 2: Specific NA locations with list rep
#' # Create a complete matrix for demonstration
#' obj_complete <- obj
#' obj_complete[is.na(obj_complete)] <- 0
#'
#' # Define specific positions to test
#' na_positions <- list(
#'   sample(1:length(obj_complete), 20, replace = FALSE),
#'   sample(1:length(obj_complete), 20, replace = FALSE),
#'   sample(1:length(obj_complete), 20, replace = FALSE)
#' )
#'
#' # Tune with predefined NA locations (useful for reproducible benchmarking)
#' results_fixed <- tune_imp(
#'   obj_complete,
#'   parameters,
#'   .f = "SlideKnn",
#'   rep = na_positions  # No num_na needed
#' )
#'
#' # # Install {yardstick} or calculate any other metrics using the result
#' # library(yardstick)
#' # met_set <- metric_set(mae, rmse, rsq)
#' # results$metrics <- lapply(
#' #   results$result,
#' #   function(x) {
#' #     met_set(x, truth = truth, estimate = estimate)
#' #   }
#' # )
#' # # Unnest the metrics
#' # tidyr::unnest(dplyr::select(results, -result), cols = "metrics")
#'
#' # Example with a custom imputation function where missing values are filled with random values
#' custom_imp <- function(obj, mean = 0, sd = 1) {
#'   na_pos <- is.na(obj)
#'   obj[na_pos] <- rnorm(sum(na_pos), mean = mean, sd = sd)
#'   return(obj)
#' }
#'
#' parameters_custom <- data.frame(
#'   mean = c(0, 0, 1),
#'   sd = c(1, 2, 1)
#' )
#'
#' set.seed(1234)
#' # Reuse the same obj
#' results_custom <- tune_imp(obj, parameters_custom, .f = custom_imp, rep = 1, num_na = 20)
#'
#' # # Similarly, compute metrics
#' # results_custom$metrics <- lapply(
#' #   results_custom$result,
#' #   function(x) {
#' #     met_set(x, truth = truth, estimate = estimate)
#' #   }
#' # )
#' # tidyr::unnest(dplyr::select(results_custom, -result), cols = "metrics")
#'
#' @export
tune_imp <- function(
    obj,
    parameters,
    .f = "SlideKnn",
    rep = 1,
    num_na = NULL,
    max_iter = 1000,
    .progress = FALSE,
    rowmax = 0.9,
    colmax = 0.9,
    cores = 1,
    strip_dimnames = FALSE) {
  fun <- NULL
  # Input validation
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
    "`.f` must be a function or 'SlideKnn' or 'knn_imp'." = (
      is.function(.f) || (is.character(.f) && (.f %in% c("SlideKnn", "knn_imp")) && length(.f) == 1)
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
  checkmate::assert_flag(strip_dimnames, .var.name = "strip_dimnames")
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

  # Process parameters based on function type
  if (is.character(.f)) {
    parameters <- unique(parameters)

    # Remove any nboot column if present and force nboot = 1 also disable bigmemory support
    parameters$nboot <- NULL

    # Add fixed parameters
    parameters$rowmax <- rowmax
    parameters$colmax <- colmax

    # Set defaults for optional parameters if not provided
    if (!"method" %in% names(parameters)) {
      parameters$method <- "euclidean"
    }
    if (!"weighted" %in% names(parameters)) {
      parameters$weighted <- FALSE
    }
    if (!"dist_pow" %in% names(parameters)) {
      parameters$dist_pow <- 1
    }
    if (!"post_imp" %in% names(parameters)) {
      parameters$post_imp <- FALSE
    }

    if (.f == "SlideKnn") {
      # Validate required parameters
      required_params <- c("n_feat", "k", "n_overlap")
      missing_params <- setdiff(required_params, names(parameters))
      if (length(missing_params) > 0) {
        stop(sprintf(
          "`SlideKnn` requires %s in parameters",
          paste(missing_params, collapse = ", ")
        ))
      }

      parameters$.progress <- FALSE
      # For SlideKnn, cores passed to SlideKnn instead
      parameters$cores <- cores
      cores <- 1

      # Select only relevant columns for SlideKnn
      valid_cols <- c(
        "n_feat", "k", "n_overlap", "rowmax", "colmax",
        "post_imp", "method", ".progress", "weighted",
        "dist_pow", "cores", "nboot"
      )
      parameters <- parameters[, intersect(names(parameters), valid_cols), drop = FALSE]
    } else if (.f == "knn_imp") {
      # Validate required parameters
      if (!"k" %in% names(parameters)) {
        stop("`knn_imp` requires `k` in parameters")
      }

      # Optimize core allocation
      total_work <- nrow(parameters) * n_reps
      if (cores <= total_work) {
        # Parallelize over iterations
        parameters$cores <- 1
      } else {
        # Parallelize within each knn_imp call
        parameters$cores <- cores
        cores <- 1
      }

      # Select only relevant columns for knn_imp
      valid_cols <- c(
        "k", "rowmax", "colmax", "post_imp", "method",
        "cores", "weighted", "dist_pow", "nboot"
      )
      parameters <- parameters[, intersect(names(parameters), valid_cols), drop = FALSE]
    }
  } else {
    # For custom functions, ensure nboot is not in parameters or set to 1
    if ("nboot" %in% names(parameters)) {
      warning("Removing 'nboot' from parameters for custom function - tuning always uses nboot=1")
      parameters$nboot <- NULL
    }
  }

  # Remove duplicates after processing
  parameters <- unique(parameters)

  # Create parameter sets and repetition indices
  .rowid <- seq_len(nrow(parameters))
  parameters_list <- lapply(split(parameters, f = as.factor(.rowid)), as.list)
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

  # Strip dimnames to reduce object size
  if (strip_dimnames) {
    rn <- rownames(obj)
    cn <- colnames(obj)
    rownames(obj) <- NULL
    colnames(obj) <- NULL
  }
  # Setup parallelization
  if (cores > 1) {
    tryCatch(
      mirai::require_daemons(),
      error = function(e) {
        stop(sprintf(
          "%d cores requested, but no mirai daemon is setup. Call mirai::daemons(%d) to set up the parallelization",
          cores, cores
        ))
      }
    )
    fn <- purrr::in_parallel
  } else {
    fn <- carrier::crate
  }

  # Create the crated function based on the type of imputation
  if (is.character(.f) && .f == "SlideKnn") {
    # SlideKnn must run sequentially (but can parallelize internally)
    # Note: we don't use fn here since SlideKnn handles its own parallelization
    crated_fn <- function(i) {
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

          # Run SlideKnn (returns a list, extract first element)
          imputed_result <- SlideKnn(
            obj = pre,
            n_feat = param_vec$n_feat,
            subset = NULL,
            n_overlap = param_vec$n_overlap,
            k = param_vec$k,
            rowmax = param_vec$rowmax,
            colmax = param_vec$colmax,
            cores = param_vec$cores,
            method = param_vec$method,
            post_imp = param_vec$post_imp,
            weighted = param_vec$weighted,
            dist_pow = param_vec$dist_pow,
            nboot = 1L, # Always use nboot = 1 for tuning
            .progress = FALSE
          )

          # Extract imputed values from first (and only) list element
          estimate_vec <- imputed_result[[1]][na_positions]

          tibble::tibble(truth = truth_vec, estimate = estimate_vec)
        },
        error = function(e) {
          warning(sprintf("Iteration %d failed: %s", i, e$message))
          tibble::tibble(truth = numeric(), estimate = numeric())
        }
      )
    }
  } else if (is.character(.f) && .f == "knn_imp") {
    # knn_imp can be parallelized
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

            # Run knn_imp (returns a list, extract first element)
            imputed_result <- do.call(
              knn_imp,
              args = c(list(obj = pre), param_vec)
            )

            # Extract imputed values from first (and only) list element
            estimate_vec <- imputed_result[[1]][na_positions]

            tibble::tibble(truth = truth_vec, estimate = estimate_vec)
          },
          error = function(e) {
            message(e)
            tibble::tibble(truth = numeric(), estimate = numeric())
          }
        )
      },
      obj = obj,
      na_loc = na_loc,
      indices = indices,
      parameters_list = parameters_list,
      knn_imp = knn_imp # Explicitly pass knn_imp to avoid environment issues
    )
  } else {
    # For custom functions
    nrow_obj <- nrow(obj)
    ncol_obj <- ncol(obj)
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

            # Run imputation function (expects matrix return)
            imputed_result <- do.call(
              fun,
              args = c(list(obj = pre), param_vec)
            )
            checkmate::assert_matrix(
              imputed_result,
              mode = "numeric",
              nrows = nrow_obj,
              ncols = ncol_obj,
              null.ok = FALSE,
              .var.name = "imputed_result"
            )
            checkmate::assert_true(sum(is.infinite(imputed_result)) == 0, .var.name = "imputed_result")

            # Extract imputed values directly from matrix
            estimate_vec <- imputed_result[na_positions]

            tibble::tibble(truth = truth_vec, estimate = estimate_vec)
          },
          error = function(e) {
            message(e)
            tibble::tibble(truth = numeric(), estimate = numeric())
          }
        )
      },
      fun = .f,
      obj = obj,
      na_loc = na_loc,
      indices = indices,
      parameters_list = parameters_list,
      nrow_obj = nrow_obj,
      ncol_obj = ncol_obj
    )
  }

  # Execute the mapping with the crated function
  result_list <- purrr::map(
    seq_len(nrow(indices)),
    crated_fn,
    .progress = .progress
  )

  # Check for failed iterations (empty results or all NA estimates)
  failed_iterations <- vapply(result_list, function(res) {
    nrow(res) == 0 || all(is.na(res$estimate))
  }, logical(1))

  # Combine parameters with results
  result_df <- tibble::as_tibble(cbind(
    parameters[indices$param_set, , drop = FALSE],
    indices,
    tibble::tibble(result = result_list)
  ))

  # Restore dimnames
  if (strip_dimnames) {
    rownames(obj) <- rn
    colnames(obj) <- cn
  }

  return(result_df)
}
