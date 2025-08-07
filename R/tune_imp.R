#' Inject NA Values into a Matrix
#'
#' This function randomly selects positions in a matrix to inject a specified number of NA values,
#' ensuring that the injection does not exceed specified missingness thresholds for rows and columns.
#' It attempts to find a valid set of positions within a maximum number of iterations.
#'
#' @inheritParams SlideKnn
#' @param num_na The number of missing values used to estimate prediction quality. Must be a positive integer.
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
#' mat <- matrix(1:100, nrow = 10, ncol = 10)
#' # Inject 10 NAs
#' na_positions <- inject_na(mat, num_na = 10)
#' mat[na_positions] <- NA
#'
#' @export
inject_na <- function(
    obj,
    num_na = 100,
    rowmax = 0.9,
    colmax = 0.9,
    max_iter = 1000) {
  checkmate::assert_number(colmax, lower = 0, upper = 1)
  checkmate::assert_number(rowmax, lower = 0, upper = 1)
  checkmate::assert_count(max_iter, positive = TRUE)
  checkmate::assert_count(num_na, positive = TRUE)

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

#' Tune Parameters for \code{SlideKnn}/\code{knn_imp} Imputation
#'
#' This function tunes the parameters for the [SlideKnn()]/[knn_imp()] imputation method by injecting missing values
#' into the dataset multiple times and evaluating the imputation performance for different parameter
#' combinations. Can also tune custom imputation functions. See details.
#'
#' @details
#' This function allows tuning of hyperparameters for matrix imputation methods, including the built-in 'SlideKnn' and 'knn_imp',
#' or a custom function provided to `.f`.
#'
#' For a custom function in `.f`, the `parameters` data.frame must have columns whose names match the argument names of `.f`
#' (excluding `obj`). The custom function must take `obj` as its first input argument and output an imputed numeric matrix
#' of the same dimensions as `obj`.
#'
#' For the built-in methods ('SlideKnn' or 'knn_imp'), certain parameters are required in `parameters` (e.g., `n_feat`, `k`, `n_overlap` for 'SlideKnn';
#' `k` for 'knn_imp'), and defaults are set if not provided (e.g., `method = "euclidean"`, `post_imp = FALSE`).
#'
#' `tune_imp` can be parallelize over iterations (number of rows of `parameters` * number of repetitions per row) with
#' the argument `cores`. Note that for `SlideKnn`, `tune_imp` will run sequentially but cores will parallelize within
#' `SlideKnn`.
#'
#' The output is a tibble with columns for each parameter combination, a `param_set` identifier, the repetition number (`rep`),
#' and a nested `results` column containing a tibble of `truth` (original values) and `estimate` (imputed values) for the injected NAs.
#' Metrics can be computed from these results, as shown in the examples.
#'
#' @inheritParams SlideKnn
#' @param obj A numeric matrix with \strong{samples in rows} and \strong{features in columns}. Note: keep `obj` small since
#' this function doesn't support `bigmemory`.
#' @param parameters A data frame specifying the parameter combinations to tune. Duplicated rows are removed.
#'   The required columns depend on `.f`; see [knn_imp()] or [SlideKnn()] for details about the parameters, and `@details` for custom functions.
#' @param .f The imputation function to tune. Can be the string "SlideKnn" (default), "knn_imp", or a custom function.
#' @param rep The number of repetitions for injecting missing values to evaluate each combination of parameters. Default is 1.
#'
#' @inheritParams inject_na
#'
#' @seealso [knn_imp()], [SlideKnn()]
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
#' # Tune `SlideKnn` function on a subset of khanmiss1
#' obj <- t(khanmiss1)[1:30, sample.int(nrow(khanmiss1), size = 200)]
#' anyNA(obj)
#' results <- tune_imp(obj, parameters, rep = 1)
#'
#' # # Install `{yardstick}` or calculate any other metrics using the result
#' # library(yardstick)
#' # met_set <- metric_set(mae, rmse, rsq)
#' # results$metrics <- lapply(
#' # results$results,
#' # function(x) {
#' # met_set(x, truth = truth, estimate = estimate)
#' # }
#' # )
#' # # Unnest the metrics
#' # tidyr::unnest(dplyr::select(results, -results), cols = "metrics")
#'
#' @export
tune_imp <- function(
    obj,
    parameters,
    .f = "SlideKnn",
    rep = 1,
    num_na = 100,
    max_iter = 1000,
    .progress = FALSE,
    rowmax = 0.9,
    colmax = 0.9,
    cores = 1) {
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
    "`.f` must be a function or 'SlideKnn' or 'knn_imp'." = (is.function(.f) ||
      is.character(.f) &&
        (.f %in% c("SlideKnn", "knn_imp")) &&
        length(.f) == 1)
  )
  checkmate::assert_count(rep, positive = TRUE, .var.name = "rep")
  checkmate::assert_count(num_na, positive = TRUE, .var.name = "num_na")
  checkmate::assert_count(max_iter, positive = TRUE, .var.name = "max_iter")
  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_number(
    rowmax,
    lower = 0,
    upper = 1,
    null.ok = FALSE,
    .var.name = "rowmax"
  )
  checkmate::assert_number(
    colmax,
    lower = 0,
    upper = 1,
    null.ok = FALSE,
    .var.name = "colmax"
  )
  checkmate::assert_integerish(
    cores,
    lower = 1,
    len = 1,
    null.ok = FALSE,
    .var.name = "cores"
  )
  checkmate::assert_data_frame(
    parameters,
    any.missing = FALSE,
    all.missing = FALSE,
    min.rows = 1,
    col.names = "unique",
    .var.name = "parameters",
    null.ok = FALSE
  )
  if (is.character(.f)) {
    parameters <- unique(parameters)
    parameters$rowmax <- rowmax
    parameters$colmax <- colmax
    if (!"method" %in% names(parameters)) {
      parameters$method <- "euclidean"
    }
    if (!"weighted" %in% names(parameters)) {
      parameters$weighted <- FALSE
    }
    if (!"dist_pow" %in% names(parameters)) {
      parameters$dist_pow <- 1
    }
    if (.f == "SlideKnn") {
      stopifnot(
        "`SlideKnn` requires `n_feat`, `k`, and `n_overlap` in parameters" = c(
          "n_feat",
          "k",
          "n_overlap"
        ) %in%
          names(parameters)
      )
      parameters$.progress <- FALSE
      # For `SlideKnn`, cores passed to `SlideKnn` instead
      parameters$cores <- cores
      cores <- 1
      parameters <- subset(
        parameters,
        select = c(
          "n_feat",
          "k",
          "n_overlap",
          "rowmax",
          "colmax",
          "post_imp",
          "method",
          ".progress",
          "weighted",
          "dist_pow",
          "cores"
        )
      )
    } else if (.f == "knn_imp") {
      stopifnot(
        "`knn_imp` requires `k` in parameters" = "k" %in% names(parameters)
      )
      # If cores is less than number of work total, then parallel the work and sequential in each task
      if (cores <= nrow(parameters) * rep) {
        parameters$cores <- 1
      } else {
        # Else sequential the work, and parallel in each task
        parameters$cores <- cores
        cores <- 1
      }
      parameters <- subset(
        parameters,
        select = c("k", "rowmax", "colmax", "post_imp", "method", "cores", "weighted", "dist_pow")
      )
      .f <- knn_imp
    }
  }
  parameters <- unique(parameters)
  .rowid <- seq_len(nrow(parameters))
  parameters_list <- lapply(split(parameters, f = as.factor(.rowid)), as.list)
  indices <- tibble::as_tibble(expand.grid(
    param_set = .rowid,
    rep = seq_len(rep)
  ))
  # Generate na_loc
  na_loc <- replicate(
    n = rep,
    inject_na(
      obj = obj,
      num_na = num_na,
      rowmax = rowmax,
      colmax = colmax,
      max_iter = max_iter
    ),
    simplify = FALSE
  )
  # strip dimnames to make object slightly smaller
  rn <- rownames(obj)
  cn <- colnames(obj)
  rownames(obj) <- NULL
  colnames(obj) <- NULL
  # setup parallelization
  if (cores > 1) {
    tryCatch(
      mirai::require_daemons(),
      error = function(e) {
        stop(sprintf(
          "%d cores requested, but no mirai daemon is setup. Call mirai::daemons(%d) to set up the parallelization",
          cores,
          cores
        ))
      }
    )
    fn <- purrr::in_parallel
  } else {
    fn <- carrier::crate
  }
  if (is.character(.f) && .f == "SlideKnn") {
    # Can't crate SlideKnn, has to run in sequential mode. But can parallelize within SlideKnn
    indices$result <- purrr::map(
      seq_len(nrow(indices)),
      function(i) {
        tryCatch(
          {
            pre <- obj
            # Inject NA based on pre-calculated position
            pre[na_loc[[indices[i, "rep", drop = TRUE]]]] <- NA
            truth_vec <- obj[na_loc[[indices[i, "rep", drop = TRUE]]]]
            param_vec <- parameters_list[[indices[i, "param_set", drop = TRUE]]]
            estimate_vec <- SlideKnn(
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
              .progress = FALSE
            )[na_loc[[indices[i, "rep", drop = TRUE]]]]
            tibble::tibble(truth = truth_vec, estimate = estimate_vec)
          },
          error = function(e) {
            message(e)
            tibble::tibble(truth = numeric(), estimate = numeric())
          }
        )
      },
      .progress = .progress
    )
  } else {
    fun <- function() {}
    indices$result <- purrr::map(
      seq_len(nrow(indices)),
      fn(
        function(i) {
          tryCatch(
            {
              pre <- obj
              # Inject NA based on pre-calculated position
              pre[na_loc[[indices[i, "rep", drop = TRUE]]]] <- NA
              truth_vec <- obj[na_loc[[indices[i, "rep", drop = TRUE]]]]
              estimate_vec <- do.call(
                fun,
                args = c(
                  list(obj = pre),
                  parameters_list[[indices[i, "param_set", drop = TRUE]]]
                )
              )[na_loc[[indices[i, "rep", drop = TRUE]]]]
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
        parameters_list = parameters_list
      ),
      .progress = .progress
    )
  }
  if (any(vapply(indices$result, nrow, double(1)) == 0)) {
    warning("Some iteration(s) failed. Check the results carefully.")
  }
  indices <- tibble::as_tibble(cbind(parameters[indices$param_set, ], indices))
  rownames(obj) <- rn
  colnames(obj) <- cn
  return(indices)
}
