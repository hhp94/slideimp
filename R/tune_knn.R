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
      stop("NA injection failed. Adjust 'num_na' or increase 'colmax' and 'rowmax'.")
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

#' Tune Parameters for [SlideKnn()]/[knn_imp()] Imputation
#'
#' This function tunes the parameters for the [SlideKnn()]/[knn_imp()] imputation method by injecting missing values
#' into the dataset multiple times and evaluating the imputation performance for different parameter
#' combinations. Set \code{n_feat = ncol(obj)} and \code{n_overlap = 0} to tune [knn_imp()]
#'
#' @inheritParams SlideKnn
#' @param parameters A data frame specifying the parameter combinations to tune. Must include columns
#'   `n_feat`, `k`, `n_overlap` and `method`. Duplicate rows are automatically removed.
#' @param rep The number of repetitions for injecting missing values and evaluating parameters.
#'   Default is 1.
#' @param verbose Print out progress? Default is TRUE.
#' @inheritParams inject_na
#'
#' @return A tibble containing the tuning results. Each row corresponds to a specific repetition and
#'   parameter combination, with columns:
#'   \itemize{
#'     \item `rep`: The repetition number.
#'     \item `param_id`: The ID of the parameter combination.
#'     \item `n_feat`, `k`, `n_overlap`, `method`: The parameter values used.
#'     \item `result`: A nested tibble with columns `truth` (original values) and `estimate`
#'       (imputed values) for the injected NAs.
#'   }
#'
#' @examples
#' \dontrun{
#' data(khanmiss1)
#'
#' parameters <- dplyr::tibble(
#'   n_feat = c(100, 100, 100),
#'   k = c(5, 10, 10),
#'   n_overlap = c(10, 10, 10),
#'   method = "euclidean"
#' )
#'
#' set.seed(1234)
#'
#' results <- tune_knn(
#'   khanmiss1,
#'   parameters,
#'   rep = 5
#' )
#'
#' # Compute multiple metrics using yardstick
#' library(yardstick)
#' met_set <- metric_set(mae, rmse, rsq)
#' results$metrics <- lapply(
#'   results$result,
#'   function(x) {
#'     met_set(x, truth = truth, estimate = estimate)
#'   }
#' )
#'
#' # Unnest the metrics
#' tidyr::unnest(dplyr::select(results, -result), cols = "metrics")
#'
#' # View metrics for the first result
#' results$metrics[[1]]
#' }
#'
#' @export
tune_knn <- function(
    obj,
    parameters,
    rep = 1,
    num_na = 100,
    max_iter = 1000,
    rowmax = 0.9,
    colmax = 0.9,
    verbose = TRUE,
    cores = 1,
    .parallel = FALSE) {
  checkmate::assert_data_frame(
    parameters,
    any.missing = FALSE,
    all.missing = FALSE,
    min.rows = 1
  )
  checkmate::assert_count(rep, positive = TRUE)
  checkmate::assert_logical(verbose, any.missing = FALSE, len = 1)
  stopifnot(all(c("n_feat", "k", "n_overlap", "method") %in% names(parameters)))
  # de-dup the parameter
  parameters <- unique(subset(parameters, select = c("n_feat", "k", "n_overlap", "method")))
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
  results_list <- vector("list", length = rep * nrow(parameters))
  z <- 1
  # Main loop over repetitions
  .progress <- if (verbose) {
    TRUE
  } else {
    FALSE
  }
  for (i in seq_along(na_loc)) {
    if (verbose) {
      message("Iter ", i, " of ", rep)
    }
    pre <- obj
    # Inject NA based on pre-calculated position
    pre[na_loc[[i]]] <- NA
    truth_vec <- obj[na_loc[[i]]]

    # Inner loop over parameter combinations
    for (j in seq_len(nrow(parameters))) {
      current_params <- parameters[j, ]

      # Imputation
      imputed_vec <- SlideKnn(
        obj = pre,
        n_feat = current_params$n_feat,
        k = current_params$k,
        n_overlap = current_params$n_overlap,
        method = current_params$method,
        rowmax = rowmax,
        colmax = colmax,
        .progress = .progress,
        .parallel = .parallel,
        cores = cores
      )[na_loc[[i]]]

      # Result as nested tibble
      run_result <- tibble::tibble(
        rep = i,
        param_id = j,
        n_feat = current_params$n_feat,
        k = current_params$k,
        n_overlap = current_params$n_overlap,
        result = list(tibble::tibble(truth = truth_vec, estimate = imputed_vec))
      )
      results_list[[z]] <- run_result
      z <- z + 1
    }
  }
  return(do.call(rbind, results_list))
}
