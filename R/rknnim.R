#' Rolling KNN Imputation
#'
#' Performs rolling window KNN imputation on a numeric matrix to handle missing values.
#' The matrix is divided into overlapping windows, and imputation is applied to each window.
#' Overlapping regions are averaged to produce the final imputed matrix.
#'
#' @details
#' See [knn_imp()] for details about the implementation.
#'
#' @param obj A numeric matrix with \strong{samples in rows} and \strong{features in columns}. See details.
#' @param n_feat Integer specifying the number of features (columns) in each window. Must be between 2 and the number of columns in \code{obj}.
#' @param n_overlap Integer specifying the overlap between consecutive windows. Default is 10. Must be between 0 and \code{n_feat - 1}.
#' @param k An integer specifying the number of nearest neighbors to use for
#' imputation. Must be between 1 and the number of columns.
#' @param colmax A numeric value between 0 and 1. This is the threshold for the
#' proportion of missing values in a column. Columns exceeding this
#' threshold will be imputed using the column mean instead of k-NN.
#' @param rowmax A numeric value between 0 and 1. This is the maximum
#' allowable proportion of missing values in any single row. If a row
#' exceeds this threshold, the function will stop with an error.
#' @param cores Number of cores to parallelize calculations of distances over. Note: if \code{.parallel} is TRUE and cores = n, then each [mirai::daemons()] process will spawn n cores.
#' @param method A character string specifying the distance metric for k-NN.
#' Acceptable values are `"euclidean"`, `"manhattan"`, or `"impute.knn"`.
#' Defaults to `"euclidean"`. See details.
#' @param post_imp KNN impute can fail. Retry with mean imputation or not? Default is TRUE.
#' @param .progress Logical indicating whether to show a progress bar. Default is FALSE.
#' @param .parallel Logical indicating whether to use parallel processing for multiple windows. Default is FALSE.
#'
#' @return A numeric matrix of the same dimensions as \code{obj} with missing values imputed.
#'
#' @examples
#' \dontrun{
#' data(khanmiss1)
#' imputed <- rknnim(t(khanmiss1), k = 10, n_feat = 100, n_overlap = 10)
#' }
#'
#' @export
rknnim <- function(
    obj,
    n_feat,
    n_overlap = 10,
    k = 10,
    rowmax = 0.9,
    colmax = 0.9,
    cores = 1,
    method = c("euclidean", "manhattan", "impute.knn"),
    post_imp = TRUE,
    .progress = FALSE,
    .parallel = FALSE) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_matrix(obj, mode = "numeric", null.ok = FALSE)
  checkmate::assert_true(sum(is.infinite(obj)) == 0)
  checkmate::assert_integerish(n_feat, lower = 2, upper = ncol(obj), len = 1, null.ok = FALSE)
  checkmate::assert_integerish(n_overlap, lower = 0, upper = n_feat - 1, len = 1, null.ok = FALSE)
  # Warning for parallel core usage with mirai.
  if (cores > 1 && .parallel) {
    message("Each parallel core used by mirai daemon will spawn additional cores. If this is not intended, set cores = 1 or .parallel = FALSE and cores = N.")
  }

  # Windowing Logic
  # Calculate the total number of steps/windows needed.
  idx <- 1 # R index at 1
  max_step <- ceiling((ncol(obj) - idx) / (n_feat - n_overlap))
  step <- 0:max_step
  start <- idx + (step * n_feat) - (step * n_overlap)
  end <- start + n_feat - 1
  # Edge case, the end might overshoot the number of rows of the obj
  n_overshoot <- sum(end > ncol(obj))

  # In which case trim off the runs that overshoot
  corrected_length <- length(end) - n_overshoot
  start <- start[idx:corrected_length]
  end <- end[idx:corrected_length]

  # And make the last window extra wide to cover the full end
  end[corrected_length] <- ncol(obj)
  # width <- end - start + 1

  # Rolling Imputation
  # Set up the function for mapping (either parallel or sequential).
  if (.parallel) {
    mirai::require_daemons()
    fn <- purrr::in_parallel
  } else {
    fn <- function(x, ...) {
      x
    }
  }

  # Perform imputation for each window using purrr::map.
  # The anonymous function (function(i) { ... }) is applied to each window index.
  # Necessary variables are explicitly passed to the parallel environment.
  imputed <- purrr::map(
    seq_along(start), # Iterate over the indices of the windows
    fn(
      function(i) {
        impute_knn(
          obj = obj[, start[i]:end[i], drop = FALSE],
          k = k,
          rowmax = rowmax,
          colmax = colmax,
          cores = cores,
          method = method,
          post_imp = post_imp
        )
      },
      # Pass all required objects to the parallel function's environment.
      impute_knn = impute_knn,
      obj = obj,
      start = start,
      end = end,
      k = k,
      rowmax = rowmax,
      colmax = colmax,
      cores = cores,
      method = method,
      post_imp = post_imp
    ),
    .progress = .progress
  )

  # Combine Imputed Windows
  # Initialize a zero matrix to accumulate the imputed values.
  final_imputed <- matrix(0, nrow = nrow(obj), ncol = ncol(obj))
  # Initialize a matrix to count how many times each cell has been imputed (for averaging).
  counts <- matrix(0, nrow = nrow(obj), ncol = ncol(obj))
  # Loop through each imputed window to combine the results.
  for (i in seq_along(imputed)) {
    # Define the column range for the current window.
    cols_to_update <- start[i]:end[i]

    # Add the imputed values from the current window to the accumulator matrix.
    final_imputed[, cols_to_update] <- final_imputed[, cols_to_update] + imputed[[i]]

    # Increment the count for these columns, indicating they've been imputed.
    counts[, cols_to_update] <- counts[, cols_to_update] + 1
  }

  # Average the values by dividing the accumulated sums by their respective counts.
  # This handles overlapping regions by averaging their contributions.
  final_imputed <- final_imputed / counts

  if (post_imp) {
    if (anyNA(final_imputed)) {
      final_imputed <- mean_impute_col(final_imputed)
    }
  }
  # Restore original column and row names to the final imputed matrix.
  colnames(final_imputed) <- colnames(obj)
  rownames(final_imputed) <- rownames(obj)

  return(final_imputed)
}

#' KNN Imputation Wrapper
#'
#' A wrapper function for KNN imputation that filters rows based on missing value proportions before imputing.
#' Rows with missing proportions exceeding \code{rowmax} are not imputed in this step to avoid
#' throwing exceptions in [knn_imp()].
#'
#' @inheritParams rknnim
#'
#' @return The imputed matrix, with only qualifying rows imputed via KNN; others remain unchanged.
#'
#' @keywords internal
#' @noRd
impute_knn <- function(obj, k, rowmax, colmax, cores, method, post_imp = TRUE) {
  na_mat <- is.na(obj)

  # Determine 'good_rows': rows where the proportion of NAs is less than 'rowmax'.
  # Rows with too many NAs (as defined by rowmax) are excluded from imputation.
  good_rows <- rowSums(na_mat) / ncol(na_mat) < rowmax

  # If no rows meet the criteria for imputation, return the original object as is.
  if (sum(good_rows) == 0) {
    return(obj)
  }

  imputed_good <- knn_imp(
    obj = obj[good_rows, , drop = FALSE],
    k = k,
    rowmax = rowmax,
    colmax = colmax,
    method = method,
    cores = cores,
    post_imp = post_imp
  )

  # Initialize the result matrix with the original object's values.
  result <- obj
  # Place the imputed values from 'imputed_good' back into the corresponding 'good_rows'
  # in the result matrix.
  result[good_rows, ] <- imputed_good
  return(result)
}

#' K-Nearest Neighbor (k-NN) Imputation
#'
#' This function imputes missing values in a numeric matrix using the k-Nearest
#' Neighbors algorithm. It follows a two-stage process. First, it imputes
#' columns with a proportion of missing values below `colmax` using k-NN.
#' Second, if requested, any remaining missing values are imputed using the column mean.
#'
#' @details
#' This implementation calculates the distances for neighbors column-wise. This is an
#' \strong{extremely} important detail. Outside of microarray data, most datasets have
#' people in columns and features (e.g., weight, height, etc.) in rows for imputation
#' However, in microarray data, genes or CpG sites for the same sample that are
#' spatially closer together carry mutual information, so you can place genes in columns
#' and samples in rows; the algorithm will then impute values based on nearby genes
#' for the same sample.
#'
#' The distance calculation between columns for identifying nearest neighbors is
#' scaled based on the number of non-missing value pairs. Specifically, the
#' raw distance is penalized by scaling it up for columns that have fewer
#' overlapping observations. This penalizes distances for columns with very few
#' shared observations used for distance calculations. The
#' \code{impute.knn} method averages the distances over the number of matching positions,
#' so a column with only one matching value to calculate distance from might have a lower
#' raw distance than a column with many matched values. See also [stats::dist()].
#'
#' @inheritParams rknnim
#'
#' @inherit rknnim return
#'
#' @export
#'
#' @examples
#' # See ?khanmiss1
#' data(khanmiss1)
#' sum(is.na(khanmiss1))
#'
#' # Perform k-NN imputation. `khanmiss1` stores genes in the row so we have to t().
#' # set method to "impute.knn" to mimic how distant is scaled in impute::impute.knn.
#' imputed <- knn_imp(obj = t(khanmiss1), k = 3, colmax = 0.5, rowmax = 0.8, method = "euclidean")
#' imputed[1:5, 1:20]
#' sum(is.na(imputed))
knn_imp <- function(
    obj,
    k,
    colmax = 0.9,
    rowmax = 0.9,
    method = c("euclidean", "manhattan", "impute.knn"),
    cores = 1,
    post_imp = TRUE) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, col.names = "unique")
  checkmate::assert_integerish(k, lower = 1, upper = ncol(obj), len = 1)
  checkmate::assert_integerish(cores, lower = 1, len = 1)
  checkmate::assert_numeric(colmax, lower = 0, upper = 1, len = 1)
  checkmate::assert_numeric(rowmax, lower = 0, upper = 1, len = 1)
  checkmate::assert_logical(post_imp, len = 1, null.ok = FALSE, any.missing = FALSE)

  if (cores >= 16) {
    warning("Setting cores too high can slow down runtime. Benchmark your data first.")
  }

  miss <- is.na(obj)
  if (sum(miss) == 0) {
    message("No missing data")
    return(obj)
  }

  rmiss <- rowSums(miss) / ncol(obj)
  if (any(rmiss >= rowmax)) {
    stop("Row(s) missing exceeded rowmax. Remove rows(s) with too high NA %")
  }

  cmiss <- colSums(miss)
  if (any(cmiss / nrow(obj) == 1)) {
    stop("Col(s) with all missing detected. Remove before proceed")
  }

  # Partition
  knn_imp_cols <- (cmiss / nrow(obj)) < colmax
  pre_imp_cols <- obj[, knn_imp_cols, drop = FALSE]
  pre_imp_miss <- miss[, knn_imp_cols, drop = FALSE]
  pre_imp_cmiss <- cmiss[knn_imp_cols]
  p_pre <- ncol(pre_imp_cols)

  if (any(rowSums(pre_imp_miss) / ncol(pre_imp_cols) == 1)) {
    stop("Row(s) missing exceeded rowmax. Remove rows(s) with too high NA %")
  }

  # Impute
  ## mean imp cols
  obj[, !knn_imp_cols] <- mean_impute_col(obj[, !knn_imp_cols])
  ## knn imp cols
  post_imp_cols <- impute_knn_naive(
    obj = pre_imp_cols,
    miss = pre_imp_miss,
    k = k,
    n_col_miss = pre_imp_cmiss,
    method = switch(method,
      "euclidean" = 0L,
      "manhattan" = 1L,
      "impute.knn" = 2L
    ),
    cores = cores
  )
  colnames(post_imp_cols) <- colnames(pre_imp_cols)
  if (post_imp) {
    if (anyNA(post_imp_cols)) {
      post_imp_cols <- mean_impute_col(post_imp_cols)
    }
  }
  ## Reconstruct obj
  obj[, knn_imp_cols] <- post_imp_cols

  return(obj)
}

#' @title Row Mean Imputation
#'
#' @description
#' Imputes missing values (NA) in a matrix by replacing them with the mean of their
#' respective rows.
#'
#' @inheritParams rknnim
#'
#' @return The matrix with NA values replaced by row means.
#'
#' @export
mean_impute_row <- function(obj) {
  na_indices <- which(is.na(obj), arr.ind = TRUE)
  row_means <- rowMeans(obj, na.rm = TRUE)
  obj[na_indices] <- row_means[na_indices[, 1]]
  return(obj)
}

#' @title Col Mean Imputation
#'
#' @description
#' Imputes missing values (NA) in a matrix by replacing them with the mean of their
#' respective cols.
#'
#' @inheritParams rknnim
#'
#' @return The matrix with NA values replaced by col means.
#'
#' @export
mean_impute_col <- function(obj) {
  na_indices <- which(is.na(obj), arr.ind = TRUE)
  column_means <- colMeans(obj, na.rm = TRUE)
  obj[na_indices] <- column_means[na_indices[, 2]]
  return(obj)
}

dummy <- function() {
  Rcpp::evalCpp()
  carrier::crate()
}
