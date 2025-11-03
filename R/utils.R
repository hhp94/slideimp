#' Simulate Methylation Beta Values with Metadata
#'
#' @description
#' This function generates a matrix of random normal data, scaled between 0 and 1
#' per column. It also creates corresponding data frames for feature and sample
#' metadata and can optionally introduce `NA` values into a specified proportion
#' of rows.
#'
#' @param n An integer specifying the number of rows (features). Default is `100`.
#' @param m An integer specifying the number of columns (samples). Default is `100`.
#' @param nchr An integer for the number of chromosome groups to assign to features (e.g., `nchr = 22` for human autosomes). Default is `2`.
#' @param ngrp An integer for the number of groups to assign to samples. Default is `1`.
#' @param perc_NA A numeric value between 0 and 1 indicating the proportion of values to set to `NA` within each selected row. Default is `0.5`.
#' @param perc_col_NA A numeric value between 0 and 1 indicating the proportion of rows to select for `NA` introduction. Default is `0.5`.
#' @param beta If TRUE (default) then simulate beta values by scaling the values between 0 and 1.
#'
#' @return A `list` containing three elements:
#' \itemize{
#'   \item `input`: The simulated `n` x `m` numeric matrix with values between 0 and 1.
#'   \item `group_feature`: A `data.frame` with feature IDs and their assigned chromosome group.
#'   \item `group_sample`: A `data.frame` with sample IDs and their assigned group.
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' sim_data <- sim_mat(n = 50, m = 10)
#'
#' # Metadata of each features
#' sim_data$group_feature[1:5, ]
#' sim_data$group_sample[1:5, ]
#'
#' # View the first few rows and columns of the matrix
#' sim_data$input[1:5, 1:5]
#'
#' # Generate a dataset with no missing values
#' sim_data_complete <- sim_mat(n = 50, m = 10, perc_NA = 0, perc_col_NA = 0)
#' sum(is.na(sim_data_complete$input))
sim_mat <- function(
  n = 100,
  m = 100,
  nchr = 2,
  ngrp = 1,
  perc_NA = 0.5,
  perc_col_NA = 0.5,
  beta = TRUE
) {
  checkmate::assert_int(n, lower = 2)
  checkmate::assert_int(m, lower = 2)
  checkmate::assert_int(nchr, lower = 1, upper = 25)
  checkmate::assert_number(perc_NA, lower = 0, upper = 1)
  checkmate::assert_number(perc_col_NA, lower = 0, upper = 1)
  checkmate::assert_flag(beta)
  # Create and scale the matrix to between 0 and 1 per column
  d_length <- n * m
  d <- matrix(stats::rnorm(d_length), nrow = n, ncol = m)
  if (beta) {
    mins <- colMMs(d, 0)[1, ]
    maxs <- colMMs(d, 1)[1, ]
    ranges <- maxs - mins
    d <- sweep(d, 2, mins, "-")
    d <- sweep(d, 2, ranges, "/")
  }
  # Generate realistic feature and sample names
  feature <- seq_len(n)
  chr <- sample(paste0("chr", seq_len(nchr)), size = n, replace = TRUE)
  group_feature <- data.frame(feature_id = paste0("feat", feature), group = chr)
  colnames(d) <- paste0("s", seq_len(m))
  grp <- sample(paste0("grp", seq_len(ngrp)), size = m, replace = TRUE)
  group_sample <- data.frame(sample_id = colnames(d), group = grp)
  rownames(d) <- group_feature$feature_id

  # Introduce missing values in selected features (rows)
  feature_miss_size <- floor(perc_col_NA * n)
  na_size <- floor(perc_NA * m)

  if (feature_miss_size > 0 && na_size > 0) {
    feature_miss <- sample.int(n, size = feature_miss_size)
    col_idx <- c(replicate(feature_miss_size, sample.int(m, na_size)))
    row_idx <- rep(feature_miss, each = na_size)
    d[cbind(row_idx, col_idx)] <- NA
  }

  return(list(
    input = d,
    group_feature = group_feature,
    group_sample = group_sample
  ))
}
