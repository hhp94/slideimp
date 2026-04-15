#' Simulate Matrix with Metadata
#'
#' @description
#' Generates a matrix of random normal data then optionally scale to between 0
#' and 1 between columns. It also creates corresponding data frames for feature
#' (column) and sample (row) metadata and can optionally introduce `NA` values
#' into a specified proportion of rows. A correlation between columns `rho`
#' (before scaling) can be added.
#'
#' @param n An integer specifying the number of rows (samples). Default is `100`.
#' @param p An integer specifying the number of columns (features). Default is `100`.
#' @param rho Columns correlation before scaling (compound symmetry). Default is `0.5`.
#' @param n_col_groups An integer for the number of groups to assign to features/columns. Default is `2`.
#' @param n_row_groups An integer for the number of groups to assign to samples/rows. Default is `1`.
#' @param perc_total_na Proportion of all cells to set to NA. Default is `0.1`.
#' @param perc_col_na Proportion of columns across which those NAs are spread. Default is `0.5`.
#' @param beta If TRUE (default) scale values between 0 and 1 column wise.
#'
#' @returns An object of class `slideimp_sim`. This is a list containing:
#' * `input`: A numeric matrix of dimension \eqn{n \times p} containing
#'   the simulated values and injected `NA`s.
#' * `col_group`: A data frame with $p$ rows mapping each `feature`
#'   to a `group`.
#' * `row_group`: A data frame with $n$ rows mapping each `sample`
#'   to a `group`.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' sim_data <- sim_mat(n = 50, p = 10, rho = 0.5)
#' sim_data
sim_mat <- function(
  n = 100,
  p = 100,
  rho = 0.5,
  n_col_groups = 2,
  n_row_groups = 1,
  perc_total_na = 0.1,
  perc_col_na = 0.5,
  beta = TRUE
) {
  checkmate::assert_int(n, lower = 2)
  checkmate::assert_int(p, lower = 2)
  checkmate::assert_number(rho, lower = 0, upper = 1 - .Machine$double.eps^2)
  checkmate::assert_int(n_col_groups, lower = 1, upper = p)
  checkmate::assert_int(n_row_groups, lower = 1, upper = n)
  checkmate::assert_number(perc_total_na, lower = 0, upper = 1)
  checkmate::assert_number(perc_col_na, lower = 0, upper = 1)
  checkmate::assert_flag(beta)

  # correlated columns: var = rho + (1 - rho) = 1, cor = rho
  if (rho > 0) {
    common <- stats::rnorm(n) * sqrt(rho)
    Z <- matrix(stats::rnorm(n * p), n, p) * sqrt(1 - rho)
    d <- common %*% matrix(1, 1, p) + Z
  } else {
    d <- matrix(stats::rnorm(n * p), n, p)
  }

  # per-column [0,1] rescale
  if (beta) {
    mins <- apply(d, 2, min, na.rm = TRUE)
    maxs <- apply(d, 2, max, na.rm = TRUE)
    ranges <- maxs - mins
    ranges[ranges == 0] <- 1
    d <- sweep(d, 2, mins, "-")
    d <- sweep(d, 2, ranges, "/")
  }

  # metadata - columns/features
  col_groups <- sample(
    paste0("group", seq_len(n_col_groups)),
    size = p,
    replace = (n_col_groups < p)
  )
  col_group <- data.frame(
    feature = paste0("feature", seq_len(p)),
    group = col_groups
  )

  # metadata - rows/samples
  colnames(d) <- col_group$feature
  rownames(d) <- paste0("sample", seq_len(n))

  row_groups <- sample(
    paste0("group", seq_len(n_row_groups)),
    size = n,
    replace = (n_row_groups < n)
  )

  row_group <- data.frame(
    sample = rownames(d),
    group = row_groups
  )

  # inject NAs
  n_total_na <- ceiling(perc_total_na * length(d))
  n_cols_na <- max(1, floor(perc_col_na * p))
  if (n_total_na > 0) {
    cols_na <- sample.int(p, n_cols_na)
    per_col <- floor(n_total_na / n_cols_na)
    remainder <- n_total_na - per_col * n_cols_na
    counts <- rep(per_col, n_cols_na) + c(rep(1, remainder), rep(0, n_cols_na - remainder))
    counts <- pmin(counts, n)
    row_idx <- unlist(lapply(counts, function(k) sample.int(n, k)))
    col_idx <- rep(cols_na, times = counts)
    d[cbind(row_idx, col_idx)] <- NA
  }

  structure(
    list(
      input = d,
      col_group = col_group,
      row_group = row_group,
      n_col_groups = n_col_groups,
      n_row_groups = n_row_groups
    ),
    class = "slideimp_sim"
  )
}
