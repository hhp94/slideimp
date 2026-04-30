#' Simulate a Matrix with Metadata
#'
#' Generate a numeric matrix with optional row and column metadata and added
#' missing values.
#'
#' @description
#' `sim_mat()` generates random normal data with optional compound-symmetric
#' column correlation. Values can optionally be scaled to the interval
#' `[0, 1]` column-wise. The function also creates feature metadata for columns
#' and sample metadata for rows, and can inject `NA` values into a specified
#' proportion of matrix cells across a specified proportion of columns.
#'
#' @param n Integer. Number of rows, interpreted as samples. Defaults to `100`.
#' @param p Integer. Number of columns, interpreted as features. Defaults to
#'   `100`.
#' @param rho Numeric. Compound-symmetric column correlation before optional
#'   scaling. Defaults to `0.5`.
#' @param n_col_groups Integer. Number of groups to assign to features.
#'   Defaults to `2`.
#' @param n_row_groups Integer. Number of groups to assign to samples.
#'   Defaults to `1`.
#' @param perc_total_na Numeric scalar between `0` and `1`. Proportion of all
#'   matrix cells to set to `NA`. Defaults to `0.1`.
#' @param perc_col_na Numeric scalar between `0` and `1`. Proportion of columns
#'   across which injected `NA` values are spread. Defaults to `0.5`.
#' @param beta Logical. If `TRUE`, scale values to the interval `[0, 1]`
#'   column-wise.
#'
#' @returns An object of class `slideimp_sim`. This is a list containing:
#' * `input`: a numeric matrix of dimension \eqn{n \times p} containing the
#'   simulated values and injected missing values.
#' * `col_group`: a data frame with \eqn{p} rows mapping each `feature` to a
#'   `group`.
#' * `row_group`: a data frame with \eqn{n} rows mapping each `sample` to a
#'   `group`.
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
  checkmate::assert_number(perc_total_na, lower = 0, upper = (n - 1) / n)
  checkmate::assert_number(perc_col_na, lower = 0, upper = 1)
  checkmate::assert_flag(beta)

  # NA feasibility checks
  n_total_na <- ceiling(perc_total_na * n * p)
  n_cols_na <- max(1L, floor(perc_col_na * p))
  max_per_col <- n - 1L

  if (n_total_na > n_cols_na * max_per_col) {
    needed <- ceiling(n_total_na / max_per_col)
    cli::cli_abort(
      c(
        "{.arg perc_col_na} is too small for the requested {.arg perc_total_na}.",
        "i" = "Requested {n_total_na} NA{cli::qty(n_total_na)}{?s} across {n_cols_na} column{cli::qty(n_cols_na)}{?s} (max {max_per_col} per column).",
        "x" = "Need at least {needed} NA-bearing column{cli::qty(needed)}{?s}.",
        ">" = "Increase {.arg perc_col_na} so that at least {needed} of {p} columns can carry NAs (currently {n_cols_na}), or decrease {.arg perc_total_na}."
      )
    )
  }

  # correlated columns: var = rho + (1 - rho) = 1, cor = rho
  if (rho > 0) {
    common <- stats::rnorm(n) * sqrt(rho)
    Z <- matrix(stats::rnorm(n * p), n, p) * sqrt(1 - rho)
    d <- common %*% matrix(1, 1, p) + Z
  } else {
    d <- matrix(stats::rnorm(n * p), n, p)
  }

  # re scale
  if (beta) {
    mm <- col_min_max(d)
    mins <- mm[1, ]
    ranges <- mm[2, ] - mins
    ranges[ranges == 0] <- 1
    d <- sweep(d, 2, mins, "-")
    d <- sweep(d, 2, ranges, "/")
  }

  # metadata
  col_groups <- sample(
    paste0("group", seq_len(n_col_groups)),
    size = p,
    replace = (n_col_groups < p)
  )
  col_group <- data.frame(
    feature = paste0("feature", seq_len(p)),
    group = col_groups
  )
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

  # inject NAs: evenly divide across selected columns, cap at max_per_col
  if (n_total_na > 0) {
    cols_na <- sample.int(p, n_cols_na)
    per_col <- floor(n_total_na / n_cols_na)
    remainder <- n_total_na - per_col * n_cols_na
    counts <- rep(per_col, n_cols_na) + c(rep(1L, remainder), rep(0L, n_cols_na - remainder))
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
