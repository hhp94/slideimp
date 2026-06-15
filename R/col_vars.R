#' Calculate Matrix Column Variances
#'
#' Compute the sample variance for each column of a numeric matrix.
#'
#' @param obj A numeric matrix.
#' @param cores Integer. Number of cores to use for parallel computation.
#'   Defaults to `1`.
#'
#' @details
#' Columns with fewer than two non-missing values are assigned `NA`.
#'
#' `NA`/`NaN` are treated as missing and dropped (equivalent to
#' `var(x, na.rm = TRUE)`). `Inf`/`-Inf` are **not** missing. They enter the
#' arithmetic, so a column's variance can be `NaN`, matching base R.
#'
#' @returns A numeric vector of column variances, named if `obj` has column
#'   names.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- matrix(rnorm(7 * 10), ncol = 7)
#' obj[1, 1] <- Inf
#' obj[1, 2] <- NA
#' obj[1:8, 3] <- NA
#' obj[8, 3] <- Inf
#' obj[1:8, 4] <- NA
#' obj[1:8, 5] <- NA
#' obj[9, 5] <- obj[10, 5]
#' obj[1:9, 6] <- NA
#' obj[, 7] <- NA
#' obj
#'
#' col_vars(obj)
#' apply(obj, 2, var, na.rm = TRUE)
col_vars <- function(obj, cores = 1) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    null.ok = FALSE,
    min.rows = 1,
    min.cols = 1,
    .var.name = "obj"
  )
  checkmate::assert_int(cores, lower = 1)
  vars <- col_vars_internal(mat = obj, cores = cores)[1, ]
  names(vars) <- colnames(obj)
  return(vars)
}
