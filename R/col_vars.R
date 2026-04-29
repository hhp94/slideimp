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
#' @returns A numeric vector of column variances, named when `obj` has column
#' names.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' obj <- matrix(rnorm(4 * 10), ncol = 4)
#' obj[1, 1] <- NA
#' obj[1:8, 2] <- NA
#' obj[1:9, 3] <- NA
#' obj[, 4] <- NA
#' obj
#'
#' col_vars(obj)
#' apply(obj, 2, var, na.rm = TRUE)
col_vars <- function(obj, cores = 1) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric", null.ok = FALSE, min.rows = 1, min.cols = 1, .var.name = "obj"
  )
  checkmate::assert_int(cores, lower = 1)
  vars <- col_vars_internal(mat = obj, cores = cores)[1, ]
  vars[is.nan(vars)] <- NA
  names(vars) <- colnames(obj)
  return(vars)
}
