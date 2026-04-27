#' Calculate Matrix Column Variances
#'
#' Computes the sample variance for each column of a numeric matrix.
#'
#' @param mat A numeric matrix.
#' @param cores Integer. Number of cores to use for parallel computation.
#'   Defaults to `1`.
#'
#' @details
#' Columns with fewer than two distinct non-missing values are assigned `NA`.
#'
#' @returns A numeric vector of column variances, named when `mat` has column
#' names.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' mat <- matrix(rnorm(4 * 10), ncol = 4)
#' mat[1, 1] <- NA
#' mat[1:8, 2] <- NA
#' mat[1:9, 3] <- NA
#' mat[, 4] <- NA
#' mat
#'
#' col_vars(mat)
#' apply(mat, 2, var, na.rm = TRUE)
col_vars <- function(mat, cores = 1) {
  checkmate::assert_matrix(
    mat,
    mode = "numeric", null.ok = FALSE, min.rows = 1, min.cols = 1, .var.name = "mat"
  )
  checkmate::assert_int(cores, lower = 1)
  vars <- col_vars_internal(mat = mat, cores = cores)[1, ]
  vars[is.nan(vars)] <- NA
  names(vars) <- colnames(mat)
  return(vars)
}
