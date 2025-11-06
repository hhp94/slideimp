#' Calculate Matrix Column Variance
#'
#' Computes the sample variance for each column of a numeric matrix
#'
#' @param mat A numeric matrix.
#' @param cores Number of cores to use for parallel computation. Defaults to 1.
#'
#' @return A named numeric vector of column variances. Variances for columns with
#'   insufficient variation (e.g., all identical values) are set to \code{NA}.
#'
#' @export
#'
#' @examples
#' col_vars(t(khanmiss1))
col_vars <- function(mat, cores = 1) {
  checkmate::assert_matrix(
    mat,
    mode = "numeric", null.ok = FALSE, min.rows = 1, min.cols = 1, .var.name = "mat"
  )
  vars <- col_vars_internal(mat = mat, cores = cores)[1, ]
  vars[is.nan(vars)] <- NA
  names(vars) <- colnames(mat)
  return(vars)
}
