#' Matrix Column Variances
#'
#' Compute the sample variance for each column of a numeric matrix.
#'
#' @param obj A double matrix.
#'
#' @details
#' `NA` and `NaN` are treated as missing and dropped, as [stats::var()] does
#' with `na.rm = TRUE`. A column with fewer than two non-missing values has
#' variance `NA`.
#'
#' `Inf` and `-Inf` are not missing. A column that contains one and has at
#' least two non-missing values has variance `NaN`, as in base R.
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
col_vars <- function(obj) {
  checkmate::assert_matrix(
    obj,
    mode = "double",
    null.ok = FALSE,
    min.rows = 1,
    min.cols = 1,
    .var.name = "obj"
  )
  vars <- col_scan(obj, var = TRUE)$var
  names(vars) <- colnames(obj)
  return(vars)
}
