#' Column or Row Missing Counts and Proportions
#'
#' Calculate the number or proportion of missing values per column or per row
#' of a numeric matrix without allocating a full logical mask matrix.
#'
#' @param obj A numeric matrix.
#' @param col Logical. If `TRUE`, compute column-wise. If `FALSE`,
#'   compute row-wise.
#' @param prop Logical. If `FALSE`, return missing-value counts. If `TRUE`,
#'   return missing-value proportions.
#'
#' @returns A numeric vector containing missing-value counts or proportions for
#'   columns or rows, named when the corresponding dimension names are present.
#'
#' @examples
#' obj <- matrix(c(1, NA, 3, 4, NA, 6, NA, 8, 9), nrow = 3)
#' obj
#'
#' # column missing counts
#' mat_miss(obj)
#'
#' # row missing counts
#' mat_miss(obj, col = FALSE)
#'
#' # column missing proportions
#' mat_miss(obj, prop = TRUE)
#'
#' @export
mat_miss <- function(obj, col = TRUE, prop = FALSE) {
  checkmate::assert_matrix(obj, mode = "numeric", null.ok = FALSE, .var.name = "obj")
  checkmate::assert_flag(col, .var.name = "col")
  checkmate::assert_flag(prop, .var.name = "prop")
  if (col) {
    vec_miss <- as.numeric(col_miss_internal(obj))
    names(vec_miss) <- colnames(obj)
    denom <- nrow(obj)
  } else {
    vec_miss <- as.numeric(row_miss_internal(obj))
    names(vec_miss) <- rownames(obj)
    denom <- ncol(obj)
  }
  if (prop) {
    vec_miss <- vec_miss / denom
  }
  vec_miss
}
