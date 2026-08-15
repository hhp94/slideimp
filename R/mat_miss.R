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
#' @returns When `prop = FALSE`, an integer vector of missing-value counts.
#'   When `prop = TRUE`, a double vector of missing-value proportions. Either
#'   is per column when `col = TRUE` and per row otherwise, and is named when
#'   the corresponding dimension names are present.
#'
#'   `NA` and `NaN` both count as missing. `Inf` and `-Inf` do not.
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
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    null.ok = FALSE,
    .var.name = "obj"
  )
  checkmate::assert_flag(col, .var.name = "col")
  checkmate::assert_flag(prop, .var.name = "prop")
  # counts are bounded by nrow/ncol, so they fit in an integer
  if (col) {
    vec_miss <- as.integer(col_miss_internal(obj))
    names(vec_miss) <- colnames(obj)
    denom <- nrow(obj)
  } else {
    vec_miss <- as.integer(row_miss_internal(obj))
    names(vec_miss) <- rownames(obj)
    denom <- ncol(obj)
  }
  if (prop) {
    vec_miss <- vec_miss / denom
  }
  vec_miss
}
