#' Column or Row Missing Counts and Proportions
#'
#' Calculate the number or proportion of missing values per column or per row
#' of a numeric matrix, without allocating a logical matrix the size of `obj`.
#'
#' @param obj A double matrix.
#' @param col Logical. Whether to compute per column (`TRUE`) or per row
#'   (`FALSE`).
#' @param prop Logical. Whether to return missing-value proportions (`TRUE`)
#'   or counts (`FALSE`).
#'
#' @details `NA` and `NaN` both count as missing. `Inf` and `-Inf` do not.
#'
#' @returns When `prop = FALSE`, an integer vector of missing-value counts.
#'   When `prop = TRUE`, a double vector of missing-value proportions. The
#'   vector has one element per column when `col = TRUE`, or per row
#'   otherwise, named from `colnames(obj)` or `rownames(obj)` when set.
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
    mode = "double",
    null.ok = FALSE,
    .var.name = "obj"
  )
  checkmate::assert_flag(col, .var.name = "col")
  checkmate::assert_flag(prop, .var.name = "prop")
  if (col) {
    vec_miss <- col_scan(obj)$n_miss
    names(vec_miss) <- colnames(obj)
    denom <- nrow(obj)
  } else {
    vec_miss <- col_scan(obj, rows = TRUE)$row_miss
    names(vec_miss) <- rownames(obj)
    denom <- ncol(obj)
  }
  if (prop) {
    vec_miss <- vec_miss / denom
  }
  vec_miss
}
