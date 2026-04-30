#' Column Mean Imputation
#'
#' Impute missing values in a matrix by replacing them with the mean of their
#' respective columns.
#'
#' @param obj A numeric matrix.
#' @param subset Optional character or integer vector specifying columns to
#'   impute. If `NULL`, all columns are imputed.
#' @param cores Integer. Number of cores to use for parallel computation.
#'   Defaults to `1`.
#'
#' @returns A numeric matrix of the same dimensions as `obj`, with missing
#' values in the selected columns replaced by column means.
#'
#' @details
#' Columns with no observed values cannot be imputed by their column mean and
#' are left unchanged.
#'
#' @examples
#' obj <- matrix(c(1, 2, NA, 4, NA, 6, NA, 8, 9, NA, NA, NA), nrow = 3)
#' colnames(obj) <- c("A", "B", "C", "D")
#' obj
#'
#' # impute missing values with column means
#' mean_imp_col(obj)
#'
#' # impute only specific columns by name
#' mean_imp_col(obj, subset = c("A", "C"))
#'
#' # impute only specific columns by index
#' mean_imp_col(obj, subset = c(1, 3))
#'
#' @export
mean_imp_col <- function(obj, subset = NULL, cores = 1) {
  # pre-cond
  checkmate::assert_matrix(obj, mode = "numeric", .var.name = "obj")
  checkmate::assert_count(cores, positive = TRUE)

  # subset resolution
  subset <- resolve_subset(subset, obj)
  if (is.null(subset)) {
    return(obj)
  }

  # imputation
  res <- mean_imp_col_internal(mat = obj, col_idx = as.integer(subset) - 1L, cores = cores)

  dimnames(res) <- dimnames(obj)

  return(res)
}
