#' Column Mean Imputation
#'
#' Impute missing values in a matrix by replacing them with the mean of their
#' respective columns.
#'
#' @inheritParams knn_imp
#'
#' @returns A numeric matrix of the same dimensions as `obj` with missing values
#' in the specified columns replaced by column means.
#'
#' @examples
#' # Create example matrix with missing values
#' mat <- matrix(c(1, 2, NA, 4, NA, 6, NA, 8, 9), nrow = 3)
#' colnames(mat) <- c("A", "B", "C")
#' mat
#'
#' # Impute missing values with column means
#' imputed_mat <- mean_imp_col(mat)
#' imputed_mat
#'
#' # Impute only specific columns by name
#' imputed_subset <- mean_imp_col(mat, subset = c("A", "C"))
#' imputed_subset
#'
#' # Impute only specific columns by index
#' imputed_idx <- mean_imp_col(mat, subset = c(1, 3))
#' imputed_idx
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
