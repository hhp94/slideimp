#' Column Mean Imputation
#'
#' @description
#' Imputes missing values in a matrix by replacing them with the mean of their
#' respective columns.
#'
#' @details
#' This function calculates the mean for each column excluding missing values and
#' replaces all missing values in that column with the computed mean.
#'
#' The `subset` parameter allows for imputation of only specific columns.
#'
#' @inheritParams knn_imp
#'
#' @return A numeric matrix of the same dimensions as `obj` with missing values
#'   in the specified columns replaced by column means.
#'
#' @examples
#' # Create example matrix with missing values
#' mat <- matrix(c(1, 2, NA, 4, 5, 6, NA, 8, 9), nrow = 3)
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
mean_imp_col <- function(obj, subset = NULL) {
  checkmate::assert_matrix(obj, mode = "numeric", .var.name = "obj")
  checkmate::assert(
    checkmate::check_character(subset, min.len = 0, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
    checkmate::check_integerish(subset, lower = 1, upper = ncol(obj), min.len = 0, any.missing = FALSE, null.ok = TRUE, unique = TRUE),
    combine = "or",
    .var.name = "subset"
  )

  if (!is.null(subset)) {
    if (length(subset) == 0) {
      return(obj)
    }
    if (is.character(subset)) {
      if (is.null(colnames(obj))) {
        stop("`subset` contains characters but `obj` doesn't have column names")
      }
      matched <- match(subset, colnames(obj), nomatch = NA)
      if (any(is.na(matched))) {
        warning("Some subset names not found in column names and will be ignored")
      }
      subset <- matched[!is.na(matched)]
      if (length(subset) == 0) {
        warning("No valid column names found in subset")
        return(obj)
      }
    }
    obj_subset <- obj[, subset, drop = FALSE]
    na_indices <- which(is.na(obj_subset), arr.ind = TRUE)
    if (nrow(na_indices) > 0) {
      column_means <- colMeans(obj_subset, na.rm = TRUE)
      obj_subset[na_indices] <- column_means[na_indices[, 2]]
      obj[, subset] <- obj_subset
    }
  } else {
    na_indices <- which(is.na(obj), arr.ind = TRUE)
    if (nrow(na_indices) > 0) {
      column_means <- colMeans(obj, na.rm = TRUE)
      obj[na_indices] <- column_means[na_indices[, 2]]
    }
  }
  return(obj)
}
