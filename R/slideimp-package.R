#' slideimp: Numeric Matrices K-NN and PCA Imputation
#'
#' Provides K-nearest-neighbor, PCA, grouped, and sliding-window imputation
#' methods for numeric matrices, along with utilities for simulation and
#' parameter tuning.
#'
#' @section Column mean fill:
#'
#' A column with a proportion of missing values greater than `colmax` takes no
#' part in the imputation. The method neither uses the column to impute other
#' columns nor imputes it.
#'
#' With `post_imp = TRUE`, column means replace the missing values that remain
#' after the method has run. Two kinds of missing value can remain: the missing
#' values of the columns that took no part, and the cells for which the method
#' found no value. Each column mean is the mean of the column's observed
#' values. With `post_imp = FALSE`, both kinds stay missing.
#'
#' The result records what the column means replaced, in two attributes:
#'
#' * `mean_filled_cols`: an integer vector with the index, in the result, of
#'   each column that took no part and was filled.
#' * `mean_filled_cells`: an integer matrix with the columns `row` and `col`,
#'   one row for each cell that the method could not impute and that was
#'   filled.
#'
#' Both attributes are empty when no column mean was used. `print()` reports
#' them when they are not empty.
#'
#' @section Missing values and non-finite input:
#'
#' `obj` must be a double matrix. An integer or logical matrix is rejected,
#' not converted. Convert one with `storage.mode(obj) <- "double"`.
#'
#' `NA` and `NaN` are both treated as missing. The package does not
#' distinguish between them.
#'
#' `Inf` and `-Inf` are rejected on entry, and so is a column with no
#' observed values (every cell `NA` or `NaN`). Locate them with
#' `which(is.infinite(obj), arr.ind = TRUE)` and
#' `which(colSums(!is.na(obj)) == 0)`.
#'
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang abort
#' @importFrom stats update
#' @useDynLib slideimp, .registration = TRUE
"_PACKAGE"
