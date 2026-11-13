#' slideimp: Numeric Matrices K-NN and PCA Imputation
#'
#' Provides K-nearest-neighbor, PCA, grouped, and sliding-window imputation
#' methods for numeric matrices, along with utilities for simulation and
#' parameter tuning.
#'
#' @section Column mean fill:
#'
#' A column with a proportion of missing values greater than `colmax` takes no
#' part in the imputation. The method does not use the column to impute other
#' columns, and the method does not impute the column.
#'
#' With `post_imp = TRUE`, column means replace the missing values that remain
#' after the method has run. Two kinds of missing value remain: the missing
#' values of the columns that took no part, and the cells for which the method
#' found no value. For such a cell, the column mean includes the values that
#' the method imputed in the same column. With `post_imp = FALSE`, both kinds
#' stay missing.
#'
#' The result records what the column means replaced, in two attributes:
#'
#' * `mean_filled_cols`: an integer vector with the index of each column that
#'   took no part and was filled.
#' * `mean_filled_cells`: an integer matrix with the columns `row` and `col`,
#'   one row for each cell that the method could not impute and that was
#'   filled.
#'
#' Both attributes are empty when no column mean was used. `print()` reports
#' them when they are not empty.
#'
#' @section Missing values and non-finite input:
#'
#' `NA` and `NaN` are both treated as missing. No function in the package
#' distinguishes between them.
#'
#' `Inf` and `-Inf` are rejected on entry and are never imputed. Locate them
#' with `which(is.infinite(x), arr.ind = TRUE)`.
#'
#' A column with no observed values (every cell `NA` or `NaN`) is also
#' rejected.
#'
#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @importFrom stats update
#' @useDynLib slideimp, .registration = TRUE
"_PACKAGE"
