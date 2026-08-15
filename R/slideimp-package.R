#' slideimp: Numeric Matrices K-NN and PCA Imputation
#'
#' Provides K-nearest-neighbor, PCA, grouped, and sliding-window imputation
#' methods for numeric matrices, along with utilities for simulation and
#' parameter tuning.
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
#' @useDynLib slideimp, .registration = TRUE
"_PACKAGE"
