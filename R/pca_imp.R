#' Impute Numeric Matrix with PCA Imputation
#'
#' Impute missing values in a numeric matrix using (regularized) iterative PCA.
#'
#' @details
#' This algorithm is based on the original `missMDA::imputePCA` function and is
#' optimized for tall/wide numeric matrices.
#'
#' @inheritParams knn_imp
#' @param ncp Integer. Number of components used to predict the missing entries.
#' @param scale Logical. If `TRUE` (default), variables are scaled to have
#'   unit variance.
#' @param method Character. Either `"regularized"` (default) or `"EM"`.
#' @param coeff.ridge Numeric. Ridge regularization coefficient (default is 1).
#'   Only used if `method = "regularized"`. Values < 1 regularize less (closer
#'   to EM); values > 1 regularize more (closer to mean imputation).
#' @param row.w Row weights (internally normalized to sum to 1). Can be one of:
#'   * `NULL` (default): All rows weighted equally.
#'   * A numeric vector: Custom positive weights of length `nrow(obj)`.
#'   * `"n_miss"`: Rows with more missing values receive lower weight.
#' @param threshold Numeric. The threshold for assessing convergence.
#' @param seed Integer. Random number generator seed.
#' @param nb.init Integer. Number of random initializations. The first
#'   initialization is always mean imputation.
#' @param maxiter Integer. Maximum number of iterations for the algorithm.
#' @param miniter Integer. Minimum number of iterations for the algorithm.
#'
#' @inherit knn_imp return
#'
#' @references
#' Josse, J. & Husson, F. (2013). Handling missing values in exploratory
#' multivariate data analysis methods. Journal de la SFdS. 153 (2), pp. 79-99.
#'
#' Josse, J. and Husson, F. (2016). missMDA: A Package for Handling Missing
#' Values in Multivariate Data Analysis. Journal of Statistical Software,
#' 70 (1), pp 1-31. \doi{10.18637/jss.v070.i01}
#'
#' @author Francois Husson and Julie Josse (original `missMDA` implementation).
#'
#' @examples
#' data("khanmiss1")
#'
#' # Transpose to put genes on columns.
#' # Randomly initialize missing values 5 times (1st time is mean).
#' pca_imp(t(khanmiss1), ncp = 2, nb.init = 5)
#'
#' @export
pca_imp <- function(
  obj, ncp = 2, scale = TRUE, method = c("regularized", "EM"),
  coeff.ridge = 1, row.w = NULL, threshold = 1e-6, seed = NULL,
  nb.init = 1, maxiter = 1000, miniter = 5
) {
  #### Main program
  checkmate::assert_matrix(obj, mode = "numeric", null.ok = FALSE, .var.name = "obj")
  if (!anyNA(obj)) {
    return(obj)
  }
  checkmate::assert_int(ncp, lower = 1L, upper = min(nrow(obj) - 2L, ncol(obj) - 1L))
  method <- match.arg(method)
  checkmate::assert_flag(scale, .var.name = "scale")
  checkmate::assert_number(coeff.ridge, lower = 0, .var.name = "coeff.ridge")
  checkmate::assert(
    checkmate::check_numeric(
      row.w,
      finite = TRUE,
      any.missing = FALSE,
      len = nrow(obj),
      lower = 1e-10
    ),
    checkmate::check_choice(row.w, choices = "n_miss"),
    checkmate::check_null(row.w),
    .var.name = "row.w"
  )
  checkmate::assert_number(threshold, lower = 0, .var.name = "threshold")
  checkmate::assert_int(nb.init, lower = 1, .var.name = "nb.init")
  checkmate::assert_int(seed, null.ok = TRUE, lower = 0, .var.name = "seed")
  if (!is.null(seed) && nb.init > 1L && seed > 2147483647L / (nb.init - 1L)) {
    stop("`seed` too large")
  }
  checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
  checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")

  obj_vars <- col_vars(obj)
  if (any(obj_vars < .Machine$double.eps | is.na(obj_vars))) {
    stop("Features with zero variance after na.rm not permitted for PCA Imputation. Try 'col_vars(obj)'")
  }
  miss <- is.na(obj)

  if (any((colSums(miss) / nrow(obj)) == 1)) {
    stop("Col(s) with all missing detected. Remove before proceed")
  }
  init_obj <- Inf

  if (is.null(row.w)) {
    row.w <- rep(1, nrow(obj))
  } else if (is.character(row.w) && row.w == "n_miss") {
    n_miss_per_row <- rowSums(miss)
    row.w <- 1 - (n_miss_per_row / ncol(obj))
    row.w[row.w < 1e-10] <- 1e-10
  }

  # These pre-fill and scale weight are now handled in C++
  # obj[miss] <- 0
  # row.w <- row.w / sum(row.w)
  for (i in seq_len(nb.init)) {
    if (!is.null(seed)) {
      set.seed(seed * (i - 1L)) # exactly as missMDA does
    }
    res.impute <- pca_imp_internal_cpp(
      X = obj,
      miss = miss,
      ncp = ncp,
      scale = scale,
      regularized = (method == "regularized"),
      threshold = threshold,
      init = if (i == 1) 0 else i,
      maxiter = maxiter,
      miniter = miniter,
      row_w = row.w,
      coeff_ridge = coeff.ridge
    )
    cur_obj <- res.impute$mse
    if (cur_obj < init_obj) {
      best_imputed <- res.impute$imputed_vals
      init_obj <- cur_obj
    }
  }
  # apply best imputation
  obj[miss] <- best_imputed

  class(obj) <- c("slideimp_results", class(obj))
  attr(obj, "imp_method") <- "pca"
  return(obj)
}
