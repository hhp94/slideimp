#' Impute dataset with PCA
#'
#' (From the missMDA package on CRAN) Impute the missing values of a dataset with the Principal Components Analysis model. Can be used as a preliminary step before performing a PCA on an completed dataset.
#'
#' @details
#' Impute the missing entries of a mixed data using the iterative PCA algorithm (method="EM") or the regularised iterative PCA algorithm (method="Regularized"). The (regularized) iterative PCA algorithm first consists imputing missing values with initial values such as the mean of the variable. If the argument seed is set to a specific value, a random initialization is performed: the initial values are drawn from a gaussian distribution
#' with mean and standard deviation calculated from the observed values. nb.init different random initialization can be drawn. In such a situation, the solution giving the smallest objective function (the mean square error between the fitted matrix and the observed one) is kept. The second step of the (regularized) iterative PCA algorithm is to perform PCA on the completed dataset. Then, it imputes the missing values with the (regularized) reconstruction formulae of order ncp (the fitted matrix computed with ncp components for the (regularized) scores and loadings). These steps of estimation of the parameters via PCA and imputation of the missing values using the (regularized) fitted matrix are iterate until convergence. The iterative PCA algorithm is also known as the EM-PCA algorithm since it corresponds to an EM algorithm of the fixed effect model where the data are generated as a fixed structure (with a low rank representation) corrupted by noise. The number of components used in the algorithm can be found using cross-validation criteria implemented in the function estim_ncpPCA.\cr
#' We advice to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. In the regularized algorithm, the singular values of the PCA are shrinked.\cr
#' The output of the algorithm can be used as an input of the PCA function of the FactoMineR package in order to perform PCA on an incomplete dataset.
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param ncp integer corresponding to the number of components used to to predict the missing entries
#' @param scale boolean. By default TRUE leading to a same weight for each variable
#' @param method "regularized" by default or "EM"
#' @param coeff.ridge 1 by default to perform the regularized pca_imp (imputePCA) algorithm; useful only if method="Regularized". Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of the EM method) or more than 1 to regularized more (to get closer to the results of the mean imputation)
#' @param row.w row weights (by default, a vector of 1 for uniform row weights)
#' @param threshold the threshold for assessing convergence
#' @param seed integer, by default seed = NULL implies that missing values are initially imputed by the mean of each variable. Other values leads to a random initialization
#' @param nb.init integer corresponding to the number of random initializations; the first initialization is the initialization with the mean imputation
#' @param maxiter integer, maximum number of iteration for the algorithm
#' @param miniter integer, minimum number of iteration for the algorithm
#'
#' @inherit knn_imp return
#'
#' @references
#' Josse, J & Husson, F. (2013). Handling missing values in exploratory multivariate data analysis methods. Journal de la SFdS. 153 (2), pp. 79-99.
#'
#' Josse, J. and Husson, F. missMDA (2016). A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70 (1), pp 1-31 \doi{doi:10.18637/jss.v070.i01}.
#'
#' @author Francois Husson \email{francois.husson@institut-agro.fr}
#'
#' @author Julie Josse \email{julie.josse@polytechnique.edu}
#'
#' @examples
#' data("khanmiss1")
#'
#' # Transpose to put genes on columns
#' pca_imp(t(khanmiss1), ncp = 2)
#'
#' @export
pca_imp <- function(
  obj, ncp = 2, scale = TRUE, method = c("regularized", "EM"),
  coeff.ridge = 1, row.w = NULL, threshold = 1e-6, seed = NULL,
  nb.init = 1, maxiter = 1000, miniter = 5
) {
  #### Main program
  checkmate::assert_matrix(obj, mode = "numeric", row.names = "named", col.names = "unique", null.ok = FALSE, .var.name = "obj")
  cn <- colnames(obj)
  method <- match.arg(method)
  checkmate::assert_flag(scale, .var.name = "scale")
  checkmate::assert_int(ncp, lower = 1, upper = min(ncol(obj), nrow(obj) - 1), .var.name = "ncp")
  checkmate::assert_number(coeff.ridge, lower = 0, .var.name = "coeff.ridge")
  checkmate::assert_number(seed, null.ok = TRUE, .var.name = "seed")
  checkmate::assert_numeric(row.w, finite = TRUE, any.missing = FALSE, len = nrow(obj), null.ok = TRUE, .var.name = "row.w")
  checkmate::assert_number(threshold, lower = 0, .var.name = "threshold")
  checkmate::assert_int(nb.init, lower = 1, .var.name = "nb.init")
  checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
  checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")
  obj_vars <- col_vars(obj)
  if (
    any(abs(obj_vars) < .Machine$double.eps^0.5 | is.na(obj_vars))
  ) {
    stop("Features with zero variance after na.rm not permitted for PCA Imputation. Try `col_vars(obj)`")
  }
  if (!anyNA(obj)) {
    return(obj)
  }

  init_obj <- Inf
  method <- tolower(method)

  nrX <- nrow(obj)
  ncX <- ncol(obj)

  if (is.null(row.w)) {
    row.w <- rep(1, nrX)
  }
  row.w <- row.w / sum(row.w)

  miss <- is.na(obj)
  cmiss <- colSums(miss)
  if (any(cmiss / nrX == 1)) {
    stop("Col(s) with all missing detected. Remove before proceed")
  }
  # pre-fill obj with 0, important. See comments in pca_imp_internal_cpp
  obj[miss] <- 0
  for (i in seq_len(nb.init)) {
    if (!is.null(seed)) {
      set.seed(seed * (i - 1))
    }
    res.impute <- pca_imp_internal_cpp(
      X = obj,
      miss = miss,
      ncp = ncp,
      scale = scale,
      regularized = (method == "regularized"),
      threshold = threshold,
      init = i,
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

  class(obj) <- c("ImputedMatrix", class(obj))
  attr(obj, "imp_method") <- "pca"
  return(obj)
}
