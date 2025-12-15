#' @importFrom bootSVD fastSVD
fast.svd.triplet <- function(X, ncp) {
  svd.usuelle <- suppressWarnings(fastSVD(X))
  U <- svd.usuelle$u[, 1:ncp, drop = FALSE]
  V <- svd.usuelle$v[, 1:ncp, drop = FALSE]
  if (ncp > 1) {
    mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
    mult[mult == 0] <- 1
    U <- t(t(U) * mult)
    V <- t(t(V) * mult)
  }
  vs <- svd.usuelle$d[1:min(ncol(X), nrow(X) - 1)]
  # num <- which(vs[1:ncp] < 1e-15)
  #   if (length(num) == 1) {
  #       U[, num] <- U[, num, drop = FALSE] * vs[num]
  #       V[, num] <- V[, num, drop = FALSE] * vs[num]
  #   }
  #   if (length(num) > 1) {
  #       U[, num] <- t(t(U[, num]) * vs[num])
  #       V[, num] <- t(t(V[, num]) * vs[num])
  #   }
  res <- list(vs = vs, U = U, V = V)
  return(res)
}

#' @importFrom stats rnorm
pca_imp_internal <- function(
  X, miss, ncp, scale, method, ind.sup, quanti.sup, threshold, seed, init, maxiter,
  miniter, row.w, coeff.ridge, nrX, ncX
) {
  nb.iter <- 1
  old <- Inf
  objective <- 0
  if (!is.null(seed)) {
    set.seed(seed)
  }
  ncp <- min(ncp, ncol(X), nrow(X) - 1)
  missing <- which(miss)
  # X has missing value, so we temporary replace the missing value with zero to
  # allow for vectorized calculation instead of apply + moy.p/ec
  X[missing] <- 0
  denom <- as.vector(row.w %*% !miss)
  weighted_sum <- as.vector(row.w %*% X)
  sum_sq <- as.vector(row.w %*% (X * X))
  mean.p <- weighted_sum / denom
  et <- sqrt(sum_sq / denom - mean.p^2)
  # Then we refill back the missing value with NA to calculate Xhat. This operation
  # is cheap
  X[missing] <- NA
  Xhat <- t(t(X) - mean.p)
  if (scale) {
    Xhat <- t(t(Xhat) / et)
  }
  if (!is.null(quanti.sup)) {
    Xhat[, quanti.sup] <- Xhat[, quanti.sup] * 1e-08
  }
  # Init with 0 for the first iteration, which is just mean init. But from here on,
  # Xhat have no more missing values because it is now either filled with zero or
  # rnorm value.
  Xhat[missing] <- 0
  if (init > 1) {
    Xhat[missing] <- rnorm(length(missing)) # random initialization
  }
  fittedX <- Xhat
  total_w <- sum(row.w)
  if (ncp == 0) {
    nb.iter <- 0
  }
  while (nb.iter > 0) {
    Xhat[missing] <- fittedX[missing]
    if (!is.null(quanti.sup)) {
      Xhat[, quanti.sup] <- Xhat[, quanti.sup] * 1e+08
    }
    if (scale) {
      Xhat <- t(t(Xhat) * et)
    }
    Xhat <- t(t(Xhat) + mean.p)
    mean.p <- as.vector(row.w %*% Xhat / total_w)
    Xhat <- t(t(Xhat) - mean.p)
    et <- sqrt(as.vector(row.w %*% (Xhat * Xhat) / total_w))
    if (scale) {
      Xhat <- t(t(Xhat) / et)
    }
    if (!is.null(quanti.sup)) {
      Xhat[, quanti.sup] <- Xhat[, quanti.sup] * 1e-08
    }

    svd.res <- fast.svd.triplet(Xhat, ncp = ncp)
    sigma2 <- nrX * ncX / min(ncX, nrX - 1) *
      sum((svd.res$vs[-c(seq_len(ncp))]^2) /
        ((nrX - 1) * ncX - (nrX - 1) * ncp - ncX * ncp + ncp^2))
    sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp + 1]^2)
    if (method == "em") {
      sigma2 <- 0
    }
    lambda.shrinked <- (svd.res$vs[seq_len(ncp)]^2 - sigma2) / svd.res$vs[seq_len(ncp)]
    fittedX <- tcrossprod(
      t(t(svd.res$U[, seq_len(ncp), drop = FALSE] * row.w) * lambda.shrinked),
      svd.res$V[, seq_len(ncp), drop = FALSE]
    )
    fittedX <- fittedX / row.w
    diff <- Xhat - fittedX
    diff[missing] <- 0
    objective <- sum(diff^2 * row.w)
    criterion <- abs(1 - objective / old)
    old <- objective
    nb.iter <- nb.iter + 1
    if (!is.nan(criterion)) {
      if ((criterion < threshold) && (nb.iter > miniter)) {
        nb.iter <- 0
      }
      if ((objective < threshold) && (nb.iter > miniter)) {
        nb.iter <- 0
      }
    }
    if (nb.iter > maxiter) {
      nb.iter <- 0
      warning(paste("Stopped after ", maxiter, " iterations"))
    }
  }

  if (!is.null(quanti.sup)) {
    Xhat[, quanti.sup] <- Xhat[, quanti.sup] * 1e+08
  }

  if (scale) {
    Xhat <- t(t(Xhat) * et)
  }
  Xhat <- t(t(Xhat) + mean.p)
  completeObs <- X
  completeObs[missing] <- Xhat[missing]
  if (!is.null(quanti.sup)) {
    fittedX[, quanti.sup] <- fittedX[, quanti.sup] * 1e+08
  }
  if (scale) {
    fittedX <- t(t(fittedX) * et)
  }
  fittedX <- t(t(fittedX) + mean.p)
  result <- list()
  result$completeObs <- completeObs
  result$fittedX <- fittedX
  return(result)
}

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
#' @param method "Regularized" by default or "EM"
#' @param coeff.ridge 1 by default to perform the regularized pca_imp (imputePCA) algorithm; useful only if method="Regularized". Other regularization terms can be implemented by setting the value to less than 1 in order to regularized less (to get closer to the results of the EM method) or more than 1 to regularized more (to get closer to the results of the mean imputation)
#' @param threshold the threshold for assessing convergence
#' @param seed integer, by default seed = NULL implies that missing values are initially imputed by the mean of each variable. Other values leads to a random initialization
#' @param nb.init integer corresponding to the number of random initializations; the first initialization is the initialization with the mean imputation
#' @param maxiter integer, maximum number of iteration for the algorithm
#' @param miniter integer, minimum number of iteration for the algorithm
#'
#' @return A `dim(obj)` matrix with missing values imputed.
#'
#' @references
#' Josse, J & Husson, F. (2013). Handling missing values in exploratory multivariate data analysis methods. Journal de la SFdS. 153 (2), pp. 79-99.
#'
#' Josse, J. and Husson, F. missMDA (2016). A Package for Handling Missing Values in Multivariate Data Analysis. Journal of Statistical Software, 70 (1), pp 1-31 \doi{doi:10.18637/jss.v070.i01}.
#'
#' @author Francois Husson \email{francois.husson@institut-agro.fr} and Julie Josse \email{julie.josse@polytechnique.edu}
#'
#' \href{https://www.youtube.com/watch?v=YDbx2pk9xNY&list=PLnZgp6epRBbQzxFnQrcxg09kRt-PA66T_&index=2}{Video showing how to perform PCA on an incomplete dataset}
#'
#' @examples
#' data("khanmiss1")
#'
#' # Transpose to put genes on columns
#' pca_imp(t(khanmiss1), ncp = 2)
#'
#' @export
pca_imp <- function(
  obj, ncp = 2, scale = TRUE, method = c("Regularized", "EM"),
  coeff.ridge = 1, threshold = 1e-6, seed = NULL,
  nb.init = 1, maxiter = 1000, miniter = 5
) {
  #### Main program
  row.w <- NULL
  ind.sup <- NULL
  quanti.sup <- NULL

  checkmate::assert_matrix(obj, mode = "numeric", row.names = "named", col.names = "unique", null.ok = FALSE, .var.name = "obj")
  cn <- colnames(obj)

  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"),
    several.ok = T
  )[1]

  checkmate::assert_flag(scale, .var.name = "scale")
  checkmate::assert_int(ncp, lower = 1, upper = min(ncol(obj), nrow(obj)), .var.name = "ncp")
  checkmate::assert_number(coeff.ridge, .var.name = "coeff.ridge")
  checkmate::assert_number(seed, null.ok = TRUE, .var.name = "seed")
  # checkmate::assert_numeric(row.w, lower = 0, upper = 1, any.missing = FALSE, len = nrow(obj), null.ok = TRUE, .var.name = "row.w")
  # checkmate::assert_integerish(ind.sup, lower = 1, upper = nrow(obj), any.missing = FALSE, len = nrow(obj), null.ok = TRUE, .var.name = "ind.sup")
  checkmate::assert_int(nb.init, lower = 1, .var.name = "nb.init")
  checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
  checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")
  obj_vars <- col_vars(obj)
  if (
    any(abs(obj_vars) < .Machine$double.eps^0.5 | is.nan(obj_vars) | is.na(obj_vars))
  ) {
    stop("Features with zero variance after na.rm not permitted for PCA Imputation. Try `col_vars(obj)`")
  }
  if (!anyNA(obj)) {
    return(obj)
  }

  init_obj <- Inf
  method <- tolower(method)

  nrX <- nrow(obj) - length(ind.sup)
  ncX <- ncol(obj) - length(quanti.sup)

  if (ncp > min(nrow(obj) - 2, ncol(obj) - 1)) {
    stop("ncp is too large")
  }

  if (is.null(row.w)) {
    row.w <- rep(1, nrow(obj)) / nrow(obj)
  }

  if (!is.null(ind.sup)) {
    row.w[ind.sup] <- row.w[ind.sup] * 1e-08
  }

  miss <- is.na(obj)
  cmiss <- colSums(miss)
  if (any(cmiss / nrow(obj) == 1)) {
    stop("Col(s) with all missing detected. Remove before proceed")
  }
  for (i in seq_len(nb.init)) {
    res.impute <- pca_imp_internal(
      X = obj,
      miss = miss, ncp = ncp, scale = scale, method = method,
      ind.sup = ind.sup, quanti.sup = quanti.sup, threshold = threshold,
      seed = if (!is.null(seed)) {
        (seed * (i - 1))
      } else {
        NULL
      },
      init = i, maxiter = maxiter, miniter = miniter, row.w = row.w,
      coeff.ridge = coeff.ridge, nrX = nrX, ncX = ncX
    )

    cur_obj <- mean((res.impute$fittedX[!miss] - obj[!miss])^2)

    if (cur_obj < init_obj) {
      res <- res.impute
      init_obj <- cur_obj
    }
  }

  class(res$completeObs) <- c("ImputedMatrix", class(res$completeObs))
  attr(res$completeObs, "imp_method") <- "pca"
  return(res$completeObs)
}
