pca_imp_internal <- function(
    X, miss, ncp, scale, method, ind.sup, quanti.sup, threshold, seed, init, maxiter,
    row.w, coeff.ridge, nrX, ncX) {
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

    svd.res <- FactoMineR::svd.triplet(Xhat, row.w = row.w, ncp = ncp)
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
      if ((criterion < threshold) && (nb.iter > 5)) {
        nb.iter <- 0
      }
      if ((objective < threshold) && (nb.iter > 5)) {
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

pca_imp <- function(X, ncp = 2, scale = TRUE,
                    method = c("Regularized", "EM"), row.w = NULL,
                    ind.sup = NULL, quanti.sup = NULL,
                    coeff.ridge = 1, threshold = 1e-6, seed = NULL,
                    nb.init = 1, maxiter = 1000, ...) {
  #### Main program
  if (!anyNA(X)) {
    return(X)
  }

  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"),
    several.ok = T
  )[1]
  obj <- Inf
  method <- tolower(method)

  nrX <- nrow(X) - length(ind.sup)
  ncX <- ncol(X) - length(quanti.sup)

  if (ncp > min(nrow(X) - 2, ncol(X) - 1)) {
    stop("ncp is too large")
  }

  if (is.null(row.w)) {
    row.w <- rep(1, nrow(X)) / nrow(X)
  }

  if (!is.null(ind.sup)) {
    row.w[ind.sup] <- row.w[ind.sup] * 1e-08
  }

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  miss <- is.na(X)
  for (i in seq_len(nb.init)) {
    res.impute <- pca_imp_internal(
      X,
      miss = miss, ncp = ncp, scale = scale, method = method,
      ind.sup = ind.sup, quanti.sup = quanti.sup, threshold = threshold,
      seed = if (!is.null(seed)) {
        (seed * (i - 1))
      } else {
        NULL
      },
      init = i, maxiter = maxiter, row.w = row.w,
      coeff.ridge = coeff.ridge, nrX = nrX, ncX = ncX
    )

    cur_obj <- mean((res.impute$fittedX[!miss] - X[!miss])^2)

    if (cur_obj < obj) {
      res <- res.impute
      obj <- cur_obj
    }
  }

  return(res)
}
