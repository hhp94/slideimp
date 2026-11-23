.resolve_transform <- function(transform, n) {
  msg <- paste0(
    "'transform' must be \"identity\", \"logit_biwhiten\", or a transform ",
    "object: list(name, forward, backward), or list(name, fit) for one that ",
    "is fitted first"
  )
  if (is.character(transform)) {
    if (length(transform) != 1L || is.na(transform)) {
      stop(msg, call. = FALSE)
    }
    tf <- switch(
      transform,
      identity = tf_identity(),
      logit_biwhiten = tf_biwhiten(tf_sv(n), rows = TRUE),
      NULL
    )
    if (is.null(tf)) {
      stop(msg, call. = FALSE)
    }
    return(tf)
  }
  ok <- is.list(transform) &&
    is.character(transform$name) &&
    length(transform$name) == 1L &&
    (is.function(transform$fit) ||
      (is.function(transform$forward) && is.function(transform$backward)))
  if (!ok) {
    stop(msg, call. = FALSE)
  }
  transform
}

#' Transforms for Bayesian PCA Imputation
#'
#' Transform `obj` prior to the Bayesian PCA imputation. Each draw is mapped
#' back to the scale of `obj` before it is returned.
#'
#' @details
#' **The two names.** `transform = "identity"` is `tf_identity()`, and
#' `transform = "logit_biwhiten"` is
#' `tf_biwhiten(tf_sv(nrow(obj)), rows = TRUE)`. Build a transform with these
#' functions to change one part of either, such as the logit without the
#' rescaling, or the columns rescaled and the rows not.
#'
#' **`tf_identity()`** fits the values as is.
#'
#' **`tf_sv()`** is for proportions such as DNAm beta values. Each value `y` is
#' squeezed to `(y * (n - 1) + 0.5) / n`, which moves `0` and `1` inside the
#' interval (Smithson and Verkuilen, 2006), and then taken to the logit scale.
#' Draws are mapped back and clamped to `[0, 1]`. Every observed value must lie
#' in `[0, 1]`.
#'
#' **`tf_biwhiten()`** rescales the columns, and with `rows = TRUE` the rows, of
#' the matrix on the scale of `base` to a common noise level before the fit, and
#' takes the scales back out of each draw. The scales are estimated once from
#' the observed values, with the column variances shrunk as in Smyth (2004). The
#' rescaling follows the idea of Landa, Zhang and Kluger (2022), with the scales
#' estimated differently.
#'
#' **A transform of your own** is a list with `name`, a string, and two
#' functions, `forward(y, ...)` and `backward(z, ...)`, where `backward()`
#' undoes `forward()` value by value. Both must take `...`, since the fit passes
#' them column and row positions. It can be the `base` of `tf_biwhiten()`.
#'
#' @returns A transform object, to pass as `transform`.
#'
#' @references Smithson M, Verkuilen J (2006). A better lemon squeezer?
#'   Maximum-likelihood regression with beta-distributed dependent variables.
#'   *Psychological Methods*, 11(1), 54-71. \doi{10.1037/1082-989X.11.1.54}
#'
#' Smyth GK (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. *Statistical Applications
#' in Genetics and Molecular Biology*, 3(1), 1-25. \doi{10.2202/1544-6115.1027}
#'
#' Landa B, Zhang TTCK, Kluger Y (2022). Biwhitening reveals the rank of a count
#' matrix. *SIAM Journal on Mathematics of Data Science*, 4(4), 1420-1446.
#' \doi{10.1137/21M1456807}
#'
#' @seealso [mipca()], whose `transform` argument these build.
#'
#' @examples
#' set.seed(1)
#' sim <- sim_mat(n = 50, p = 20, beta = TRUE)
#' obj <- sim$input
#'
#' # The same model as transform = "logit_biwhiten"
#' tf_logit_biwhiten <- tf_biwhiten(tf_sv(nrow(obj)), rows = TRUE)
#'
#' # The logit without the rescaling
#' tf_logit <- tf_sv(nrow(obj))
#'
#' # The columns rescaled, the rows not
#' tf_cols <- tf_biwhiten(tf_sv(nrow(obj)), rows = FALSE)
#' fit <- mipca(obj, rank = 2, transform = tf_cols, seed = 1)
#' fit
#' @export
tf_identity <- function() {
  list(name = "identity", forward = identity, backward = function(z, ...) z)
}

#' @rdname tf_identity
#' @param n Integer. The sample size of the squeeze, usually `nrow(obj)`. At
#'   least `2`.
#' @export
tf_sv <- function(n) {
  n <- .chk_count(n, "n")
  if (n < 2L) {
    stop("'n' must be at least 2", call. = FALSE)
  }
  nd <- as.double(n)
  list(
    name = paste0("sv-logit(n=", n, ")"),
    forward = function(y, ...) {
      if (!checkmate::test_numeric(y, lower = 0, upper = 1)) {
        stop(
          "this transform reads proportions, so every observed cell of ",
          "'obj' must lie in [0, 1]. For data on another scale, use ",
          "'transform = \"identity\"'",
          call. = FALSE
        )
      }
      stats::qlogis((y * (nd - 1) + 0.5) / nd)
    },
    backward = function(z, ...) {
      pmin(pmax((stats::plogis(z) * nd - 0.5) / (nd - 1), 0), 1)
    }
  )
}

.trigamma_inverse <- function(x) {
  y <- 0.5 + 1 / x
  for (iter in seq_len(50L)) {
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2L)
    y <- y + dif
    if (max(-dif / y) < 1e-8) {
      return(y)
    }
  }
  stop("internal error: the trigamma inverse did not converge")
}

.fit_fdist <- function(s2, df) {
  e <- log(s2) - digamma(df / 2) + log(df / 2)
  emean <- mean(e)
  evar <- if (length(e) > 1L) {
    sum((e - emean)^2) / (length(e) - 1) - mean(trigamma(df / 2))
  } else {
    -Inf
  }
  if (evar > 0) {
    d0 <- 2 * .trigamma_inverse(evar)
    s20 <- exp(emean + digamma(d0 / 2) - log(d0 / 2))
  } else {
    d0 <- Inf
    s20 <- exp(emean)
  }
  list(d0 = d0, s20 = s20)
}

.squeeze_var <- function(s2, df, d0, s20) {
  if (is.infinite(d0)) {
    return(rep(s20, length(s2)))
  }
  (d0 * s20 + df * s2) / (d0 + df)
}

.BIWHITEN_D_MIN <- 2L

#' @rdname tf_identity
#' @param base A transform object. The scale the rescaling is applied on:
#'   [tf_identity()], [tf_sv()] or a transform of your own, not another
#'   `tf_biwhiten()`. Defaults to `tf_identity()`.
#' @param rows Logical. If `TRUE`, rescale the rows as well as the columns.
#'   Defaults to `TRUE`.
#' @export
tf_biwhiten <- function(base = tf_identity(), rows = TRUE) {
  if (
    !is.list(base) || !is.function(base$forward) || !is.function(base$backward)
  ) {
    stop(
      "'base' must be a transform: list(name, forward, backward)",
      call. = FALSE
    )
  }
  if (is.function(base$fit)) {
    stop(
      "'base' must be a transform with no 'fit' step, such as tf_identity() ",
      "or tf_sv()",
      call. = FALSE
    )
  }
  if (!isTRUE(rows) && !isFALSE(rows)) {
    stop("'rows' must be TRUE or FALSE", call. = FALSE)
  }
  nm <- paste0(if (rows) "biwhiten(" else "whiten(", base$name, ")")
  list(
    name = nm,
    fit = function(X, rank) {
      rank <- .chk_count(rank, "rank")
      Z <- base$forward(X)
      obs <- !is.na(Z)
      n <- nrow(Z)
      p <- ncol(Z)
      n_obs <- colSums(obs)
      mu <- colSums(Z, na.rm = TRUE) / pmax(n_obs, 1L)
      Zf <- Z - rep(mu, each = n)
      rm(Z)
      Zf[!obs] <- 0
      e <- eigen(tcrossprod(Zf), symmetric = TRUE)
      W <- e$vectors[, seq_len(rank), drop = FALSE]
      rm(e)
      R <- Zf - W %*% crossprod(W, Zf)
      R[!obs] <- 0
      k <- pmax(0L, pmin(rank, n_obs - 1L - .BIWHITEN_D_MIN))
      full <- n_obs == n & k == rank
      lev_row <- drop(obs %*% (1 / pmax(n_obs, 1L))) + rowSums(W^2) * sum(full)
      for (j in which(!full)) {
        O <- which(obs[, j])
        if (k[j] == 0L) {
          R[O, j] <- Zf[O, j]
          next
        }
        Wo <- W[O, seq_len(k[j]), drop = FALSE]
        Gi <- solve(crossprod(Wo))
        b <- Gi %*% crossprod(Wo, Zf[O, j])
        R[O, j] <- Zf[O, j] - Wo %*% b
        lev_row[O] <- lev_row[O] + rowSums((Wo %*% Gi) * Wo)
      }
      rm(Zf)
      R <- R^2
      df <- pmax(n_obs - 1L, 0L) - k
      rss <- colSums(R)
      s2 <- numeric(p)
      pos <- df > 0L
      s2[pos] <- rss[pos] / df[pos]
      use <- pos & rss > 0
      pr <- .fit_fdist(s2[use], df[use])
      s <- sqrt(.squeeze_var(s2, df, pr$d0, pr$s20))
      if (!all(is.finite(s)) || any(s <= 0)) {
        stop(
          "tf_biwhiten() could not estimate a finite, positive noise scale ",
          "for every column of 'obj'. Use 'transform = \"identity\"' instead",
          call. = FALSE
        )
      }
      d0 <- pr$d0
      s0 <- sqrt(pr$s20)
      t <- rep(1, n)
      if (rows) {
        q <- drop(R %*% (1 / s^2))
        df_i <- rowSums(obs) - lev_row
        has <- df_i > 0 & q > 0
        t[has] <- sqrt(q[has] / df_i[has])
        if (!all(is.finite(t)) || any(t <= 0)) {
          stop(
            "internal error: the row scale is not finite and positive in every row"
          )
        }
        t[has] <- t[has] / exp(mean(log(t[has])))
      }
      rm(
        list = setdiff(ls(all.names = TRUE), c("s", "t", "n", "p", "d0", "s0"))
      )
      list(
        name = nm,
        scale = s,
        rowscale = t,
        d0 = d0,
        s0 = s0,
        forward = function(y, cols) {
          if (length(cols) != ncol(y)) {
            stop(
              "whiten forward: 'cols' must index every column of the block",
              call. = FALSE
            )
          }
          if (NROW(y) != n) {
            stop(
              "whiten forward: the block must carry every row of the fitted matrix",
              call. = FALSE
            )
          }
          z <- base$forward(y) / rep(s[cols], each = NROW(y))
          if (rows) {
            z <- z / t
          }
          z
        },
        backward = function(z, col, row = NULL) {
          if (length(col) != NROW(z)) {
            stop(
              "whiten backward: one column index per cell is required",
              call. = FALSE
            )
          }
          if (any(col < 1L | col > p)) {
            stop(
              "whiten backward: column index outside the fitted matrix",
              call. = FALSE
            )
          }
          if (rows) {
            if (is.null(row)) {
              stop(
                "whiten backward: this transform has a row scale; 'row' is required",
                call. = FALSE
              )
            }
            if (length(row) != NROW(z)) {
              stop(
                "whiten backward: one row index per cell is required",
                call. = FALSE
              )
            }
            if (any(row < 1L | row > n)) {
              stop(
                "whiten backward: row index outside the fitted matrix",
                call. = FALSE
              )
            }
          }
          base$backward(if (rows) z * (s[col] * t[row]) else z * s[col])
        }
      )
    }
  )
}

.forward_at <- function(tf, blk, cols) {
  if (identical(tf$forward, identity)) {
    blk
  } else {
    tf$forward(blk, cols)
  }
}

.backward_at <- function(tf, d, map, rows = NULL) {
  if (is.null(rows)) {
    tf$backward(d, map[, "col"], map[, "row"])
  } else {
    tf$backward(d, map[rows, "col"], map[rows, "row"])
  }
}

.tf_unfitted <- function(name) {
  force(name)
  list(name = name, backward = function(z, ...) z)
}
