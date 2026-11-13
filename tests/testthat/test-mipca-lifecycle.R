# Contracts of the environment-backed fit: a torn fit refuses, and a fit does
# not carry the frame it was born in. si_mi_pca is environment-backed; si_mi
# is a value.

## The tear must happen INSIDE update()'s mutation region, so the test breaks
## the stored draws (nothing there can be mocked to throw).
lifecycle_fit <- function() {
  set.seed(20260819L)
  n <- 30L
  p <- 60L
  Z <- matrix(rnorm(n * 3L), n) %*%
    matrix(rnorm(3L * p), 3L) +
    0.4 * matrix(rnorm(n * p), n)
  X <- plogis(Z / 2)
  dimnames(X) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  X[cbind(sample.int(n, 60L, TRUE), sample.int(40L, 60L, TRUE))] <- NA_real_
  mipca(
    X,
    rank = 3L,
    transform = "identity",
    keep_state = TRUE,
    chains = 2L,
    iter_sampling = 8L,
    iter_warmup = 10L,
    seed = 1L
  )
}

test_that("a fit torn by a failed update() refuses every method but print()", {
  # Catches an assignment moved outside the failed-marked region, or a method
  # added without .life_check().
  fit <- lifecycle_fit()
  expect_identical(fit$lifecycle, "ok")

  # Chain 2 one row short: update() grows chain 1, then throws on chain 2.
  fit$draws[[2L]] <- fit$draws[[2L]][-1L, , drop = FALSE]
  expect_error(update(fit, 4L))
  # Non-vacuity: chain 1 grew and chain 2 did not, so the fit really is torn.
  expect_identical(
    c(ncol(fit$draws[[1L]]), ncol(fit$draws[[2L]])),
    c(12L, 8L)
  )

  expect_identical(fit$lifecycle, "failed")
  expect_error(update(fit, 4L), "marked failed")
  expect_error(mipca_finalize(fit), "marked failed")
  expect_error(as.matrix(fit), "marked failed")
  expect_error(summary(fit), "marked failed")
  expect_error(inspect_completion(fit), "marked failed")
  # Regression: clamp_frac() read the torn draws and returned a number.
  expect_error(clamp_frac(fit), "marked failed")

  # print() is exempt: a broken object must be visible AS broken.
  expect_output(print(fit), "FAILED")
})

test_that("a fit does not serialize the frame it was built in", {
  # Pins parent = emptyenv(): a fit must not carry its caller's matrix to disk.
  # runif, not zeros: saveRDS gzips, and leaked zeros would pass the bar.
  ballast_mb <- 8
  # Base tempfile: withr is not a dependency.
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  local({
    ballast <- matrix(runif(1e6L), 1000L) # 8 MB, in the fit's birth frame
    expect_gt(as.numeric(object.size(ballast)) / 1e6, ballast_mb - 0.1)
    saveRDS(lifecycle_fit(), path)
  })

  # Falsified: clean ~0.3 MB, one field closing over the frame ~5.6 MB.
  expect_lt(file.size(path) / 1e6, 2)
})
