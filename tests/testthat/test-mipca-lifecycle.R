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
  fit <- lifecycle_fit()
  expect_identical(fit$lifecycle, "ok")

  fit$draws[[2L]] <- fit$draws[[2L]][-1L, , drop = FALSE]
  expect_error(update(fit, 4L))
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
  expect_error(clamp_frac(fit), "marked failed")

  expect_output(print(fit), "FAILED")
})

test_that("a fit does not serialize the frame it was built in", {
  ballast_mb <- 8
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  local({
    ballast <- matrix(runif(1e6L), 1000L)
    expect_gt(as.numeric(object.size(ballast)) / 1e6, ballast_mb - 0.1)
    saveRDS(lifecycle_fit(), path)
  })

  expect_lt(file.size(path) / 1e6, 2)
})
