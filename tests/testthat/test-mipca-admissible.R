# Regression: a fully missing column ran silently with every draw of that
# column NaN.

test_that("a fully missing column is refused, not silently returned as NaN", {
  set.seed(42)
  n <- 30L
  p <- 12L
  X <- matrix(stats::rnorm(n * p), n, p)
  X[sample.int(n * p, 40L)] <- NA_real_

  run <- function(Z) {
    mipca_core(
      Z,
      rank = 2L,
      iter_sampling = 3L,
      iter_warmup = 2L,
      thin = 1L,
      chains = 1L
    )
  }

  # Non-vacuity: the fixture must fit with finite draws before the refusal.
  fit <- run(X)
  ## [[1L]]: on the per-chain list these two checks would pass vacuously.
  d <- fit$draws[[1L]]
  expect_gt(length(d), 0L)
  expect_false(anyNA(d))

  X[, 5L] <- NA_real_
  expect_error(run(X), "no observed cells")
})
