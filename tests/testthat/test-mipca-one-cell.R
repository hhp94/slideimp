# Regression: a one-cell diagnostic sample was refused as an internal error,
# because a one-row map read comes back NAMED and .diag_sample()'s identical()
# guard failed on the name. Common for a flank target with one or two holes.

test_that("a fit whose diagnostic sample is one cell finalizes", {
  set.seed(20261001L)
  n <- 30L
  p <- 12L
  X <- matrix(stats::rnorm(n * p), n, p)
  dimnames(X) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  X[4L, 7L] <- NA_real_

  fit <- mipca(
    X,
    rank = 2L,
    transform = "identity",
    chains = 2L,
    iter_sampling = 8L,
    iter_warmup = 4L,
    seed = 1L
  )

  # Non-vacuity: the fit must really be the one-cell shape.
  expect_identical(nrow(fit$chunks[[1L]]$map), 1L)
  cells <- summary(fit)$cells
  expect_identical(nrow(cells), 1L)
  expect_identical(cells$name, "c07")
  expect_identical(cells$row, 4L)
  expect_false(anyNA(as.matrix(fit, 1L)))
})
