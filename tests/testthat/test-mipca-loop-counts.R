# Regression: iter_sampling = 0 with keep_state = TRUE made summary() segfault
# R. Pinned on the message.

test_that("iter_sampling = 0 is refused rather than segfaulting summary()", {
  set.seed(1)
  Y <- matrix(stats::rnorm(320), 16, 20)
  Y[sample(seq_along(Y), 30L)] <- NA
  dimnames(Y) <- list(paste0("r", 1:16), paste0("c", 1:20))

  # Non-vacuity: the fixture must fit at all.
  ok <- mipca(
    Y,
    rank = 2L,
    transform = "identity",
    chains = 1L,
    iter_warmup = 4L,
    iter_sampling = 1L,
    keep_state = TRUE
  )
  expect_s3_class(ok, "si_mi_pca")

  expect_error(
    mipca(
      Y,
      rank = 2L,
      transform = "identity",
      chains = 1L,
      iter_warmup = 4L,
      iter_sampling = 0L,
      keep_state = TRUE
    ),
    "iter_sampling must be a single integer >= 1",
    fixed = TRUE
  )
})
