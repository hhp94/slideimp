test_that("ess_ok is NA, not NaN, when no ESS can be estimated", {
  set.seed(2)
  n <- 20L
  p <- 40L
  Y <- matrix(stats::rnorm(n * 3L), n) %*%
    matrix(stats::rnorm(3L * p), 3L) +
    0.3 * matrix(stats::rnorm(n * p), n)
  dimnames(Y) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  Y[cbind(sample.int(n, 40L, TRUE), sample.int(p, 40L, TRUE))] <- NA_real_
  nl <- slideimp::sample_na_loc(Y, num_na = 12L, n_rows = 3L, n_reps = 1L)

  tw <- mipca_tune(
    Y,
    2L,
    na_loc = nl,
    transform = "identity",
    n_diag = 20L,
    chains = 2L,
    iter_warmup = 5L,
    iter_sampling = 10L,
    seed = 1L
  )

  expect_true(is.na(tw$error[1L]))
  expect_true(is.na(tw$ess_min[1L]))

  expect_true(is.na(tw$ess_ok[1L]))
  expect_false(is.nan(tw$ess_ok[1L]))
})
