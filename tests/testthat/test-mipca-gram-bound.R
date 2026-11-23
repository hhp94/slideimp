skip_on_cran()

test_that("rank min(n, p) - 2 fits and min(n, p) - 1 refuses", {
  set.seed(20261002L)
  n <- 10L
  p <- 30L
  X <- matrix(stats::runif(n * p), n, p)
  dimnames(X) <- list(
    sprintf("r%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  holes <- cbind(1:6, c(2L, 5L, 9L, 14L, 20L, 27L))
  X[holes] <- NA_real_
  top <- min(n, p) - 2L

  expect_gt(2L * (top + 1L), min(n, p))

  fit <- mipca(
    X,
    rank = top,
    transform = "identity",
    chains = 2L,
    iter_sampling = 10L,
    iter_warmup = 10L,
    seed = 1L
  )
  expect_identical(fit$chunks[[1L]]$prov$rank, top)
  m <- as.matrix(fit, 1L)
  expect_true(all(is.finite(m[holes])))

  expect_error(
    mipca(X, rank = top + 1L, transform = "identity", seed = 1L),
    "'rank' must be at most the smaller dimension minus two (8 here",
    fixed = TRUE
  )
})
