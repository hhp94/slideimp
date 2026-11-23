skip_on_cran()

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

  fit <- run(X)
  d <- fit$draws[[1L]]
  expect_gt(length(d), 0L)
  expect_false(anyNA(d))

  X[, 5L] <- NA_real_
  expect_error(
    run(X),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )
})

test_that("abb_cols() counts each tail over the observed cells, strictly past the bound", {
  b <- 0.02
  mid <- 0.5
  cols <- list(
    lo_19 = c(rep(b / 2, 19), mid), # 19 of 20 below: listed
    lo_18 = c(rep(b / 2, 18), mid, mid), # 18 of 20: not listed
    at_lo = c(rep(b, 19), b / 2), # a cell at the bound is not below it
    hi_19 = c(rep(1 - b / 2, 19), mid), # listed
    at_hi = c(rep(1 - b, 19), 1), # a cell at 1 - bound is not above it
    na_10 = c(rep(NA, 10), rep(b / 2, 10)), # 10 of 10 observed: listed
    nan_10 = c(rep(NaN, 10), rep(b / 2, 9), mid), # 9 of 10: not listed
    empty = rep(NA_real_, 20), # never listed
    split = rep(c(0, 1), each = 10), # half in each tail: not listed
    zeros = rep(0, 20) # listed
  )
  X <- do.call(cbind, cols)
  rownames(X) <- sprintf("r%02d", 1:20)
  expect_identical(
    abb_cols(X, bound = b),
    match(
      c("lo_19", "hi_19", "na_10", "zeros"),
      colnames(X)
    )
  )
  # q = 0.9 needs 18 of 20 and 9 of 10, so the columns with 18 and 9 join
  expect_identical(
    abb_cols(X, bound = b, q = 0.9),
    match(
      c("lo_19", "lo_18", "hi_19", "na_10", "nan_10", "zeros"),
      colnames(X)
    )
  )
  # a bound past every cell at the old bound counts those cells
  expect_true(all(
    match(c("at_lo", "at_hi"), colnames(X)) %in% abb_cols(X, bound = 2 * b)
  ))
})

test_that("abb_cols() and tf_sv() refuse a cell outside [0, 1], and skip NA and NaN", {
  set.seed(5)
  X <- matrix(stats::runif(20 * 4), 20, 4)
  dimnames(X) <- list(sprintf("r%02d", 1:20), sprintf("c%d", 1:4))
  X[1, 1] <- NA_real_
  X[2, 2] <- NaN
  X[3, 3] <- 0
  X[4, 3] <- 1
  fwd <- tf_sv(nrow(X))$forward
  expect_no_error(abb_cols(X))
  expect_no_error(fwd(X))
  # a matrix with no observed cell has no cell to check
  empty <- X
  empty[] <- NA_real_
  expect_identical(abb_cols(empty), integer(0))
  expect_no_error(fwd(empty))

  for (bad in c(-1e-12, 1 + 1e-12, -Inf, Inf)) {
    Y <- X
    Y[5, 4] <- bad
    if (is.finite(bad)) {
      r <- range(Y, na.rm = TRUE)
      expect_error(
        abb_cols(Y),
        sprintf(
          "every observed cell must lie in [0, 1] (this one spans %s to %s)",
          format(r[1L]),
          format(r[2L])
        ),
        fixed = TRUE,
        info = format(bad)
      )
    }
    expect_error(
      fwd(Y),
      "every observed cell of 'obj' must lie in [0, 1]",
      fixed = TRUE,
      info = format(bad)
    )
  }
})
