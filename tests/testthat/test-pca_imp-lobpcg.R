skip_on_cran()

# lobpcg ----
test_that("a tampered lobpcg_control object is re-validated, not taken as-is", {
  set.seed(1234)
  x <- sim_mat(30, 40, perc_total_na = 0.2, perc_col_na = 1)$input

  run <- function(ctrl) {
    pca_imp(
      x,
      ncp = 2,
      solver = "lobpcg",
      lobpcg_control = ctrl,
      miniter = 2,
      maxiter = 5,
      seed = 1
    )
  }

  # the fixture is a VALID control object, so anything below that errors does
  # so because of the field that was assigned into, not because of the object.
  ok <- lobpcg_control()
  expect_s3_class(ok, "slideimp_lobpcg_control")
  expect_no_error(suppressWarnings(run(ok)))

  # Regression: new_lobpcg_control() took any object carrying the class with
  # `out <- x` and no field checks. A control object is an ordinary list, so
  # each of these survived construction. Before the fix: `tol = -1` ran to
  # completion silently, `maxiter = NULL` raised `missing value where
  # TRUE/FALSE needed` from the `out$maxiter == 0L` test, and `maxiter = -5L`
  # reached C++ and came back as an int-range error naming no argument.
  bad_tol <- ok
  bad_tol$tol <- -1
  expect_error(run(bad_tol), "lobpcg_control\\$tol")

  bad_maxiter <- ok
  bad_maxiter$maxiter <- -5L
  expect_error(run(bad_maxiter), "lobpcg_control\\$maxiter")

  bad_warmup <- ok
  bad_warmup$warmup_iters <- -1L
  expect_error(run(bad_warmup), "lobpcg_control\\$warmup_iters")

  # above INT_MAX the value does not fit the `int` the C++ backend takes; R
  # refuses it first
  big <- ok
  big$maxiter <- 2^40
  expect_error(run(big), "lobpcg_control\\$maxiter")

  # a dropped field is an error, not a silent fall back to the default
  dropped <- ok
  dropped$maxiter <- NULL
  expect_error(run(dropped), "maxiter")

  # a field replaced by something that is not a number at all
  wrong_type <- ok
  wrong_type$tol <- "small"
  expect_error(run(wrong_type), "lobpcg_control")

  # and the constructor itself still refuses the same values
  expect_error(lobpcg_control(tol = -1), "tol")
  expect_error(lobpcg_control(maxiter = -5L), "maxiter")
  expect_error(lobpcg_control(warmup_iters = 2^40), "warmup_iters")
})

test_that("LOBPCG mode during warmup matches forced exact path", {
  set.seed(1234)
  x <- sim_mat(20, 80)$input
  pca_iters <- 12L

  # `run_pca_fixed_iters()` drives pca_imp_internal_cpp() directly, and the
  # kernel is silent on non-convergence by design - pca_imp() raises that
  # warning from R. See "pca_imp warns from R when the EM loop hits maxiter".
  ref <- run_pca_fixed_iters(
    x,
    solver = "exact",
    pca_iters = pca_iters
  )

  got <- run_pca_fixed_iters(
    x,
    solver = "lobpcg",
    ctrl = lobpcg_control(
      maxiter = 50L,
      warmup_iters = pca_iters + 1L,
      tol = 1e-10
    ),
    pca_iters = pca_iters
  )

  expect_false(anyNA(ref$mat))
  expect_false(anyNA(got$mat))
  expect_equal(dim(got$mat), dim(ref$mat))

  # same numerical path: warmup_iters > pca_iters means LOBPCG never triggers.
  # Measured bitwise equal (worst abs and rel error both 0); the absolute bound
  # is the one this block had, and the relative bound matches it.
  expect_worst_err(got$mat, ref$mat, abs_tol = 1e-10, rel_tol = 1e-10)

  expect_equal(ref$n_lobpcg_ok, 0L)
  expect_equal(ref$n_lobpcg_bad, 0L)
  expect_equal(ref$n_exact, pca_iters)

  expect_equal(got$n_lobpcg_ok, 0L)
  expect_equal(got$n_lobpcg_bad, 0L)
  expect_equal(got$n_exact, pca_iters)
})

test_that("LOBPCG enabled agrees with exact eigensolver branch", {
  cases <- list(
    tall = c(80L, 25L),
    wide = c(25L, 80L)
  )

  set.seed(1234)

  for (case in cases) {
    n <- case[[1L]]
    p <- case[[2L]]

    x <- sim_mat(n, p)$input
    pca_iters <- 14L
    warmup <- 2L

    # kernel-direct, so no non-convergence warning. See the note above.
    ref <- run_pca_fixed_iters(
      x,
      solver = "exact",
      pca_iters = pca_iters
    )

    got <- run_pca_fixed_iters(
      x,
      solver = "lobpcg",
      ctrl = lobpcg_control(
        maxiter = 100L,
        warmup_iters = warmup,
        tol = 1e-8
      ),
      pca_iters = pca_iters
    )

    expect_false(anyNA(ref$mat))
    expect_false(anyNA(got$mat))
    expect_equal(dim(got$mat), dim(ref$mat))
    # worst cell over both shapes, measured: 1.64e-8 absolute, 1.67e-7
    # relative (over the nonzero cells), where one component fewer misses by
    # 0.107 or more; the bounds leave about 60 times the measured error
    expect_worst_err(got$mat, ref$mat, abs_tol = 1e-6, rel_tol = 1e-5)

    # forced exact: every iteration uses dsyevr.
    expect_equal(ref$n_lobpcg_ok, 0L)
    expect_equal(ref$n_lobpcg_bad, 0L)
    expect_equal(ref$n_exact, pca_iters)

    # LOBPCG: warmup iterations use exact, then LOBPCG should converge.
    # the structural invariant always holds
    expect_equal(got$n_exact + got$n_lobpcg_ok + got$n_lobpcg_bad, pca_iters)
    # warmup did use exact
    expect_gte(got$n_exact, warmup)
    # the vast majority of post-warmup iters used LOBPCG successfully
    expect_gte(got$n_lobpcg_ok, pca_iters - warmup - 2L)
  }
})

test_that("LOBPCG tolerance stays relative when the Gram norm is below 1", {
  # `A` is the row-weight-scaled Gram and `row.w` sums to 1, so A's diagonal is
  # a weighted column variance. Under `scale = TRUE` that diagonal is exactly 1
  # and ||A||_inf >= 1; under `scale = FALSE` on bounded data it lands well
  # below 1 (0.16 for the matrix below). Regression: lobpcg_solve() floored
  # normA at 1.0, so for those Grams the convergence test became ABSOLUTE and
  # the error grew in lock step with the data scale - 6.2e-09 at factor 1 but
  # 3.0e-04 at factor 0.01, against a requested tol of 1e-9. The invariant is
  # that the RELATIVE error does not move when the data is rescaled.
  set.seed(11)
  n <- 60L
  p <- 20L
  latent <- matrix(rnorm(n * 3L), n, 3L)
  loadings <- matrix(rnorm(3L * p), 3L, p)
  base <- plogis(scale(latent %*% loadings) * 0.6)
  na_idx <- sample(seq_len(n * p), round(0.2 * n * p))

  ctrl <- lobpcg_control(warmup_iters = 3L, tol = 1e-9, maxiter = 40L)
  # miniter = maxiter fixes the iteration count, so both solvers take the
  # same number of steps and any difference is the eigensolvers'. The
  # threshold, far above the criterion at that count, decides only
  # `converged`.
  run <- function(x, solver) {
    pca_imp(
      x,
      ncp = 3,
      scale = FALSE,
      solver = solver,
      seed = 1,
      nb.init = 1,
      maxiter = 40,
      miniter = 40,
      threshold = 1e-8,
      post_imp = FALSE,
      lobpcg_control = ctrl
    )
  }

  factors <- c(1, 0.01)
  max_abs <- numeric(length(factors))
  max_rel <- numeric(length(factors))

  for (i in seq_along(factors)) {
    x <- base * factors[[i]]
    x[na_idx] <- NA_real_

    ref <- run(x, "exact")
    got <- run(x, "lobpcg")
    # LOBPCG ran, or the comparison is the exact fit against itself
    expect_gt(attr(got, "n_lobpcg_ok"), 0L)

    a <- ref[is.na(x)]
    b <- got[is.na(x)]
    err <- abs(a - b)

    # the relative figure is taken over the imputed cells whose exact value is
    # not effectively zero; near zero it says nothing.
    keep <- abs(a) > 1e-12
    expect_gt(sum(keep), 0L)

    max_abs[[i]] <- max(err)
    max_rel[[i]] <- max(err[keep] / abs(a[keep]))
  }

  # worst case, both figures. The absolute bound tracks the data scale; the
  # relative bound does not, which is the whole point.
  expect_lt(max_abs[[1L]], 1e-8 * factors[[1L]])
  expect_lt(max_abs[[2L]], 1e-8 * factors[[2L]])
  expect_lt(max_rel[[1L]], 1e-8)
  expect_lt(max_rel[[2L]], 1e-8)

  # scale invariance: rescaling the data by 100x must not move the relative
  # error. Before the fix this ratio was ~5e4.
  expect_lt(max(max_rel) / min(max_rel), 10)
})
