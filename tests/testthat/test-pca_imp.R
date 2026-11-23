skip_on_cran()

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- sim_mat(
    20,
    50,
    perc_total_na = 0.25,
    perc_col_na = 1,
    rho = 0.75
  )$input
  # a row with no observed value is imputed with the rest
  to_test[1, ] <- NA
  res <- pca_imp(to_test, ncp = 2, seed = 1234)
  expect_false(anyNA(res))
  # a column with no observed value is refused before the fit
  to_test[, 1] <- NA
  expect_error(
    pca_imp(to_test, ncp = 2, seed = 1234),
    "'obj' holds infinite values or columns with no observed value",
    fixed = TRUE
  )
})

test_that("row.w = 'n_miss' floors near-zero weights", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("row", 1:10)
  colnames(mat) <- paste0("col", 1:10)
  mat[1, -1] <- NA # row 1 has 9/10 missing -> weight = 0.1
  mat[2, ] <- NA # row 2 has 10/10 missing -> weight = 0, under the floor
  mat[3, 1:5] <- NA # row 3 has 5/10 missing -> weight = 0.5

  # the weights "n_miss" stands for, floored at 1e-10 as pca_imp() floors
  # them; the premise is that one of them is below the floor
  w <- unname(1 - rowSums(is.na(mat)) / ncol(mat))
  expect_identical(w[2], 0)
  w[w < 1e-10] <- 1e-10

  # the shortcut is those weights, floor included: the same fit. miniter =
  # maxiter fixes the iteration count, so both fits take the same steps, and
  # the threshold, far above the criterion at that count, decides only
  # `converged`.
  run <- function(row.w) {
    pca_imp(
      mat,
      ncp = 2,
      row.w = row.w,
      seed = 123,
      threshold = 1e-6,
      maxiter = 30,
      miniter = 30
    )
  }
  shortcut <- run("n_miss")
  by_hand <- run(w)
  expect_worst_err(shortcut[,], by_hand[,])
  # the record is the same but for `criterion_final`, a ratio of two
  # near-equal objectives that rounding moves
  same <- setdiff(names(attributes(shortcut)), "criterion_final")
  expect_identical(attributes(shortcut)[same], attributes(by_hand)[same])
  message(sprintf(
    "row.w = 'n_miss' against its weights bitwise: %s",
    identical(shortcut[,], by_hand[,])
  ))
})

test_that("row.w rejects an unknown string and a vector of the wrong length", {
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("row", 1:10)
  colnames(mat) <- paste0("col", 1:10)
  mat[1, 1] <- NA

  expect_error(pca_imp(mat, ncp = 2, row.w = "invalid"), regexp = "row.w")
  expect_error(pca_imp(mat, ncp = 2, row.w = c(67, 69)), regexp = "row.w")
})

# eligibility resolution ----
test_that("pca_imp handles ineligible columns (high miss rate / zero variance) correctly", {
  set.seed(1234)
  to_test <- sim_mat(40, 12, perc_total_na = 0.25, perc_col_na = 0.6)$input

  # force ineligible columns:
  # - Column 1: miss_rate > colmax (0.925 > 0.9)
  to_test[1:37, 1] <- NA
  mean_1 <- mean(to_test[, 1], na.rm = TRUE)
  # - Column 2: constant column -> variance = 0
  to_test[, 2] <- 69
  # - Column 3: near-zero variance with a few NAs
  to_test[, 3] <- 3.14 + rnorm(40, sd = 1e-10)
  to_test[1:3, 3] <- NA
  mean_3 <- mean(to_test[, 3], na.rm = TRUE)
  expect_true(anyNA(to_test))
  expect_true(col_vars(to_test[, 3, drop = F]) < .Machine$double.eps)

  # 1. post_imp = TRUE: ineligible columns are mean-imputed
  res <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 5,
    seed = 1234,
    colmax = 0.9,
    scale = FALSE
  )

  expect_false(anyNA(res))

  # ineligible high-miss column becomes constant (mean imputation)
  expect_true(all(res[1:37, 1] == mean_1))
  expect_equal(unname(res[1:37, 1]), rep(mean_1, times = 37))

  # constant column untouched
  expect_equal(length(unique(res[, 2])), 1L)

  # near-zero variance column: NAs filled with column mean
  expect_equal(unname(res[1:3, 3]), rep(mean_3, 3))

  # 2. post_imp = FALSE: only eligible columns are PCA-imputed;
  # ineligible columns keep their original NAs
  res_no_post <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 5,
    seed = 1234,
    colmax = 0.9,
    post_imp = FALSE,
    scale = FALSE
  )

  expect_true(anyNA(res_no_post))
  # the ineligible columns keep their ORIGINAL NAs, in the original positions,
  # and nothing else is touched. `expect_gt(mean(is.na(res_no_post[, 1])),
  # mean_1)` used to stand in for this, comparing an NA PROPORTION (0.95)
  # against the column's DATA mean (0.448) - two quantities with nothing in
  # common. It passed on the coincidence that sim_mat() draws near zero:
  # measured on the same fixture shifted by +10, a CORRECT result fails it.
  # rows 1:37 were forced NA, and sim_mat() left one more NA further down, so
  # the set is taken from the input rather than written out.
  na_1 <- unname(which(is.na(to_test[, 1])))
  expect_gt(length(na_1), 37L)
  expect_identical(unname(which(is.na(res_no_post[, 1]))), na_1)
  expect_identical(unname(res_no_post[-na_1, 1]), unname(to_test[-na_1, 1]))
  expect_equal(unique(res_no_post[, 2]), 69)
  expect_identical(unname(which(is.na(res_no_post[, 3]))), 1:3)
  expect_identical(unname(res_no_post[4:40, 3]), unname(to_test[4:40, 3]))
  expect_false(anyNA(res_no_post[, 4:12]))
})

test_that("pca_imp refuses an ncp above the components its eligible columns give", {
  set.seed(1234)
  to_test <- sim_mat(30, 8, perc_total_na = 0.1, perc_col_na = 0.3)$input
  # this column will exceed colmax
  to_test[1:29, 1] <- NA
  expect_no_error(
    res <- pca_imp(
      to_test,
      ncp = 3,
      nb.init = 3,
      seed = 1234,
      colmax = 0.9,
      post_imp = FALSE
    )
  )
  # make most columns ineligible (all-NA)
  to_test[, 1:6] <- NA
  for (i in 1:6) {
    to_test[sample.int(30, size = 1), i] <- rnorm(1)
  }

  # only 2 eligible columns left -> ncp = 3 > min(28, 1) -> error (1 usable component)
  expect_error(
    pca_imp(
      to_test,
      ncp = 3,
      nb.init = 3,
      seed = 1234,
      colmax = 0.9,
      post_imp = TRUE
    ),
    "must be at most the number of usable components",
    fixed = TRUE
  )
})

test_that("pca_imp clamps imputed values to specified bounds", {
  set.seed(1234)
  to_test <- sim_mat(30, 8, perc_total_na = 0.1)$input
  na_pos <- which(is.na(to_test), arr.ind = TRUE)
  # the premise: there are imputed cells to clamp
  expect_gt(nrow(na_pos), 0L)
  observed_mask <- !is.na(to_test)
  observed_before <- to_test[observed_mask]

  # lower bound far above any plausible imputation -> all imputed == 999
  res_lo <- pca_imp(
    to_test,
    ncp = 2,
    seed = 1234,
    clamp = c(999, Inf),
    post_imp = FALSE
  )
  expect_true(all(res_lo[na_pos] == 999))
  expect_equal(res_lo[observed_mask], observed_before)

  # upper bound far below any plausible imputation -> all imputed == -999
  res_hi <- pca_imp(
    to_test,
    ncp = 2,
    seed = 1234,
    clamp = c(-Inf, -999),
    post_imp = FALSE
  )
  expect_true(all(res_hi[na_pos] == -999))
  expect_equal(res_hi[observed_mask], observed_before)
})

test_that("pca_imp doesn't mess up the original object", {
  set.seed(1234)
  to_test <- sim_mat(30, 30, perc_total_na = 0.1, perc_col_na = 1)$input
  expect_true(anyNA(to_test))
  # a copy with memory of its own: `passed_obj <- to_test` would share the
  # vector, and the kernel reads `obj` in place, so a write into it from C++
  # would change both and still compare equal
  passed_obj <- to_test + 0
  res <- pca_imp(
    to_test,
    ncp = 3,
    nb.init = 3,
    seed = 1234,
    colmax = 0.9,
    post_imp = FALSE
  )
  expect_identical(to_test, passed_obj)
  expect_equal(is.na(to_test), is.na(passed_obj))
})

test_that("pca_imp restores object even on bad input", {
  set.seed(1234)
  to_test <- sim_mat(30, 30, perc_total_na = 0.1, perc_col_na = 1)$input
  # a copy with memory of its own, as in the test above
  passed_obj <- to_test + 0
  expect_true(anyNA(to_test))
  # refused in R, before the kernel
  expect_error(
    pca_imp(
      to_test,
      ncp = 9999,
      nb.init = 3,
      seed = 1234,
      colmax = 0.9,
      post_imp = FALSE
    ),
    "Assertion on 'ncp' failed",
    fixed = TRUE
  )
  expect_identical(to_test, passed_obj)

  # refused inside the kernel, which has read the input by then: values
  # whose squares overflow stop it at the Gram
  big <- to_test * 1e200
  passed_big <- big + 0
  expect_error(
    pca_imp(big, ncp = 3, scale = FALSE, seed = 1234, post_imp = FALSE),
    "the cross-product of the data is not finite",
    fixed = TRUE,
    class = "Rcpp::exception"
  )
  expect_identical(big, passed_big)
})

test_that("Throw on Inf", {
  set.seed(1234)
  to_test <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  to_test[1, 1] <- Inf
  expect_error(
    pca_imp(to_test, ncp = 3),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )
  to_test[1, 1] <- -Inf
  expect_error(
    pca_imp(to_test, ncp = 3),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )
})

test_that("Inf is refused ahead of the paths that used to hide it", {
  set.seed(1234)

  # no missing cells at all: the "No missing values" early return used to hand
  # the matrix back with the Inf still in it.
  complete <- sim_mat(20, 20, perc_total_na = 0)$input
  complete[3, 4] <- Inf
  expect_error(
    pca_imp(complete, ncp = 2),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )

  # otherwise infeasible: the infeasibility abort used to win the race and
  # report eligible columns without ever mentioning the Inf.
  infeasible <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  infeasible[3, 4] <- Inf
  infeasible[, 6:20] <- 1
  expect_error(
    pca_imp(infeasible, ncp = 8),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )

  # and the input is scanned once per call, not once per `nb.init`: the
  # kernel has no scan of its own, so the scan at the door is the gate,
  # ahead of the restarts. Counted on an input that reaches every restart.
  # The one other scan of the call reads the result, to check the write.
  holed <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  expect_true(anyNA(holed))
  real_scan <- col_scan
  n_scan <- 0L
  local_mocked_bindings(
    col_scan = function(obj, ...) {
      n_scan <<- n_scan + identical(obj, holed)
      real_scan(obj, ...)
    }
  )
  res <- pca_imp(holed, ncp = 2, nb.init = 3)
  expect_false(anyNA(res))
  expect_identical(n_scan, 1L)
})

test_that("the no-missing-values return has the shape of the main path", {
  set.seed(1234)
  complete <- sim_mat(20, 10, perc_total_na = 0)$input
  expect_false(anyNA(complete))

  msgs <- testthat::capture_messages(
    r_early <- pca_imp(complete, ncp = 2)
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "No missing values in input")
  expect_identical(r_early[,, drop = FALSE], complete)

  incomplete <- complete
  incomplete[2, 3] <- NA
  r_main <- pca_imp(incomplete, ncp = 2)

  # the same class, and the same attributes with the same types, so a caller
  # cannot tell from the result which path ran
  expect_identical(class(r_early), class(r_main))
  expect_identical(
    vapply(attributes(r_early), typeof, character(1)),
    vapply(attributes(r_main), typeof, character(1))
  )
  expect_identical(attr(r_early, "imp_method"), "pca")
  expect_identical(attr(r_early, "fallback"), FALSE)
  expect_identical(attr(r_early, "solver_requested"), "auto")

  # and the solver attributes say that no solve ran
  expect_identical(attr(r_early, "solver_chosen"), NA_character_)
  expect_false(attr(r_early, "solver_probed"))
  expect_identical(attr(r_early, "n_iter"), 0L)
  expect_identical(attr(r_early, "criterion_final"), NA_real_)
  expect_identical(attr(r_early, "converged"), NA)
  expect_identical(attr(r_early, "n_exact"), 0L)
})

test_that("pca_imp reports whether the auto verdict came from a probe", {
  set.seed(1234)
  # n_gram = min(nrow, n_elig) = 60, below the 250 threshold, so the R-side
  # heuristic demotes this to exact before the kernel is ever called and no
  # probe runs. `solver_chosen` still says "exact", because that is what ran.
  small <- sim_mat(60, 60, perc_total_na = 0.2, perc_col_na = 1)$input

  res <- pca_imp(small, ncp = 2, solver = "auto", seed = 1)
  expect_identical(attr(res, "solver_chosen"), "exact")
  expect_false(attr(res, "solver_probed"))

  # a forced solver is not a verdict either: auto never ran
  for (s in c("exact", "lobpcg")) {
    r <- pca_imp(small, ncp = 2, solver = s, seed = 1)
    expect_identical(attr(r, "solver_chosen"), s)
    expect_false(attr(r, "solver_probed"))
  }

  # and the TRUE side, so this is not one-sided. n_gram = 260 clears the
  # threshold, and threshold = 0 forces 60 EM iterations, comfortably more than
  # the auto_min_exact_iter warmup plus the 5 LOBPCG probe iterations the
  # decision needs, so the probe always completes here. WHICH solver it picks
  # is a wall-clock comparison and is deliberately not asserted.
  big <- sim_mat(260, 400, perc_total_na = 0.1, perc_col_na = 1)$input
  probed <- suppressWarnings(pca_imp(
    big,
    ncp = 2,
    solver = "auto",
    threshold = 0,
    maxiter = 60,
    miniter = 60,
    seed = 1
  ))
  expect_true(attr(probed, "solver_probed"))
  expect_true(attr(probed, "solver_chosen") %in% c("exact", "lobpcg"))
})

test_that("pca_imp warns from R when the EM loop hits maxiter", {
  set.seed(1234)
  x <- sim_mat(30, 40, perc_total_na = 0.2, perc_col_na = 1)$input

  # threshold = 0 is unreachable (the criterion and the objective are both
  # non-negative), so every restart runs out at `maxiter`.
  run <- function(nb.init = 1L) {
    pca_imp(
      x,
      ncp = 2,
      solver = "exact",
      threshold = 0,
      maxiter = 3,
      miniter = 3,
      nb.init = nb.init,
      seed = 1,
      post_imp = FALSE
    )
  }

  w <- expect_warning(run(), "Stopped after")

  # Regression: this warning used to be raised by Rcpp::warning() inside the EM
  # loop, which for a std::string is a bare Rf_warning. Under options(warn = 2),
  # or any calling handler that invokes a restart, R longjmps out of the C++
  # frame; END_RCPP catches exceptions, not longjmps, so the destructors of the
  # kernel's working matrices, the Rcout precision guard and the wrapper's
  # RNGScope were all skipped. A cli/rlang condition is the proof that R raises
  # it now.
  expect_s3_class(w, "rlang_warning")
  expect_false(inherits(w, "simpleWarning"))

  # one warning per call, not one per restart: only the winning restart's
  # `converged` flag reaches the caller, so only it may warn.
  msgs <- character(0)
  withCallingHandlers(
    run(nb.init = 3L),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(grep("Stopped after", msgs), 1L)

  # the warning has to agree with the flag it is derived from
  res <- suppressWarnings(run())
  expect_false(attr(res, "converged"))

  # and a converging call stays silent
  expect_no_warning(
    pca_imp(
      x,
      ncp = 2,
      solver = "exact",
      threshold = 1e-6,
      maxiter = 1000,
      seed = 1,
      post_imp = FALSE
    )
  )
})

# argument validation ----
test_that("pca_imp validates seed and the miniter/maxiter ordering from R", {
  set.seed(1234)
  x <- sim_mat(20, 20, perc_total_na = 0.1, perc_col_na = 1)$input

  # maxiter is deliberately small everywhere below, so the non-convergence
  # warning is expected and not what any of these assertions are about.
  run <- function(...) {
    suppressWarnings(pca_imp(
      x,
      ncp = 2,
      solver = "exact",
      ...
    ))
  }

  # Regression: this was `stop("`seed` too large")` - base stop, no condition
  # class, a hard-coded 2147483647L, and no mention of `nb.init`, which is the
  # argument that makes the bound what it is. `seed * (nb.init - 1)` is what
  # reaches set.seed(), so the bound moves with `nb.init`.
  cap <- .Machine$integer.max %/% 4L
  expect_error(run(seed = cap + 1L, nb.init = 5), "Assertion on 'seed'")

  # the bound is not simply rejecting everything: the largest admissible value
  # is accepted, and the SAME value is rejected only once `nb.init` raises the
  # multiplier past what that value leaves room for.
  expect_no_error(run(seed = cap, nb.init = 5, miniter = 2, maxiter = 5))
  expect_no_error(run(
    seed = .Machine$integer.max,
    nb.init = 1L,
    miniter = 2,
    maxiter = 5
  ))
  expect_error(
    run(seed = .Machine$integer.max, nb.init = 3L),
    "Assertion on 'seed'"
  )

  # Regression: `miniter > maxiter` was unchecked on the R side and fell
  # through to a bare `Rcpp::stop("miniter must be <= maxiter")`, which is not
  # an rlang condition and names neither value.
  err <- tryCatch(run(miniter = 10, maxiter = 5), error = function(e) e)
  expect_s3_class(err, "rlang_error")
  expect_match(conditionMessage(err), "miniter")
  expect_match(conditionMessage(err), "maxiter")
  expect_no_error(run(miniter = 5, maxiter = 5))
})

# return path ----
test_that("post_imp mean-imputes only the columns that can still hold NA", {
  set.seed(1234)
  x <- sim_mat(40, 30, perc_total_na = 0.1, perc_col_na = 1)$input

  # miniter = maxiter fixes the iteration count, so every fit below takes the
  # same steps; 5 iterations do not converge, hence suppressWarnings()
  run <- function(m, post_imp) {
    suppressWarnings(pca_imp(
      m,
      ncp = 2,
      solver = "exact",
      miniter = 5,
      maxiter = 5,
      seed = 1,
      post_imp = post_imp
    ))
  }

  # Regression: pca_imp() ran mean_imp_col(obj) over every column whenever
  # post_imp = TRUE. On an all-eligible input that copied the matrix twice
  # more (the C++ result, then Rcpp's wrap of it) and scanned every column,
  # to change nothing. The mean part now covers the ineligible columns only,
  # so fallback_part() is asked for exactly those, and for none on an
  # all-eligible input.
  real_fallback_part <- fallback_part
  calls <- list()
  local_mocked_bindings(
    fallback_part = function(obj, cols, method = c("mean", "median")) {
      calls <<- c(calls, list(cols))
      real_fallback_part(obj, cols, method)
    }
  )

  # every column of x is eligible, so the kernel fills every NA and the post
  # step has nothing to do. Check that premise before asserting on it.
  res_no_post <- run(x, post_imp = FALSE)
  expect_false(anyNA(res_no_post))
  expect_length(calls, 0L)
  res_post <- run(x, post_imp = TRUE)
  expect_length(calls, 1L)
  expect_identical(calls[[1L]], integer(0))
  # two fits, so the values agree to rounding
  expect_worst_err(as.vector(res_post), as.vector(res_no_post))
  expect_true(attr(res_post, "post_imp"))

  # column 1 above colmax, column 2 zero variance with NA, column 3 zero
  # variance without NA. The first two are the only columns that can still
  # hold NA after the kernel; the third is ineligible but has nothing to fill,
  # so it is not copied into the fallback either.
  y <- x
  y[sample(40, 38), 1] <- NA
  y[, 2] <- 7
  y[1:3, 2] <- NA
  y[, 3] <- 7
  # sim_mat() may already have left NA in column 1, so take the count from y
  n_na_1 <- sum(is.na(y[, 1]))
  expect_gte(n_na_1, 38L)
  calls <- list()
  res_y <- run(y, post_imp = TRUE)
  expect_length(calls, 1L)
  expect_identical(calls[[1L]], c(1L, 2L))
  expect_false(anyNA(y[, 3]))
  expect_false(anyNA(res_y))
  expect_equal(
    unname(res_y[is.na(y[, 1]), 1]),
    rep(mean(y[, 1], na.rm = TRUE), n_na_1)
  )
  expect_identical(unname(res_y[, 2]), rep(7, 40L))

  # and the result is what the unconditional call used to produce: every
  # column with a hole left after the kernel filled with its mean, one mean
  # part over those columns, as mean_imp_col() built it. `no_post` is a
  # second fit, so the values agree to rounding.
  no_post <- matrix(as.vector(run(y, post_imp = FALSE)), 40L)
  holed <- which(colSums(is.na(no_post)) > 0L)
  expect_identical(holed, c(1L, 2L))
  old_way <- fill_parts(
    no_post,
    list(real_fallback_part(no_post, holed, "mean"))
  )
  expect_worst_err(as.vector(res_y), as.vector(old_way))
})

test_that("pca_imp records what post_imp set to a column mean", {
  set.seed(1234)
  x <- sim_mat(40, 30, perc_total_na = 0.1, perc_col_na = 1)$input
  # column 1 above colmax, column 2 zero variance with NA, column 3 zero
  # variance without NA: ineligible, but with nothing to fill
  x[sample(40, 38), 1] <- NA
  x[, 2] <- 7
  x[1:3, 2] <- NA
  x[, 3] <- 7
  no_cells <- cbind(row = integer(0), col = integer(0))

  run <- function(m, post_imp) {
    suppressWarnings(pca_imp(
      m,
      ncp = 2,
      solver = "exact",
      miniter = 2,
      maxiter = 5,
      seed = 1,
      post_imp = post_imp
    ))
  }

  res <- run(x, post_imp = TRUE)
  expect_false(anyNA(res[,]))
  expect_identical(attr(res, "mean_filled_cols"), c(1L, 2L))
  expect_identical(attr(res, "mean_filled_cells"), no_cells)

  # the record is exactly what PCA did not fill: without post_imp, those are
  # the only columns still holding a missing value
  res_off <- run(x, post_imp = FALSE)
  expect_identical(unname(which(colSums(is.na(res_off[,])) > 0L)), c(1L, 2L))
  expect_identical(attr(res_off, "mean_filled_cols"), integer(0))
  expect_identical(attr(res_off, "mean_filled_cells"), no_cells)

  # the no-missing-values exit fills nothing and still carries both attributes
  full <- matrix(rnorm(40), 8)
  r_early <- suppressMessages(pca_imp(full, ncp = 1))
  expect_identical(attr(r_early, "mean_filled_cols"), integer(0))
  expect_identical(attr(r_early, "mean_filled_cells"), no_cells)
})

test_that("pca_imp refuses a fit whose mse is not finite, as a defect", {
  # The kernel's mse sums the residual over the holes too, so a NaN or Inf
  # imputed value makes it non-finite. Finite input never gets there, so the
  # kernel is stubbed: the real fit, then a broken mse or value planted in it.
  set.seed(1234)
  x <- sim_mat(40, 30, perc_total_na = 0.1, perc_col_na = 1)$input
  expect_true(anyNA(x))
  real_kernel <- pca_imp_internal_cpp
  plant <- NULL
  local_mocked_bindings(
    pca_imp_internal_cpp = function(...) {
      res <- real_kernel(...)
      plant(res)
    }
  )
  run <- function(nb.init = 1) {
    suppressWarnings(pca_imp(
      x,
      ncp = 2,
      solver = "exact",
      miniter = 2,
      maxiter = 5,
      seed = 1,
      nb.init = nb.init
    ))
  }

  # the stub itself changes nothing: the premise of every case below
  plant <- function(res) res
  expect_false(anyNA(run()))

  for (bad in c(NaN, NA, Inf, -Inf)) {
    plant <- function(res) {
      res$mse <- bad
      res
    }
    expect_error(
      run(),
      "^slideimp internal, pca_imp: the fit's mse is not finite \\(",
      class = "simpleError"
    )
  }

  # checked on every initialization, not only the first
  n_call <- 0L
  plant <- function(res) {
    n_call <<- n_call + 1L
    if (n_call == 2L) {
      res$mse <- NaN
    }
    res
  }
  expect_error(run(nb.init = 3), "at initialization 2$")
  expect_identical(n_call, 2L)

  # a NaN value under a finite mse cannot come from the kernel, but if it
  # did, the part's validator refuses it before it reaches the result
  plant <- function(res) {
    expect_gt(length(res$value), 0L)
    res$value[1L] <- NaN
    res
  }
  expect_error(run(), "slideimp internal, part: 'values' holds NA")
})

test_that("pca_imp returns on replicate samples with ncp above their rank", {
  # Regression test. 30 samples of which 6 are distinct, holes included,
  # leave the centered completed data with rank 5 at most, so each call below
  # asks dsyevr for a range of eigenvalues that cuts a group of equal ones.
  # Up to 1.2.0 the chain gave dsyevr one eigenvalue cell per eigenpair asked
  # for, where LAPACK documents the order of the matrix, and dsyevr wrote
  # past them: each of these calls ended the R process (measured on Windows;
  # dev/review_p2.md at d5a34e5, finding 1). A build with that defect fails
  # here by taking the test process down, not by an expectation.
  set.seed(21)
  base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  x <- base6[rep(seq_len(6), times = 5), ]
  rownames(x) <- sprintf("s%02d", seq_len(nrow(x)))
  expect_true(anyNA(x))
  expect_identical(nrow(unique(x)), 6L)

  cases <- list(
    list(ncp = 5L, method = "regularized", k_eig = 6L),
    list(ncp = 7L, method = "regularized", k_eig = 8L),
    list(ncp = 12L, method = "regularized", k_eig = 13L),
    list(ncp = 7L, method = "EM", k_eig = 7L)
  )
  for (case in cases) {
    label <- sprintf("ncp = %d, method = %s", case$ncp, case$method)
    set.seed(1)
    res <- pca_imp(x, ncp = case$ncp, method = case$method)
    expect_identical(dim(res), dim(x), label = label)
    expect_false(anyNA(res), label = label)
    # the premise: more eigenpairs asked for than the completed data has rank
    m <- matrix(as.vector(res), nrow(res))
    centered <- sweep(m, 2, colMeans(m))
    expect_lt(qr(centered)$rank, case$k_eig, label = label)
  }
})

test_that("criterion_final is NA when the objective alone stopped the fit", {
  # Regression test. With ncp at or above the rank of the data the fit is
  # exact, and its objective is rounding residue below `threshold` by the
  # first iteration the stop is tested at, `miniter`. The criterion,
  # abs(1 - objective / old), is then a ratio of two such residues: it was
  # NaN, Inf or 0.5 by rounding, and moved with the build. Rounding can also
  # put that ratio below `threshold`, and then it is reported, so what holds
  # on every build is that the fit stops at `miniter` and reports NA or a
  # criterion below `threshold`, never NaN or a larger number.
  set.seed(21)
  base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  x <- base6[rep(seq_len(6), times = 5), ]
  rownames(x) <- sprintf("s%02d", seq_len(nrow(x)))
  for (ncp in c(5L, 6L, 7L, 12L)) {
    set.seed(1)
    res <- pca_imp(x, ncp = ncp)
    crit <- attr(res, "criterion_final")
    expect_true(attr(res, "converged"), label = sprintf("ncp = %d", ncp))
    expect_identical(attr(res, "n_iter"), 5L, label = sprintf("ncp = %d", ncp))
    expect_true(
      identical(crit, NA_real_) || isTRUE(crit < 1e-6),
      label = sprintf("ncp = %d, criterion_final %s", ncp, format(crit))
    )
  }

  # The NA branch, on an input where it does not rest on rounding: full rank,
  # but scaled down so far that the objective is below `threshold` from the
  # first iteration while the criterion is a well-defined ratio far above it.
  # The fit stops at `miniter` on the objective alone and reports NA.
  set.seed(4)
  y <- sim_mat(n = 30, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  tiny <- y * 1e-5
  set.seed(1)
  res <- pca_imp(tiny, ncp = 2, scale = FALSE)
  expect_true(attr(res, "converged"))
  expect_identical(attr(res, "n_iter"), 5L)
  expect_identical(attr(res, "criterion_final"), NA_real_)
  # the premise: at that count the criterion did not stop the fit. The same
  # five iterations under threshold = 0, which nothing reaches, report it.
  set.seed(1)
  expect_warning(
    at_miniter <- pca_imp(
      tiny,
      ncp = 2,
      scale = FALSE,
      threshold = 0,
      maxiter = 5,
      miniter = 5
    ),
    "Stopped after"
  )
  expect_false(attr(at_miniter, "converged"))
  expect_gt(attr(at_miniter, "criterion_final"), 100 * 1e-6)

  # A fit the criterion stopped still reports it, below the threshold.
  set.seed(4)
  y <- sim_mat(n = 30, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  set.seed(1)
  res <- pca_imp(y, ncp = 2, threshold = 1e-6)
  expect_true(attr(res, "converged"))
  expect_lt(attr(res, "criterion_final"), 1e-6)

  # So does a fit that ran out of iterations, above it.
  set.seed(1)
  res <- suppressWarnings(pca_imp(y, ncp = 2, threshold = 1e-12, maxiter = 6))
  expect_false(attr(res, "converged"))
  expect_gt(attr(res, "criterion_final"), 1e-12)
})

test_that("finite input whose squares overflow stops before the eigensolver", {
  # Values near 1e200 are finite, so the scan passes them, and their
  # cross-products are not: the Gram holds Inf. The kernel tests the Gram
  # where it is formed and stops before either solver is entered. This pins
  # the stop a user with such data sees, class and wording.
  set.seed(1234)
  x0 <- sim_mat(30, 10, perc_total_na = 0.1, perc_col_na = 1)$input
  x <- x0 * 1e200
  expect_true(anyNA(x))
  expect_true(all(is.finite(x[!is.na(x)])))

  # the premise: the same data at its own scale is fitted
  expect_false(anyNA(pca_imp(x0, ncp = 2, scale = FALSE)))

  for (solver in c("exact", "lobpcg")) {
    expect_error(
      pca_imp(x, ncp = 2, scale = FALSE, solver = solver),
      paste0(
        "the cross-product of the data is not finite; values this large ",
        "overflow when squared, rescale the input"
      ),
      fixed = TRUE,
      class = "Rcpp::exception",
      label = solver
    )
  }
})

test_that("pca_imp() takes unnamed input", {
  set.seed(1)
  obj <- sim_mat(n = 40, p = 12, perc_total_na = 0.1, perc_col_na = 0.5)$input
  obj <- unname(obj)
  expect_true(anyNA(obj))
  out <- pca_imp(obj, ncp = 2, seed = 1)
  expect_s3_class(out, "slideimp_results")
  expect_false(anyNA(out))
  expect_null(dimnames(out))
})
