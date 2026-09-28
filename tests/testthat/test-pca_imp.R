test_that("same results as imputePCA", {
  skip_if_not_manual()
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(
    20,
    50,
    perc_total_na = 0.25,
    perc_col_na = 1,
    rho = 0.75
  )$input
  expect_true(anyNA(to_test))
  # expected orientation (wide)
  r1 <- missMDA::imputePCA(to_test, ncp = 2, nb.init = 1, seed = 1234)
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 1,
    seed = 1234,
    solver = "exact",
    colmax = 1
  )
  expect_equal(r1$completeObs, r2[,])

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missMDA::imputePCA(
    to_test,
    ncp = 2,
    row.w = row.w,
    nb.init = 5,
    seed = 1234
  )
  set.seed(1234)
  r4 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = row.w,
    seed = 1234,
    solver = "exact"
  )
  expect_equal(r3$completeObs, r4[,])

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missMDA::imputePCA(to_test_t, ncp = 2, nb.init = 10, seed = 1234)
  set.seed(1234)
  r2_t <- pca_imp(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    solver = "exact"
  )
  expect_equal(r1_t$completeObs, r2_t[,])

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missMDA::imputePCA(
    to_test_t,
    ncp = 2,
    row.w = row.w_t,
    nb.init = 5,
    seed = 1234
  )
  set.seed(1234)
  r4_t <- pca_imp(
    to_test_t,
    ncp = 2,
    nb.init = 5,
    row.w = row.w_t,
    seed = 1234,
    solver = "exact"
  )
  expect_equal(r3_t$completeObs, r4_t[,])
})

test_that("same results as imputePCA, method = 'EM'", {
  skip_if_not_manual()
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(
    20,
    50,
    perc_total_na = 0.25,
    perc_col_na = 1,
    rho = 0.75
  )$input
  expect_true(anyNA(to_test))
  # expected orientation (wide)
  r1 <- missMDA::imputePCA(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_equal(r1$completeObs, r2[,])

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missMDA::imputePCA(
    to_test,
    ncp = 2,
    row.w = row.w,
    nb.init = 5,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r4 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = row.w,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_equal(r3$completeObs, r4[,])

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missMDA::imputePCA(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r2_t <- pca_imp(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_equal(r1_t$completeObs, r2_t[,])

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missMDA::imputePCA(
    to_test_t,
    ncp = 2,
    row.w = row.w_t,
    nb.init = 5,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r4_t <- pca_imp(
    to_test_t,
    ncp = 2,
    nb.init = 5,
    row.w = row.w_t,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_equal(r3_t$completeObs, r4_t[,])
})

test_that("same results as imputePCA, scale = FALSE", {
  skip_if_not_manual()
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(
    20,
    50,
    perc_total_na = 0.25,
    perc_col_na = 1,
    rho = 0.75
  )$input
  expect_true(anyNA(to_test))

  # expected orientation (wide)
  r1 <- missMDA::imputePCA(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_equal(r1$completeObs, r2[,])

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missMDA::imputePCA(
    to_test,
    ncp = 2,
    row.w = row.w,
    nb.init = 5,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r4 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = row.w,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_equal(r3$completeObs, r4[,])

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missMDA::imputePCA(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r2_t <- pca_imp(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_equal(r1_t$completeObs, r2_t[,])

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missMDA::imputePCA(
    to_test_t,
    ncp = 2,
    row.w = row.w_t,
    nb.init = 5,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r4_t <- pca_imp(
    to_test_t,
    ncp = 2,
    nb.init = 5,
    row.w = row.w_t,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_equal(r3_t$completeObs, r4_t[,])
})

test_that("row.w = 'n_miss' matches missMDA::imputePCA with equivalent weights", {
  skip_if_not_manual()
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(
    20,
    50,
    perc_total_na = 0.25,
    perc_col_na = 1,
    rho = 0.75
  )$input

  # compute expected weights manually
  miss <- is.na(to_test)
  n_miss_per_row <- rowSums(miss)
  expected_w <- 1 - (n_miss_per_row / ncol(to_test))
  expected_w[expected_w < 1e-8] <- 1e-8
  expected_w <- expected_w / sum(expected_w)

  # compare "n_miss" shortcut against missMDA with explicit weights
  set.seed(1234)
  r1 <- missMDA::imputePCA(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = expected_w,
    seed = 1234
  )
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = "n_miss",
    seed = 1234,
    solver = "exact"
  )
  expect_equal(r1$completeObs, r2[,])

  set.seed(1234)
  r3 <- pca_imp(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = expected_w,
    seed = 1234,
    solver = "exact"
  )
  expect_equal(r2, r3)
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- sim_mat(
    20,
    50,
    perc_total_na = 0.25,
    perc_col_na = 1,
    rho = 0.75
  )$input
  to_test[1, ] <- NA
  expect_no_error(pca_imp(to_test, ncp = 2, seed = 1234))
  to_test[, 1] <- NA
  expect_error(pca_imp(to_test, ncp = 2, seed = 1234))
  expect_true(all(is.na(to_test[, 1])))
})

test_that("row.w = 'n_miss' floors near-zero weights", {
  set.seed(42)
  # create matrix where one row has almost all missing
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(mat) <- paste0("row", 1:10)
  colnames(mat) <- paste0("col", 1:10)
  mat[1, -1] <- NA # row 1 has 9/10 missing -> weight = 0.1
  mat[2, ] <- NA
  mat[2, 1] <- rnorm(1) # row 2 has 9/10 missing -> weight = 0.1
  mat[3, 1:5] <- NA # row 3 has 5/10 missing -> weight = 0.5

  expect_no_error(pca_imp(mat, ncp = 2, row.w = "n_miss", seed = 123))
})

test_that("row.w rejects invalid strings", {
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

test_that("pca_imp falls back to mean imputation when ncp > usable eligible columns", {
  set.seed(1234)
  to_test <- sim_mat(30, 8, perc_total_na = 0.1, perc_col_na = 0.3)$input
  # this column will excceed colmax
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
    "exceeds the maximum usable components"
  )
})

test_that("pca_imp clamps imputed values to specified bounds", {
  set.seed(1234)
  to_test <- sim_mat(30, 8, perc_total_na = 0.1)$input
  na_pos <- which(is.na(to_test), arr.ind = TRUE)
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

test_that("post_imp fills NAs in ineligible columns with column mean", {
  set.seed(1234)
  to_test <- sim_mat(30, 8, perc_total_na = 0.05)$input
  to_test[1:28, 1] <- NA # column 1 exceeds colmax = 0.9

  res_with <- pca_imp(
    to_test,
    ncp = 2,
    seed = 1234,
    colmax = 0.9,
    post_imp = TRUE
  )
  res_without <- pca_imp(
    to_test,
    ncp = 2,
    seed = 1234,
    colmax = 0.9,
    post_imp = FALSE
  )

  expect_false(anyNA(res_with))
  expect_true(anyNA(res_without[, 1]))
  # filled values should equal column mean of remaining observed
  filled <- res_with[is.na(to_test[, 1]), 1]
  expect_true(all(filled == mean(to_test[, 1], na.rm = TRUE)))
})

test_that("pca_imp doesn't mess up the original object", {
  set.seed(1234)
  to_test <- sim_mat(30, 30, perc_total_na = 0.1, perc_col_na = 1)$input
  expect_true(anyNA(to_test))
  passed_obj <- to_test
  res <- pca_imp(
    to_test,
    ncp = 3,
    nb.init = 3,
    seed = 1234,
    colmax = 0.9,
    post_imp = FALSE
  )
  expect_equal(passed_obj, to_test)
  expect_equal(is.na(to_test), is.na(passed_obj))
})

test_that("pca_imp restores object even on bad input", {
  set.seed(1234)
  to_test <- sim_mat(30, 30, perc_total_na = 0.1, perc_col_na = 1)$input
  passed_obj <- to_test
  expect_true(anyNA(to_test))
  expect_error(pca_imp(
    to_test,
    ncp = 9999,
    nb.init = 3,
    seed = 1234,
    colmax = 0.9,
    post_imp = FALSE
  ))

  expect_equal(to_test, passed_obj)
})

test_that("Throw on Inf", {
  set.seed(1234)
  to_test <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  to_test[1, 1] <- Inf
  expect_error(pca_imp(to_test, ncp = 3), "Infinite")
  to_test[1, 1] <- -Inf
  expect_error(pca_imp(to_test, ncp = 3), "Infinite")
})

test_that("Inf is refused ahead of the paths that used to hide it", {
  set.seed(1234)

  # no missing cells at all: the "No missing values" early return used to hand
  # the matrix back with the Inf still in it.
  complete <- sim_mat(20, 20, perc_total_na = 0)$input
  complete[3, 4] <- Inf
  expect_error(pca_imp(complete, ncp = 2), "Infinite")

  # otherwise infeasible: the infeasibility abort used to win the race and
  # report eligible columns without ever mentioning the Inf.
  infeasible <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  infeasible[3, 4] <- Inf
  infeasible[, 6:20] <- 1
  expect_error(pca_imp(infeasible, ncp = 8), "Infinite")

  # and the Inf scan runs once, not once per `nb.init`
  expect_error(pca_imp(complete, ncp = 2, nb.init = 3), "Infinite")
})

test_that("pca_imp reports whether the auto verdict came from a probe", {
  set.seed(1234)
  # n_gram = min(nrow, n_elig) = 60, below the 250 threshold, so the R-side
  # heuristic demotes this to exact before the kernel is ever called and no
  # probe runs. `solver_chosen` still says "exact", because that is what ran.
  small <- sim_mat(60, 60, perc_total_na = 0.2, perc_col_na = 1)$input

  res <- pca_imp(small, ncp = 2, solver = "auto", seed = 1, na_check = FALSE)
  expect_identical(attr(res, "solver_chosen"), "exact")
  expect_false(attr(res, "solver_probed"))

  # a forced solver is not a verdict either: auto never ran
  for (s in c("exact", "lobpcg")) {
    r <- pca_imp(small, ncp = 2, solver = s, seed = 1, na_check = FALSE)
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
    seed = 1,
    na_check = FALSE
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
      post_imp = FALSE,
      na_check = FALSE
    )
  }

  w <- expect_warning(run(), "Stopped after")

  # Regression: this warning used to be raised by Rcpp::warning() inside the EM
  # loop, which for a std::string is a bare Rf_warning. Under options(warn = 2),
  # or any calling handler that invokes a restart, R longjmps out of the C++
  # frame; END_RCPP catches exceptions, not longjmps, so the Armadillo working
  # matrices, the Rcout precision guard and the wrapper's RNGScope were all
  # skipped. A cli/rlang condition is the proof that R raises it now.
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
      post_imp = FALSE,
      na_check = FALSE
    )
  )
})

# argument validation ----
test_that("a tampered lobpcg_control object is re-validated, not taken as-is", {
  set.seed(1234)
  x <- sim_mat(30, 40, perc_total_na = 0.2, perc_col_na = 1)$input

  run <- function(ctrl) {
    pca_imp(
      x, ncp = 2, solver = "lobpcg", lobpcg_control = ctrl,
      miniter = 2, maxiter = 5, seed = 1, na_check = FALSE
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

  # above INT_MAX the C++ backend aborts on the uword conversion; catch it here
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

test_that("pca_imp validates seed and the miniter/maxiter ordering from R", {
  set.seed(1234)
  x <- sim_mat(20, 20, perc_total_na = 0.1, perc_col_na = 1)$input

  # maxiter is deliberately small everywhere below, so the non-convergence
  # warning is expected and not what any of these assertions are about.
  run <- function(...) {
    suppressWarnings(pca_imp(x, ncp = 2, solver = "exact", na_check = FALSE, ...))
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
    seed = .Machine$integer.max, nb.init = 1L, miniter = 2, maxiter = 5
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

  run <- function(m, post_imp) {
    suppressWarnings(pca_imp(
      m, ncp = 2, solver = "exact", miniter = 2, maxiter = 5, seed = 1,
      post_imp = post_imp
    ))
  }

  # Regression: pca_imp() ran mean_imp_col(obj) over every column whenever
  # post_imp = TRUE. On an all-eligible input that copied the matrix twice
  # more (the C++ result, then Rcpp's wrap of it) and scanned every column,
  # to change nothing. knn_imp() already restricted the call to the columns
  # that could still hold NA and skipped it when there were none.
  real_mean_imp_col <- mean_imp_col
  calls <- list()
  local_mocked_bindings(
    mean_imp_col = function(obj, subset = NULL, cores = 1) {
      # wrapped in list() so that a NULL subset - the whole-matrix call this
      # test exists to catch - is recorded rather than dropped by `[[<-`
      calls <<- c(calls, list(subset))
      real_mean_imp_col(obj, subset = subset, cores = cores)
    }
  )

  # every column of x is eligible, so the kernel fills every NA and the post
  # step has nothing to do. Check that premise before asserting on it.
  res_no_post <- run(x, post_imp = FALSE)
  expect_false(anyNA(res_no_post))
  res_post <- run(x, post_imp = TRUE)
  expect_length(calls, 0L)
  expect_identical(as.vector(res_post), as.vector(res_no_post))
  expect_true(attr(res_post, "post_imp"))

  # column 1 above colmax, column 2 zero variance with NA, column 3 zero
  # variance without NA. The first two are the only columns that can still
  # hold NA after the kernel; the third is ineligible but has nothing to fill.
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
  expect_false(anyNA(res_y))
  expect_equal(
    unname(res_y[is.na(y[, 1]), 1]),
    rep(mean(y[, 1], na.rm = TRUE), n_na_1)
  )
  expect_identical(unname(res_y[, 2]), rep(7, 40L))

  # and the result is what the unconditional call used to produce
  old_way <- real_mean_imp_col(matrix(as.vector(run(y, post_imp = FALSE)), 40L))
  expect_identical(as.vector(res_y), as.vector(old_way))
})

# lobpcg ----
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
  expect_lt(max_abs_diff(got$mat, ref$mat), 1e-10)

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
    expect_lt(max_abs_diff(got$mat, ref$mat), 1e-4)

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
  run <- function(x, solver) {
    pca_imp(
      x,
      ncp = 3,
      scale = FALSE,
      solver = solver,
      seed = 1,
      nb.init = 1,
      maxiter = 200,
      miniter = 5,
      threshold = 1e-8,
      post_imp = FALSE,
      na_check = FALSE,
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

test_that("LOBPCG fallback to exact still produces correct result", {
  set.seed(1234)
  x <- sim_mat(60, 30)$input
  pca_iters <- 10L
  warmup <- 2L

  # kernel-direct, so no non-convergence warning. See the note above.
  ref <- run_pca_fixed_iters(
    x,
    solver = "exact",
    pca_iters = pca_iters
  )

  # tol below machine precision + maxiter = 1 forces LOBPCG to fail every
  # post-warmup iteration, exercising the exact fallback path.
  got <- run_pca_fixed_iters(
    x,
    solver = "lobpcg",
    ctrl = lobpcg_control(
      maxiter = 1L,
      warmup_iters = warmup,
      tol = 1e-20
    ),
    pca_iters = pca_iters
  )

  expect_false(anyNA(got$mat))
  expect_lt(max_abs_diff(got$mat, ref$mat), 1e-4)

  # every post-warmup iteration should attempt LOBPCG, fail, then fall back to exact.
  expect_equal(got$n_lobpcg_ok, 0L)
  expect_equal(got$n_lobpcg_bad, pca_iters - warmup)
  expect_equal(got$n_exact, pca_iters) # warmup exact + fallback exact
})
