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
      post_imp = FALSE
    )
  )
})

# argument validation ----
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

  # Regression: pca_imp() ran mean_imp_col(obj) over every column whenever
  # post_imp = TRUE. On an all-eligible input that copied the matrix twice
  # more (the C++ result, then Rcpp's wrap of it) and scanned every column,
  # to change nothing. The mean part now covers the ineligible columns only,
  # so fallback_part() is asked for exactly those, and for none on an
  # all-eligible input.
  real_mean_imp_col <- mean_imp_col
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
  expect_identical(as.vector(res_post), as.vector(res_no_post))
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

  # and the result is what the unconditional call used to produce
  old_way <- real_mean_imp_col(matrix(as.vector(run(y, post_imp = FALSE)), 40L))
  expect_identical(as.vector(res_y), as.vector(old_way))
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
  expect_error(run(), "slideimp internal, part: `values` holds NA")
})

test_that("pca_imp returns on replicate samples with ncp above their rank", {
  # Regression test. 30 samples of which 6 are distinct, holes included,
  # leave the centered completed data with rank 5 at most, so each call below
  # asks dsyevr for a range of eigenvalues that cuts a group of equal ones.
  # Up to 1.2.0 the chain gave dsyevr one eigenvalue cell per eigenpair asked
  # for, where LAPACK documents the order of the matrix, and dsyevr wrote
  # past them: each of these calls ended the R process (measured on Windows;
  # dev/review_p2.md, finding 1). A build with that defect fails here by
  # taking the test process down, not by an expectation.
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
  # exact and stops because its objective is below `threshold`. The
  # criterion, abs(1 - objective / old), is then a ratio of two values near
  # zero: it was NaN, Inf or 0.5 by rounding, and moved with the build.
  # Rounding can also put that ratio below `threshold`, and then it is
  # reported, so what holds on every build is that a converged fit reports
  # NA or a criterion below `threshold`, never NaN or a larger number.
  set.seed(21)
  base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  x <- base6[rep(seq_len(6), times = 5), ]
  rownames(x) <- sprintf("s%02d", seq_len(nrow(x)))
  for (ncp in c(5L, 6L, 7L, 12L)) {
    set.seed(1)
    res <- pca_imp(x, ncp = ncp)
    crit <- attr(res, "criterion_final")
    expect_true(attr(res, "converged"), label = sprintf("ncp = %d", ncp))
    expect_true(
      identical(crit, NA_real_) || isTRUE(crit < 1e-6),
      label = sprintf("ncp = %d, criterion_final %s", ncp, format(crit))
    )
  }

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
  # Values near 1e200 are finite, so check_finite() passes them, and their
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

test_that("wide recover_factors block agrees with the stored reference", {
  # SVD_triplet()'s wide branch forms V = Xhat' * Vn only for the changing
  # (NA-containing) columns inside the EM loop, and the fixed columns' rows
  # once after it (src/svd_triplet.h, src/armaSVD.cpp postprocessing). The
  # fixed rows reach R only through `mse`, the nb.init selection statistic;
  # the imputed values never read them. So both are compared, for a wide
  # input with fixed columns (Gram cache active), a wide input with an NA in
  # every column (no fixed block, cache inactive), and a tall input, whose
  # path is unchanged.
  #
  # The reference was recorded from the code before the change, at 96dec4f,
  # by a generator script committed alongside it in 7d8508c
  # (`git show 7d8508c:dev/gen_pca_recover_block_ref.R`). It is bytes from
  # another build: the comparison is on a tolerance, and bitwise agreement is
  # reported, not asserted.
  ref <- readRDS(test_path("fixtures", "pca_recover_block_ref.rds"))
  claims <- list(
    wide_fixed = list(wide = TRUE, fixed = TRUE),
    wide_all = list(wide = TRUE, fixed = FALSE),
    tall_fixed = list(wide = FALSE, fixed = TRUE)
  )
  expect_setequal(names(ref$cases), names(claims))

  for (nm in names(claims)) {
    case <- ref$cases[[nm]]
    x <- case$x
    n_fixed <- sum(colSums(is.na(x)) == 0L)

    # the input exercises the branch the case is named for
    expect_true(anyNA(x), label = nm)
    expect_identical(ncol(x) > nrow(x), claims[[nm]]$wide, label = nm)
    expect_identical(n_fixed > 0L, claims[[nm]]$fixed, label = nm)

    got <- pca_imp_internal_cpp(
      obj = x,
      eligible_idx = seq_len(ncol(x)) - 1L,
      ncp = ref$ncp,
      scale = TRUE,
      regularized = TRUE,
      threshold = 0,
      init = 0L,
      maxiter = ref$iters,
      miniter = ref$iters,
      row_w = case$row_w,
      coeff_ridge = 1,
      solver = 0L,
      warmup_iters = 0L,
      lobpcg_tol = 1e-8,
      lobpcg_maxiter = 1L
    )
    want <- case$imputed_values

    # structural: the same cells, in the same order. The reference holds
    # each cell as a row of doubles (row, col, value), as the kernel returned
    # them when it was recorded; the kernel now returns the two indices in
    # the integer matrix `loc` and the values in `value`.
    expect_identical(dim(got$loc), c(nrow(want), 2L), label = nm)
    expect_identical(colnames(got$loc), c("row", "col"), label = nm)
    expect_identical(length(got$value), nrow(want), label = nm)
    expect_identical(unname(got$loc[, 1]), as.integer(want[, 1]), label = nm)
    expect_identical(unname(got$loc[, 2]), as.integer(want[, 2]), label = nm)

    a <- got$value
    b <- want[, 3]
    expect_gt(length(b), 0L)
    err <- abs(a - b)
    # the relative figure is taken over the imputed cells whose reference
    # value is not near zero
    keep <- abs(b) >= 1e-3
    expect_gt(sum(keep), 0L)
    max_abs <- max(err)
    max_rel <- max(err[keep] / abs(b[keep]))
    expect_gt(case$mse, 0)
    mse_abs <- abs(got$mse - case$mse)
    mse_rel <- mse_abs / abs(case$mse)

    cat(sprintf(
      paste0(
        "\n%s: imputed max abs %.3g, max rel %.3g over %d of %d cells; ",
        "mse abs %.3g, rel %.3g; bitwise imputed %s, mse %s\n"
      ),
      nm,
      max_abs,
      max_rel,
      sum(keep),
      length(b),
      mse_abs,
      mse_rel,
      identical(a, b),
      identical(got$mse, case$mse)
    ))

    expect_lt(max_abs, 1e-8, label = nm)
    expect_lt(max_rel, 1e-8, label = nm)
    expect_lt(mse_abs, 1e-8, label = nm)
    expect_lt(mse_rel, 1e-8, label = nm)
  }
})

test_that("the PCA kernel returns integers and refuses what it cannot count with", {
  set.seed(1234)
  x <- sim_mat(20, 12, perc_total_na = 0.1, perc_col_na = 1)$input
  holes <- which(is.na(x), arr.ind = TRUE)
  expect_gt(nrow(holes), 0L)
  args <- list(
    obj = x,
    eligible_idx = seq_len(ncol(x)) - 1L,
    ncp = 2L,
    scale = TRUE,
    regularized = TRUE,
    threshold = 1e-6,
    init = 0L,
    maxiter = 20L,
    miniter = 5L,
    row_w = rep(1, nrow(x)),
    coeff_ridge = 1,
    solver = 0L,
    warmup_iters = 0L,
    lobpcg_tol = 1e-8,
    lobpcg_maxiter = 1L
  )
  got <- do.call(pca_imp_internal_cpp, args)

  # `loc` has one row per hole, 1-based, in column-major order, which is the
  # order which() lists them in
  expect_true(is.integer(got$loc))
  expect_identical(colnames(got$loc), c("row", "col"))
  expect_identical(got$loc[, "row"], unname(holes[, 1]))
  expect_identical(got$loc[, "col"], unname(holes[, 2]))
  expect_true(is.double(got$value))
  expect_null(dim(got$value))
  expect_length(got$value, nrow(holes))
  expect_false(anyNA(got$value))
  expect_true(is.integer(got$n_iter))
  expect_gte(got$n_iter, 5L)

  # a count below its range, NA included; as an unsigned index each of
  # these wrapped to a huge value
  bad <- function(...) {
    do.call(pca_imp_internal_cpp, utils::modifyList(args, list(...)))
  }
  expect_error(bad(maxiter = 0L), "maxiter must be >= 1")
  expect_error(bad(maxiter = NA_integer_), "maxiter must be >= 1")
  expect_error(bad(miniter = -1L), "miniter must be >= 0")
  expect_error(bad(miniter = 21L), "miniter must be <= maxiter")
  expect_error(bad(ncp = 0L), "ncp must be >= 1")
  expect_error(bad(ncp = NA_integer_), "ncp must be >= 1")
  expect_error(bad(init = -1L), "init must be >= 0")
  expect_error(bad(trace_iter = -1L), "trace_iter must be >= 0")
  expect_error(
    bad(solver = 1L, lobpcg_maxiter = 0L),
    "lobpcg_maxiter must be >= 1"
  )
  expect_error(
    bad(solver = 1L, warmup_iters = -1L),
    "warmup_iters must be >= 0"
  )

  # a column index outside the matrix, NA included
  expect_error(
    bad(eligible_idx = c(0L, ncol(x))),
    "eligible_idx holds 12 at position 2"
  )
  expect_error(
    bad(eligible_idx = c(0L, NA_integer_)),
    "eligible_idx holds NA at position 2"
  )
  expect_error(
    bad(eligible_idx = c(-1L, 0L)),
    "eligible_idx holds -1 at position 1"
  )
  expect_error(bad(row_w = rep(1, nrow(x) - 1L)), "row_w length must equal")
})
