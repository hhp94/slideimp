skip_if_not_manual()

# pca_imp() against missMDA::imputePCA(), cell by cell on the worst case.
# Over every comparison below the two differ by at most 5.85e-14 absolute
# and 6.25e-13 relative (over the cells whose missMDA value is nonzero), so
# these bounds leave about 17 and 160 times that. One iteration more on the
# missMDA side moves them by 2e-5 or more.
missmda_abs_tol <- 1e-12
missmda_rel_tol <- 1e-10

# Both run a fixed number of EM iterations, so they take the same steps
# whatever rounding does to the criterion near convergence. imputePCA() has
# no `miniter`: under threshold = 0, which no criterion reaches, it stops at
# `maxiter` and warns that it did, and only that warning is muffled.
# pca_imp() fixes the count with miniter = maxiter, and its threshold, far
# above the criterion at that count, decides only `converged`.
fixed_iters <- 30
missmda_fixed <- function(...) {
  withCallingHandlers(
    missMDA::imputePCA(..., threshold = 0, maxiter = fixed_iters),
    warning = function(w) {
      if (startsWith(conditionMessage(w), "Stopped after ")) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
pca_fixed <- function(...) {
  pca_imp(..., threshold = 0.1, maxiter = fixed_iters, miniter = fixed_iters)
}

test_that("same results as imputePCA", {
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
  r1 <- missmda_fixed(to_test, ncp = 2, nb.init = 1, seed = 1234)
  set.seed(1234)
  r2 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 1,
    seed = 1234,
    solver = "exact",
    colmax = 1
  )
  expect_worst_err(r2[,], r1$completeObs, missmda_abs_tol, missmda_rel_tol)

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missmda_fixed(
    to_test,
    ncp = 2,
    row.w = row.w,
    nb.init = 5,
    seed = 1234
  )
  set.seed(1234)
  r4 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = row.w,
    seed = 1234,
    solver = "exact"
  )
  expect_worst_err(r4[,], r3$completeObs, missmda_abs_tol, missmda_rel_tol)

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missmda_fixed(to_test_t, ncp = 2, nb.init = 10, seed = 1234)
  set.seed(1234)
  r2_t <- pca_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    solver = "exact"
  )
  expect_worst_err(r2_t[,], r1_t$completeObs, missmda_abs_tol, missmda_rel_tol)

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missmda_fixed(
    to_test_t,
    ncp = 2,
    row.w = row.w_t,
    nb.init = 5,
    seed = 1234
  )
  set.seed(1234)
  r4_t <- pca_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 5,
    row.w = row.w_t,
    seed = 1234,
    solver = "exact"
  )
  expect_worst_err(r4_t[,], r3_t$completeObs, missmda_abs_tol, missmda_rel_tol)
})

test_that("same results as imputePCA, method = 'EM'", {
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
  r1 <- missmda_fixed(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r2 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_worst_err(r2[,], r1$completeObs, missmda_abs_tol, missmda_rel_tol)

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missmda_fixed(
    to_test,
    ncp = 2,
    row.w = row.w,
    nb.init = 5,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r4 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = row.w,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_worst_err(r4[,], r3$completeObs, missmda_abs_tol, missmda_rel_tol)

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missmda_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r2_t <- pca_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_worst_err(r2_t[,], r1_t$completeObs, missmda_abs_tol, missmda_rel_tol)

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missmda_fixed(
    to_test_t,
    ncp = 2,
    row.w = row.w_t,
    nb.init = 5,
    seed = 1234,
    method = "EM"
  )
  set.seed(1234)
  r4_t <- pca_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 5,
    row.w = row.w_t,
    seed = 1234,
    method = "EM",
    solver = "exact"
  )
  expect_worst_err(r4_t[,], r3_t$completeObs, missmda_abs_tol, missmda_rel_tol)
})

test_that("same results as imputePCA, scale = FALSE", {
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
  r1 <- missmda_fixed(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r2 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_worst_err(r2[,], r1$completeObs, missmda_abs_tol, missmda_rel_tol)

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missmda_fixed(
    to_test,
    ncp = 2,
    row.w = row.w,
    nb.init = 5,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r4 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = row.w,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_worst_err(r4[,], r3$completeObs, missmda_abs_tol, missmda_rel_tol)

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missmda_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r2_t <- pca_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 10,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_worst_err(r2_t[,], r1_t$completeObs, missmda_abs_tol, missmda_rel_tol)

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missmda_fixed(
    to_test_t,
    ncp = 2,
    row.w = row.w_t,
    nb.init = 5,
    seed = 1234,
    scale = FALSE
  )
  set.seed(1234)
  r4_t <- pca_fixed(
    to_test_t,
    ncp = 2,
    nb.init = 5,
    row.w = row.w_t,
    seed = 1234,
    scale = FALSE,
    solver = "exact"
  )
  expect_worst_err(r4_t[,], r3_t$completeObs, missmda_abs_tol, missmda_rel_tol)
})

test_that("row.w = 'n_miss' matches missMDA::imputePCA with equivalent weights", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(
    20,
    50,
    perc_total_na = 0.25,
    perc_col_na = 1,
    rho = 0.75
  )$input

  # compute expected weights manually, with the floor pca_imp() applies
  miss <- is.na(to_test)
  n_miss_per_row <- rowSums(miss)
  expected_w <- 1 - (n_miss_per_row / ncol(to_test))
  expected_w[expected_w < 1e-10] <- 1e-10
  expected_w <- expected_w / sum(expected_w)

  # compare "n_miss" shortcut against missMDA with explicit weights
  set.seed(1234)
  r1 <- missmda_fixed(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = expected_w,
    seed = 1234
  )
  set.seed(1234)
  r2 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = "n_miss",
    seed = 1234,
    solver = "exact"
  )
  expect_worst_err(r2[,], r1$completeObs, missmda_abs_tol, missmda_rel_tol)

  set.seed(1234)
  r3 <- pca_fixed(
    to_test,
    ncp = 2,
    nb.init = 5,
    row.w = expected_w,
    seed = 1234,
    solver = "exact"
  )
  # the shortcut and its weights given by hand: the same values to rounding,
  # and the same record but for `criterion_final`, a ratio of two near-equal
  # objectives that rounding moves
  expect_worst_err(r2[,], r3[,], missmda_abs_tol, missmda_rel_tol)
  same <- setdiff(names(attributes(r2)), "criterion_final")
  expect_identical(attributes(r2)[same], attributes(r3)[same])
})
