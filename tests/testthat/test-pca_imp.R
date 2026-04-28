test_that("same results as imputePCA", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.25, perc_col_na = 1, rho = 0.75)$input
  expect_true(anyNA(to_test))
  # expected orientation (wide)
  r1 <- missMDA::imputePCA(to_test, ncp = 2, nb.init = 1, seed = 1234)
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 1, seed = 1234,
    lobpcg_control = lobpcg_control(maxiter = 0),
    colmax = 1
  )
  expect_equal(r1$completeObs, r2[, ])

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missMDA::imputePCA(to_test, ncp = 2, row.w = row.w, nb.init = 5, seed = 1234)
  set.seed(1234)
  r4 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 5, row.w = row.w, seed = 1234,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r3$completeObs, r4[, ])

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missMDA::imputePCA(to_test_t, ncp = 2, nb.init = 10, seed = 1234)
  set.seed(1234)
  r2_t <- pca_imp(
    to_test_t,
    ncp = 2, nb.init = 10, seed = 1234,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r1_t$completeObs, r2_t[, ])

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missMDA::imputePCA(to_test_t, ncp = 2, row.w = row.w_t, nb.init = 5, seed = 1234)
  set.seed(1234)
  r4_t <- pca_imp(
    to_test_t,
    ncp = 2, nb.init = 5, row.w = row.w_t, seed = 1234,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r3_t$completeObs, r4_t[, ])
})

test_that("same results as imputePCA, method = 'EM'", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.25, perc_col_na = 1, rho = 0.75)$input
  expect_true(anyNA(to_test))
  # expected orientation (wide)
  r1 <- missMDA::imputePCA(to_test, ncp = 2, nb.init = 10, seed = 1234, method = "EM")
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 10, seed = 1234, method = "EM",
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r1$completeObs, r2[, ])

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missMDA::imputePCA(
    to_test,
    ncp = 2, row.w = row.w, nb.init = 5, seed = 1234, method = "EM"
  )
  set.seed(1234)
  r4 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 5, row.w = row.w, seed = 1234, method = "EM",
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r3$completeObs, r4[, ])

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missMDA::imputePCA(to_test_t, ncp = 2, nb.init = 10, seed = 1234, method = "EM")
  set.seed(1234)
  r2_t <- pca_imp(
    to_test_t,
    ncp = 2, nb.init = 10, seed = 1234, method = "EM",
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r1_t$completeObs, r2_t[, ])

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missMDA::imputePCA(to_test_t, ncp = 2, row.w = row.w_t, nb.init = 5, seed = 1234, method = "EM")
  set.seed(1234)
  r4_t <- pca_imp(
    to_test_t,
    ncp = 2, nb.init = 5, row.w = row.w_t, seed = 1234, method = "EM",
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r3_t$completeObs, r4_t[, ])
})

test_that("same results as imputePCA, scale = FALSE", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.25, perc_col_na = 1, rho = 0.75)$input
  expect_true(anyNA(to_test))

  # expected orientation (wide)
  r1 <- missMDA::imputePCA(to_test, ncp = 2, nb.init = 10, seed = 1234, scale = FALSE)
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 10, seed = 1234, scale = FALSE,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r1$completeObs, r2[, ])

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)
  set.seed(1234)
  r3 <- missMDA::imputePCA(to_test, ncp = 2, row.w = row.w, nb.init = 5, seed = 1234, scale = FALSE)
  set.seed(1234)
  r4 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 5, row.w = row.w, seed = 1234, scale = FALSE,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r3$completeObs, r4[, ])

  # transposed input also gives identical results
  set.seed(1234)
  to_test_t <- t(to_test)
  r1_t <- missMDA::imputePCA(to_test_t, ncp = 2, nb.init = 10, seed = 1234, scale = FALSE)
  set.seed(1234)
  r2_t <- pca_imp(
    to_test_t,
    ncp = 2, nb.init = 10, seed = 1234, scale = FALSE,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r1_t$completeObs, r2_t[, ])

  row.w_t <- runif(nrow(to_test_t))
  row.w_t <- row.w_t / sum(row.w_t)
  set.seed(1234)
  r3_t <- missMDA::imputePCA(to_test_t, ncp = 2, row.w = row.w_t, nb.init = 5, seed = 1234, scale = FALSE)
  set.seed(1234)
  r4_t <- pca_imp(
    to_test_t,
    ncp = 2, nb.init = 5, row.w = row.w_t, seed = 1234, scale = FALSE,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r3_t$completeObs, r4_t[, ])
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.25, perc_col_na = 1, rho = 0.75)$input
  to_test[1, ] <- NA
  expect_no_error(pca_imp(to_test, ncp = 2, seed = 1234))
  to_test[, 1] <- NA
  expect_error(pca_imp(to_test, ncp = 2, seed = 1234))
  expect_true(all(is.na(to_test[, 1])))
})

test_that("row.w = 'n_miss' matches missMDA::imputePCA with equivalent weights", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.25, perc_col_na = 1, rho = 0.75)$input

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
    ncp = 2, nb.init = 5, row.w = expected_w, seed = 1234
  )
  set.seed(1234)
  r2 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 5, row.w = "n_miss", seed = 1234,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r1$completeObs, r2[, ])

  set.seed(1234)
  r3 <- pca_imp(
    to_test,
    ncp = 2, nb.init = 5, row.w = expected_w, seed = 1234,
    lobpcg_control = lobpcg_control(maxiter = 0)
  )
  expect_equal(r2, r3)
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

  # Force ineligible columns:
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
  expect_gt(mean(is.na(res_no_post[, 1])), mean_1)
  expect_equal(unique(res_no_post[, 2]), 69)
  expect_equal(sum(is.na(res_no_post[, 3])), 3L)
  expect_false(anyNA(res_no_post[, 4:12]))
})

test_that("pca_imp falls back to mean imputation when ncp > usable eligible columns", {
  set.seed(1234)
  to_test <- sim_mat(30, 8, perc_total_na = 0.1, perc_col_na = 0.3)$input
  # This column will excceed colmax
  to_test[1:29, 1] <- NA
  expect_no_error(res <- pca_imp(
    to_test,
    ncp = 3,
    nb.init = 3,
    seed = 1234,
    colmax = 0.9,
    post_imp = FALSE
  ))
  # Make most columns ineligible (all-NA)
  to_test[, 1:6] <- NA
  for (i in 1:6) {
    to_test[sample.int(30, size = 1), i] <- rnorm(1)
  }

  # Only 2 eligible columns left -> ncp = 3 > min(28, 1) -> error (1 usable component)
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

test_that("LOBPCG enabled during warmup matches disabled exact path", {
  set.seed(1234)
  x <- sim_mat(20, 80)$input
  pca_iters <- 12L
  expect_warning(
    ref <- run_pca_fixed_iters(
      x,
      ctrl = lobpcg_control(maxiter = 0L, warmup_iters = 0L),
      pca_iters = pca_iters
    )
  )
  expect_warning(
    got <- run_pca_fixed_iters(
      x,
      ctrl = lobpcg_control(
        maxiter = 50L,
        warmup_iters = pca_iters + 1L,
        tol = 1e-10
      ),
      pca_iters = pca_iters
    )
  )
  expect_false(anyNA(ref$mat))
  expect_false(anyNA(got$mat))
  expect_equal(dim(got$mat), dim(ref$mat))

  # same numerical path (LOBPCG never triggered), so this should be tight.
  expect_lt(max_abs_diff(got$mat, ref$mat), 1e-10)

  # warmup_iters > pca_iters means LOBPCG is never attempted
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
    expect_warning(
      ref <- run_pca_fixed_iters(
        x,
        ctrl = lobpcg_control(maxiter = 0L, warmup_iters = 0L),
        pca_iters = pca_iters
      )
    )
    expect_warning(
      got <- run_pca_fixed_iters(
        x,
        ctrl = lobpcg_control(
          maxiter = 100L,
          warmup_iters = warmup,
          tol = 1e-8
        ),
        pca_iters = pca_iters
      )
    )
    expect_false(anyNA(ref$mat))
    expect_false(anyNA(got$mat))
    expect_equal(dim(got$mat), dim(ref$mat))
    expect_lt(max_abs_diff(got$mat, ref$mat), 1e-4)

    # ref: LOBPCG disabled, every iter is exact
    expect_equal(ref$n_lobpcg_ok, 0L)
    expect_equal(ref$n_lobpcg_bad, 0L)
    expect_equal(ref$n_exact, pca_iters)

    # got: warmup iters use exact, rest should converge via LOBPCG
    expect_equal(got$n_lobpcg_bad, 0L)
    expect_equal(got$n_exact, warmup)
    expect_equal(got$n_lobpcg_ok, pca_iters - warmup)
  }
})

test_that("LOBPCG fallback to exact still produces correct result", {
  set.seed(1234)
  x <- sim_mat(60, 30)$input
  pca_iters <- 10L
  warmup <- 2L

  expect_warning(
    ref <- run_pca_fixed_iters(
      x,
      ctrl = lobpcg_control(maxiter = 0L, warmup_iters = 0L),
      pca_iters = pca_iters
    )
  )
  # tol below machine precision + maxiter = 1 forces LOBPCG to fail every
  # post-warmup iteration, exercising the exact fallback path.
  expect_warning(
    got <- run_pca_fixed_iters(
      x,
      ctrl = lobpcg_control(
        maxiter = 1L,
        warmup_iters = warmup,
        tol = 1e-20
      ),
      pca_iters = pca_iters
    )
  )

  expect_false(anyNA(got$mat))
  expect_lt(max_abs_diff(got$mat, ref$mat), 1e-4)

  # every post-warmup iter should attempt LOBPCG, fail, fall back to exact
  expect_equal(got$n_lobpcg_ok, 0L)
  expect_equal(got$n_lobpcg_bad, pca_iters - warmup)
  expect_equal(got$n_exact, pca_iters) # warmup + fallbacks
})
