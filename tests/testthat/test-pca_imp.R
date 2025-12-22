test_that("same results as imputePCA", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)
  r1 <- missMDA::imputePCA(to_test, ncp = 2, nb.init = 10, seed = 1234)
  set.seed(1234)
  r2 <- pca_imp(to_test, ncp = 2, nb.init = 10, seed = 1234)
  expect_equal(r1$completeObs, r2[, ])

  row.w <- runif(nrow(to_test))
  row.w <- row.w / sum(row.w)

  set.seed(1234)
  r3 <- missMDA::imputePCA(to_test, ncp = 2, row.w = row.w, nb.init = 5, seed = 1234)
  set.seed(1234)
  r4 <- pca_imp(to_test, ncp = 2, nb.init = 5, row.w = row.w, seed = 1234)
  expect_equal(r3$completeObs, r4[, ])
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)
  to_test[1, ] <- NA
  expect_no_error(pca_imp(to_test, ncp = 2, seed = 1234))
  to_test[, 1] <- NA
  expect_error(pca_imp(to_test, ncp = 2, seed = 1234))
})

test_that("row.w = 'n_miss' matches missMDA::imputePCA with equivalent weights", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)

  # compute expected weights manually
  miss <- is.na(to_test)
  n_miss_per_row <- rowSums(miss)
  expected_w <- 1 - (n_miss_per_row / ncol(to_test))
  expected_w[expected_w < 1e-8] <- 1e-8
  expected_w <- expected_w / sum(expected_w)

  # compare "n_miss" shortcut against missMDA with explicit weights
  set.seed(1234)
  r1 <- missMDA::imputePCA(to_test, ncp = 2, nb.init = 5, row.w = expected_w, seed = 1234)
  set.seed(1234)
  r2 <- pca_imp(to_test, ncp = 2, nb.init = 5, row.w = "n_miss", seed = 1234)
  expect_equal(r1$completeObs, r2[, ])
  r3 <- pca_imp(to_test, ncp = 2, nb.init = 5, row.w = expected_w, seed = 1234)
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

  expect_error(pca_imp(mat, ncp = 2, row.w = "invalid"))
})
