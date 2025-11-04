test_that("same results as imputePCA", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)
  r1 <- missMDA::imputePCA(to_test, quanti.sup = c(1, 2, 3), ncp = 2, nb.init = 10, seed = 1234)
  set.seed(1234)
  r2 <- pca_imp(to_test, quanti.sup = c(1, 2, 3), ncp = 2, nb.init = 10, seed = 1234)
  expect_equal(r1, r2)
  # print(
  #   bench::mark(
  #     missMDA::imputePCA(to_test, ncp = 2, nb.init = 10, seed = 1234),
  #     pca_imp(to_test, ncp = 2, nb.init = 10, seed = 1234),
  #     memory = FALSE,
  #     min_iterations = 10
  #   )
  # )
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 20, perc_NA = 0.2, perc_col_NA = 1)$input)
  to_test[1, ] <- NA
  expect_no_error(pca_imp(to_test, ncp = 2, seed = 1234))
  to_test[, 1] <- NA
  expect_error(pca_imp(to_test, ncp = 2, seed = 1234))
})
