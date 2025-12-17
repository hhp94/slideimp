test_that("same results as imputePCA", {
  skip_if_not_installed("missMDA")
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)
  r1 <- missMDA::imputePCA(to_test, ncp = 2, nb.init = 10, seed = 1234)
  set.seed(1234)
  r2 <- pca_imp(to_test, ncp = 2, nb.init = 10, seed = 1234)
  expect_equal(r1$completeObs, r2[, ])
  # print(
  #   bench::mark(
  #     missMDA::imputePCA(to_test, ncp = 2, nb.init = 10, seed = 1234),
  #     pca_imp(to_test, ncp = 2, nb.init = 10, seed = 1234),
  #     memory = FALSE,
  #     check = F,
  #     min_iterations = 10
  #   )
  # )
})

test_that("fastSVD_triplet is the same as svd.triplet", {
  skip_if_not_installed("FactoMineR")
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0)$input)
  row_weight <- runif(nrow(to_test))
  row_weight <- row_weight / sum(row_weight)
  r1 <- FactoMineR::svd.triplet(to_test, row.w = row_weight, ncp = 3)
  r2 <- fastSVD_triplet(to_test, row_w = row_weight, ncp = 3)
  expect_equal(r1$U, r2$U)
  expect_equal(r1$V, r2$V)
  expect_equal(r1$vs, r2$vs[, 1])
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)
  to_test[1, ] <- NA
  expect_no_error(pca_imp(to_test, ncp = 2, seed = 1234))
  to_test[, 1] <- NA
  expect_error(pca_imp(to_test, ncp = 2, seed = 1234))
})
