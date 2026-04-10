test_that("`mean_imp_col` works", {
  set.seed(4321)
  to_test <- sim_mat(10, 20, perc_total_na = 0.5, perc_col_na = 1)$input
  c_manual <- to_test
  r_manual <- to_test
  na_indices <- which(is.na(to_test), arr.ind = TRUE)
  column_means <- colMeans(to_test, na.rm = TRUE)
  row_means <- rowMeans(to_test, na.rm = TRUE)

  c_manual[na_indices] <- column_means[na_indices[, 2]]
  expect_equal(mean_imp_col(to_test), c_manual)

  ## Test subset feature
  c_subset <- to_test
  for (i in c(1, 5, 10)) {
    c_subset[is.na(c_subset[, i]), i] <- mean(c_subset[, i], na.rm = TRUE)
  }
  expect_equal(mean_imp_col(to_test, subset = c(1, 5, 10)), c_subset)
})
