test_that("col_vars works", {
  obj <- sim_mat(10, 50, perc_col_na = 1)$input
  expect_equal(apply(obj, 2, var, na.rm = TRUE), col_vars(obj))
})
