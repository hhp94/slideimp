test_that("col_vars matches apply(var, na.rm = TRUE)", {
  set.seed(1234)
  obj <- sim_mat(10, 50, perc_col_na = 1)$input
  obj[1, 1] <- Inf
  expect_equal(col_vars(obj), apply(obj, 2, var, na.rm = TRUE))
})
