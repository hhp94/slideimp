test_that("col_vars works", {
  expect_equal(apply(khanmiss1, 2, var, na.rm = TRUE), col_vars(khanmiss1))
})
