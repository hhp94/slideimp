skip_on_cran()

test_that("col_vars matches apply(var, na.rm = TRUE)", {
  set.seed(1234)
  obj <- sim_mat(10, 50, perc_col_na = 1)$input
  # the two documented non-values: a column with Inf and >= 2 values is NaN,
  # a column with one value is NA
  obj[1, 1] <- Inf
  obj[, 2] <- NA
  obj[1, 2] <- 0.5
  got <- col_vars(obj)
  want <- apply(obj, 2, var, na.rm = TRUE)
  expect_identical(names(got), names(want))
  # is.na() is TRUE for NaN and waldo takes NaN as equal to NA, so the two are
  # told apart by is.nan()
  expect_true(is.nan(got[[1]]))
  expect_true(is.na(got[[2]]) && !is.nan(got[[2]]))
  expect_identical(is.nan(got), is.nan(want))
  # measured at 2.8e-17 absolute and 3.5e-16 relative over the 48 finite
  # columns
  expect_worst_err(got, want, abs_tol = 1e-14, rel_tol = 1e-13)
})
