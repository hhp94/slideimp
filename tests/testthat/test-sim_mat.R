skip_on_cran()

test_that("sim_mat works", {
  expect_no_error(sim_mat())
  norm_mat <- as.vector(
    sim_mat(
      n = 10,
      p = 10,
      beta = TRUE,
      perc_total_na = 0,
      perc_col_na = 0
    )$input
  )
  expect_true(!any(norm_mat > 1 | norm_mat < 0))
})

test_that("NA count matches perc_total_na and every column keeps >=1 obs", {
  set.seed(1)
  # 50 NAs over 6 columns of 10 rows is close to the cap of n - 1 = 9 per
  # column, so some columns are filled to it and keep exactly one observation
  res <- sim_mat(n = 10, p = 10, perc_total_na = 0.5, perc_col_na = 0.6)
  expected_na <- ceiling(0.5 * 10 * 10)
  expect_equal(sum(is.na(res$input)), expected_na)
  n_na <- colSums(is.na(res$input))
  expect_identical(sum(n_na > 0), 6L)
  # the cap is reached, so the check below is made at its edge
  expect_identical(max(n_na), 9)
  # no column fully NA
  expect_true(all(colSums(!is.na(res$input)) >= 1))
})

test_that("perc_total_na = 0 produces no NAs", {
  res <- sim_mat(n = 20, p = 10, perc_total_na = 0, perc_col_na = 0.5)
  expect_equal(sum(is.na(res$input)), 0)
})

test_that("metadata has correct shape and names align with matrix", {
  res <- sim_mat(n = 15, p = 8, n_col_groups = 3, n_row_groups = 2)
  expect_equal(nrow(res$col_group), 8)
  expect_equal(nrow(res$row_group), 15)
  expect_setequal(res$col_group$feature, colnames(res$input))
  expect_setequal(res$row_group$sample, rownames(res$input))
  expect_true(all(res$col_group$group %in% paste0("group", 1:3)))
  expect_true(all(res$row_group$group %in% paste0("group", 1:2)))
})

test_that("print.slideimp_sim() truncates rows to n and input columns to p", {
  res <- sim_mat(n = 15, p = 8, n_col_groups = 3, n_row_groups = 2)
  out <- capture.output(v <- withVisible(print(res, n = 2, p = 3)))
  expect_false(v$visible)
  expect_identical(v$value, res)

  want <- c(
    "$col_group (3 column groups)",
    sprintf("# Showing 2 of %d rows", nrow(res$col_group)),
    "$row_group (2 row groups)",
    sprintf("# Showing 2 of %d rows", nrow(res$row_group)),
    sprintf("$input (%d x %d)", nrow(res$input), ncol(res$input)),
    sprintf(
      "# Showing 2 of %d rows and 3 of %d columns",
      nrow(res$input),
      ncol(res$input)
    )
  )
  expect_true(all(want %in% out))
  # columns 1 to 3 of input are shown, the 4th is not; no 3rd row anywhere
  has <- function(s) any(grepl(s, out, fixed = TRUE))
  cn <- colnames(res$input)
  expect_true(has(cn[3]))
  expect_false(has(cn[4]))
  expect_false(has(rownames(res$input)[3]))
})

test_that("aborts when perc_col_na is too small to fit requested NAs", {
  # n=10, p=10 -> max_per_col = 9. perc_col_na = 0.1 -> 1 NA-bearing column.
  # perc_total_na = 0.5 -> need 50 NAs, but 1 column holds at most 9.
  expect_error(
    sim_mat(n = 10, p = 10, perc_total_na = 0.5, perc_col_na = 0.1),
    "perc_col_na"
  )
})
