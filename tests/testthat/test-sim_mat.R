test_that("sim_mat works", {
  expect_no_error(sim_mat())
  norm_mat <- as.vector(sim_mat(
    n = 10,
    p = 10,
    beta = TRUE,
    perc_total_na = 0,
    perc_col_na = 0
  )$input)
  expect_true(!any(norm_mat > 1 | norm_mat < 0))
})

test_that("NA count matches perc_total_na and every column keeps >=1 obs", {
  set.seed(1)
  res <- sim_mat(n = 50, p = 20, perc_total_na = 0.1, perc_col_na = 0.5)
  expected_na <- ceiling(0.1 * 50 * 20)
  expect_equal(sum(is.na(res$input)), expected_na)
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

test_that("aborts when perc_col_na is too small to fit requested NAs", {
  # n=10, p=10 -> max_per_col = 9. perc_col_na = 0.1 -> 1 NA-bearing column.
  # perc_total_na = 0.5 -> need 50 NAs, but 1 column holds at most 9.
  expect_error(
    sim_mat(n = 10, p = 10, perc_total_na = 0.5, perc_col_na = 0.1),
    "perc_col_na"
  )
})
