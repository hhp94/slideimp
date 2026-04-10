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
