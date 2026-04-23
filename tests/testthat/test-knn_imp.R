test_that("`impute_knn_brute` and `impute_knn_mlpack` calculate the missing location correctly", {
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  miss <- is.na(to_test)
  cmiss <- colSums(miss)
  miss_rate <- cmiss / nrow(to_test)
  # same preprocessing knn_imp() does before calling the low-level functions
  colmax <- 0.9
  eligible <- miss_rate < colmax
  pre_imp_cols <- to_test[, eligible, drop = FALSE]
  pre_imp_miss <- miss[, eligible, drop = FALSE]

  local_has_miss <- which(cmiss[eligible] > 0L)
  grp_impute <- as.integer(local_has_miss - 1L)

  expected_local <- unname(which(pre_imp_miss, arr.ind = TRUE))

  imputed_brute <- impute_knn_brute(
    obj = pre_imp_cols,
    k = 5,
    grp_impute = grp_impute,
    grp_miss_no_imp = integer(0L),
    grp_complete = integer(0L),
    method = 0L,
    dist_pow = 1,
    cores = 1
  )

  imputed_mlpack <- impute_knn_mlpack(
    obj = pre_imp_cols,
    k = 5,
    grp_impute = grp_impute,
    grp_miss_no_imp = integer(0L),
    grp_complete = integer(0L),
    method = 0L,
    dist_pow = 1,
    cores = 1
  )

  # extract only the location columns (row, local_col_1based) that the C++ already returns
  idx_brute <- imputed_brute[, 1:2, drop = FALSE]
  idx_mlpack <- imputed_mlpack[, 1:2, drop = FALSE]

  # sort so the comparison is order-independent
  expected_sorted <- expected_local[order(expected_local[, 1], expected_local[, 2]), , drop = FALSE]
  brute_sorted <- idx_brute[order(idx_brute[, 1], idx_brute[, 2]), , drop = FALSE]
  mlpack_sorted <- idx_mlpack[order(idx_mlpack[, 1], idx_mlpack[, 2]), , drop = FALSE]

  expect_equal(brute_sorted, expected_sorted)
  expect_equal(mlpack_sorted, expected_sorted)
})

test_that("`knn_imp` works", {
  obj <- sim_mat(50, 100)$input
  expect_no_error(knn_imp(
    obj,
    k = 3,
    method = "euclidean"
  ))

  expect_no_error(knn_imp(
    obj,
    k = 3,
    method = "manhattan",
    tree = TRUE
  ))
})

test_that("`knn_imp` tree and brute is the same for few missing values", {
  set.seed(1234)
  to_test <- sim_mat(20, n = 1000, perc_total_na = 0, perc_col_na = 0)$input
  to_test[1, 1] <- NA
  to_test[2, 2] <- NA

  expect_identical(
    knn_imp(to_test, k = 3, method = "euclidean"),
    knn_imp(to_test, k = 3, method = "euclidean", tree = TRUE)
  )
})

test_that("Exactly replicate `impute.knn`", {
  skip("Manual Testing Only")
  # library(impute)
  # set.seed(1234)
  # # post_imp behavior can cause differences
  # obj <- sim_mat(100, 100, perc_total_na = 0.05)$input
  #
  # # Check if the 'impute' package is installed
  #
  # # Perform imputation using knn_imp with method "impute.knn" on transposed data
  # r1 <- knn_imp(obj, k = 3, method = "euclidean", post_imp = FALSE)
  #
  # # Perform imputation using the original impute.knn function
  # # Transpose the result to match the orientation
  # r2 <- t(
  #   impute.knn(
  #     t(obj),
  #     k = 3,
  #     maxp = ncol(obj)
  #   )$data
  # )
  #
  # # Verify that the results from knn_imp match exactly with impute.knn
  # expect_equal(r1[, ], r2[, ])
})

test_that("`subset` feature of `knn_imp` works with post_imp = FALSE/TRUE", {
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.2, perc_col_na = 1)$input
  # Impute just 3 columns
  ## Check subset using numeric index
  r1 <- knn_imp(to_test, k = 3, post_imp = FALSE, subset = c(1, 3, 5))
  expect_true(!anyNA(r1[, c(1, 3, 5)]))
  expect_equal(is.na(r1[, -c(1, 3, 5)]), is.na(to_test[, -c(1, 3, 5)]))
  ## Check subset using character vector
  r2 <- knn_imp(
    to_test,
    k = 3,
    post_imp = FALSE,
    subset = paste0("feature", c(1, 3, 5))
  )
  expect_equal(r1, r2)

  # Test with post_imp = TRUE and a column requiring post imputation
  to_test_post <- to_test
  # Column 5 will be colMeans if post_imp is TRUE
  to_test_post[2:nrow(to_test_post), 5] <- NA
  r3 <- knn_imp(to_test_post, k = 3, post_imp = TRUE, subset = c(1, 3, 5))
  expect_true(!anyNA(r3[, c(1, 3, 5)]))
  # Expect that only the subset columns are imputed. The rests are untouched
  expect_equal(is.na(r3[, -c(1, 3, 5)]), is.na(to_test_post[, -c(1, 3, 5)]))
  # Verify post_imp on column 5
  col5_mean <- mean(to_test_post[, 5], na.rm = TRUE)
  expect_equal(unname(r3[, 5]), rep(col5_mean, nrow(to_test_post)))
  r4 <- knn_imp(
    to_test_post,
    k = 3,
    post_imp = TRUE,
    subset = paste0("feature", c(1, 3, 5))
  )
  expect_equal(r3, r4)
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  # row 1 is all NA
  to_test[1, ] <- NA
  expect_no_error(knn_imp(to_test, k = 3, post_imp = FALSE))
  expect_true(!anyNA(knn_imp(to_test, k = 3, post_imp = TRUE)))

  to_test[, 1] <- NA
  expect_error(knn_imp(to_test, k = 3, post_imp = FALSE), "All NA/Inf")

  # not admissible
  mat <- matrix(NA, nrow = 20, ncol = 20)
  diag(mat) <- rnorm(20)
  expect_error(knn_imp(mat, k = 2), "exceeds usable columns")
})
