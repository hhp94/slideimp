test_that("`impute_knn_brute` and `impute_knn_mlpack` calculate the missing location correctly", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)
  miss <- is.na(to_test)
  cmiss <- colSums(miss)
  miss_rate <- cmiss / nrow(to_test)
  # same preprocessing knn_imp() does before calling the low-level functions
  colmax <- 0.9
  eligible <- miss_rate < colmax
  pre_imp_cols <- to_test[, eligible, drop = FALSE]
  pre_imp_miss <- miss[, eligible, drop = FALSE]
  pre_imp_cols[pre_imp_miss] <- 0.0

  # local groups (0-based indices into pre_imp_cols)
  local_has_miss <- which(cmiss[eligible] > 0L)
  grp_impute <- as.integer(local_has_miss - 1L)

  # Expected missing positions *local to the submatrix pre_imp_cols* (1-based)
  expected_local <- unname(which(pre_imp_miss, arr.ind = TRUE))

  # brute-force path
  imputed_brute <- impute_knn_brute(
    obj = pre_imp_cols,
    nmiss = !pre_imp_miss,
    k = 5,
    grp_impute = grp_impute,
    grp_miss_no_imp = integer(0L),
    grp_complete = integer(0L),
    method = 0L,
    dist_pow = 1,
    cache = FALSE,
    cores = 1
  )

  # mlpack path
  imputed_mlpack <- impute_knn_mlpack(
    obj = mean_imp_col(pre_imp_cols),
    nmiss = !pre_imp_miss,
    k = 5,
    grp_impute = grp_impute,
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
  data("khanmiss1")
  expect_no_error(knn_imp(
    t(khanmiss1),
    k = 3,
    method = "euclidean"
  ))

  expect_no_error(knn_imp(
    t(khanmiss1),
    k = 3,
    method = "manhattan",
    tree = TRUE
  ))
})

test_that("`knn_imp` cache and non cache path is the same", {
  data("khanmiss1")
  expect_identical(
    knn_imp(t(khanmiss1), k = 10, method = "euclidean"),
    knn_imp(t(khanmiss1), k = 10, method = "euclidean", max_cache = 0)
  )
})

test_that("`knn_imp` tree and brute is the same for few missing values", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 1000, n = 20, perc_NA = 0, perc_col_NA = 0)$input)
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
  # data("khanmiss1")
  #
  # # Check if the 'impute' package is installed
  #
  # # Perform imputation using knn_imp with method "impute.knn" on transposed data
  # r1 <- knn_imp(t(khanmiss1), k = 3, method = "euclidean")
  #
  # # Perform imputation using the original impute.knn function
  # # Transpose the result to match the orientation
  # r2 <- t(
  #   impute.knn(
  #     khanmiss1,
  #     k = 3,
  #     maxp = nrow(khanmiss1)
  #   )$data
  # )
  #
  # # Verify that the results from knn_imp match exactly with impute.knn
  # expect_equal(r1[, ], r2[, ])
  #
  # # Test to see if the post_imp strategy would replicate the results completely
  # # Set seed for reproducibility in simulation
  # set.seed(1234)
  #
  # # Generate a simulated matrix with missing values (500 rows, 30 columns, 50% NA, 80% columns with NA)
  # to_test <- t(sim_mat(n = 500, m = 30, perc_NA = 0.5, perc_col_NA = 0.8)$input)
  #
  # # Pre-compute row means before imputation (ignoring NAs)
  # pre_impute <- rowMeans(to_test, na.rm = TRUE)
  #
  # # Impute using knn_imp without post-imputation step; expect some NAs to remain
  # r1.1 <- knn_imp(to_test, k = 5, method = "euclidean", post_imp = FALSE)
  # expect_true(anyNA(r1.1))
  #
  # # impute.knn uses the pre-imputation row means to impute the data.
  # # After knn_imp, we row impute the data with pre-calculated row_means
  # # Identify indices of remaining NAs
  # indices <- which(is.na(r1.1), arr.ind = TRUE)
  #
  # # Fill remaining NAs with pre-computed row means
  # r1.1[indices] <- pre_impute[indices[, 1]]
  #
  # # Verify no NAs remain after manual post-imputation
  # expect_true(!anyNA(r1.1))
  #
  # # Perform imputation using impute.knn on the transposed simulated data
  # r2.1 <- t(
  #   impute.knn(
  #     t(to_test),
  #     k = 5,
  #     maxp = ncol(to_test)
  #   )$data
  # )
  #
  # # Verify that the manually post-imputed knn_imp matches impute.knn
  # expect_equal(r1.1[, ], r2.1[, ])
  #
  # # Test subset. strategy is to use subset, then impute.knn on the same data
  # # and pull out the same subset then compare the two matrices
  # # Set seed for reproducibility in subset selection
  # set.seed(2345)
  #
  # # Generate another simulated matrix (100 rows, 200 columns, 10% NA, all columns with NA)
  # to_test_subset <- t(sim_mat(n = 100, m = 200, perc_NA = 0.1, perc_col_NA = 1)$input)
  #
  # # Randomly select 10 subset columns
  # subset_cols <- sample(colnames(to_test_subset), size = 10)
  #
  # # Verify that the subset has NAs before imputation
  # expect_true(anyNA(to_test_subset[, subset_cols]))
  #
  # # Impute only the subset columns using knn_imp without post_imp
  # r1_subset <- knn_imp(
  #   to_test_subset,
  #   k = 10,
  #   method = "euclidean",
  #   post_imp = FALSE,
  #   subset = subset_cols
  # )[, subset_cols]
  #
  # # Verify no NAs remain in the imputed subset
  # expect_true(!anyNA(r1_subset))
  #
  # # Perform full imputation using impute.knn and extract the subset
  # r2_subset <- t(
  #   impute.knn(
  #     t(to_test_subset),
  #     k = 10,
  #     maxp = ncol(to_test_subset)
  #   )$data
  # )[, subset_cols]
  #
  # # Verify no NAs in the extracted subset from full imputation
  # expect_true(!anyNA(r2_subset))
  #
  # # Verify that the subset imputation matches the extracted subset from full imputation
  # expect_equal(r1_subset, r2_subset)
})

test_that("`subset` feature of `knn_imp` works with post_imp = FALSE/TRUE", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)
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
    subset = paste0("feat", c(1, 3, 5))
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
    subset = paste0("feat", c(1, 3, 5))
  )
  expect_equal(r3, r4)
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 20, perc_NA = 0.2, perc_col_NA = 1)$input)
  # row 1 is all NA
  to_test[1, ] <- NA
  expect_no_error(knn_imp(to_test, k = 3, post_imp = FALSE))
  expect_true(!anyNA(knn_imp(to_test, k = 3, post_imp = TRUE)))

  to_test[, 1] <- NA
  expect_error(knn_imp(to_test, k = 3, post_imp = FALSE))
})
