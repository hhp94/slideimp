test_that("grouped result is correct with aux columns", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)

  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

  # impute only first 3 values of group 1, the rest are aux. Group 2 do 4 features.
  group_df <- tibble::tibble(
    features = list(group_1[1:3], group_2[1:4]),
    aux = list(group_1, group_2)
  )

  # run grouped imputation
  obj <- t(to_test$input)
  grouped_results <- group_knn_imp(obj, group = group_df, k = 5)[[1]]

  # manual imputation for comparison
  sub1_cols <- unique(c(group_df$features[[1]], group_df$aux[[1]]))
  sub1 <- knn_imp(obj[, sub1_cols], k = 5, subset = group_df$features[[1]])[[1]]
  sub2_cols <- unique(c(group_df$features[[2]], group_df$aux[[2]]))
  sub2 <- knn_imp(obj[, sub2_cols], k = 5, subset = group_df$features[[2]])[[1]]

  expected_results <- cbind(sub1, sub2)[, colnames(obj)]
  # Compare results
  expect_identical(grouped_results, expected_results)
})

test_that("grouped result is correct with n_imp = 2 and n_pmm = 0", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

  # impute only first 3 values of group 1, the rest are aux. Group 2 do 4 features.
  group_df <- tibble::tibble(
    features = list(group_1[1:3], group_2[1:4]),
    aux = list(group_1, group_2)
  )

  # run grouped imputation with n_imp = 2, n_pmm = 0
  obj <- t(to_test$input)
  grouped_results <- group_knn_imp(
    obj,
    group = group_df,
    k = 5,
    n_imp = 2,
    n_pmm = 0,
    seed = 1234
  )
  expect_length(grouped_results, 2)
  expect_identical(dim(grouped_results[[1]]), dim(obj))
  expect_identical(dim(grouped_results[[2]]), dim(obj))

  # ,anual imputation for comparison (for both imputations)
  sub1_cols <- unique(c(group_df$features[[1]], group_df$aux[[1]]))
  sub1_imp <- knn_imp(
    obj[, sub1_cols],
    k = 5,
    subset = group_df$features[[1]],
    n_imp = 2,
    n_pmm = 0,
    seed = 1234
  )

  sub2_cols <- unique(c(group_df$features[[2]], group_df$aux[[2]]))
  sub2_imp <- knn_imp(
    obj[, sub2_cols],
    k = 5,
    subset = group_df$features[[2]],
    n_imp = 2,
    n_pmm = 0,
    seed = 1234
  )

  # build expected results for both imputations
  expected_results_1 <- cbind(sub1_imp[[1]], sub2_imp[[1]])[, colnames(obj)]
  expected_results_2 <- cbind(sub1_imp[[2]], sub2_imp[[2]])[, colnames(obj)]

  # compare results for both imputations
  expect_identical(grouped_results[[1]], expected_results_1)
  expect_identical(grouped_results[[2]], expected_results_2)

  # verify that the two imputations are different
  imputed_cols <- c(group_df$features[[1]], group_df$features[[2]])
  expect_false(identical(
    grouped_results[[1]][, imputed_cols],
    grouped_results[[2]][, imputed_cols]
  ))

  # non-imputed columns remain unchanged across imputations
  non_imputed_cols <- setdiff(colnames(obj), imputed_cols)
  if (length(non_imputed_cols) > 0) {
    expect_identical(
      grouped_results[[1]][, non_imputed_cols],
      grouped_results[[2]][, non_imputed_cols]
    )
    expect_identical(
      grouped_results[[1]][, non_imputed_cols],
      obj[, non_imputed_cols]
    )
  }
})

test_that("group-specific parameters work correctly", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

  # Different k values for each group
  group_df <- tibble::tibble(
    features = list(group_1[1:3], group_2[1:4]),
    aux = list(group_1, group_2),
    parameters = list(
      list(k = 3, weighted = TRUE),
      list(k = 7, weighted = FALSE)
    )
  )

  obj <- t(to_test$input)
  grouped_results <- group_knn_imp(obj, group = group_df, k = 5)[[1]]

  # Manual verification with different parameters
  sub1 <- knn_imp(obj[, group_1], k = 3, subset = group_1[1:3], weighted = TRUE)[[1]]
  sub2 <- knn_imp(obj[, group_2], k = 7, subset = group_2[1:4], weighted = FALSE)[[1]]
  expected_results <- cbind(sub1, sub2)[, colnames(obj)]

  expect_identical(grouped_results, expected_results)
})

test_that("duplicate features across groups throws error", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

  group_df <- tibble::tibble(
    features = list(group_1[1:5], c(group_1[5], group_2[1:3])), # group_1[5] in both
    aux = list(group_1, group_2)
  )

  obj <- t(to_test$input)
  expect_error(
    group_knn_imp(obj, group = group_df, k = 5),
    "Same features can't be in more than 1 groups"
  )
})

test_that("grouped imputation works without aux columns", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

  # no aux columns, only features
  group_df <- tibble::tibble(
    features = list(group_1[1:20], group_2[1:10]),
    aux = list(NULL, NULL)
  )

  obj <- t(to_test$input)
  expect_no_error(grouped_results <- group_knn_imp(obj, group = group_df, k = 5))

  # Build expected results: start with original and update only imputed columns
  sub1 <- knn_imp(obj[, group_1[1:20]], k = 5)[[1]]
  sub2 <- knn_imp(obj[, group_2[1:10]], k = 5)[[1]]

  expected_results <- obj
  expected_results[, group_1[1:20]] <- sub1
  expected_results[, group_2[1:10]] <- sub2

  expect_identical(grouped_results[[1]], expected_results)
})
