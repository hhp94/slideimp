test_that("grouped result is correct with aux columns, knn", {
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
  grouped_results <- group_imp(obj, group = group_df, k = 5)

  # manual imputation for comparison
  sub1_cols <- unique(c(group_df$features[[1]], group_df$aux[[1]]))
  sub1 <- knn_imp(obj[, sub1_cols], k = 5, subset = group_df$features[[1]])
  sub2_cols <- unique(c(group_df$features[[2]], group_df$aux[[2]]))
  sub2 <- knn_imp(obj[, sub2_cols], k = 5, subset = group_df$features[[2]])

  expected_results <- cbind(sub1, sub2)[, colnames(obj)]
  # Compare results
  expect_identical(grouped_results, expected_results)
})

test_that("grouped result is correct with aux columns, pca", {
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
  grouped_results <- group_imp(obj, group = group_df, ncp = 2)

  # manual imputation for comparison
  sub1_cols <- unique(c(group_df$features[[1]], group_df$aux[[1]]))
  sub1 <- pca_imp(obj[, sub1_cols], ncp = 2)[, group_df$features[[1]]]
  sub2_cols <- unique(c(group_df$features[[2]], group_df$aux[[2]]))
  sub2 <- pca_imp(obj[, sub2_cols], ncp = 2)[, group_df$features[[2]]]

  expected_results <- cbind(sub1, sub2)
  # Compare results
  expect_identical(grouped_results[, colnames(expected_results)], expected_results)
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
      list(k = 3, dist_pow = 0),
      list(k = 7, dist_pow = 1)
    )
  )

  obj <- t(to_test$input)
  grouped_results <- group_imp(obj, group = group_df, k = 5)

  # Manual verification with different parameters
  sub1 <- knn_imp(obj[, group_1], k = 3, subset = group_1[1:3], dist_pow = 0)
  sub2 <- knn_imp(obj[, group_2], k = 7, subset = group_2[1:4], dist_pow = 1)
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
    group_imp(obj, group = group_df, k = 5),
    "Same features can't be in more than 1 groups"
  )
})

test_that("grouped imputation works without aux columns, knn", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

  # no aux columns, only features
  group_df <- tibble::tibble(
    features = list(group_1[1:20], group_2[1:10])
  )

  obj <- t(to_test$input)
  expect_no_error(grouped_results <- group_imp(obj, group = group_df, k = 5))

  # Build expected results: start with original and update only imputed columns
  sub1 <- knn_imp(obj[, group_1[1:20]], k = 5)
  sub2 <- knn_imp(obj[, group_2[1:10]], k = 5)

  expected_results <- obj
  expected_results[, group_1[1:20]] <- sub1
  expected_results[, group_2[1:10]] <- sub2

  expect_identical(grouped_results, expected_results)
})

test_that("grouped imputation works without aux columns, pca", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

  # no aux columns, only features
  group_df <- tibble::tibble(
    features = list(group_1[1:20], group_2[1:10])
  )

  obj <- t(to_test$input)
  expect_no_error(grouped_results <- group_imp(obj, group = group_df, ncp = 2, seed = 1234))

  # Build expected results: start with original and update only imputed columns
  sub1 <- pca_imp(obj[, group_1[1:20]], ncp = 2, seed = 1234)
  sub2 <- pca_imp(obj[, group_2[1:10]], ncp = 2, seed = 1234)

  expected_results <- obj
  expected_results[, group_1[1:20]] <- sub1
  expected_results[, group_2[1:10]] <- sub2

  expect_identical(grouped_results, expected_results)
})
