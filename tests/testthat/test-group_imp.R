# group_imp ----
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
  expect_identical(grouped_results[, ], expected_results)
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
  grouped_results <- group_imp(obj, group = group_df, ncp = 2, seed = 1234)

  # manual imputation for comparison
  sub1_cols <- unique(c(group_df$features[[1]], group_df$aux[[1]]))
  sub1 <- pca_imp(obj[, sub1_cols], ncp = 2, seed = 1234)[, group_df$features[[1]]]
  sub2_cols <- unique(c(group_df$features[[2]], group_df$aux[[2]]))
  sub2 <- pca_imp(obj[, sub2_cols], ncp = 2, seed = 1234)[, group_df$features[[2]]]

  expected_results <- cbind(sub1, sub2)
  expect_equal(grouped_results[, colnames(expected_results)], expected_results)
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
  grouped_results <- group_imp(obj, group = group_df)

  # Manual verification with different parameters
  sub1 <- knn_imp(obj[, group_1], k = 3, subset = group_1[1:3], dist_pow = 0)
  sub2 <- knn_imp(obj[, group_2], k = 7, subset = group_2[1:4], dist_pow = 1)
  expected_results <- cbind(sub1, sub2)[, colnames(obj)]

  expect_identical(grouped_results[, ], expected_results)
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

  expect_identical(grouped_results[, ], expected_results)
})

test_that("group-specific parameters work correctly, pca", {
  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id
  # Different ncp and coeff.ridge values for each group
  group_df <- tibble::tibble(
    features = list(group_1[1:3], group_2[1:4]),
    aux = list(group_1, group_2),
    parameters = list(
      list(ncp = 2, coeff.ridge = 1, seed = 1234),
      list(ncp = 3, coeff.ridge = 2, seed = 1234)
    )
  )
  obj <- t(to_test$input)
  grouped_results <- group_imp(obj, group = group_df)

  # Manual verification with different parameters
  sub1_full <- pca_imp(obj[, group_1], ncp = 2, coeff.ridge = 1, seed = 1234)
  sub1 <- obj[, group_1]
  sub1[, group_1[1:3]] <- sub1_full[, group_1[1:3]]
  sub2_full <- pca_imp(obj[, group_2], ncp = 3, coeff.ridge = 2, seed = 1234)
  sub2 <- obj[, group_2]
  sub2[, group_2[1:4]] <- sub2_full[, group_2[1:4]]
  expected_results <- cbind(sub1, sub2)[, colnames(obj)]
  expect_equal(grouped_results[, ], expected_results)
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

  expect_identical(grouped_results[, ], expected_results)
})

test_that("group-specific parameters work correctly in parallel, pca", {
  skip_on_cran()
  skip_on_ci()

  set.seed(1234)
  to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)
  group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
  group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id
  # Different ncp and coeff.ridge values for each group
  group_df <- tibble::tibble(
    features = list(group_1[1:3], group_2[1:4]),
    aux = list(group_1, group_2),
    parameters = list(
      list(ncp = 2, coeff.ridge = 1),
      list(ncp = 3, coeff.ridge = 2)
    )
  )
  obj <- t(to_test$input)
  mirai::daemons(2, seed = 1234)
  grouped_results <- group_imp(obj, group = group_df, cores = 2, seed = 1234, nb.init = 10)
  mirai::daemons(0)
  # Manual verification with different parameters
  sub1_full <- pca_imp(obj[, group_1], ncp = 2, coeff.ridge = 1, seed = 1234, nb.init = 10)
  sub1 <- obj[, group_1]
  sub1[, group_1[1:3]] <- sub1_full[, group_1[1:3]]
  sub2_full <- pca_imp(obj[, group_2], ncp = 3, coeff.ridge = 2, seed = 1234, nb.init = 10)
  sub2 <- obj[, group_2]
  sub2[, group_2[1:4]] <- sub2_full[, group_2[1:4]]
  expected_results <- cbind(sub1, sub2)[, colnames(obj)]

  imputed_cols <- c(group_1[1:3], group_2[1:4])
  obj_orig <- obj[, imputed_cols]
  grouped_values <- grouped_results[, imputed_cols][is.na(obj_orig)]
  expected_values <- expected_results[, imputed_cols][is.na(obj_orig)]
  # seeding in parallel is hard to reproduce correctly
  expect_true(
    cor(grouped_values, expected_values) > 0.999
  )
})

# group_features ----
test_that("group_features returns correct structure without k/ncp", {
  obj <- matrix(1:10, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e")))
  features_df <- data.frame(
    feature_id = c("a", "b", "c", "d"),
    group = c("g1", "g1", "g2", "g2")
  )

  result <- group_features(obj, features_df)

  expect_true(tibble::is_tibble(result))
  expect_true("group" %in% names(result))
  expect_true("features" %in% names(result))
  expect_equal(sort(result$group), c("g1", "g2"))
})

test_that("group_features handles subset correctly", {
  obj <- matrix(1:12, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e", "f")))
  features_df <- data.frame(
    feature_id = c("a", "b", "c", "d", "e", "f"),
    group = c("g1", "g1", "g1", "g2", "g2", "g2")
  )

  result <- group_features(obj, features_df, subset = c("a", "b", "d", "e"))

  expect_true("aux" %in% names(result))

  g1_row <- result[result$group == "g1", ]
  g2_row <- result[result$group == "g2", ]

  # features in subset go to features column
  expect_setequal(g1_row$features[[1]], c("a", "b"))
  expect_setequal(g2_row$features[[1]], c("d", "e"))

  # features not in subset go to aux column
  expect_equal(g1_row$aux[[1]], "c")
  expect_equal(g2_row$aux[[1]], "f")
})

test_that("group_features errors when no subset element found", {
  obj <- matrix(1:6, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
  features_df <- data.frame(
    feature_id = c("a", "b", "c"),
    group = c("g1", "g1", "g2")
  )

  expect_error(
    suppressWarnings(group_features(obj, features_df, subset = c("x", "y", "z"))),
    "No element"
  )
})

test_that("group_features adds parameters column with k", {
  obj <- matrix(1:10, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e")))
  features_df <- data.frame(
    feature_id = c("a", "b", "c", "d", "e"),
    group = c("g1", "g1", "g1", "g2", "g2")
  )

  result <- group_features(obj, features_df, k = 5)

  expect_true("parameters" %in% names(result))

  g1_row <- result[result$group == "g1", ]
  g2_row <- result[result$group == "g2", ]

  # k is capped at group_size - 1
  expect_equal(g1_row$parameters[[1]]$k, 2) # 3 features, so k = min(2, 5) = 2
  expect_equal(g2_row$parameters[[1]]$k, 1) # 2 features, so k = min(1, 5) = 1
})

test_that("group_features adds parameters column with ncp", {
  obj <- matrix(1:10, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e")))
  features_df <- data.frame(
    feature_id = c("a", "b", "c", "d", "e"),
    group = c("g1", "g1", "g1", "g2", "g2")
  )

  result <- group_features(obj, features_df, ncp = 5)

  expect_true("parameters" %in% names(result))

  g1_row <- result[result$group == "g1", ]
  g2_row <- result[result$group == "g2", ]

  # ncp is capped at min(ncol(obj), group_size - 1, ncp)
  expect_equal(g1_row$parameters[[1]]$ncp, 2) # min(5, 2, 5) = 2
  expect_equal(g2_row$parameters[[1]]$ncp, 1) # min(5, 1, 5) = 1
})

test_that("group_features pads groups to min_group_size", {
  obj <- matrix(1:12, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e", "f")))
  features_df <- data.frame(
    feature_id = c("a", "b", "c"),
    group = c("g1", "g1", "g2")
  )

  result <- group_features(obj, features_df, min_group_size = 4, seed = 123)

  expect_true("aux" %in% names(result))

  g1_row <- result[result$group == "g1", ]
  g2_row <- result[result$group == "g2", ]

  # g1 has 2 features, needs 2 more to reach min_group_size of 4
  expect_equal(length(g1_row$features[[1]]) + length(g1_row$aux[[1]]), 4)

  # g2 has 1 feature, needs 3 more
  expect_equal(length(g2_row$features[[1]]) + length(g2_row$aux[[1]]), 4)
})

test_that("group_features errors when min_group_size too large", {
  obj <- matrix(1:6, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
  features_df <- data.frame(
    feature_id = c("a", "b", "c"),
    group = c("g1", "g1", "g2")
  )

  expect_error(
    group_features(obj, features_df, min_group_size = 100),
    "too large"
  )
})

test_that("group_features errors when no colnames match features_df", {
  obj <- matrix(1:6, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
  features_df <- data.frame(
    feature_id = c("x", "y", "z"),
    group = c("g1", "g1", "g2")
  )

  expect_error(
    group_features(obj, features_df),
    "matched with no group"
  )
})
