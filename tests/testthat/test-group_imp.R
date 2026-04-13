# group_imp ----
test_that("group column + feature column API works correctly", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1, rho = 0.75)
  obj <- to_test$input

  # prep_groups to match the expected API
  group_long <- data.frame(
    group = to_test$col_group$group,
    feature = to_test$col_group$feature
  )

  # equivalent list-column form for comparison
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  group_list <- data.frame(
    feature = I(list(group_1, group_2))
  )

  expect_identical(
    group_imp(obj, group = group_long, k = 3),
    group_imp(obj, group = group_list, k = 3)
  )

  expect_identical(
    group_imp(obj, group = group_long, ncp = 5, nb.init = 10, seed = 1234),
    group_imp(obj, group = group_list, ncp = 5, nb.init = 10, seed = 1234)
  )
})

test_that("group column API collapses duplicate groups correctly", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  obj <- to_test$input

  # duplicate some rows so the same group label appears in separate blocks
  gf <- to_test$col_group
  group_doubled <- data.frame(
    group = c(gf$group, gf$group),
    feature = c(gf$feature, gf$feature)
  )

  group_single <- data.frame(
    group = gf$group,
    feature = gf$feature
  )

  expect_identical(
    group_imp(obj, group = group_doubled, k = 3),
    group_imp(obj, group = group_single, k = 3)
  )
})

test_that("group_imp() handles aux columns present in only some groups (padded)", {
  set.seed(1234)
  to_test <- sim_mat(n = 20, p = 50, n_col_groups = 2, perc_total_na = 0.3, perc_col_na = 1)
  obj <- to_test$input
  meta <- to_test$col_group

  # move feature 1 into its own tiny group so padding is needed for that
  # group only — the other groups will have zero aux columns.
  meta[1, "group"] <- "group3"

  prepped <- prep_groups(colnames(obj), group = meta, min_group_size = 10)

  # sanity: exactly one group has aux, the others have none. This is the
  # configuration that previously broke the rep(iter, aux_lengths) split.
  aux_lens <- lengths(prepped$aux)
  expect_true(sum(aux_lens > 0) == 1L)
  expect_true(any(aux_lens == 0L))

  # should run without "subscript out of bounds" and return a full matrix.
  expect_no_error(
    res <- group_imp(obj, group = prepped, k = 3)
  )
  expect_identical(dim(res), dim(obj))
  expect_identical(colnames(res), colnames(obj))
  expect_false(anyNA(res))
})

test_that("group column API accepts factor group column", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  obj <- to_test$input

  group_df <- data.frame(
    group = factor(to_test$col_group$group),
    feature = to_test$col_group$feature
  )

  expect_no_error(group_imp(obj, group = group_df, k = 2))
})

# --- Preconditioning tests for the `group` parameter ---
test_that("group must be a data.frame", {
  obj <- matrix(1:12, nrow = 3, dimnames = list(NULL, paste0("f", 1:4)))
  expect_error(
    group_imp(obj, group = list(feature = list("f1")), k = 2),
    "group"
  )
})

test_that("group must contain a 'feature' column", {
  obj <- matrix(1:12, nrow = 3, dimnames = list(NULL, paste0("f", 1:4)))
  group_df <- data.frame(stuff = I(list(c("f1", "f2"))))
  expect_error(group_imp(obj, group = group_df, k = 2), "feature")
})

test_that("group errors when feature is character without group column", {
  obj <- matrix(1:12, nrow = 3, dimnames = list(NULL, paste0("f", 1:4)))
  group_df <- data.frame(feature = c("f1", "f2"))
  expect_error(
    group_imp(obj, group = group_df, k = 2),
    "no group column"
  )
})

test_that("group errors on NA in group column", {
  obj <- matrix(1:12, nrow = 3, dimnames = list(NULL, paste0("f", 1:4)))
  group_df <- data.frame(
    group = c("a", NA_character_),
    feature = c("f1", "f2")
  )
  expect_error(group_imp(obj, group = group_df, k = 2), "NA")
})

test_that("group feature column must not have NAs", {
  obj <- matrix(1:12, nrow = 3, dimnames = list(NULL, paste0("f", 1:4)))
  group_df <- data.frame(
    group = c("a", "a"),
    feature = c("f1", NA_character_)
  )
  expect_error(group_imp(obj, group = group_df, k = 2), "feature")
})

test_that("grouped result is correct with aux columns, knn", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)

  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # impute only first 3 values of group 1, the rest are aux. Group 2 do 4 feature.
  group_df <- data.frame(
    feature = I(list(group_1[1:3], group_2[1:4])),
    aux = I(list(group_1, group_2))
  )

  # run grouped imputation
  obj <- to_test$input
  grouped_results <- group_imp(obj, group = group_df, k = 3)

  # manual imputation for comparison
  sub1_cols <- unique(c(group_df$feature[[1]], group_df$aux[[1]]))
  sub1 <- knn_imp(obj[, sub1_cols], k = 3, subset = group_df$feature[[1]])
  sub2_cols <- unique(c(group_df$feature[[2]], group_df$aux[[2]]))
  sub2 <- knn_imp(obj[, sub2_cols], k = 3, subset = group_df$feature[[2]])

  expected_results <- cbind(sub1, sub2)[, colnames(obj)]
  # Compare results
  expect_identical(grouped_results[, ], expected_results)
})

test_that("grouped result is correct with aux columns, pca", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)

  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # impute only first 3 values of group 1, the rest are aux. Group 2 do 4 feature.
  group_df <- data.frame(
    feature = I(list(group_1[1:3], group_2[1:4])),
    aux = I(list(group_1, group_2))
  )

  # run grouped imputation
  obj <- to_test$input
  grouped_results <- group_imp(obj, group = group_df, ncp = 2, seed = 1234)

  # manual imputation for comparison
  sub1_cols <- unique(c(group_df$feature[[1]], group_df$aux[[1]]))
  sub1 <- pca_imp(obj[, sub1_cols], ncp = 2, seed = 1234)[, group_df$feature[[1]]]
  sub2_cols <- unique(c(group_df$feature[[2]], group_df$aux[[2]]))
  sub2 <- pca_imp(obj[, sub2_cols], ncp = 2, seed = 1234)[, group_df$feature[[2]]]

  expected_results <- cbind(sub1, sub2)
  expect_equal(grouped_results[, colnames(expected_results)], expected_results)
})

test_that("group-specific parameters work correctly", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # Different k values for each group
  group_df <- data.frame(
    feature = I(list(group_1[1:3], group_2[1:4])),
    aux = I(list(group_1, group_2)),
    parameters = I(list(
      list(k = 3, dist_pow = 0),
      list(k = 7, dist_pow = 1)
    ))
  )

  obj <- to_test$input
  grouped_results <- group_imp(obj, group = group_df)

  # Manual verification with different parameters
  sub1 <- knn_imp(obj[, group_1], k = 3, subset = group_1[1:3], dist_pow = 0)
  sub2 <- knn_imp(obj[, group_2], k = 7, subset = group_2[1:4], dist_pow = 1)
  expected_results <- cbind(sub1, sub2)[, colnames(obj)]

  expect_identical(grouped_results[, ], expected_results)
})

test_that("duplicate feature across groups throws error", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  group_df <- data.frame(
    feature = I(list(group_1[1:5], c(group_1[5], group_2[1:3]))), # group_1[5] in both
    aux = I(list(group_1, group_2))
  )

  obj <- to_test$input
  expect_error(
    group_imp(obj, group = group_df, k = 3),
    "appear in more than one group"
  )
})

test_that("grouped imputation works without aux columns, knn", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # no aux columns, only feature
  group_df <- data.frame(
    feature = I(list(group_1[1:5], group_2[6:10]))
  )
  obj <- to_test$input

  grouped_results <- group_imp(
    obj,
    group = group_df,
    k = 3,
    allow_unmapped = TRUE
  )

  # Build expected results: start with original and update only imputed columns
  sub1 <- knn_imp(obj[, group_1[1:5]], k = 3)
  sub2 <- knn_imp(obj[, group_2[6:10]], k = 3)

  expected_results <- obj
  expected_results[, group_1[1:5]] <- sub1
  expected_results[, group_2[6:10]] <- sub2

  expect_identical(grouped_results[, ], expected_results)
})

test_that("group-specific parameters work correctly, pca", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # Different ncp and coeff.ridge values for each group
  group_df <- data.frame(
    feature = I(list(group_1[1:3], group_2[1:4])),
    aux = I(list(group_1, group_2)),
    parameters = I(
      list(
        list(ncp = 2, coeff.ridge = 1, seed = 1234),
        list(ncp = 3, coeff.ridge = 2, seed = 1234)
      )
    )
  )
  obj <- to_test$input
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
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # no aux columns, only feature
  group_df <- data.frame(
    feature = I(list(group_1[1:5], group_2[6:10]))
  )
  obj <- to_test$input

  grouped_results <- group_imp(
    obj,
    group = group_df,
    ncp = 2,
    seed = 1234,
    allow_unmapped = TRUE
  )

  # Build expected results: start with original and update only imputed columns
  sub1 <- pca_imp(obj[, group_1[1:5]], ncp = 2, seed = 1234)
  sub2 <- pca_imp(obj[, group_2[6:10]], ncp = 2, seed = 1234)

  expected_results <- obj
  expected_results[, group_1[1:5]] <- sub1
  expected_results[, group_2[6:10]] <- sub2

  expect_identical(grouped_results[, ], expected_results)
})

test_that("group-specific parameters work correctly in parallel, pca", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("carrier")
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature
  # Different ncp and coeff.ridge values for each group
  group_df <- data.frame(
    feature = I(list(group_1[1:3], group_2[1:4])),
    aux = I(list(group_1, group_2)),
    parameters = I(list(
      list(ncp = 2, coeff.ridge = 1),
      list(ncp = 3, coeff.ridge = 2)
    ))
  )
  obj <- to_test$input
  mirai::daemons(2, seed = 1234)
  grouped_results <- group_imp(obj, group = group_df, cores = 1, seed = 1234, nb.init = 10)
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

# prep_groups ----
test_that("prep_groups returns correct structure without k/ncp", {
  obj <- matrix(rnorm(2 * 5), nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e")))
  features_df <- data.frame(
    feature = c("a", "b", "c", "d"),
    group = c("g1", "g1", "g2", "g2")
  )

  result <- prep_groups(colnames(obj), features_df, allow_unmapped = TRUE)
  expect_true(inherits(result, "slideimp_tbl"))
  expect_true("group" %in% names(result))
  expect_true("feature" %in% names(result))
  expect_equal(sort(result$group), c("g1", "g2"))
})

test_that("prep_groups handles subset correctly", {
  obj <- matrix(rnorm(2 * 6), nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e", "f")))
  features_df <- data.frame(
    feature = c("a", "b", "c", "d", "e", "f"),
    group = c("g1", "g1", "g1", "g2", "g2", "g2")
  )

  result <- prep_groups(
    colnames(obj), features_df,
    subset = c("a", "b", "d", "e")
  )

  expect_true("aux" %in% names(result))

  g1_row <- result[result$group == "g1", ]
  g2_row <- result[result$group == "g2", ]

  # features in subset go to features column
  expect_setequal(g1_row$feature[[1]], c("a", "b"))
  expect_setequal(g2_row$feature[[1]], c("d", "e"))

  # features not in subset go to aux column
  expect_equal(g1_row$aux[[1]], "c")
  expect_equal(g2_row$aux[[1]], "f")
})

test_that("prep_groups errors when no subset element found", {
  obj <- matrix(rnorm(2 * 3), nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b", "c")))
  features_df <- data.frame(
    feature = c("a", "b", "c"),
    group = c("g1", "g1", "g2")
  )

  expect_error(
    suppressWarnings(prep_groups(colnames(obj), features_df, subset = c("x", "y", "z"))),
    "x, y, z"
  )
})

test_that("prep_groups pads groups to min_group_size", {
  n <- 2
  p <- 6
  obj <- matrix(rnorm(n * p), nrow = n)
  rownames(obj) <- paste0("r", seq_len(n))
  colnames(obj) <- letters[1:p]

  features_df <- data.frame(
    feature = c("a", "b", "c"),
    group = c("g1", "g1", "g2")
  )

  result <- prep_groups(
    colnames(obj),
    features_df,
    min_group_size = 4,
    seed = 123,
    allow_unmapped = TRUE
  )
  expect_true("aux" %in% names(result))
  g1_row <- result[result$group == "g1", ]
  g2_row <- result[result$group == "g2", ]
  # g1 has 2 features, needs 2 more to reach min_group_size of 4
  expect_equal(length(g1_row$feature[[1]]) + length(g1_row$aux[[1]]), 4)
  # g2 has 1 feature, needs 3 more
  expect_equal(length(g2_row$feature[[1]]) + length(g2_row$aux[[1]]), 4)
})

test_that("prep_groups errors when min_group_size too large", {
  n <- 2
  p <- 3
  obj <- matrix(rnorm(n * p), nrow = n)
  rownames(obj) <- paste0("r", seq_len(n))
  colnames(obj) <- letters[1:p]

  features_df <- data.frame(
    feature = c("a", "b", "c"),
    group = c("g1", "g1", "g2")
  )

  expect_error(
    prep_groups(colnames(obj), features_df, min_group_size = 100),
    "too large"
  )
})

test_that("prep_groups errors when no colnames match features_df", {
  n <- 2
  p <- 3
  obj <- matrix(rnorm(n * p), nrow = n)
  rownames(obj) <- paste0("r", seq_len(n))
  colnames(obj) <- letters[1:p]

  features_df <- data.frame(
    feature = c("x", "y", "z"),
    group = c("g1", "g1", "g2")
  )

  expect_error(
    prep_groups(colnames(obj), features_df),
    "No groups remain after pruning"
  )
})

# slideimp.extra ----
test_that("slideimp_extra_manifests works with prep_groups", {
  skip_if_not_installed("slideimp.extra")
  skip_on_cran()
  slideimp.extra::set_slideimp_path("dev")
  msa <- slideimp.extra::ilmn_manifest("MSA", deduped = TRUE, rawdir = "dev")
  n_feat <- length(msa$feature)
  sim_mat <- matrix(rnorm(1 * n_feat), nrow = 1, dimnames = list(NULL, msa$feature))
  expect_no_error(prep_groups(colnames(sim_mat), group = msa))
  expect_no_error(prep_groups(colnames(sim_mat), group = "MSA_deduped"))
  slideimp.extra::set_slideimp_path(NULL)
})
