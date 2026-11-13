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

  expect_equal(
    group_imp(obj, group = group_long, k = 3),
    group_imp(obj, group = group_list, k = 3),
    ignore_attr = "fallback"
  )

  expect_equal(
    group_imp(obj, group = group_long, ncp = 5, nb.init = 10, seed = 1234),
    group_imp(obj, group = group_list, ncp = 5, nb.init = 10, seed = 1234),
    ignore_attr = "fallback"
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

  expect_equal(
    group_imp(obj, group = group_doubled, k = 3),
    group_imp(obj, group = group_single, k = 3)
  )
})

test_that("group_imp() handles aux columns present in only some groups (padded)", {
  set.seed(1234)
  to_test <- sim_mat(
    n = 20,
    p = 50,
    n_col_groups = 2,
    perc_total_na = 0.3,
    perc_col_na = 1
  )
  obj <- to_test$input
  meta <- to_test$col_group

  # move feature 1 into its own tiny group so padding is needed for that
  # group only - the other groups will have zero aux columns.
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
  expect_equal(dim(res), dim(obj))
  expect_equal(colnames(res), colnames(obj))
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
  expect_equal(grouped_results[,], expected_results)
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
  sub1 <- pca_imp(obj[, sub1_cols], ncp = 2, seed = 1234)[, group_df$feature[[
    1
  ]]]
  sub2_cols <- unique(c(group_df$feature[[2]], group_df$aux[[2]]))
  sub2 <- pca_imp(obj[, sub2_cols], ncp = 2, seed = 1234)[, group_df$feature[[
    2
  ]]]

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

  expect_equal(grouped_results[,], expected_results)
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

  expect_equal(grouped_results[,], expected_results)
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
  expect_equal(grouped_results[,], expected_results)
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

  expect_equal(grouped_results[,], expected_results)
})

test_that("group-specific parameters work correctly in parallel, pca", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature
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
  withr::defer(mirai::daemons(0))
  expect_message(
    grouped_results <- group_imp(
      obj,
      group = group_df,
      cores = 1,
      seed = 1234,
      nb.init = 10
    ),
    regexp = "Running mode.*mirai"
  )
  sub1_full <- pca_imp(
    obj[, group_1],
    ncp = 2,
    coeff.ridge = 1,
    seed = 1234,
    nb.init = 10
  )
  sub1 <- obj[, group_1]
  sub1[, group_1[1:3]] <- sub1_full[, group_1[1:3]]
  sub2_full <- pca_imp(
    obj[, group_2],
    ncp = 3,
    coeff.ridge = 2,
    seed = 1234,
    nb.init = 10
  )
  sub2 <- obj[, group_2]
  sub2[, group_2[1:4]] <- sub2_full[, group_2[1:4]]
  expected_results <- cbind(sub1, sub2)[, colnames(obj)]
  imputed_cols <- c(group_1[1:3], group_2[1:4])
  obj_orig <- obj[, imputed_cols]
  grouped_values <- grouped_results[, imputed_cols][is.na(obj_orig)]
  expected_values <- expected_results[, imputed_cols][is.na(obj_orig)]
  expect_true(cor(grouped_values, expected_values) > 0.999)
})

test_that("group_imp aborts when a group fails inside a mirai daemon", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1234)
  to_test <- sim_mat(30, 12, perc_total_na = 0.2, perc_col_na = 1)
  obj <- to_test$input
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # `coeff.ridge = -1` passes group_imp()'s own validation untouched and is
  # rejected by pca_imp() as an ordinary error, not a `slideimp_infeasible`
  # condition, so the worker's tryCatch never sees it. Regression: the mirai
  # path used to collect that failure silently - the flags were coerced to
  # character, `big_out` was never written for the group, and group 2 came
  # back all NA with observed values destroyed and no condition raised.
  group_df <- data.frame(
    feature = I(list(group_1, group_2)),
    aux = I(list(character(0), character(0))),
    parameters = I(list(
      list(ncp = 2, coeff.ridge = 1),
      list(ncp = 2, coeff.ridge = -1)
    ))
  )

  # the sequential path lets the error propagate; the mirai path must match
  expect_error(
    suppressMessages(group_imp(obj, group = group_df, seed = 1234)),
    regexp = "coeff.ridge"
  )

  mirai::daemons(2, seed = 1234)
  withr::defer(mirai::daemons(0))
  expect_error(
    suppressMessages(group_imp(obj, group = group_df, seed = 1234)),
    regexp = "Imputation failed in group 2"
  )
})

test_that("group_imp under mirai matches the sequential run group by group", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(5)
  obj <- matrix(
    rnorm(40 * 24),
    40,
    24,
    dimnames = list(NULL, paste0("f", 1:24))
  )
  obj[sample(40, 6), 1:8] <- NA # group "a": imputed
  obj[1:3, 17:23] <- NA # group "c": imputed
  obj[5, 24] <- NA # group "d": one column, too small
  grp <- data.frame(
    feature = colnames(obj),
    group = rep(c("a", "b", "c", "d"), c(8, 8, 7, 1))
  )
  # group "b" is complete, so it never reaches a daemon. The tasks that do
  # run then write to output columns that are not their own row positions.
  expect_false(anyNA(obj[, 9:16]))
  run <- function() {
    suppressMessages(group_imp(
      obj,
      group = grp,
      k = 3,
      on_infeasible = "mean",
      .progress = FALSE
    ))
  }
  seq_res <- run()
  expect_identical(attr(seq_res, "fallback"), "d")

  mirai::daemons(2, seed = 1234)
  withr::defer(mirai::daemons(0))
  par_res <- run()
  expect_identical(attr(par_res, "fallback"), "d")
  err_abs <- abs(unclass(par_res) - unclass(seq_res))
  nz <- seq_res != 0
  max_abs <- max(err_abs)
  max_rel <- max(err_abs[nz] / abs(seq_res[nz]))
  expect_lte(
    max_abs,
    1e-12,
    label = sprintf(
      "max abs error %.3g (max rel %.3g over nonzero cells)",
      max_abs,
      max_rel
    )
  )
  expect_lte(max_rel, 1e-12, label = sprintf("max rel error %.3g", max_rel))

  # an infeasible group keeps its class on the way back from a daemon. Group
  # "d" is made complete, or the host would refuse it before any daemon ran.
  bad <- obj
  bad[5, 24] <- 0
  bad[1:39, 1:8] <- NA
  got <- tryCatch(
    {
      suppressMessages(group_imp(bad, group = grp, k = 3, .progress = FALSE))
      "no error"
    },
    slideimp_infeasible = function(e) conditionMessage(e),
    error = function(e) "some other error"
  )
  expect_match(got, "Group a is infeasible")
})

# prep_groups ----
test_that("prep_groups returns correct structure without k/ncp", {
  obj <- matrix(
    rnorm(2 * 5),
    nrow = 2,
    dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e"))
  )
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
  obj <- matrix(
    rnorm(2 * 6),
    nrow = 2,
    dimnames = list(c("r1", "r2"), c("a", "b", "c", "d", "e", "f"))
  )
  features_df <- data.frame(
    feature = c("a", "b", "c", "d", "e", "f"),
    group = c("g1", "g1", "g1", "g2", "g2", "g2")
  )

  result <- prep_groups(
    colnames(obj),
    features_df,
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
  obj <- matrix(
    rnorm(2 * 3),
    nrow = 2,
    dimnames = list(c("r1", "r2"), c("a", "b", "c"))
  )
  features_df <- data.frame(
    feature = c("a", "b", "c"),
    group = c("g1", "g1", "g2")
  )

  expect_error(
    suppressWarnings(prep_groups(
      colnames(obj),
      features_df,
      subset = c("x", "y", "z")
    )),
    "x, y, z"
  )
})

test_that("prep_groups pads groups to min_group_size", {
  n <- 2
  p <- 6
  obj <- matrix(rnorm(n * p), nrow = n)
  rownames(obj) <- paste0("r", seq_len(n))
  colnames(obj) <- letters[1:p]

  # every column is mapped: padding draws from the columns `group` maps
  features_df <- data.frame(
    feature = c("a", "b", "c", "d", "e", "f"),
    group = c("g1", "g1", "g2", "g3", "g3", "g3")
  )

  result <- prep_groups(
    colnames(obj),
    features_df,
    min_group_size = 4,
    seed = 123
  )
  expect_true("aux" %in% names(result))
  g1_row <- result[result$group == "g1", ]
  g2_row <- result[result$group == "g2", ]
  # g1 has 2 features, needs 2 more to reach min_group_size of 4
  expect_equal(length(g1_row$feature[[1]]) + length(g1_row$aux[[1]]), 4)
  # g2 has 1 feature, needs 3 more
  expect_equal(length(g2_row$feature[[1]]) + length(g2_row$aux[[1]]), 4)
})

test_that("prep_groups pads only from columns that group maps", {
  cn <- paste0("c", 1:20)
  g <- data.frame(feature = cn[1:8], group = rep(c("g1", "g2"), c(3, 5)))
  unmapped <- cn[9:20]

  # Regression: the padding pool was every column of `obj`, so under
  # `allow_unmapped = TRUE` the columns documented as "not used as auxiliary
  # columns" were drawn as aux. With 12 unmapped columns against 8 mapped
  # ones, most draws came from them.
  res <- suppressMessages(prep_groups(
    cn,
    g,
    allow_unmapped = TRUE,
    min_group_size = 6,
    seed = 1
  ))
  aux <- unlist(res$aux)
  expect_length(aux, 4L) # g1 needs 3 more columns, g2 needs 1
  expect_false(any(aux %in% unmapped))
  expect_identical(lengths(res$feature) + lengths(res$aux), c(6L, 6L))

  # g1 has 3 columns and only the 5 of g2 to draw from, so 9 is out of reach.
  # The unmapped columns do not make up the difference.
  expect_error(
    suppressMessages(prep_groups(
      cn,
      g,
      allow_unmapped = TRUE,
      min_group_size = 9
    )),
    "too large"
  )
})

test_that("prep_groups demotes non-subset features and drops emptied groups", {
  cn <- letters[1:8]
  g <- data.frame(
    group = c("g1", "g2", "g3"),
    feature = I(list(c("a", "b", "c"), c("d", "e"), c("f", "g"))),
    aux = I(list("h", character(0), "a"))
  )

  # g1 keeps a and c and demotes b behind its own aux. g2 has no feature in
  # `subset` and is dropped, named by its label. g3 keeps f and demotes g.
  expect_message(
    res <- prep_groups(cn, g, subset = c("c", "a", "f")),
    "Dropped group: g2"
  )
  expect_identical(res$group, c("g1", "g3"))
  expect_identical(res$feature, list(c("a", "c"), "f"))
  expect_identical(res$aux, list(c("h", "b"), c("a", "g")))
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
    "columns in `obj` have no matching entry"
  )
})

# infeasible ----
test_that("group_imp: on_infeasible = 'error' rethrows slideimp_infeasible", {
  set.seed(1234)
  sim <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)
  mat <- sim$input
  grp <- sim$col_group

  mat[1:19, "feature6"] <- NA
  mat[20, ] <- rnorm(20)

  bad_grp <- data.frame(
    feature = colnames(mat)[1:5],
    group = "bad_group",
    stringsAsFactors = FALSE
  )
  mat_bad <- sim$input
  mat_bad[1:19, 1:5] <- NA # all subset cols fully NA

  expect_error(
    group_imp(
      mat_bad,
      group = bad_grp,
      k = 3,
      colmax = 0.9,
      on_infeasible = "error",
      .progress = FALSE,
      allow_unmapped = TRUE
    ),
    class = "slideimp_infeasible"
  )
})

test_that("group_imp: on_infeasible = 'skip' retains original NAs and flags group", {
  set.seed(1234)
  sim <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)
  mat <- sim$input
  mat[1:19, 1:5] <- NA # subset cols of the infeasible group all NA
  mat[20, ] <- rnorm(20)

  bad_grp <- data.frame(
    feature = colnames(mat)[1:5],
    group = "bad_group"
  )

  res <- suppressMessages(
    suppressWarnings(
      group_imp(
        mat,
        group = bad_grp,
        k = 3,
        colmax = 0.9,
        on_infeasible = "skip",
        .progress = FALSE,
        allow_unmapped = TRUE
      )
    )
  )

  # columns should still be fully NA (skip = original retained)
  expect_true(all(is.na(res[1:19, 1:5])))
  # group should be recorded in fallback attribute
  expect_true("bad_group" %in% attr(res, "fallback"))
  expect_equal(attr(res, "fallback_action"), "skip")
})

test_that("group_imp: on_infeasible = 'mean' fills with column means and flags group", {
  set.seed(1234)
  sim <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)
  mat <- sim$input

  mat[1:19, 1:5] <- NA # 95% missing, exceeds default colmax = 0.9
  # leave row 20 with values so mean_imp_col produces something
  mat[20, 1:5] <- c(0.1, 0.2, 0.3, 0.4, 0.5)

  bad_grp <- data.frame(
    feature = colnames(mat)[1:5],
    group = "bad_group"
  )

  res <- suppressMessages(
    suppressWarnings(
      group_imp(
        mat,
        group = bad_grp,
        k = 3,
        colmax = 0.9,
        on_infeasible = "mean",
        allow_unmapped = TRUE,
        .progress = FALSE
      )
    )
  )

  expect_true("bad_group" %in% attr(res, "fallback"))
  expect_equal(attr(res, "fallback_action"), "mean")
  # no remaining NA in the group's subset columns (means were computable)
  expect_false(anyNA(res[, 1:5]))
})

test_that("group_imp: mixed feasible + infeasible groups - only bad group falls back", {
  set.seed(1234)
  sim <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)
  mat <- sim$input

  # force group2 (features 1:5) to be infeasible; leave group1 (feature6) alone
  # But group1 has only 1 feature -> k would cap to 0. Expand group1 artificially.
  good_grp <- data.frame(
    feature = colnames(mat)[6:15],
    group = "good_group"
  )
  bad_grp <- data.frame(
    feature = colnames(mat)[1:5],
    group = "bad_group"
  )
  grp <- rbind(good_grp, bad_grp)

  mat[1:19, 1:5] <- NA
  mat[20, ] <- rnorm(20)

  res <- suppressMessages(
    suppressWarnings(
      group_imp(
        mat,
        group = grp,
        k = 3,
        colmax = 0.9,
        on_infeasible = "skip",
        .progress = FALSE,
        allow_unmapped = TRUE
      )
    )
  )

  expect_setequal(attr(res, "fallback"), "bad_group")
  # good_group columns should be fully imputed
  expect_false(anyNA(res[, 6:15]))
  # bad_group columns should remain NA
  expect_true(all(is.na(res[1:19, 1:5])))
})

test_that("group_imp: on_infeasible = 'error' is the default and propagates", {
  set.seed(1234)
  sim <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)
  mat <- sim$input
  mat[1:19, 1:5] <- NA
  mat[20, ] <- rnorm(20)

  bad_grp <- data.frame(
    feature = colnames(mat)[1:5],
    group = "bad_group",
    stringsAsFactors = FALSE
  )

  # no on_infeasible specified -> default "error"
  expect_error(
    suppressMessages(
      group_imp(
        mat,
        group = bad_grp,
        k = 3,
        colmax = 0.9,
        allow_unmapped = TRUE,
        .progress = FALSE
      )
    ),
    class = "slideimp_infeasible"
  )
})

test_that("knn_imp: slideimp_infeasible class is set on both abort sites", {
  # Site 1: k > n_elig - 1
  mat <- matrix(NA_real_, nrow = 20, ncol = 20)
  diag(mat) <- rnorm(20)
  expect_error(knn_imp(mat, k = 2), class = "slideimp_infeasible")

  # Site 2: all subset cols exceed colmax
  set.seed(1)
  mat2 <- matrix(rnorm(400), 20, 20)
  mat2[1:19, 1:5] <- NA
  mat2[20, ] <- rnorm(20)
  expect_error(
    knn_imp(mat2, k = 3, colmax = 0.5, subset = 1:5),
    class = "slideimp_infeasible"
  )
})

test_that("group_imp: the infeasible error reaches a slideimp_infeasible handler", {
  set.seed(1234)
  sim <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)
  mat_bad <- sim$input
  mat_bad[1:19, 1:5] <- NA
  bad_grp <- data.frame(feature = colnames(mat_bad)[1:5], group = "bad_group")

  # Regression: the rethrow carried the class only on its parent.
  # `expect_error(class = )` searches the parent chain and passed all along;
  # a caller's `tryCatch(slideimp_infeasible = )` matches the condition itself
  # and did not catch it.
  for (param in c("k", "ncp")) {
    args <- list(
      mat_bad,
      group = bad_grp,
      on_infeasible = "error",
      .progress = FALSE,
      allow_unmapped = TRUE
    )
    args[[param]] <- 3
    got <- tryCatch(
      {
        suppressMessages(suppressWarnings(do.call(group_imp, args)))
        "no error"
      },
      slideimp_infeasible = function(e) "caught",
      error = function(e) "some other error"
    )
    expect_identical(got, "caught", info = param)
  }
})

test_that("group_imp routes a one-column group through on_infeasible", {
  set.seed(1234)
  sim <- sim_mat(30, 20, perc_total_na = 0.1, perc_col_na = 1)
  obj <- sim$input
  grp <- data.frame(
    feature = colnames(obj)[1:11],
    group = c(rep("a", 10), "lonely")
  )
  # the lone column has something to impute, or the group would be left alone
  expect_true(anyNA(obj[, 11]))
  run <- function(...) {
    suppressMessages(group_imp(
      obj,
      group = grp,
      allow_unmapped = TRUE,
      .progress = FALSE,
      ...
    ))
  }

  # Regression: the group was capped to `k = 0` with a warning and then
  # aborted the run as an input error, whatever `on_infeasible` said.
  err <- tryCatch(run(k = 3), slideimp_infeasible = function(e) e)
  expect_s3_class(err, "slideimp_infeasible")
  expect_match(conditionMessage(err), "Group lonely is infeasible")

  expect_no_warning(res <- run(k = 3, on_infeasible = "mean"))
  expect_identical(attr(res, "fallback"), "lonely")
  expect_false(anyNA(res[, 1:10]))
  ref <- obj[, 11]
  ref[is.na(ref)] <- mean(ref, na.rm = TRUE)
  err_abs <- abs(res[, 11] - ref)
  nz <- ref != 0
  expect_gt(sum(nz), 0L)
  max_abs <- max(err_abs)
  max_rel <- max(err_abs[nz] / abs(ref[nz]))
  expect_lte(
    max_abs,
    1e-12,
    label = sprintf(
      "max abs error %.3g (max rel %.3g over nonzero cells)",
      max_abs,
      max_rel
    )
  )
  expect_lte(max_rel, 1e-12, label = sprintf("max rel error %.3g", max_rel))

  expect_no_warning(res <- run(k = 3, on_infeasible = "skip"))
  expect_identical(attr(res, "fallback"), "lonely")
  expect_identical(res[, 11], obj[, 11])
  expect_false(anyNA(res[, 1:10]))

  expect_no_warning(res <- run(ncp = 2, on_infeasible = "skip"))
  expect_identical(attr(res, "fallback"), "lonely")
  expect_identical(res[, 11], obj[, 11])

  # a one-column group with nothing to impute is not infeasible
  obj_full <- obj
  obj_full[, 11] <- seq_len(30) / 30
  res <- suppressMessages(group_imp(
    obj_full,
    group = grp,
    k = 3,
    allow_unmapped = TRUE,
    .progress = FALSE
  ))
  expect_length(attr(res, "fallback"), 0L)
  expect_identical(res[, 11], obj_full[, 11])
})

test_that("group_imp does not impute a group whose features are complete", {
  set.seed(2)
  obj <- matrix(
    rnorm(30 * 16),
    30,
    16,
    dimnames = list(NULL, paste0("f", 1:16))
  )
  obj[1:29, 6] <- NA # an aux column of group 1, above colmax = 0.9
  obj[1:4, 11:13] <- NA # group 2 has features to impute
  g <- data.frame(
    group = c("complete", "missing"),
    feature = I(list(colnames(obj)[1:5], colnames(obj)[11:16])),
    aux = I(list(colnames(obj)[6:10], character(0)))
  )
  expect_false(anyNA(obj[, 1:5]))
  # pca_imp() itself refuses group 1's columns, so calling it aborts the run
  expect_error(
    suppressMessages(pca_imp(obj[, 1:10], ncp = 2)),
    class = "slideimp_infeasible"
  )

  # Regression: pca_imp() has no `subset`, so the whole submatrix of every
  # group was handed to it. Group 1 needed nothing and still aborted the run
  # under the default `on_infeasible = "error"`.
  real_pca_imp <- pca_imp
  seen_ncol <- integer()
  local_mocked_bindings(
    pca_imp = function(obj, ...) {
      seen_ncol <<- c(seen_ncol, ncol(obj))
      real_pca_imp(obj, ...)
    }
  )
  res <- suppressMessages(group_imp(obj, group = g, ncp = 2, .progress = FALSE))
  expect_identical(seen_ncol, 6L) # group 2 only
  expect_length(attr(res, "fallback"), 0L)
  expect_identical(res[, 1:10], obj[, 1:10])
  expect_false(anyNA(res[, 11:16]))
})

# parameters ----
test_that("group_imp: a global seed works with K-NN and seeds the padding", {
  set.seed(1234)
  sim <- sim_mat(30, 20, perc_total_na = 0.1, perc_col_na = 1)
  obj <- sim$input
  grp <- sim$col_group
  # one group is below 12 columns, so the padding has columns to draw
  expect_true(any(table(grp$group) < 12))

  # Regression: `seed` was copied into every group's parameters and K-NN
  # rejected it as an unknown parameter in `group$parameters`, so the one
  # argument that makes the padding reproducible could not be used with `k`.
  run <- function(seed) {
    suppressMessages(group_imp(
      obj,
      group = grp,
      k = 3,
      seed = seed,
      min_group_size = 12,
      .progress = FALSE
    ))
  }
  first <- run(1)
  expect_identical(first, run(1))

  # it is the padding prep_groups() draws with that seed
  pp1 <- prep_groups(colnames(obj), grp, min_group_size = 12, seed = 1)
  expect_identical(
    first,
    suppressMessages(group_imp(obj, group = pp1, k = 3, .progress = FALSE))
  )

  # and the seed is used, not ignored: another one pads with other columns
  pp2 <- prep_groups(colnames(obj), grp, min_group_size = 12, seed = 2)
  expect_false(identical(pp1$aux, pp2$aux))
  expect_false(identical(first, run(2)))
})

test_that("group_imp names a global argument of the other backend", {
  set.seed(1234)
  sim <- sim_mat(30, 20, perc_total_na = 0.1, perc_col_na = 1)
  obj <- sim$input
  grp <- sim$col_group

  # Regression: these were reported as unknown names in `group$parameters`,
  # which the caller never set.
  expect_error(
    group_imp(obj, group = grp, k = 3, scale = TRUE),
    "scale.* does not apply to KNN"
  )
  expect_error(
    group_imp(obj, group = grp, k = 3, scale = TRUE, maxiter = 10),
    "scale.*maxiter.* do not apply to KNN"
  )
  expect_error(
    group_imp(obj, group = grp, ncp = 2, dist_pow = 1),
    "dist_pow.* does not apply to PCA"
  )

  # a name a group sets for itself is still reported against the table
  pp <- prep_groups(colnames(obj), grp)
  pp$parameters <- rep(list(list(k = 3, seed = 1)), nrow(pp))
  expect_error(
    group_imp(obj, group = pp),
    "Unknown parameter in .group\\$parameters. for KNN"
  )
})

test_that("group_imp validates a per-group k before it compares it", {
  set.seed(3)
  sim <- sim_mat(30, 40, perc_total_na = 0.1, perc_col_na = 1)
  obj <- sim$input
  pp <- prep_groups(
    colnames(obj),
    data.frame(feature = colnames(obj), group = rep(1:2, each = 20))
  )
  with_k <- function(k) {
    pp$parameters <- rep(list(list(k = k)), nrow(pp))
    pp
  }

  # Regression: the cap compared the value first. A character "3" is greater
  # than an integer 19, so both groups ran with `k = 19` and a "capped from 3
  # to 19" warning.
  expect_no_warning(expect_error(
    group_imp(obj, group = with_k("3")),
    "group$parameters[[1]]$k",
    fixed = TRUE
  ))
  expect_error(
    group_imp(obj, group = with_k(NA)),
    "group$parameters[[1]]$k",
    fixed = TRUE
  )
  expect_error(
    group_imp(obj, group = with_k(c(2, 3))),
    "group$parameters[[1]]$k",
    fixed = TRUE
  )
  # the global argument is named as itself
  expect_error(
    group_imp(obj, group = pp, k = "3"),
    "Assertion on 'k'",
    fixed = TRUE
  )
})

test_that("group_imp validates its scalar arguments", {
  set.seed(1234)
  sim <- sim_mat(30, 20, perc_total_na = 0.1, perc_col_na = 1)
  obj <- sim$input
  grp <- sim$col_group

  expect_error(group_imp(obj, group = grp, k = 3, cores = NULL), "'cores'")
  expect_error(
    group_imp(obj, group = grp, k = 3, .progress = NA),
    "'.progress'",
    fixed = TRUE
  )
  expect_error(
    group_imp(obj, group = grp, k = 3, allow_unmapped = NA),
    "'allow_unmapped'"
  )
})

test_that("group_imp names a group by its label, not its row", {
  set.seed(1234)
  sim <- sim_mat(30, 20, perc_total_na = 0.1, perc_col_na = 1)
  obj <- sim$input
  # labels 20 and 10: prep_groups() sorts them, so label 10 is row 1
  grp <- data.frame(
    feature = colnames(obj)[1:12],
    group = c(rep(20, 8), rep(10, 4))
  )
  expect_warning(
    suppressMessages(group_imp(
      obj,
      group = grp,
      k = 5,
      allow_unmapped = TRUE,
      .progress = FALSE
    )),
    "Group 10: .k. capped from 5 to 3"
  )
})

test_that("a group_imp() result can be imputed again", {
  set.seed(1234)
  sim <- sim_mat(30, 20, perc_total_na = 0.1, perc_col_na = 1)
  first <- suppressMessages(group_imp(
    sim$input,
    group = sim$col_group,
    k = 3,
    .progress = FALSE
  ))
  expect_identical(attr(first, "metacaller"), "group_imp")

  # Regression: the second pass attached the class a second time and kept
  # `metacaller` and `fallback_action` from the first, so print() reported a
  # group_imp() run with "1 group used a fallback" and "See group: FALSE".
  again <- suppressMessages(knn_imp(first, k = 3))
  expect_identical(class(again), c("slideimp_results", "matrix", "array"))
  expect_null(attr(again, "metacaller"))
  expect_null(attr(again, "fallback_action"))
  printed <- capture.output(print(again))
  expect_match(printed[[1]], "^Method: KNN imputation$")
  expect_false(any(grepl("fallback", printed)))

  slid <- suppressMessages(slide_imp(
    first,
    location = seq_len(ncol(first)),
    window_size = 10,
    overlap_size = 5,
    min_window_n = 5,
    k = 3,
    .progress = FALSE
  ))
  expect_identical(class(slid), c("slideimp_results", "matrix", "array"))
  expect_identical(attr(slid, "metacaller"), "slide_imp")
})

# One matrix for the post_imp tests below. Group "b" is columns 7:12, so its
# column over colmax, f9, is local column 3 and column 9 of `obj`.
post_imp_case <- function() {
  set.seed(1234)
  obj <- matrix(
    rnorm(20 * 12),
    20,
    12,
    dimnames = list(NULL, paste0("f", 1:12))
  )
  obj[1:3, 2] <- NA # group "a": ordinary missingness
  obj[2:20, 9] <- NA # group "b": miss_rate 0.95, above the default colmax
  obj[4:5, 11] <- NA # group "b": ordinary missingness
  list(
    obj = obj,
    group = data.frame(
      feature = colnames(obj),
      group = rep(c("a", "b"), c(6, 6))
    )
  )
}

test_that("group_imp records post_imp as the groups applied it", {
  case <- post_imp_case()
  run <- function(group = case$group, ...) {
    suppressMessages(group_imp(
      case$obj,
      group = group,
      k = 3,
      .progress = FALSE,
      ...
    ))
  }

  # Regression: the default `post_imp = NULL` hands each group the backend
  # default, TRUE, and the attribute recorded the NULL that was passed.
  res <- run()
  expect_false(anyNA(res[, 9]))
  expect_identical(attr(res, "post_imp"), TRUE)
  expect_identical(attr(run(post_imp = TRUE), "post_imp"), TRUE)

  res_off <- run(post_imp = FALSE)
  expect_true(anyNA(res_off[, 9]))
  expect_identical(attr(res_off, "post_imp"), FALSE)

  # set group by group, and not the same in both: no one value describes it
  mixed <- data.frame(
    feature = I(split(case$group$feature, case$group$group)),
    parameters = I(list(list(post_imp = FALSE), list()))
  )
  res_mixed <- run(group = mixed)
  expect_identical(attr(res_mixed, "post_imp"), NA)
  # group "b" kept the default, so its column over colmax is filled
  expect_false(anyNA(res_mixed[, 9]))

  # no group has anything to impute, so no value was applied: the attribute
  # is the value passed, and NA when none was
  full <- case$obj
  full[is.na(full)] <- 0
  run_full <- function(...) {
    suppressMessages(group_imp(
      full,
      group = case$group,
      k = 3,
      .progress = FALSE,
      ...
    ))
  }
  expect_identical(attr(run_full(), "post_imp"), NA)
  expect_identical(attr(run_full(post_imp = FALSE), "post_imp"), FALSE)
})

test_that("group_imp records mean-filled columns by their index in obj", {
  case <- post_imp_case()
  no_cells <- cbind(row = integer(0), col = integer(0))
  run <- function(...) {
    suppressMessages(group_imp(
      case$obj,
      group = case$group,
      .progress = FALSE,
      ...
    ))
  }

  knn <- run(k = 3)
  expect_false(anyNA(knn[,]))
  expect_identical(attr(knn, "mean_filled_cols"), 9L)
  expect_identical(attr(knn, "mean_filled_cells"), no_cells)
  expect_equal(unname(knn[2:20, 9]), rep(case$obj[[1, 9]], 19L))

  knn_off <- run(k = 3, post_imp = FALSE)
  expect_true(anyNA(knn_off[, 9]))
  expect_identical(attr(knn_off, "mean_filled_cols"), integer(0))
  expect_identical(attr(knn_off, "mean_filled_cells"), no_cells)

  pca <- run(ncp = 2)
  expect_identical(attr(pca, "mean_filled_cols"), 9L)
  expect_identical(attr(pca, "mean_filled_cells"), no_cells)

  printed <- capture.output(print(knn))
  expect_true(any(grepl("1 column was excluded from imputation", printed)))
})

test_that("group_imp under mirai records the same mean-filled columns", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  case <- post_imp_case()
  run <- function() {
    suppressMessages(group_imp(
      case$obj,
      group = case$group,
      k = 3,
      .progress = FALSE
    ))
  }
  seq_res <- run()
  expect_identical(attr(seq_res, "mean_filled_cols"), 9L)

  mirai::daemons(2, seed = 1234)
  withr::defer(mirai::daemons(0))
  par_res <- run()
  expect_identical(attr(par_res, "mean_filled_cols"), 9L)
  expect_identical(
    attr(par_res, "mean_filled_cells"),
    attr(seq_res, "mean_filled_cells")
  )
  expect_identical(attr(par_res, "post_imp"), attr(seq_res, "post_imp"))
})
