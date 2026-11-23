skip_on_cran()

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

# --- Preconditioning tests for the `group` parameter ---
test_that("group must be a data.frame", {
  obj <- matrix(
    as.double(1:12),
    nrow = 3,
    dimnames = list(NULL, paste0("f", 1:4))
  )
  # a list has no slideimp_resolve_group() method, so the default one refuses
  # it before any check on its columns
  expect_error(
    group_imp(obj, group = list(feature = list("f1")), k = 2),
    "No .slideimp_resolve_group\\(\\). method for objects of class .list.",
    class = "rlang_error"
  )
})

test_that("group must contain a 'feature' column", {
  obj <- matrix(
    as.double(1:12),
    nrow = 3,
    dimnames = list(NULL, paste0("f", 1:4))
  )
  group_df <- data.frame(stuff = I(list(c("f1", "f2"))))
  expect_error(group_imp(obj, group = group_df, k = 2), "feature")
})

test_that("group errors when feature is character without group column", {
  obj <- matrix(
    as.double(1:12),
    nrow = 3,
    dimnames = list(NULL, paste0("f", 1:4))
  )
  group_df <- data.frame(feature = c("f1", "f2"))
  expect_error(
    group_imp(obj, group = group_df, k = 2),
    "no group column"
  )
})

test_that("group errors on NA in group column", {
  obj <- matrix(
    as.double(1:12),
    nrow = 3,
    dimnames = list(NULL, paste0("f", 1:4))
  )
  group_df <- data.frame(
    group = c("a", NA_character_),
    feature = c("f1", "f2")
  )
  expect_error(group_imp(obj, group = group_df, k = 2), "NA")
})

test_that("group feature column must not have NAs", {
  obj <- matrix(
    as.double(1:12),
    nrow = 3,
    dimnames = list(NULL, paste0("f", 1:4))
  )
  group_df <- data.frame(
    group = c("a", "a"),
    feature = c("f1", NA_character_)
  )
  expect_error(group_imp(obj, group = group_df, k = 2), "feature")
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
    "in more than one group",
    fixed = TRUE
  )
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
    prep_groups(
      colnames(obj),
      features_df,
      subset = c("x", "y", "z")
    ),
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
    "columns have no matching entry"
  )
})

test_that("prep_groups refuses a group table with no feature among the columns", {
  # under allow_unmapped = TRUE the unmapped columns are let through, and the
  # group with no feature left is dropped, so no group remains
  expect_error(
    suppressMessages(prep_groups(
      c("a", "b"),
      data.frame(feature = c("x", "y"), group = 1),
      allow_unmapped = TRUE
    )),
    "^No group has a feature that matches a column name\\.\n. Check that .group\\$feature. holds the column names of the data\\.$",
    class = "rlang_error"
  )
})

test_that("prep_groups refuses a subset feature that no group has", {
  # "c" is a column, and an aux column of the group, but no group's feature
  g <- data.frame(group = 1, feature = I(list(c("a", "b"))), aux = I(list("c")))
  expect_error(
    prep_groups(c("a", "b", "c"), g, subset = "c"),
    "^1 feature in .subset. is not assigned to any group:\n. c\n. Add it to a group, or remove it from .subset.\\.$",
    class = "rlang_error"
  )
})

test_that("prep_groups warns on non-integer numeric group labels", {
  g <- data.frame(feature = c("a", "b", "c"), group = c(1.5, 1.5, 2))
  expect_warning(
    res <- prep_groups(c("a", "b", "c"), g),
    "^.group\\$group. has non-integer numeric values\\.\n. Each distinct value forms one group\\.\n. Use integer or character labels\\.$",
    class = "rlang_warning"
  )
  # as the warning says, each distinct value forms one group
  expect_identical(res$group, c(1.5, 2))
  expect_identical(res$feature, list(c("a", "b"), "c"))
})
