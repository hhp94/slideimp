skip_on_cran()

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

  # the PCA fits run a fixed iteration count, miniter = maxiter, so both
  # calls take the same steps; the threshold, far above the criterion at that
  # count, decides only `converged`. Two runs through BLAS, so the values
  # agree to rounding.
  run_pca <- function(group) {
    group_imp(
      obj,
      group = group,
      ncp = 5,
      nb.init = 10,
      seed = 1234,
      threshold = 0.1,
      maxiter = 30,
      miniter = 30
    )
  }
  pca_long <- run_pca(group_long)
  pca_list <- run_pca(group_list)
  expect_identical(dimnames(pca_long), dimnames(pca_list))
  expect_worst_err(pca_long[,], pca_list[,])
  same <- setdiff(names(attributes(pca_long)), "fallback")
  expect_identical(attributes(pca_long)[same], attributes(pca_list)[same])
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

  # manual imputation for comparison. knn_imp() returns only the `subset`
  # columns, which are the features; every other column comes back as it is.
  sub1_cols <- unique(c(group_df$feature[[1]], group_df$aux[[1]]))
  sub1 <- knn_imp(obj[, sub1_cols], k = 3, subset = group_df$feature[[1]])
  sub2_cols <- unique(c(group_df$feature[[2]], group_df$aux[[2]]))
  sub2 <- knn_imp(obj[, sub2_cols], k = 3, subset = group_df$feature[[2]])
  expect_identical(colnames(sub1), group_df$feature[[1]])
  expect_identical(colnames(sub2), group_df$feature[[2]])

  expected_results <- obj
  expected_results[, colnames(sub1)] <- sub1
  expected_results[, colnames(sub2)] <- sub2
  # Compare results
  expect_identical(dimnames(grouped_results), dimnames(expected_results))
  expect_worst_err(grouped_results[,], expected_results)
})

test_that("group_imp places each group by column, whatever its feature order", {
  # Each group's parts come back in its slice's numbering, features first,
  # and group_imp() moves them to obj's (rebase_parts()). Here the features
  # are named out of obj's order and aux columns sit between them, so the
  # move has to sort; the result is each group's own backend result written
  # back by column, in one process.
  #
  # Every group reads obj as it was passed. Group "b"'s aux columns are
  # features of group "a" with holes, and "b" sees those holes, not "a"'s
  # imputed values, in either running mode. Until the groups returned parts,
  # the sequential path wrote each group into obj before the next one read
  # its slice, so "b" saw "a"'s values there and the mirai path did not.
  set.seed(11)
  obj <- sim_mat(40, 24, perc_total_na = 0.2, perc_col_na = 1)$input
  cn <- colnames(obj)
  g <- data.frame(
    group = c("a", "b"),
    feature = I(list(
      cn[c(9, 2, 5, 12, 1, 7)],
      cn[c(20, 14, 23, 17, 13, 24, 15, 18, 16, 19, 21, 22)]
    )),
    aux = I(list(cn[c(3, 4, 6, 8, 10, 11)], cn[c(1, 2)]))
  )
  for (f in g$feature) {
    expect_true(is.unsorted(match(f, cn)))
    expect_gt(sum(is.na(obj[, f])), 0L)
  }
  expect_true(all(g$aux[[2]] %in% g$feature[[1]]))
  expect_gt(sum(is.na(obj[, g$aux[[2]]])), 0L)

  by_hand <- function(imp) {
    out <- obj
    for (i in seq_len(nrow(g))) {
      f <- g$feature[[i]]
      sub <- obj[, c(f, g$aux[[i]]), drop = FALSE]
      out[, f] <- suppressMessages(imp(sub, seq_along(f)))
    }
    out
  }

  knn <- suppressMessages(group_imp(obj, group = g, k = 3, .progress = FALSE))
  knn_hand <- by_hand(function(sub, s) knn_imp(sub, k = 3, subset = s))
  expect_identical(knn[,], knn_hand)

  # the PCA fits run a fixed iteration count, miniter = maxiter, so both
  # paths take the same steps; the threshold, far above the criterion at that
  # count, decides only `converged`. They are separate fits through BLAS, so
  # the values agree to rounding.
  fixed <- list(threshold = 1e-6, maxiter = 30, miniter = 30)
  pca <- suppressMessages(do.call(
    group_imp,
    c(list(obj, group = g, ncp = 2, seed = 1, .progress = FALSE), fixed)
  ))
  pca_hand <- by_hand(function(sub, s) {
    do.call(pca_imp, c(list(sub, ncp = 2, seed = 1, subset = s), fixed))
  })
  expect_identical(dimnames(pca), dimnames(pca_hand))
  expect_worst_err(pca[,], pca_hand)
})

test_that("group_imp with subset returns only the subset columns, in obj's order", {
  # prep_groups() keeps a group's `subset` features as its features, in the
  # group's order, and makes its other features aux columns after them, so
  # each group's slice is those columns in that order. The result holds the
  # `subset` columns only, in obj's order, whatever order `subset` names
  # them in, and the record counts its columns.
  set.seed(5)
  obj <- sim_mat(40, 24, perc_total_na = 0.2, perc_col_na = 1)$input
  cn <- colnames(obj)
  obj[1:38, 9] <- NA # over colmax: a mean part, in `mean_filled_cols`
  obj[5, 1:12] <- NA # K-NN finds no neighbor observed at row 5 in group a
  g <- data.frame(
    group = c("a", "b"),
    feature = I(list(cn[1:12], cn[13:24]))
  )
  sel <- cn[c(20, 9, 3, 14)]
  ord <- c(3L, 9L, 14L, 20L)
  expect_identical(sort(match(sel, cn)), ord)

  by_hand <- function(imp) {
    out <- obj
    for (f in g$feature) {
      kept <- f[f %in% sel]
      sub <- obj[, c(kept, setdiff(f, kept)), drop = FALSE]
      out[, kept] <- suppressMessages(imp(sub, seq_along(kept)))
    }
    out[, ord]
  }
  run <- function(...) {
    suppressMessages(group_imp(
      obj,
      group = g,
      subset = sel,
      .progress = FALSE,
      ...
    ))
  }

  knn <- run(k = 3)
  expect_identical(colnames(knn), cn[ord])
  expect_identical(
    knn[,],
    by_hand(function(sub, s) knn_imp(sub, k = 3, subset = s))
  )
  # column 9 is the result's 2nd, and K-NN left row 5 of column 3, the 1st
  expect_identical(attr(knn, "mean_filled_cols"), 2L)
  expect_identical(attr(knn, "mean_filled_cells"), cbind(row = 5L, col = 1L))

  # the PCA fits run a fixed iteration count, as in the test above, and
  # agree to rounding
  fixed <- list(threshold = 1e-3, maxiter = 30, miniter = 30)
  pca <- do.call(run, c(list(ncp = 2, seed = 1), fixed))
  pca_hand <- by_hand(function(sub, s) {
    do.call(pca_imp, c(list(sub, ncp = 2, seed = 1, subset = s), fixed))
  })
  expect_identical(dimnames(pca), dimnames(pca_hand))
  expect_worst_err(pca[,], pca_hand)
  expect_identical(attr(pca, "mean_filled_cols"), 2L)

  # a `subset` whose columns have no hole is returned as it is, though other
  # columns have holes
  seq_mat <- matrix(
    as.double(seq_len(40L * 24L)),
    40L,
    24L,
    dimnames = dimnames(obj)
  )
  seq_mat[cbind(1:10, 1:10)] <- NA
  complete <- cn[c(22, 15)]
  expect_false(anyNA(seq_mat[, complete]))
  out <- suppressMessages(
    group_imp(seq_mat, group = g, subset = complete, k = 3, .progress = FALSE)
  )
  expect_identical(out[,], seq_mat[, c(15L, 22L)])
  expect_identical(attr(out, "mean_filled_cols"), integer(0))
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

  # Manual verification with different parameters. knn_imp() returns only
  # the `subset` columns; every other column comes back as it is.
  sub1 <- knn_imp(obj[, group_1], k = 3, subset = group_1[1:3], dist_pow = 0)
  sub2 <- knn_imp(obj[, group_2], k = 7, subset = group_2[1:4], dist_pow = 1)
  expected_results <- obj
  expected_results[, colnames(sub1)] <- sub1
  expected_results[, colnames(sub2)] <- sub2

  expect_identical(dimnames(grouped_results), dimnames(expected_results))
  expect_worst_err(grouped_results[,], expected_results)
})

test_that("grouped imputation works without aux columns, knn", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # no aux columns, only feature. group_2[6:10] needs sim_mat() to have drawn
  # at least 10 features into group2
  expect_gte(length(group_2), 10L)
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

  expect_identical(dimnames(grouped_results), dimnames(expected_results))
  expect_worst_err(grouped_results[,], expected_results)
})

test_that("group-specific parameters work correctly, pca", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # Different ncp and coeff.ridge values for each group. Every fit runs a
  # fixed iteration count, miniter = maxiter, so both paths take the same
  # steps; the threshold, far above the criterion at that count, decides only
  # `converged`.
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
  grouped_results <- group_imp(
    obj,
    group = group_df,
    threshold = 1e-4,
    maxiter = 30,
    miniter = 30
  )

  # Manual verification with different parameters
  fit <- function(x, ...) {
    pca_imp(x, seed = 1234, threshold = 1e-4, maxiter = 30, miniter = 30, ...)
  }
  sub1_full <- fit(obj[, group_1], ncp = 2, coeff.ridge = 1)
  sub1 <- obj[, group_1]
  sub1[, group_1[1:3]] <- sub1_full[, group_1[1:3]]
  sub2_full <- fit(obj[, group_2], ncp = 3, coeff.ridge = 2)
  sub2 <- obj[, group_2]
  sub2[, group_2[1:4]] <- sub2_full[, group_2[1:4]]
  expected_results <- cbind(sub1, sub2)[, colnames(obj)]
  expect_identical(dimnames(grouped_results), dimnames(expected_results))
  expect_worst_err(grouped_results[,], expected_results)
})

test_that("grouped imputation works without aux columns, pca", {
  set.seed(1234)
  to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
  group_1 <- subset(to_test$col_group, group == "group1")$feature
  group_2 <- subset(to_test$col_group, group == "group2")$feature

  # no aux columns, only feature. group_2[6:10] needs sim_mat() to have drawn
  # at least 10 features into group2
  expect_gte(length(group_2), 10L)
  group_df <- data.frame(
    feature = I(list(group_1[1:5], group_2[6:10]))
  )
  obj <- to_test$input

  # every fit runs a fixed iteration count, as in the test above
  fixed <- list(threshold = 1e-4, maxiter = 30, miniter = 30)
  grouped_results <- do.call(
    group_imp,
    c(
      list(obj, group = group_df, ncp = 2, seed = 1234, allow_unmapped = TRUE),
      fixed
    )
  )

  # Build expected results: start with original and update only imputed columns
  fit <- function(x) do.call(pca_imp, c(list(x, ncp = 2, seed = 1234), fixed))
  sub1 <- fit(obj[, group_1[1:5]])
  sub2 <- fit(obj[, group_2[6:10]])

  expected_results <- obj
  expected_results[, group_1[1:5]] <- sub1
  expected_results[, group_2[6:10]] <- sub2

  expect_identical(dimnames(grouped_results), dimnames(expected_results))
  expect_worst_err(grouped_results[,], expected_results)
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
    "scale.* does not apply to K-NN"
  )
  expect_error(
    group_imp(obj, group = grp, k = 3, scale = TRUE, maxiter = 10),
    "scale.*maxiter.* do not apply to K-NN"
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
    "1 unknown parameter in .group\\$parameters. for K-NN"
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

# Two groups, "a" and "b", each with a hole, for the tests of group_imp()'s
# refusals and running-mode messages below.
msg_case <- function() {
  set.seed(1)
  obj <- matrix(rnorm(20 * 6), 20, 6, dimnames = list(NULL, paste0("c", 1:6)))
  obj[1:2, 1] <- NA
  obj[3, 4] <- NA
  list(
    obj = obj,
    group = data.frame(
      group = c("a", "b"),
      feature = I(list(paste0("c", 1:3), paste0("c", 4:6)))
    )
  )
}

test_that("group_imp refuses k and ncp passed together", {
  case <- msg_case()
  expect_error(
    group_imp(case$obj, group = case$group, k = 2, ncp = 2, .progress = FALSE),
    "^Supply at most one of .k. \\(K-NN imputation\\) or .ncp. \\(PCA imputation\\)\\.",
    class = "rlang_error"
  )
})

test_that("group_imp refuses a group that has both k and ncp", {
  case <- msg_case()
  g <- case$group
  g$parameters <- I(list(list(ncp = 2), list()))
  expect_error(
    group_imp(case$obj, group = g, k = 2, .progress = FALSE),
    "^1 group has both .k. and .ncp.:\n. a\n",
    class = "rlang_error"
  )
})

test_that("group_imp refuses groups that have neither k nor ncp", {
  case <- msg_case()
  expect_error(
    group_imp(case$obj, group = case$group, .progress = FALSE),
    "^2 groups have neither .k. nor .ncp.:\n. a, b\n",
    class = "rlang_error"
  )
})

test_that("group_imp refuses groups that mix k and ncp", {
  case <- msg_case()
  g <- case$group
  g$parameters <- I(list(list(k = 2), list(ncp = 2)))
  expect_error(
    group_imp(case$obj, group = g, .progress = FALSE),
    "^Some groups use .k. and others use .ncp.\\.",
    class = "rlang_error"
  )
})

test_that("group_imp refuses a per-group method of the other backend", {
  case <- msg_case()
  g <- case$group
  g$parameters <- I(list(list(method = "EM"), list()))
  expect_error(
    group_imp(case$obj, group = g, k = 2, .progress = FALSE),
    paste0(
      "^1 group has an invalid .method. for K-NN imputation:\n. a\n",
      ".*Set .method. to \"euclidean\" or \"manhattan\"\\."
    ),
    class = "rlang_error"
  )
  g$parameters <- I(list(list(method = "euclidean"), list()))
  expect_error(
    group_imp(case$obj, group = g, ncp = 2, .progress = FALSE),
    paste0(
      "^1 group has an invalid .method. for PCA imputation:\n. a\n",
      ".*Set .method. to \"regularized\" or \"EM\"\\."
    ),
    class = "rlang_error"
  )
})

test_that("group_imp says cores is ignored for PCA without daemons", {
  case <- msg_case()
  msgs <- rlang_messages(
    group_imp(
      case$obj,
      group = case$group,
      ncp = 1,
      cores = 2,
      .progress = FALSE
    )
  )
  expect_match(
    msgs,
    "^.cores. = 2 is ignored for PCA imputation\\.\n. Use .mirai::daemons\\(\\). to impute the groups in parallel\\.$",
    all = FALSE
  )
  expect_true("Running mode: sequential" %in% msgs)
})

test_that("group_imp says cores is ignored for K-NN under mirai daemons", {
  case <- msg_case()
  full <- case$obj
  full[is.na(full)] <- 0
  # daemons are reported set, but no group needs imputing, so no task is
  # sent to one
  local_mocked_bindings(mirai_daemons_set = function() TRUE)
  msgs <- rlang_messages(
    group_imp(full, group = case$group, k = 2, cores = 2, .progress = FALSE)
  )
  expect_match(
    msgs,
    "^.cores. = 2 is ignored, to avoid nested parallelism with the mirai daemons\\.\n. To use 2 threads in place of the daemons, call .mirai::daemons\\(0\\). first\\.$",
    all = FALSE
  )
  expect_true("Running mode: mirai" %in% msgs)
})

test_that("group_imp reports the threaded running mode for K-NN", {
  case <- msg_case()
  msgs <- rlang_messages(
    group_imp(case$obj, group = case$group, k = 2, cores = 2, .progress = FALSE)
  )
  expect_true("Running mode: threaded (2 cores)" %in% msgs)
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
  expect_match(printed[[1]], "^Method: K-NN imputation$")
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
  expect_worst_err(unname(knn[2:20, 9]), rep(case$obj[[1, 9]], 19L))

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

test_that("a daemon missing a slideimp internal gets the build hint", {
  # the daemons' results are mocked: what is pinned is the host's reading of
  # a failed task, with and without a missing function in its message
  x <- matrix(rnorm(20 * 6), 20, 6, dimnames = list(NULL, paste0("c", 1:6)))
  x[1, 1] <- NA
  g <- data.frame(feature = colnames(x), group = 1)
  fail_with <- function(msg) {
    local_mocked_bindings(
      mirai_daemons_set = function() TRUE,
      collect_mirai_map = function(mapped, .progress) {
        list(structure(list(), class = c("miraiError", "errorValue")))
      },
      mirai_error_text = function(x) msg
    )
    local_mocked_bindings(
      mirai_map = function(...) NULL,
      .package = "mirai"
    )
    tryCatch(
      suppressMessages(group_imp(x, group = g, k = 2, .progress = FALSE)),
      error = conditionMessage
    )
  }
  missing <- fail_with("could not find function \"mean_filled_record\"")
  expect_match(missing, "Imputation failed in 1 group", fixed = TRUE)
  expect_match(missing, "another build of slideimp", fixed = TRUE)
  other <- fail_with("subscript out of bounds")
  expect_match(other, "Imputation failed in 1 group", fixed = TRUE)
  expect_no_match(other, "another build", fixed = TRUE)
})

test_that("group_imp keeps named dimnames, as knn_imp does", {
  set.seed(11)
  x <- matrix(
    rnorm(20 * 12),
    20,
    12,
    dimnames = list(sample = paste0("s", 1:20), probe = paste0("c", 1:12))
  )
  x[sample(length(x), 20)] <- NA
  attr(x, "extra") <- "dropped"
  g <- data.frame(feature = colnames(x), group = rep(1:2, each = 6))
  # with holes, with `subset`, and with neither: the fill, the cut and the
  # matrix as passed reach the result by three paths
  complete <- x
  complete[is.na(complete)] <- 0
  sub <- colnames(x)[c(2, 9)]
  for (case in list(
    list(obj = x, subset = NULL),
    list(obj = x, subset = sub),
    list(obj = complete, subset = NULL),
    list(obj = complete, subset = sub)
  )) {
    res <- suppressMessages(group_imp(
      case$obj,
      group = g,
      k = 3,
      subset = case$subset,
      .progress = FALSE
    ))
    keep <- if (is.null(case$subset)) colnames(x) else case$subset
    expect_identical(dimnames(res), dimnames(case$obj[, keep, drop = FALSE]))
    expect_identical(names(dimnames(res)), c("sample", "probe"))
    expect_null(attr(res, "extra"))
    ref <- knn_imp(case$obj, k = 3, subset = case$subset, .progress = FALSE)
    expect_identical(dimnames(res), dimnames(ref))
  }
})
