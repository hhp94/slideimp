skip_if_not_manual()

test_that("group-specific parameters work correctly in parallel, pca", {
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
  # every fit runs a fixed iteration count, miniter = maxiter, so the daemons
  # and the fits by hand take the same steps; the threshold, far above the
  # criterion at that count, decides only `converged`
  mirai::daemons(2, seed = 1234)
  withr::defer(mirai::daemons(0))
  expect_message(
    grouped_results <- group_imp(
      obj,
      group = group_df,
      cores = 1,
      seed = 1234,
      nb.init = 10,
      threshold = 1e-3,
      maxiter = 30,
      miniter = 30
    ),
    regexp = "Running mode.*mirai"
  )
  # `daemons(seed = )` runs the daemons on L'Ecuyer-CMRG, and `seed` reseeds
  # within that kind, so the nb.init starts on a daemon are those of the
  # same seed under L'Ecuyer-CMRG here, not under the session's own kind
  by_hand <- function() {
    sub1_full <- pca_imp(
      obj[, group_1],
      ncp = 2,
      coeff.ridge = 1,
      seed = 1234,
      nb.init = 10,
      threshold = 1e-3,
      maxiter = 30,
      miniter = 30
    )
    sub1 <- obj[, group_1]
    sub1[, group_1[1:3]] <- sub1_full[, group_1[1:3]]
    sub2_full <- pca_imp(
      obj[, group_2],
      ncp = 3,
      coeff.ridge = 2,
      seed = 1234,
      nb.init = 10,
      threshold = 1e-3,
      maxiter = 30,
      miniter = 30
    )
    sub2 <- obj[, group_2]
    sub2[, group_2[1:4]] <- sub2_full[, group_2[1:4]]
    cbind(sub1, sub2)[, colnames(obj)]
  }
  expected_results <- withr::with_seed(
    1234,
    by_hand(),
    .rng_kind = "L'Ecuyer-CMRG"
  )
  expect_identical(dimnames(grouped_results), dimnames(expected_results))
  expect_worst_err(grouped_results[,], expected_results)
})

test_that("group_imp aborts when a group fails inside a mirai daemon", {
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
    regexp = "Imputation failed in 1 group: 2",
    fixed = TRUE
  )
})

test_that("group_imp under mirai matches the sequential run group by group", {
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
  expect_gt(sum(nz), 0L)
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

test_that("group_imp under mirai records the same mean-filled columns", {
  skip_if_not_installed("withr")
  case <- post_imp_case()
  run <- function(...) {
    suppressMessages(group_imp(
      case$obj,
      group = case$group,
      k = 3,
      .progress = FALSE,
      ...
    ))
  }
  seq_res <- run()
  expect_identical(attr(seq_res, "mean_filled_cols"), 9L)
  # with `subset`, column 9 is the result's 2nd: the result holds columns 4,
  # 9 and 11 of `obj`
  sel <- colnames(case$obj)[c(11, 4, 9)]
  seq_sub <- run(subset = sel)
  expect_identical(attr(seq_sub, "mean_filled_cols"), 2L)

  mirai::daemons(2, seed = 1234)
  withr::defer(mirai::daemons(0))
  par_res <- run()
  expect_identical(attr(par_res, "mean_filled_cols"), 9L)
  expect_identical(
    attr(par_res, "mean_filled_cells"),
    attr(seq_res, "mean_filled_cells")
  )
  expect_identical(attr(par_res, "post_imp"), attr(seq_res, "post_imp"))
  par_sub <- run(subset = sel)
  expect_identical(par_sub, seq_sub)
})
