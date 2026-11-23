skip_on_cran()

# infeasible ----
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
  # leave row 20 with values so each column has a mean
  mat[20, 1:5] <- c(0.1, 0.2, 0.3, 0.4, 0.5)

  bad_grp <- data.frame(
    feature = colnames(mat)[1:5],
    group = "bad_group"
  )

  res <- suppressMessages(
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

  expect_true("bad_group" %in% attr(res, "fallback"))
  expect_equal(attr(res, "fallback_action"), "mean")
  # no remaining NA in the group's subset columns (means were computable)
  expect_false(anyNA(res[, 1:5]))
})

test_that("group_imp: mixed feasible + infeasible groups - only bad group falls back", {
  set.seed(1234)
  sim <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)
  mat <- sim$input

  # bad_group (features 1:5) is made infeasible below; good_group (features
  # 6:15) is left as it is
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
        suppressMessages(do.call(group_imp, args))
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

  # Regression: before groups passed their features to pca_imp() as
  # `subset`, the whole submatrix of every group was fitted. Group 1 needed
  # nothing and still aborted the run under the default
  # `on_infeasible = "error"`. It is still not sent to the backend, which
  # group_imp() reaches as pca_imp_parts(), at all.
  real_pca_imp_parts <- pca_imp_parts
  seen_ncol <- integer()
  local_mocked_bindings(
    pca_imp_parts = function(obj, ...) {
      seen_ncol <<- c(seen_ncol, ncol(obj))
      real_pca_imp_parts(obj, ...)
    }
  )
  res <- suppressMessages(group_imp(obj, group = g, ncp = 2, .progress = FALSE))
  expect_identical(seen_ncol, 6L) # group 2 only
  expect_length(attr(res, "fallback"), 0L)
  expect_identical(res[, 1:10], obj[, 1:10])
  expect_false(anyNA(res[, 11:16]))
})

test_that("a group whose feature holes are all over colmax is infeasible", {
  # Each group passes its features to the backend as `subset`, and both
  # backends refuse a `subset` whose holes are all in columns over colmax,
  # even when the group's aux columns have holes they could fit. So such a
  # group goes through on_infeasible under PCA as under K-NN. Before groups
  # passed `subset` to pca_imp(), it fitted the whole group and mean-filled
  # the feature.
  set.seed(3)
  obj <- matrix(
    rnorm(30 * 8),
    30,
    8,
    dimnames = list(NULL, paste0("f", 1:8))
  )
  obj[1:28, 1] <- NA # a feature over colmax = 0.9
  obj[1:3, 5:8] <- NA # aux columns with holes the backends could fit
  g <- data.frame(
    group = "g",
    feature = I(list(c("f1", "f2"))),
    aux = I(list(paste0("f", 3:8)))
  )
  # non-vacuity: the feature holes are all in f1, and the aux columns are
  # eligible with holes
  expect_false(anyNA(obj[, 2]))
  expect_true(all(colMeans(is.na(obj[, 5:8])) <= 0.9))
  expect_true(anyNA(obj[, 5:8]))

  for (backend in list(list(k = 2), list(ncp = 2))) {
    info <- names(backend)
    expect_error(
      suppressMessages(do.call(
        group_imp,
        c(list(obj, group = g, .progress = FALSE), backend)
      )),
      class = "slideimp_infeasible"
    )
    res <- suppressMessages(do.call(
      group_imp,
      c(
        list(obj, group = g, on_infeasible = "mean", .progress = FALSE),
        backend
      )
    ))
    expect_identical(attr(res, "fallback"), "g", info = info)
    expect_worst_err(unname(res[1:28, 1]), rep(mean(obj[29:30, 1]), 28L))
    # the aux columns come back as they are
    expect_identical(res[, 3:8], obj[, 3:8], info = info)
  }
})

test_that("group_imp's mean fallback is handed only the features with a hole", {
  # Regression, moved from mean_imp_col(): the fallback part over an
  # infeasible group took its complete features too, so fallback_part()
  # copied them and their is.na() to fill nothing. A complete column in a
  # filling part still validates, so only the columns it is handed show it.
  x <- matrix(
    rnorm(20 * 10),
    20,
    10,
    dimnames = list(NULL, paste0("c", 1:10))
  )
  x[1:19, 1:9] <- NA # over colmax: the group is infeasible
  expect_false(anyNA(x[, 10]))
  g <- data.frame(group = "a", feature = I(list(colnames(x))))
  real_fallback_part <- fallback_part
  calls <- list()
  local_mocked_bindings(
    fallback_part = function(obj, cols, method = c("mean", "median")) {
      calls <<- c(calls, list(cols))
      real_fallback_part(obj, cols, method)
    }
  )
  res <- suppressMessages(group_imp(
    x,
    group = g,
    k = 3,
    on_infeasible = "mean",
    .progress = FALSE
  ))
  expect_identical(as.character(attr(res, "fallback")), "a")
  expect_identical(calls, list(1:9))
  expect_false(anyNA(res[,]))
  expect_identical(unname(res[, 10]), unname(x[, 10]))
})
