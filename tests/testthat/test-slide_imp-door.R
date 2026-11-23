skip_on_cran()

test_that("slide_imp rejects ncp above the row cap pca_imp enforces", {
  set.seed(1234)
  # 10 rows, so pca_imp() allows at most nrow - 2 = 8 components. The column
  # half of its cap is n_elig - 1, well above 8 for a 50-column window.
  to_test <- sim_mat(10, 300, perc_total_na = 0.2, perc_col_na = 1)$input
  location <- seq_len(ncol(to_test))

  run <- function(ncp) {
    slide_imp(
      to_test,
      location = location,
      window_size = 50,
      overlap_size = 0,
      min_window_n = 10,
      ncp = ncp,
      miniter = 2,
      colmax = 1,
      seed = 1234,
      .progress = FALSE
    )
  }

  # Regression: slide_imp() bounded `ncp` at min(nrow, ncol) - 1 = 9, one above
  # pca_imp()'s row cap. `ncp = 9` passed validation and then aborted EVERY
  # window with `slideimp_infeasible`, which the default on_infeasible = "skip"
  # turned into the input returned unimputed - no error, no warning, every NA
  # still there.
  expect_error(run(9), "ncp")

  # the value one below is genuinely feasible and imputes, so the bound is not
  # simply rejecting everything.
  res <- suppressMessages(run(8))
  expect_true(anyNA(to_test))
  expect_false(anyNA(res[,]))
})

test_that("slide_imp scans the union of covered columns, not window by window", {
  set.seed(1234)
  mat <- sim_mat(20, 30, perc_total_na = 0.2)$input
  location <- 1:30

  args <- list(
    location = location,
    k = 3,
    window_size = 10,
    overlap_size = 5,
    min_window_n = 5,
    .progress = FALSE
  )

  # The windows are [1, 10], [6, 15], [11, 20], ... One scan of `obj` is
  # taken before any window runs, so a dry run, which runs none, refuses
  # what a full run does. A scan window by window would let the dry run
  # through.
  for (dry_run in c(FALSE, TRUE)) {
    run_args <- c(args, list(dry_run = dry_run))

    # an Inf anywhere under a window is refused
    inf_mat <- mat
    inf_mat[3, 12] <- Inf
    expect_error(
      suppressMessages(do.call(slide_imp, c(list(inf_mat), run_args))),
      "holds infinite values or columns with no observed value",
      fixed = TRUE
    )

    # so is an all-NA column, including one sitting in an overlap region
    na_mat <- mat
    na_mat[, 8] <- NA
    expect_error(
      suppressMessages(do.call(slide_imp, c(list(na_mat), run_args))),
      "holds infinite values or columns with no observed value",
      fixed = TRUE
    )
  }
})

test_that("slide_imp scopes the all-NA check to covered columns but Inf to all", {
  set.seed(1234)
  mat <- sim_mat(20, 60, perc_total_na = 0.2)$input
  location <- 1:60

  # subset = 25:35 keeps two windows covering columns 21:40, so 1:20 and 41:60
  # are covered by no window.
  args <- list(
    location = location,
    k = 3,
    window_size = 10,
    overlap_size = 0,
    min_window_n = 5,
    subset = 25:35,
    .progress = FALSE
  )

  # an all-NA column outside every window is never imputed, so it is allowed
  uncovered_na <- mat
  uncovered_na[, 5] <- NA
  expect_no_error(
    suppressMessages(do.call(slide_imp, c(list(uncovered_na), args)))
  )

  # the same column inside a window is refused
  covered_na <- mat
  covered_na[, 25] <- NA
  expect_error(
    suppressMessages(do.call(slide_imp, c(list(covered_na), args))),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )

  # Inf is refused wherever it sits, under a window or not
  uncovered_inf <- mat
  uncovered_inf[3, 5] <- Inf
  expect_error(
    suppressMessages(do.call(slide_imp, c(list(uncovered_inf), args))),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )

  covered_inf <- mat
  covered_inf[4, 25] <- Inf
  expect_error(
    suppressMessages(do.call(slide_imp, c(list(covered_inf), args))),
    "holds infinite values or columns with no observed value",
    fixed = TRUE
  )
})

test_that("slide_imp refuses a subset it cannot return in full", {
  set.seed(1234)
  x <- sim_mat(20, 30, perc_total_na = 0.2)$input
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  run <- function(...) {
    slide_imp(
      x,
      location = seq_len(30),
      window_size = 10,
      min_window_n = 5,
      k = 3,
      .progress = FALSE,
      ...
    )
  }

  # The result holds the `subset` columns, so a name it would drop is an
  # error, where it was skipped with a message, and so is an empty `subset`,
  # where the input came back unimputed. Both refusals name slide_imp(), and
  # come before the windows, so a dry run gives the same one.
  for (dry_run in c(FALSE, TRUE)) {
    err <- expect_error(
      run(subset = c("c2", "not_a_column"), dry_run = dry_run),
      "not in `obj`: not_a_column",
      fixed = TRUE
    )
    expect_match(deparse(conditionCall(err))[[1]], "^slide_imp\\(")
    expect_error(
      run(subset = character(0), dry_run = dry_run),
      "selects no columns",
      fixed = TRUE
    )
    expect_error(
      run(subset = integer(0), flank = TRUE, dry_run = dry_run),
      "selects no columns",
      fixed = TRUE
    )
  }
  expect_error(run(subset = c("c2", "c2")), "more than once", fixed = TRUE)
  expect_error(run(subset = 31L), "subset", fixed = TRUE)
  unnamed <- unname(x)
  expect_error(
    slide_imp(
      unnamed,
      location = seq_len(30),
      window_size = 10,
      min_window_n = 5,
      k = 3,
      subset = "c2",
      .progress = FALSE
    ),
    "has no column names",
    fixed = TRUE
  )
})

test_that("slide_imp ignores overlap_size when flank = TRUE", {
  set.seed(1234)
  x <- sim_mat(20, 30, perc_total_na = 0.2)$input
  args <- list(
    x,
    location = seq_len(30),
    window_size = 5,
    overlap_size = 8,
    min_window_n = 2,
    dry_run = TRUE,
    .progress = FALSE
  )

  # Regression: the comparison against `window_size` ran before the flank
  # branch discarded `overlap_size`, so a documented no-op argument errored.
  expect_no_error(do.call(slide_imp, c(args, list(flank = TRUE, subset = 10L))))
  # it is still enforced where it means something
  expect_error(do.call(slide_imp, args), "strictly less than")
})
