skip_on_cran()

test_that("slide_imp: on_infeasible = 'mean' fills skipped windows with column means", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
  # Make cols 1:10 infeasible (0.95 > colmax) but not fully NA, so a column
  # mean exists. colmax is inclusive: 18/20 = 0.9 would still be eligible.
  mat[1:19, 1:10] <- NA
  mat[20, 1:10] <- rnorm(10)

  expect_warning(
    res <- suppressMessages(slide_imp(
      mat,
      location = location,
      k = 3,
      window_size = 10,
      overlap_size = 0,
      min_window_n = 5,
      colmax = 0.9,
      on_infeasible = "mean",
      .progress = FALSE
    )),
    "filled with column means",
    class = "slideimp_fallback"
  )

  expect_gt(length(attr(res, "fallback")), 0)
  expect_equal(attr(res, "fallback_action"), "mean")
  # No remaining NA because column means were available
  expect_false(anyNA(res[, 1:10]))
})

test_that("slide_imp: mixed feasible + infeasible windows - skip isolates correctly", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
  mat[1:19, 11:20] <- NA
  mat[20, 11:20] <- rnorm(10)

  expect_warning(
    res <- suppressMessages(slide_imp(
      mat,
      location = location,
      k = 3,
      window_size = 10,
      overlap_size = 0,
      min_window_n = 5,
      colmax = 0.9,
      on_infeasible = "skip",
      .progress = FALSE
    )),
    class = "slideimp_fallback"
  )

  expect_true(all(is.na(res[1:19, 11:20])))
  expect_false(anyNA(res[, 1:10]))
  expect_false(anyNA(res[, 21:30]))
  # the skipped window should be recorded
  expect_gt(length(attr(res, "fallback")), 0)
  expect_equal(attr(res, "fallback_action"), "skip")
})

test_that("slide_imp: flank mode - infeasible flank window skips only its target", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
  # kill a region around target 15 so its flanking window is infeasible
  mat[1:19, 10:20] <- NA
  mat[20, 10:20] <- rnorm(11)

  expect_warning(
    res <- suppressMessages(slide_imp(
      mat,
      location = location,
      k = 2,
      window_size = 10,
      flank = TRUE,
      subset = c(15, 60),
      min_window_n = 5,
      colmax = 0.9,
      on_infeasible = "skip",
      .progress = FALSE
    )),
    class = "slideimp_fallback"
  )

  # the result holds the two targets, columns 15 and 60 of `mat`
  expect_identical(dim(res), c(nrow(mat), 2L))
  # target 15's window was skipped -> col 15 retains original NA
  expect_true(anyNA(mat[, 15]))
  expect_identical(res[, 1], mat[, 15])
  # target 60's window was feasible -> col 60 is imputed
  expect_true(anyNA(mat[, 60]))
  expect_false(anyNA(res[, 2]))
})

test_that("slide_imp: overlapping windows - a skipped window is not counted in the average", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
  # Windows are [1, 20], [11, 30], [21, 40], ... Columns 1:18 are 95% missing,
  # so window [1, 20] has 2 usable columns against k = 3 and is infeasible,
  # while [11, 30] still has 12 and is not.
  mat[1:19, 1:18] <- NA
  mat[20, 1:18] <- rnorm(18)

  expect_warning(
    res <- suppressMessages(slide_imp(
      mat,
      location = location,
      k = 3,
      window_size = 20,
      overlap_size = 10,
      min_window_n = 10,
      colmax = 0.9,
      on_infeasible = "skip",
      .progress = FALSE
    )),
    class = "slideimp_fallback"
  )

  # the fixture really skips window 1 and nothing else. An earlier version of
  # this test left 10 usable columns in window 1, so no window was skipped and
  # the skip path was never reached.
  expect_identical(attr(res, "fallback"), 1L)

  # columns 11:20 sit under [1, 20] (skipped) and [11, 30] (feasible). Only
  # windows that ran are counted, so their count is 1 and they must equal
  # window [11, 30] alone, with no division. Counting the skipped window too
  # halves every imputed value here.
  alone <- suppressMessages(knn_imp(
    mat[, 11:30],
    k = 3,
    colmax = 0.9,
    post_imp = TRUE,
    .progress = FALSE
  ))
  a <- res[, 11:20]
  b <- alone[, 1:10]
  expect_false(anyNA(b))
  expect_identical(is.na(a), is.na(b))
  expect_equal(max(abs(a - b)), 0)
  # relative error over the cells of columns 11:20 whose reference is not
  # exactly 0, which sim_mat() can produce
  nz <- b != 0
  expect_gt(sum(nz), 100L)
  expect_equal(max(abs(a[nz] - b[nz]) / abs(b[nz])), 0)

  # columns 1:10 sit under the skipped window only, so they come back unchanged
  expect_identical(res[, 1:10], mat[, 1:10])
})

test_that("slide_imp: all windows infeasible under 'skip' returns original matrix", {
  set.seed(1234)
  mat <- matrix(NA_real_, nrow = 20, ncol = 100)
  colnames(mat) <- sample.int(100, size = 100)
  for (j in seq_len(ncol(mat))) {
    mat[1, j] <- rnorm(1)
  }
  location <- 1:100

  # one warning for the whole run, not one per window, and it is raised with
  # `.progress = FALSE`. Regression: this run used to return the input
  # unimputed with no condition of any kind.
  warns <- list()
  res <- withCallingHandlers(
    suppressMessages(slide_imp(
      mat,
      location = location,
      k = 3,
      window_size = 10,
      overlap_size = 0,
      min_window_n = 5,
      colmax = 0.5,
      on_infeasible = "skip",
      .progress = FALSE
    )),
    warning = function(w) {
      warns[[length(warns) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warns, 1L)
  expect_s3_class(warns[[1L]], "slideimp_fallback")
  expect_match(conditionMessage(warns[[1L]]), "10 of 10 windows are infeasible")
  expect_identical(attr(res, "fallback"), 1:10)

  # every window was skipped, so the result is the input, cell for cell, NA
  # included
  expect_identical(res[,], mat)
  expect_equal(attr(res, "fallback_action"), "skip")
})

test_that("slide_imp: on_infeasible = 'error' names the window that failed", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  # window 3 is columns 21:30
  mat[1:19, 21:30] <- NA
  mat[20, 21:30] <- rnorm(10)

  err <- expect_error(
    suppressMessages(slide_imp(
      mat,
      location = 1:100,
      k = 3,
      window_size = 10,
      overlap_size = 0,
      min_window_n = 5,
      on_infeasible = "error",
      .progress = FALSE
    )),
    class = "slideimp_infeasible"
  )
  expect_match(
    conditionMessage(err),
    "Window 3 (columns 21 to 30)",
    fixed = TRUE
  )
  # the backend's own condition travels with it
  expect_s3_class(err$parent, "slideimp_infeasible")

  err_flank <- expect_error(
    suppressMessages(slide_imp(
      mat,
      location = 1:100,
      k = 3,
      window_size = 4,
      flank = TRUE,
      subset = 25L,
      min_window_n = 5,
      on_infeasible = "error",
      .progress = FALSE
    )),
    class = "slideimp_infeasible"
  )
  expect_match(conditionMessage(err_flank), "target column 25", fixed = TRUE)
})

test_that("slide_imp's mean fallback is handed only the targets with a hole", {
  # Regression, moved from mean_imp_col(): the fallback part over a window's
  # targets took the complete ones too, so fallback_part() copied them and
  # their is.na() to fill nothing. A complete column in a filling part still
  # validates, so only the columns it is handed show the defect.
  x <- matrix(rnorm(20 * 10), 20, 10)
  x[1:19, 1:9] <- NA # over colmax: the one window is infeasible
  expect_false(anyNA(x[, 10]))
  real_fallback_part <- fallback_part
  calls <- list()
  local_mocked_bindings(
    fallback_part = function(obj, cols, method = c("mean", "median")) {
      calls <<- c(calls, list(cols))
      real_fallback_part(obj, cols, method)
    }
  )
  res <- suppressWarnings(slide_imp(
    x,
    location = seq_len(10),
    window_size = 10,
    min_window_n = 5,
    k = 3,
    on_infeasible = "mean",
    .progress = FALSE
  ))
  expect_identical(attr(res, "fallback"), 1L)
  expect_identical(calls, list(1:9))
  expect_false(anyNA(res[,]))
  expect_identical(res[, 10], x[, 10])
})
