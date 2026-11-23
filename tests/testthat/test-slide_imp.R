skip_on_cran()

test_that("slide_imp knn mode works", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## knn_imp, which we test correctness elsewhere
  to_test <- sim_mat(100, 280, rho = 0.5, perc_col_na = 1)$input
  # Init
  counts <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  final_imputed <- counts
  # 1 to 100 is the first window;
  final_imputed[, 1:100] <- final_imputed[, 1:100] +
    knn_imp(
      obj = to_test[, 1:100],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 91 to 190 is the second window;
  final_imputed[, 91:190] <- final_imputed[, 91:190] +
    knn_imp(
      obj = to_test[, 91:190],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )
  counts[, 91:190] <- counts[, 91:190] + 1
  # 181 to 280 is the last window
  final_imputed[, 181:280] <- final_imputed[, 181:280] +
    knn_imp(
      obj = to_test[, 181:280],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )
  counts[, 181:280] <- counts[, 181:280] + 1
  final_imputed <- final_imputed / counts

  # slide_imp should exactly replicate this result
  location <- 1:ncol(to_test)
  simple_mean <- slide_imp(
    to_test,
    location = location,
    window_size = 100,
    overlap_size = 10,
    k = 3,
    min_window_n = 10,
    colmax = 0.9,
    post_imp = TRUE
  )

  expect_worst_err(simple_mean[,], final_imputed)
  expect_identical(dimnames(simple_mean), dimnames(final_imputed))

  # slide_imp weighted should be different than simple mean
  weighted_1 <- slide_imp(
    to_test,
    location = location,
    window_size = 100,
    overlap_size = 10,
    k = 3,
    min_window_n = 10,
    colmax = 0.9,
    post_imp = TRUE,
    dist_pow = 1
  )
  weighted_2 <- slide_imp(
    to_test,
    location = location,
    window_size = 100,
    overlap_size = 10,
    k = 3,
    min_window_n = 10,
    colmax = 0.9,
    post_imp = TRUE,
    dist_pow = 2
  )
  expect_true(sum((simple_mean[,] - weighted_1[,])^2) > 0)
  expect_true(sum((weighted_2[,] - weighted_1[,])^2) > 0)
})

test_that("slide_imp subset works", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## knn_imp, which we test correctness elsewhere
  to_test <- sim_mat(10, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  subset <- c(1, 6, 10, 50)
  # Init
  counts <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  final_imputed <- counts
  # 1 to 20 is the first window;
  window_cols <- 1:20
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols[local_subset]] <- final_imputed[, window_cols[
    local_subset
  ]] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # 16 to 35 is the second window. It covers no subset column, so slide_imp()
  # drops it, and knn_imp() refuses an empty subset; it would add nothing to
  # the subset columns compared below.
  expect_length(which(16:35 %in% subset), 0L)
  # 31 to 50 is the last window
  window_cols <- 31:50
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols[local_subset]] <- final_imputed[, window_cols[
    local_subset
  ]] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  final_imputed <- final_imputed / counts
  # slide_imp should exactly replicate this result, and return only the
  # subset columns
  location <- 1:ncol(to_test)
  res <- slide_imp(
    to_test,
    location = location,
    window_size = 20,
    overlap_size = 5,
    k = 3,
    min_window_n = 10,
    colmax = 0.9,
    post_imp = TRUE,
    subset = subset
  )
  want <- final_imputed[, subset, drop = F]
  expect_worst_err(res[,], want)
  expect_identical(dimnames(res), dimnames(want))
})

test_that("slide_imp edge case no overlap", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## knn_imp, which we test correctness elsewhere
  to_test <- sim_mat(100, 300, perc_total_na = 0.5, perc_col_na = 1)$input
  # Init
  counts <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  final_imputed <- counts

  # 1 to 100 is the first window;
  final_imputed[, 1:100] <- final_imputed[, 1:100] +
    knn_imp(
      obj = to_test[, 1:100],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 101 to 200 is the second window;
  final_imputed[, 101:200] <- final_imputed[, 101:200] +
    knn_imp(
      obj = to_test[, 101:200],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )
  counts[, 101:200] <- counts[, 101:200] + 1
  # 201 to 300 is the last window
  final_imputed[, 201:300] <- final_imputed[, 201:300] +
    knn_imp(
      obj = to_test[, 201:300],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )
  counts[, 201:300] <- counts[, 201:300] + 1
  final_imputed <- final_imputed / counts
  # slide_imp should exactly replicate this result
  location <- 1:ncol(to_test)
  res <- slide_imp(
    to_test,
    location = location,
    window_size = 100,
    overlap_size = 0,
    min_window_n = 10,
    k = 3,
    colmax = 0.9,
    post_imp = TRUE
  )
  expect_worst_err(res[,], final_imputed)
  expect_identical(dimnames(res), dimnames(final_imputed))
})

test_that("slide_imp pca mode works", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## pca_imp, which we test correctness elsewhere
  # Every fit here runs a fixed iteration count, miniter = maxiter, so the
  # windows and the hand-built reference take the same steps; the threshold,
  # far above the criterion at that count, decides only `converged`.
  to_test <- sim_mat(20, 280, perc_total_na = 0.5, perc_col_na = 1)$input
  # Init
  counts <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  final_imputed <- counts
  # 1 to 100 is the first window;
  final_imputed[, 1:100] <- final_imputed[, 1:100] +
    pca_imp(
      obj = to_test[, 1:100],
      ncp = 2,
      threshold = 0.1,
      maxiter = 30,
      miniter = 30,
      seed = 1234
    )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 91 to 190 is the second window;
  final_imputed[, 91:190] <- final_imputed[, 91:190] +
    pca_imp(
      obj = to_test[, 91:190],
      ncp = 2,
      threshold = 0.1,
      maxiter = 30,
      miniter = 30,
      seed = 1234
    )
  counts[, 91:190] <- counts[, 91:190] + 1
  # 181 to 280 is the last window
  final_imputed[, 181:280] <- final_imputed[, 181:280] +
    pca_imp(
      obj = to_test[, 181:280],
      ncp = 2,
      threshold = 0.1,
      maxiter = 30,
      miniter = 30,
      seed = 1234
    )
  counts[, 181:280] <- counts[, 181:280] + 1
  final_imputed <- final_imputed / counts
  set.seed(1234)
  # slide_imp should exactly replicate this result
  location <- 1:ncol(to_test)
  simple_mean <- slide_imp(
    to_test,
    location = location,
    window_size = 100,
    overlap_size = 10,
    min_window_n = 10,
    ncp = 2,
    threshold = 0.1,
    maxiter = 30,
    miniter = 30,
    seed = 1234
  )
  expect_worst_err(simple_mean[,], final_imputed)
  expect_identical(dimnames(simple_mean), dimnames(final_imputed))
})

test_that("slide_imp handling of errors on zero-variance features in PCA mode", {
  set.seed(1234)
  to_test <- sim_mat(10, 200, perc_total_na = 0.5, perc_col_na = 1)$input
  to_test[, 1] <- 1
  location <- 1:ncol(to_test)
  expect_no_error(
    slide_imp(
      to_test,
      location = location,
      window_size = 100,
      overlap_size = 10,
      min_window_n = 10,
      ncp = 2,
      miniter = 2
    )
  )
})

test_that("slide_imp does not pin the solver on a window that was never probed", {
  set.seed(1234)
  # 260 rows, so a full-width flank window clears pca_imp()'s n_gram >= 250
  # threshold and gets probed. The window centered on column 1 is truncated to
  # 201 columns, so it is demoted to exact WITHOUT a probe.
  x <- sim_mat(260, 1000, perc_total_na = 0.1, perc_col_na = 1)$input
  location <- seq_len(ncol(x))
  targets <- c(1L, 400L, 800L)

  args <- list(
    x,
    location = location,
    window_size = 200,
    flank = TRUE,
    subset = targets,
    min_window_n = 50,
    .progress = FALSE
  )

  windows <- do.call(slide_imp, c(args, list(dry_run = TRUE)))
  first_cols <- windows$start[[1L]]:windows$end[[1L]]
  expect_lt(length(first_cols), 250L)
  expect_gt(nrow(windows), 1L)

  # the fixture really is the case in question: window 1 reports "exact"
  # without ever probing.
  first <- pca_imp(
    x[, first_cols, drop = FALSE],
    ncp = 2,
    solver = "auto",
    seed = 1
  )
  expect_identical(attr(first, "solver_chosen"), "exact")
  expect_false(attr(first, "solver_probed"))

  res <- suppressMessages(do.call(
    slide_imp,
    c(args, list(ncp = 2, solver = "auto", seed = 1))
  ))

  # Regression: slide_imp pinned every later window to window 1's reported
  # `solver_chosen`, so an unprobed "exact" locked at window 1 and the wider
  # windows, which auto would have probed, never got the chance. Measured on
  # 400 x 4000 with eight windows: 1.66s pinned against 1.25s unpinned, the
  # pinned figure being bit-for-bit the cost of forcing solver = "exact".
  # Which solver a probed window picks is timing-based and not assertable;
  # that the lock did not land on window 1 is.
  expect_false(identical(attr(res, "solver_lock_window"), 1L))
})

test_that("slide_imp flank works with knn", {
  set.seed(1234)
  to_test <- sim_mat(10, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  location <- 1:ncol(to_test)
  subset <- c(5, 25, 45)
  window_size <- 20
  min_window_n <- 10

  # the windows, written out by hand: window_size columns either side of each
  # target, cut at the ends of `location`. All three are at least
  # min_window_n wide, so none is dropped.
  start <- c(1L, 5L, 25L)
  end <- c(25L, 45L, 50L)
  subset_local <- c(5L, 21L, 21L)
  windows <- slide_imp(
    to_test,
    location = location,
    window_size = window_size,
    flank = TRUE,
    min_window_n = min_window_n,
    subset = subset,
    dry_run = TRUE,
    .progress = FALSE
  )
  expect_identical(windows$start, start)
  expect_identical(windows$end, end)

  result <- to_test
  for (i in seq_along(start)) {
    window_cols <- start[i]:end[i]
    imputed_window <- knn_imp(
      obj = to_test[, window_cols, drop = FALSE],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = subset_local[i]
    )
    local_idx <- subset_local[i]
    result[, window_cols[local_idx]] <- imputed_window
  }

  res <- slide_imp(
    to_test,
    location = location,
    window_size = window_size,
    flank = TRUE,
    k = 3,
    min_window_n = min_window_n,
    colmax = 0.9,
    post_imp = TRUE,
    subset = subset,
    .progress = FALSE
  )
  expect_worst_err(res[,], result[, subset])
  expect_identical(dimnames(res), dimnames(result[, subset]))
})

test_that("slide_imp PCA skips windows not covering any subset features", {
  set.seed(1234)
  to_test <- sim_mat(100, 65, perc_total_na = 0.5, perc_col_na = 1)$input
  location <- 1:ncol(to_test)
  # Windows of 20 with overlap 5 are [1, 20], [16, 35], [31, 50] and
  # [46, 65]. [16, 35] holds no subset column and is dropped. Column 48 sits
  # under both [31, 50] and [46, 65], so its holes are the mean of the two.
  subset <- c(1, 6, 48, 65)
  windows <- slide_imp(
    to_test,
    location = location,
    window_size = 20,
    overlap_size = 5,
    min_window_n = 10,
    subset = subset,
    dry_run = TRUE,
    .progress = FALSE
  )
  expect_identical(windows$start, c(1L, 31L, 46L))
  expect_identical(windows$end, c(20L, 50L, 65L))
  expect_true(anyNA(to_test[, 48]))

  # every fit runs a fixed iteration count, as in "slide_imp pca mode works"
  counts <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  final_imputed <- counts
  # window 1: 1-20 - covers subset cols 1, 6
  window_cols <- 1:20
  imp <- pca_imp(
    obj = to_test[, window_cols],
    ncp = 2,
    scale = TRUE,
    seed = 1,
    threshold = 1e-2,
    maxiter = 30,
    miniter = 30
  )
  local_sub <- intersect(subset, window_cols)
  final_imputed[, local_sub] <- final_imputed[, local_sub] +
    imp[, match(local_sub, window_cols), drop = FALSE]
  counts[, local_sub] <- counts[, local_sub] + 1

  # window 2: 16-35 - no subset features, SKIPPED

  # window 3: 31-50 - covers subset col 48
  window_cols <- 31:50
  imp <- pca_imp(
    obj = to_test[, window_cols],
    ncp = 2,
    scale = TRUE,
    seed = 1,
    threshold = 1e-2,
    maxiter = 30,
    miniter = 30
  )
  local_sub <- intersect(subset, window_cols)
  final_imputed[, local_sub] <- final_imputed[, local_sub] +
    imp[, match(local_sub, window_cols), drop = FALSE]
  counts[, local_sub] <- counts[, local_sub] + 1

  # window 4: 46-65 - covers subset cols 48, 65
  window_cols <- 46:65
  imp <- pca_imp(
    obj = to_test[, window_cols],
    ncp = 2,
    scale = TRUE,
    seed = 1,
    threshold = 1e-2,
    maxiter = 30,
    miniter = 30
  )
  local_sub <- intersect(subset, window_cols)
  final_imputed[, local_sub] <- final_imputed[, local_sub] +
    imp[, match(local_sub, window_cols), drop = FALSE]
  counts[, local_sub] <- counts[, local_sub] + 1

  # divide each subset column by the number of windows holding it. `counts`
  # holds that number in every row, so read it from row 1.
  expect_identical(unname(counts[1, subset]), c(1, 1, 2, 1))
  final_imputed[, subset] <- final_imputed[, subset] / counts[, subset]
  res <- slide_imp(
    to_test,
    location = location,
    window_size = 20,
    overlap_size = 5,
    ncp = 2,
    min_window_n = 10,
    scale = TRUE,
    seed = 1,
    threshold = 1e-2,
    maxiter = 30,
    miniter = 30,
    subset = subset
  )
  expect_worst_err(res[,], final_imputed[, subset])
  expect_identical(dimnames(res), dimnames(final_imputed[, subset]))
})

test_that("slide_imp flank works with pca", {
  set.seed(1234)
  to_test <- sim_mat(10, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  location <- 1:ncol(to_test)
  subset <- c(5, 25, 45)
  window_size <- 20
  min_window_n <- 10

  # the windows, written out by hand, as in the K-NN flank test above
  start <- c(1L, 5L, 25L)
  end <- c(25L, 45L, 50L)
  subset_local <- c(5L, 21L, 21L)
  windows <- slide_imp(
    to_test,
    location = location,
    window_size = window_size,
    flank = TRUE,
    min_window_n = min_window_n,
    subset = subset,
    dry_run = TRUE,
    .progress = FALSE
  )
  expect_identical(windows$start, start)
  expect_identical(windows$end, end)

  # every fit runs a fixed iteration count, as in "slide_imp pca mode works"
  result <- to_test
  for (i in seq_along(start)) {
    window_cols <- start[i]:end[i]
    imputed_window <- pca_imp(
      obj = to_test[, window_cols, drop = FALSE],
      ncp = 2,
      scale = TRUE,
      method = "regularized",
      seed = 1234,
      threshold = 0.1,
      maxiter = 30,
      miniter = 30
    )
    local_idx <- subset_local[i]
    result[, window_cols[local_idx]] <- imputed_window[, local_idx]
  }

  res <- slide_imp(
    to_test,
    location = location,
    window_size = window_size,
    flank = TRUE,
    ncp = 2,
    min_window_n = min_window_n,
    scale = TRUE,
    subset = subset,
    seed = 1234,
    threshold = 0.1,
    maxiter = 30,
    miniter = 30,
    .progress = FALSE
  )
  expect_worst_err(res[,], result[, subset])
  expect_identical(dimnames(res), dimnames(result[, subset]))
})

test_that("slide_imp matches a hand-built reference on gappy locations", {
  set.seed(1234)
  # three blocks of columns separated by gaps wider than a window. Blocks 1 and
  # 2 hold 25 columns each, block 3 holds 10 that are 96.7% missing.
  location <- c(1:25, 41:65, 81:90)
  x <- matrix(rnorm(30 * 60, mean = 5), 30, 60)
  x[sample(30 * 50, 150)] <- NA
  x[1:29, 51:60] <- NA

  args <- list(
    x,
    location = location,
    window_size = 10,
    overlap_size = 5,
    min_window_n = 5,
    .progress = FALSE
  )

  # the windows, written out by hand. Each block ends on a window that the
  # trailing overlap cannot extend, which is where a nested window used to be.
  starts <- c(1L, 6L, 11L, 16L, 26L, 31L, 36L, 41L, 51L)
  ends <- c(10L, 15L, 20L, 25L, 35L, 40L, 45L, 50L, 60L)
  windows <- do.call(slide_imp, c(args, list(dry_run = TRUE)))
  expect_identical(windows$start, starts)
  expect_identical(windows$end, ends)

  # reference: impute each window alone, add up, divide by the cover count
  sums <- matrix(0, nrow(x), ncol(x))
  counts <- integer(ncol(x))
  for (i in 1:8) {
    cols <- starts[i]:ends[i]
    sums[, cols] <- sums[, cols] +
      suppressMessages(knn_imp(x[, cols], k = 2, .progress = FALSE))[,]
    counts[cols] <- counts[cols] + 1L
  }
  ref <- sweep(sums, 2, pmax(counts, 1L), "/")
  # window 9 is infeasible, so its columns keep the input
  expect_error(
    suppressMessages(knn_imp(x[, 51:60], k = 2, .progress = FALSE)),
    class = "slideimp_infeasible"
  )
  ref[, 51:60] <- x[, 51:60]
  # the fixture has both singly and doubly covered columns
  expect_true(any(counts == 1L))
  expect_true(any(counts == 2L))

  expect_warning(
    res <- suppressMessages(do.call(slide_imp, c(args, list(k = 2)))),
    class = "slideimp_fallback"
  )
  expect_identical(attr(res, "fallback"), 9L)

  out <- res[,]
  dimnames(out) <- NULL
  expect_identical(is.na(out), is.na(ref))
  ok <- !is.na(ref)
  # the 290 missing cells of block 3 stay NA. Both errors are taken over every
  # other cell, 1510 of them.
  expect_identical(sum(!ok), 290L)
  expect_lte(max(abs(out[ok] - ref[ok])), 1e-12)
  expect_lte(max(abs(out[ok] - ref[ok]) / abs(ref[ok])), 1e-12)
})

test_that("slide_imp PCA does not call pca_imp on a window with complete targets", {
  set.seed(2)
  y <- matrix(rnorm(20 * 30), 20, 30)
  y[1:19, 5] <- NA # column 5 is 95% missing, above colmax = 0.9
  expect_false(anyNA(y[, 15]))

  # the flank window of target 15 is columns 5:25. Its only missing values sit
  # in an ineligible column, so pca_imp() itself refuses it.
  expect_error(
    suppressMessages(pca_imp(y[, 5:25], ncp = 2, .progress = FALSE)),
    class = "slideimp_infeasible"
  )

  # Regression: slide_imp() handed that window to pca_imp() although the target
  # needed nothing, so `on_infeasible = "error"` aborted the whole run and
  # "skip" recorded a fallback. K-NN on the same input did neither.
  res <- suppressMessages(slide_imp(
    y,
    location = seq_len(30),
    window_size = 10,
    flank = TRUE,
    subset = 15L,
    min_window_n = 5,
    ncp = 2,
    on_infeasible = "error",
    .progress = FALSE
  ))
  expect_length(attr(res, "fallback"), 0L)
  expect_identical(res[, 1], y[, 15])
  # no window was probed, and the attributes keep one type either way
  expect_identical(attr(res, "solver_chosen"), NA_character_)
  expect_identical(attr(res, "solver_lock_window"), NA_integer_)
})

test_that("slide_imp with subset returns the full result's subset columns", {
  # The windows impute their targets from the same columns whatever the
  # targets are, so narrowing changes what is returned, not the values: the
  # narrowed result is the full result's columns, the record renumbered to
  # them. A window that covers no subset column is dropped, which adds
  # nothing to the subset columns.
  set.seed(42)
  x <- sim_mat(30, 80, perc_total_na = 0.2, perc_col_na = 1)$input
  x[1:28, 27] <- NA # over colmax, filled with its mean
  x[7, 41:60] <- NA # a row no neighbor holds in the window of columns 41:60
  # Column 28 shares every window with 27. A window whose only targets are
  # over colmax is infeasible, so 27 alone would leave some of its windows
  # out of the narrowed run and not of the full one.
  sel <- c(27L, 3L, 44L, 28L, 79L, 50L)
  ord <- sort(sel)
  run <- function(...) {
    suppressMessages(slide_imp(
      x,
      location = seq_len(ncol(x)),
      window_size = 20,
      overlap_size = 12,
      min_window_n = 10,
      .progress = FALSE,
      ...
    ))
  }
  # the PCA fits run a fixed iteration count, miniter = maxiter, so the three
  # runs take the same steps; the threshold, far above the criterion at that
  # count, decides only `converged`
  pca_args <- list(
    ncp = 2,
    solver = "exact",
    threshold = 1e-3,
    maxiter = 30,
    miniter = 30
  )
  for (args in list(list(k = 4), pca_args)) {
    full <- do.call(run, args)
    narrow <- do.call(run, c(args, list(subset = sel)))
    by_name <- do.call(run, c(args, list(subset = colnames(x)[sel])))
    if ("k" %in% names(args)) {
      # K-NN reaches no BLAS, so one process gives the same bits
      expect_identical(narrow[,], full[, ord])
      expect_identical(by_name, narrow)
    } else {
      # PCA goes through BLAS, which need not return the same bits for the
      # same call twice: the values are asserted on a tolerance and bitwise
      # reported, the attributes, all structural, stay identical
      expect_worst_err(narrow[,], full[, ord])
      expect_worst_err(by_name[,], narrow[,])
      expect_identical(attributes(by_name), attributes(narrow))
      message(sprintf(
        "subset PCA bitwise: narrow vs full %s, by_name vs narrow %s",
        identical(narrow[,], full[, ord]),
        identical(by_name, narrow)
      ))
    }
    expect_identical(colnames(narrow), colnames(x)[ord])

    # the record, in the result's numbering: column 27 is the result's 2nd
    full_cols <- attr(full, "mean_filled_cols")
    expect_true(27L %in% full_cols)
    expect_identical(
      attr(narrow, "mean_filled_cols"),
      match(intersect(full_cols, ord), ord)
    )
    expect_identical(attr(narrow, "mean_filled_cols"), 2L)
    full_cells <- attr(full, "mean_filled_cells")
    want_cells <- full_cells[full_cells[, "col"] %in% ord, , drop = FALSE]
    want_cells[, "col"] <- match(want_cells[, "col"], ord)
    expect_identical(attr(narrow, "mean_filled_cells"), want_cells)
    if ("k" %in% names(args)) {
      # K-NN leaves row 7 of the columns 41:60, column 44 and 50 among them
      expect_true(all(
        c(44L, 50L) %in% full_cells[full_cells[, "row"] == 7L, "col"]
      ))
      expect_gt(nrow(want_cells), 0L)
    }
    nms <- setdiff(
      names(attributes(full)),
      c("dim", "dimnames", "mean_filled_cols", "mean_filled_cells")
    )
    expect_identical(attributes(narrow)[nms], attributes(full)[nms])
  }
})

test_that("slide_imp leaves observed cells as they are under any overlap", {
  # Regression: slide_imp() added every window's copy of a column, observed
  # cells included, and divided by the number of windows. Under three
  # windows ((x + x) + x) / 3 can round away from x, so observed values came
  # back changed in their last bit. Only the holes are averaged now.
  set.seed(1)
  x <- sim_mat(n = 60, p = 120, perc_total_na = 0.1, perc_col_na = 0.5)$input
  args <- list(
    x,
    location = seq_len(ncol(x)),
    window_size = 30,
    overlap_size = 21,
    min_window_n = 10,
    .progress = FALSE
  )
  w <- do.call(slide_imp, c(args, list(dry_run = TRUE)))
  cover <- count_window_cover(w$start, w$end, ncol(x))
  cols3 <- which(cover >= 3L)
  expect_gt(length(cols3), 0L)
  obs <- !is.na(x)
  obs3 <- obs
  obs3[, -cols3] <- FALSE
  # the fixture holds observed values that ((x + x) + x) / 3 changes
  x3 <- x[obs3]
  expect_gt(sum(((x3 + x3) + x3) / 3 != x3), 0L)

  for (extra in list(list(k = 5), list(ncp = 2, solver = "exact"))) {
    res <- suppressMessages(do.call(slide_imp, c(args, extra)))
    expect_identical(unclass(res)[obs], x[obs])
    expect_false(anyNA(res[,]))
  }
})

test_that("slide_imp averages a hole over the windows, NA when one leaves it", {
  # With post_imp = FALSE, K-NN can leave a hole: here row 7 of columns
  # 21:30, where every neighbor inside the window of columns 21:40 is
  # missing. Those columns are also in the window of columns 11:30, which
  # fills them: with k = 10, at least one neighbor of each is in columns
  # 11:20, all observed at row 7. So the hole stays NA only because one
  # window left it, as it did when the windows' matrices were summed. The
  # other holes of those columns are the mean over the two windows.
  set.seed(7)
  x <- sim_mat(30, 50, perc_total_na = 0.1, perc_col_na = 1)$input
  row7 <- x[7, 1:20]
  row7[is.na(row7)] <- 0.5
  x[7, 1:20] <- row7
  x[7, 21:40] <- NA
  run_args <- list(
    location = seq_len(ncol(x)),
    window_size = 20,
    overlap_size = 10,
    min_window_n = 12,
    k = 10,
    post_imp = FALSE,
    .progress = FALSE
  )
  w <- do.call(slide_imp, c(list(x, dry_run = TRUE), run_args))
  expect_identical(w$start[2:3], c(11L, 21L))
  expect_identical(w$end[2:3], c(30L, 40L))
  # no other window covers columns 21:30
  expect_identical(sum(w$start <= 21L & w$end >= 30L), 2L)
  expect_identical(sum(w$start <= 30L & w$end >= 21L), 2L)
  res <- suppressMessages(do.call(slide_imp, c(list(x), run_args)))

  win <- function(cols) {
    suppressMessages(knn_imp(
      x[, cols],
      k = 10,
      post_imp = FALSE,
      .progress = FALSE
    ))[,]
  }
  a <- win(11:30)[, 11:20]
  b <- win(21:40)[, 1:10]
  # window 21:40 leaves row 7, window 11:30 fills it
  expect_true(all(is.na(b[7, ])))
  expect_false(anyNA(a[7, ]))
  expect_true(all(is.na(res[7, 21:30])))
  # every other hole of columns 21:30 is the mean of the two windows
  holes <- is.na(x[, 21:30])
  holes[7, ] <- FALSE
  expect_gt(sum(holes), 0L)
  want <- (a + b) / 2
  expect_identical(res[, 21:30][holes], want[holes])
  # and no observed cell moved
  obs <- !is.na(x)
  expect_identical(unclass(res)[obs], x[obs])
})

test_that("slide_imp records mean-filled columns and cells by their column in the result", {
  set.seed(1234)
  mat <- matrix(rnorm(20 * 60), 20, 60)
  windows <- slide_imp(
    mat,
    location = 1:60,
    window_size = 20,
    overlap_size = 5,
    min_window_n = 5,
    dry_run = TRUE,
    .progress = FALSE
  )
  # the premises: a window that does not start at column 1, so a local index
  # differs from the index in `mat`, and a later window that holds column 52
  expect_identical(c(windows$start[3], windows$end[3]), c(31L, 50L))
  w52 <- which(windows$start <= 52L & windows$end >= 52L)
  expect_length(w52, 1L)
  expect_gt(windows$start[w52], 1L)

  mat[1:19, 52] <- NA # miss_rate 0.95, above the default colmax
  mat[2:3, 10] <- NA # ordinary missingness, the method fills it
  # row 5 is missing across the whole of window 3, so inside that window no
  # neighbor has a value there
  mat[5, 31:50] <- NA
  no_cells <- cbind(row = integer(0), col = integer(0))

  run <- function(...) {
    suppressMessages(slide_imp(
      mat,
      location = 1:60,
      window_size = 20,
      overlap_size = 5,
      min_window_n = 5,
      .progress = FALSE,
      ...
    ))
  }

  knn <- run(k = 3)
  expect_false(anyNA(knn[,]))
  expect_identical(attr(knn, "mean_filled_cols"), 52L)
  expect_identical(
    attr(knn, "mean_filled_cells"),
    cbind(row = rep(5L, 20L), col = 31:50)
  )
  expect_equal(unname(knn[1:19, 52]), rep(mat[20, 52], 19L))

  # nothing is filled without post_imp, so nothing is recorded
  knn_off <- run(k = 3, post_imp = FALSE)
  expect_true(anyNA(knn_off[, 52]))
  expect_identical(attr(knn_off, "mean_filled_cols"), integer(0))
  expect_identical(attr(knn_off, "mean_filled_cells"), no_cells)

  pca <- run(ncp = 2)
  expect_identical(attr(pca, "mean_filled_cols"), 52L)
  # with `subset`, the result holds only those columns. Column 52 is not one
  # of them, so nothing of it is recorded.
  pca_sub <- run(ncp = 2, subset = setdiff(1:60, 52L))
  expect_identical(ncol(pca_sub), 59L)
  expect_identical(attr(pca_sub, "mean_filled_cols"), integer(0))

  # with `subset = 41:60`, column 52 is the result's 12th and columns 41:50
  # its 1st to 10th, so the record counts them so
  knn_sub <- run(k = 3, subset = 41:60)
  expect_identical(dim(knn_sub), c(20L, 20L))
  expect_identical(attr(knn_sub, "mean_filled_cols"), 12L)
  expect_identical(
    attr(knn_sub, "mean_filled_cells"),
    cbind(row = rep(5L, 10L), col = 1:10)
  )
  expect_identical(knn_sub[,], knn[, 41:60])
})

test_that("hole_row_index() returns integer rows from double positions", {
  # Regression: the rows of a window's holes were taken from which() with
  # integer arithmetic, so when which() returned double positions, as it does
  # for a matrix past 2^31 - 1 cells, the rows stayed double and the part
  # built from them failed validation at the final write.
  at <- c(1, 4, 7, 12)
  expect_true(is.double(at))
  rows <- hole_row_index(at, 5L)
  expect_identical(rows, c(1L, 4L, 2L, 2L))
  expect_identical(hole_row_index(as.integer(at), 5L), rows)
  # the door refuses double rows, so the conversion is what makes them pass
  x <- matrix(as.double(c(NA, 1:4, NA, 6:9)), 5, 2)
  rows_of <- function(r) {
    slide_parts(
      subset = 1:2,
      hole_n = c(1L, 1L),
      acc = c(0.5, 0.25),
      n_ran = c(1L, 1L),
      seen = c(TRUE, TRUE),
      rows = r,
      kind = "knn"
    )
  }
  expect_error(
    fill_parts(x, rows_of(c(1, 1))),
    "'map' is not an integer matrix",
    fixed = TRUE
  )
  out <- fill_parts(x, rows_of(hole_row_index(c(1, 6), 5L)))
  expect_identical(out[c(1, 6)], c(0.5, 0.25))
})
