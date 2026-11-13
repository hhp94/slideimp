test_that("slide_imp knn mode works", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## knn_imp, which we test correctness elsewhere
  to_test <- sim_mat(100, 280, 0.5, perc_col_na = 1)$input
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

  expect_equal(simple_mean[,], final_imputed)

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
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # 16 to 35 is the second window;
  window_cols <- 16:35
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # 31 to 50 is the last window
  window_cols <- 31:50
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  final_imputed <- final_imputed / counts
  # slide_imp should exactly replicate this result
  location <- 1:ncol(to_test)
  expect_equal(
    slide_imp(
      to_test,
      location = location,
      window_size = 20,
      overlap_size = 5,
      k = 3,
      min_window_n = 10,
      colmax = 0.9,
      post_imp = TRUE,
      subset = subset
    )[, subset, drop = F],
    final_imputed[, subset, drop = F]
  )
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
  expect_equal(
    slide_imp(
      to_test,
      location = location,
      window_size = 100,
      overlap_size = 0,
      min_window_n = 10,
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )[,],
    final_imputed
  )
})

test_that("slide_imp pca mode works", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## pca_imp, which we test correctness elsewhere
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
      miniter = 2,
      seed = 1234
    )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 91 to 190 is the second window;
  final_imputed[, 91:190] <- final_imputed[, 91:190] +
    pca_imp(
      obj = to_test[, 91:190],
      ncp = 2,
      miniter = 2,
      seed = 1234
    )
  counts[, 91:190] <- counts[, 91:190] + 1
  # 181 to 280 is the last window
  final_imputed[, 181:280] <- final_imputed[, 181:280] +
    pca_imp(
      obj = to_test[, 181:280],
      ncp = 2,
      miniter = 2,
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
    miniter = 2,
    seed = 1234
  )
  expect_equal(simple_mean[,], final_imputed)
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

test_that("slide_imp flank works with knn", {
  set.seed(1234)
  to_test <- sim_mat(10, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  location <- 1:ncol(to_test)
  subset <- c(5, 25, 45)
  window_size <- 20
  min_window_n <- 10

  fw <- find_windows_flank(location, subset, window_size)
  start <- fw$start
  end <- fw$end
  subset_local <- fw$subset_local

  window_n <- end - start + 1L
  keep <- window_n >= min_window_n
  start <- start[keep]
  end <- end[keep]
  subset_local <- subset_local[keep]
  target_cols <- subset[keep]

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
    result[, window_cols[local_idx]] <- imputed_window[, local_idx]
  }

  expect_equal(
    slide_imp(
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
    )[,],
    result
  )
})

test_that("slide_imp K-NN skips windows not covering any subset features", {
  set.seed(1234)
  to_test <- sim_mat(10, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  subset <- c(1, 6, 45, 50)
  counts <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  final_imputed <- counts
  # window 1: 1 to 20 - covers subset cols 1, 6
  window_cols <- 1:20
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # window 2: 16 to 35 - no subset features, SKIPPED
  # window 3: 31 to 50 - covers subset cols 45, 50
  window_cols <- 31:50
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # average overlaps, restore originals where uncovered
  for (j in which(colSums(counts) > 1)) {
    final_imputed[, j] <- final_imputed[, j] / counts[, j]
  }
  uncovered <- which(colSums(counts) == 0)
  final_imputed[, uncovered] <- to_test[, uncovered]
  location <- 1:ncol(to_test)
  expect_equal(
    slide_imp(
      to_test,
      location = location,
      window_size = 20,
      overlap_size = 5,
      k = 3,
      min_window_n = 10,
      colmax = 0.9,
      post_imp = TRUE,
      subset = subset
    )[,],
    final_imputed[,]
  )
})

test_that("slide_imp PCA skips windows not covering any subset features", {
  set.seed(1234)
  to_test <- sim_mat(100, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  subset <- c(1, 6, 45, 50)
  counts <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  final_imputed <- counts
  # window 1: 1-20 - covers subset cols 1, 6
  window_cols <- 1:20
  imp <- pca_imp(obj = to_test[, window_cols], ncp = 2, scale = TRUE, seed = 1)
  local_sub <- intersect(subset, window_cols)
  final_imputed[, local_sub] <- final_imputed[, local_sub] +
    imp[, match(local_sub, window_cols), drop = FALSE]
  counts[, local_sub] <- counts[, local_sub] + 1

  # window 2: 16-35 - no subset features, SKIPPED

  # window 3: 31-50 - covers subset cols 45, 50
  window_cols <- 31:50
  imp <- pca_imp(obj = to_test[, window_cols], ncp = 2, scale = TRUE, seed = 1)
  local_sub <- intersect(subset, window_cols)
  final_imputed[, local_sub] <- final_imputed[, local_sub] +
    imp[, match(local_sub, window_cols), drop = FALSE]
  counts[, local_sub] <- counts[, local_sub] + 1
  # average overlaps, restore originals where uncovered
  for (j in which(colSums(counts) > 1)) {
    final_imputed[, j] <- final_imputed[, j] / counts[, j]
  }
  uncovered <- which(colSums(counts) == 0)
  final_imputed[, uncovered] <- to_test[, uncovered]
  location <- 1:ncol(to_test)
  expect_equal(
    slide_imp(
      to_test,
      location = location,
      window_size = 20,
      overlap_size = 5,
      ncp = 2,
      min_window_n = 10,
      scale = TRUE,
      seed = 1,
      subset = subset
    )[,],
    final_imputed[,]
  )
})

test_that("slide_imp flank works with pca", {
  set.seed(1234)
  to_test <- sim_mat(10, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  location <- 1:ncol(to_test)
  subset <- c(5, 25, 45)
  window_size <- 20
  min_window_n <- 10

  fw <- find_windows_flank(location, subset, window_size)
  start <- fw$start
  end <- fw$end
  subset_local <- fw$subset_local

  window_n <- end - start + 1L
  keep <- window_n >= min_window_n
  start <- start[keep]
  end <- end[keep]
  subset_local <- subset_local[keep]
  target_cols <- subset[keep]

  result <- to_test
  for (i in seq_along(start)) {
    window_cols <- start[i]:end[i]
    imputed_window <- pca_imp(
      obj = to_test[, window_cols, drop = FALSE],
      ncp = 2,
      scale = TRUE,
      method = "regularized",
      seed = 1234
    )
    local_idx <- subset_local[i]
    result[, window_cols[local_idx]] <- imputed_window[, local_idx]
  }

  expect_equal(
    slide_imp(
      to_test,
      location = location,
      window_size = window_size,
      flank = TRUE,
      ncp = 2,
      min_window_n = min_window_n,
      scale = TRUE,
      subset = subset,
      seed = 1234,
      .progress = FALSE
    )[,],
    result
  )
})

test_that("slide_imp: on_infeasible = 'error' rethrows slideimp_infeasible", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100

  # force one window to be infeasible
  mat[1:19, 1:10] <- NA
  mat[20, 1:10] <- rnorm(10)

  expect_error(
    suppressMessages(slide_imp(
      mat,
      location = location,
      k = 3,
      window_size = 10,
      overlap_size = 0,
      min_window_n = 5,
      colmax = 0.9,
      on_infeasible = "error",
      .progress = FALSE
    )),
    class = "slideimp_infeasible"
  )
})

test_that("slide_imp: on_infeasible = 'skip' marks windows and retains originals", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
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
      on_infeasible = "skip",
      .progress = FALSE
    )),
    class = "slideimp_fallback"
  )

  # at least one window was skipped
  expect_gt(length(attr(res, "fallback")), 0)
  expect_equal(attr(res, "fallback_action"), "skip")

  # skipped window columns should retain their original (NA) values
  expect_true(all(is.na(res[1:19, 1:10])))

  # downstream windows (that are feasible) should have imputed values
  expect_false(anyNA(res[, 50:60]))
})

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

  # target 15's window was skipped -> col 15 retains original NA
  expect_true(is.na(res[1, 15]) || all(is.na(res[which(is.na(mat[, 15])), 15])))
  # target 60's window was feasible -> col 60 is imputed
  expect_false(anyNA(res[, 60]))
})

test_that("slide_imp: overlapping windows - skip decrements overlap counts correctly", {
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
  # the skip-decrement was never reached.
  expect_identical(attr(res, "fallback"), 1L)

  # columns 11:20 sit under [1, 20] (skipped) and [11, 30] (feasible). After the
  # decrement their count is 1, so they must equal window [11, 30] alone, with
  # no division. A count left at 2 halves every value here.
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

test_that("slide_imp: all windows infeasible under 'error' fails with slideimp_infeasible", {
  set.seed(1234)
  mat <- matrix(NA_real_, nrow = 20, ncol = 100)
  colnames(mat) <- sample.int(100, size = 100)
  # Leave a handful of values so check_finite passes per column, but colmax rejects
  for (j in seq_len(ncol(mat))) {
    mat[1, j] <- rnorm(1)
  }
  location <- 1:100

  expect_error(
    suppressMessages(slide_imp(
      mat,
      location = location,
      k = 3,
      window_size = 10,
      overlap_size = 0,
      min_window_n = 5,
      colmax = 0.5,
      on_infeasible = "error",
      .progress = FALSE
    )),
    class = "slideimp_infeasible"
  )
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

  # all windows skipped -> every value either original or 0-filled then overwritten
  # subset cols should equal obj (NAs preserved)
  expect_equal(sum(is.na(res)), sum(is.na(mat)))
  expect_length(attr(res, "fallback"), length(attr(res, "fallback"))) # sanity
  expect_equal(attr(res, "fallback_action"), "skip")
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

  # an Inf anywhere under a window is refused
  inf_mat <- mat
  inf_mat[3, 12] <- Inf
  expect_error(
    suppressMessages(do.call(slide_imp, c(list(inf_mat), args))),
    "Infinite"
  )

  # so is an all-NA column, including one sitting in an overlap region
  na_mat <- mat
  na_mat[, 8] <- NA
  expect_error(
    suppressMessages(do.call(slide_imp, c(list(na_mat), args))),
    "All NA/NaN"
  )
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
    "All NA/NaN"
  )

  # Inf is refused wherever it sits, and the reported position refers to the
  # caller's matrix rather than to the window slice
  uncovered_inf <- mat
  uncovered_inf[3, 5] <- Inf
  expect_error(
    suppressMessages(do.call(slide_imp, c(list(uncovered_inf), args))),
    "row 3, column 5"
  )

  covered_inf <- mat
  covered_inf[4, 25] <- Inf
  expect_error(
    suppressMessages(do.call(slide_imp, c(list(covered_inf), args))),
    "row 4, column 25"
  )
})

test_that("count_window_cover counts every window, duplicate ends included", {
  # windows 1-2 and 3-4 end on the same column
  start <- c(1L, 6L, 11L, 16L, 21L)
  end <- c(10L, 10L, 20L, 20L, 30L)
  expect_gt(anyDuplicated(end), 0L)

  expected <- integer(32L)
  for (i in seq_along(start)) {
    cols <- start[i]:end[i]
    expected[cols] <- expected[cols] + 1L
  }

  # Regression: the counts came from `delta[end + 1] <- delta[end + 1] - 1`,
  # which applies a repeated index once. Two windows ending on column 10
  # lowered the running count by 1 instead of 2, so every later column was
  # counted one too high. slide_imp() then halved singly covered columns and
  # left zeros in columns it should have restored.
  expect_identical(count_window_cover(start, end, 32L), expected)
  # columns past the last window are covered by nothing
  expect_identical(count_window_cover(start, end, 32L)[31:32], c(0L, 0L))
  # no windows at all
  expect_identical(count_window_cover(integer(0), integer(0), 4L), integer(4L))
})

test_that("find_windows emits no window nested inside the one before it", {
  # five clusters of 10 columns, 0.1 apart, each followed by a gap of 0.6
  loc <- unlist(lapply(0:4, function(g) g * 1.5 + seq(0, 0.9, by = 0.1)))
  w <- find_windows(loc, window = 1, overlap = 0.5)

  # Regression: the start found in the trailing overlap region was used even
  # when the next coordinate was a full window away, which produced [1, 10]
  # followed by [6, 10]: a window that adds no column and shares its end.
  expect_identical(w$start, c(1L, 11L, 21L, 31L, 41L))
  expect_identical(w$end, c(10L, 20L, 30L, 40L, 50L))

  # the property on irregular coordinates, at several overlaps
  set.seed(1)
  gappy <- sort(c(runif(200, 0, 10), runif(200, 30, 40), runif(50, 100, 101)))
  for (ov in c(0, 1, 2.5, 4.9)) {
    w <- find_windows(gappy, window = 5, overlap = ov)
    expect_true(all(diff(w$start) > 0))
    expect_true(all(diff(w$end) > 0))
    # every column is still inside some window
    covered <- sort(unique(sequence(w$end - w$start + 1L, from = w$start)))
    expect_identical(covered, seq_along(gappy))
    # and a positive overlap still produces overlapping windows, so the fix
    # did not turn every window disjoint
    if (ov > 0) {
      expect_true(any(w$start[-1L] <= w$end[-length(w$end)]))
    }
  }
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
  expect_identical(res[, 15], y[, 15])
  # no window was probed, and the attributes keep one type either way
  expect_identical(attr(res, "solver_chosen"), NA_character_)
  expect_identical(attr(res, "solver_lock_window"), NA_integer_)
})

test_that("slide_imp keeps its return shape when subset resolves to nothing", {
  set.seed(1234)
  x <- sim_mat(20, 30, perc_total_na = 0.2)$input
  colnames(x) <- paste0("c", seq_len(ncol(x)))
  args <- list(
    x,
    location = seq_len(30),
    window_size = 10,
    min_window_n = 5,
    subset = "not_a_column",
    .progress = FALSE
  )

  # Regression: both of these returned the bare input matrix.
  dry <- suppressMessages(do.call(slide_imp, c(args, list(dry_run = TRUE))))
  expect_s3_class(dry, "slideimp_tbl")
  expect_identical(nrow(dry), 0L)
  expect_named(dry, c("start", "end", "window_n", "subset_local"))

  dry_flank <- suppressMessages(do.call(
    slide_imp,
    c(args, list(dry_run = TRUE, flank = TRUE))
  ))
  expect_s3_class(dry_flank, "slideimp_tbl")
  expect_named(
    dry_flank,
    c("start", "end", "window_n", "target", "subset_local")
  )

  res <- suppressMessages(do.call(slide_imp, c(args, list(k = 3))))
  expect_s3_class(res, "slideimp_results")
  expect_identical(attr(res, "metacaller"), "slide_imp")
  expect_identical(attr(res, "imp_method"), "knn")
  expect_length(attr(res, "fallback"), 0L)
  expect_identical(res[,], x)
  expect_no_error(utils::capture.output(print(res), print(dry)))

  # same attribute set as a run that imputes
  full <- suppressMessages(slide_imp(
    x,
    location = seq_len(30),
    window_size = 10,
    min_window_n = 5,
    k = 3,
    .progress = FALSE
  ))
  expect_setequal(names(attributes(res)), names(attributes(full)))
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

test_that("slide_imp records mean-filled columns and cells by their index in obj", {
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

  # pca_imp() has no `subset` and records every ineligible column of a window.
  # Only the requested ones belong in the result.
  pca <- run(ncp = 2)
  expect_identical(attr(pca, "mean_filled_cols"), 52L)
  pca_sub <- run(ncp = 2, subset = setdiff(1:60, 52L))
  expect_identical(attr(pca_sub, "mean_filled_cols"), integer(0))
  expect_true(anyNA(pca_sub[, 52]))
})
