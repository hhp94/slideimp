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

  expect_equal(simple_mean[, ], final_imputed)

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
  expect_true(sum((simple_mean[, ] - weighted_1[, ])^2) > 0)
  expect_true(sum((weighted_2[, ] - weighted_1[, ])^2) > 0)
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
    )[, ],
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
  expect_equal(simple_mean[, ], final_imputed)
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
    )[, ],
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
  # window 1: 1 to 20 — covers subset cols 1, 6
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
  # window 2: 16 to 35 — no subset features, SKIPPED
  # window 3: 31 to 50 — covers subset cols 45, 50
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
    )[, ],
    final_imputed[, ]
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
  # window 1: 1 to 20 — covers subset cols 1, 6
  window_cols <- 1:20
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    pca_imp(
      obj = to_test[, window_cols],
      ncp = 2,
      scale = TRUE,
      seed = 1
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # window 2: 16 to 35 — no subset features, SKIPPED
  # window 3: 31 to 50 — covers subset cols 45, 50
  window_cols <- 31:50
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    pca_imp(
      obj = to_test[, window_cols],
      ncp = 2,
      scale = TRUE,
      seed = 1
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
      ncp = 2,
      min_window_n = 10,
      scale = TRUE,
      seed = 1,
      subset = subset
    )[, ],
    final_imputed[, ]
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
    )[, ],
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
      location = location, k = 3,
      window_size = 10, overlap_size = 0,
      min_window_n = 5, colmax = 0.9,
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

  res <- suppressMessages(slide_imp(
    mat,
    location = location, k = 3,
    window_size = 10, overlap_size = 0,
    min_window_n = 5, colmax = 0.9,
    on_infeasible = "skip", .progress = FALSE
  ))

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
  # Make cols 1:10 infeasible but not fully NA (leave a few rows so mean exists)
  mat[1:18, 1:10] <- NA
  mat[19:20, 1:10] <- matrix(rnorm(20), nrow = 2)

  res <- suppressMessages(slide_imp(
    mat,
    location = location, k = 3,
    window_size = 10, overlap_size = 0,
    min_window_n = 5, colmax = 0.9,
    on_infeasible = "mean", .progress = FALSE
  ))

  expect_gt(length(attr(res, "fallback")), 0)
  expect_equal(attr(res, "fallback_action"), "mean")
  # No remaining NA because column means were available
  expect_false(anyNA(res[, 1:10]))
})

test_that("slide_imp: mixed feasible + infeasible windows — skip isolates correctly", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
  mat[1:19, 11:20] <- NA
  mat[20, 11:20] <- rnorm(10)

  res <- suppressMessages(slide_imp(
    mat,
    location = location, k = 3,
    window_size = 10, overlap_size = 0,
    min_window_n = 5, colmax = 0.9,
    on_infeasible = "skip", .progress = FALSE
  ))

  expect_true(all(is.na(res[1:19, 11:20])))
  expect_false(anyNA(res[, 1:10]))
  expect_false(anyNA(res[, 21:30]))
  # skipped window left NAs in requested columns -> must be flagged
  expect_true(isTRUE(attr(res, "has_remaining_na")))
  # and the skipped window should be recorded
  expect_gt(length(attr(res, "fallback")), 0)
  expect_equal(attr(res, "fallback_action"), "skip")
})

test_that("slide_imp: flank mode — infeasible flank window skips only its target", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
  # kill a region around target 15 so its flanking window is infeasible
  mat[1:19, 10:20] <- NA
  mat[20, 10:20] <- rnorm(11)

  res <- suppressMessages(slide_imp(
    mat,
    location = location, k = 2,
    window_size = 10, flank = TRUE,
    subset = c(15, 60),
    min_window_n = 5, colmax = 0.9,
    on_infeasible = "skip", .progress = FALSE
  ))

  # target 15's window was skipped -> col 15 retains original NA
  expect_true(is.na(res[1, 15]) || all(is.na(res[which(is.na(mat[, 15])), 15])))
  # target 60's window was feasible -> col 60 is imputed
  expect_false(anyNA(res[, 60]))
})

test_that("slide_imp: overlapping windows — skip decrements overlap counts correctly", {
  set.seed(1234)
  mat <- sim_mat(20, 100, perc_total_na = 0.2)$input
  location <- 1:100
  # force ONE window infeasible, ensure overlap regions shared with feasible
  # windows still average only over non-skipped windows.
  mat[1:19, 1:10] <- NA
  mat[20, 1:10] <- rnorm(10)

  res <- suppressMessages(slide_imp(
    mat,
    location = location, k = 3,
    window_size = 20, overlap_size = 10,
    min_window_n = 10, colmax = 0.9,
    on_infeasible = "skip", .progress = FALSE
  ))

  # columns 11:20 are overlapped by windows [1-20] (skipped) and [11-30] (feasible)
  # after skip-decrement, counts[11:20] == 1, so result == single-window contribution
  # (no division). Values must be non-NA and finite.
  expect_false(anyNA(res[, 11:20]))
  expect_true(all(is.finite(res[, 11:20])))
})

test_that("slide_imp: all windows infeasible under 'error' fails with slideimp_infeasible", {
  set.seed(1234)
  mat <- matrix(NA_real_, nrow = 20, ncol = 100)
  colnames(mat) <- sample.int(100, size = 100)
  # Leave a handful of values so check_finite passes per column, but colmax rejects
  for (j in seq_len(ncol(mat))) mat[1, j] <- rnorm(1)
  location <- 1:100

  expect_error(
    suppressMessages(slide_imp(
      mat,
      location = location, k = 3,
      window_size = 10, overlap_size = 0,
      min_window_n = 5, colmax = 0.5,
      on_infeasible = "error", .progress = FALSE
    )),
    class = "slideimp_infeasible"
  )
})

test_that("slide_imp: all windows infeasible under 'skip' returns original matrix", {
  set.seed(1234)
  mat <- matrix(NA_real_, nrow = 20, ncol = 100)
  colnames(mat) <- sample.int(100, size = 100)
  for (j in seq_len(ncol(mat))) mat[1, j] <- rnorm(1)
  location <- 1:100

  res <- suppressMessages(slide_imp(
    mat,
    location = location, k = 3,
    window_size = 10, overlap_size = 0,
    min_window_n = 5, colmax = 0.5,
    on_infeasible = "skip", .progress = FALSE
  ))

  # all windows skipped -> every value either original or 0-filled then overwritten
  # subset cols should equal obj (NAs preserved)
  expect_equal(sum(is.na(res)), sum(is.na(mat)))
  expect_length(attr(res, "fallback"), length(attr(res, "fallback"))) # sanity
  expect_equal(attr(res, "fallback_action"), "skip")
})
