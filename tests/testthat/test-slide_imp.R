test_that("`slide_imp` knn mode works", {
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

  expect_identical(simple_mean[, ], final_imputed)

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

test_that("`slide_imp` subset works", {
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

test_that("`slide_imp` edge case no overlap", {
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

test_that("`slide_imp` pca mode works", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## pca_imp, which we test correctness elsewhere
  to_test <- sim_mat(100, 280, perc_total_na = 0.5, perc_col_na = 1)$input
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
  expect_identical(simple_mean[, ], final_imputed)
})

test_that("`slide_imp` don't errors on zero-variance features in PCA mode", {
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

test_that("`slide_imp` flank works with knn", {
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

test_that("`slide_imp` KNN skips windows not covering any subset features", {
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
  # Window 1: 1 to 20 — covers subset cols 1, 6
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
  # Window 2: 16 to 35 — no subset features, SKIPPED
  # Window 3: 31 to 50 — covers subset cols 45, 50
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
  # Average overlaps, restore originals where uncovered
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

test_that("`slide_imp` PCA skips windows not covering any subset features", {
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
  # Window 1: 1 to 20 — covers subset cols 1, 6
  window_cols <- 1:20
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    pca_imp(
      obj = to_test[, window_cols],
      ncp = 2,
      scale = TRUE,
      seed = 1
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # Window 2: 16 to 35 — no subset features, SKIPPED
  # Window 3: 31 to 50 — covers subset cols 45, 50
  window_cols <- 31:50
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    pca_imp(
      obj = to_test[, window_cols],
      ncp = 2,
      scale = TRUE,
      seed = 1
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  # Average overlaps, restore originals where uncovered
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

test_that("`slide_imp` flank works with pca", {
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
