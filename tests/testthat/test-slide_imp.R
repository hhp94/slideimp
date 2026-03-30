test_that("`slide_imp` knn mode works", {
  set.seed(1234)
  ## Manual minimal implementation to test slide_imp functionality by using
  ## knn_imp, which we test correctness elsewhere
  to_test <- t(
    sim_mat(
      n = 280,
      m = 100,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
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
  to_test <- t(
    sim_mat(
      n = 50,
      m = 10,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
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
  to_test <- t(
    sim_mat(
      n = 300,
      m = 100,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
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
  to_test <- t(
    sim_mat(
      n = 280,
      m = 100,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
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

test_that("`slide_imp` errors on zero-variance features in PCA mode", {
  set.seed(1234)
  to_test <- t(
    sim_mat(
      n = 200,
      m = 10,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
  to_test[, 1] <- 1
  location <- 1:ncol(to_test)
  expect_error(
    slide_imp(
      to_test,
      location = location,
      window_size = 100,
      overlap_size = 10,
      min_window_n = 10,
      ncp = 2,
      miniter = 2
    ),
    regexp = "Features with zero variance after na.rm not permitted for PCA Imputation"
  )
})

test_that("`slide_imp` flank works with knn", {
  set.seed(1234)
  to_test <- t(
    sim_mat(
      n = 50,
      m = 10,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
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

test_that("`slide_imp` flank works with pca", {
  set.seed(1234)
  to_test <- t(
    sim_mat(
      n = 50,
      m = 10,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
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

test_that("`compute_windows` returns expected structure", {
  result <- compute_windows(1:100, window_size = 50, overlap_size = 10)
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("start", "end", "window_n", "keep"))
  expect_true(all(result$window_n == (result$end - result$start + 1L)))
})

test_that("`compute_windows` flags small windows with `min_window_n`", {
  result <- compute_windows(1:100, window_size = 50, overlap_size = 10, min_window_n = 60)
  expect_true(all(result$keep == (result$window_n >= 60)))
  expect_true(any(!result$keep))
})

test_that("`compute_windows` produces more windows with overlap", {
  no_overlap <- compute_windows(1:100, window_size = 50, overlap_size = 0)
  with_overlap <- compute_windows(1:100, window_size = 50, overlap_size = 25)
  expect_gte(nrow(with_overlap), nrow(no_overlap))
})

test_that("`compute_windows` rejects invalid inputs", {
  expect_error(compute_windows(1:100, window_size = 50, overlap_size = 50))
  expect_error(compute_windows(1:100, window_size = -1))
  expect_error(compute_windows(c(3, 1, 2), window_size = 50))
})

test_that("`compute_windows` returns correct windows for known input", {
  result <- compute_windows(1:10, window_size = 5, overlap_size = 0)
  expect_equal(result$start, c(1L, 6L))
  expect_equal(result$end, c(5L, 10L))
})

test_that("`compute_windows` flank returns expected structure", {
  result <- compute_windows(
    1:50,
    window_size = 10,
    subset = c(5, 25, 45),
    flank = TRUE
  )
  expect_s3_class(result, "tbl_df")
  expect_named(result, c("start", "end", "target", "subset_local", "window_n", "keep"))
  expect_equal(nrow(result), 3L)
  expect_true(all(result$window_n == (result$end - result$start + 1L)))
})

test_that("`compute_windows` flank produces correct windows for known input", {
  # location = 1:50, window_size = 10, subset = c(5, 25, 45)
  # center=5:  [1, 15]  -> local idx 5
  # center=25: [15, 35] -> local idx 11
  # center=45: [35, 50] -> local idx 11
  result <- compute_windows(
    1:50,
    window_size = 10,
    subset = c(5, 25, 45),
    flank = TRUE
  )
  expect_equal(result$start, c(1L, 15L, 35L))
  expect_equal(result$end, c(15L, 35L, 50L))
  expect_equal(result$target, c(5L, 25L, 45L))
  expect_equal(result$subset_local, c(5L, 11L, 11L))
})

test_that("`compute_windows` flank flags small windows with `min_window_n`", {
  # Use a small window_size so edge targets produce small windows
  # location = 1:20, subset = c(1, 10, 20), window_size = 3
  # center=1:  [1, 4]  -> n=4
  # center=10: [7, 13] -> n=7
  # center=20: [17, 20] -> n=4
  result <- compute_windows(
    1:20,
    window_size = 3,
    subset = c(1, 10, 20),
    flank = TRUE,
    min_window_n = 5
  )
  expect_equal(result$keep, c(FALSE, TRUE, FALSE))
})

test_that("`compute_windows` flank requires subset", {
  expect_error(
    compute_windows(1:50, window_size = 10, flank = TRUE),
    "subset"
  )
})
