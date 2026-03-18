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
