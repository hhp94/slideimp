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
  simple_mean <- slide_imp(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    colmax = 0.9,
    post_imp = TRUE
  )
  expect_identical(simple_mean[, ], final_imputed)

  # slide_imp weighted should be different than simple mean
  weighted_1 <- slide_imp(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    colmax = 0.9,
    post_imp = TRUE,
    dist_pow = 1
  )
  weighted_2 <- slide_imp(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
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
  expect_equal(
    slide_imp(
      to_test,
      n_feat = 20,
      n_overlap = 5,
      k = 3,
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
  expect_equal(
    slide_imp(
      to_test,
      n_feat = 100,
      n_overlap = 0,
      k = 3,
      colmax = 0.9,
      post_imp = TRUE
    )[, ],
    final_imputed
  )
})

test_that("`weighted_row_means` works", {
  # Generate test data
  set.seed(123)
  to_test <- t(sim_mat(n = 10, m = 10, perc_NA = 0.5, perc_col_NA = 1)$input)

  # Create miss matrix (1 = missing, 0 = observed)
  miss <- is.na(to_test)

  # All columns with equal weights should match rowMeans
  n_cols <- ncol(to_test)
  nn_columns <- 0:(n_cols - 1) # 0-indexed for C++
  nn_weights <- rep(1, n_cols) # Equal weights

  # All Columns
  r1 <- weighted_row_means(to_test, miss, nn_columns, nn_weights)[, 1]
  e1 <- unname(rowMeans(to_test, na.rm = TRUE))

  expect_equal(r1, e1)

  # Selected Columns
  selected_cols <- c(1, 3, 4) # R indexing
  r2 <- weighted_row_means(to_test, miss, selected_cols - 1, nn_weights)[, 1]
  e2 <- unname(rowMeans(to_test[, selected_cols, drop = FALSE], na.rm = TRUE))

  expect_equal(r2, e2)

  # Weighted all cols
  set.seed(1234)
  r_weights <- runif(n_cols, min = 0.1, max = 2)
  r3 <- weighted_row_means(to_test, miss, nn_columns, r_weights)[, 1]

  # Manual Calculation
  weighted_mat <- sweep(to_test, MARGIN = 2, r_weights, FUN = "*")
  weight_mat <- sweep(!is.na(to_test), MARGIN = 2, r_weights, FUN = "*")
  e3 <- rowSums(weighted_mat, na.rm = TRUE) / rowSums(weight_mat, na.rm = TRUE)

  expect_equal(r3, unname(e3))

  # Weighted selected cols
  r4 <- weighted_row_means(to_test, miss, selected_cols - 1, r_weights[selected_cols])[, 1]
  sel_mat <- to_test[, selected_cols, drop = FALSE]
  weighted_mat_sel <- sweep(sel_mat, MARGIN = 2, r_weights[selected_cols], FUN = "*")
  weight_mat_sel <- sweep(!is.na(sel_mat), MARGIN = 2, r_weights[selected_cols], FUN = "*")
  e4 <- rowSums(weighted_mat_sel, na.rm = TRUE) / rowSums(weight_mat_sel, na.rm = TRUE)

  expect_equal(r4, unname(e4))
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
  simple_mean <- slide_imp(
    to_test,
    n_feat = 100,
    n_overlap = 10,
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
  expect_error(
    slide_imp(
      to_test,
      n_feat = 100,
      n_overlap = 10,
      ncp = 2,
      miniter = 2
    ),
    regexp = "Features with zero variance after na.rm not permitted for PCA Imputation"
  )
})
