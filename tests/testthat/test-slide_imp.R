test_that("`slide_imp` in-memory matrix mode works", {
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
      rowmax = 0.9,
      post_imp = TRUE
    )[[1]]
  counts[, 1:100] <- counts[, 1:100] + 1
  # 91 to 190 is the second window;
  final_imputed[, 91:190] <- final_imputed[, 91:190] +
    knn_imp(
      obj = to_test[, 91:190],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )[[1]]
  counts[, 91:190] <- counts[, 91:190] + 1
  # 181 to 280 is the last window
  final_imputed[, 181:280] <- final_imputed[, 181:280] +
    knn_imp(
      obj = to_test[, 181:280],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )[[1]]
  counts[, 181:280] <- counts[, 181:280] + 1
  final_imputed <- final_imputed / counts

  # slide_imp should exactly replicate this result
  simple_mean <- slide_imp(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    rowmax = 0.9,
    colmax = 0.9,
    post_imp = TRUE
  )[[1]]
  expect_equal(simple_mean[, ], final_imputed)

  # slide_imp weighted should be different than simple mean
  weighted_1 <- slide_imp(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    rowmax = 0.9,
    colmax = 0.9,
    post_imp = TRUE,
    weighted = TRUE
  )[[1]]
  weighted_2 <- slide_imp(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    rowmax = 0.9,
    colmax = 0.9,
    post_imp = TRUE,
    weighted = TRUE,
    dist_pow = 2
  )[[1]]
  expect_true(sum((simple_mean[, ] - weighted_1[, ])^2) > 0)
  expect_true(sum((weighted_2[, ] - weighted_1[, ])^2) > 0)
})

test_that("`slide_imp` in-memory subset works", {
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
      rowmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )[[1]]
  counts[, window_cols] <- counts[, window_cols] + 1
  # 16 to 35 is the second window;
  window_cols <- 16:35
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )[[1]]
  counts[, window_cols] <- counts[, window_cols] + 1
  # 31 to 50 is the last window
  window_cols <- 31:50
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] +
    knn_imp(
      obj = to_test[, window_cols],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )[[1]]
  counts[, window_cols] <- counts[, window_cols] + 1
  final_imputed <- final_imputed / counts
  # slide_imp should exactly replicate this result
  expect_equal(
    slide_imp(
      to_test,
      n_feat = 20,
      n_overlap = 5,
      k = 3,
      rowmax = 0.9,
      colmax = 0.9,
      post_imp = TRUE,
      subset = subset
    )[[1]][, ],
    final_imputed[, subset, drop = F]
  )
})

test_that("`slide_imp` in-memory edge case no overlap", {
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
      rowmax = 0.9,
      post_imp = TRUE
    )[[1]]
  counts[, 1:100] <- counts[, 1:100] + 1
  # 101 to 200 is the second window;
  final_imputed[, 101:200] <- final_imputed[, 101:200] +
    knn_imp(
      obj = to_test[, 101:200],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )[[1]]
  counts[, 101:200] <- counts[, 101:200] + 1
  # 201 to 300 is the last window
  final_imputed[, 201:300] <- final_imputed[, 201:300] +
    knn_imp(
      obj = to_test[, 201:300],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )[[1]]
  counts[, 201:300] <- counts[, 201:300] + 1
  final_imputed <- final_imputed / counts
  # slide_imp should exactly replicate this result
  expect_equal(
    slide_imp(
      to_test,
      n_feat = 100,
      n_overlap = 0,
      k = 3,
      rowmax = 0.9,
      colmax = 0.9,
      post_imp = TRUE
    )[[1]][, ],
    final_imputed
  )
})

test_that("`slide_imp` bigmemory matrix mode and parallelization works", {
  set.seed(1234)
  # Simulated data
  sim <- t(
    sim_mat(
      n = 200,
      m = 50,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
  dimnames(sim) <- NULL

  # Create explicit temp dir for sim_bm
  temp_dir <- withr::local_tempdir(pattern = "sim_bm")

  # {bigmemory} version
  sim_bm <- bigmemory::as.big.matrix(
    x = sim,
    type = "double",
    backingfile = "sim_bm.bim",
    descriptorfile = "sim_bm.desc",
    backingpath = temp_dir
  )
  expect_true(bigmemory::is.big.matrix(sim_bm))

  # in-memory version
  ram <- slide_imp(
    sim,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 1,
    post_imp = TRUE
  )[[1]]

  # Quickly check that treed version works
  expect_no_error(
    slide_imp(
      sim,
      n_feat = 100,
      n_overlap = 10,
      k = 5,
      cores = 1,
      post_imp = TRUE,
      tree = "ball",
      method = "euclidean"
    )
  )
  # Explicit temp file path for bm output
  temp_bm <- withr::local_tempfile(pattern = "bm")

  # {bigmemory} version
  bm <- slide_imp(
    sim_bm,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 1,
    post_imp = TRUE,
    output = temp_bm,
    overwrite = TRUE
  )[[1]]
  expect_equal(bm[, ], ram[, ])

  # parallel version
  # in-memory version
  skip_if_not(interactive())
  mirai::daemons(4)
  on.exit(mirai::daemons(0))
  ram_4 <- slide_imp(
    sim,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 4,
    post_imp = TRUE
  )[[1]]

  # Explicit temp file path for bm_4 output
  temp_bm4 <- withr::local_tempfile(pattern = "bm_4")

  # {bigmemory} version
  bm_4 <- slide_imp(
    sim_bm,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 4,
    post_imp = TRUE,
    output = temp_bm4,
    overwrite = TRUE
  )[[1]]
  expect_equal(ram[, ], ram_4[, ])
  expect_equal(bm[, ], bm_4[, ])
  mirai::daemons(0)
})

test_that("`impute_knn` ignore all na rows", {
  set.seed(123)

  mat <- matrix(rnorm(100), nrow = 10, ncol = 10, dimnames = list(1:10, 1:10))
  # Introduce all NA rows
  mat[3, ] <- NA
  mat[7, ] <- NA
  # Introduce some specific NA values that should be imputed
  mat[2, 2] <- NA
  mat[4, 4] <- NA

  # Call impute_knn. For all-NA rows to be ignored, their NA proportion (100%)
  # must be >= colmax. Setting colmax to 1 ensures this.
  imputed <- impute_knn(
    obj = mat,
    k = 5,
    rowmax = 0.5,
    colmax = 1,
    method = "euclidean",
    post_imp = TRUE,
    cores = 1,
    weighted = TRUE,
    dist_pow = 1,
    subset = NULL,
    knn_imp = knn_imp,
    impute_knn_brute = impute_knn_brute,
    mean_imp_col = mean_imp_col,
    tree = NULL,
    n_imp = 1,
    n_pmm = 1,
    seed = 42
  )[[1]]

  # Expect that the all-NA rows remain all NA
  expect_true(all(is.na(imputed[3, ])) && all(is.na(imputed[7, ])))
  # Expect that the specific NA values are imputed
  expect_true(!is.na(imputed[2, 2]) && !is.na(imputed[4, 4]))
})


test_that("`mean_imp_col` works", {
  set.seed(4321)
  to_test <- t(
    sim_mat(
      n = 20,
      m = 10,
      perc_NA = 0.5,
      perc_col_NA = 1
    )$input
  )
  c_manual <- to_test
  r_manual <- to_test
  na_indices <- which(is.na(to_test), arr.ind = TRUE)
  column_means <- colMeans(to_test, na.rm = TRUE)
  row_means <- rowMeans(to_test, na.rm = TRUE)

  c_manual[na_indices] <- column_means[na_indices[, 2]]
  expect_equal(mean_imp_col(to_test), c_manual)

  ## Test subset feature
  c_subset <- to_test
  for (i in c(1, 5, 10)) {
    c_subset[is.na(c_subset[, i]), i] <- mean(c_subset[, i], na.rm = TRUE)
  }
  expect_equal(mean_imp_col(to_test, subset = c(1, 5, 10)), c_subset)
})

test_that("`slide_imp` with n_imp = 2 and subset = first 2 columns works", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 100, m = 50, perc_NA = 0.3, perc_col_NA = 1)$input)
  to_test_bm <- bigmemory::as.big.matrix(to_test)
  expect_warning(
    results_in_memory <- slide_imp(
      to_test,
      n_feat = 30,
      n_overlap = 5,
      k = 5,
      n_imp = 2,
      n_pmm = 2,
      subset = c(1, 2),
      post_imp = TRUE,
      overwrite = TRUE,
      seed = 123
    )
  )
  # Test with n_imp = 2 and subset = first 2 columns
  results <- slide_imp(
    to_test_bm,
    n_feat = 30,
    n_overlap = 5,
    k = 5,
    n_imp = 2,
    n_pmm = 2,
    subset = c(1, 2),
    post_imp = TRUE,
    output = withr::local_tempfile(),
    overwrite = TRUE,
    seed = 123
  )
  for (i in 1:2) {
    expect_equal(results_in_memory[[i]][, ], results[[i]][, ])
  }
  # Should return a list of length 2
  expect_equal(length(results), 2)

  # Each result should is a matrix with same dimensions as input
  expect_true(all(vapply(results, function(x) all(dim(x[, ]) == dim(to_test_bm[, c(1, 2)])), logical(1))))

  # Only first 2 columns should be imputed (no NAs)
  expect_true(all(vapply(results, function(x) !anyNA(x[, 1:2]), logical(1))))

  # Results should be proper matrices with correct attributes
  expect_true(all(vapply(results, \(x) is.matrix(x[, ]), logical(1))))
  expect_true(all(vapply(results, \(x) is.numeric(x[, ]), logical(1))))
})

test_that("`slide_imp` MI reproducibility with seeds", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 10, n = 30, perc_NA = 0.2, perc_col_NA = 1)$input)
  to_test_bm <- bigmemory::as.big.matrix(to_test)

  # Test slide_imp reproducibility
  result1_slide <- slide_imp(
    to_test_bm,
    n_feat = 15,
    n_overlap = 3,
    k = 3,
    n_imp = 2,
    n_pmm = 0,
    subset = c(1, 2),
    output = withr::local_tempfile(),
    seed = 456
  )
  result2_slide <- slide_imp(
    to_test_bm,
    n_feat = 15,
    n_overlap = 3,
    k = 3,
    n_imp = 2,
    n_pmm = 0,
    subset = c(1, 2),
    output = withr::local_tempfile(),
    seed = 456
  )

  for (i in 1:2) {
    expect_equal(result1_slide[[i]][, ], result2_slide[[i]][, ])
  }
})

test_that("`find_knn_brute` returns correct neighbors as manual implementation", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 10, m = 30, perc_NA = 0.5, perc_col_NA = 1)$input)
  miss <- is.na(to_test)
  # Ensure columns 1 and 2 have at least one NA
  expect_true(anyNA(to_test[, 1]))
  expect_true(anyNA(to_test[, 2]))
  n_col_miss <- colSums(miss)
  n_col_name <- colnames(to_test)
  # Test only column 1 and 2
  k <- 3
  n_col_miss[3:length(n_col_miss)] <- 0
  result <- find_knn_brute(
    obj = to_test,
    miss = miss,
    k = k,
    n_col_miss = n_col_miss,
    n_col_name = n_col_name,
    method = 0,
    cores = 1
  )
  r_dist <- as.matrix(dist(t(to_test)))

  # For each of the two columns, compute expected neighbors
  for (col_idx in 1:2) {
    col_name <- paste0("feat", col_idx)

    # Get distances for this column, excluding self (diagonal)
    distances <- r_dist[, col_idx]
    distances[col_idx] <- Inf # Exclude self-distance

    # Find k nearest neighbors
    k_nearest_indices <- order(distances)[1:k]
    k_nearest_distances <- unname(distances[k_nearest_indices])

    # Check that the function returned the correct indices
    expect_equal(result[[col_name]]$indices, k_nearest_indices)

    # Check that the function returned the correct distances
    expect_equal(
      sqrt(result[[col_name]]$distances * nrow(to_test)),
      k_nearest_distances
    )

    # Check that k neighbors were returned
    expect_equal(result[[col_name]]$n_neighbors, k)
  }
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

test_that("PMM works", {
  # Manually use implement pmm in R
  # Simulate a 30*6 matrix and column 0 has missing data. Column 1, 2 are the
  # nearest neighbors. Let's do 5 donors, 2 nearest neighbors
  set.seed(123)
  nrow_obj <- 30
  ncol_obj <- 6
  n_miss <- 5
  target_col <- 1 # Column with missing is column 1 (R based)
  nn_col <- c(2, 3) # nearest neighbor is 2, 3
  weights <- c(1, 0.5)
  l <- nrow_obj * ncol_obj
  # Just make sure generated values are unique for easy debug
  while (TRUE) {
    values <- rnorm(l)
    if (length(values) == length(unique(values))) {
      break
    }
  }
  obj <- matrix(values, ncol = ncol_obj)
  # 4 values that requires pmm imputation
  obj[, target_col][sample.int(nrow_obj, size = n_miss)] <- NA
  miss <- is.na(obj)

  # Calculate the donor pool
  y_hat <- weighted_row_means(
    obj[, nn_col],
    miss = miss[, nn_col],
    seq_along(nn_col) - 1, # All columns of the subset matrix.
    weights # Weights
  )[, 1]
  # miss rows
  miss_rows <- which(miss[, target_col])
  # Set values with missing to NA
  y_hat_pool <- y_hat
  y_hat_pool[miss_rows] <- NA
  # Add names for easy tracking
  names(y_hat_pool) <- seq_along(y_hat_pool)

  n_pmm <- 5
  donor_pool <- lapply(
    miss_rows,
    \(x) {
      miss_hat <- y_hat[x]
      # distance to all other rows
      dist_vect <- abs(y_hat_pool - miss_hat)
      # get the actual imputed value pool from the donors back
      obj[as.numeric(names(sort(dist_vect)[seq_len(n_pmm)])), target_col]
    }
  )

  # Initialize result matrix
  n_imp <- 50
  result <- matrix(0, nrow = n_miss, ncol = 2 + n_imp)

  # Modify the data in place
  impute_column_values_pmm(
    result = result,
    obj = obj,
    miss = miss,
    col_offset = 0, # first column
    target_col_idx = target_col - 1, # C++ index is 0 based
    nn_columns = nn_col - 1,
    nn_weights = weights, # weight vec
    n_imp = n_imp,
    n_pmm = n_pmm,
    seed = 42
  )

  # Simulated result matrix. Where column 2:2+n_imp holds each draw of the pmm
  miss_rows_pred <- asplit(result[, c(3:(2 + n_imp)), drop = F], MARGIN = 1)
  for (i in seq_along(miss_rows_pred)) {
    expect_true(all(miss_rows_pred[[i]] %in% donor_pool[[i]]))
  }
})

test_that("`bigmem_impute_colmeans` works", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)
  to_test_bm <- bigmemory::as.big.matrix(
    to_test,
    backingpath = withr::local_tempdir(),
    backingfile = "temp.bin",
    descriptorfile = "temp.desc",
    type = "double"
  )

  to_test_sub <- mean_imp_col(to_test, subset = c(1, 2))
  bigmem_impute_colmeans(to_test_bm@address, col_indices = c(1, 2), cores = 1)
  expect_equal(to_test_bm[, ], to_test_sub)

  to_test_all <- mean_imp_col(to_test_sub)
  bigmem_impute_colmeans(to_test_bm@address, col_indices = seq_len(ncol(to_test)), cores = 4)

  expect_equal(to_test_bm[, ], to_test_all)
  rm(to_test_bm)
  invisible(gc())
})

test_that("`bigmem_add_windows` works", {
  set.seed(1234)
  right <- matrix(rnorm(20 * 50), nrow = 20, ncol = 50)
  temp_dir <- withr::local_tempdir()
  right_bm <- bigmemory::as.big.matrix(
    right,
    backingpath = temp_dir,
    backingfile = "right.bin",
    descriptorfile = "right.desc",
    type = "double"
  )

  # Left matrix starts with zeros
  left_bm <- bigmemory::as.big.matrix(
    matrix(0, nrow = 20, ncol = 10),
    backingpath = temp_dir,
    backingfile = "left.bin",
    descriptorfile = "left.desc",
    type = "double"
  )

  # Add one window (left cols 1-3 += right cols 1-3)
  bigmem_add_windows(
    left_bm@address, right_bm@address,
    start_l = 1, end_l = 3,
    start_r = 1, end_r = 3
  )
  expected <- matrix(0, nrow = 20, ncol = 10)
  expected[, 1:3] <- expected[, 1:3] + right[, 1:3]
  expect_equal(left_bm[, ], expected)

  # Re add to existing values (test addition, not replacement)
  bigmem_add_windows(
    left_bm@address, right_bm@address,
    start_l = 1, end_l = 3,
    start_r = 1, end_r = 3
  )
  expected[, 1:3] <- expected[, 1:3] + right[, 1:3]
  expect_equal(left_bm[, ], expected)

  # Overlapping windows
  left_bm[, ] <- 0
  # Windows that overlap at column 5
  start_l <- c(4, 5)
  end_l <- c(6, 7)
  start_r <- c(1, 10)
  end_r <- c(3, 12)
  bigmem_add_windows(
    left_bm@address, right_bm@address,
    start_l = start_l, end_l = end_l,
    start_r = start_r, end_r = end_r
  )
  expected <- matrix(0, nrow = 20, ncol = 10)
  expected[, 4:6] <- expected[, 4:6] + right[, 1:3]
  expected[, 5:7] <- expected[, 5:7] + right[, 10:12]
  expect_equal(left_bm[, ], expected)

  rm(right_bm, left_bm)
  invisible(gc())
})

test_that("`bigmem_avg` works", {
  set.seed(1234)
  n_rows <- 5

  # Case 1: Overlaps with constant count >1 in regions (cores=1)
  n_cols <- 8
  to_test <- matrix(1, nrow = n_rows, ncol = n_cols)
  to_test_bm <- bigmemory::as.big.matrix(
    to_test,
    backingpath = withr::local_tempdir(),
    backingfile = "temp1.bin",
    descriptorfile = "temp1.desc",
    type = "double"
  )
  interval_start <- c(1, 3, 6)
  interval_end <- c(4, 7, 8)
  res <- find_overlap_regions(interval_start, interval_end)
  start_vec <- res$region[, "start"]
  end_vec <- res$region[, "end"]
  denom_vec <- res$counts_vec
  expected <- sweep(to_test, MARGIN = 2, STATS = denom_vec, FUN = "/")

  bigmem_avg(to_test_bm@address, start_vec, end_vec, denom_vec, cores = 1)
  expect_equal(to_test_bm[, ], expected)
  rm(to_test_bm)
  invisible(gc())

  # Case 2: No overlaps (empty regions, no changes expected)
  n_cols <- 6
  to_test <- matrix(1, nrow = n_rows, ncol = n_cols)
  to_test_bm <- bigmemory::as.big.matrix(
    to_test,
    backingpath = withr::local_tempdir(),
    backingfile = "temp3.bin",
    descriptorfile = "temp3.desc",
    type = "double"
  )
  interval_start <- c(1, 4)
  interval_end <- c(3, 6)
  res <- find_overlap_regions(interval_start, interval_end)
  start_vec <- res$region[, "start"]
  end_vec <- res$region[, "end"]
  denom_vec <- res$counts_vec
  expected <- sweep(to_test, MARGIN = 2, STATS = denom_vec, FUN = "/")
  bigmem_avg(to_test_bm@address, start_vec, end_vec, denom_vec, cores = 1)
  expect_equal(to_test_bm[, ], expected)
  rm(to_test_bm)
  invisible(gc())
})

test_that("`bigmem_copy` works", {
  set.seed(1234)
  right <- matrix(rnorm(20 * 50), nrow = 20, ncol = 50)
  temp_dir <- withr::local_tempdir()
  right_bm <- bigmemory::as.big.matrix(
    right,
    backingpath = temp_dir,
    backingfile = "right.bin",
    descriptorfile = "right.desc",
    type = "double"
  )
  # Copy just one column
  left_sub <- matrix(0, nrow = 20, ncol = 2)
  left_sub_bm <- bigmemory::as.big.matrix(
    left_sub,
    backingpath = temp_dir,
    backingfile = "left_sub.bin",
    descriptorfile = "left_sub.desc",
    type = "double"
  )
  bigmem_copy(left_sub_bm@address, right_bm@address, col_idx_r = c(1, 2), cores = 1)
  expect_equal(left_sub_bm[, ], right[, c(1, 2)])
  left_all <- matrix(0, nrow = 20, ncol = 50)
  left_all_bm <- bigmemory::as.big.matrix(
    left_all,
    backingpath = temp_dir,
    backingfile = "left_all.bin",
    descriptorfile = "left_all.desc",
    type = "double"
  )
  # Copy all columns, shuffled
  col_idx_all <- sample(seq_len(50))
  expected_all <- right[, col_idx_all]
  bigmem_copy(left_all_bm@address, right_bm@address, col_idx_r = col_idx_all, cores = 4)
  expect_equal(left_all_bm[, ], expected_all)
  rm(right_bm, left_sub_bm, left_all_bm)
  invisible(gc())
})
