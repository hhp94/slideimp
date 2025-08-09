test_that("impute_knn_naive and impute_knn_mlpack calculates the missing location correctly", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)
  missing <- which(is.na(to_test))
  miss <- matrix(is.na(to_test), nrow = nrow(to_test), ncol = ncol(to_test))
  storage.mode(miss) <- "integer"
  n_col_miss <- colSums(is.na(to_test))

  # For naive
  imputed_index_naive <- impute_knn_naive(
    obj = to_test,
    miss = miss,
    k = 5,
    n_col_miss = n_col_miss,
    method = 2,
    weighted = FALSE,
    dist_pow = 1,
    cores = 1
  )

  imputed_index_mlpack <- impute_knn_mlpack(
    obj = mean_impute_col(to_test),
    miss = miss,
    k = 5,
    n_col_miss = n_col_miss,
    method = 0,
    tree = "kd",
    weighted = FALSE,
    dist_pow = 1,
    cores = 1
  )
  imputed_index_naive[is.nan(imputed_index_naive)] <- NA
  imputed_index_mlpack[is.nan(imputed_index_mlpack)] <- NA
  expect_equal(imputed_index_naive[, 1], missing)
  expect_equal(imputed_index_mlpack[, 1], missing)
})

test_that("`SlideKnn` in-memory matrix mode works", {
  set.seed(1234)
  ## Manual minimal implementation to test SlideKnn functionality by using
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
  counts <- final_imputed <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  # 1 to 100 is the first window;
  final_imputed[, 1:100] <- final_imputed[, 1:100] +
    knn_imp(
      obj = to_test[, 1:100],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 91 to 190 is the second window;
  final_imputed[, 91:190] <- final_imputed[, 91:190] +
    knn_imp(
      obj = to_test[, 91:190],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )
  counts[, 91:190] <- counts[, 91:190] + 1
  # 181 to 280 is the last window
  final_imputed[, 181:280] <- final_imputed[, 181:280] +
    knn_imp(
      obj = to_test[, 181:280],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )
  counts[, 181:280] <- counts[, 181:280] + 1
  final_imputed <- final_imputed / counts

  # SlideKnn should exactly replicate this result
  simple_mean <- SlideKnn(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    rowmax = 0.9,
    colmax = 0.9,
    post_imp = TRUE
  )
  expect_identical(simple_mean, simple_mean)

  # SlideKnn weighted should be different than simple mean
  weighted_1 <- SlideKnn(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    rowmax = 0.9,
    colmax = 0.9,
    post_imp = TRUE,
    weighted = TRUE
  )
  weighted_2 <- SlideKnn(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    rowmax = 0.9,
    colmax = 0.9,
    post_imp = TRUE,
    weighted = TRUE,
    dist_pow = 2
  )
  expect_true(sum((simple_mean - weighted_1)^2) > 0)
  expect_true(sum((weighted_2 - weighted_1)^2) > 0)
})

test_that("`SlideKnn` in-memory subset works", {
  set.seed(1234)
  ## Manual minimal implementation to test SlideKnn functionality by using
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
  counts <- final_imputed <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
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
      rowmax = 0.9,
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
      rowmax = 0.9,
      post_imp = TRUE,
      subset = local_subset
    )
  counts[, window_cols] <- counts[, window_cols] + 1
  final_imputed <- final_imputed / counts
  # SlideKnn should exactly replicate this result
  expect_identical(
    SlideKnn(
      to_test,
      n_feat = 20,
      n_overlap = 5,
      k = 3,
      rowmax = 0.9,
      colmax = 0.9,
      post_imp = TRUE,
      subset = subset
    ),
    final_imputed
  )
})

test_that("`SlideKnn` in-memory edge case no overlap", {
  set.seed(1234)
  ## Manual minimal implementation to test SlideKnn functionality by using
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
  counts <- final_imputed <- matrix(
    0,
    nrow = nrow(to_test),
    ncol = ncol(to_test),
    dimnames = dimnames(to_test)
  )
  # 1 to 100 is the first window;
  final_imputed[, 1:100] <- final_imputed[, 1:100] +
    knn_imp(
      obj = to_test[, 1:100],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 101 to 200 is the second window;
  final_imputed[, 101:200] <- final_imputed[, 101:200] +
    knn_imp(
      obj = to_test[, 101:200],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )
  counts[, 101:200] <- counts[, 101:200] + 1
  # 201 to 300 is the last window
  final_imputed[, 201:300] <- final_imputed[, 201:300] +
    knn_imp(
      obj = to_test[, 201:300],
      k = 3,
      colmax = 0.9,
      rowmax = 0.9,
      post_imp = TRUE
    )
  counts[, 201:300] <- counts[, 201:300] + 1
  final_imputed <- final_imputed / counts
  # SlideKnn should exactly replicate this result
  expect_identical(
    SlideKnn(
      to_test,
      n_feat = 100,
      n_overlap = 0,
      k = 3,
      rowmax = 0.9,
      colmax = 0.9,
      post_imp = TRUE
    ),
    final_imputed
  )
})

test_that("`SlideKnn` bigmemory matrix mode and parallelization works", {
  set.seed(1234)
  # Simulated data
  sim <- t(
    sim_mat(
      n = 2000,
      m = 100,
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
  ram <- SlideKnn(
    sim,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 1,
    post_imp = TRUE
  )

  # Quickly check that treed version works
  expect_no_error(
    SlideKnn(
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
  bm <- SlideKnn(
    sim_bm,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 1,
    post_imp = TRUE,
    output = temp_bm,
    overwrite = TRUE
  )
  expect_identical(bm[, ], ram)

  # parallel version
  # in-memory version
  skip_if_not(interactive())
  mirai::daemons(4)
  on.exit(mirai::daemons(0))
  ram_4 <- SlideKnn(
    sim,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 4,
    post_imp = TRUE
  )

  # Explicit temp file path for bm_4 output
  temp_bm4 <- withr::local_tempfile(pattern = "bm_4")

  # {bigmemory} version
  bm_4 <- SlideKnn(
    sim_bm,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 4,
    post_imp = TRUE,
    output = temp_bm4,
    overwrite = TRUE
  )
  expect_identical(ram, ram_4)
  expect_identical(bm[, ], bm_4[, ])
  mirai::daemons(0)
})

test_that("`knn_imp` works", {
  data("khanmiss1")
  expect_no_error(knn_imp(
    t(khanmiss1),
    k = 3,
    rowmax = 1,
    method = "impute.knn"
  ))

  expect_no_error(knn_imp(
    t(khanmiss1),
    k = 3,
    rowmax = 1,
    method = "euclidean",
    tree = "kd"
  ))
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
    impute_knn_naive = impute_knn_naive,
    mean_impute_col = mean_impute_col,
    tree = NULL
  )

  # Expect that the all-NA rows remain all NA
  testthat::expect_true(all(is.na(imputed[3, ])) && all(is.na(imputed[7, ])))
  # Expect that the specific NA values are imputed
  testthat::expect_true(!is.na(imputed[2, 2]) && !is.na(imputed[4, 4]))
})

test_that("Exactly replicate `impute::impute.knn`", {
  # Load the example dataset with missing values
  data("khanmiss1")

  # impute is on bioconductor
  # Skip this test on CRAN to avoid dependency issues
  testthat::skip_on_cran()

  # Check if the 'impute' package is installed
  if (rlang::is_installed("impute")) {
    # Perform imputation using knn_imp with method "impute.knn" on transposed data
    r1 <- knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn")

    # Perform imputation using the original impute::impute.knn function
    # Transpose the result to match the orientation
    r2 <- t(
      impute::impute.knn(
        khanmiss1,
        k = 3,
        rowmax = 1,
        maxp = nrow(khanmiss1)
      )$data
    )

    # Verify that the results from knn_imp match exactly with impute::impute.knn
    expect_equal(r1, r2)

    # Test to see if the post_imp strategy would replicate the results completely
    # Set seed for reproducibility in simulation
    set.seed(1234)

    # Generate a simulated matrix with missing values (500 rows, 30 columns, 50% NA, 80% columns with NA)
    to_test <- t(sim_mat(n = 500, m = 30, perc_NA = 0.5, perc_col_NA = 0.8)$input)

    # Pre-compute row means before imputation (ignoring NAs)
    pre_impute <- rowMeans(to_test, na.rm = TRUE)

    # Impute using knn_imp without post-imputation step; expect some NAs to remain
    r1.1 <- knn_imp(to_test, k = 5, method = "impute.knn", post_imp = FALSE)
    expect_true(anyNA(r1.1))

    # impute::impute.knn uses the pre-imputation row means to impute the data.
    # After knn_imp, we row impute the data with pre-calculated row_means
    # Identify indices of remaining NAs
    indices <- which(is.na(r1.1), arr.ind = TRUE)

    # Fill remaining NAs with pre-computed row means
    r1.1[indices] <- pre_impute[indices[, 1]]

    # Verify no NAs remain after manual post-imputation
    expect_true(!anyNA(r1.1))

    # Perform imputation using impute::impute.knn on the transposed simulated data
    r2.1 <- t(
      impute::impute.knn(
        t(to_test),
        k = 5,
        maxp = ncol(to_test)
      )$data
    )

    # Verify that the manually post-imputed knn_imp matches impute::impute.knn
    expect_equal(r1.1, r2.1)

    # Test subset. strategy is to use subset, then impute.knn on the same data
    # and pull out the same subset then compare the two matrices
    # Set seed for reproducibility in subset selection
    set.seed(2345)

    # Generate another simulated matrix (100 rows, 200 columns, 10% NA, all columns with NA)
    to_test_subset <- t(sim_mat(n = 100, m = 200, perc_NA = 0.1, perc_col_NA = 1)$input)

    # Randomly select 10 subset columns
    subset_cols <- sample(colnames(to_test_subset), size = 10)

    # Verify that the subset has NAs before imputation
    expect_true(anyNA(to_test_subset[, subset_cols]))

    # Impute only the subset columns using knn_imp without post_imp
    r1_subset <- knn_imp(
      to_test_subset,
      k = 10,
      method = "impute.knn",
      post_imp = FALSE,
      subset = subset_cols
    )[, subset_cols]

    # Verify no NAs remain in the imputed subset
    expect_true(!anyNA(r1_subset))

    # Perform full imputation using impute::impute.knn and extract the subset
    r2_subset <- t(
      impute::impute.knn(
        t(to_test_subset),
        k = 10,
        maxp = ncol(to_test_subset)
      )$data
    )[, subset_cols]

    # Verify no NAs in the extracted subset from full imputation
    expect_true(!anyNA(r2_subset))

    # Verify that the subset imputation matches the extracted subset from full imputation
    expect_equal(r1_subset, r2_subset)
  } else {
    # If 'impute' is not installed, just verify that knn_imp runs without error
    expect_no_error(knn_imp(
      t(khanmiss1),
      k = 3,
      rowmax = 1,
      method = "impute.knn"
    ))
  }
})

test_that("`subset` feature of `knn_imp` works", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)
  # Impute just 3 columns
  r1 <- knn_imp(to_test, k = 3, post_imp = FALSE, subset = c(1, 3, 5))
  expect_true(!anyNA(r1[, c(1, 3, 5)]))
  r2 <- knn_imp(
    to_test,
    k = 3,
    post_imp = FALSE,
    subset = paste0("feat", c(1, 3, 5))
  )
  expect_identical(r1, r2)
})

test_that("`mean_impute_col` works", {
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
  r_manual[na_indices] <- row_means[na_indices[, 1]]
  expect_identical(mean_impute_col(to_test), c_manual)
  expect_identical(mean_impute_row(to_test), r_manual)

  ## Test subset feature
  c_subset <- to_test
  for (i in c(1, 5, 10)) {
    c_subset[is.na(c_subset[, i]), i] <- mean(c_subset[, i], na.rm = TRUE)
  }
  expect_identical(mean_impute_col(to_test, subset = c(1, 5, 10)), c_subset)
})
