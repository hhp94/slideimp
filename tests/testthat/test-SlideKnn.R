test_that("impute_knn_brute and impute_knn_mlpack calculates the missing location correctly", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.5, perc_col_NA = 1)$input)
  missing <- unname(which(is.na(to_test), arr.ind = TRUE))
  miss <- matrix(is.na(to_test), nrow = nrow(to_test), ncol = ncol(to_test))
  storage.mode(miss) <- "integer"
  n_col_miss <- colSums(is.na(to_test))

  # For brute
  imputed_index_brute <- impute_knn_brute(
    obj = to_test,
    miss = miss,
    k = 5,
    n_col_miss = n_col_miss,
    method = 0,
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

  imputed_index_brute[is.nan(imputed_index_brute)] <- NA
  imputed_index_mlpack[is.nan(imputed_index_mlpack)] <- NA
  expect_equal(imputed_index_brute[, c(1, 2)], missing)
  expect_equal(imputed_index_mlpack[, c(1, 2)], missing)
})

test_that("impute_knn_brute with n_imp > 1 produces correct bootstrap results", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 50, m = 20, perc_NA = 0.3, perc_col_NA = 1)$input)
  missing <- unname(which(is.na(to_test), arr.ind = TRUE))
  miss <- matrix(is.na(to_test), nrow = nrow(to_test), ncol = ncol(to_test))
  storage.mode(miss) <- "integer"
  n_col_miss <- colSums(is.na(to_test))
  n_imp <- 5
  k <- 5
  # Test with bootstrap
  imputed_bootstrap <- impute_knn_brute(
    obj = to_test,
    miss = miss,
    k = k,
    n_col_miss = n_col_miss,
    method = 0,
    weighted = TRUE, # Should be forced to FALSE when n_imp > 1
    dist_pow = 1,
    n_imp = n_imp,
    seed = 42,
    cores = 1
  )
  # Convert NaN to NA for easier testing
  imputed_bootstrap[is.nan(imputed_bootstrap)] <- NA
  # Correct dimensions
  expect_equal(ncol(imputed_bootstrap), n_imp + 2)
  expect_equal(imputed_bootstrap[, 1:2], missing)
  # All bootstrap columns should contain imputed values
  # (some may be NA if no valid neighbors, but structure should be consistent)
  for (b in 3:(2 + n_imp)) {
    expect_true(is.numeric(imputed_bootstrap[, b]))
  }
  # Bootstrap should produce some variability
  # Count how many missing values have at least 2 different imputed values across bootstraps
  n_with_variability <- 0
  for (i in 1:nrow(imputed_bootstrap)) {
    bootstrap_vals <- imputed_bootstrap[i, 3:(2 + n_imp)]
    bootstrap_vals <- bootstrap_vals[!is.na(bootstrap_vals)]
    if (length(bootstrap_vals) > 1 && length(unique(bootstrap_vals)) > 1) {
      n_with_variability <- n_with_variability + 1
    }
  }
  # At least some missing values should show bootstrap variability
  expect_true(n_with_variability > 0)
})

test_that("impute_knn_brute with n_imp > 1 produces reproducible results with same seed", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 30, m = 15, perc_NA = 0.4, perc_col_NA = 1)$input)
  miss <- matrix(is.na(to_test), nrow = nrow(to_test), ncol = ncol(to_test))
  storage.mode(miss) <- "integer"
  n_col_miss <- colSums(is.na(to_test))
  n_imp <- 3
  seed_val <- 123
  # Run twice with same seed
  result1 <- impute_knn_brute(
    obj = to_test,
    miss = miss,
    k = 5,
    n_col_miss = n_col_miss,
    method = 0,
    weighted = FALSE,
    dist_pow = 1,
    n_imp = n_imp,
    seed = seed_val,
    cores = 1
  )
  result2 <- impute_knn_brute(
    obj = to_test,
    miss = miss,
    k = 5,
    n_col_miss = n_col_miss,
    method = 0,
    weighted = FALSE,
    dist_pow = 1,
    n_imp = n_imp,
    seed = seed_val,
    cores = 1
  )
  # Results should be identical
  expect_equal(result1, result2)
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

  # SlideKnn should exactly replicate this result
  simple_mean <- SlideKnn(
    to_test,
    n_feat = 100,
    n_overlap = 10,
    k = 3,
    rowmax = 0.9,
    colmax = 0.9,
    post_imp = TRUE
  )[[1]]
  expect_identical(simple_mean, final_imputed)

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
  )[[1]]
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
  )[[1]]
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
    )[[1]],
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
    )[[1]],
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
  )[[1]]

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
  )[[1]]
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
  )[[1]]

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
  )[[1]]
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
    method = "euclidean"
  ))

  expect_no_error(knn_imp(
    t(khanmiss1),
    k = 3,
    rowmax = 1,
    method = "manhattan",
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
    impute_knn_brute = impute_knn_brute,
    mean_impute_col = mean_impute_col,
    tree = NULL,
    n_imp = 1,
    n_pmm = 1,
    seed = 42
  )[[1]]

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
    r1 <- knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "euclidean")[[1]]

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
    r1.1 <- knn_imp(to_test, k = 5, method = "euclidean", post_imp = FALSE)[[1]]
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
      method = "euclidean",
      post_imp = FALSE,
      subset = subset_cols
    )[[1]][, subset_cols]

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
      method = "euclidean"
    ))
  }
})

test_that("bigmemory functionality in knn_imp works correctly", {
  # Save current option and set bigmemory.allow.dimnames
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = TRUE)
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)
  skip_if_not(isTRUE(getOption("bigmemory.allow.dimnames")))

  # Create test data
  set.seed(1234)
  test_data <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)

  # Create temporary directory for output
  temp_dir <- withr::local_tempdir()
  output_path <- fs::path(temp_dir, "test_imputed")

  # Test 1: Run with output, expect no error
  expect_no_error({
    result_bigmem <- knn_imp(
      test_data,
      k = 3,
      output = output_path,
      overwrite = TRUE,
      n_imp = 3,
      seed = 42
    )
  })

  # Verify that files were created
  base_no_ext <- fs::path_ext_remove(output_path)
  for (i in 1:3) {
    suffix <- paste0("_boot", i)
    bin_file <- paste0(base_no_ext, suffix, ".bin")
    desc_file <- paste0(base_no_ext, suffix, ".desc")
    expect_true(fs::file_exists(bin_file),
      info = paste("Binary file should exist:", bin_file)
    )
    expect_true(fs::file_exists(desc_file),
      info = paste("Descriptor file should exist:", desc_file)
    )
  }

  # Test 2: Rerun with overwrite = FALSE, expect error
  expect_error(
    {
      knn_imp(
        test_data,
        k = 3,
        output = output_path,
        overwrite = FALSE,
        n_imp = 3,
        seed = 42
      )
    },
    "Output files already exist"
  )

  # Test 3: Rerun with overwrite = TRUE, check files exist
  expect_no_error({
    result_bigmem_overwrite <- knn_imp(
      test_data,
      k = 3,
      output = output_path,
      overwrite = TRUE,
      n_imp = 3,
      seed = 42
    )
  })

  if (.Platform$OS.type == "windows") {
    # Force garbage collection to release memory mapped files
    invisible(gc(verbose = FALSE, full = TRUE))
    Sys.sleep(0.1)
  }

  # Verify files still exist after overwrite
  for (i in 1:3) {
    suffix <- paste0("_boot", i)
    bin_file <- paste0(base_no_ext, suffix, ".bin")
    desc_file <- paste0(base_no_ext, suffix, ".desc")
    expect_true(fs::file_exists(bin_file))
    expect_true(fs::file_exists(desc_file))
  }

  # Test 4: Compare in-memory vs file-backed results
  # Run in-memory version with same seed
  result_memory <- knn_imp(
    test_data,
    k = 3,
    output = NULL, # Memory version
    n_imp = 3,
    seed = 42
  )

  # Compare each bootstrap iteration
  for (i in 1:3) {
    # Get the realized matrix from big.matrix
    bigmem_mat <- result_bigmem_overwrite[[i]][, ]
    mem_mat <- result_memory[[i]]

    # Check values are identical (or very close due to potential floating point differences)
    expect_identical(bigmem_mat, mem_mat,
      info = paste("Values should match for bootstrap", i)
    )

    # Check dimnames are preserved
    expect_equal(dimnames(bigmem_mat), dimnames(mem_mat),
      info = paste("Dimnames should match for bootstrap", i)
    )
  }

  # Test that big.matrix objects are properly created
  for (i in 1:3) {
    expect_s4_class(result_bigmem_overwrite[[i]], "big.matrix")
  }
})

test_that("bigmemory with single bootstrap (n_imp=1) works correctly", {
  # Set option
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = TRUE)
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)

  # Create test data
  set.seed(1234)
  test_data <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)

  # Create temporary directory for output
  temp_dir <- withr::local_tempdir()
  output_path <- fs::path(temp_dir, "test_single_boot")

  # Run with n_imp=1 (should not add _boot suffix)
  result_single <- knn_imp(
    test_data,
    k = 3,
    output = output_path,
    overwrite = TRUE,
    n_imp = 1,
    seed = 42
  )

  # Check that files exist without _boot suffix
  base_no_ext <- fs::path_ext_remove(output_path)
  bin_file <- paste0(base_no_ext, ".bin")
  desc_file <- paste0(base_no_ext, ".desc")

  expect_true(fs::file_exists(bin_file))
  expect_true(fs::file_exists(desc_file))

  # Compare with in-memory version
  result_memory <- knn_imp(
    test_data,
    k = 3,
    output = NULL,
    n_imp = 1,
    seed = 42
  )

  expect_identical(result_single[[1]][, ], result_memory[[1]])
})

test_that("bigmemory with subset parameter works correctly", {
  # Set option
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = TRUE)
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)

  # Create test data
  set.seed(1234)
  test_data <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)

  # Create temporary directory for output
  temp_dir <- withr::local_tempdir()
  output_path <- fs::path(temp_dir, "test_subset")

  # Test with subset using numeric index
  result_bigmem_subset <- knn_imp(
    test_data,
    k = 3,
    output = output_path,
    overwrite = TRUE,
    subset = c(1, 3, 5),
    post_imp = FALSE,
    n_imp = 2,
    seed = 42
  )

  # Compare with in-memory version
  result_memory_subset <- knn_imp(
    test_data,
    k = 3,
    output = NULL,
    subset = c(1, 3, 5),
    post_imp = FALSE,
    n_imp = 2,
    seed = 42
  )

  # Check both bootstrap iterations
  for (i in 1:2) {
    bigmem_mat <- result_bigmem_subset[[i]][, ]
    mem_mat <- result_memory_subset[[i]]

    # Check that only specified columns are imputed
    expect_true(!anyNA(bigmem_mat[, c(1, 3, 5)]))
    expect_equal(is.na(bigmem_mat[, -c(1, 3, 5)]), is.na(test_data[, -c(1, 3, 5)]))

    # Compare with memory version
    expect_identical(bigmem_mat, mem_mat)
  }
})

test_that("bigmemory error handling works correctly", {
  # Create test data
  set.seed(1234)
  test_data <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)

  # Test error when bigmemory.allow.dimnames is FALSE
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = FALSE)
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)

  expect_error(
    {
      knn_imp(
        test_data,
        k = 3,
        output = tempfile()
      )
    },
    "`bigmemory.allow.dimnames` must be TRUE"
  )

  # Reset option for next tests
  options(bigmemory.allow.dimnames = TRUE)

  # Test error with invalid output path
  expect_error(
    {
      knn_imp(
        test_data,
        k = 3,
        output = "/invalid/path/that/does/not/exist/output.bin"
      )
    },
    "doesn't exist or is not readable/writable"
  )
})

test_that("`subset` feature of `knn_imp` works with post_imp = FALSE/TRUE", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)
  # Impute just 3 columns
  ## Check subset using numeric index
  r1 <- knn_imp(to_test, k = 3, post_imp = FALSE, subset = c(1, 3, 5))[[1]]
  expect_true(!anyNA(r1[, c(1, 3, 5)]))
  expect_equal(is.na(r1[, -c(1, 3, 5)]), is.na(to_test[, -c(1, 3, 5)]))
  ## Check subset using character vector
  r2 <- knn_imp(
    to_test,
    k = 3,
    post_imp = FALSE,
    subset = paste0("feat", c(1, 3, 5))
  )[[1]]
  expect_identical(r1, r2)

  # Test with post_imp = TRUE and a column requiring post imputation
  to_test_post <- to_test
  # Column 5 will be colMeans if post_imp is TRUE
  to_test_post[2:nrow(to_test_post), 5] <- NA
  r3 <- knn_imp(to_test_post, k = 3, post_imp = TRUE, subset = c(1, 3, 5))[[1]]
  expect_true(!anyNA(r3[, c(1, 3, 5)]))
  # Expect that only the subset columns are imputed. The rests are untouched
  expect_equal(is.na(r3[, -c(1, 3, 5)]), is.na(to_test_post[, -c(1, 3, 5)]))
  # Verify post_imp on column 5
  col5_mean <- mean(to_test_post[, 5], na.rm = TRUE)
  expect_equal(unname(r3[, 5]), rep(col5_mean, nrow(to_test_post)))
  r4 <- knn_imp(
    to_test_post,
    k = 3,
    post_imp = TRUE,
    subset = paste0("feat", c(1, 3, 5))
  )[[1]]
  expect_identical(r3, r4)
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
  expect_identical(mean_impute_col(to_test), c_manual)

  ## Test subset feature
  c_subset <- to_test
  for (i in c(1, 5, 10)) {
    c_subset[is.na(c_subset[, i]), i] <- mean(c_subset[, i], na.rm = TRUE)
  }
  expect_identical(mean_impute_col(to_test, subset = c(1, 5, 10)), c_subset)
})

test_that("`knn_imp` with n_imp = 2 and subset = first 2 columns works", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1)$input)

  # Test with n_imp = 2 and subset = first 2 columns
  results <- knn_imp(
    to_test,
    k = 5,
    n_imp = 2,
    subset = c(1, 2),
    post_imp = TRUE,
    seed = 123
  )

  # Should return a list of length 2
  expect_equal(length(results), 2)

  # Each result should be a matrix with same dimensions as input
  expect_true(all(sapply(results, function(x) all(dim(x) == dim(to_test)))))

  # Only first 2 columns should be imputed (no NAs)
  expect_true(all(sapply(results, function(x) !anyNA(x[, 1:2]))))

  # Other columns should have same NA pattern as original
  expect_true(all(sapply(results, function(x) {
    identical(is.na(x[, -(1:2)]), is.na(to_test[, -(1:2)]))
  })))

  # Bootstrap results should potentially differ
  # (though with small data and specific seed, they might be identical)
  expect_true(is.list(results) && length(results) == 2)
})

test_that("`SlideKnn` with n_imp = 2 and subset = first 2 columns works", {
  set.seed(1234)
  to_test <- t(sim_mat(n = 100, m = 50, perc_NA = 0.3, perc_col_NA = 1)$input)
  to_test_bm <- bigmemory::as.big.matrix(to_test)
  expect_warning(
    results_in_memory <- SlideKnn(
      to_test,
      n_feat = 30,
      n_overlap = 5,
      k = 5,
      n_imp = 2,
      subset = c(1, 2),
      post_imp = TRUE,
      overwrite = TRUE,
      seed = 123
    )
  )
  # Test with n_imp = 2 and subset = first 2 columns
  results <- SlideKnn(
    to_test_bm,
    n_feat = 30,
    n_overlap = 5,
    k = 5,
    n_imp = 2,
    subset = c(1, 2),
    post_imp = TRUE,
    output = withr::local_tempdir(),
    overwrite = TRUE,
    seed = 123
  )
  for (i in 1:2) {
    expect_identical(unname(results_in_memory[[i]]), results[[i]][, ])
  }
  # Should return a list of length 2
  expect_equal(length(results), 2)

  # Each result should be a matrix with same dimensions as input
  expect_true(all(sapply(results, function(x) all(dim(x) == dim(to_test_bm[, ])))))

  # Only first 2 columns should be imputed (no NAs)
  expect_true(all(sapply(results, function(x) !anyNA(x[, 1:2]))))

  # Other columns should have same NA pattern as original
  expect_true(all(sapply(results, function(x) {
    identical(is.na(x[, ][, -(1:2)]), unname(is.na(to_test[, -(1:2)])))
  })))

  # Results should be proper matrices with correct attributes
  expect_true(all(sapply(results, \(x) is.matrix(x[, ]))))
  expect_true(all(sapply(results, \(x) is.numeric(x[, ]))))
})

test_that("`knn_imp` and `SlideKnn` bootstrap reproducibility with seeds", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 10, n = 30, perc_NA = 0.2, perc_col_NA = 1)$input)
  to_test_bm <- bigmemory::as.big.matrix(to_test)

  # Test knn_imp reproducibility
  result1_knn <- knn_imp(to_test, k = 3, n_imp = 2, subset = c(1, 2), seed = 456)
  result2_knn <- knn_imp(to_test, k = 3, n_imp = 2, subset = c(1, 2), seed = 456)

  expect_identical(result1_knn, result2_knn)

  # Test SlideKnn reproducibility
  result1_slide <- SlideKnn(
    to_test_bm,
    n_feat = 15,
    n_overlap = 3,
    k = 3,
    n_imp = 2,
    subset = c(1, 2),
    output = withr::local_tempdir(),
    seed = 456
  )
  result2_slide <- SlideKnn(
    to_test_bm,
    n_feat = 15,
    n_overlap = 3,
    k = 3,
    n_imp = 2,
    subset = c(1, 2),
    output = withr::local_tempdir(),
    seed = 456
  )

  for (i in 1:2) {
    expect_identical(result1_slide[[i]][, ], result2_slide[[i]][, ])
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

test_that("weighted_row_means matches rowMeans with equal weights and all columns", {
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
