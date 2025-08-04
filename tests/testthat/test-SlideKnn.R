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
  counts <- final_imputed <- matrix(0, nrow = nrow(to_test), ncol = ncol(to_test), dimnames = dimnames(to_test))
  # 1 to 100 is the first window;
  final_imputed[, 1:100] <- final_imputed[, 1:100] + knn_imp(
    obj = to_test[, 1:100], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE
  )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 91 to 190 is the second window;
  final_imputed[, 91:190] <- final_imputed[, 91:190] + knn_imp(
    obj = to_test[, 91:190], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE
  )
  counts[, 91:190] <- counts[, 91:190] + 1
  # 181 to 280 is the last window
  final_imputed[, 181:280] <- final_imputed[, 181:280] + knn_imp(
    obj = to_test[, 181:280], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE
  )
  counts[, 181:280] <- counts[, 181:280] + 1
  final_imputed <- final_imputed / counts

  # SlideKnn should exactly replicate this result
  expect_identical(
    SlideKnn(
      to_test,
      n_feat = 100,
      n_overlap = 10,
      k = 3,
      rowmax = 0.9,
      colmax = 0.9,
      post_imp = TRUE
    ),
    final_imputed
  )
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
  counts <- final_imputed <- matrix(0, nrow = nrow(to_test), ncol = ncol(to_test), dimnames = dimnames(to_test))
  # 1 to 20 is the first window;
  window_cols <- 1:20
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] + knn_imp(
    obj = to_test[, window_cols], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE, subset = local_subset
  )
  counts[, window_cols] <- counts[, window_cols] + 1
  # 16 to 35 is the second window;
  window_cols <- 16:35
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] + knn_imp(
    obj = to_test[, window_cols], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE, subset = local_subset
  )
  counts[, window_cols] <- counts[, window_cols] + 1
  # 31 to 50 is the last window
  window_cols <- 31:50
  local_subset <- which(window_cols %in% subset)
  final_imputed[, window_cols] <- final_imputed[, window_cols] + knn_imp(
    obj = to_test[, window_cols], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE, subset = local_subset
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
  counts <- final_imputed <- matrix(0, nrow = nrow(to_test), ncol = ncol(to_test), dimnames = dimnames(to_test))
  # 1 to 100 is the first window;
  final_imputed[, 1:100] <- final_imputed[, 1:100] + knn_imp(
    obj = to_test[, 1:100], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE
  )
  counts[, 1:100] <- counts[, 1:100] + 1
  # 101 to 200 is the second window;
  final_imputed[, 101:200] <- final_imputed[, 101:200] + knn_imp(
    obj = to_test[, 101:200], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE
  )
  counts[, 101:200] <- counts[, 101:200] + 1
  # 201 to 300 is the last window
  final_imputed[, 201:300] <- final_imputed[, 201:300] + knn_imp(
    obj = to_test[, 201:300], k = 3, colmax = 0.9, rowmax = 0.9, post_imp = TRUE
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
  ram <- SlideKnn(sim, n_feat = 100, n_overlap = 10, k = 5, cores = 1, post_imp = TRUE)

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
  ram_4 <- SlideKnn(sim, n_feat = 100, n_overlap = 10, k = 5, cores = 4, post_imp = TRUE)

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
  expect_no_error(knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn"))
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
    knn_imp = knn_imp,
    subset = NULL
  )

  # Expect that the all-NA rows remain all NA
  testthat::expect_true(all(is.na(imputed[3, ])) && all(is.na(imputed[7, ])))
  # Expect that the specific NA values are imputed
  testthat::expect_true(!is.na(imputed[2, 2]) && !is.na(imputed[4, 4]))
})

test_that("Exactly replicate `impute::impute.knn`", {
  data("khanmiss1")
  # impute is on bioconductor
  testthat::skip_on_cran()
  if (rlang::is_installed("impute")) {
    r1 <- knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn")
    r2 <- t(impute::impute.knn(khanmiss1, k = 3, rowmax = 1, maxp = nrow(khanmiss1))$data)
    expect_equal(r1, r2)
  } else {
    expect_no_error(knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn"))
  }
})

test_that("`subset` feature of `knn_imp` works", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)
  # Impute just 3 columns
  r1 <- knn_imp(to_test, k = 3, post_imp = FALSE, subset = c(1, 3, 5))
  expect_true(!anyNA(r1[, c(1, 3, 5)]))
  r2 <- knn_imp(to_test, k = 3, post_imp = FALSE, subset = paste0("feat", c(1, 3, 5)))
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
