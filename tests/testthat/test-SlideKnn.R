test_that("SlideKnn in-memory matrix mode works", {
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

test_that("SlideKnn bigmemory matrix mode and parallelization works", {
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
  # {bigmemory} version
  sim_bm <- bigmemory::as.big.matrix(
    x = sim,
    type = "double",
    backingfile = "sim_bm.bim",
    descriptorfile = "sim_bm.desc",
    backingpath = withr::local_tempdir(pattern = "sim_bm")
  )
  expect_true(bigmemory::is.big.matrix(sim_bm))
  # in-memory version
  ram <- SlideKnn(sim, n_feat = 100, n_overlap = 10, k = 5, cores = 1, post_imp = TRUE)
  # {bigmemory} version
  bm <- SlideKnn(
    sim_bm,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 1,
    post_imp = TRUE,
    output = withr::local_tempfile(pattern = "bm"),
    overwrite = TRUE
  )
  expect_identical(bm[, ], ram)
  # parallel version
  skip_on_cran()
  # in-memory version
  mirai::daemons(4)
  on.exit(mirai::daemons(0))
  ram_4 <- SlideKnn(sim, n_feat = 100, n_overlap = 10, k = 5, cores = 4, post_imp = TRUE)
  # {bigmemory} version
  bm_4 <- SlideKnn(
    sim_bm,
    n_feat = 100,
    n_overlap = 10,
    k = 5,
    cores = 4,
    post_imp = TRUE,
    output = withr::local_tempfile(pattern = "bm_4"),
    overwrite = TRUE
  )
  expect_identical(ram, ram_4)
  expect_identical(bm[, ], bm_4[, ])
  mirai::daemons(0)
})

test_that("knn_imp works", {
  data("khanmiss1")
  expect_no_error(knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn"))
})

test_that("impute_knn ignore all na rows", {
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
  imputed <- impute_knn(obj = mat, k = 5, rowmax = 0.5, colmax = 1, method = "euclidean", cores = 1)

  # Expect that the all-NA rows remain all NA
  testthat::expect_true(all(is.na(imputed[3, ])) && all(is.na(imputed[7, ])))
  # Expect that the specific NA values are imputed
  testthat::expect_true(!is.na(imputed[2, 2]) && !is.na(imputed[4, 4]))
})

test_that("Exactly replicate impute::impute.knn", {
  data("khanmiss1")
  testthat::skip_on_cran()
  if (rlang::is_installed("impute")) {
    r1 <- knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn")
    r2 <- t(impute::impute.knn(khanmiss1, k = 3, rowmax = 1, maxp = nrow(khanmiss1))$data)
    expect_equal(r1, r2)
  } else {
    expect_no_error(knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn"))
  }
})
