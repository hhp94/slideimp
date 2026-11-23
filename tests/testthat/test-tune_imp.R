skip_on_cran()

# resolve_na_loc ----
test_that("resolve_na_loc() rounds near-integer positions", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:90), nrow = 30)

  # Regression: 0.29 * 100 is 28.999999999999996. It passes the integerish
  # check, and as.integer() truncated it to 28: the wrong cell was masked.
  p <- 0.29 * 100
  expect_identical(as.integer(p), 28L)
  expect_identical(unname(rna(obj, p)[[1]][1, ]), c(29L, 1L))
  expect_identical(
    unname(rna(obj, matrix(c(p, 2), ncol = 2))[[1]][1, ]),
    c(29L, 2L)
  )
})

test_that("resolve_na_loc() rejects repeated and missing matrix positions", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:50), nrow = 10)

  # Regression: the matrix form took (1, 1) twice, and tune_imp() scored the
  # cell twice. The vector form already refused the same repeat.
  dup <- matrix(c(1L, 1L, 2L, 1L, 1L, 2L), ncol = 2)
  expect_error(rna(obj, list(dup)), "positions must be unique")
  expect_error(rna(obj, list(c(1L, 1L, 12L))), "duplicated")
  # a repeat that only appears after rounding
  expect_error(rna(obj, list(c(0.29 * 100, 29))), "positions must be unique")
  fill0 <- function(obj) {
    obj[is.na(obj)] <- 0
    obj
  }
  expect_error(
    suppressMessages(tune_imp(
      obj,
      .f = fill0,
      na_loc = list(dup),
      .progress = FALSE
    )),
    "positions must be unique"
  )

  # Regression: an NA inside a matrix read as an out-of-bounds index
  # ("'na_loc[[1]] row indices' failed: Must be TRUE")
  expect_error(rna(obj, list(matrix(c(1, NA, 1, 1), 2))), "missing values")
})

test_that("resolve_na_loc() reads a matrix named col, row by name", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:50), nrow = 10)

  # Regression: cbind(col = 2, row = 5) was read by position, as row 2,
  # column 5.
  expect_identical(
    rna(obj, list(cbind(col = 2L, row = 5L)))[[1]],
    cbind(row = 5L, col = 2L)
  )
  # sample_na_loc() output, and unnamed matrices, are read as before
  expect_identical(
    rna(obj, list(cbind(row = 5L, col = 2L)))[[1]],
    cbind(row = 5L, col = 2L)
  )
  expect_identical(
    rna(obj, list(cbind(5L, 2L)))[[1]],
    cbind(row = 5L, col = 2L)
  )
})

test_that("resolve_na_loc() argument errors are cli conditions", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:50), nrow = 10)

  # Regression: both were base stop() calls
  err <- expect_error(
    rna(obj, "a"),
    "must be a two-column matrix",
    fixed = TRUE
  )
  expect_s3_class(err, "rlang_error")
  err <- expect_error(rna(obj, list(1:2, 1:3)), "same number of NA positions")
  expect_s3_class(err, "rlang_error")
})

# tune imp ----
test_that("tune_imp works", {
  slide_imp_par <- data.frame(
    window_size = c(100, 100),
    k = c(5, 10),
    overlap_size = c(10, 10),
    min_window_n = 20,
    method = "euclidean",
    post_imp = FALSE
  )
  set.seed(1234)
  obj <- sim_mat(50, 1000, perc_col_na = 0.5)$input
  expect_true(anyNA(obj))
  location <- 1:ncol(obj)
  # a failed iteration has numeric columns too, with 0 rows: each run also
  # has no error row and scores every masked cell
  expect_clean_run <- function(res, n) {
    expect_true(all(is.na(res$error)))
    expect_true(all(vapply(res$result, nrow, integer(1)) == n))
  }
  # Check `slide_imp`
  expect_no_error({
    slide_imp_imp_res <- tune_imp(
      obj,
      slide_imp_par,
      .f = "slide_imp",
      location = location,
      n_reps = 1,
      num_na = 200
    )
  })

  # `slide_imp` requires parameters: the window ones are checked first
  expect_error(
    {
      slide_imp_imp_res <- tune_imp(
        obj,
        .f = "slide_imp",
        location = location,
        n_reps = 1,
        num_na = 200
      )
    },
    regexp = "requires `window_size`, `overlap_size`, and `min_window_n`",
    fixed = TRUE
  )

  expect_true(
    all(
      vapply(
        slide_imp_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )
  expect_clean_run(slide_imp_imp_res, 200L)

  # Check `knn_imp`
  knn_imp_par <- data.frame(
    k = c(5, 10),
    method = "euclidean",
    post_imp = TRUE
  )
  expect_no_error({
    knn_imp_res <- tune_imp(
      obj,
      knn_imp_par,
      .f = "knn_imp",
      n_reps = 1,
      num_na = 100
    )
  })

  expect_true(
    all(
      vapply(
        knn_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )
  expect_clean_run(knn_imp_res, 100L)

  # Check `pca_imp`
  pca_imp_par <- data.frame(ncp = 2, miniter = 2)
  expect_no_error({
    pca_imp_res <- tune_imp(
      obj,
      pca_imp_par,
      .f = "pca_imp",
      n_reps = 1,
      num_na = 100
    )
  })

  expect_true(
    all(
      vapply(
        pca_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )
  expect_clean_run(pca_imp_res, 100L)

  # Check custom function
  custom_fun <- function(obj, value) {
    obj[is.na(obj)] <- value
    return(obj)
  }
  custom_par <- data.frame(
    value = c(0, 1)
  )
  expect_no_error({
    custom_imp_res <- tune_imp(
      obj,
      custom_par,
      n_reps = 1,
      num_na = 100,
      .f = custom_fun
    )
  })

  expect_true(
    all(
      vapply(
        custom_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )
  expect_clean_run(custom_imp_res, 100L)
})

test_that("tune_imp works when na_loc is a list of NA locations", {
  # Create a complete matrix (no NAs) for testing
  obj <- sim_mat(50, 200)$input
  obj[is.na(obj)] <- 0 # Fill any existing NAs

  # Create predefined NA location sets
  # Each set has 10 locations, all within matrix bounds
  set.seed(42)
  na_loc_list <- list(
    sample(1:length(obj), 10, replace = FALSE),
    sample(1:length(obj), 10, replace = FALSE),
    sample(1:length(obj), 10, replace = FALSE)
  )

  # Test with slide_imp
  slide_imp_par <- data.frame(
    window_size = 100,
    k = 5,
    overlap_size = 10,
    method = "euclidean",
    min_window_n = 10,
    post_imp = FALSE
  )

  location <- 1:ncol(obj)
  expect_no_error({
    slide_imp_res <- tune_imp(
      location = location,
      obj,
      slide_imp_par,
      .f = "slide_imp",
      na_loc = na_loc_list, # Using list instead of integer
    )
  })

  # Check that we get 3 results (one for each NA location set)
  expect_equal(nrow(slide_imp_res), 3)

  # Check that each result has the correct number of estimates (10 each)
  expect_true(
    all(vapply(slide_imp_res$result, function(x) nrow(x) == 10, logical(1)))
  )

  # Verify the truth values match the original matrix values at those locations
  for (i in 1:3) {
    truth_values <- slide_imp_res$result[[i]]$truth
    expected_truth <- obj[na_loc_list[[i]]]
    expect_equal(truth_values, expected_truth)
  }

  # Test with knn_imp
  knn_imp_par <- data.frame(
    k = c(5, 10),
    method = "euclidean",
    post_imp = FALSE
  )

  expect_no_error({
    knn_imp_res <- tune_imp(
      obj,
      knn_imp_par,
      .f = "knn_imp",
      na_loc = na_loc_list
    )
  })

  # Should have 2 parameters x 3 repetitions = 6 rows
  expect_equal(nrow(knn_imp_res), 6)

  # Check that results contain numeric estimates
  expect_true(
    all(vapply(
      knn_imp_res$result,
      function(x) {
        is.numeric(x$estimate) && is.numeric(x$truth)
      },
      logical(1)
    ))
  )
  # which a failed, 0-row iteration also passes: no error, all 10 scored
  expect_true(all(is.na(knn_imp_res$error)))
  expect_true(all(vapply(knn_imp_res$result, nrow, integer(1)) == 10L))

  # Test with custom function
  custom_fun <- function(obj, value) {
    obj[is.na(obj)] <- value
    return(obj)
  }

  custom_par <- data.frame(value = c(0.5, 1.5))

  expect_no_error({
    custom_res <- tune_imp(
      obj,
      custom_par,
      na_loc = na_loc_list,
      .f = custom_fun
    )
  })

  # Should have 2 parameters x 3 repetitions = 6 rows
  expect_equal(nrow(custom_res), 6)
  # all() below is TRUE over a failed, 0-row iteration: no error, all scored
  expect_true(all(is.na(custom_res$error)))
  expect_true(all(vapply(custom_res$result, nrow, integer(1)) == 10L))

  # Verify custom function fills with the specified values
  for (i in 1:nrow(custom_res)) {
    expected_value <- custom_res$value[i]
    estimates <- custom_res$result[[i]]$estimate
    expect_true(all(estimates == expected_value))
  }

  # Test with another list: 5 positions per element instead of 10
  varied_na_locs <- list(
    sample(1:length(obj), 5, replace = FALSE),
    sample(1:length(obj), 5, replace = FALSE)
  )

  location <- 1:ncol(obj)
  expect_no_error({
    varied_res <- tune_imp(
      obj,
      location = location,
      slide_imp_par,
      .f = "slide_imp",
      na_loc = varied_na_locs
    )
  })

  expect_equal(nrow(varied_res), 2)
  expect_equal(nrow(varied_res$result[[1]]), 5)
  expect_equal(nrow(varied_res$result[[2]]), 5)
})

test_that("tune_imp handles mixed linear and 2D positions in list", {
  set.seed(789)
  obj <- matrix(as.double(1:100), nrow = 10, ncol = 10)

  # mix of linear and 2D positions
  na_locations_mixed <- list(
    c(1, 11, 21), # linear
    matrix(c(10, 10, 10, 1, 2, 3), ncol = 2), # 2D, row 10, column 1, 2, 3
    c(45, 55, 65) # linear
  )

  simple_imp <- function(obj, fill_value) {
    obj[is.na(obj)] <- fill_value
    return(obj)
  }

  params <- data.frame(fill_value = 67)

  result <- tune_imp(
    obj,
    params,
    na_loc = na_locations_mixed,
    .f = simple_imp
  )

  expected_linear <- list(
    c(1, 11, 21),
    c(10, 20, 30),
    c(45, 55, 65)
  )

  for (i in 1:3) {
    res <- result$result[[i]]
    expected_truth <- obj[expected_linear[[i]]]
    expect_equal(res$truth, expected_truth)
    expect_true(all(res$estimate == 67))
  }
})

test_that("compute_metrics works with slideimp_tune and data.frame", {
  set.seed(123)
  obj <- matrix(as.double(1:100), nrow = 10, ncol = 10)

  simple_imp <- function(obj, mu) {
    miss <- is.na(obj)
    obj[miss] <- stats::rnorm(n = sum(miss), mean = mu)
    return(obj)
  }

  params <- data.frame(mu = 42)
  result_tune <- tune_imp(
    obj,
    params,
    n_reps = 2,
    num_na = 10,
    .f = simple_imp
  )

  # slideimp_tune object
  out_tune <- compute_metrics(result_tune)
  expect_s3_class(out_tune, "data.frame")
  expect_true(all(
    c(".metric", ".estimator", ".estimate", "n", "n_miss") %in% names(out_tune)
  ))
  expect_equal(sort(unique(out_tune$.metric)), c("mae", "rmse"))

  # Plain data.frame
  result_df <- as.data.frame(result_tune)
  class(result_df) <- "data.frame"
  out_df <- compute_metrics(result_df)
  expect_s3_class(out_df, "data.frame")
  expect_equal(out_df$.estimate, out_tune$.estimate)
})

test_that("compute_metrics correctly computes n and n_miss with NA estimates", {
  set.seed(456)
  obj <- matrix(as.double(1:100), nrow = 10, ncol = 10)

  simple_imp <- function(obj, mu) {
    miss <- is.na(obj)
    obj[miss] <- rnorm(n = sum(miss), mean = mu)
    return(obj)
  }

  params <- data.frame(mu = 42)
  result <- tune_imp(
    obj,
    params,
    n_reps = 2,
    num_na = 10,
    .f = simple_imp
  )

  # No NAs case: all estimates should be present
  out_clean <- compute_metrics(result)
  # mae and rmse for each of the 2 reps
  expect_identical(nrow(out_clean), 4L)
  expect_true(all(out_clean$n == 10))
  expect_true(all(out_clean$n_miss == 0))

  # Inject NAs into the estimate column of each result element
  result$result[[1]]$estimate[c(1, 3)] <- NA
  result$result[[2]]$estimate[c(2, 5, 7)] <- NA

  out_na <- compute_metrics(
    result,
    metrics = c("mae", "rmse", "mape", "bias", "rsq", "rsq_trad")
  )

  # Rep 1: 10 rows, 2 missing
  rows_rep1 <- out_na[out_na$rep_id == 1, ]
  # all() is TRUE over no rows: one row per metric
  expect_identical(nrow(rows_rep1), 6L)
  expect_true(all(rows_rep1$n == 10))
  expect_true(all(rows_rep1$n_miss == 2))

  # Rep 2: 10 rows, 3 missing
  rows_rep2 <- out_na[out_na$rep_id == 2, ]
  expect_identical(nrow(rows_rep2), 6L)
  expect_true(all(rows_rep2$n == 10))
  expect_true(all(rows_rep2$n_miss == 3))

  # n and n_miss are consistent across metrics within the same rep
  for (r in unique(out_na$rep_id)) {
    subset <- out_na[out_na$rep_id == r, ]
    expect_length(unique(subset$n), 1)
    expect_length(unique(subset$n_miss), 1)
  }
})

test_that("tune_imp works with custom function and list-column parameters", {
  set.seed(42)
  obj <- matrix(rnorm(200), nrow = 10, ncol = 20)

  # custom function that takes a vector of weights per column and fills NAs
  # with a weighted column mean
  weighted_fill <- function(obj, weights) {
    stopifnot(length(weights) == ncol(obj))
    for (j in seq_len(ncol(obj))) {
      col_mean <- mean(obj[, j], na.rm = TRUE)
      obj[is.na(obj[, j]), j] <- col_mean * weights[j]
    }
    return(obj)
  }

  # parameters with a list column: each row holds a different weight vector
  custom_par <- data.frame(
    weights = I(list(
      rep(1, 20),
      rep(0.5, 20),
      seq(0.1, 2, length.out = 20)
    ))
  )

  expect_no_error({
    res <- tune_imp(
      obj,
      custom_par,
      .f = weighted_fill,
      n_reps = 2,
      num_na = 15
    )
  })

  # Should have 3 param sets * 2 reps = 6 rows
  expect_equal(nrow(res), 6)
  expect_true("result" %in% names(res))
  expect_true(
    all(
      vapply(
        res$result,
        \(x) {
          is.numeric(x$estimate) && nrow(x) > 0
        },
        logical(1)
      )
    )
  )
  # The weights list column should be preserved in the output
  expect_true("weights" %in% names(res))
  expect_true(is.list(res$weights))
})

test_that("tune_imp works with custom function and NULL parameters", {
  set.seed(99)
  obj <- matrix(rnorm(150), nrow = 10, ncol = 15)

  # a function with only `obj` - fills NAs with 0
  zero_fill <- function(obj) {
    obj[is.na(obj)] <- 0
    return(obj)
  }

  expect_no_error({
    res <- tune_imp(
      obj,
      parameters = NULL,
      .f = zero_fill,
      n_reps = 3,
      num_na = 10
    )
  })

  # 1 param set * 3 reps = 3 rows
  expect_equal(nrow(res), 3)
  expect_true("result" %in% names(res))
  # placeholder column should be stripped
  expect_false(".placeholder" %in% names(res))
  expect_true(
    all(
      vapply(
        res$result,
        \(x) {
          is.numeric(x$estimate) && is.numeric(x$truth) && nrow(x) > 0
        },
        logical(1)
      )
    )
  )
})

test_that("tune_imp with NULL parameters and a function that has defaults", {
  set.seed(7)
  obj <- matrix(rnorm(100), nrow = 5, ncol = 20)

  fill_with_default <- function(obj, value = -999, scale = 1.0) {
    obj[is.na(obj)] <- value * scale
    return(obj)
  }

  # NULL parameters should run the function using its defaults
  expect_no_error({
    res_null <- tune_imp(
      obj,
      parameters = NULL,
      .f = fill_with_default,
      n_reps = 1,
      num_na = 10
    )
  })

  expect_equal(nrow(res_null), 1)
  # all imputed values should be -999 (the defaults)
  expect_true(all(res_null$result[[1]]$estimate == -999))

  # compare with explicit parameters to make sure NULL truly uses defaults
  explicit_par <- data.frame(value = -999, scale = 1.0)
  expect_no_error({
    res_explicit <- tune_imp(
      obj,
      parameters = explicit_par,
      .f = fill_with_default,
      n_reps = 1,
      num_na = 10
    )
  })

  expect_equal(res_null$result[[1]]$estimate, res_explicit$result[[1]]$estimate)
})

test_that("serial branch reuses `pre` correctly across param_sets within a rep", {
  set.seed(1234)
  obj <- matrix(rnorm(50 * 100), nrow = 50, ncol = 100) # clean, no pre-existing NAs

  seen <- list()
  capture_fun <- function(obj, tag) {
    seen[[length(seen) + 1L]] <<- which(is.na(obj), arr.ind = TRUE)
    obj[is.na(obj)] <- 0
    obj
  }

  params <- data.frame(tag = c("a", "b", "c")) # 3 param sets

  res <- tune_imp(
    obj,
    params,
    .f = capture_fun,
    n_reps = 2,
    num_na = 20,
    .progress = FALSE
  )

  # 1 probe call + 2 reps * 3 param_sets = 7
  expect_length(seen, 7L)

  rep1 <- seen[2:4]
  rep2 <- seen[5:7]

  # within-rep: NA mask identical across param_sets (pre is reused correctly)
  expect_equal(rep1[[1]], rep1[[2]])
  expect_equal(rep1[[2]], rep1[[3]])
  expect_equal(rep2[[1]], rep2[[2]])
  expect_equal(rep2[[2]], rep2[[3]])

  # across-rep: masks differ (pre is rebuilt with the new rep's NAs)
  expect_false(identical(rep1[[1]], rep2[[1]]))

  # truth_list is indexed by the correct rep_id.
  # sort both sides because na_loc ordering (from C++) != which() column-major order.
  for (i in seq_len(nrow(res))) {
    mask <- if (res$rep_id[i] == 1L) rep1[[1]] else rep2[[1]]
    expect_equal(sort(res$result[[i]]$truth), sort(obj[mask]))
  }
})

test_that("slide_imp flank mode is tuned correctly", {
  set.seed(1234)
  ncols <- 100
  obj <- sim_mat(n = 50, p = ncols, perc_col_na = 1)$input
  # we test imputing subset feature 5, 20, and 95
  testing_features <- paste0("feature", c(5, 20, 95))
  na_loc <- sample_na_loc(
    obj = obj,
    n_cols = 3,
    na_col_subset = testing_features,
    n_rows = 2,
    n_reps = 1
  )

  # we manually use slide imp flank mode by manual NA injection into obj
  slide_imp_pre <- obj
  truth <- slide_imp_pre[na_loc[[1]]]
  slide_imp_pre[na_loc[[1]]] <- NA
  location <- seq_len(ncols)

  ## slide_imp on the amputed object. Windows are
  #  # slideimp table: 3 x 5
  # start end window_n target  subset_local
  #     1  10       10      5 <integer [1]>
  #    15  25       11     20 <integer [1]>
  #    90 100       11     95 <integer [1]>

  slide_imp_knn <- slide_imp(
    obj = slide_imp_pre,
    location = location,
    k = 2,
    subset = testing_features,
    min_window_n = 3,
    window_size = 5,
    overlap_size = 0,
    flank = TRUE
  )

  # the PCA fits run a fixed iteration count, miniter = maxiter, so slide_imp()
  # here and under tune_imp() take the same steps; the threshold, far above
  # the criterion at that count, decides only `converged`
  slide_imp_pca <- slide_imp(
    obj = slide_imp_pre,
    location = location,
    ncp = 2,
    threshold = 1e-6,
    maxiter = 30,
    miniter = 30,
    subset = testing_features,
    min_window_n = 3,
    window_size = 5,
    overlap_size = 0,
    flank = TRUE
  )

  ## now we setup tune_imp
  params <- data.frame(
    min_window_n = 3,
    window_size = 5,
    overlap_size = 0,
    flank = TRUE
  )
  knn_params <- params
  knn_params$k <- 2
  pca_params <- params
  pca_params$ncp <- 2
  pca_params$threshold <- 1e-6
  pca_params$maxiter <- 30
  pca_params$miniter <- 30

  tune_imp_knn <- tune_imp(
    obj,
    parameters = knn_params,
    na_loc = na_loc,
    .f = "slide_imp",
    location = location
  )

  tune_imp_pca <- tune_imp(
    obj,
    parameters = pca_params,
    na_loc = na_loc,
    .f = "slide_imp",
    location = location
  )

  # slide_imp() returns only the `subset` columns, 5, 20 and 95, so a masked
  # cell is read at its column's position among them
  loc <- na_loc[[1]]
  expect_identical(ncol(slide_imp_knn), 3L)
  read_at <- cbind(loc[, 1L], match(loc[, 2L], c(5L, 20L, 95L)))
  expect_false(anyNA(read_at))
  expect_equal(
    tune_imp_knn$result[[1]]$estimate,
    slide_imp_knn[read_at]
  )
  expect_equal(
    tune_imp_knn$result[[1]]$truth,
    truth
  )
  # two runs through BLAS, so the values agree to rounding
  expect_worst_err(
    tune_imp_pca$result[[1]]$estimate,
    slide_imp_pca[read_at]
  )
  expect_equal(
    tune_imp_pca$result[[1]]$truth,
    truth
  )
})

test_that("supplied na_loc positions must point at observed cells", {
  set.seed(1234)
  obj <- sim_mat(20, 20, perc_total_na = 0.2)$input
  obj[5, 5] <- NA
  obj[6, 6] <- NaN
  observed <- which(!is.na(obj), arr.ind = TRUE)

  fill0 <- function(obj, dummy = 1) {
    obj[is.na(obj)] <- 0
    obj
  }
  params <- data.frame(dummy = 1)

  # a position sitting on an existing hole has no truth behind it
  on_hole <- rbind(observed[1:3, ], c(5L, 5L))
  expect_error(
    tune_imp(
      obj,
      params,
      .f = fill0,
      na_loc = list(on_hole),
      .progress = FALSE
    ),
    "observed cells"
  )

  # NaN is a hole too
  on_nan <- rbind(observed[1:3, ], c(6L, 6L))
  expect_error(
    tune_imp(obj, params, .f = fill0, na_loc = list(on_nan), .progress = FALSE),
    "observed cells"
  )

  # linear positions go through the same check
  expect_error(
    tune_imp(
      obj,
      params,
      .f = fill0,
      na_loc = list(c(which(is.na(obj))[1], which(!is.na(obj))[1])),
      .progress = FALSE
    ),
    "observed cells"
  )

  # and fully observed positions are scored in full, with no silent shrink
  clean <- observed[1:6, ]
  res <- tune_imp(
    obj,
    params,
    .f = fill0,
    na_loc = list(clean),
    .progress = FALSE
  )
  d <- res$result[[1]]
  expect_equal(nrow(d), 6L)
  expect_false(anyNA(d$truth))
  expect_equal(sum(is.finite(d$truth - d$estimate)), 6L)
})

# tune_imp() driver ----
test_that("a grid built with expand.grid() tunes a built-in method", {
  # expand.grid() makes a character column a factor, and match.arg() in the
  # methods rejects a factor: every iteration used to fail
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  grid <- expand.grid(k = c(3, 5), method = c("euclidean", "manhattan"))
  expect_true(is.factor(grid$method))

  expect_no_warning(
    res <- suppressMessages(
      tune_imp(obj, grid, .f = "knn_imp", num_na = 10, .progress = FALSE)
    )
  )
  expect_identical(nrow(res), 4L)
  expect_true(all(is.na(res$error)))
  expect_true(all(vapply(res$result, nrow, integer(1)) == 10L))
  # the caller's column comes back as it was given
  expect_true(is.factor(res$method))
})

test_that("the output of a built-in method holds only the caller's parameters", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  res <- suppressMessages(
    tune_imp(
      obj,
      data.frame(k = 3),
      .f = "knn_imp",
      num_na = 10,
      .progress = FALSE
    )
  )
  expect_identical(names(res), c("k", "param_set", "rep_id", "result", "error"))

  # a `.progress` column is dropped, with a note, and the call still runs
  expect_message(
    res2 <- tune_imp(
      obj,
      data.frame(k = c(3, 3), .progress = c(TRUE, FALSE)),
      .f = "knn_imp",
      num_na = 10,
      .progress = FALSE
    ),
    "is ignored",
    fixed = TRUE
  )
  expect_identical(names(res2), names(res))
  expect_identical(nrow(res2), 1L)
  expect_true(all(is.na(res2$error)))
})

test_that("a misplaced k or ncp gets the message that names the other method", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  run <- function(parameters, .f) {
    suppressMessages(
      tune_imp(obj, parameters, .f = .f, num_na = 10, .progress = FALSE)
    )
  }
  expect_error(
    run(data.frame(k = 3, ncp = 2), "knn_imp"),
    "To tune .*pca_imp"
  )
  expect_error(
    run(data.frame(ncp = 2, k = 3), "pca_imp"),
    "To tune .*knn_imp"
  )
  # any other unknown name still gets the general message
  expect_error(
    run(data.frame(k = 3, subset = 1), "knn_imp"),
    "unknown parameter"
  )
  expect_error(
    run(data.frame(ncp = 2, subset = 2L), "pca_imp"),
    "unknown parameter"
  )
  expect_error(
    run(data.frame(k = 3, bogus = 1), "slide_imp"),
    "1 unknown parameter in",
    fixed = TRUE
  )
})

test_that("every tunable slide_imp() argument is accepted in parameters", {
  # what tune_imp() supplies itself, or what breaks the matrix return
  own <- c("obj", "location", "subset", "cores", "dry_run")
  tunable <- setdiff(names(formals(slide_imp)), own)
  expect_gt(length(tunable), 5L)
  expect_setequal(slide_tune_params(), tunable)
  expect_true(all(.param_registry$knn$allowed %in% slide_tune_params()))
  expect_true(all(.param_registry$pca$allowed %in% slide_tune_params()))
})

test_that("parameters cannot use a name the tuning output needs", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  f_n <- function(obj, n) {
    obj[is.na(obj)] <- n
    obj
  }
  err <- expect_error(
    tune_imp(
      obj,
      data.frame(n = c(3, 7)),
      .f = f_n,
      num_na = 10,
      .progress = FALSE
    ),
    "Name in use: \"n\""
  )
  expect_s3_class(err, "rlang_error")

  f_two <- function(obj, error, result) {
    obj[is.na(obj)] <- 0
    obj
  }
  expect_error(
    tune_imp(
      obj,
      data.frame(error = 1, result = 2),
      .f = f_two,
      num_na = 10,
      .progress = FALSE
    ),
    "Names in use: \"error\" and \"result\""
  )
})

test_that("a custom .f return is checked on every parameter set", {
  # the probe runs parameter set 1 only. A later set that returned another
  # shape was scored: a matrix index still finds cells in a taller matrix.
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  f_shape <- function(obj, mode) {
    obj[is.na(obj)] <- if (mode == 4) Inf else 0.5
    if (mode == 2) {
      return(rbind(obj, obj))
    }
    if (mode == 3) {
      return(as.vector(obj))
    }
    obj
  }
  expect_warning(
    res <- suppressMessages(
      tune_imp(
        obj,
        data.frame(mode = c(1, 2, 3, 4)),
        .f = f_shape,
        num_na = 10,
        .progress = FALSE
      )
    ),
    "tuning iteration",
    fixed = TRUE
  )
  expect_identical(res$mode, c(1, 2, 3, 4))
  expect_true(is.na(res$error[1]))
  expect_identical(nrow(res$result[[1]]), 10L)
  expect_match(
    res$error[2],
    "returned a 40 x 30 matrix, but 'obj' is 20 x 30",
    fixed = TRUE
  )
  expect_match(res$error[3], "did not return a numeric matrix")
  expect_match(res$error[4], "returned Inf for a masked cell", fixed = TRUE)
  expect_identical(vapply(res$result[2:4], nrow, integer(1)), c(0L, 0L, 0L))
  # the metrics do not stop on the failed, 0-row iterations: they count no
  # cell there (mae and rmse for each set)
  cm <- compute_metrics(res)
  expect_identical(cm$n[cm$mode == 1], c(10L, 10L))
  expect_identical(cm$n[cm$mode != 1], rep(0L, 6L))
})

test_that("the custom .f probe checks Inf in the masked cells only", {
  # Regression: the probe refused Inf anywhere in the return, while every
  # later parameter set refuses it in the masked cells only, so the same .f
  # failed on set 1 and passed on set 2.
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  orig_na <- which(is.na(obj))
  expect_gt(length(orig_na), 0L)
  f_outside <- function(obj, mode) {
    obj[is.na(obj)] <- 0.5
    obj[orig_na[1]] <- Inf
    obj
  }
  res <- suppressMessages(
    tune_imp(
      obj,
      data.frame(mode = 1),
      .f = f_outside,
      num_na = 10,
      .progress = FALSE
    )
  )
  expect_true(is.na(res$error[1]))
  expect_identical(nrow(res$result[[1]]), 10L)
  f_inside <- function(obj, mode) {
    obj[is.na(obj)] <- Inf
    obj
  }
  expect_error(
    suppressMessages(
      tune_imp(
        obj,
        data.frame(mode = 1),
        .f = f_inside,
        num_na = 10,
        .progress = FALSE
      )
    ),
    "failed on the first parameter set",
    fixed = TRUE
  )
})

test_that("cores and metrics are validated by name", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  expect_error(
    suppressMessages(tune_imp(
      obj,
      data.frame(k = 3),
      .f = "knn_imp",
      num_na = 10,
      cores = NA,
      .progress = FALSE
    )),
    "'cores'"
  )
  res <- suppressMessages(
    tune_imp(
      obj,
      data.frame(k = 3),
      .f = "knn_imp",
      num_na = 10,
      .progress = FALSE
    )
  )
  expect_error(compute_metrics(res, metrics = character(0)), "'metrics'")
  expect_error(compute_metrics(res, metrics = c("mae", NA)), "'metrics'")
})

test_that("compute_metrics() errors are cli conditions", {
  err <- expect_error(compute_metrics(data.frame(x = 1)), "result")
  expect_s3_class(err, "rlang_error")

  bad_result <- data.frame(result = I(list(data.frame(a = 1))))
  err <- expect_error(compute_metrics(bad_result), "truth.*estimate")
  expect_s3_class(err, "rlang_error")

  ok <- data.frame(result = I(list(data.frame(truth = 1:3, estimate = 3:1))))
  err <- expect_error(
    compute_metrics(ok, metrics = c("mae", "nope", "zzz")),
    "Unknown metrics in .metrics.: nope, zzz"
  )
  expect_s3_class(err, "rlang_error")
})

test_that("mape does not drop a cell whose true value is 0", {
  # true value 0, estimate not 0: the ratio is Inf
  expect_identical(calc_mape(c(0, 1), c(0.5, 1)), Inf)
  # true value 0, estimate 0: 0 / 0. Dropping that cell gave 50, the error of
  # the other cell alone.
  expect_true(is.nan(calc_mape(c(0, 2), c(0, 1))))
  # a cell with no estimate is still dropped
  expect_equal(calc_mape(c(2, 4, 1), c(1, NA, 1)), 25)
  # no zero: as before
  expect_equal(calc_mape(c(2, 4), c(1, 5)), 37.5)
})

test_that("rsq is NA, without a warning, against a constant", {
  expect_no_warning(const_est <- calc_rsq(c(0.1, 0.5, 0.9), c(0.3, 0.3, 0.3)))
  expect_identical(const_est, NA_real_)
  expect_no_warning(const_truth <- calc_rsq(c(0.3, 0.3, 0.3), c(0.1, 0.5, 0.9)))
  expect_identical(const_truth, NA_real_)
  # and the usual value otherwise
  expect_equal(
    calc_rsq(c(1, 2, 3, 5), c(2, 4, 7, 8)),
    cor(c(1, 2, 3, 5), c(2, 4, 7, 8))^2
  )
})

test_that("tune_imp() reads knn_imp() estimates at the masked cells", {
  # knn_imp() returns only the rep's `subset` columns, so each estimate is
  # read at its column's position in that subset. Compared with a direct
  # call, in one process, cell by cell.
  set.seed(1)
  obj <- sim_mat(40, 30, perc_total_na = 0.05, perc_col_na = 0.5)$input
  na_loc <- sample_na_loc(obj, num_na = 20, n_reps = 2)
  # non-vacuity: the masked columns are not a prefix of obj's, so reading
  # na_loc as is would hit other columns or fall outside the result
  for (m in na_loc) {
    sub <- sort(unique(m[, 2]))
    expect_false(identical(as.integer(sub), seq_along(sub)))
  }
  res <- suppressMessages(tune_imp(
    obj,
    data.frame(k = 3),
    .f = "knn_imp",
    na_loc = na_loc,
    .progress = FALSE
  ))
  expect_true(all(is.na(res$error)))
  for (r in seq_along(na_loc)) {
    m <- na_loc[[r]]
    pre <- obj
    pre[m] <- NA
    sub <- sort(unique(as.integer(m[, 2])))
    direct <- suppressMessages(knn_imp(pre, k = 3, subset = sub))
    want <- direct[cbind(m[, 1], match(m[, 2], sub))]
    expect_length(want, 20L)
    expect_identical(res$result[[r]]$estimate, want, info = r)
    expect_identical(res$result[[r]]$truth, obj[m], info = r)
  }
})

test_that("tune_imp fails an iteration whose slide_imp result is not the subset's width", {
  # Regression: slide_imp() returns only the rep's `subset` columns and
  # tune_imp() reads its estimates there, so a result holding every column,
  # as a daemon running an older slideimp returns, was read at the wrong
  # cells with no error. It now fails the iteration.
  set.seed(1234)
  obj <- sim_mat(n = 30, p = 40, perc_col_na = 1)$input
  na_loc <- sample_na_loc(obj, n_cols = 3, n_rows = 2, n_reps = 1)
  # non-vacuity: the subset is narrower than the matrix
  expect_lt(length(unique(na_loc[[1]][, 2L])), ncol(obj))
  local_mocked_bindings(slide_imp = function(obj, ...) obj)
  suppressMessages(expect_warning(
    res <- tune_imp(
      obj,
      parameters = data.frame(
        k = 2,
        window_size = 10,
        overlap_size = 0,
        min_window_n = 5
      ),
      na_loc = na_loc,
      .f = "slide_imp",
      location = seq_len(ncol(obj)),
      .progress = FALSE
    ),
    "failed"
  ))
  expect_match(
    res$error[[1]],
    "returned 40 columns, but the repetition's subset has 3"
  )
  expect_identical(nrow(res$result[[1]]), 0L)
})
