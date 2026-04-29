# sample_na_loc ----
test_that("returns a list of length n_reps, each a 2-col row/col matrix", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 3, n_rows = 2, n_reps = 4)

  expect_type(out, "list")
  expect_length(out, 4)
  for (rep in out) {
    expect_equal(nrow(rep), 3 * 2)
  }
})

# n_cols / n_rows semantics
test_that("each rep uses exactly n_cols distinct columns, each with n_rows NAs", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 3, n_reps = 5)
  for (rep in out) {
    tab <- table(rep[, "col"])
    expect_length(tab, 4) # exactly n_cols distinct columns
    expect_true(all(tab == 3)) # each with exactly n_rows NAs
  }
})

test_that("no duplicated (row, col) pairs within a rep", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 3, n_reps = 5)
  for (rep in out) {
    key <- paste(rep[, "row"], rep[, "col"], sep = ":")
    expect_false(anyDuplicated(key) > 0)
  }
})

# num_na distribution
test_that("num_na divisible by n_rows distributes evenly", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, num_na = 12, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 12)
  tab <- table(rep[, "col"])
  expect_length(tab, 4) # 12 / 3 = 4 columns
  expect_true(all(tab == 3))
})

test_that("num_na not divisible by n_rows bumps last buckets by at least 1", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  # num_na = 13, n_rows = 3 -> n_cols = 4, na_per_col = c(4, 3, 3, 3)
  # (remainder = 1 column gets a +1)
  out <- sample_na_loc(m, num_na = 13, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 13)
  tab <- sort(as.integer(table(rep[, "col"])))
  expect_equal(tab, c(3L, 3L, 3L, 4L))
})

# na_col_subset handling
test_that("numeric na_col_subset restricts columns to the pool", {
  m <- sim_mat(p = 10, perc_total_na = 0)$input
  pool <- c(2L, 4L, 6L, 8L)
  out <- sample_na_loc(m, n_cols = 3, n_rows = 2, na_col_subset = pool, n_reps = 10)
  used <- unique(unlist(lapply(out, function(r) r[, "col"])))
  expect_true(all(used %in% pool))
})

test_that("character na_col_subset resolves via colnames", {
  m <- sim_mat(p = 6, perc_total_na = 0)$input # colnames feature1 ... feature6
  out <- sample_na_loc(m,
    n_cols = 2, n_rows = 1,
    na_col_subset = c("feature2", "feature5"), n_reps = 5
  )
  used <- unique(unlist(lapply(out, function(r) r[, "col"])))
  expect_setequal(used, c(2L, 5L))
})

# zero-variance cols protection
test_that("columns keep >= 2 distinct observed values after injection", {
  m <- sim_mat(20, p = 6, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 5, n_reps = 10)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    touched <- unique(rep[, "col"])
    for (j in touched) {
      obs <- m2[, j][!is.na(m2[, j])]
      expect_gte(length(unique(obs)), 2L)
    }
  }
})

test_that("columns without enough eligible rows (after keeping 2 uniques) are skipped", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  # Column 3 has only 4 observed values with exactly 2 uniques.
  # After keeping 2 (one of each), left_over = 2 < needed = 3 -> must be skipped.
  m[1:4, 3] <- c(1, 1, 2, 2)
  m[5:20, 3] <- NA

  out <- sample_na_loc(
    m,
    n_cols = 4,
    n_rows = 3,
    na_col_subset = 1:5,
    n_reps = 10
  )

  for (rep in out) {
    expect_false(
      3L %in% rep[, "col"],
      label = "column 3 should never be selected (not enough sacrificable rows)"
    )
  }
})

test_that("pre-check aborts on columns with zero variance", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  m[, 3] <- 7 # truly zero variance

  expect_error(
    sample_na_loc(m, n_cols = 4, n_rows = 1, na_col_subset = 1:5),
    "Some columns already have zero"
  )
})

test_that("aborts when requested n_cols exceeds available pool", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  expect_error(
    sample_na_loc(m, n_cols = 10, n_rows = 1, na_col_subset = 1:5),
    "Cannot place"
  )
})

# row / col budget enforcement
test_that("resulting matrix respects colmax", {
  m <- sim_mat(20, 8, perc_total_na = 0)$input
  colmax <- 0.5
  out <- sample_na_loc(m,
    n_cols = 4, n_rows = 5,
    colmax = colmax, n_reps = 10
  )
  cap <- floor(nrow(m) * colmax)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    col_miss <- colSums(is.na(m2))
    expect_true(all(col_miss <= cap))
  }
})

test_that("resulting matrix respects rowmax", {
  m <- sim_mat(30, 10, perc_total_na = 0)$input
  rowmax <- 0.4
  out <- sample_na_loc(m,
    n_cols = 5, n_rows = 4,
    rowmax = rowmax, n_reps = 10
  )
  cap <- floor(ncol(m) * rowmax)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    row_miss <- rowSums(is.na(m2))
    expect_true(all(row_miss <= cap))
  }
})

# rep independence / reproducibility
test_that("reps are independently sampled (not identical)", {
  m <- sim_mat(30, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 3, n_reps = 5)
  # Hash each rep; at least 2 distinct with high probability.
  sigs <- vapply(out, function(r) {
    paste(r[, "row"], r[, "col"], collapse = ",")
  }, character(1))
  expect_gt(length(unique(sigs)), 1)
})

test_that("set.seed makes sample_na_loc reproducible", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  set.seed(42)
  a <- sample_na_loc(m, n_cols = 3, n_rows = 2, n_reps = 3)
  set.seed(42)
  b <- sample_na_loc(m, n_cols = 3, n_rows = 2, n_reps = 3)
  expect_equal(a, b)
})

# failure path
test_that("aborts when budgets make sampling infeasible", {
  # 10 rows allows n_rows=3 under the "keep >=2 values" hard bound,
  # but colmax=0.1 -> floor(10*0.1)=1 NA per column, so needing 3 is impossible.
  m <- sim_mat(10, 6, perc_total_na = 0)$input
  expect_error(
    sample_na_loc(m,
      n_cols = 2, n_rows = 3, colmax = 0.1,
      max_attempts = 3
    ),
    "Failed to sample NA locations"
  )
})

# pre-existing NAs
test_that("sampled positions never collide with pre-existing NAs", {
  set.seed(1)
  m <- sim_mat(30, 10, perc_total_na = 0)$input
  # scatter some NAs, keeping columns healthy
  preset <- cbind(row = c(1, 2, 3, 4, 5), col = c(1, 2, 3, 4, 5))
  m[preset] <- NA

  out <- sample_na_loc(m, n_cols = 5, n_rows = 3, n_reps = 20)
  for (rep in out) {
    # none of the sampled (row, col) pairs should already be NA
    expect_true(all(!is.na(m[rep])))
  }
})

test_that("final colmax / rowmax account for pre-existing NAs", {
  set.seed(2)
  m <- sim_mat(20, 8, perc_total_na = 0)$input
  # preload column 1 with 4 NAs and row 1 with 2 NAs (avoiding overlap with col 1)
  m[1:4, 1] <- NA
  m[1, 3:4] <- NA

  colmax <- 0.5 # cap = floor(20 * 0.5) = 10
  rowmax <- 0.5 # cap = floor(8 * 0.5)  = 4

  out <- sample_na_loc(
    m,
    n_cols = 4, n_rows = 3,
    colmax = colmax, rowmax = rowmax,
    n_reps = 20
  )
  col_cap <- floor(nrow(m) * colmax)
  row_cap <- floor(ncol(m) * rowmax)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    expect_true(all(colSums(is.na(m2)) <= col_cap))
    expect_true(all(rowSums(is.na(m2)) <= row_cap))
  }
})

test_that("columns with exhausted col_room are skipped even if individually healthy", {
  set.seed(3)
  m <- sim_mat(20, 6, perc_total_na = 0)$input
  colmax <- 0.5 # cap = 10
  # column 2 already has 9 NAs -> col_room = 1, so needed = 3 can't fit
  na_rows <- sample.int(20, 9)
  m[na_rows, 2] <- NA

  out <- sample_na_loc(
    m,
    n_cols = 4, n_rows = 3,
    colmax = colmax, n_reps = 20
  )
  for (rep in out) {
    expect_false(2L %in% rep[, "col"])
  }
})

test_that("num_na with remainder > n_cols distributes via larger bumps", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  # num_na = 5, n_rows = 3 -> n_cols = 1, one column takes all 5
  out <- sample_na_loc(m, num_na = 5, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 5)
  tab <- as.integer(table(rep[, "col"]))
  expect_equal(tab, 5L)
})

test_that("num_na = 11 with n_rows = 3 yields c(3, 4, 4)", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, num_na = 11, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 11)
  expect_equal(sort(as.integer(table(rep[, "col"]))), c(3L, 4L, 4L))
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

  # `slide_imp` requires parameters
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
    regexp = "requires"
  )

  expect_true(
    all(
      vapply(
        slide_imp_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) == "numeric"
    )
  )

  # Check `knn_imp`
  knn_imp_par <- data.frame(
    k = c(5, 10),
    method = "euclidean",
    post_imp = TRUE
  )
  expect_no_error({
    knn_imp_res <- tune_imp(obj, knn_imp_par, .f = "knn_imp", n_reps = 1, num_na = 100)
  })

  expect_true(
    all(
      vapply(
        knn_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) == "numeric"
    )
  )

  # Check `pca_imp`
  pca_imp_par <- data.frame(ncp = 2, miniter = 2)
  expect_no_error({
    pca_imp_res <- tune_imp(obj, pca_imp_par, .f = "pca_imp", n_reps = 1, num_na = 100)
  })

  expect_true(
    all(
      vapply(
        pca_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) == "numeric"
    )
  )


  # Check custom function
  f1 <- function() {}
  custom_fun <- function(obj, value) {
    obj[is.na(obj)] <- value
    f1()
    return(obj)
  }
  custom_par <- data.frame(
    value = c(0, 1)
  )
  expect_no_error({
    custom_imp_res <- tune_imp(obj, custom_par, n_reps = 1, num_na = 100, .f = custom_fun)
  })

  expect_true(
    all(
      vapply(custom_imp_res$result, \(x) {
        class(x$estimate)
      }, character(1)) == "numeric"
    )
  )
})

test_that("tune_imp works when n_reps is a list of NA locations", {
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

  # Should have 2 parameters × 3 repetitions = 6 rows
  expect_equal(nrow(knn_imp_res), 6)

  # Check that results contain numeric estimates
  expect_true(
    all(vapply(knn_imp_res$result, function(x) {
      is.numeric(x$estimate) && is.numeric(x$truth)
    }, logical(1)))
  )

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

  # Should have 2 parameters × 3 repetitions = 6 rows
  expect_equal(nrow(custom_res), 6)

  # Verify custom function fills with the specified values
  for (i in 1:nrow(custom_res)) {
    expected_value <- custom_res$value[i]
    estimates <- custom_res$result[[i]]$estimate
    expect_true(all(estimates == expected_value))
  }

  # Test with different length NA location sets
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

test_that("tune_imp correctly uses provided NA locations from list", {
  # Create a simple matrix for easier verification
  set.seed(123)
  obj <- matrix(1:100, nrow = 10, ncol = 10)

  # Define specific NA locations
  na_locations <- list(
    c(1, 11, 21), # First column positions
    c(10, 20, 30), # Last position of first 3 rows
    c(50, 60, 70) # Middle positions
  )

  simple_imp <- function(obj, fill_value) {
    obj[is.na(obj)] <- fill_value
    return(obj)
  }

  params <- data.frame(fill_value = 42)

  result <- tune_imp(
    obj,
    params,
    na_loc = na_locations,
    .f = simple_imp
  )

  # Verify each repetition used the correct NA locations
  for (i in 1:3) {
    res <- result$result[[i]]

    # Check truth values match original matrix at specified locations
    expected_truth <- obj[na_locations[[i]]]
    expect_equal(res$truth, expected_truth)

    # Check all estimates are the fill value
    expect_true(all(res$estimate == 42))

    # Check we have the right number of values
    expect_equal(length(res$truth), length(na_locations[[i]]))
  }
})

test_that("tune_imp handles mixed linear and 2D positions in list", {
  set.seed(789)
  obj <- matrix(1:100, nrow = 10, ncol = 10)

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
  obj <- matrix(1:100, nrow = 10, ncol = 10)

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
  expect_true(all(c(".metric", ".estimator", ".estimate", "n", "n_miss") %in% names(out_tune)))
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
  obj <- matrix(1:100, nrow = 10, ncol = 10)

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
  expect_true(all(rows_rep1$n == 10))
  expect_true(all(rows_rep1$n_miss == 2))

  # Rep 2: 10 rows, 3 missing
  rows_rep2 <- out_na[out_na$rep_id == 2, ]
  expect_true(all(rows_rep2$n == 10))
  expect_true(all(rows_rep2$n_miss == 3))

  # n and n_miss are consistent across metrics within the same rep
  for (r in unique(out_na$rep_id)) {
    subset <- out_na[out_na$rep_id == r, ]
    expect_length(unique(subset$n), 1)
    expect_length(unique(subset$n_miss), 1)
  }
})

test_that("compute_metrics.data.frame errors without required columns", {
  bad_df <- data.frame(x = 1:3)
  expect_error(compute_metrics(bad_df), "result")

  bad_result <- data.frame(result = I(list(data.frame(a = 1, b = 2))))
  expect_error(compute_metrics(bad_result), "truth.*estimate")
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
      vapply(res$result, \(x) {
        is.numeric(x$estimate) && nrow(x) > 0
      }, logical(1))
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
      vapply(res$result, \(x) {
        is.numeric(x$estimate) && is.numeric(x$truth) && nrow(x) > 0
      }, logical(1))
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
    obj, params,
    .f = capture_fun,
    n_reps = 2, num_na = 20, .progress = FALSE
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

  slide_imp_pca <- slide_imp(
    obj = slide_imp_pre,
    location = location,
    ncp = 2,
    subset = testing_features,
    min_window_n = 3,
    window_size = 5,
    overlap_size = 0,
    flank = TRUE
  )

  ## now we setup tune_imp
  params <- data.frame(
    min_window_n = 3, window_size = 5, overlap_size = 0, flank = TRUE
  )
  knn_params <- params
  knn_params$k <- 2
  pca_params <- params
  pca_params$ncp <- 2

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

  expect_equal(
    tune_imp_knn$result[[1]]$estimate,
    slide_imp_knn[na_loc[[1]]]
  )
  expect_equal(
    tune_imp_knn$result[[1]]$truth,
    truth
  )
  expect_equal(
    tune_imp_pca$result[[1]]$estimate,
    slide_imp_pca[na_loc[[1]]]
  )
  expect_equal(
    tune_imp_pca$result[[1]]$truth,
    truth
  )
})
