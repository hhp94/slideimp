test_that("grid_to_linear correctly converts 2D positions to linear indices", {
  n <- 10
  m <- 10

  pos_2d <- matrix(c(1, 1, 1, 2, 2, 1, 10, 10), ncol = 2, byrow = TRUE)
  pos_1d <- grid_to_linear(pos_2d, n, m)
  sim_dat <- matrix(rnorm(n * m), ncol = n, nrow = m)
  expect_identical(sim_dat[pos_2d], sim_dat[pos_1d])
})

test_that("tune_imp works", {
  data(khanmiss1)
  slide_imp_par <- data.frame(
    window_size = c(100, 100),
    k = c(5, 10),
    overlap_size = c(10, 10),
    min_window_n = 20,
    method = "euclidean",
    post_imp = FALSE
  )
  set.seed(1234)
  # Tune `slide_imp` function on a subset of khanmiss1
  obj <- t(khanmiss1)[1:30, sample.int(nrow(khanmiss1), size = 200)]
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
  data(khanmiss1)

  # Create a complete matrix (no NAs) for testing
  obj <- t(khanmiss1)[1:30, sample.int(nrow(khanmiss1), size = 200)]
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

  # Custom function that takes a vector of weights per column and fills NAs
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
      num_na = 30
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

  # a function with only `obj` — fills NAs with 0
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
      num_na = 20
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
