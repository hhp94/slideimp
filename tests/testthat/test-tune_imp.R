test_that("tune_imp .f = `slide_imp`, `knn_imp`, `pca_imp` or function works", {
  data(khanmiss1)
  slide_imp_par <- data.frame(
    n_feat = c(100, 100),
    k = c(5, 10),
    n_overlap = c(10, 10),
    knn_method = "euclidean",
    post_imp = FALSE
  )
  set.seed(1234)
  # Tune `slide_imp` function on a subset of khanmiss1
  obj <- t(khanmiss1)[1:30, sample.int(nrow(khanmiss1), size = 200)]
  expect_true(anyNA(obj))

  # Check `slide_imp`
  expect_no_error({
    slide_imp_imp_res <- tune_imp(obj, slide_imp_par, rep = 1, num_na = 200, .f = "slide_imp")
  })

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
    knn_imp_res <- tune_imp(obj, knn_imp_par, rep = 1, num_na = 100, .f = "knn_imp")
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
    pca_imp_res <- tune_imp(obj, pca_imp_par, rep = 1, num_na = 100, .f = "pca_imp")
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
    custom_imp_res <- tune_imp(obj, custom_par, rep = 1, num_na = 100, .f = custom_fun)
  })

  expect_true(
    all(
      vapply(custom_imp_res$result, \(x) {
        class(x$estimate)
      }, character(1)) == "numeric"
    )
  )
})

test_that("tune_imp works when rep is a list of NA locations", {
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
    n_feat = 100,
    k = 5,
    n_overlap = 10,
    knn_method = "euclidean",
    post_imp = FALSE
  )

  expect_no_error({
    slide_imp_res <- tune_imp(
      obj,
      slide_imp_par,
      rep = na_loc_list, # Using list instead of integer
      .f = "slide_imp"
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
      rep = na_loc_list,
      .f = "knn_imp"
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
      rep = na_loc_list,
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

  expect_no_error({
    varied_res <- tune_imp(
      obj,
      slide_imp_par,
      rep = varied_na_locs,
      .f = "slide_imp"
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

  # Simple imputation function that replaces NA with a fixed value
  simple_imp <- function(obj, fill_value) {
    obj[is.na(obj)] <- fill_value
    return(obj)
  }

  params <- data.frame(fill_value = 42)

  result <- tune_imp(
    obj,
    params,
    rep = na_locations,
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
