test_that("`impute_knn_brute` and impute_knn_mlpack calculates the missing location correctly", {
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
    dist_pow = 1,
    cores = 1
  )

  imputed_index_mlpack <- impute_knn_mlpack(
    obj = mean_imp_col(to_test),
    miss = miss,
    k = 5,
    n_col_miss = n_col_miss,
    method = 0,
    tree = "kd",
    dist_pow = 1,
    cores = 1
  )

  imputed_index_brute[is.nan(imputed_index_brute)] <- NA
  imputed_index_mlpack[is.nan(imputed_index_mlpack)] <- NA
  expect_equal(imputed_index_brute[, c(1, 2)], missing)
  expect_equal(imputed_index_mlpack[, c(1, 2)], missing)
})

test_that("`knn_imp` works", {
  data("khanmiss1")
  expect_no_error(knn_imp(
    t(khanmiss1),
    k = 3,
    method = "euclidean"
  ))

  expect_no_error(knn_imp(
    t(khanmiss1),
    k = 3,
    method = "manhattan",
    tree = "kd"
  ))
})

test_that("Exactly replicate `impute::impute.knn`", {
  # impute is on bioconductor
  # Skip this test on CRAN to avoid dependency issues
  skip("Uncomment and test manually")
  # data("khanmiss1")
  #
  # # Check if the 'impute' package is installed
  #
  # # Perform imputation using knn_imp with method "impute.knn" on transposed data
  # r1 <- knn_imp(t(khanmiss1), k = 3, method = "euclidean")
  #
  # # Perform imputation using the original impute::impute.knn function
  # # Transpose the result to match the orientation
  # r2 <- t(
  #   impute::impute.knn(
  #     khanmiss1,
  #     k = 3,
  #     maxp = nrow(khanmiss1)
  #   )$data
  # )
  #
  # # Verify that the results from knn_imp match exactly with impute::impute.knn
  # expect_equal(r1[, ], r2[, ])
  #
  # # Test to see if the post_imp strategy would replicate the results completely
  # # Set seed for reproducibility in simulation
  # set.seed(1234)
  #
  # # Generate a simulated matrix with missing values (500 rows, 30 columns, 50% NA, 80% columns with NA)
  # to_test <- t(sim_mat(n = 500, m = 30, perc_NA = 0.5, perc_col_NA = 0.8)$input)
  #
  # # Pre-compute row means before imputation (ignoring NAs)
  # pre_impute <- rowMeans(to_test, na.rm = TRUE)
  #
  # # Impute using knn_imp without post-imputation step; expect some NAs to remain
  # r1.1 <- knn_imp(to_test, k = 5, method = "euclidean", post_imp = FALSE)
  # expect_true(anyNA(r1.1))
  #
  # # impute::impute.knn uses the pre-imputation row means to impute the data.
  # # After knn_imp, we row impute the data with pre-calculated row_means
  # # Identify indices of remaining NAs
  # indices <- which(is.na(r1.1), arr.ind = TRUE)
  #
  # # Fill remaining NAs with pre-computed row means
  # r1.1[indices] <- pre_impute[indices[, 1]]
  #
  # # Verify no NAs remain after manual post-imputation
  # expect_true(!anyNA(r1.1))
  #
  # # Perform imputation using impute::impute.knn on the transposed simulated data
  # r2.1 <- t(
  #   impute::impute.knn(
  #     t(to_test),
  #     k = 5,
  #     maxp = ncol(to_test)
  #   )$data
  # )
  #
  # # Verify that the manually post-imputed knn_imp matches impute::impute.knn
  # expect_equal(r1.1[, ], r2.1[, ])
  #
  # # Test subset. strategy is to use subset, then impute.knn on the same data
  # # and pull out the same subset then compare the two matrices
  # # Set seed for reproducibility in subset selection
  # set.seed(2345)
  #
  # # Generate another simulated matrix (100 rows, 200 columns, 10% NA, all columns with NA)
  # to_test_subset <- t(sim_mat(n = 100, m = 200, perc_NA = 0.1, perc_col_NA = 1)$input)
  #
  # # Randomly select 10 subset columns
  # subset_cols <- sample(colnames(to_test_subset), size = 10)
  #
  # # Verify that the subset has NAs before imputation
  # expect_true(anyNA(to_test_subset[, subset_cols]))
  #
  # # Impute only the subset columns using knn_imp without post_imp
  # r1_subset <- knn_imp(
  #   to_test_subset,
  #   k = 10,
  #   method = "euclidean",
  #   post_imp = FALSE,
  #   subset = subset_cols
  # )[, subset_cols]
  #
  # # Verify no NAs remain in the imputed subset
  # expect_true(!anyNA(r1_subset))
  #
  # # Perform full imputation using impute::impute.knn and extract the subset
  # r2_subset <- t(
  #   impute::impute.knn(
  #     t(to_test_subset),
  #     k = 10,
  #     maxp = ncol(to_test_subset)
  #   )$data
  # )[, subset_cols]
  #
  # # Verify no NAs in the extracted subset from full imputation
  # expect_true(!anyNA(r2_subset))
  #
  # # Verify that the subset imputation matches the extracted subset from full imputation
  # expect_equal(r1_subset, r2_subset)
})

test_that("`subset` feature of `knn_imp` works with post_imp = FALSE/TRUE", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 50, perc_NA = 0.2, perc_col_NA = 1)$input)
  # Impute just 3 columns
  ## Check subset using numeric index
  r1 <- knn_imp(to_test, k = 3, post_imp = FALSE, subset = c(1, 3, 5))
  expect_true(!anyNA(r1[, c(1, 3, 5)]))
  expect_equal(is.na(r1[, -c(1, 3, 5)]), is.na(to_test[, -c(1, 3, 5)]))
  ## Check subset using character vector
  r2 <- knn_imp(
    to_test,
    k = 3,
    post_imp = FALSE,
    subset = paste0("feat", c(1, 3, 5))
  )
  expect_equal(r1, r2)

  # Test with post_imp = TRUE and a column requiring post imputation
  to_test_post <- to_test
  # Column 5 will be colMeans if post_imp is TRUE
  to_test_post[2:nrow(to_test_post), 5] <- NA
  r3 <- knn_imp(to_test_post, k = 3, post_imp = TRUE, subset = c(1, 3, 5))
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
  )
  expect_equal(r3, r4)
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- t(sim_mat(m = 20, n = 20, perc_NA = 0.2, perc_col_NA = 1)$input)
  # row 1 is all NA
  to_test[1, ] <- NA
  expect_no_error(knn_imp(to_test, k = 3, post_imp = FALSE))
  expect_true(!anyNA(knn_imp(to_test, k = 3, post_imp = TRUE)))

  to_test[, 1] <- NA
  expect_error(knn_imp(to_test, k = 3, post_imp = FALSE))
})

# test_that("`find_knn_brute` returns correct neighbors as manual implementation", {
#   set.seed(1234)
#   to_test <- t(sim_mat(n = 10, m = 30, perc_NA = 0.5, perc_col_NA = 1)$input)
#   miss <- is.na(to_test)
#   # Ensure columns 1 and 2 have at least one NA
#   expect_true(anyNA(to_test[, 1]))
#   expect_true(anyNA(to_test[, 2]))
#   n_col_miss <- colSums(miss)
#   n_col_name <- colnames(to_test)
#   # Test only column 1 and 2
#   k <- 3
#   n_col_miss[3:length(n_col_miss)] <- 0
#   result <- find_knn_brute(
#     obj = to_test,
#     miss = miss,
#     k = k,
#     n_col_miss = n_col_miss,
#     n_col_name = n_col_name,
#     method = 0,
#     cores = 1
#   )
#   r_dist <- as.matrix(dist(t(to_test)))
#
#   # For each of the two columns, compute expected neighbors
#   for (col_idx in 1:2) {
#     col_name <- paste0("feat", col_idx)
#
#     # Get distances for this column, excluding self (diagonal)
#     distances <- r_dist[, col_idx]
#     distances[col_idx] <- Inf # Exclude self-distance
#
#     # Find k nearest neighbors
#     k_nearest_indices <- order(distances)[1:k]
#     k_nearest_distances <- unname(distances[k_nearest_indices])
#
#     # Check that the function returned the correct indices
#     expect_equal(result[[col_name]]$indices, k_nearest_indices)
#
#     # Check that the function returned the correct distances
#     expect_equal(
#       sqrt(result[[col_name]]$distances * nrow(to_test)),
#       k_nearest_distances
#     )
#
#     # Check that k neighbors were returned
#     expect_equal(result[[col_name]]$n_neighbors, k)
#   }
# })
