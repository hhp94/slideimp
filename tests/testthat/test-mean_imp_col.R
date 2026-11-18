test_that("`mean_imp_col` works", {
  set.seed(4321)
  to_test <- sim_mat(10, 20, perc_total_na = 0.5, perc_col_na = 1)$input
  c_manual <- to_test
  r_manual <- to_test
  na_indices <- which(is.na(to_test), arr.ind = TRUE)
  column_means <- colMeans(to_test, na.rm = TRUE)
  row_means <- rowMeans(to_test, na.rm = TRUE)

  c_manual[na_indices] <- column_means[na_indices[, 2]]
  # each hole takes its column's colMeans(), the same arithmetic in the same
  # process, so the result is identical, not just close
  expect_identical(mean_imp_col(to_test), c_manual)

  ## Test subset feature
  c_subset <- to_test
  for (i in c(1, 5, 10)) {
    c_subset[is.na(c_subset[, i]), i] <- mean(c_subset[, i], na.rm = TRUE)
  }
  expect_equal(mean_imp_col(to_test, subset = c(1, 5, 10)), c_subset)
})

test_that("`mean_imp_col` refuses Inf, including when subset resolves to NULL", {
  set.seed(1234)
  to_test <- sim_mat(20, 10, perc_total_na = 0.2)$input
  to_test[3, 4] <- Inf

  expect_error(mean_imp_col(to_test), "Infinite")
  expect_error(mean_imp_col(to_test, subset = 4), "Infinite")

  # the Inf sits in a column that is not selected, and the selected columns
  # have no missingness, so `subset` resolves to NULL and the old code path
  # returned the input unchanged with the Inf still in it.
  complete <- sim_mat(20, 10, perc_total_na = 0)$input
  complete[3, 4] <- Inf
  expect_error(mean_imp_col(complete, subset = 1), "Infinite")
})

# A 6 x 5 input: column 1 with two holes, column 2 complete, column 3 all NA,
# column 4 with a NaN and an NA, column 5 with one hole. Each column with a
# hole and an observed value is skewed, so its mean and median differ.
mean_input <- function() {
  m <- matrix(
    c(
      c(1, NA, 3, 4, NA, 10),
      c(2, 4, 6, 8, 10, 12),
      rep(NA, 6),
      c(NaN, 1.5, 2.5, NA, 3.5, 8.5),
      c(10, 20, 30, 40, 90, NA)
    ),
    6,
    5,
    dimnames = list(paste0("r", 1:6), paste0("c", 1:5))
  )
  m
}

test_that("fallback_part() maps each hole of its columns to the column statistic", {
  m <- mean_input()
  cols <- c(1L, 2L, 4L, 5L)
  # the holes in column-major order; the complete column 2 has none, and the
  # NaN in column 4 is a hole like the NA
  expected_map <- cbind(
    row = c(2L, 5L, 1L, 4L, 6L),
    col = c(1L, 1L, 4L, 4L, 5L)
  )
  stats <- list(
    mean = colMeans(m, na.rm = TRUE),
    median = apply(m, 2L, stats::median, na.rm = TRUE)
  )
  # the two statistics differ at every filled column, so a part that took
  # the wrong one fails
  expect_true(all(stats$mean[c(1, 4, 5)] != stats$median[c(1, 4, 5)]))

  for (method in names(stats)) {
    part <- fallback_part(m, cols, method)
    expect_s3_class(part, "slideimp_part")
    expect_identical(part$kind, method)
    expect_false(part$left_na)
    expect_identical(part$cols, cols)
    expect_identical(part$map, expected_map)
    expect_identical(
      part$values,
      unname(stats[[method]][expected_map[, "col"]])
    )
    # validate_part() passes it, returning the holes' linear indices
    expect_identical(
      validate_part(part, nrow(m), ncol(m)),
      c(2, 5, 19, 22, 30)
    )

    # no columns, and columns without a hole: an empty map
    none <- fallback_part(m, integer(0), method)
    expect_identical(none$kind, method)
    expect_identical(nrow(none$map), 0L)
    expect_identical(none$values, numeric(0))
    complete <- fallback_part(m, 2L, method)
    expect_identical(complete$cols, 2L)
    expect_identical(nrow(complete$map), 0L)
    expect_identical(complete$values, numeric(0))

    # a column with no observed value has no statistic, and the door
    # refuses it
    expect_error(
      fill_parts(m, list(fallback_part(m, 3L, method))),
      "`values` holds NA",
      fixed = TRUE
    )
  }

  # with no dimnames, the same part
  expect_identical(fallback_part(unname(m), cols), fallback_part(m, cols))
})

test_that("mean_imp_col() leaves a column with no observed value as it is", {
  m <- mean_input()
  expect_true(all(is.na(m[, 3])))
  for (subset in list(NULL, 3L, c(1L, 3L))) {
    out <- mean_imp_col(m, subset = subset)
    expect_identical(out[, 3], m[, 3])
  }
  out <- mean_imp_col(m)
  expect_identical(sum(is.na(out)), 6L)
  expect_identical(out[, -3], {
    e <- m[, -3]
    e[is.na(e)] <- colMeans(m[, -3], na.rm = TRUE)[col(e)[is.na(e)]]
    e
  })
})

test_that("mean_imp_col() hands the fallback only the columns it fills", {
  # Regression: every selected column with an observed value went into the
  # fallback part, complete ones included, so fallback_part() copied them and
  # their is.na() to fill nothing. The slide_imp() and group_imp() fallbacks
  # select every target column of a window or group.
  m <- mean_input()
  # non-vacuity: column 2 is complete and column 3 has no observed value
  expect_false(anyNA(m[, 2]))
  expect_true(all(is.na(m[, 3])))
  real_fallback_part <- fallback_part
  calls <- list()
  local_mocked_bindings(
    fallback_part = function(obj, cols, method = c("mean", "median")) {
      calls <<- c(calls, list(cols))
      real_fallback_part(obj, cols, method)
    }
  )
  out <- mean_imp_col(m)
  expect_identical(calls, list(c(1L, 4L, 5L)))
  expect_identical(out[, 2], m[, 2])
  expect_false(anyNA(out[, c(1, 4, 5)]))

  # complete columns only: an empty part, and the input's values back
  calls <- list()
  expect_identical(mean_imp_col(m, subset = 2L), m)
  expect_identical(calls, list(integer(0)))
})

test_that("mean_imp_col() returns a double matrix with dim and dimnames only", {
  m <- mean_input()
  plain <- mean_imp_col(m, subset = c(1L, 5L))

  # the order of `subset` does not matter
  expect_identical(mean_imp_col(m, subset = c(5L, 1L)), plain)
  expect_identical(mean_imp_col(m, subset = c("c5", "c1")), plain)
  # columns outside `subset` keep their holes
  expect_identical(plain[, c(2, 3, 4)], m[, c(2, 3, 4)])
  expect_false(anyNA(plain[, c(1, 5)]))

  # any other attribute of the input is dropped, whether the door writes or not
  tagged <- structure(m, class = c("tagged", "matrix"), note = "x")
  expect_identical(mean_imp_col(tagged, subset = c(1L, 5L)), plain)
  expect_identical(mean_imp_col(tagged, subset = 2L), m)

  # no dimnames in, none out
  bare <- unname(m)
  expect_identical(mean_imp_col(bare, subset = c(1L, 5L)), unname(plain))

  # integer storage comes back double
  mi <- matrix(c(1L, NA, 3L, 4L, 5L, NA), 3, 2)
  out <- mean_imp_col(mi)
  expect_true(is.double(out))
  expect_identical(out, matrix(c(1, 2, 3, 4, 5, 4.5), 3, 2))

  # a subset that resolves to nothing returns the input itself, attributes
  # and storage included
  expect_identical(
    suppressMessages(mean_imp_col(tagged, subset = integer(0))),
    tagged
  )
})
