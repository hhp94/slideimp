test_that("mat_miss counts NAs correctly", {
  m <- matrix(
    c(1, NA, 3, NA, NA, 6),
    nrow = 2,
    dimnames = list(c("r1", "r2"), c("a", "b", "c"))
  )

  expect_equal(mat_miss(m), c(a = 1, b = 1, c = 1))
  expect_equal(mat_miss(m, col = FALSE), c(r1 = 1, r2 = 2))
  expect_equal(mat_miss(m, prop = TRUE), c(a = 0.5, b = 0.5, c = 0.5))
  expect_equal(mat_miss(m, col = FALSE, prop = TRUE), c(r1 = 1 / 3, r2 = 2 / 3))
})

test_that("mat_miss handles edge cases", {
  # no NAs
  m1 <- matrix(1:6, 2, 3)
  expect_equal(unname(mat_miss(m1)), c(0, 0, 0))

  # all NAs
  m2 <- matrix(NA_real_, 2, 3)
  expect_equal(unname(mat_miss(m2, prop = TRUE)), c(1, 1, 1))

  # NaN should also count
  m3 <- matrix(c(1, NaN, NA, 4), 2, 2)
  expect_equal(unname(mat_miss(m3)), c(1, 1))
})

test_that("mat_miss returns integer counts and double proportions", {
  m <- matrix(c(1, NA, 3, NA, NA, 6), nrow = 2)

  expect_type(mat_miss(m), "integer")
  expect_type(mat_miss(m, col = FALSE), "integer")
  expect_type(mat_miss(m, prop = TRUE), "double")
  expect_type(mat_miss(m, col = FALSE, prop = TRUE), "double")
})

test_that("mat_miss does not count Inf as missing", {
  m <- matrix(c(1, Inf, NA, NaN, -Inf, 6), nrow = 3)
  expect_equal(unname(mat_miss(m)), c(1L, 1L))
})

test_that("mat_miss validates inputs", {
  expect_error(mat_miss(data.frame(x = 1)), "matrix")
  expect_error(mat_miss(matrix("a", 1, 1)), "numeric")
  expect_error(mat_miss(matrix(1, 1, 1), col = "yes"))
})

test_that("mat_miss agrees with base R on sim_mat output", {
  set.seed(1)
  s <- sim_mat(n = 50, p = 20, perc_total_na = 0.2)
  expect_equal(unname(mat_miss(s$input)), unname(colSums(is.na(s$input))))
  expect_equal(
    unname(mat_miss(s$input, col = FALSE)),
    unname(rowSums(is.na(s$input)))
  )
})

test_that("the package is built with 64-bit Armadillo indices", {
  # -DARMA_64BIT_WORD=1 in src/Makevars. Without it arma indices are 32 bits,
  # capping a matrix at 2^32-1 elements; the zero-copy constructor used for
  # `const arma::mat&` does not check that bound, so n_elem and colptr()
  # arithmetic wrap silently instead of erroring. At 850k columns the cap lands
  # near 5,052 rows, which real data reaches.
  expect_equal(arma_uword_bytes(), 8L)
})
