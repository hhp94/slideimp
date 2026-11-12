# new_slideimp_results ----
test_that("new_slideimp_results() starts from a bare matrix", {
  m <- matrix(c(1, NA, 3, 4), 2, dimnames = list(c("a", "b"), c("x", "y")))
  first <- new_slideimp_results(m, "knn", FALSE, TRUE, na_check = TRUE)
  # what group_imp() and slide_imp() add on top
  attr(first, "metacaller") <- "group_imp"
  attr(first, "fallback_action") <- "skip"

  # Regression: a result passed back in kept its class, which was then
  # attached a second time, and kept every attribute this constructor does
  # not set.
  second <- new_slideimp_results(first, "pca", FALSE, FALSE, na_check = TRUE)
  expect_identical(class(second), c("slideimp_results", "matrix", "array"))
  expect_setequal(
    names(attributes(second)),
    c(
      "dim",
      "dimnames",
      "class",
      "imp_method",
      "fallback",
      "post_imp",
      "has_remaining_na"
    )
  )
  expect_identical(attr(second, "imp_method"), "pca")
  expect_identical(dimnames(second), dimnames(m))
  expect_identical(unclass(second)[,], m)
  expect_true(attr(second, "has_remaining_na"))

  # no dimnames in, none out
  expect_null(dimnames(
    new_slideimp_results(unname(m), "knn", FALSE, TRUE, na_check = FALSE)
  ))
})

test_that("new_slideimp_results() computes has_remaining_na before attaching the class", {
  skip_if_not(capabilities("profmem"), "R built without memory profiling")

  n <- 200L
  p <- 500L
  # no NA anywhere, so anyNA() has to reach the last cell on either path
  m <- matrix(rnorm(n * p), n)

  # Regression: `has_remaining_na = if (na_check) anyNA(obj) else NULL` is a
  # lazy default and was first forced at the attr() call, after
  # `class(obj) <-`. anyNA() on an object carrying a class is any(is.na(x)),
  # which allocates an n x p logical - 400,048 bytes here - and scans the whole
  # matrix; on the bare matrix it allocates nothing and stops at the first NA.
  f <- tempfile(fileext = ".out")
  Rprofmem(f, threshold = 1e5)
  out <- new_slideimp_results(
    m,
    "pca",
    fallback = FALSE,
    post_imp = FALSE,
    na_check = TRUE
  )
  Rprofmem(NULL)
  big <- grep("^[0-9]+ :", readLines(f), value = TRUE)
  unlink(f)
  expect_length(big, 0L)

  expect_s3_class(out, "slideimp_results")
  expect_false(attr(out, "has_remaining_na"))

  # the flag itself is unchanged
  m[1L] <- NA
  expect_true(attr(
    new_slideimp_results(m, "pca", FALSE, FALSE, na_check = TRUE),
    "has_remaining_na"
  ))
  expect_null(attr(
    new_slideimp_results(m, "pca", FALSE, FALSE, na_check = FALSE),
    "has_remaining_na"
  ))
})
