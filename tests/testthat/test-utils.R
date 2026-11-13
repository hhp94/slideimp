# new_slideimp_results ----
test_that("new_slideimp_results() starts from a bare matrix", {
  m <- matrix(c(1, NA, 3, 4), 2, dimnames = list(c("a", "b"), c("x", "y")))
  first <- new_slideimp_results(m, "knn", FALSE, TRUE)
  # what group_imp() and slide_imp() add on top
  attr(first, "metacaller") <- "group_imp"
  attr(first, "fallback_action") <- "skip"

  # Regression: a result passed back in kept its class, which was then
  # attached a second time, and kept every attribute this constructor does
  # not set.
  second <- new_slideimp_results(first, "pca", FALSE, FALSE)
  expect_identical(class(second), c("slideimp_results", "matrix", "array"))
  expect_setequal(
    names(attributes(second)),
    c(
      "dim",
      "dimnames",
      "class",
      "imp_method",
      "fallback",
      "post_imp"
    )
  )
  expect_identical(attr(second, "imp_method"), "pca")
  expect_identical(dimnames(second), dimnames(m))
  expect_identical(unclass(second)[,], m)

  # no dimnames in, none out
  expect_null(dimnames(
    new_slideimp_results(unname(m), "knn", FALSE, TRUE)
  ))
})
