skip_on_cran()

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
      "post_imp",
      "mean_filled_cols",
      "mean_filled_cells"
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

test_that("mean_filled_cells_matrix orders by column and drops duplicates", {
  expect_identical(
    mean_filled_cells_matrix(),
    cbind(row = integer(0), col = integer(0))
  )
  # doubles in, as the kernels return them, and one cell given twice
  expect_identical(
    mean_filled_cells_matrix(
      row = c(2, 4, 1, 2),
      col = c(3, 1, 3, 3),
      n_row = 5L
    ),
    cbind(row = c(4L, 1L, 2L), col = c(1L, 3L, 3L))
  )
  # `col * n_row` is beyond integer range here
  expect_identical(
    mean_filled_cells_matrix(
      row = c(1L, 1L),
      col = c(3000000L, 2999999L),
      n_row = 1000L
    ),
    cbind(row = c(1L, 1L), col = c(2999999L, 3000000L))
  )
})

test_that("combine_mean_filled merges the records of overlapping windows", {
  # two windows that both filled cell (1, 3) and column 2
  expect_identical(
    combine_mean_filled(
      list(
        list(cols = 2L, row = c(1L, 4L), col = c(3L, 3L)),
        NULL,
        list(cols = c(5L, 2L), row = c(2L, 1L), col = c(1L, 3L))
      ),
      n_row = 4L
    ),
    list(
      cols = c(2L, 5L),
      cells = cbind(row = c(2L, 1L, 4L), col = c(1L, 3L, 3L))
    )
  )
  # no record at all still combines to the two empty attributes
  expect_identical(
    combine_mean_filled(list(NULL, NULL), n_row = 4L),
    list(cols = integer(0), cells = mean_filled_cells_matrix())
  )
})

test_that("print reports what column means replaced", {
  m <- matrix(as.numeric(1:6), 2)
  note_lines <- function(x) {
    out <- capture.output(print(x))
    out[grepl("^Note: |attribute", out)]
  }

  # nothing filled, nothing reported
  expect_length(note_lines(new_slideimp_results(m, "knn", FALSE, TRUE)), 0L)

  one <- new_slideimp_results(
    m,
    "knn",
    FALSE,
    TRUE,
    mean_filled_cols = 2L,
    mean_filled_cells = mean_filled_cells_matrix(1L, 3L, 2L)
  )
  expect_identical(
    note_lines(one),
    c(
      "Note: 1 column was excluded from imputation. Column means replaced its missing values.",
      "  See attribute \"mean_filled_cols\".",
      "Note: 1 cell could not be imputed. Column means replaced it.",
      "  See attribute \"mean_filled_cells\"."
    )
  )

  many <- new_slideimp_results(
    m,
    "knn",
    FALSE,
    TRUE,
    mean_filled_cols = 1:2,
    mean_filled_cells = mean_filled_cells_matrix(1:2, c(3L, 3L), 2L)
  )
  expect_identical(
    note_lines(many),
    c(
      "Note: 2 columns were excluded from imputation. Column means replaced their missing values.",
      "  See attribute \"mean_filled_cols\".",
      "Note: 2 cells could not be imputed. Column means replaced them.",
      "  See attribute \"mean_filled_cells\"."
    )
  )

  # an object saved before the attributes existed prints without them
  old <- one
  attr(old, "mean_filled_cols") <- NULL
  attr(old, "mean_filled_cells") <- NULL
  expect_length(note_lines(old), 0L)
})

# mirai_daemons_set ----
test_that("mirai_daemons_set() answers without raising", {
  skip_if(mirai::daemons_set(), "daemons are set in this session")
  expect_false(mirai_daemons_set())
})
