# pca_imp(subset = ): the result narrowed to the selected columns ----

# A named 40 x 12 input. Column 6 is over the default colmax (0.9), column 9
# is constant with no hole, and about half the columns have holes.
subset_input <- function() {
  set.seed(1)
  obj <- sim_mat(n = 40, p = 12, perc_total_na = 0.1, perc_col_na = 0.5)$input
  obj[1:37, 6] <- NA
  obj[, 9] <- 0.5
  obj
}

# The values, dims and dimnames only, as as.matrix() of a fit returns them.
bare_matrix <- function(x) {
  attributes(x) <- list(dim = dim(x), dimnames = dimnames(x))
  x
}

# Every attribute of a matrix output but dim, dimnames and mean_filled_cols:
# the record that `subset` must not change.
record_but_cols <- function(x) {
  a <- attributes(x)
  a[setdiff(names(a), c("dim", "dimnames", "mean_filled_cols"))]
}

test_that("subset returns the full output's columns, in input order", {
  obj <- subset_input()
  holed <- which(colSums(is.na(obj)) > 0)
  # non-vacuity: eligible columns with holes both inside and outside the
  # subsets below, so the fit uses columns the result does not hold
  expect_true(all(c(2L, 6L) %in% holed))
  expect_gt(length(setdiff(holed, c(2L, 6L, 9L, 11L))), 1L)

  subsets <- list(
    sorted = c(2L, 6L, 9L, 11L),
    unsorted = c(11L, 2L, 9L, 6L),
    names = colnames(obj)[c(11, 2, 9, 6)],
    double = c(2, 6, 9, 11),
    one = 2L
  )
  for (post_imp in c(TRUE, FALSE)) {
    full <- pca_imp(obj, ncp = 2, seed = 1, post_imp = post_imp)
    for (nm in names(subsets)) {
      info <- paste("post_imp", post_imp, nm)
      sel <- sort(match(
        colnames(obj[, subsets[[nm]], drop = FALSE]),
        colnames(obj)
      ))
      out <- pca_imp(
        obj,
        ncp = 2,
        seed = 1,
        post_imp = post_imp,
        subset = subsets[[nm]]
      )
      expect_s3_class(out, "slideimp_results")
      expect_identical(
        bare_matrix(out),
        bare_matrix(full)[, sel, drop = FALSE],
        info = info
      )
      expect_identical(record_but_cols(out), record_but_cols(full), info = info)
    }
  }
})

test_that("a subset of every column is the call without subset", {
  obj <- subset_input()
  full <- pca_imp(obj, ncp = 2, seed = 1)
  expect_identical(pca_imp(obj, ncp = 2, seed = 1, subset = 12:1), full)
  expect_identical(
    pca_imp(obj, ncp = 2, seed = 1, subset = rev(colnames(obj))),
    full
  )
})

test_that("mean_filled_cols counts the result's columns", {
  obj <- subset_input()
  full <- pca_imp(obj, ncp = 2, seed = 1)
  out <- pca_imp(obj, ncp = 2, seed = 1, subset = c(2L, 6L, 9L, 11L))
  # column 6 is the one column over colmax with a hole: the input's 6th
  # column, and the result's 2nd
  expect_identical(attr(full, "mean_filled_cols"), 6L)
  expect_identical(attr(out, "mean_filled_cols"), 2L)
  expect_identical(colnames(out)[2], colnames(obj)[6])

  # a subset without column 6 records none
  out <- pca_imp(obj, ncp = 2, seed = 1, subset = c(2L, 11L))
  expect_identical(attr(out, "mean_filled_cols"), integer(0))
})

test_that("the mipca output with subset is the narrowed matrix output", {
  obj <- subset_input()
  sel <- c(2L, 6L, 9L, 11L)
  for (post_imp in c(TRUE, FALSE)) {
    m <- pca_imp(obj, ncp = 2, seed = 1, post_imp = post_imp, subset = sel)
    fit <- pca_imp(
      obj,
      ncp = 2,
      seed = 1,
      post_imp = post_imp,
      subset = sel,
      output = "mipca"
    )
    expect_identical(
      as.matrix(fit),
      bare_matrix(m),
      info = paste("post_imp", post_imp)
    )
    expect_identical(fit$rownames, rownames(obj))
    # the record on every chunk is the matrix output's, in its numbering
    rec <- attributes(m)[setdiff(
      names(attributes(m)),
      c("dim", "dimnames", "class")
    )]
    for (ch in fit$chunks) {
      expect_identical(ch$info, rec)
    }
  }

  # without post_imp, column 6 is the one column left NA
  fit <- pca_imp(
    obj,
    ncp = 2,
    seed = 1,
    post_imp = FALSE,
    subset = sel,
    output = "mipca"
  )
  expect_output(print(fit), "left NA: 1 column not imputed", fixed = TRUE)
  out <- as.matrix(fit)
  expect_identical(unname(which(colSums(is.na(out)) > 0)), 2L)
})

test_that("a subset without holes returns its columns without a fit", {
  obj <- subset_input()
  done <- unname(which(colSums(is.na(obj)) == 0))
  expect_gt(length(done), 1L)
  expect_true(9L %in% done)

  # the kernel is never reached; a subset with a hole shows the stub is in
  # place
  local_mocked_bindings(
    pca_imp_internal_cpp = function(...) stop("the kernel ran")
  )
  expect_error(pca_imp(obj, ncp = 2, subset = 2L), "the kernel ran")

  expect_message(
    out <- pca_imp(obj, ncp = 2, subset = done),
    "No missing values in the `subset` columns",
    fixed = TRUE
  )
  expect_s3_class(out, "slideimp_results")
  expect_identical(bare_matrix(out), obj[, done, drop = FALSE])
  expect_identical(attr(out, "n_iter"), 0L)

  expect_message(
    fit <- pca_imp(obj, ncp = 2, subset = done, output = "mipca"),
    "No missing values in the `subset` columns",
    fixed = TRUE
  )
  expect_length(fit$chunks, 1L)
  expect_identical(nrow(fit$chunks[[1]]$map), 0L)
  expect_identical(as.matrix(fit), obj[, done, drop = FALSE])
})

test_that("a subset whose holes are all in ineligible columns is infeasible", {
  obj <- subset_input()
  expect_true(anyNA(obj[, 6]))
  expect_false(anyNA(obj[, 9]))
  expect_error(
    pca_imp(obj, ncp = 2, subset = c(6L, 9L)),
    class = "slideimp_infeasible"
  )
  expect_error(
    pca_imp(obj, ncp = 2, subset = c(6L, 9L)),
    "All `subset` columns with missing values are ineligible",
    fixed = TRUE
  )
})

test_that("subset refuses what it cannot select, before the fit", {
  obj <- subset_input()
  local_mocked_bindings(
    pca_imp_internal_cpp = function(...) stop("the kernel ran")
  )
  refuse <- function(subset, msg, x = obj) {
    cnd <- expect_error(pca_imp(x, ncp = 2, subset = subset), msg, fixed = TRUE)
    expect_identical(conditionCall(cnd)[[1]], as.name("pca_imp"))
  }
  refuse(integer(0), "`subset` selects no columns.")
  refuse(character(0), "`subset` selects no columns.")
  refuse(
    c(colnames(obj)[1], "nope", "nada"),
    "`subset` names 2 columns not in `obj`: nope, nada."
  )
  refuse(
    colnames(obj)[c(1, 2, 1)],
    "`subset` names 1 column more than once"
  )
  refuse(
    colnames(obj)[1:2],
    "`subset` names columns, but `obj` has no column names.",
    x = unname(obj)
  )

  # numbers go through checkmate, which names the argument
  for (bad in list(13L, 0L, 2.5, c(2L, 2L), c(2L, NA))) {
    expect_error(
      pca_imp(obj, ncp = 2, subset = bad),
      "subset",
      info = deparse(bad)
    )
  }
  expect_error(pca_imp(obj, ncp = 2, subset = TRUE), "subset")
})

test_that("subset is not a parameter tune_imp() takes", {
  obj <- subset_input()
  expect_error(
    tune_imp(
      obj,
      parameters = data.frame(ncp = 2, subset = 2L),
      .f = "pca_imp",
      n_reps = 1
    ),
    "Unknown parameter"
  )
})
