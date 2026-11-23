skip_on_cran()

# knn_imp() on parts: `subset` narrows the result ----

# A named 40 x 12 input. Column 6 is over the default colmax (0.9), and about
# half the columns have holes.
knn_subset_input <- function() {
  set.seed(1)
  obj <- sim_mat(n = 40, p = 12, perc_total_na = 0.1, perc_col_na = 0.5)$input
  obj[1:37, 6] <- NA
  obj
}

# The record of a result: its attributes after dim, dimnames and class.
knn_record <- function(x) {
  a <- attributes(x)
  a[setdiff(names(a), c("dim", "dimnames", "class"))]
}

test_that("knn_left_input() leaves the cells it claims to", {
  obj <- knn_left_input()
  r <- knn_imp(obj, k = 1, post_imp = FALSE)
  expect_identical(
    unname(which(is.na(r[, -7]), arr.ind = TRUE)),
    cbind(c(1L, 1L), c(2L, 3L))
  )
  expect_false(is.na(r[5, 2]))
  expect_false(anyNA(r[, 5]))
  expect_identical(sum(is.na(r[, 7])), 11L)
})

test_that("subset returns the full output's columns, in input order", {
  obj <- knn_subset_input()
  holed <- which(colSums(is.na(obj)) > 0)
  # non-vacuity: columns with holes both inside and outside the subsets
  # below, so the kernel has candidates the result does not hold
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
    # every column is imputed from the same candidates with or without
    # subset, so within one process the values are the same bits
    full <- knn_imp(obj, k = 3, post_imp = post_imp)
    for (nm in names(subsets)) {
      info <- paste("post_imp", post_imp, nm)
      sel <- sort(match(
        colnames(obj[, subsets[[nm]], drop = FALSE]),
        colnames(obj)
      ))
      out <- knn_imp(obj, k = 3, post_imp = post_imp, subset = subsets[[nm]])
      expect_s3_class(out, "slideimp_results")
      expect_identical(
        bare_knn(out),
        bare_knn(full)[, sel, drop = FALSE],
        info = info
      )
      # the record is the full call's, cut to the selected columns and
      # numbered as the result's columns. This input leaves no cell for
      # `mean_filled_cells` (the record test below covers one that does), so
      # only the cut of `mean_filled_cols` is exercised here.
      rec <- knn_record(out)
      rec_full <- knn_record(full)
      cols_full <- rec_full$mean_filled_cols
      rec_full$mean_filled_cols <- match(cols_full[cols_full %in% sel], sel)
      cells <- rec_full$mean_filled_cells
      cells <- cells[cells[, "col"] %in% sel, , drop = FALSE]
      cells[, "col"] <- match(cells[, "col"], sel)
      rec_full$mean_filled_cells <- cells
      expect_identical(rec, rec_full, info = info)
    }
  }

  # a subset of every column, in any order, is the call without it
  full <- knn_imp(obj, k = 3)
  expect_identical(knn_imp(obj, k = 3, subset = 12:1), full)
  expect_identical(knn_imp(obj, k = 3, subset = rev(colnames(obj))), full)
})

test_that("the record counts the result's columns", {
  obj <- knn_left_input()
  full <- knn_imp(obj, k = 1)
  # c7 is over colmax, and K-NN left row 1 of c2 and c3
  expect_identical(attr(full, "mean_filled_cols"), 7L)
  expect_identical(
    attr(full, "mean_filled_cells"),
    cbind(row = c(1L, 1L), col = c(2L, 3L))
  )

  # in the result of c3, c5 and c7 they are its columns 3 and 1
  out <- knn_imp(obj, k = 1, subset = c("c7", "c3", "c5"))
  expect_identical(colnames(out), c("c3", "c5", "c7"))
  expect_identical(attr(out, "mean_filled_cols"), 3L)
  expect_identical(attr(out, "mean_filled_cells"), cbind(row = 1L, col = 1L))
  expect_identical(bare_knn(out), bare_knn(full)[, c(3, 5, 7)])
})

test_that("a subset without holes returns its columns without the kernel", {
  obj <- knn_subset_input()
  done <- unname(which(colSums(is.na(obj)) == 0))
  expect_gt(length(done), 1L)

  # the kernel is never reached; a subset with a hole shows the stub is in
  # place
  local_mocked_bindings(
    impute_knn_brute = function(...) stop("the kernel ran")
  )
  expect_error(knn_imp(obj, k = 3, subset = 2L), "the kernel ran")

  expect_message(
    out <- knn_imp(obj, k = 3, subset = done),
    "No missing values in the `subset` columns",
    fixed = TRUE
  )
  expect_s3_class(out, "slideimp_results")
  expect_identical(bare_knn(out), obj[, done, drop = FALSE])
})

test_that("subset refuses what it cannot select, before the kernel", {
  obj <- knn_subset_input()
  local_mocked_bindings(
    impute_knn_brute = function(...) stop("the kernel ran")
  )
  refuse <- function(subset, msg, x = obj) {
    cnd <- expect_error(knn_imp(x, k = 3, subset = subset), msg, fixed = TRUE)
    expect_identical(conditionCall(cnd)[[1]], as.name("knn_imp"))
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
      knn_imp(obj, k = 3, subset = bad),
      "subset",
      info = deparse(bad)
    )
  }
})
