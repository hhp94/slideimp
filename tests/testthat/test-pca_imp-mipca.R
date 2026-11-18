# pca_imp(output = "mipca"): the point estimate as a one-draw si_mi ----

# A named 40 x 12 input. Column 6 is over the default colmax (0.9), in the
# middle, so the parts' columns interleave and the fit needs a placement.
# Column 9 is constant with no hole: ineligible, but with nothing to fill.
mipca_input <- function() {
  set.seed(1)
  obj <- sim_mat(n = 40, p = 12, perc_total_na = 0.1, perc_col_na = 0.5)$input
  obj[1:37, 6] <- NA
  obj[, 9] <- 0.5
  obj
}

# The matrix output reduced to what as.matrix() of a fit returns: the values,
# the dims and the dimnames, with no class or other attribute.
bare <- function(x) {
  attributes(x) <- list(dim = dim(x), dimnames = dimnames(x))
  x
}

test_that("as.matrix() of the mipca output is the matrix output", {
  obj <- mipca_input()
  # non-vacuity: holes in eligible columns, an ineligible column with holes,
  # and an ineligible column without
  expect_gt(sum(is.na(obj[, -c(6, 9)])), 0L)
  expect_gt(sum(is.na(obj[, 6])), 0.9 * nrow(obj))
  expect_false(anyNA(obj[, 9]))

  for (post_imp in c(TRUE, FALSE)) {
    m <- pca_imp(obj, ncp = 2, seed = 1, post_imp = post_imp)
    fit <- pca_imp(
      obj,
      ncp = 2,
      seed = 1,
      post_imp = post_imp,
      output = "mipca"
    )
    expect_s3_class(fit, "si_mi")
    expect_identical(fit$nkeep, 1L)
    out <- as.matrix(fit)
    expect_identical(out, bare(m), info = paste("post_imp", post_imp))
    # columns in input order, through a placement that is not the identity
    expect_identical(colnames(out), colnames(obj))
    expect_false(is.null(fit$placement))
    # the same cells left NA: none with post_imp, column 6 without
    expect_identical(is.na(out), is.na(m))
    expect_identical(anyNA(out), !post_imp)
    if (!post_imp) {
      expect_identical(unname(which(colSums(is.na(out)) > 0)), 6L)
    }
  }
})

test_that("the mipca output's chunks carry the matrix output's record", {
  obj <- mipca_input()
  m <- pca_imp(obj, ncp = 2, seed = 1, post_imp = FALSE)
  fit <- pca_imp(obj, ncp = 2, seed = 1, post_imp = FALSE, output = "mipca")

  # three parts: the eligible columns, the ineligible column left NA, and the
  # constant column with nothing to fill
  expect_length(fit$chunks, 3L)
  cols <- lapply(fit$chunks, function(ch) colnames(ch$obj))
  expect_identical(cols[[2]], colnames(obj)[6])
  expect_identical(cols[[3]], colnames(obj)[9])
  expect_identical(
    vapply(fit$chunks, function(ch) isTRUE(ch$prov$left_na), TRUE),
    c(FALSE, TRUE, FALSE)
  )

  rec <- attributes(m)[setdiff(
    names(attributes(m)),
    c("dim", "dimnames", "class")
  )]
  expect_length(rec, 14L)
  for (ch in fit$chunks) {
    expect_identical(ch$info, rec)
  }
  expect_identical(fit$chunks[[1]]$prov$rank, 2L)
  expect_identical(fit$chunks[[1]]$prov$iter, attr(m, "n_iter"))
  expect_null(fit$chunks[[2]]$prov$rank)

  # the maps are in each chunk's own numbering, integer, and the chunks that
  # are slices hold no row names: the object holds them once
  for (ch in fit$chunks) {
    expect_true(is.integer(ch$map))
    expect_null(rownames(ch$obj))
  }
  expect_identical(fit$rownames, rownames(obj))

  # print() counts only the column that keeps its holes as left NA
  expect_output(print(fit), "1 draw", fixed = TRUE)
  expect_output(print(fit), "left NA: 1 column not imputed", fixed = TRUE)

  # the call is stored without the matrix, even when the matrix is passed by
  # value
  fit2 <- do.call(pca_imp, list(obj = obj, ncp = 2, output = "mipca"))
  expect_identical(fit2$call[[1]], as.name("pca_imp"))
  expect_identical(fit2$call$obj, "<matrix>")
})

test_that("an all-eligible input is one chunk holding the input itself", {
  obj <- mipca_input()[, -c(6, 9)]
  fit <- pca_imp(obj, ncp = 2, seed = 1, output = "mipca")
  expect_length(fit$chunks, 1L)
  expect_null(fit$placement)
  # not a slice: the row names are kept on the chunk too
  expect_identical(fit$chunks[[1]]$obj, obj)
  expect_identical(
    as.matrix(fit),
    bare(pca_imp(obj, ncp = 2, seed = 1))
  )
})

test_that("a one-chunk fit drops the input's class and attributes", {
  # Regression: the one-chunk fit held the input itself, and as.matrix()
  # starts from that matrix, so a `slideimp_results` from an earlier pass
  # came back with its class and that pass's record (`imp_method = "knn"`
  # after a knn_imp()), and a plain matrix with its own attributes. A fit of
  # more chunks is built from slices, which carry neither.
  obj <- mipca_input()[, -c(6, 9)]
  earlier <- pca_imp(obj, ncp = 2, seed = 1)
  again <- earlier
  again[cbind(1:5, 1:5)] <- NA
  tagged <- obj
  attr(tagged, "foo") <- "bar"
  inputs <- list(
    results_with_holes = again,
    results_without_holes = earlier,
    plain_with_attribute = tagged
  )
  # non-vacuity: each input carries more than dims and dimnames, and the
  # first one has holes, so it takes the main path, not the no-holes exit
  for (x in inputs) {
    expect_gt(length(setdiff(names(attributes(x)), c("dim", "dimnames"))), 0L)
  }
  expect_true(anyNA(again))
  expect_false(anyNA(earlier))

  for (nm in names(inputs)) {
    x <- inputs[[nm]]
    fit <- suppressMessages(pca_imp(x, ncp = 2, seed = 1, output = "mipca"))
    # the path that reuses the input rather than slicing it
    expect_length(fit$chunks, 1L)
    expect_null(fit$placement)
    out <- as.matrix(fit)
    expect_identical(names(attributes(out)), c("dim", "dimnames"), info = nm)
    expect_identical(
      out,
      bare(suppressMessages(pca_imp(x, ncp = 2, seed = 1))),
      info = nm
    )
    expect_identical(mipca_complete(fit, m = 1)[[1]], out, info = nm)
  }
  # the inputs keep what they carry
  expect_s3_class(again, "slideimp_results")
  expect_identical(attr(tagged, "foo"), "bar")
})

test_that("an input with no holes is one part with an empty map", {
  obj <- mipca_input()
  obj[is.na(obj)] <- 0.5
  expect_message(
    fit <- pca_imp(obj, ncp = 2, output = "mipca"),
    "No missing values"
  )
  m <- suppressMessages(pca_imp(obj, ncp = 2))
  expect_length(fit$chunks, 1L)
  expect_identical(nrow(fit$chunks[[1]]$map), 0L)
  expect_identical(as.matrix(fit), bare(m))
  expect_identical(fit$chunks[[1]]$info$n_iter, 0)

  # Integer input with no holes: the matrix output keeps it integer, while
  # as.matrix() of the fit returns double, as it does for every input with a
  # hole. The fit holds the integer input itself; as.matrix.si_mi() assigns
  # its (empty) double draws into it, and R coerces on that assignment.
  oi <- matrix(
    1:60 %% 7L,
    12,
    5,
    dimnames = list(paste0("r", 1:12), paste0("c", 1:5))
  )
  fi <- suppressMessages(pca_imp(oi, ncp = 2, output = "mipca"))
  mi <- suppressMessages(pca_imp(oi, ncp = 2))
  expect_true(is.integer(mi))
  expect_true(is.integer(fi$chunks[[1]]$obj))
  expected <- bare(mi)
  storage.mode(expected) <- "double"
  expect_identical(as.matrix(fi), expected)
})

test_that("output = 'mipca' refuses an input without unique names, before the fit", {
  obj <- mipca_input()
  bad <- list(
    no_rownames = `rownames<-`(obj, NULL),
    no_colnames = `colnames<-`(obj, NULL),
    no_dimnames = unname(obj),
    repeated_row = `rownames<-`(obj, rep(c("a", "b"), 20)),
    repeated_col = `colnames<-`(obj, c("x", "x", colnames(obj)[-(1:2)])),
    empty_col = `colnames<-`(obj, c("", colnames(obj)[-1]))
  )
  # the kernel is never reached; the named input shows the stub is in place
  local_mocked_bindings(
    pca_imp_internal_cpp = function(...) stop("the kernel ran")
  )
  expect_error(pca_imp(obj, ncp = 2, output = "mipca"), "the kernel ran")
  for (nm in names(bad)) {
    expect_error(
      pca_imp(bad[[nm]], ncp = 2, output = "mipca"),
      "must have unique row and column names",
      info = nm
    )
    # the rule mipca() applies
    expect_error(.chk_matrix(bad[[nm]]), info = nm)
  }
  expect_error(
    pca_imp(bad$repeated_row, ncp = 2, output = "mipca"),
    "repeats 2 row names: a, b.",
    fixed = TRUE
  )
  cnd <- expect_error(
    pca_imp(bad$no_dimnames, ncp = 2, output = "mipca"),
    "has no column names",
    fixed = TRUE
  )
  # the error names the function the user called
  expect_identical(conditionCall(cnd)[[1]], as.name("pca_imp"))
  # the no-holes exit refuses too
  full <- bad$no_colnames
  full[is.na(full)] <- 0.5
  expect_error(
    pca_imp(full, ncp = 2, output = "mipca"),
    "must have unique row and column names"
  )
  expect_error(pca_imp(obj, ncp = 2, output = "fit"), "should be one of")
})

test_that("the matrix output still takes unnamed input", {
  obj <- unname(mipca_input())
  out <- pca_imp(obj, ncp = 2, seed = 1)
  expect_s3_class(out, "slideimp_results")
  expect_false(anyNA(out))
  expect_null(dimnames(out))
})

test_that("a one-draw fit completes once, summarizes nothing, and joins one-draw fits only", {
  obj <- mipca_input()
  fit <- pca_imp(obj, ncp = 2, seed = 1, output = "mipca")

  done <- mipca_complete(fit, m = 1)
  expect_length(done, 1L)
  expect_identical(done[[1]], as.matrix(fit))
  expect_error(mipca_complete(fit, m = 2), "kept draws", fixed = TRUE)

  s <- summary(fit)
  expect_identical(nrow(s$cells), 0L)

  # a mipca() fit over other columns of the same rows
  set.seed(2)
  y <- sim_mat(n = 40, p = 8, perc_total_na = 0.1, perc_col_na = 0.5)$input
  dimnames(y) <- list(rownames(obj), paste0("y", 1:8))
  expect_true(anyNA(y))
  many <- mipca(
    y,
    rank = 2L,
    transform = "identity",
    chains = 1L,
    iter_sampling = 100L,
    iter_warmup = 10L,
    seed = 1L
  )
  expect_identical(many$nkeep, 100L)
  expect_error(cbind(fit, many), "different numbers of kept draws")

  one <- mipca(
    y,
    rank = 2L,
    transform = "identity",
    chains = 1L,
    iter_sampling = 1L,
    iter_warmup = 10L,
    seed = 1L
  )
  joined <- cbind(fit, one)
  expect_identical(colnames(as.matrix(joined)), c(colnames(obj), colnames(y)))
  expect_identical(as.matrix(joined)[, colnames(obj)], as.matrix(fit))
  expect_identical(as.matrix(joined)[, colnames(y)], as.matrix(one))
})

test_that("output is not a parameter group_imp() or tune_imp() takes", {
  obj <- mipca_input()
  expect_error(
    tune_imp(
      obj,
      parameters = data.frame(ncp = 2, output = "mipca"),
      .f = "pca_imp",
      n_reps = 1
    ),
    "Unknown parameter"
  )
})
