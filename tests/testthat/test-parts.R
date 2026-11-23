skip_on_cran()

# parts and the matrix door ----

# A 5 x 6 input with holes in columns 1, 2, 4 and 6; columns 3 and 5 are
# complete.
part_input <- function() {
  m <- matrix(
    as.numeric(1:30),
    5,
    6,
    dimnames = list(paste0("r", 1:5), paste0("c", 1:6))
  )
  m[cbind(c(1, 4, 2, 2, 5, 3), c(1, 1, 2, 4, 4, 6))] <- NA
  m
}

# The holes of `m` in `cols`, as a part map: which() walks the matrix in
# column-major order, which is the order a map must have.
holes_map <- function(m, cols) {
  idx <- which(is.na(m), arr.ind = TRUE)
  idx <- idx[idx[, "col"] %in% cols, , drop = FALSE]
  cbind(row = as.integer(idx[, "row"]), col = as.integer(idx[, "col"]))
}

# A part that fills every hole in `cols`, with values from `start` upward so
# that each written cell is told apart.
filled_part <- function(m, cols, kind = "pca", start = 100) {
  map <- holes_map(m, cols)
  new_part(
    cols = as.integer(cols),
    map = map,
    values = start + seq_len(nrow(map)),
    kind = kind
  )
}

test_that("fill_parts() writes each part's values at its holes", {
  m <- part_input()
  m0 <- m
  a <- filled_part(m, 1:3, "pca", 100)
  b <- filled_part(m, 4:6, "mean", 200)
  expect_gt(nrow(a$map), 0L)
  expect_gt(nrow(b$map), 0L)

  expected <- m
  expected[a$map] <- a$values
  expected[b$map] <- b$values
  out <- fill_parts(m, list(a, b))
  expect_identical(out, expected)
  expect_false(anyNA(out))
  expect_identical(dimnames(out), dimnames(m))
  # observed cells untouched
  expect_identical(out[!is.na(m0)], m0[!is.na(m0)])
  # the caller's input is not written to
  expect_identical(m, m0)

  # parts given out of column order
  expect_identical(fill_parts(m, list(b, a)), expected)
})

test_that("fill_parts() takes a part over no columns, one, or every column", {
  m <- part_input()

  # no columns: writes nothing, alone or beside another part
  none <- new_part(cols = integer(0), kind = "mean")
  expect_identical(fill_parts(m, list(none)), m)
  every <- filled_part(m, 1:6)
  expect_identical(
    fill_parts(m, list(none, every)),
    fill_parts(m, list(every))
  )

  # one column: the columns no part covers keep their holes
  one <- filled_part(m, 4L)
  expect_gt(nrow(one$map), 0L)
  out <- fill_parts(m, list(one))
  expect_false(anyNA(out[, 4]))
  expect_identical(out[, -4], m[, -4])
  expect_true(anyNA(out[, -4]))

  # every column
  expect_false(anyNA(fill_parts(m, list(every))))

  # a part over complete columns only: an empty map, and nothing to fill
  done <- filled_part(m, c(3L, 5L))
  expect_identical(nrow(done$map), 0L)
  expect_identical(fill_parts(m, list(done)), m)

  # no parts at all
  expect_identical(fill_parts(m, list()), m)
})

test_that("fill_parts() leaves the holes of a left_na part NA", {
  m <- part_input()
  a <- filled_part(m, 1:3)
  b <- new_part(cols = c(4L, 6L), kind = "mean", left_na = TRUE)
  done <- new_part(cols = 5L, kind = "mean")
  expect_true(anyNA(m[, 4:6]))
  out <- fill_parts(m, list(a, b, done))
  expect_false(anyNA(out[, 1:3]))
  expect_identical(out[, 4:6], m[, 4:6])
})

# A left_na part whose map names some of the holes of its columns, as K-NN
# leaves a cell no neighbor is observed at: columns 1 and 4 have two holes
# each, and the map names one of each, (4, 1) and (5, 4).
left_part <- function(start = 300) {
  new_part(
    cols = c(1L, 4L),
    map = cbind(row = c(4L, 5L), col = c(1L, 4L)),
    values = start + 1:2,
    kind = "knn",
    left_na = TRUE
  )
}

test_that("a left_na part with a map writes its cells and leaves the rest", {
  m <- part_input()
  lp <- left_part()
  rest <- filled_part(m, c(2L, 3L, 5L, 6L), "knn", 100)
  expect_identical(colSums(is.na(m))[c(1, 4)], c(c1 = 2, c4 = 2))

  out <- fill_parts(m, list(lp, rest))
  expect_identical(out[4, 1], 301)
  expect_identical(out[5, 4], 302)
  # the holes the map does not name stay NA, and only they
  expect_identical(
    unname(which(is.na(out), arr.ind = TRUE)),
    cbind(c(1L, 2L), c(1L, 4L))
  )

  # cut by keep
  expect_identical(
    fill_parts(m, list(lp, rest), keep = c(2L, 4L)),
    out[, c(2L, 4L)]
  )
})

test_that("hole counts past the integer range do not overflow", {
  m <- part_input()
  # col_scan() counts in integer, and the door sums the counts over a
  # part's columns. A matrix with that many holes is too large to build here,
  # so the count is raised for columns 1 and 2 instead, by the same amount
  # before and after the write: each then counts .Machine$integer.max - 1 or
  # more, and the two sum past the integer range.
  real_scan <- col_scan
  big <- .Machine$integer.max - 2L
  local_mocked_bindings(
    col_scan = function(obj, ...) {
      scan <- real_scan(obj, ...)
      scan$n_miss[1:2] <- scan$n_miss[1:2] + big
      scan
    }
  )
  n_miss <- col_scan(m)$n_miss
  expect_identical(n_miss[1:2], c(big + 2L, big + 1L))
  expect_gt(sum(as.double(n_miss[1:2])), .Machine$integer.max)

  # a left_na part that names one hole of each column: a valid part, and the
  # write goes through
  lp <- new_part(
    cols = 1:2,
    map = cbind(row = c(4L, 2L), col = 1:2),
    values = c(301, 302),
    kind = "knn",
    left_na = TRUE
  )
  out <- fill_parts(m, list(lp))
  expect_identical(out[4, 1], 301)
  expect_identical(out[2, 2], 302)
  expect_identical(out[, 3:6], m[, 3:6])

  # a filling part that names 3 of that many holes is the count defect
  expect_error(
    fill_parts(m, list(filled_part(m, 1:2))),
    "one row per hole in its columns (rows 3, holes 4294967293",
    fixed = TRUE
  )
})

test_that("fill_parts() refuses a left_na part that keeps no hole in a column", {
  m <- part_input()
  rest <- filled_part(m, c(3L, 5L, 6L))
  # column 2 has one hole, and the map names it
  whole_col <- new_part(
    cols = c(1L, 2L, 4L),
    map = cbind(row = c(4L, 2L, 5L), col = c(1L, 2L, 4L)),
    values = c(1, 2, 3),
    kind = "knn",
    left_na = TRUE
  )
  expect_error(
    fill_parts(m, list(whole_col, rest)),
    "a 'left_na' part keeps no hole in column 2",
    fixed = TRUE
  )
  # a complete column has no hole to keep, even under an empty map
  complete_col <- new_part(cols = c(4L, 5L), kind = "mean", left_na = TRUE)
  expect_false(anyNA(m[, 5]))
  expect_error(
    fill_parts(m, list(complete_col)),
    "a 'left_na' part keeps no hole in column 5",
    fixed = TRUE
  )
})

test_that("fill_parts() refuses a left_na map cell over an observed value", {
  m <- part_input()
  lp <- left_part()
  # (4, 1) is a hole; (5, 1) is observed. Column 1 still keeps a hole, so
  # only the observed-cell checks can catch it.
  expect_false(is.na(m[5, 1]))
  lp$map[1, "row"] <- 5L
  rest <- filled_part(m, c(2L, 3L, 5L, 6L))
  expect_error(
    fill_parts(m, list(lp, rest)),
    "writes over an observed cell (row 5, col 1)",
    fixed = TRUE
  )
})

test_that("fill_parts() refuses a hole written twice", {
  m <- part_input()

  # twice in one part: the map lists a hole twice
  a <- filled_part(m, 1:3)
  dup <- a
  dup$map <- a$map[c(1, 1:nrow(a$map)), ]
  dup$values <- c(0, a$values)
  expect_identical(dup$map[1, ], dup$map[2, ])
  expect_error(
    fill_parts(m, list(dup)),
    "not strictly increasing in column-major order",
    fixed = TRUE
  )

  # once in each of two parts over the same column
  b <- filled_part(m, 2:4)
  expect_true(anyNA(m[, 2]))
  expect_error(
    fill_parts(m, list(a, b)),
    "two parts cover column 2",
    fixed = TRUE
  )
})

test_that("fill_parts() refuses a write over an observed cell", {
  m <- part_input()
  a <- filled_part(m, 1L)
  # swap the hole at (4, 1) for the observed cell (5, 1): the hole count
  # still matches, so only the observed-cell check can catch it
  expect_identical(a$map[2, ], c(row = 4L, col = 1L))
  expect_false(is.na(m[5, 1]))
  a$map[2, "row"] <- 5L
  expect_error(
    fill_parts(m, list(a)),
    "writes over an observed cell (row 5, col 1)",
    fixed = TRUE
  )
})

test_that("fill_parts() refuses a part that leaves a hole in its columns", {
  m <- part_input()
  a <- filled_part(m, 1:2)
  short <- new_part(
    cols = a$cols,
    map = a$map[-1, , drop = FALSE],
    values = a$values[-1],
    kind = "pca"
  )
  expect_error(
    fill_parts(m, list(short)),
    "one row per hole in its columns (rows 2, holes 3",
    fixed = TRUE
  )

  # and one whose map has more rows than its columns have holes: column 1
  # has holes at rows 1 and 4, and the map adds the observed cell (2, 1)
  long <- new_part(
    cols = 1L,
    map = cbind(row = c(1L, 2L, 4L), col = 1L),
    values = c(1, 2, 3),
    kind = "pca"
  )
  expect_false(is.na(m[2, 1]))
  expect_error(
    fill_parts(m, list(long)),
    "one row per hole in its columns (rows 3, holes 2",
    fixed = TRUE
  )
})

test_that("a part defect names the part's position in the list", {
  m <- part_input()
  a <- filled_part(m, 1:3)
  b <- filled_part(m, 4:6)
  expect_identical(nrow(b$map), 3L)

  # a field check, from validate_part()
  bad <- b
  bad$left_na <- NA
  expect_error(
    fill_parts(m, list(a, bad)),
    "'left_na' is not TRUE or FALSE, in part 2",
    fixed = TRUE
  )
  expect_error(
    fill_parts(m, list(bad, a)),
    "'left_na' is not TRUE or FALSE, in part 1",
    fixed = TRUE
  )

  # the count check
  short <- b
  short$map <- b$map[-1, , drop = FALSE]
  short$values <- b$values[-1]
  expect_error(
    fill_parts(m, list(a, short)),
    "(rows 2, holes 3), in part 2",
    fixed = TRUE
  )

  # the check after the write: the hole (3, 6) swapped for the observed
  # cell (4, 6)
  over <- b
  expect_identical(over$map[3, ], c(row = 3L, col = 6L))
  expect_false(is.na(m[4, 6]))
  over$map[3, "row"] <- 4L
  expect_error(
    fill_parts(m, list(a, over)),
    "writes over an observed cell (row 4, col 6), in part 2",
    fixed = TRUE
  )

  # validate_part() called with no position names none
  msg <- tryCatch(
    validate_part(bad, nrow(m), ncol(m)),
    error = conditionMessage
  )
  expect_identical(
    msg,
    "slideimp internal, part: 'left_na' is not TRUE or FALSE"
  )
})

test_that("validate_part() refuses malformed parts", {
  m <- part_input()
  a <- filled_part(m, c(1L, 2L, 4L))
  expect_gte(nrow(a$map), 3L)
  refuse <- function(part, msg) {
    expect_error(fill_parts(m, list(part)), msg, fixed = TRUE)
  }

  p <- a
  p$map <- a$map[c(2, 1, 3:nrow(a$map)), ]
  refuse(p, "not strictly increasing in column-major order, at map row 2")

  p <- a
  storage.mode(p$map) <- "double"
  refuse(p, "'map' is not an integer matrix")

  p <- a
  colnames(p$map) <- c("i", "j")
  refuse(p, "'map' is not an integer matrix")

  p <- a
  p$map[1, "row"] <- 6L
  refuse(p, "'map' has a row outside the input")

  p <- a
  p$cols <- c(1L, 2L)
  refuse(p, "'map' has a column outside 'cols'")

  p <- a
  p$cols <- c(2L, 1L, 4L)
  refuse(p, "'cols' is not strictly increasing")

  p <- a
  p$cols <- c(1, 2, 4)
  refuse(p, "'cols' is not integer columns of the input")

  p <- a
  p$cols <- c(1L, 2L, 7L)
  refuse(p, "'cols' is not integer columns of the input")

  p <- a
  p$values[2] <- NA
  refuse(p, "'values' holds NA")

  # Inf is no hole to the counts, so only this check can catch it
  p <- a
  p$values[2] <- Inf
  refuse(p, "'values' holds an infinite value")

  p <- a
  p$values[2] <- -Inf
  refuse(p, "'values' holds an infinite value")

  p <- a
  p$values <- p$values[-1]
  refuse(p, "'values' is not double with one value per map row")

  p <- a
  p$values <- as.integer(p$values)
  refuse(p, "'values' is not double with one value per map row")

  p <- a
  p$kind <- "kmeans"
  refuse(p, "'kind' is not one of \"pca\", \"knn\", \"mean\", \"median\"")

  p <- a
  p$info <- "note"
  refuse(p, "'info' is not a list")

  p <- a
  p$left_na <- NA
  refuse(p, "'left_na' is not TRUE or FALSE")

  # a's map names every hole of its columns, so none keeps one
  p <- a
  p$left_na <- TRUE
  refuse(p, "a 'left_na' part keeps no hole in column 1")

  p <- unclass(a)
  refuse(p, "not a \"slideimp_part\"")

  p <- a
  p$left_na <- NULL
  refuse(p, "fields are not 'cols', 'map', 'values'")

  # the door's own arguments
  expect_error(
    fill_parts(m, a),
    "'parts' is not a list of parts",
    fixed = TRUE
  )
  mi <- m
  storage.mode(mi) <- "integer"
  expect_error(
    fill_parts(mi, list(a)),
    "the input is not a double matrix",
    fixed = TRUE
  )
})

# keep: the door narrowed to some columns ----

# The two parts interleave: a over columns 1, 2, 5 and b over 3, 4, 6. Each
# `keep` below is checked against the full fill cut to it afterwards.
keep_cases <- list(
  across_both = c(2L, 4L), # each part crosses the edge of keep
  interleaved = c(1L, 3L, 5L, 6L),
  one_part = c(1L, 2L, 5L), # a wholly inside, b wholly outside
  one_col = 4L,
  complete_only = c(3L, 5L) # no hole kept
)

test_that("fill_parts() with keep is the full fill cut to keep", {
  m <- part_input()
  m0 <- m
  a <- filled_part(m, c(1L, 2L, 5L), "pca", 100)
  b <- filled_part(m, c(3L, 4L, 6L), "mean", 200)
  full <- fill_parts(m, list(a, b))
  expect_true(anyNA(m[, keep_cases$across_both]))
  expect_false(anyNA(m[, keep_cases$complete_only]))

  for (nm in names(keep_cases)) {
    keep <- keep_cases[[nm]]
    expect_identical(
      fill_parts(m, list(a, b), keep = keep),
      full[, keep, drop = FALSE],
      label = nm
    )
    expect_identical(
      fill_parts(m, list(b, a), keep = keep),
      full[, keep, drop = FALSE],
      label = nm
    )
  }
  # every column is the same as none: no cut, the input's own copy
  expect_identical(fill_parts(m, list(a, b), keep = 1:6), full)
  expect_identical(m, m0)

  # a left_na part keeps its holes in the kept columns
  la <- filled_part(m, c(1L, 2L, 3L, 5L))
  lb <- new_part(cols = c(4L, 6L), kind = "mean", left_na = TRUE)
  out <- fill_parts(m, list(la, lb), keep = c(2L, 4L))
  expect_identical(colnames(out), c("c2", "c4"))
  expect_false(anyNA(out[, "c2"]))
  expect_true(anyNA(out[, "c4"]))
  expect_identical(out[, "c4"], m[, "c4"])

  # a kept column no part covers comes back as it is, as without keep
  out <- fill_parts(m, list(la), keep = c(2L, 4L))
  expect_identical(out[, "c4"], m[, "c4"])
})

test_that("fill_parts() checks the kept columns of a cut part", {
  m <- part_input()

  # a part over columns 1 and 2 that leaves out the first hole of column 1:
  # refused when column 1 is kept, cut away with it when it is not
  a <- filled_part(m, 1:2)
  short <- new_part(
    cols = a$cols,
    map = a$map[-1, , drop = FALSE],
    values = a$values[-1],
    kind = "pca"
  )
  expect_identical(unname(a$map[1, ]), c(1L, 1L))
  expect_error(
    fill_parts(m, list(short), keep = c(1L, 3L)),
    "one row per hole in its columns (rows 1, holes 2",
    fixed = TRUE
  )
  expect_false(anyNA(fill_parts(m, list(short), keep = 2L)))

  # a value over an observed cell, found in the narrowed copy and reported in
  # the input's numbering
  over <- filled_part(m, 1:2)
  expect_identical(over$map[2, ], c(row = 4L, col = 1L))
  expect_false(is.na(m[5, 1]))
  over$map[2, "row"] <- 5L
  expect_error(
    fill_parts(m, list(over), keep = c(1L, 4L)),
    "writes over an observed cell (row 5, col 1)",
    fixed = TRUE
  )
})

test_that("fill_parts() refuses a malformed keep", {
  m <- part_input()
  a <- filled_part(m, 1:3)
  b <- filled_part(m, 4:6)
  bad <- list(
    "'keep' is not strictly increasing" = c(4L, 2L),
    "'keep' is not strictly increasing" = c(2L, 2L),
    "'keep' is not one or more integer columns of the input" = c(2, 4),
    "'keep' is not one or more integer columns of the input" = 7L,
    "'keep' is not one or more integer columns of the input" = 0L,
    "'keep' is not one or more integer columns of the input" = c(2L, NA),
    "'keep' is not one or more integer columns of the input" = integer(0)
  )
  for (i in seq_along(bad)) {
    expect_error(
      fill_parts(m, list(a, b), keep = bad[[i]]),
      names(bad)[[i]],
      fixed = TRUE
    )
  }
})

test_that("the parts internals take their wrappers' formals", {
  # group_imp() and slide_imp() call the internals with what they would pass
  # the wrappers, so a default that differs between the two would change
  # their groups and windows silently. The internals take one more, the
  # scan, second and with no default.
  for (pair in list(
    list(knn_imp_parts, knn_imp),
    list(pca_imp_parts, pca_imp)
  )) {
    inner <- formals(pair[[1L]])
    expect_identical(names(inner)[2L], "scan")
    expect_identical(inner[["scan"]], quote(expr = ))
    expect_identical(as.list(inner)[-2L], as.list(formals(pair[[2L]])))
  }
})

test_that("rebase_parts() moves a slice's parts into the input's numbering", {
  m <- part_input()
  # the slice takes input columns out of order. It keeps its first three,
  # input columns 6, 2 and 4, so the moved map is out of column order until
  # rebase_parts() sorts it.
  col_idx <- c(6L, 2L, 4L, 1L, 3L)
  keep <- 1:3
  expect_true(is.unsorted(col_idx[keep]))
  s <- m[, col_idx]
  dimnames(s) <- NULL
  # `a` crosses the cut: slice column 4 (input column 1) holds two holes.
  # `b` lies wholly outside it.
  a <- filled_part(s, 1:4, "pca", 100)
  b <- filled_part(s, 5L, "mean", 200)
  expect_gt(sum(a$map[, "col"] == 4L), 0L)

  moved <- rebase_parts(list(a, b), col_idx, keep)
  expect_length(moved, 1L)
  expect_identical(moved[[1]]$cols, c(2L, 4L, 6L))
  expect_true(is.unsorted(col_idx[a$map[a$map[, "col"] <= 3L, "col"]]))

  # filling the moved part into the input is filling the slice and writing
  # its kept columns back
  expected <- m
  expected[, col_idx[keep]] <- fill_parts(s, list(a, b))[, keep]
  expect_identical(fill_parts(m, moved), expected)

  # without `keep`, every slice column moves
  all_moved <- rebase_parts(list(a, b), col_idx)
  expect_identical(
    lapply(all_moved, `[[`, "cols"),
    list(c(1L, 2L, 4L, 6L), 3L)
  )
  expected <- m
  expected[, col_idx] <- fill_parts(s, list(a, b))
  expect_identical(fill_parts(m, all_moved), expected)
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
      "'values' holds NA",
      fixed = TRUE
    )
  }

  # with no dimnames, the same part
  expect_identical(fallback_part(unname(m), cols), fallback_part(m, cols))
})
