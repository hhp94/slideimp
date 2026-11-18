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
  b <- new_part(cols = 4:6, kind = "mean", left_na = TRUE)
  expect_true(anyNA(m[, 4:6]))
  out <- fill_parts(m, list(a, b))
  expect_false(anyNA(out[, 1:3]))
  expect_identical(out[, 4:6], m[, 4:6])
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
    "fills 2 of the 3 holes in its columns",
    fixed = TRUE
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
  refuse(p, "`map` is not an integer matrix")

  p <- a
  colnames(p$map) <- c("i", "j")
  refuse(p, "`map` is not an integer matrix")

  p <- a
  p$map[1, "row"] <- 6L
  refuse(p, "`map` has a row outside the input")

  p <- a
  p$cols <- c(1L, 2L)
  refuse(p, "`map` has a column outside `cols`")

  p <- a
  p$cols <- c(2L, 1L, 4L)
  refuse(p, "`cols` is not strictly increasing")

  p <- a
  p$cols <- c(1, 2, 4)
  refuse(p, "`cols` is not integer columns of the input")

  p <- a
  p$cols <- c(1L, 2L, 7L)
  refuse(p, "`cols` is not integer columns of the input")

  p <- a
  p$values[2] <- NA
  refuse(p, "`values` holds NA")

  p <- a
  p$values <- p$values[-1]
  refuse(p, "`values` is not double with one value per map row")

  p <- a
  p$values <- as.integer(p$values)
  refuse(p, "`values` is not double with one value per map row")

  p <- a
  p$kind <- "knn"
  refuse(p, "`kind` is not one of pca, mean, median")

  p <- a
  p$info <- "note"
  refuse(p, "`info` is not a list")

  p <- a
  p$left_na <- NA
  refuse(p, "`left_na` is not TRUE or FALSE")

  p <- a
  p$left_na <- TRUE
  refuse(p, "a `left_na` part has a non-empty `map`")

  p <- unclass(a)
  refuse(p, "not a slideimp_part")

  # the door's own arguments
  expect_error(
    fill_parts(m, a),
    "`parts` is not a list of parts",
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

# the si_mi door ----

to_si_mi <- function(m, parts, prov_extra = NULL, keep = NULL) {
  parts_to_si_mi(
    m,
    parts,
    call = quote(probe(m)),
    fn = "probe",
    prov_extra = prov_extra,
    keep = keep
  )
}

test_that("parts_to_si_mi() completes to what fill_parts() writes", {
  m <- part_input()
  a <- filled_part(m, c(1L, 2L, 5L), "pca", 100)
  b <- filled_part(m, c(3L, 4L, 6L), "mean", 200)
  expect_gt(nrow(a$map), 0L)
  expect_gt(nrow(b$map), 0L)
  filled <- fill_parts(m, list(a, b))

  fit <- to_si_mi(m, list(a, b))
  expect_s3_class(fit, "si_mi")
  expect_identical(as.matrix(fit), filled)
  # interleaved columns need a placement
  expect_identical(fit$placement, c(1L, 2L, 5L, 3L, 4L, 6L))
  # parts given out of column order
  expect_identical(as.matrix(to_si_mi(m, list(b, a))), filled)

  # each chunk: its slice without row names, its map in the slice's own
  # numbering as an integer matrix, the values as the one draw
  ch <- fit$chunks[[2]]
  slice <- m[, c(3, 4, 6)]
  rownames(slice) <- NULL
  expect_identical(ch$obj, slice)
  expect_true(is.integer(ch$map))
  # holes at (2, 4), (5, 4) and (3, 6): input columns 4 and 6 are the
  # slice's 2 and 3
  expect_identical(unname(ch$map[, "col"]), c(2L, 2L, 3L))
  expect_identical(ch$draws, matrix(b$values, ncol = 1L))
  expect_identical(ch$prov$seam_state, "never")
  expect_identical(ch$prov$thin_effective, 1L)
  expect_identical(nrow(ch$diag$cells), 0L)
  expect_null(ch$info)
  # the call is rewritten to the caller's name, the matrix not stored
  expect_identical(fit$call[[1]], as.name("probe"))
})

test_that("parts_to_si_mi() takes one part, a part over no columns, and left_na", {
  m <- part_input()

  # one part over every column: no slice, the input itself
  every <- filled_part(m, 1:6)
  fit <- to_si_mi(m, list(every))
  expect_length(fit$chunks, 1L)
  expect_identical(fit$chunks[[1]]$obj, m)
  expect_null(fit$placement)
  expect_identical(as.matrix(fit), fill_parts(m, list(every)))
  # an input with more than dims and dimnames: the chunk holds it without
  # them, since as.matrix() starts from the chunk's matrix on this path. The
  # input itself keeps them.
  tagged <- m
  attr(tagged, "foo") <- "bar"
  class(tagged) <- c("slideimp_results", class(m))
  fit <- to_si_mi(tagged, list(every))
  expect_length(fit$chunks, 1L)
  expect_identical(fit$chunks[[1]]$obj, m)
  expect_identical(as.matrix(fit), fill_parts(m, list(every)))
  expect_identical(attr(tagged, "foo"), "bar")
  expect_s3_class(tagged, "slideimp_results")

  # a part over no columns is dropped
  none <- new_part(cols = integer(0), kind = "mean")
  expect_length(to_si_mi(m, list(none, every))$chunks, 1L)

  # a left_na part: no map, no draws, its holes stay NA
  a <- filled_part(m, 1:3)
  b <- new_part(cols = 4:6, kind = "mean", left_na = TRUE)
  expect_true(anyNA(m[, 4:6]))
  fit <- to_si_mi(m, list(a, b))
  ch <- fit$chunks[[2]]
  expect_true(isTRUE(ch$prov$left_na))
  expect_identical(dim(ch$draws), c(0L, 1L))
  expect_null(fit$chunks[[1]]$prov$left_na)
  expect_identical(as.matrix(fit), fill_parts(m, list(a, b)))

  # info and the imputer's prov fields reach their chunk
  a$info <- list(note = "a")
  fit <- to_si_mi(m, list(a, b), prov_extra = list(list(rank = 2L), NULL))
  expect_identical(fit$chunks[[1]]$info, list(note = "a"))
  expect_identical(fit$chunks[[1]]$prov$rank, 2L)
  expect_null(fit$chunks[[2]]$info)
})

test_that("parts_to_si_mi() refuses what it cannot build", {
  m <- part_input()
  a <- filled_part(m, 1:3)
  b <- filled_part(m, 4:6)

  expect_error(
    to_si_mi(m, list(a)),
    "the parts cover 3 of the 6 columns",
    fixed = TRUE
  )
  expect_error(
    to_si_mi(m, list(a, filled_part(m, 3:6))),
    "two parts cover column 3",
    fixed = TRUE
  )
  # the observed-cell swap fill_parts() is tested with: the hole count still
  # matches, so only the check at birth catches it
  over <- filled_part(m, 1L)
  over$map[2, "row"] <- 5L
  rest <- filled_part(m, 2:6)
  expect_error(
    to_si_mi(m, list(over, rest)),
    "writes over an observed cell (row 5, col 1)",
    fixed = TRUE
  )
  expect_error(
    to_si_mi(unname(m), list(a, b)),
    "the si_mi door needs unique row and column names",
    fixed = TRUE
  )
  expect_error(
    to_si_mi(m, list(a, b), prov_extra = list(list(left_na = TRUE), NULL)),
    "`prov_extra` is not one named list",
    fixed = TRUE
  )
  expect_error(
    to_si_mi(m, list(a, b), prov_extra = list(NULL)),
    "`prov_extra` is not one named list",
    fixed = TRUE
  )
  expect_error(
    to_si_mi(as.character(m), list(a, b)),
    "the input is not a numeric matrix",
    fixed = TRUE
  )
})

# keep: the doors narrowed to some columns ----

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
  la <- filled_part(m, 1:3)
  lb <- new_part(cols = 4:6, kind = "mean", left_na = TRUE)
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
    "fills 1 of the 2 holes in its columns",
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

test_that("both doors refuse a malformed keep", {
  m <- part_input()
  a <- filled_part(m, 1:3)
  b <- filled_part(m, 4:6)
  bad <- list(
    "`keep` is not strictly increasing" = c(4L, 2L),
    "`keep` is not strictly increasing" = c(2L, 2L),
    "`keep` is not one or more integer columns of the input" = c(2, 4),
    "`keep` is not one or more integer columns of the input" = 7L,
    "`keep` is not one or more integer columns of the input" = 0L,
    "`keep` is not one or more integer columns of the input" = c(2L, NA),
    "`keep` is not one or more integer columns of the input" = integer(0)
  )
  for (i in seq_along(bad)) {
    expect_error(
      fill_parts(m, list(a, b), keep = bad[[i]]),
      names(bad)[[i]],
      fixed = TRUE
    )
    expect_error(
      to_si_mi(m, list(a, b), keep = bad[[i]]),
      names(bad)[[i]],
      fixed = TRUE
    )
  }
})

test_that("parts_to_si_mi() with keep completes to the narrowed fill", {
  m <- part_input()
  a <- filled_part(m, c(1L, 2L, 5L), "pca", 100)
  b <- filled_part(m, c(3L, 4L, 6L), "mean", 200)
  full <- fill_parts(m, list(a, b))

  for (nm in names(keep_cases)) {
    keep <- keep_cases[[nm]]
    expect_identical(
      as.matrix(to_si_mi(m, list(a, b), keep = keep)),
      full[, keep, drop = FALSE],
      label = nm
    )
    expect_identical(
      as.matrix(to_si_mi(m, list(b, a), keep = keep)),
      full[, keep, drop = FALSE],
      label = nm
    )
  }

  # interleaved: a keeps input columns 1 and 5, b keeps 3 and 6, which are
  # the fit's columns 1, 3 and 2, 4
  fit <- to_si_mi(m, list(a, b), keep = keep_cases$interleaved)
  expect_length(fit$chunks, 2L)
  expect_identical(fit$placement, c(1L, 3L, 2L, 4L))
  ch <- fit$chunks[[2]]
  slice <- m[, c(3, 6)]
  rownames(slice) <- NULL
  expect_identical(ch$obj, slice)
  # b's kept hole is (3, 6), the slice's column 2; its holes in column 4 are
  # cut away with their values
  expect_identical(b$map[3, ], c(row = 3L, col = 6L))
  expect_identical(unname(ch$map[, "col"]), 2L)
  expect_true(is.integer(ch$map))
  expect_identical(ch$draws, matrix(b$values[[3]], ncol = 1L))

  # a part wholly outside keep is dropped
  expect_length(
    to_si_mi(m, list(a, b), keep = keep_cases$one_part)$chunks,
    1L
  )

  # one part over every column, narrowed: a slice, not the input
  every <- filled_part(m, 1:6)
  fit <- to_si_mi(m, list(every), keep = c(2L, 4L))
  expect_length(fit$chunks, 1L)
  expect_identical(colnames(fit$chunks[[1]]$obj), c("c2", "c4"))
  expect_identical(
    as.matrix(fit),
    fill_parts(m, list(every))[, c(2L, 4L)]
  )
  # every column is the same as none: the input itself, no slice
  fit <- to_si_mi(m, list(every), keep = 1:6)
  expect_identical(fit$chunks[[1]]$obj, m)

  # a left_na part, cut
  la <- filled_part(m, 1:3)
  lb <- new_part(cols = 4:6, kind = "mean", left_na = TRUE)
  fit <- to_si_mi(m, list(la, lb), keep = c(2L, 4L))
  expect_true(isTRUE(fit$chunks[[2]]$prov$left_na))
  expect_identical(
    as.matrix(fit),
    fill_parts(m, list(la, lb), keep = c(2L, 4L))
  )

  # the parts must cover every kept column
  expect_error(
    to_si_mi(m, list(la), keep = c(2L, 4L)),
    "the parts cover 1 of the 2 kept columns",
    fixed = TRUE
  )
})
