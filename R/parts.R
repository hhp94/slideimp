# Parts: what the imputers pass to each other, and the matrix door that turns
# them into a result.
#
# A part records what one imputer or fallback did to some columns of an input:
# the columns it covers, the holes it filled and the value of each, all in the
# INPUT's numbering. It holds no copy of the input, so passing parts between
# functions copies nothing of matrix size. fill_parts() copies the input once
# and writes every part's values at its holes. A second door, which builds an
# `si_mi` from the same parts, comes later.
#
# The fields, in this order:
# - `cols`: integer, the input columns the part covers, strictly increasing.
#   May be empty.
# - `map`: integer matrix with columns `row` and `col`, one row per hole the
#   part fills, in column-major order (by column, then by row), so each hole
#   appears once. Every `col` is one of `cols`.
# - `values`: double, the value of each hole, in `map` order. None is NA.
# - `kind`: which imputer or fallback filled the part, "pca", "mean" or
#   "median".
# - `info`: list, what that imputer reports about the fill.
# - `left_na`: TRUE when the part covers its columns but fills none of their
#   holes, as `post_imp = FALSE` leaves them. Its `map` is empty.
#
# new_part() only assembles the fields. validate_part() checks them, given the
# input's dimensions, and fill_parts() validates every part it is given, so no
# part reaches a result unchecked. A failed check is a package defect, never a
# user error, so it is a plain stop() with a fixed prefix (dev/WRITING.md,
# section 1).
#
# A part holds one map row and one value per hole, so its checks and its write
# cost memory in proportion to the holes. validate_part() forms the holes'
# linear index once, for its order check, and returns it, and fill_parts()
# writes with it: writing with the two-column map would convert it to a
# linear index again.

part_kinds <- c("pca", "mean", "median")

new_part <- function(
  cols,
  map = cbind(row = integer(0), col = integer(0)),
  values = numeric(0),
  kind,
  info = list(),
  left_na = FALSE
) {
  structure(
    list(
      cols = cols,
      map = map,
      values = values,
      kind = kind,
      info = info,
      left_na = left_na
    ),
    class = "slideimp_part"
  )
}

part_defect <- function(...) {
  stop("slideimp internal, part: ", ..., call. = FALSE)
}

validate_part <- function(part, n_row, n_col) {
  if (!inherits(part, "slideimp_part")) {
    part_defect("not a slideimp_part")
  }
  fields <- c("cols", "map", "values", "kind", "info", "left_na")
  if (!identical(names(part), fields)) {
    part_defect("fields are not ", paste(fields, collapse = ", "))
  }

  cols <- part$cols
  if (!is.integer(cols) || anyNA(cols) || any(cols < 1L | cols > n_col)) {
    part_defect("`cols` is not integer columns of the input")
  }
  if (length(cols) > 1L && any(diff(cols) <= 0L)) {
    part_defect("`cols` is not strictly increasing")
  }

  map <- part$map
  if (
    !is.matrix(map) ||
      !is.integer(map) ||
      !identical(colnames(map), c("row", "col"))
  ) {
    part_defect("`map` is not an integer matrix with columns row and col")
  }
  # min(), max(), tabulate() and is.unsorted() each make one pass over the
  # map and build nothing of its length; the linear index is the one vector
  # of its length these checks keep.
  n_map <- nrow(map)
  lin <- numeric(0)
  if (n_map > 0L) {
    row <- map[, "row"]
    col <- map[, "col"]
    if (anyNA(map) || min(row) < 1L || max(row) > n_row) {
      part_defect("`map` has a row outside the input")
    }
    # `cols` is unique, so its columns' counts sum to the map's length only
    # when every map column is one of them
    if (
      min(col) < 1L ||
        max(col) > n_col ||
        sum(tabulate(col, n_col)[cols]) != n_map
    ) {
      part_defect("`map` has a column outside `cols`")
    }
    # Strictly increasing linear index. It is formed in double, exact below
    # 2^53 cells: in integer, (col - 1) * n_row + row overflows past 2^31 - 1.
    # A hole listed twice fails here too.
    lin <- (col - 1) * as.double(n_row) + row
    if (is.unsorted(lin, strictly = TRUE)) {
      i <- which(diff(lin) <= 0)[[1]] + 1L
      part_defect(
        "`map` is not strictly increasing in column-major order, at map row ",
        i,
        " (row ",
        row[[i]],
        ", col ",
        col[[i]],
        ")"
      )
    }
  }

  values <- part$values
  if (!is.double(values) || length(values) != n_map) {
    part_defect("`values` is not double with one value per map row")
  }
  if (anyNA(values)) {
    part_defect("`values` holds NA")
  }

  kind <- part$kind
  if (!is.character(kind) || length(kind) != 1L || !(kind %in% part_kinds)) {
    part_defect("`kind` is not one of ", paste(part_kinds, collapse = ", "))
  }
  if (!is.list(part$info)) {
    part_defect("`info` is not a list")
  }
  left_na <- part$left_na
  if (!isTRUE(left_na) && !isFALSE(left_na)) {
    part_defect("`left_na` is not TRUE or FALSE")
  }
  if (left_na && n_map > 0L) {
    part_defect("a `left_na` part has a non-empty `map`")
  }
  invisible(lin)
}

# The matrix door. `obj` is the input the parts were built from, a double
# matrix; `parts` is a list of parts over disjoint columns. Columns no part
# covers are returned as they are, holes included, so a caller that means to
# fill every column passes parts that cover every column. In each part that is
# not `left_na`, the map must be exactly the holes of the part's columns:
# every cell it writes is a hole (a value over an observed cell is a defect),
# and no hole in its columns is left out.
#
# The result is `obj` with the values written, its attributes untouched; the
# caller sets the class. It is one copy of `obj`, made by the first write, and
# no copy at all when no part writes anything.
fill_parts <- function(obj, parts) {
  if (!is.matrix(obj) || !is.double(obj)) {
    part_defect("the input is not a double matrix")
  }
  if (!is.list(parts) || inherits(parts, "slideimp_part")) {
    part_defect("`parts` is not a list of parts")
  }
  n_row <- nrow(obj)
  n_col <- ncol(obj)
  lins <- lapply(parts, validate_part, n_row = n_row, n_col = n_col)

  all_cols <- unlist(lapply(parts, `[[`, "cols"), use.names = FALSE)
  if (anyDuplicated(all_cols)) {
    part_defect(
      "two parts cover column ",
      all_cols[[anyDuplicated(all_cols)]]
    )
  }

  # Holes per column, counted by C++ reading the input in place. A part that
  # is not `left_na` lists as many cells as its columns have holes.
  n_miss <- col_miss_internal(obj)
  filling <- !vapply(parts, `[[`, logical(1), "left_na")
  for (i in which(filling)) {
    part <- parts[[i]]
    if (sum(n_miss[part$cols]) != nrow(part$map)) {
      part_defect(
        "fills ",
        nrow(part$map),
        " of the ",
        sum(n_miss[part$cols]),
        " holes in its columns"
      )
    }
  }

  out <- obj
  wrote <- FALSE
  for (i in seq_along(parts)) {
    if (length(lins[[i]])) {
      out[lins[[i]]] <- parts[[i]]$values
      wrote <- TRUE
    }
  }

  # Each filling part lists as many distinct cells of its columns as they have
  # holes, and writes no NA. So a hole still open in those columns means one
  # of its cells was not a hole: a value went over an observed cell. Counting
  # after the write costs a pass over the result and nothing the size of the
  # map, which reading the input at every cell before the write would.
  fill_cols <- unlist(lapply(parts[filling], `[[`, "cols"), use.names = FALSE)
  if (wrote && any(col_miss_internal(out)[fill_cols] > 0)) {
    for (i in which(filling)) {
      observed <- which(!is.na(obj[lins[[i]]]))
      if (length(observed)) {
        cell <- parts[[i]]$map[observed[[1]], ]
        part_defect(
          "writes over an observed cell (row ",
          cell[["row"]],
          ", col ",
          cell[["col"]],
          ")"
        )
      }
    }
    part_defect("leaves a hole open in a column it fills")
  }
  out
}
