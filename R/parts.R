# Parts: what the imputers pass to each other, and the matrix door that turns
# them into a result.
#
# A part records what one imputer or fallback did to some columns of an input:
# the columns it covers, the holes it filled and the value of each, all in the
# INPUT's numbering. It holds no copy of the input, so passing parts between
# functions copies nothing of matrix size. fill_parts() copies the input once
# and writes every part's values at its holes.
#
# The fields, in this order:
# - `cols`: integer, the input columns the part covers, strictly increasing.
#   May be empty.
# - `map`: integer matrix with columns `row` and `col`, one row per hole the
#   part fills, in column-major order (by column, then by row), so each hole
#   appears once. Every `col` is one of `cols`.
# - `values`: double, the value of each hole, in `map` order. Each is finite:
#   none is NA, and none is infinite, as the imputers refuse an input that
#   holds an infinite value.
# - `kind`: which imputer or fallback filled the part, "pca", "knn", "mean"
#   or "median".
# - `info`: list, what that imputer reports about the fill.
# - `left_na`: TRUE when every column of the part keeps a hole: a hole in its
#   columns that `map` does not name stays NA. Under `post_imp = FALSE` the
#   columns `colmax` holds out are such a part with an empty map, and the
#   K-NN columns with a cell no neighbor could fill are one whose map names
#   the cells K-NN did fill. FALSE when `map` names every hole of its columns.
#
# fallback_part() builds the part that sets each hole of some columns to its
# column's mean: the imputers' `post_imp` fill of the columns they leave out,
# and the `on_infeasible = "mean"` fill of group_imp() and slide_imp().
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
#
# The door takes `keep`, the input columns to return, so a call given
# `subset` returns only those columns. It cuts the parts to `keep` after
# validating them, so an imputer may hand over a part that also covers columns
# outside it, as the PCA kernel's part does.
#
# knn_imp() and pca_imp() are each a door that takes the scan (R/col_scan.R),
# an internal that returns parts (knn_imp_parts(), pca_imp_parts()), and
# fill_imp_parts(). group_imp() calls the internals on each group's slice,
# moves the parts into its input's numbering with rebase_parts(), and fills
# them all in one fill_parts().
# slide_imp() calls them on each window; its windows overlap, so it averages
# their values at the holes and fills the averages, as parts, in one
# fill_parts() (R/slide_imp.R: `slide_parts <- function(`).

part_kinds <- c("pca", "knn", "mean", "median")

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

# A fallback part over `cols` of `obj`: each hole in those columns set to the
# mean or the median of its column's observed values, and `kind` set to
# `method`. `cols` is integer and strictly increasing, and none of its columns
# is without an observed value: such a column has neither statistic, and a
# part over one fails validation, its values holding NA.
#
# Vectorized over the columns: their submatrix is copied out once, with its
# is.na(). So its callers (the imputers, group_imp() and slide_imp()) pass
# only columns with a hole: a complete column would be copied to fill
# nothing. The fallback covers few columns beside those the imputer filled,
# so a copy of them is small. The mean is .colMeans(), a long double sum
# divided in long double. The median is collapse::fmedian(), the mean of the
# two middle values for an even count, as stats::median() takes it.
fallback_part <- function(obj, cols, method = c("mean", "median")) {
  method <- match.arg(method)
  if (length(cols) == 0L) {
    return(new_part(cols = cols, kind = method))
  }
  sub <- obj[, cols, drop = FALSE]
  map <- which(is.na(sub), arr.ind = TRUE, useNames = FALSE)
  colnames(map) <- c("row", "col")
  stat <- switch(
    method,
    mean = .colMeans(sub, nrow(sub), ncol(sub), na.rm = TRUE),
    median = collapse::fmedian(sub, na.rm = TRUE, use.g.names = FALSE)
  )
  values <- as.vector(stat)[map[, "col"]]
  map[, "col"] <- cols[map[, "col"]]
  new_part(cols = cols, map = map, values = values, kind = method)
}

# `at`, when given, is the position of the part at fault in the list the
# door was handed, appended after the wording, so that a defect in one of
# several parts can be traced to it.
part_defect <- function(..., at = NULL) {
  if (is.null(at)) {
    stop("slideimp internal, part: ", ..., call. = FALSE)
  }
  stop("slideimp internal, part: ", ..., ", in part ", at, call. = FALSE)
}

validate_part <- function(part, n_row, n_col, at = NULL) {
  defect <- function(...) part_defect(..., at = at)
  if (!inherits(part, "slideimp_part")) {
    defect("not a \"slideimp_part\"")
  }
  fields <- c("cols", "map", "values", "kind", "info", "left_na")
  if (!identical(names(part), fields)) {
    defect("fields are not ", paste0("'", fields, "'", collapse = ", "))
  }

  cols <- part$cols
  if (!is.integer(cols) || anyNA(cols) || any(cols < 1L | cols > n_col)) {
    defect("'cols' is not integer columns of the input")
  }
  if (length(cols) > 1L && any(diff(cols) <= 0L)) {
    defect("'cols' is not strictly increasing")
  }

  map <- part$map
  if (
    !is.matrix(map) ||
      !is.integer(map) ||
      !identical(colnames(map), c("row", "col"))
  ) {
    defect(
      "'map' is not an integer matrix with columns \"row\" and \"col\""
    )
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
      defect("'map' has a row outside the input")
    }
    # `cols` is unique, so its columns' counts sum to the map's length only
    # when every map column is one of them
    if (
      min(col) < 1L ||
        max(col) > n_col ||
        sum(tabulate(col, n_col)[cols]) != n_map
    ) {
      defect("'map' has a column outside 'cols'")
    }
    # Strictly increasing linear index. It is formed in double, exact below
    # 2^53 cells: in integer, (col - 1) * n_row + row overflows past 2^31 - 1.
    # A hole listed twice fails here too.
    lin <- (col - 1) * as.double(n_row) + row
    if (is.unsorted(lin, strictly = TRUE)) {
      i <- which(diff(lin) <= 0)[[1]] + 1L
      defect(
        "'map' is not strictly increasing in column-major order, at map row ",
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
    defect("'values' is not double with one value per map row")
  }
  if (anyNA(values)) {
    defect("'values' holds NA")
  }
  # Inf is no hole to col_scan(), so the counts after the write would not show
  # one. min() and max() build nothing of the values' length.
  if (n_map > 0L && (min(values) == -Inf || max(values) == Inf)) {
    defect("'values' holds an infinite value")
  }

  kind <- part$kind
  if (!is.character(kind) || length(kind) != 1L || !(kind %in% part_kinds)) {
    defect(
      "'kind' is not one of ",
      paste0("\"", part_kinds, "\"", collapse = ", ")
    )
  }
  if (!is.list(part$info)) {
    defect("'info' is not a list")
  }
  left_na <- part$left_na
  if (!isTRUE(left_na) && !isFALSE(left_na)) {
    defect("'left_na' is not TRUE or FALSE")
  }
  invisible(lin)
}

# The door's `keep`: NULL, or the input columns to return, strictly increasing,
# at least one. Returns NULL for every column, so a `keep` that names every
# column takes the same path as none, with no copy where none is needed.
check_keep <- function(keep, n_col) {
  if (is.null(keep)) {
    return(NULL)
  }
  if (
    !is.integer(keep) ||
      length(keep) == 0L ||
      anyNA(keep) ||
      min(keep) < 1L ||
      max(keep) > n_col
  ) {
    part_defect("'keep' is not one or more integer columns of the input")
  }
  if (length(keep) > 1L && any(diff(keep) <= 0L)) {
    part_defect("'keep' is not strictly increasing")
  }
  if (length(keep) == n_col) NULL else keep
}

# Cuts validated parts to the columns in `keep`, with their linear indices.
# A part wholly inside `keep` is passed on as it is, so the cut allocates only
# for parts that cross its edge: their kept map rows and values. A cut part
# is still valid, since a subset of a strictly increasing map stays so.
cut_parts <- function(parts, lins, keep, n_col) {
  kept <- logical(n_col)
  kept[keep] <- TRUE
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    in_keep <- kept[part$cols]
    if (all(in_keep)) {
      next
    }
    part$cols <- part$cols[in_keep]
    if (nrow(part$map) > 0L) {
      rows <- kept[part$map[, "col"]]
      part$map <- part$map[rows, , drop = FALSE]
      part$values <- part$values[rows]
      lins[[i]] <- lins[[i]][rows]
    }
    parts[[i]] <- part
  }
  list(parts = parts, lins = lins)
}

# The holes in `cols`, from the holes per column. The counts are integer
# (col_scan()), and sum() of integers is NA past 2^31 - 1, which the
# holes of a part's columns can pass in a large input, so they are summed in
# double, exact below 2^53.
sum_holes <- function(n_miss, cols) {
  sum(as.double(n_miss[cols]))
}

# The checks the door makes on its parts, given the input: every part
# valid, no column in two parts, each part that is not `left_na` lists as
# many cells as its columns have holes, and each column of a `left_na` part
# has more holes than its map lists, so it keeps one. With `keep` (from
# check_keep()), the parts are cut to it after the first two checks, so the
# hole counts hold each part to its kept columns only; a defect in a column
# not kept is not looked for, since nothing of that column is returned.
# `n_miss` is the holes of each input column, from the scan the door took;
# NULL counts them here, over every column. Returns the parts, cut or not,
# each part's linear index in the input (from validate_part()), which parts
# fill, and the holes per input column.
check_parts <- function(obj, parts, keep = NULL, n_miss = NULL) {
  if (!is.list(parts) || inherits(parts, "slideimp_part")) {
    part_defect("'parts' is not a list of parts")
  }
  n_row <- nrow(obj)
  n_col <- ncol(obj)
  lins <- lapply(seq_along(parts), function(i) {
    validate_part(parts[[i]], n_row, n_col, at = i)
  })

  all_cols <- unlist(lapply(parts, `[[`, "cols"), use.names = FALSE)
  if (anyDuplicated(all_cols)) {
    part_defect(
      "two parts cover column ",
      all_cols[[anyDuplicated(all_cols)]]
    )
  }

  if (!is.null(keep)) {
    cut <- cut_parts(parts, lins, keep, n_col)
    parts <- cut$parts
    lins <- cut$lins
  }

  # Holes per column, the door's or counted by the scan reading the input
  # in place. A part that is not `left_na` lists as many cells as its columns
  # have holes. A `left_na` part lists fewer than each of its columns has, one
  # tabulate() of its map's columns, so the flag holds of every column the
  # part covers: a column it fills whole belongs in a part that is not
  # `left_na`.
  if (is.null(n_miss)) {
    n_miss <- col_scan(obj)$n_miss
  } else if (!is.integer(n_miss) || length(n_miss) != n_col || anyNA(n_miss)) {
    part_defect("'n_miss' is not one integer count per input column")
  }
  filling <- !vapply(parts, `[[`, logical(1), "left_na")
  for (i in seq_along(parts)) {
    part <- parts[[i]]
    if (filling[[i]]) {
      holes <- sum_holes(n_miss, part$cols)
      if (holes != nrow(part$map)) {
        # The map may have fewer rows than the holes or more, so the wording
        # takes no side. "%.0f": paste() would print a round double count as
        # 1e+05.
        part_defect(
          "'map' does not have one row per hole in its columns (rows ",
          nrow(part$map),
          ", holes ",
          sprintf("%.0f", holes),
          ")",
          at = i
        )
      }
    } else if (length(part$cols)) {
      listed <- tabulate(part$map[, "col"], n_col)[part$cols]
      full <- which(listed >= n_miss[part$cols])
      if (length(full)) {
        part_defect(
          "a 'left_na' part keeps no hole in column ",
          part$cols[[full[[1]]]],
          at = i
        )
      }
    }
  }
  list(parts = parts, lins = lins, filling = filling, n_miss = n_miss)
}

# The matrix door. `obj` is the input the parts were built from, a double
# matrix; `parts` is a list of parts over disjoint columns. Columns no part
# covers are returned as they are, holes included, so a caller that means to
# fill every column passes parts that cover every column. Every cell a map
# names must be a hole (a value over an observed cell is a defect). In each
# part that is not `left_na`, the map must be exactly the holes of the
# part's kept columns, so no hole in them is left out; in a `left_na` part,
# the holes it does not name stay NA.
#
# Without `keep` (or with every column), the result is `obj` with the values
# written, its attributes untouched; the caller sets the class. It is one copy
# of `obj`, made by the first write, and no copy at all when no part writes
# anything. With `keep`, it is one copy of the kept columns, made by `[`, so
# it keeps the input's dimnames for them and no other attribute, and the
# values go to their cells in the copy's own numbering.
#
# `n_miss`, the holes of each column of `obj` before the write, is the
# door's scan, so the input is not read again; NULL counts them here.
# Nothing writes into the input between the scan and this call.
fill_parts <- function(obj, parts, keep = NULL, n_miss = NULL) {
  if (!is.matrix(obj) || !is.double(obj)) {
    part_defect("the input is not a double matrix")
  }
  n_row <- nrow(obj)
  n_col <- ncol(obj)
  keep <- check_keep(keep, n_col)
  chk <- check_parts(obj, parts, keep, n_miss)
  parts <- chk$parts
  lins <- chk$lins

  if (is.null(keep)) {
    out <- obj
    out_col <- seq_len(n_col)
  } else {
    out <- obj[, keep, drop = FALSE]
    # input column -> its column in `out`, 0 for a column not kept
    out_col <- integer(n_col)
    out_col[keep] <- seq_along(keep)
  }
  wrote <- FALSE
  for (i in seq_along(parts)) {
    if (length(lins[[i]])) {
      at <- if (is.null(keep)) {
        lins[[i]]
      } else {
        map <- parts[[i]]$map
        (out_col[map[, "col"]] - 1) * as.double(n_row) + map[, "row"]
      }
      out[at] <- parts[[i]]$values
      wrote <- TRUE
    }
  }

  # Each map lists distinct cells of its part's columns and writes no NA, so
  # the holes left in those columns after the write are the holes before less
  # the map's rows, unless a map cell was not a hole: a value went over an
  # observed cell. For a filling part that leaves none. Counting after the
  # write costs a pass over the result and nothing the size of the maps,
  # which reading the input at every cell before the write would.
  if (wrote) {
    n_miss <- chk$n_miss
    after <- col_scan(out)$n_miss
    for (i in which(lengths(lins) > 0L)) {
      cols <- parts[[i]]$cols
      if (
        sum_holes(after, out_col[cols]) !=
          sum_holes(n_miss, cols) - nrow(parts[[i]]$map)
      ) {
        refuse_observed_cells(obj, parts, lins)
        part_defect(
          "the holes left after the write do not match the maps",
          at = i
        )
      }
    }
  }
  out
}

# The door for what knn_imp_parts() and pca_imp_parts() return: the columns
# in `res$keep` of `res$obj`, filled from `res$parts`. When `res$parts` is
# NULL no kept column had a hole, and they come back as they were.
fill_imp_parts <- function(res) {
  if (!is.null(res$parts)) {
    fill_parts(res$obj, res$parts, keep = res$keep, n_miss = res$n_miss)
  } else if (is.null(res$keep)) {
    res$obj
  } else {
    res$obj[, res$keep, drop = FALSE]
  }
}

# Moves parts built on a slice of an input into the input's numbering, for a
# driver that imputes the input slice by slice (group_imp()). `col_idx` gives
# the input column of each slice column. `keep`, NULL for every slice column,
# is the slice columns to carry: a part's other columns and their map rows
# are dropped, as fill_parts() drops what its `keep` leaves out, and a part
# left with no column is dropped. The slice's columns need not be in the
# input's order, so each part's columns are sorted again, and its map is
# ordered stably by its new columns, which keeps each column's rows in the
# order they had. It builds vectors the size of the maps and nothing of the
# input's size. The parts are not validated here: fill_parts() validates
# them in the input.
rebase_parts <- function(parts, col_idx, keep = NULL) {
  kept <- logical(length(col_idx))
  kept[if (is.null(keep)) seq_along(col_idx) else keep] <- TRUE
  out <- lapply(parts, function(part) {
    cols <- part$cols[kept[part$cols]]
    if (length(cols) == 0L) {
      return(NULL)
    }
    map <- part$map
    values <- part$values
    if (nrow(map) > 0L) {
      rows <- kept[map[, "col"]]
      if (!all(rows)) {
        map <- map[rows, , drop = FALSE]
        values <- values[rows]
      }
      map[, "col"] <- col_idx[map[, "col"]]
      if (is.unsorted(map[, "col"])) {
        o <- order(map[, "col"], method = "radix")
        map <- map[o, , drop = FALSE]
        values <- values[o]
      }
    }
    part$cols <- sort(col_idx[cols])
    part$map <- map
    part$values <- values
    part
  })
  out[!vapply(out, is.null, logical(1))]
}

# Stops at the first map cell of any part that is observed in `obj`. It reads
# `obj` at every cell of the maps, so it costs a vector the length of the
# maps.
refuse_observed_cells <- function(obj, parts, lins) {
  for (i in which(lengths(lins) > 0L)) {
    observed <- which(!is.na(obj[lins[[i]]]))
    if (length(observed)) {
      cell <- parts[[i]]$map[observed[[1]], ]
      part_defect(
        "writes over an observed cell (row ",
        cell[["row"]],
        ", col ",
        cell[["col"]],
        ")",
        at = i
      )
    }
  }
  invisible(NULL)
}
