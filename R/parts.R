# Parts: what the imputers pass to each other, and the matrix door that turns
# them into a result.
#
# A part records what one imputer or fallback did to some columns of an input:
# the columns it covers, the holes it filled and the value of each, all in the
# INPUT's numbering. It holds no copy of the input, so passing parts between
# functions copies nothing of matrix size. fill_parts() copies the input once
# and writes every part's values at its holes. parts_to_si_mi() is the second
# door: it builds a one-draw `si_mi`, the object mipca() returns, from the
# same parts.
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

# The checks both doors make on their parts, given the input: every part
# valid, no column in two parts, and each part that is not `left_na` lists as
# many cells as its columns have holes. Returns each part's linear index (from
# validate_part()) and which parts fill.
check_parts <- function(obj, parts) {
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
  list(lins = lins, filling = filling)
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
  chk <- check_parts(obj, parts)
  lins <- chk$lins
  filling <- chk$filling

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
    refuse_observed_cells(obj, parts, lins, filling)
    part_defect("leaves a hole open in a column it fills")
  }
  out
}

# Stops at the first map cell of a filling part that is observed in `obj`.
# It reads `obj` at every cell of the maps, so it costs a vector the length
# of the maps.
refuse_observed_cells <- function(obj, parts, lins, filling) {
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
  invisible(NULL)
}

# The user-facing rule of the `si_mi` door: `output = "mipca"` takes only an
# input with unique row and column names, as mipca() does (R/mi_pca_imp.R:
# `.chk_matrix <- function(obj) {`). The join matches columns and rows by
# name, so names an imputer made up would join two unnamed fits by position,
# unchecked. A name that is NA or empty counts as missing, as checkmate's
# "unique" reads it there. Called before the fit, so no fit is spent on an
# input the door refuses. `call` is the imputer's frame, so the error names
# the function the user called.
check_si_mi_names <- function(obj, call = parent.frame()) {
  # Each problem is rendered here and reaches cli_abort() below as an
  # interpolated value, which cli does not parse again: a name may hold `{`.
  problem <- function(nms, what) {
    if (is.null(nms)) {
      return(cli::format_inline("{.arg obj} has no {what} names."))
    }
    bad <- sum(is.na(nms) | !nzchar(nms))
    if (bad > 0L) {
      return(cli::format_inline(
        "{.arg obj} has {bad} missing or empty {what} {cli::qty(bad)}name{?s}."
      ))
    }
    dup <- unique(nms[duplicated(nms)])
    if (length(dup)) {
      return(cli::format_inline(
        "{.arg obj} repeats {length(dup)} {what} {cli::qty(length(dup))}name{?s}: {fmt_trunc(dup, 6)}."
      ))
    }
    NULL
  }
  found <- c(problem(rownames(obj), "row"), problem(colnames(obj), "column"))
  if (length(found)) {
    bullets <- stats::setNames(
      sprintf("{found[[%d]]}", seq_along(found)),
      rep("x", length(found))
    )
    cli::cli_abort(
      c(
        "With {.code output = 'mipca'}, {.arg obj} must have unique row and column names.",
        bullets,
        "i" = "Set {.code rownames(obj)} and {.code colnames(obj)} to unique values, or use {.code output = 'matrix'}."
      ),
      call = call
    )
  }
  invisible(NULL)
}

# The `si_mi` door. From the same parts as fill_parts(), and with the same
# checks, it builds a one-draw `si_mi`: the object mipca() returns, read with
# as.matrix(), mipca_complete(), summary() and cbind(). This file is the one
# caller of the vendored internals it uses (dev/vendored.md), so a port that
# changes them breaks here only.
#
# - `obj`: the input the parts were built from, a numeric matrix with unique
#   row and column names (check_si_mi_names() has refused any other).
# - `parts`: parts over disjoint columns that together cover EVERY column:
#   the placement that puts the columns back in input order needs each column
#   once.
# - `call`: the caller's match.call(), and `fn` the caller's name.
#   .scrub_call() rewrites it so that a matrix passed by value through
#   do.call() is not stored in the fit.
# - `prov_extra`: NULL, or one list per part, of `prov` fields the imputer
#   records on that part's chunk (the rank and iteration count of a PCA fit).
#
# Each part with columns becomes one chunk; a part over no columns is
# dropped, since a chunk needs one. A chunk holds the part's columns of `obj`
# (a slice, with its row names dropped since the object holds them once; no
# slice when one part covers every column), the map in the slice's own
# numbering, the values as the one draw, the identity backward map, and the
# empty diagnostic frame an ABB piece takes: one draw has nothing to
# diagnose. The part's `info`, when it has any, is stored on the chunk, since
# `.si_mi_new()` rebuilds the object from its own fields on every join.
parts_to_si_mi <- function(obj, parts, call, fn, prov_extra = NULL) {
  if (!is.matrix(obj) || !is.numeric(obj)) {
    part_defect("the input is not a numeric matrix")
  }
  rn <- rownames(obj)
  cn <- colnames(obj)
  if (is.null(rn) || is.null(cn) || anyDuplicated(rn) || anyDuplicated(cn)) {
    part_defect("the si_mi door needs unique row and column names")
  }
  chk <- check_parts(obj, parts)
  n_col <- ncol(obj)
  cols <- lapply(parts, `[[`, "cols")
  # check_parts() refused a column in two parts, so this many columns means
  # every column
  n_cov <- sum(lengths(cols))
  if (n_cov != n_col) {
    part_defect(
      "the parts cover ",
      n_cov,
      " of the ",
      n_col,
      " columns, and the si_mi door needs every column"
    )
  }
  # Here, at birth, rather than at the first as.matrix(), which asserts it
  # too, so a defect shows where it was made.
  refuse_observed_cells(obj, parts, chk$lins, chk$filling)

  if (is.null(prov_extra)) {
    prov_extra <- vector("list", length(parts))
  }
  door_prov <- c(
    "chains",
    "thin",
    "ndraws",
    "nmiss",
    "thin_final",
    "thin_effective",
    "seam_state",
    "left_na"
  )
  ok_extra <- function(x) {
    is.null(x) ||
      (is.list(x) &&
        length(x) > 0L &&
        !is.null(names(x)) &&
        all(nzchar(names(x))) &&
        !any(names(x) %in% door_prov))
  }
  if (
    !is.list(prov_extra) ||
      length(prov_extra) != length(parts) ||
      !all(vapply(prov_extra, ok_extra, logical(1)))
  ) {
    part_defect(
      "`prov_extra` is not one named list (or NULL) per part, of fields the door does not set"
    )
  }

  call <- .scrub_call(call, fn)
  keep <- which(lengths(cols) > 0L)
  whole <- length(keep) == 1L
  fits <- lapply(keep, function(i) {
    part_chunk_fit(obj, parts[[i]], whole, prov_extra[[i]], call)
  })
  if (length(fits) == 1L) {
    return(fits[[1L]])
  }
  .si_mi_join(
    fits,
    call = call,
    placement = unlist(cols[keep], use.names = FALSE)
  )
}

# One part as a one-chunk `si_mi`; parts_to_si_mi() checked the part.
part_chunk_fit <- function(obj, part, whole, prov_extra, call) {
  cols <- part$cols
  map <- part$map
  if (whole) {
    sub <- obj
  } else {
    sub <- obj[, cols, drop = FALSE]
    dimnames(sub) <- list(NULL, colnames(sub))
    # the map's columns in the slice's numbering, an integer matrix still
    inv <- integer(ncol(obj))
    inv[cols] <- seq_along(cols)
    map[, "col"] <- inv[map[, "col"]]
  }
  prov <- c(
    list(
      chains = 1L,
      thin = 1L,
      ndraws = 1L,
      nmiss = nrow(map),
      thin_final = 1L,
      thin_effective = 1L,
      seam_state = "never"
    ),
    prov_extra
  )
  if (part$left_na) {
    prov$left_na <- TRUE
  }
  chunk <- .si_mi_chunk(
    obj = sub,
    map = map,
    draws = matrix(part$values, ncol = 1L),
    transform = tf_identity(),
    diag = .diag_frame(
      map = map[integer(0), , drop = FALSE],
      draws = list(matrix(NA_real_, 0L, 1L)),
      p = ncol(sub),
      nms = colnames(sub),
      prop_diag = .prop_diag_default()
    ),
    prov = prov,
    call = call,
    info = if (length(part$info)) part$info
  )
  .si_mi_new(
    chunks = list(chunk),
    rownames = rownames(obj),
    nkeep = 1L,
    call = call
  )
}
