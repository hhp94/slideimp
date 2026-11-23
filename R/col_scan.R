# The column scan: what every entry point knows about its input before any real
# work, from one pass over the columns it reads. It counts and never refuses;
# each entry point applies its own policy to the result, through the helpers
# below. Every entry point takes a double matrix only, and refuses any other
# with checkmate before it gets here, so nothing converts one.

# One pass over the columns `cols` of the double matrix `obj`, every column
# when `cols` is NULL. Returns a list, each per-column field in the order of
# `cols`:
# - `n_row`: the rows of `obj`;
# - `n_miss`: integer, the missing cells (NA and NaN) of each column read;
# - `has_inf`: logical, whether each column read holds Inf or -Inf;
# - `var`: with `var = TRUE`, the sample variance of each column read, as
#   col_vars() gives it, otherwise NULL;
# - `row_miss`: with `rows = TRUE`, the missing cells of each row over the
#   columns read, otherwise NULL;
# - `n_below`, `n_above`: with `tails = c(below, above)`, integer, the cells
#   of each column read below `below` and above `above`, NA and NaN in
#   neither, otherwise NULL.
col_scan <- function(
  obj,
  cols = NULL,
  var = FALSE,
  rows = FALSE,
  tails = NULL
) {
  if (!is.matrix(obj) || !is.double(obj)) {
    stop("internal error: col_scan() takes a double matrix", call. = FALSE)
  }
  if (!is.null(tails) && (!is.double(tails) || length(tails) != 2L)) {
    stop(
      "internal error: col_scan() takes tails as two doubles",
      call. = FALSE
    )
  }
  j <- if (is.null(cols)) seq_len(ncol(obj)) else cols
  scan <- col_scan_cpp(
    obj,
    j - 1L,
    var,
    rows,
    if (is.null(tails)) double(0) else tails
  )
  c(list(n_row = nrow(obj)), scan)
}

# The scan of the columns `j` of the matrix that `scan` was taken of,
# for a block of those columns (a group of group_imp(), a window of
# slide_imp()). A column's counts and variance do not depend on the other
# columns, so they are the scan's own, indexed. The row counts of a block
# do, so they are not carried: a block that needs them counts its own.
scan_cols <- function(scan, j) {
  list(
    n_row = scan$n_row,
    n_miss = scan$n_miss[j],
    has_inf = scan$has_inf[j],
    var = if (is.null(scan$var)) NULL else scan$var[j],
    row_miss = NULL,
    n_below = if (is.null(scan$n_below)) NULL else scan$n_below[j],
    n_above = if (is.null(scan$n_above)) NULL else scan$n_above[j]
  )
}

# The refusal of the imputers: an infinite value in any column the scan
# read, or a column with no observed value among them, or among its columns
# `empty_j` when given (slide_imp() refuses one only where a window reads
# it). One message for both, naming no cell or column; the help section
# "Missing values and non-finite input" says how to find them.
refuse_nonfinite <- function(scan, empty_j = NULL) {
  empty <- scan$n_miss == scan$n_row
  if (!is.null(empty_j)) {
    empty <- empty[empty_j]
  }
  if (any(scan$has_inf) || any(empty)) {
    stop(
      "'obj' holds infinite values or columns with no observed value. ",
      "Remove them before imputation",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# The refusal of an entry point that accepts a column with no observed value:
# an infinite value in the columns the scan read, or in its columns `j`.
refuse_inf <- function(scan, j = NULL) {
  inf <- if (is.null(j)) scan$has_inf else scan$has_inf[j]
  if (any(inf)) {
    stop("'obj' holds infinite values, which are not supported", call. = FALSE)
  }
  invisible(NULL)
}
