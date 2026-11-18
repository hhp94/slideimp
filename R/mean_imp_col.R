#' Column Mean Imputation
#'
#' Set the holes of the selected columns of a matrix to the mean of each
#' column's observed values. The fallback the imputers share: it builds one
#' mean part with `fallback_part()` and fills it through the matrix door,
#' `fill_parts()`.
#'
#' @param obj A numeric matrix.
#' @param subset Optional character or integer vector of the columns to fill.
#'   If `NULL`, every column is filled.
#'
#' @returns A double matrix with the dimensions and dimnames of `obj` and no
#'   other attribute, with the holes of the selected columns set to their
#'   column means. A column with no observed value is returned as it is. When
#'   `subset` resolves to no column, `obj` itself is returned.
#'
#' @noRd
#' @keywords internal
mean_imp_col <- function(obj, subset = NULL) {
  # pre-cond
  checkmate::assert_matrix(obj, mode = "numeric", .var.name = "obj")
  # before the subset early return, which would skip this check
  check_inf(obj)

  # subset resolution
  subset <- resolve_subset(subset, obj)
  if (is.null(subset)) {
    return(obj)
  }
  if (!is.double(obj)) {
    storage.mode(obj) <- "double"
  }

  # imputation. A column with no hole has nothing to fill, and one with no
  # observed value has no mean. Both are left out of the part and come back
  # from the door as they are, so fallback_part() copies only the columns
  # it fills.
  n_miss <- col_miss_internal(obj)
  cols <- sort(subset)
  cols <- cols[n_miss[cols] > 0L & n_miss[cols] < nrow(obj)]
  out <- fill_parts(obj, list(fallback_part(obj, cols, "mean")))

  # The result carries dim and dimnames only. When the door wrote, `out` is
  # its own copy and this edits it in place.
  keep <- intersect(c("dim", "dimnames"), names(attributes(out)))
  if (length(attributes(out)) > length(keep)) {
    attributes(out) <- attributes(out)[keep]
  }
  out
}

# A fallback part over `cols` of `obj`: each hole in those columns set to the
# mean or the median of its column's observed values, and `kind` set to
# `method`. `cols` is integer and strictly increasing, and none of its columns
# is without an observed value: such a column has neither statistic, and a
# part over one fails validation, its values holding NA.
#
# Vectorized over the columns: their submatrix is copied out once, with its
# is.na(). So its callers, mean_imp_col() and pca_imp(), pass only columns
# with a hole: a complete column would be copied to fill nothing. The
# fallback covers few columns beside those the imputer filled, so a copy of
# them is small. The mean is .colMeans(), a long double sum
# divided in long double, within one ulp of the exact mean on the chr21 data
# and on a 2000 x 20000 normal matrix (dev/probe_mean_source.R). The median
# is collapse::fmedian(), the mean of the two middle values for an even
# count, as stats::median() takes it.
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
