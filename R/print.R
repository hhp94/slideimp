#' Print method for slide_impList
#' @param x A slide_impList object
#' @param m Number of rows to show in preview (default: 5)
#' @param n Number of columns to show in preview (default: 10)
#' @param ... Additional arguments (not used)
#' @export
print.slide_impList <- function(x, m = 5, n = 10, ...) {
  # Determine object class for display
  obj_class <- ifelse(inherits(x, "slide_impList"), "slide_impList", "KnnImpList")
  n_imp <- length(x)
  rn <- attr(x, "rownames")
  cn <- attr(x, "colnames")
  n_rows <- nrow(x[[1]])
  n_cols <- ncol(x[[1]])
  # Matrix type - slide_impList is always big.matrix
  mat_type <- if (obj_class == "slide_impList") {
    "big.matrix"
  } else if (bigmemory::is.big.matrix(x[[1]])) {
    "big.matrix"
  } else {
    "matrix"
  }
  # Print header with appropriate class name
  cat(
    sprintf(
      "%s: List of %d imputation%s of a %d x %d %s\n",
      obj_class,
      n_imp,
      ifelse(n_imp == 1, "", "s"),
      n_rows,
      n_cols,
      mat_type
    )
  )
  # Subset info
  subset_attr <- attr(x, "subset")
  total_cols <- attr(x, "ncol")
  if (!is.null(subset_attr) && !is.null(total_cols) && length(subset_attr) < total_cols) {
    sliding_text <- ifelse(obj_class == "slide_impList", " (sliding window)", "")
    cat(
      sprintf(
        "Imputed columns: %d of %d total columns%s\n",
        length(subset_attr),
        total_cols,
        sliding_text
      )
    )
  }
  cat("\n")
  # Show preview of first imputation
  cat("Preview of imputation 1:\n")
  if (n_imp > 0 && n_rows > 0 && n_cols > 0) {
    preview_rows <- min(m, n_rows)
    preview_cols <- min(n, n_cols)
    preview_data <- x[[1]][1:preview_rows, 1:preview_cols, drop = FALSE]
    print(preview_data, digits = 4)
    # Remaining rows/cols info
    if (n_rows > m || n_cols > n) {
      cat("\n")
      if (n_rows > m && n_cols > n) {
        cat(sprintf(
          "[ ... %d more rows, %d more columns ]\n",
          n_rows - m, n_cols - n
        ))
      } else if (n_rows > m) {
        cat(sprintf("[ ... %d more rows ]\n", n_rows - m))
      } else if (n_cols > n) {
        cat(sprintf("[ ... %d more columns ]\n", n_cols - n))
      }
    }
    # Additional imputations info
    if (n_imp > 1) {
      cat(sprintf(
        "\n[ ... and %d more imputation%s ]\n",
        n_imp - 1,
        ifelse(n_imp - 1 == 1, "", "s")
      ))
    }
  }
  invisible(x)
}

#' Print method for KnnImpList
#' @param x A KnnImpList object
#' @param m Number of rows to show in preview (default: 5)
#' @param n Number of columns to show in preview (default: 10)
#' @param ... Additional arguments (not used)
#' @export
print.KnnImpList <- function(x, m = 5, n = 10, ...) {
  # Simply call the KnnImpList print method which handles both classes
  print.slide_impList(x, m, n, ...)
}
