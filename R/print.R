#' Print SlideImpImputedMatrix
#'
#' @param x An SlideImpImputedMatrix
#' @param n Number of rows to print
#' @param m Number of cols to print
#' @param ... Not used
#'
#' @returns Invisible object of class SlideImpImputedMatrix
#' @export
#'
#' @examples
#' data(khanmiss1)
#' t_khanmiss1 <- t(khanmiss1)
#' result <- knn_imp(t_khanmiss1, k = 5)
#' print(result, n = 6, m = 6)
print.SlideImpImputedMatrix <- function(x, n = 5, m = 5, ...) {
  imp_method <- attr(x, "imp_method")

  # Cleaner header
  cat("SlideImpImputedMatrix (", toupper(imp_method), ")\n", sep = "")
  cat("Dimensions: ", nrow(x), " x ", ncol(x), "\n\n", sep = "")

  # Print subset
  subset_x <- x[seq_len(min(n, nrow(x))), seq_len(min(m, ncol(x))), drop = FALSE]
  print(subset_x, ...)

  # Better truncation message
  if (n < nrow(x) || m < ncol(x)) {
    cat("\n# Showing [1:", min(n, nrow(x)), ", 1:", min(m, ncol(x)),
      "] of full matrix\n",
      sep = ""
    )
  }

  invisible(x)
}


#' Print a SlideImpSimMat
#'
#' @param x A `SlideImpSimMat` object.
#' @param n Number of rows of each component to show.
#' @param m Number of columns of `input` to show.
#' @param ... Not used.
#'
#' @returns Invisibly, `x`.
#' @export
print.SlideImpSimMat <- function(x, n = 6L, m = 6L, ...) {
  cat("$col_group (", x$n_col_groups, " column groups)\n", sep = "")
  print(x$col_group[seq_len(min(n, nrow(x$col_group))), , drop = FALSE])
  cat("\n")

  cat("$row_group (", x$n_row_groups, " row groups)\n", sep = "")
  print(x$row_group[seq_len(min(n, nrow(x$row_group))), , drop = FALSE])
  cat("\n")

  d <- x$input
  nr <- nrow(d)
  nc <- ncol(d)
  cat("$input (", nr, " x ", nc, ")\n", sep = "")
  print(d[seq_len(min(n, nr)), seq_len(min(m, nc)), drop = FALSE])

  if (n < nr || m < nc) {
    cat("# Showing [1:", min(n, nr), ", 1:", min(m, nc),
      "] of full matrix\n",
      sep = ""
    )
  }

  invisible(x)
}
