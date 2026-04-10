#' Print a `slideimp_results` Object
#'
#' @param x An `slideimp_results` object
#' @param n Number of rows to print
#' @param p Number of cols to print
#' @param ... Not used
#'
#' @returns Invisible `x`
#'
#' @export
#'
#' @examples
#' data(khanmiss1)
#' t_khanmiss1 <- t(khanmiss1)
#' result <- knn_imp(t_khanmiss1, k = 5)
#' print(result, n = 6, p = 6)
print.slideimp_results <- function(x, n = 5, p = 5, ...) {
  imp_method <- attr(x, "imp_method")

  # Cleaner header
  cat("slideimp_results (", toupper(imp_method), ")\n", sep = "")
  cat("Dimensions: ", nrow(x), " x ", ncol(x), "\n\n", sep = "")

  # Print subset
  subset_x <- x[seq_len(min(n, nrow(x))), seq_len(min(p, ncol(x))), drop = FALSE]
  print(subset_x, ...)

  # Better truncation message
  if (n < nrow(x) || p < ncol(x)) {
    cat("\n# Showing [1:", min(n, nrow(x)), ", 1:", min(p, ncol(x)),
      "] of full matrix\n",
      sep = ""
    )
  }

  invisible(x)
}


#' Print a `slideimp_sim` Object
#'
#' @param x A `slideimp_sim` object.
#' @param n Number of rows of each component to show.
#' @param p Number of columns of `input` to show.
#' @param ... Not used.
#'
#' @returns Invisible `x`
#'
#' @export
print.slideimp_sim <- function(x, n = 6L, p = 6L, ...) {
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
  print(d[seq_len(min(n, nr)), seq_len(min(p, nc)), drop = FALSE])

  if (n < nr || p < nc) {
    cat("# Showing [1:", min(n, nr), ", 1:", min(p, nc),
      "] of full matrix\n",
      sep = ""
    )
  }

  invisible(x)
}


#' Print a `slideimp_tbl` Object
#'
#' @param x A `slideimp_tbl` object.
#' @param n Number of rows to show. Defaults to 10.
#' @param ... Not used.
#'
#' @returns Invisible `x`
#'
#' @export
print.slideimp_tbl <- function(x, n = NULL, ...) {
  if (is.null(n)) n <- 10L
  n <- min(n, nrow(x))
  cat(sprintf("# slideimp table: %d x %d\n", nrow(x), ncol(x)))
  if (nrow(x) == 0L) {
    return(invisible(x))
  }
  disp <- as.data.frame(x)
  class(disp) <- "data.frame" # prevent recursion
  for (nm in names(disp)) {
    if (is.list(disp[[nm]])) {
      disp[[nm]] <- vapply(disp[[nm]], function(elem) {
        if (is.null(elem)) {
          "<NULL>"
        } else if (is.data.frame(elem)) {
          sprintf("<df [%d x %d]>", nrow(elem), ncol(elem))
        } else if (is.list(elem)) {
          sprintf("<list [%d]>", length(elem))
        } else {
          sprintf("<%s [%d]>", typeof(elem), length(elem))
        }
      }, character(1L))
    }
  }
  print(disp[seq_len(n), , drop = FALSE], row.names = FALSE, ...)
  if (nrow(x) > n) cat(sprintf("# ... with %d more rows\n", nrow(x) - n))
  invisible(x)
}
