#' Print a `slideimp_results` Object
#'
#' @param x An `slideimp_results` object
#' @param n Number of rows to print
#' @param p Number of cols to print
#' @param ... Not used
#'
#' @returns Invisible `x`
#'
#' @examples
#' set.seed(1234)
#' mat <- sim_mat(n = 10, p = 10)
#' result <- knn_imp(mat$input, k = 5)
#' class(result)
#' print(result, n = 6, p = 6)
#'
#' @export
print.slideimp_results <- function(x, n = 6L, p = 6L, ...) {
  imp_method <- toupper(attr(x, "imp_method"))
  metacaller <- attr(x, "metacaller")
  fallback <- attr(x, "fallback")
  if (!is.null(metacaller)) {
    cat("Method: ", metacaller, " (", imp_method, " imputation",")\n", sep = "")
  } else {
    cat("Method: ", imp_method, " imputation", "\n", sep = "")
  }
  cat("Dimensions: ", nrow(x), " x ", ncol(x), "\n", sep = "")
  # fallback notes
  if (!is.null(metacaller) && length(fallback) > 0L) {
    cat(
      "Note: ", length(fallback), " group(s) fell back to mean imputation: ",
      fmt_trunc(fallback), "\n",
      sep = ""
    )
  } else if (isTRUE(fallback)) {
    if (isTRUE(attr(x, "post_imp"))) {
      cat("Note: fell back to mean imputation\n")
    } else {
      cat("Note: imputation skipped (insufficient eligible columns)\n")
    }
  }
  # remaining NA note
  if (isTRUE(attr(x, "has_remaining_na"))) {
    cat("Note: matrix still contains NA values\n")
  }

  cat("\n")
  subset_x <- x[seq_len(min(n, nrow(x))), seq_len(min(p, ncol(x))), drop = FALSE]
  print(subset_x, ...)
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
#' @examples
#' set.seed(123)
#' sim_data <- sim_mat(n = 50, p = 10, rho = 0.5)
#' class(sim_data)
#' print(sim_data)
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
#' @examples
#' mat <- sim_mat(n = 10, p = 500)
#' set.seed(1234)
#' results <- tune_imp(mat$input, parameters = data.frame(k = 5), .f = "knn_imp")
#' class(results)
#' print(results)
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
