#' Print a `slideimp_results` Object
#'
#' Print the output of [knn_imp()], [pca_imp()], [group_imp()], [slide_imp()].
#'
#' @param x A `slideimp_results` object.
#' @param n Number of rows to print.
#' @param p Number of cols to print.
#' @param ... Not used.
#'
#' @returns Invisible `x`.
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
  fallback_action <- attr(x, "fallback_action")

  # Header
  if (!is.null(metacaller)) {
    cat("Method: ", metacaller, " (", imp_method, " imputation)\n", sep = "")
  } else {
    cat("Method: ", imp_method, " imputation\n", sep = "")
  }
  cat("Dimensions: ", nrow(x), " x ", ncol(x), "\n", sep = "")

  # Fallback note
  if (!is.null(metacaller) && length(fallback) > 0L) {
    unit <- if (identical(metacaller, "slide_imp")) "window" else "group"
    action <- if (is.null(fallback_action)) {
      "used a fallback"
    } else {
      switch(fallback_action,
        mean = "used mean imputation as fallback",
        skip = "skipped imputation (insufficient eligible columns; original values retained)",
        "used a fallback"
      )
    }
    n_fb <- length(fallback)
    unit_plural <- if (n_fb == 1L) unit else paste0(unit, "s")
    cat(
      "Note: ", n_fb, " ", unit_plural, " ", action, ".\n",
      "  See ", unit_plural, ": ", fmt_trunc(fallback), "\n",
      sep = ""
    )
  }

  # Remaining NA note
  if (isTRUE(attr(x, "has_remaining_na"))) {
    cat("Note: requested columns still contain NA values\n")
  }
  cat("\n")

  # Preview
  n_show <- min(n, nrow(x))
  p_show <- min(p, ncol(x))
  print(x[seq_len(n_show), seq_len(p_show), drop = FALSE], ...)
  if (n_show < nrow(x) || p_show < ncol(x)) {
    cat(sprintf(
      "# Showing %d of %d rows and %d of %d columns\n",
      n_show, nrow(x), p_show, ncol(x)
    ))
  }
  invisible(x)
}

#' Print a `slideimp_sim` Object
#'
#' Print the output of [sim_mat()].
#'
#' @param x A `slideimp_sim` object.
#' @param n Number of rows of each component to show.
#' @param p Number of columns of `input` to show.
#' @param ... Not used.
#'
#' @returns Invisible `x`.
#'
#' @examples
#' set.seed(123)
#' sim_data <- sim_mat(n = 50, p = 10, rho = 0.5)
#' class(sim_data)
#' print(sim_data)
#'
#' @export
print.slideimp_sim <- function(x, n = 6L, p = 6L, ...) {
  # $col_group
  cg_n <- nrow(x$col_group)
  cg_show <- min(n, cg_n)
  cat("$col_group (", x$n_col_groups, " column groups)\n", sep = "")
  print(x$col_group[seq_len(cg_show), , drop = FALSE])
  if (cg_show < cg_n) {
    cat(sprintf("# Showing %d of %d rows\n", cg_show, cg_n))
  }
  cat("\n")

  # $row_group
  rg_n <- nrow(x$row_group)
  rg_show <- min(n, rg_n)
  cat("$row_group (", x$n_row_groups, " row groups)\n", sep = "")
  print(x$row_group[seq_len(rg_show), , drop = FALSE])
  if (rg_show < rg_n) {
    cat(sprintf("# Showing %d of %d rows\n", rg_show, rg_n))
  }
  cat("\n")

  # $input
  d <- x$input
  nr <- nrow(d)
  nc <- ncol(d)
  n_show <- min(n, nr)
  p_show <- min(p, nc)
  cat("$input (", nr, " x ", nc, ")\n", sep = "")
  print(d[seq_len(n_show), seq_len(p_show), drop = FALSE])
  if (n_show < nr || p_show < nc) {
    cat(sprintf(
      "# Showing %d of %d rows and %d of %d columns\n",
      n_show, nr, p_show, nc
    ))
  }

  invisible(x)
}

#' Print a `slideimp_tbl` Object
#'
#' Print `slideimp_tbl` objects (which inherit `data.frame`) with nicer looking
#' list-columns (similar to `tibble`).
#'
#' @param x A `slideimp_tbl` object.
#' @param n Number of rows to show. Defaults to 10.
#' @param ... Not used.
#'
#' @returns Invisible `x`.
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
  n_show <- min(n, nrow(x))
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
  print(disp[seq_len(n_show), , drop = FALSE], row.names = FALSE, ...)
  if (nrow(x) > n_show) {
    cat(sprintf("# ... with %d more rows\n", nrow(x) - n_show))
  }
  invisible(x)
}
