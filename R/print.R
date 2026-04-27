#' Print a `slideimp_results` Object
#'
#' Print the output of [knn_imp()], [pca_imp()], [group_imp()], or
#' [slide_imp()].
#'
#' @param x A `slideimp_results` object.
#' @param n Number of rows to print.
#' @param p Number of columns to print.
#' @param ... Not used.
#'
#' @returns `x`, invisibly.
#'
#' @examples
#' set.seed(1234)
#' mat <- sim_mat(n = 10, p = 10)
#' result <- knn_imp(mat$input, k = 5, .progress = FALSE)
#' class(result)
#' print(result, n = 6, p = 6)
#'
#' @method print slideimp_results
#' @export
print.slideimp_results <- function(x, n = 6L, p = 6L, ...) {
  imp_method <- toupper(attr(x, "imp_method"))
  metacaller <- attr(x, "metacaller")
  fallback <- attr(x, "fallback")
  fallback_action <- attr(x, "fallback_action")

  # header
  if (!is.null(metacaller)) {
    cat("Method: ", metacaller, " (", imp_method, " imputation)\n", sep = "")
  } else {
    cat("Method: ", imp_method, " imputation\n", sep = "")
  }
  cat("Dimensions: ", nrow(x), " x ", ncol(x), "\n", sep = "")

  # fallback note
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

  # remaining NA note
  if (isTRUE(attr(x, "has_remaining_na"))) {
    cat("Note: requested columns still contain NA values\n")
  }
  cat("\n")

  # preview
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
#' @returns `x`, invisibly.
#'
#' @examples
#' set.seed(123)
#' sim_data <- sim_mat(n = 50, p = 10, rho = 0.5)
#' class(sim_data)
#' print(sim_data)
#'
#' @method print slideimp_sim
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
#' Print `slideimp_tbl` objects, which inherit from `data.frame`, with compact
#' display of list-columns.
#'
#' @param x A `slideimp_tbl` object.
#' @param n Number of rows to show. If `NULL`, a default is used.
#' @param ... Not used.
#'
#' @returns `x`, invisibly.
#'
#' @examples
#' sim <- sim_mat(n = 10, p = 20)
#' tbl <- prep_groups(colnames(sim$input), sim$col_group)
#' class(tbl)
#' print(tbl)
#'
#' @method print slideimp_tbl
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
