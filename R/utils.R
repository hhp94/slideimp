#' Format Truncated Output for Messages
#'
#' Helper function to format a vector as a comma-separated string, truncating
#' it to a maximum number of elements to keep warning and error messages readable.
#' Appends a hint pointing the user to the `prep_groups` documentation.
#'
#' @param x A character to format.
#' @param n Integer scalar. The maximum number of elements to display. Default is `5`.
#'
#' @returns A length-1 character string.
#'
#' @noRd
#' @keywords internal
fmt_trunc <- function(x, n = 5) {
  n <- min(n, length(x))
  truncated <- x[seq_len(n)]
  suffix <- if (length(x) > n) ", ..." else ""
  paste0(paste(truncated, collapse = ", "), suffix)
}

#' Check and Validate BLAS Thread Pinning
#'
#' Verifies that the `RhpcBLASctl` package is installed if BLAS pinning is
#' explicitly requested (`pin_blas = TRUE`). If pinning is not requested but
#' multiple BLAS threads are detected, it emits a helpful tip suggesting
#' pinning for better parallel performance.
#'
#' @param pin_blas Logical scalar. Whether BLAS thread pinning is requested.
#'
#' @returns `NULL` invisibly. Called for its side effects (messages or errors).
#'
#' @noRd
#' @keywords internal
check_pin_blas <- function(pin_blas) {
  has_pkg <- requireNamespace("RhpcBLASctl", quietly = TRUE)
  if (pin_blas && !has_pkg) {
    cli::cli_abort(c(
      "{.code pin_blas = TRUE} requires the {.pkg RhpcBLASctl} package.",
      "i" = "Install it with {.code install.packages(\"RhpcBLASctl\")}"
    ))
  }
  if (!pin_blas && has_pkg && RhpcBLASctl::blas_get_num_procs() > 1L) {
    message("Tip: set `pin_blas = TRUE` may improve parallel performance.")
  }
}

#' Ampute NA given the output of `sample_each_rep()`
#'
#' @param obj Input
#' @param loc Output of `sample_each_rep()`
#'
#' @returns `NULL` invisibly. Called for its side effects (messages or errors).
#'
#' @noRd
#' @keywords internal
apply_na <- function(obj, loc) {
  obj[cbind(loc[, "row"], loc[, "col"])] <- NA_real_
  obj
}
