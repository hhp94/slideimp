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

#' Ampute NA Given the Output of `sample_each_rep()`
#'
#' Used in testthat only.
#'
#' @param obj Input.
#' @param loc Output of `sample_each_rep()`.
#'
#' @returns `NULL` invisibly. Called for its side effects (messages or errors).
#'
#' @noRd
#' @keywords internal
apply_na <- function(obj, loc) {
  obj[cbind(loc[, "row"], loc[, "col"])] <- NA_real_
  obj
}

#' Resolve a Subset Argument to Sorted Integer Column Indices
#'
#' Resolve the subset argument from integers, characters, to integers.
#'
#' @param subset NULL (all columns), character, or integerish vector.
#' @param obj matrix/data.frame whose columns are referenced.
#' @param sort whether to sort the resulting indices (needed by `slide_imp()`).
#'
#' @returns integer vector of column indices, or NULL to signal early return.
#'
#' @noRd
#' @keywords internal
resolve_subset <- function(subset, obj, sort = FALSE) {
  nc <- ncol(obj)
  cn <- colnames(obj)

  if (is.null(subset)) {
    subset <- seq_len(nc)
  } else if (is.character(subset)) {
    if (is.null(cn)) {
      cli::cli_abort("{.arg subset} contains characters but {.arg obj} has no column names.")
    }
    if (anyDuplicated(subset)) {
      cli::cli_abort("{.arg subset} contains duplicate feature names.")
    }
    matched <- match(subset, cn, nomatch = NA_integer_)
    if (anyNA(matched)) {
      cli::cli_inform("Feature(s) in {.arg subset} not found in {.code colnames(obj)} and dropped.")
    }
    subset <- matched[!is.na(matched)]
  } else {
    checkmate::assert_integerish(
      subset,
      lower = 1L, upper = nc, any.missing = FALSE, unique = TRUE,
      min.len = 0L, .var.name = "`subset`"
    )
    subset <- as.integer(subset)
  }

  if (length(subset) == 0L) {
    cli::cli_inform("No features in {.arg subset} detected. No imputation was performed.")
    return(NULL)
  }

  if (sort) sort(subset) else subset
}

#' @noRd
#' @keywords internal
as_slideimp_results <- function(
  obj,
  imp_method,
  fallback,
  post_imp,
  na_check,
  has_remaining_na = if (na_check) anyNA(obj) else NULL
) {
  class(obj) <- c("slideimp_results", class(obj))
  attr(obj, "imp_method") <- imp_method
  attr(obj, "fallback") <- fallback
  attr(obj, "post_imp") <- post_imp
  attr(obj, "has_remaining_na") <- has_remaining_na
  obj
}
