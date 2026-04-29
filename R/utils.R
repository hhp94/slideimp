# S3 class ---------------------------------------------------------------------

#' @noRd
#' @keywords internal
new_slideimp_results <- function(
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

# Formatting -------------------------------------------------------------------

#' Format Truncated Output for Messages
#'
#' Helper function to format a vector as a comma-separated string, truncating
#' it to a maximum number of elements to keep warning and error messages readable.
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

# Input validation -------------------------------------------------------------

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

# Param ----
#' Method-specific parameter registry
#'
#' @description
#' Internal registry defining required, allowed, and valid method parameters for
#' each imputation backend.
#'
#' @format
#' A named list with one entry per imputation mode. Each mode contains:
#'
#' \describe{
#'   \item{`required`}{The required tuning parameter for the method, such as
#'   `"k"` for K-NN or `"ncp"` for PCA.}
#'   \item{`methods`}{Valid values for the method-specific `method` argument.}
#'   \item{`allowed`}{Parameter names accepted for that imputation backend.}
#' }
#'
#' @details
#' This object is used internally by `check_unknown_params()`, `group_imp()`,
#' and `tune_imp()` to validate user-supplied parameter names and method values.
#'
#' It is intentionally not exported.
#'
#' @keywords internal
#' @noRd
.param_registry <- list(
  knn = list(
    required = "k",
    methods = c("euclidean", "manhattan"),
    allowed = c(
      "k", "method", "colmax", "post_imp",
      "dist_pow", "tree", ".progress"
    )
  ),
  pca = list(
    required = "ncp",
    methods = c("regularized", "EM"),
    allowed = c(
      "ncp", "scale", "method", "coeff.ridge", "row.w",
      "threshold", "seed", "nb.init", "maxiter", "miniter",
      "lobpcg_control", "solver", "colmax", "post_imp", "clamp"
    )
  )
)

#' Check parameter names against the method registry
#'
#' @description
#' Internal helper used by meta-functions to detect unsupported parameter names
#' supplied for K-NN or PCA imputation.
#'
#' @param param_names Character vector of parameter names to check.
#' @param mode Character scalar. Imputation mode. Must correspond to an entry in
#'   `.param_registry`, currently `"knn"` or `"pca"`.
#' @param extra Character vector of additional allowed parameter names.
#' @param arg Character scalar used in error messages to identify the source of
#'   the parameters.
#'
#' @returns
#' Invisibly returns `NULL` if all parameter names are valid. Otherwise aborts
#' with a formatted error message.
#'
#' @details
#' This helper centralizes validation of method-specific parameter names. It is
#' used by functions such as `group_imp()` and `tune_imp()` to avoid silently
#' accepting misspelled or unsupported arguments.
#'
#' @keywords internal
#' @noRd
check_unknown_params <- function(
  param_names, mode, extra = character(),
  arg = "parameters"
) {
  allowed <- unique(c(.param_registry[[mode]]$allowed, extra))
  unknown <- setdiff(param_names, allowed)
  if (length(unknown)) {
    cli::cli_abort(c(
      "{cli::qty(length(unknown))}Unknown parameter{?s} in {.arg {arg}} for {.strong {toupper(mode)}}:",
      "x" = "{fmt_trunc(unknown, 10)}",
      "i" = "Allowed: {.arg {allowed}}."
    ))
  }
}
