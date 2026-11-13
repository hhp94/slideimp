# S3 class ---------------------------------------------------------------------

#' @noRd
#' @keywords internal
new_slideimp_results <- function(
  obj,
  imp_method,
  fallback,
  post_imp,
  mean_filled_cols = integer(0),
  mean_filled_cells = mean_filled_cells_matrix()
) {
  # start from a bare matrix. The input can itself be a `slideimp_results`,
  # from an earlier imputation pass: left as it is, the class would be attached
  # a second time and the attributes this constructor does not set
  # (`metacaller`, `fallback_action`) would stay and describe that earlier call.
  #
  # Attaching anything to a matrix the caller still holds copies it, because
  # the caller's matrix must not gain the class. `attr<-`, `class<-` and
  # `attributes<-` each copy a shared 76 MB matrix once, as measured, so
  # no rewrite of this body avoids it. On an early return, where `obj` is the
  # caller's own input, that one copy is what returning the class costs.
  attributes(obj) <- list(dim = dim(obj), dimnames = dimnames(obj))
  class(obj) <- c("slideimp_results", class(obj))
  attr(obj, "imp_method") <- imp_method
  attr(obj, "fallback") <- fallback
  attr(obj, "post_imp") <- post_imp
  # what `post_imp` set to a column mean, in the two kinds it fills. Both are
  # always present, empty when nothing was filled, so that a caller reads them
  # without asking first whether they exist.
  attr(obj, "mean_filled_cols") <- mean_filled_cols
  attr(obj, "mean_filled_cells") <- mean_filled_cells
  obj
}

# The `mean_filled_cells` attribute: one row per cell, columns `row` and `col`,
# ordered by column and then by row. Duplicated cells are kept once: slide_imp()
# can meet the same cell in two overlapping windows. `n_row` is the number of
# rows of the matrix the cells index, needed only to tell duplicates apart.
mean_filled_cells_matrix <- function(
  row = integer(0),
  col = integer(0),
  n_row = 0L
) {
  if (length(row) > 1L) {
    # a double key: `col * n_row` overflows an integer on a large matrix
    key <- (as.numeric(col) - 1) * n_row + row
    o <- order(key)
    o <- o[!duplicated(key[o])]
    row <- row[o]
    col <- col[o]
  }
  cbind(row = as.integer(row), col = as.integer(col))
}

# What one backend result recorded as mean-filled, and the `post_imp` it ran
# with, for a driver that takes only some of its columns. `res` is what knn_imp() or pca_imp() returned for a
# submatrix, or a plain matrix when the backend did not run. `keep` holds the
# submatrix columns the driver takes, and `global` maps every submatrix column
# to its index in the driver's own matrix.
#
# The records of local_mean_filled(), one per window or group and NULL where
# there is none, as the two attributes of the driver's result.
combine_mean_filled <- function(records, n_row) {
  pull <- function(what) {
    unlist(lapply(records, `[[`, what), use.names = FALSE)
  }
  list(
    # as.integer(): unlist() over records that are all NULL is NULL
    cols = as.integer(sort(unique(pull("cols")))),
    cells = mean_filled_cells_matrix(
      row = pull("row"),
      col = pull("col"),
      n_row = n_row
    )
  )
}

# Base R only, and it calls nothing else from the package: group_imp() sends
# it to mirai daemons inside a crate, so it must not depend on what a daemon
# has loaded.
local_mean_filled <- function(res, keep, global) {
  cols <- attr(res, "mean_filled_cols", exact = TRUE)
  cells <- attr(res, "mean_filled_cells", exact = TRUE)
  cols <- cols[cols %in% keep]
  if (is.null(cells)) {
    cell_row <- integer(0)
    cell_col <- integer(0)
  } else {
    taken <- cells[, "col"] %in% keep
    cell_row <- unname(cells[taken, "row"])
    cell_col <- global[cells[taken, "col"]]
  }
  list(
    cols = as.integer(global[cols]),
    row = cell_row,
    col = cell_col,
    # the `post_imp` the backend ran with, NULL when it did not run
    post_imp = attr(res, "post_imp", exact = TRUE)
  )
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
      cli::cli_abort(
        "{.arg subset} contains characters but {.arg obj} has no column names."
      )
    }
    if (anyDuplicated(subset)) {
      cli::cli_abort("{.arg subset} contains duplicate feature names.")
    }
    matched <- match(subset, cn, nomatch = NA_integer_)
    if (anyNA(matched)) {
      cli::cli_inform(
        "Feature(s) in {.arg subset} not found in {.code colnames(obj)} and dropped."
      )
    }
    subset <- matched[!is.na(matched)]
  } else {
    checkmate::assert_integerish(
      subset,
      lower = 1L,
      upper = nc,
      any.missing = FALSE,
      unique = TRUE,
      min.len = 0L,
      .var.name = "`subset`"
    )
    subset <- as.integer(subset)
  }

  if (length(subset) == 0L) {
    cli::cli_inform(
      "No features in {.arg subset} detected. No imputation was performed."
    )
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
    cli::cli_inform(
      "Setting {.code pin_blas = TRUE} may improve parallel performance."
    )
  }
}

# mirai ----
#' Are mirai Daemons Set?
#'
#' The one place the drivers ask. `mirai::daemons_set()` is the predicate.
#' `mirai::require_daemons()` is the prompt built on it: it raises an error
#' when there are none, and it is newer.
#'
#' @returns A logical scalar.
#'
#' @noRd
#' @keywords internal
mirai_daemons_set <- function() {
  mirai::daemons_set()
}

#' Collect a mirai Map
#'
#' `[.mirai_map` evaluates its arguments among mirai's own collection options,
#' so `mapped[.progress = .progress]` asks for the progress bar whatever the
#' caller's flag holds. Ask for it by name, or not at all.
#'
#' @param mapped A `mirai_map` object.
#' @param .progress Logical scalar. Whether to show mirai's progress bar.
#'
#' @returns A list with one element per task. A task that failed outside its
#'   own error handling is a mirai error value: see `mirai_error_text()`.
#'
#' @noRd
#' @keywords internal
collect_mirai_map <- function(mapped, .progress) {
  if (.progress) mapped[.progress] else mapped[]
}

#' One Line of Text for a mirai Error Value
#'
#' An R error raised in a task comes back as a `miraiError`, which carries the
#' message. A task whose daemon stopped comes back as a bare integer error
#' code with no `conditionMessage()` method. Its print method names the code.
#'
#' @param x An object for which `mirai::is_error_value()` is `TRUE`.
#'
#' @returns A character scalar.
#'
#' @noRd
#' @keywords internal
mirai_error_text <- function(x) {
  if (inherits(x, "miraiError")) {
    return(conditionMessage(x))
  }
  paste(utils::capture.output(print(x)), collapse = " ")
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
      "k",
      "method",
      "colmax",
      "post_imp",
      "dist_pow",
      ".progress"
    )
  ),
  pca = list(
    required = "ncp",
    methods = c("regularized", "EM"),
    allowed = c(
      "ncp",
      "scale",
      "method",
      "coeff.ridge",
      "row.w",
      "threshold",
      "seed",
      "nb.init",
      "maxiter",
      "miniter",
      "lobpcg_control",
      "solver",
      "colmax",
      "post_imp",
      "clamp",
      ".progress"
    )
  )
)

#' Parameter names `tune_imp()` accepts for `slide_imp()`
#'
#' @description
#' The window arguments plus everything the registry allows for K-NN and PCA,
#' so an argument added to the registry is tunable through `slide_imp()` too.
#'
#' @returns A character vector.
#'
#' @keywords internal
#' @noRd
slide_tune_params <- function() {
  unique(c(
    "window_size",
    "overlap_size",
    "min_window_n",
    "k",
    "ncp",
    "flank",
    "on_infeasible",
    .param_registry$knn$allowed,
    .param_registry$pca$allowed
  ))
}

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
  param_names,
  mode,
  extra = character(),
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
