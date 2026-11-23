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
  # Attaching anything to a matrix the caller still holds copies it once, so
  # on an early return, where `obj` is the caller's own input, that copy is
  # what returning the class costs.
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

# One record of what `post_imp` set to a column mean, from a backend's
# `info`, for combine_mean_filled(). The backend numbers `mean_filled_cols`
# and the `col` of `mean_filled_cells` as its own result's columns; `to[t]`
# is the column of the driver's result that the backend's `t`-th result
# column becomes. `post_imp` is the value the backend ran with, which
# group_imp() reads back.
mean_filled_record <- function(info, to) {
  cells <- info$mean_filled_cells
  list(
    cols = to[info$mean_filled_cols],
    row = unname(cells[, "row"]),
    col = to[cells[, "col"]],
    post_imp = info$post_imp
  )
}

# The records of what `post_imp` set to a column mean, one per window or
# group and NULL where there is none, as the two attributes of the driver's
# result. Every record is built by mean_filled_record() and numbered as the
# driver's result's columns.
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

# Formatting -------------------------------------------------------------------

#' Truncated Comma-Separated List for Messages
#'
#' Joins the first `n` elements of `x` into one comma-separated string, with
#' `, ...` appended when elements were cut, so that a message stays readable.
#'
#' @param x An atomic vector. The values to list.
#' @param n Integer. The maximum number of elements to show. Defaults to `5`.
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

#' Imputation Method Label for Messages
#'
#' @param x Character. An imputation method code, such as `"knn"` or `"pca"`.
#'
#' @returns A length-1 character string: `"K-NN"` for `"knn"`, otherwise `x`
#'   in upper case.
#'
#' @details
#' Messages and printed output spell the method as the help pages do, so
#' `"knn"` is `"K-NN"` and never `"KNN"`.
#'
#' @noRd
#' @keywords internal
method_label <- function(x) {
  if (identical(x, "knn")) "K-NN" else toupper(x)
}

# Input validation -------------------------------------------------------------

#' Column Indices from a Subset Argument
#'
#' Turns `subset` into integer column indices of `obj`. Unknown names are
#' dropped with a message, and an empty selection gives `NULL` with a message.
#'
#' @param subset Optional character or integer vector. The columns to select,
#'   by name or by number, or `NULL` for every column.
#' @param obj A numeric matrix. The matrix whose columns `subset` selects.
#' @param sort Logical. Whether to sort the indices. Defaults to `FALSE`.
#'
#' @returns An integer vector of column indices, or `NULL` when `subset`
#'   selects no column.
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
        c(
          "{.arg subset} names columns, but {.arg obj} has no column names.",
          "i" = "Set {.code colnames(obj)}, or give column numbers."
        )
      )
    }
    if (anyDuplicated(subset)) {
      cli::cli_abort(
        c(
          "{.arg subset} names a column more than once.",
          "i" = "Name each column once."
        )
      )
    }
    matched <- match(subset, cn, nomatch = NA_integer_)
    if (anyNA(matched)) {
      cli::cli_inform(
        c(
          "{.arg subset} names {sum(is.na(matched))} column{?s} not in {.arg obj}, which {?is/are} skipped: {fmt_trunc(subset[is.na(matched)], 6)}.",
          "i" = "Use names from {.code colnames(obj)}."
        )
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
      .var.name = "subset"
    )
    subset <- as.integer(subset)
  }

  if (length(subset) == 0L) {
    cli::cli_inform(
      c(
        "{.arg subset} selects no columns, so nothing is imputed.",
        "i" = "Name at least one column of {.arg obj}, or use {.code subset = NULL} for every column."
      )
    )
    return(NULL)
  }

  if (sort) sort(subset) else subset
}

# `subset` for an imputer whose result holds only the selected columns.
# Unlike resolve_subset() above, an unknown name and an empty selection are
# refused, as mipca()'s .resolve_cols() refuses them: the result's width
# depends on `subset`, so a dropped name would be a column the caller asked
# for and does not get. Returns NULL for NULL (every column)
# and for a selection of every column, so both take the same path; otherwise
# the selected columns, ascending, which is the result's column order.
# `call` is the imputer's frame, so an error names the function the user
# called.
resolve_narrow_subset <- function(subset, obj, call = parent.frame()) {
  if (is.null(subset)) {
    return(NULL)
  }
  if (length(subset) == 0L) {
    cli::cli_abort(
      c(
        "{.arg subset} selects no columns.",
        "i" = "Name at least one column of {.arg obj}, or use {.code subset = NULL} for every column."
      ),
      call = call
    )
  }
  if (is.character(subset)) {
    cn <- colnames(obj)
    if (is.null(cn)) {
      cli::cli_abort(
        c(
          "{.arg subset} names columns, but {.arg obj} has no column names.",
          "i" = "Set {.code colnames(obj)}, or give column numbers."
        ),
        call = call
      )
    }
    dup <- unique(subset[duplicated(subset)])
    if (length(dup) > 0L) {
      cli::cli_abort(
        c(
          "{.arg subset} names {length(dup)} column{?s} more than once: {fmt_trunc(dup, 6)}.",
          "i" = "Name each column once."
        ),
        call = call
      )
    }
    idx <- match(subset, cn)
    bad <- subset[is.na(idx)]
    if (length(bad) > 0L) {
      cli::cli_abort(
        c(
          "{.arg subset} names {length(bad)} column{?s} not in {.arg obj}: {fmt_trunc(bad, 6)}.",
          "i" = "Use names from {.code colnames(obj)}."
        ),
        call = call
      )
    }
  } else {
    checkmate::assert_integerish(
      subset,
      lower = 1L,
      upper = ncol(obj),
      any.missing = FALSE,
      unique = TRUE,
      .var.name = "subset"
    )
    idx <- as.integer(subset)
  }
  if (length(idx) == ncol(obj)) NULL else sort(idx)
}

# mirai ----
#' Whether mirai Daemons Are Set
#'
#' `group_imp()` and `tune_imp()` check for daemons only through this
#' function. It calls `mirai::daemons_set()`, not `mirai::require_daemons()`,
#' which raises an error when no daemons are set.
#'
#' @returns A logical scalar.
#'
#' @noRd
#' @keywords internal
mirai_daemons_set <- function() {
  mirai::daemons_set()
}

#' Results of a mirai Map
#'
#' `[.mirai_map` evaluates its arguments among mirai's own collection options,
#' so `mapped[.progress = .progress]` asks for the progress bar whatever the
#' caller's flag holds. Ask for it by name, or not at all.
#'
#' @param mapped A `mirai_map` object.
#' @param .progress Logical. Whether to show mirai's progress bar.
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
#' @param x An object. A mirai error value, one for which
#'   `mirai::is_error_value()` is `TRUE`.
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
#' Method-Specific Parameter Registry
#'
#' The tuning parameter, the valid `method` values and the allowed parameter
#' names of each imputation mode.
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
#' `check_unknown_params()`, `slide_tune_params()` and `group_imp()` read it.
#' `tune_imp()` reaches it through the first two.
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

#' Parameter Names `tune_imp()` Accepts for `slide_imp()`
#'
#' The window arguments of `slide_imp()`, `k`, `ncp`, `on_infeasible` and every
#' name the registry allows for K-NN and PCA. An argument added to the registry
#' is then tunable through `slide_imp()` too.
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

#' Parameter Name Check Against the Method Registry
#'
#' Raises an error for a parameter name that the registry does not allow for
#' K-NN or PCA imputation.
#'
#' @param param_names Character. The parameter names to check.
#' @param mode Character. One of `"knn"` or `"pca"`, the modes in
#'   `.param_registry`.
#' @param extra Character. More parameter names to allow.
#' @param arg Character. The argument name the error message gives as the
#'   source of the parameters.
#'
#' @returns `NULL` invisibly. Raises an error for an unknown name.
#'
#' @details
#' `group_imp()` and `tune_imp()` call it, so a misspelled argument is an
#' error.
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
      "{length(unknown)} unknown parameter{?s} in {.arg {arg}} for {.strong {method_label(mode)}}:",
      "x" = "{fmt_trunc(unknown, 10)}",
      "i" = "Allowed: {.arg {allowed}}."
    ))
  }
}
