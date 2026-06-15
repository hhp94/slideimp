# -----------------------------------------------------------------------------
# Reduced base-R re-implementation of carrier::crate()
#
# Derived from carrier: https://github.com/r-lib/carrier
# Copyright (c) Posit Software, PBC
# Original implementation by the carrier authors
# Licensed under the MIT License (see inst/LICENSE.note)
# -----------------------------------------------------------------------------

#' Test whether an object is a crate
#' @noRd
.is_crate <- function(x) inherits(x, "crate")

#' Tag a function as a crate
#' @noRd
.new_crate <- function(fn) {
  structure(fn, class = "crate")
}

#' Create a self-contained ("crated") function
#'
#' Reduced base-R reimplementation of `carrier::crate()` for internal use.
#' Isolates `.fn` in a child of `.parent_env` so that only the explicitly
#' supplied objects are packaged with it when sent to parallel workers.
#'
#' @param .fn A function closure to isolate.
#' @param ... Named, already-evaluated objects to package with `.fn`.
#' @param .parent_env Parent environment of the crate's environment; defaults
#'   to [baseenv()] so nothing from the calling scope leaks in.
#' @return A function of class `"crate"` whose enclosing environment holds only
#'   the supplied objects.
#'
#' @keywords internal
#' @noRd
.crate <- function(.fn, ..., .parent_env = baseenv()) {
  vars <- list(...)
  nms <- names(vars)

  if (
    length(vars) &&
      (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms))
  ) {
    stop("All `...` arguments must be uniquely named.", call. = FALSE)
  }

  if (typeof(.fn) != "closure") {
    stop("`.fn` must be a function closure.", call. = FALSE)
  }

  env <- new.env(parent = .parent_env)
  if (length(vars)) {
    list2env(vars, envir = env)
  }

  for (nm in nms) {
    x <- env[[nm]]
    if (
      typeof(x) == "closure" &&
        !isNamespace(environment(x)) &&
        !.is_crate(x)
    ) {
      environment(x) <- env
      env[[nm]] <- utils::removeSource(x)
    }
  }

  environment(.fn) <- env
  .new_crate(utils::removeSource(.fn))
}
