# -----------------------------------------------------------------------------
# Reduced base-R re-implementation of carrier::crate()
#
# Derived from carrier: https://github.com/r-lib/carrier
# Copyright (c) Posit Software, PBC
# Original implementation by the carrier authors
# Licensed under the MIT License (see inst/LICENSE.note)
# -----------------------------------------------------------------------------

#' Crate Class Check
#' @noRd
.is_crate <- function(x) inherits(x, "crate")

#' Crate Class Constructor
#' @noRd
.new_crate <- function(fn) {
  structure(fn, class = "crate")
}

#' Self-Contained (Crated) Function
#'
#' Reduced base-R reimplementation of `carrier::crate()`. Isolates `.fn` in a
#' child of `.parent_env`, so that only the supplied objects travel with it to
#' parallel workers.
#'
#' @param .fn A function. The closure to isolate.
#' @param ... Named objects. Each name is unique, and each value is packaged
#'   with `.fn` under that name.
#' @param .parent_env An environment. The parent of the crate's environment.
#'   Defaults to [baseenv()].
#' @details A closure among the `...` objects is moved into the crate's
#'   environment too, unless its environment is a namespace or it is already
#'   a crate.
#' @returns A function of class `"crate"` whose enclosing environment holds only
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
    stop("'...' arguments to .crate() must be uniquely named", call. = FALSE)
  }

  if (typeof(.fn) != "closure") {
    stop("'.fn' passed to .crate() must be a function closure", call. = FALSE)
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
