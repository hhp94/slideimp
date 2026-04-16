.slideimp_env <- new.env(parent = emptyenv())

#' Register a Group Resolver
#'
#' Called by the companion package `slideimp.extra` from their `.onLoad()`
#' to register a callback that turns character manifest into a `data.frame` for
#' `group_imp()`/`prep_groups()`.
#'
#' @param resolver A function that takes the manifest character and returns
#' a `data.frame`.
#'
#' @returns Invisibly returns `NULL`. The function is called only for its
#' side effect of registering the resolver.
#'
#' @export
#' @examples
#' # Typically called from `.onLoad` by `slideimp.extra` package, not by users directly.
#' dummy_resolver <- function(x) data.frame(feature = character(), group = character())
#' register_group_resolver(dummy_resolver)
register_group_resolver <- function(resolver) {
  checkmate::assert_function(resolver, nargs = 1)
  .slideimp_env$group_resolver <- resolver
}
