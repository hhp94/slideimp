# To do

Work that is agreed on but not yet done. An entry stays here until the change
lands; then it is deleted, not marked done - git history is the record of what
happened, this file is the record of what has not.

An entry says what to change and why, enough that a session with no memory of
the conversation can pick it up. If an entry turns out to be wrong, delete it
and say so in chat.

## Hoist the early returns in `knn_imp()` and `pca_imp()`

Run `/simplify` over `R/knn_imp.R` and `R/pca_imp.R`, targeting the
early-return paths.

Both functions bail out early on inputs that need no work - an empty `subset`,
no missing values in the subset, no missing values at all - and the two
functions do it differently. `knn_imp()` returns through
`new_slideimp_results()` on both of its early paths; `pca_imp()` returns a bare
matrix at `R/pca_imp.R:383-386`, so the class and attributes of its result
depend on which path the data took. That asymmetry is the thing to remove.

What the hoisted form has to preserve:

- **Return the input unchanged, with no copy.** These paths do no imputation,
  so nothing should touch the matrix. Attaching a class and attributes is a
  modification under copy-on-modify semantics, and on a 450k-column matrix the
  duplication is not free. Check whether the constructor can be made to attach
  in place, or whether the copy is unavoidable and merely needs to be
  acknowledged; measure with `tracemem()` rather than reasoning about it.
- **One message per path, through cli.** `resolve_subset()` already reports an
  empty subset, so the caller must not report it a second time. Every message
  on these paths goes through `cli::cli_inform()` with inline markup, the same
  as the rest of the package.
- **The same class and attributes as the main path.** `slideimp_results` with
  `imp_method`, `fallback`, `post_imp` and `has_remaining_na`, so a caller
  cannot tell from the result's shape which path produced it.
  `tests/testthat/test-knn_imp.R` pins this for `knn_imp()`; `pca_imp()` needs
  the same test.

Open question to settle while doing it: `has_remaining_na` is computed over the
whole matrix, but `print.slideimp_results` describes it as covering the
requested columns. Whichever way that is resolved, both functions and both
paths have to agree, and the help page has to say which one it is.
