# Decisions

Closed questions, stated as facts about the code. Each was weighed and settled,
so a session does not reopen it on its own initiative: not as a to-do entry,
not as a review finding, not as a drive-by fix. The maintainer can reopen one;
when that happens, the entry is changed or deleted in the same change as the
code.

Work still owed lives in [to-do.md](to-do.md), never here.

## Kernels

- The K-NN kernel's `obj_masked` and `nmiss_masked`, the buffers of
  `src/pca_kernel.cpp`, and each new `Block` in `src/lobpcg_warm.h` keep
  their zero-fill by `std::vector`.
- The fallback's column means stay in R, `.colMeans()` in `fallback_part()`
  (`R/parts.R`). They are within one ulp of the exact mean, and nothing
  requires them bitwise equal to the 1.2.0 kernel's.
- The PCA kernel's own peak memory above the input is not pursued.

## Sampler

The multiple-imputation sampler is `R/mi_pca_*.R`, `src/mipca*.cpp` and
`tests/testthat/test-mipca-*.R`.

Rulings:

- The sampler's refusals stay plain `stop(..., call. = FALSE)`, and its
  messages and warnings stay base `message()` and `warning()`, in its files
  and in the `test-mipca-*` files that pin their wording.
- No bytes warning before allocating, and no refusal at `2^31` cells. The
  draws refuse in plain words when R cannot allocate them, naming the size
  in GB and what makes it up. File-backed draws are parked until a user hits
  that refusal.
- `subset` keeps its `NULL` default at `mipca()`, meaning the feature is
  off: draws for every column with a hole. In flank mode it is required, and
  `mipca_flank_plan()` and `mipca_flank_window()` refuse a missing or `NULL`
  `subset`. The `@param subset` of `mipca()` says `NULL` is "every column"
  rather than every column with a hole, and stays so.
- `cbind.si_mi()` joins finalized fits only: identical row names in the
  same order, the same `nkeep`, each name once; no positional fallback.
- Flank keeps its two calls, `mipca_flank_plan()` then `mipca_flank()`, and
  `slide_imp()` keeps its one call with `dry_run`. Neither is made to look
  like the other.
- No GCV scoring in `tune_imp()`.
- `mipca()` sets no column aside to fill outside its fit: the core refuses a
  column that holds `Inf`, is all `NA`, or has one observed value. `pca_imp()` refuses
  the first two and fills the third with its one value; in `mipca()` that
  column would take the same value in every draw, so the refusal points at
  `abb` instead, where the caller asks for it.
- With `subset`, the check for nothing to impute reads only the `subset`
  columns, so a fit with nothing to impute accepts an `Inf` outside
  `subset`. A fit that samples reads every column, since its model is fit on
  all of them. The check is `.nothing_to_impute()` (`R/mi_pca_imp.R`:
  `.nothing_to_impute <- function(obj, subset, scan = NULL) {`).
- Flank reads only the columns a window can hold, and the `abb` targets, so
  `mipca_flank_plan()` and `mipca_flank()` accept an `Inf` in a column no
  window holds. The plan's read set is built in `.flank_scan()`
  (`R/mi_pca_flank.R`:
  `span <- find_windows_flank(location, tg[!tg_abb], window_size)`).

Constraints the code keeps; a change to one is deliberate:

- The sampler starts no thread of its own and pins no BLAS, and its pages
  never tell users to pin BLAS threads. If `mipca_tune()` ever runs its
  arms in parallel, `tune_imp()`'s guard (`cores` to 1 under live mirai
  daemons) comes with it.
- `.resolve_cols()` refuses an unknown `subset` name itself, before
  `resolve_subset()`, which informs and drops one, because `subset` decides
  the fit's shape.
- The constant-column test, variance below `.Machine$double.eps` on the
  scan's variance (`R/col_scan.R`), is shared by `pca_imp()`,
  `sample_na_loc()` and the sampler's `.admissible_scan()`; they change
  together.
- The scan's counts are integer as `col_scan_cpp()` returns them, and
  nothing coerces them (`src/mat_stats.cpp`: `Rcpp::List col_scan_cpp(`).
- A change to the state's fields, or to how `g_fixed` is formed, bumps
  `.STATE_VERSION`; an older state is refused, with no migration. A resume
  starts each chain's first eigensolve from the stored `eig_start`, and
  `resuming` stays a named argument of `.mipca_da_loop_multi()`, never
  inferred from empty state vectors.
- A failing flank window stops the run; a window is seeded by the seed and
  its target's column index in the plan's matrix.
- Every BLAS call goes through `R_ext/BLAS.h` with `USE_FC_LEN_T` defined
  before any R header and one `FCONE` per character argument. Progress
  lines are formatted by `tfm::format` and sent through `Rcpp::message()`,
  so they are R messages.

## CI

The run tiers and platforms of a release are in [ci.md](ci.md).

- ubuntu-release is in the minimal run tier: it is the one full check,
  documentation included. The sanitizer and valgrind images run R CMD check
  with --extra-arch, which skips most checks.
- The full run tier runs the CRAN test tier on ubuntu-release, nold, atlas,
  mkl and blis: CRAN's extra checks leave NOT_CRAN unset, so on their
  machines only `tests/testthat/test-smoke.R` runs, and these are the
  platforms where its numbers can move.
- These of CRAN's extra checks get no run. donttest: slideimp's examples
  have no \donttest or \dontrun. noOMP: `src/Makevars` and
  `src/Makevars.win` use no OpenMP. BLAS: `src/mipca.cpp` passes FCONE, and
  any install's log shows a fault. C23, Strict and noRemap: R-devel's
  defaults, so every R-devel image covers them. Intel and 0len: CRAN no
  longer runs them.
- rchk's one error, "possible protection stack imbalance" in
  `Rcpp::Rcpp_protect` (`Rcpp/protection/Shield.h`), stays. `Shield`
  protects in its constructor and unprotects in its destructor, and Rcpp
  routes every such call through `Rcpp_protect()` so that rchk reports them
  at one place (the comment in `Shield.h` cites Rcpp's #892). Every Rcpp
  vector, `wrap()`, `as()` and `BEGIN_RCPP` uses it, the generated wrappers
  in `src/RcppExports.cpp` included, so only C++ without Rcpp would remove
  it. CRAN's filter drops it: CRAN's unfiltered rchk results for slideimp
  1.2.0 hold it, as dqrng's do, and neither package is in CRAN's rchk list.
  Every other rchk error is fixed, not kept.

## Releases

- No 1.2.x patch release on CRAN. The fix for the 1.2.0 `pca_imp()`
  segfault ships in the next release.
