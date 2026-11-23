# slideimp 1.3.0

## Breaking changes

* `Depends` is now `R (>= 4.4.0)` (was 4.3.0), for C++17.

* Every function that takes `obj` requires a double matrix. 1.2.0 converted
integer and logical matrices.

* The `na_check` argument of `knn_imp()`, `pca_imp()`, `slide_imp()` and
`group_imp()` is removed. Use `anyNA()` on the result instead.

* `mean_imp_col()` is removed.

* `col_vars()` loses its `cores` argument and runs on one thread.

* `subset` now narrows the result in `knn_imp()`, `slide_imp()` and
`group_imp()`: the result holds only the selected columns, in input order.
In `knn_imp()`, an unknown name, a repeated name, names on an input without
column names, and an empty `subset` are errors; 1.2.0 dropped an unknown
name with a message and returned the input for an empty `subset`.
`slide_imp()` refuses an unknown name and an empty `subset`, under
`dry_run = TRUE` too. `tune_imp()` with `.f = "slide_imp"` scores the same
cells as before.

* `mat_miss()` returns an integer vector when `prop = FALSE` (still a double
vector when `prop = TRUE`).

* `tune_imp()` errors when a supplied `na_loc` position points at a cell that
is already `NA` or `NaN`. Positions from `sample_na_loc()` are unaffected.

* Many errors, messages and `print()` notes are reworded. Code that matches
on message text may need updating.

## New features

* Multiple imputation by Bayesian PCA, a reimplementation of the sampler of
`missMDA::MIPCA()` (Audigier, Husson and Josse, 2016):

  * `mipca()` draws multiple imputations of the missing values.
  * `mipca_tune()` chooses the rank on held-out values, and `mipca_fit()`
    fits at the tuned rank.
  * `mipca_finalize()` and `update()` extend or finalize a resumable fit.
  * `mipca_flank()`, `mipca_flank_plan()` and `mipca_flank_window()` impute
    each column on a flanking window of its neighbors.
  * `tf_identity()`, `tf_sv()` and `tf_biwhiten()` transform the data before
    the fit; `spectra_opts()` sets the eigensolver; `abb_cols()` proposes
    columns to draw by approximate Bayesian bootstrap.
  * Results have `print()`, `summary()`, `as.matrix()`, `[` and `cbind()`
    methods.

* `pca_imp()` gains `subset`: the result holds only the selected columns,
while every eligible column still takes part in the fit.

## Changes in results

* `knn_imp()` with `post_imp = TRUE` sets a cell that K-NN leaves missing to
the mean of its column's observed values. 1.2.0 took a mean that included
the K-NN imputed values. These cells change value, in `slide_imp()` and
`group_imp()` with K-NN too.

* The column means that `post_imp` and `on_infeasible = "mean"` set are now
computed with base R's `.colMeans()`, so they can differ from 1.2.0's in the
last bits.

* `group_imp()`: every group now reads `obj` as it was passed. 1.2.0's
sequential path wrote each group's imputed values into `obj` before the next
group ran, so a later group could use them as observed values in its
auxiliary columns.

* `slide_imp()` and `group_imp()` with PCA: a window or group whose missing
target values all lie in columns that `pca_imp()` leaves out is infeasible
and goes through `on_infeasible`, as with K-NN. 1.2.0 fitted it and set the
targets to their column means.

* `slide_imp()` returns the observed cells of `obj` unchanged. Where three or
more windows overlapped, 1.2.0 averaged every window's copy of a column,
which could move an observed value in its last bit.

* `pca_imp()` now calls R's BLAS and LAPACK directly instead of going through
Armadillo, so its results differ from 1.2.0's in the last digits under the
exact solver, and by up to about 1e-8 under LOBPCG.

* `pca_imp()` on a matrix with no missing values returns a
`slideimp_results` object like any other call, where 1.2.0 returned the
input unchanged.

## Minor improvements and fixes

* `pca_imp()` no longer crashes the R session when `ncp` exceeds the rank of
the completed data, as with replicate samples.

* `pca_imp()`: finite input whose squares overflow stops with "the
cross-product of the data is not finite; values this large overflow when
squared, rescale the input", before either solver. A LAPACK failure reads
"the symmetric eigensolver failed (LAPACK dsyevr, info = <n>)". The 1.2.0
error "`eig_sym_sel` failed to converge" is gone.

* `pca_imp()` refuses eligible columns holding more than
`.Machine$integer.max` missing values.

* `pca_imp()` scans for infinities once per call rather than once per
`nb.init` initialization. When the input has more columns than rows, only
the loadings of the columns being imputed are updated inside the EM loop.

* `pin_blas` in `group_imp()` and `tune_imp()` now defaults to `TRUE`,
pinning BLAS and OpenMP threads to 1 on mirai daemons; a sequential run pins
nothing. RhpcBLASctl moves from `Suggests` to `Imports`.

* `tune_imp()` under mirai sends `obj` to each daemon once, with
`mirai::everywhere()`, so memory grows with the number of daemons.
bigmemory is no longer imported, and mirai must be `>= 2.5.0`. A new error
reports when `obj` did not reach every daemon.

* `tune_imp()` with `.f = "knn_imp"` or `"slide_imp"` fails an iteration when
the result is not the width of the repetition's `subset`. That happens when a
daemon runs an older slideimp; 1.2.0 read the wrong cells without an error.

* `sample_na_loc()`, and so `tune_imp()`, no longer draws a column whose
variance is above zero but below `.Machine$double.eps`, the threshold
`pca_imp()` uses for eligibility.

# slideimp 1.2.0

## Breaking changes

* `knn_imp()` no longer supports Ball Tree-based imputation due to lower 
accuracy compared to brute-force K-NN.

## New features

* `pca_imp()` gained a `.progress` argument to display imputation progress.

## Minor improvements and fixes

* Fixed `pca_imp()` tests that broke due to upstream changes in
`missMDA::imputePCA()` and `FactoMineR::svd.triplet()`.

* Implemented an internal `.crate()` function, a minimal version of 
`carrier::crate()`.

* `DESCRIPTION` now describes the LOBPCG eigensolver and its warm-start
behavior across EM iterations in `pca_imp()`.

* Minor wording improvements.

* Added citation information.

# slideimp 1.1.0

## Breaking changes

* Minimum R version is now 4.3.0 (was 4.1.0).
* `knn_imp()` no longer caches pair-wise distances. The `use_cache` argument
  has been removed.

## New features

* `pca_imp()` gained a `solver` argument to select the LOBPCG eigensolver,
  configured via the new `lobpcg_control()` helper.

* `pca_imp()` gained a `clamp` argument to bound imputed values, e.g.
  `c(0, 1)` for methylation beta values. Also threaded through `slide_imp()`
  and `group_imp()`.

## Minor improvements and fixes

* `knn_imp()` now uses `{RcppThread}` instead of OpenMP, enabling parallelism
  on macOS.

* `tune_imp()` now infers the subset from `na_loc`, speeding up tuning for
  `knn_imp()` and `slide_imp()`.

* `prep_groups()` is now an S3 generic, replacing the previous register-on-load
  pattern with `{slideimp.extra}`.

* New article describing how to speed up PCA imputation by choosing `solver`
  and `threshold`.

* Fixed numerical tolerance check failure on CRAN's ATLAS configuration.

# slideimp 1.0.0

## Breaking changes

* `group_imp()` enforces stricter data validation. The requested feature subset 
must be a subset the object's column names which must be a subset of the mapping 
`data.frame`. Set `allow_unmapped = TRUE` to bypass errors when intersections are
incomplete.

* `group_imp()` and `tune_imp()` now error when arguments are supplied that do 
not apply to the chosen imputation method, rather than silently ignoring them.

* `knn_imp()` now uses a logical `tree` argument to toggle between Ball tree 
(`TRUE`) and brute force (`FALSE`). KD tree is no longer supported.

* `knn_imp()` and `pca_imp()` gain more early errors and early exits.

* `pca_imp()` gains the same `colmax` and `post_imp` arguments as `knn_imp()`.

* `prep_groups()` (formerly `group_features()`) is the new name for the grouping
function. It now accepts a column name vector instead of a full matrix.

* `sample_na_loc()` (formerly `inject_na()`) is now exported. The original 
remains accessible via `slideimp:::inject_na()` for legacy code.

* `sim_mat()` now returns a matrix in sample-by-column format for immediate 
compatibility with other package functions. `perc_NA` is renamed to 
`perc_total_na`, and dimensions are now specified via `n` (rows) and `p` (columns).

* `tune_imp()` gains a unified `method` argument that applies to both 
`pca_imp()` and `knn_imp()`, replacing `pca_method` and `knn_method`. 
The `rep` argument is renamed to `n_reps`.

* `tune_imp()` results from v0.5.4 are no longer reproducible because internal
NA generation now uses `sample_na_loc()`.

* The `khanmiss1` dataset has been removed.

## New features

* `compute_metrics()` now supports data frames with a `result` list column 
containing truth and estimate columns, similar to `{yardstick}`.

* `group_imp()` and `prep_groups()` automatically look up Illumina manifests
using the register-on-load pattern for `{slideimp.extra}`.

* `knn_imp()` gains `max_cache` to control the internal cache size 
(defaults to 4GB).

* `sim_mat()` gains a `rho` argument to support compound symmetry correlation 
structures in simulated matrices.

* `sim_mat()` and `tune_imp()` gain dedicated print methods that provide concise
summaries instead of dumping raw data to the console.

* `slide_imp()` gains `location`, `flank`, and `dry_run` arguments for 
fixed-window imputation, "flank mode" for features surrounding a subset, and 
pre-computation inspection of window statistics.

* `tune_imp()` gains granular control over NA injection via `n_cols`, `n_rows`,
`num_na`, and `na_col_subset`. Pre-calculated locations can also be passed to
`na_loc` to compare methods using identical NA patterns.

## Minor improvements and fixes

* `col_vars()` and `mean_imp_col()` have been overhauled to use the faster 
`{RcppArmadillo}` backend and now support parallel computation with OpenMP.

* Dependencies are streamlined. `{tibble}` and `{purrr}` are removed as hard 
dependencies, `{cli}` is added for more informative messaging, and `{carrier}`
is added as an explicit dependency.

* Documentation is thoroughly overhauled with numerous consistency improvements 
and bug fixes.

* `{RhpcBLASctl}` is added as a suggested package to allow pinning BLAS cores
and avoid thrashing during parallel runs.

* `group_imp()` and `tune_imp()` prioritize process-level parallelization via
`{mirai}`. `knn_imp()` supports OpenMP-controlled parallelization via the `cores` 
argument when `{mirai}` daemons are not active.

* `knn_imp()` and `pca_imp()` use optimized internal Rcpp functions for better 
performance.

# slideimp 0.5.4

* CRAN resubmission.

* `group_features()` is added to help with creating the group tibble needed for
`group_imp()`.

* `pca_imp()` now allows `row.w = "n_miss"` to scale row weights by the number 
of missing values per row.

# slideimp 0.5.3

* Initial CRAN submission.
