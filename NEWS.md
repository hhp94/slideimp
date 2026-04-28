# slideimp 1.1.0

## Breaking changes

* `knn_imp()` no longer caches pair-wise distances. The `use_cache` argument has been removed.

## New features

* `pca_imp()` gains support for the LOBPCG eigensolver, selected via the new `solver` argument and tuned with `lobpcg_control()`.

## Minor improvements and fixes

* `knn_imp()` now uses `{RcppThread}` instead of OpenMP, enabling parallelism on macOS.

* `tune_imp()` now infers the subset from `na_loc`, speeding up tuning for `knn_imp()` and `slide_imp()`.

* `prep_groups()` is now an S3 generic, replacing the previous register-on-load pattern with `{slideimp.extra}`.

* New article describing how to speed up PCA imputation by choosing `solver` and `threshold`.

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
