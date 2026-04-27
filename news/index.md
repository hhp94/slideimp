# Changelog

## slideimp 1.0.1

### Breaking changes

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  no longer uses the cache; the `use_cache` argument has been removed.
- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  now runs the warm-start LOBPCG solver by default. Use
  `lobpcg_control(maxiter = 0)` to go back to the full `dsyevr` solver.

### Minor improvements and fixes

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  now uses [RcppThread](https://github.com/tnagler/RcppThread) instead
  of OpenMP for macOS support.

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  now infers the subset from `na_loc` to speed up tuning for
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  and
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md).

- [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  is now an S3 generic instead of using the register-on-load pattern
  with `{slideimp.extra}`.

- Fixed CRAN ATLAS numerical tolerance check.

## slideimp 1.0.0

CRAN release: 2026-04-16

### Breaking changes

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  enforces stricter data validation. The requested feature subset must
  be a subset the object’s column names which must be a subset of the
  mapping `data.frame`. Set `allow_unmapped = TRUE` to bypass errors
  when intersections are incomplete.

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  and
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  now error when arguments are supplied that do not apply to the chosen
  imputation method, rather than silently ignoring them.

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  now uses a logical `tree` argument to toggle between Ball tree
  (`TRUE`) and brute force (`FALSE`). KD tree is no longer supported.

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  and
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  gain more early errors and early exits.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  gains the same `colmax` and `post_imp` arguments as
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md).

- [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  (formerly `group_features()`) is the new name for the grouping
  function. It now accepts a column name vector instead of a full
  matrix.

- [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  (formerly `inject_na()`) is now exported. The original remains
  accessible via `slideimp:::inject_na()` for legacy code.

- [`sim_mat()`](https://hhp94.github.io/slideimp/reference/sim_mat.md)
  now returns a matrix in sample-by-column format for immediate
  compatibility with other package functions. `perc_NA` is renamed to
  `perc_total_na`, and dimensions are now specified via `n` (rows) and
  `p` (columns).

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  gains a unified `method` argument that applies to both
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  and
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
  replacing `pca_method` and `knn_method`. The `rep` argument is renamed
  to `n_reps`.

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  results from v0.5.4 are no longer reproducible because internal NA
  generation now uses
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md).

- The `khanmiss1` dataset has been removed.

### New features

- [`compute_metrics()`](https://hhp94.github.io/slideimp/reference/compute_metrics.md)
  now supports data frames with a `result` list column containing truth
  and estimate columns, similar to
  [yardstick](https://github.com/tidymodels/yardstick).

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  and
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  automatically look up Illumina manifests using the register-on-load
  pattern for `{slideimp.extra}`.

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  gains `max_cache` to control the internal cache size (defaults to
  4GB).

- [`sim_mat()`](https://hhp94.github.io/slideimp/reference/sim_mat.md)
  gains a `rho` argument to support compound symmetry correlation
  structures in simulated matrices.

- [`sim_mat()`](https://hhp94.github.io/slideimp/reference/sim_mat.md)
  and
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  gain dedicated print methods that provide concise summaries instead of
  dumping raw data to the console.

- [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  gains `location`, `flank`, and `dry_run` arguments for fixed-window
  imputation, “flank mode” for features surrounding a subset, and
  pre-computation inspection of window statistics.

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  gains granular control over NA injection via `n_cols`, `n_rows`,
  `num_na`, and `na_col_subset`. Pre-calculated locations can also be
  passed to `na_loc` to compare methods using identical NA patterns.

### Minor improvements and fixes

- [`col_vars()`](https://hhp94.github.io/slideimp/reference/col_vars.md)
  and
  [`mean_imp_col()`](https://hhp94.github.io/slideimp/reference/mean_imp_col.md)
  have been overhauled to use the faster
  [RcppArmadillo](https://github.com/RcppCore/RcppArmadillo) backend and
  now support parallel computation with OpenMP.

- Dependencies are streamlined. [tibble](https://tibble.tidyverse.org/)
  and [purrr](https://purrr.tidyverse.org/) are removed as hard
  dependencies, [cli](https://cli.r-lib.org) is added for more
  informative messaging, and [carrier](https://github.com/r-lib/carrier)
  is added as an explicit dependency.

- Documentation is thoroughly overhauled with numerous consistency
  improvements and bug fixes.

- [RhpcBLASctl](https://prs.ism.ac.jp/~nakama/Rhpc/) is added as a
  suggested package to allow pinning BLAS cores and avoid thrashing
  during parallel runs.

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  and
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  prioritize process-level parallelization via
  [mirai](https://mirai.r-lib.org).
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  supports OpenMP-controlled parallelization via the `cores` argument
  when [mirai](https://mirai.r-lib.org) daemons are not active.

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  and
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  use optimized internal Rcpp functions for better performance.

## slideimp 0.5.4

CRAN release: 2026-01-07

- CRAN resubmission.

- `group_features()` is added to help with creating the group tibble
  needed for
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  now allows `row.w = "n_miss"` to scale row weights by the number of
  missing values per row.

## slideimp 0.5.3

- Initial CRAN submission.
