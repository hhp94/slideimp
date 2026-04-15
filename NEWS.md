# slideimp 1.0.0

## Breaking changes

* `group_features()` is renamed to `prep_groups()` to better reflect its
  purpose. It now accepts a column name vector instead of a full matrix.

* `group_imp()` enforces stricter data validation. The requested feature
  subset must be a subset of both the object's column names and the
  mapping data frame. Set `allow_unmapped = TRUE` to bypass errors when
  intersections are incomplete.

* `group_imp()` and `tune_imp()` now error when arguments are supplied
  that do not apply to the chosen imputation method, rather than
  silently ignoring them.

* `inject_na()` is renamed to `sample_na_loc()` and is now exported. The
  original remains accessible via `slideimp:::inject_na()` for legacy
  code.

* `knn_imp()` uses a logical `tree` argument to toggle between Ball tree
  (`TRUE`) and brute force (`FALSE`). KD tree is no longer supported.

* `sim_mat()` returns a matrix in sample-by-column format for immediate
  compatibility with other package functions. `perc_NA` is renamed to
  `perc_total_na`, and dimensions are now specified via `n` (rows) and
  `p` (columns).

* `tune_imp()` gains a unified `method` argument that applies to both
  `pca_imp()` and `knn_imp()`, replacing `pca_method` and `knn_method`.
  The `rep` argument is renamed to `n_reps`.

* `tune_imp()` results from v0.5.4 are no longer reproducible because
  internal NA generation now uses `sample_na_loc()`.

## New features

* `compute_metrics()` supports data frames with a `result` list column
  containing truth and estimate columns, similar to `{yardstick}`.

* `knn_imp()` gains `max_cache` to control the internal cache size
  (defaults to 4GB).

* `sim_mat()` gains a `rho` argument to support compound symmetry
  correlation structures in simulated matrices.

* `sim_mat()` and `tune_imp()` gain dedicated print methods that provide
  concise summaries instead of dumping raw data to the console.

* `slide_imp()` gains `location`, `flank`, and `dry_run` arguments for
  fixed-window imputation, "flank mode" for features surrounding a
  subset, and pre-computation inspection of window statistics.

* `tune_imp()` gains granular control over NA injection via `n_cols`,
  `n_rows`, `num_na`, and `na_col_subset`. Pre-calculated locations can
  also be passed to `na_loc` to compare methods using identical NA
  patterns.

## Minor improvements and fixes

* Dependencies are streamlined. `{tibble}` and `{purrr}` are removed as
  hard dependencies, `{cli}` is added for more informative messaging,
  and `{carrier}` is added as an explicit dependency.

* Documentation is thoroughly overhauled with numerous consistency
  improvements and bug fixes.

* `{RhpcBLASctl}` is added as a suggested package to allow pinning BLAS
  cores and avoid thrashing during parallel runs.

* `{slideimp.extra}` is added as a suggested package, available via
  R-universe. It provides lightweight Illumina manifests for frictionless
  `group_imp()` calls.

* `group_imp()` and `tune_imp()` prioritize process-level parallelization
  via `{mirai}`. `knn_imp()` supports OpenMP-controlled parallelization
  via the `cores` argument when `{mirai}` daemons are not active.

* `knn_imp()` and `pca_imp()` use optimized internal Rcpp functions for
  better performance.
