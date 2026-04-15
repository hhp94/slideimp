# Changelog

## slideimp 1.0.0

### Breaking changes

- `group_features()` is renamed to
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  to better reflect its purpose. It now accepts a column name vector
  instead of a full matrix.

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  enforces stricter data validation. The requested feature subset must
  be a subset of both the object’s column names and the mapping data
  frame. Set `allow_unmapped = TRUE` to bypass errors when intersections
  are incomplete.

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  and
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  now error when arguments are supplied that do not apply to the chosen
  imputation method, rather than silently ignoring them.

- `inject_na()` is renamed to
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  and is now exported. The original remains accessible via
  `slideimp:::inject_na()` for legacy code.

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  uses a logical `tree` argument to toggle between Ball tree (`TRUE`)
  and brute force (`FALSE`). KD tree is no longer supported.

- [`sim_mat()`](https://hhp94.github.io/slideimp/reference/sim_mat.md)
  returns a matrix in sample-by-column format for immediate
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

### New features

- [`compute_metrics()`](https://hhp94.github.io/slideimp/reference/compute_metrics.md)
  supports data frames with a `result` list column containing truth and
  estimate columns, similar to
  [yardstick](https://github.com/tidymodels/yardstick).

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

- `{slideimp.extra}` is added as a suggested package, available via
  R-universe. It provides lightweight Illumina manifests for
  frictionless
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  calls.

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
