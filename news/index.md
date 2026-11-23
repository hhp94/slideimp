# Changelog

## slideimp 1.3.0

### Breaking changes

- `Depends` is now `R (>= 4.4.0)` (was 4.3.0), for C++17.

- Every function that takes `obj` requires a double matrix. 1.2.0
  converted integer and logical matrices.

- The `na_check` argument of
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md),
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  and
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  is removed. Use [`anyNA()`](https://rdrr.io/r/base/NA.html) on the
  result instead.

- `mean_imp_col()` is removed.

- [`col_vars()`](https://hhp94.github.io/slideimp/reference/col_vars.md)
  loses its `cores` argument and runs on one thread.

- `subset` now narrows the result in
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  and
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md):
  the result holds only the selected columns, in input order. In
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
  an unknown name, a repeated name, names on an input without column
  names, and an empty `subset` are errors; 1.2.0 dropped an unknown name
  with a message and returned the input for an empty `subset`.
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  refuses an unknown name and an empty `subset`, under `dry_run = TRUE`
  too.
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  with `.f = "slide_imp"` scores the same cells as before.

- [`mat_miss()`](https://hhp94.github.io/slideimp/reference/mat_miss.md)
  returns an integer vector when `prop = FALSE` (still a double vector
  when `prop = TRUE`).

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  errors when a supplied `na_loc` position points at a cell that is
  already `NA` or `NaN`. Positions from
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  are unaffected.

- Many errors, messages and
  [`print()`](https://rdrr.io/r/base/print.html) notes are reworded.
  Code that matches on message text may need updating.

### New features

- Multiple imputation by Bayesian PCA, a reimplementation of the sampler
  of [`missMDA::MIPCA()`](https://rdrr.io/pkg/missMDA/man/MIPCA.html)
  (Audigier, Husson and Josse, 2016):

  - [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md)
    draws multiple imputations of the missing values.
  - [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
    chooses the rank on held-out values, and
    [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
    fits at the tuned rank.
  - [`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
    and [`update()`](https://rdrr.io/r/stats/update.html) extend or
    finalize a resumable fit.
  - [`mipca_flank()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md),
    [`mipca_flank_plan()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
    and
    [`mipca_flank_window()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
    impute each column on a flanking window of its neighbors.
  - [`tf_identity()`](https://hhp94.github.io/slideimp/reference/tf_identity.md),
    [`tf_sv()`](https://hhp94.github.io/slideimp/reference/tf_identity.md)
    and
    [`tf_biwhiten()`](https://hhp94.github.io/slideimp/reference/tf_identity.md)
    transform the data before the fit;
    [`spectra_opts()`](https://hhp94.github.io/slideimp/reference/spectra_opts.md)
    sets the eigensolver;
    [`abb_cols()`](https://hhp94.github.io/slideimp/reference/abb_cols.md)
    proposes columns to draw by approximate Bayesian bootstrap.
  - Results have [`print()`](https://rdrr.io/r/base/print.html),
    [`summary()`](https://rdrr.io/r/base/summary.html),
    [`as.matrix()`](https://rdrr.io/r/base/matrix.html), `[` and
    [`cbind()`](https://rdrr.io/r/base/cbind.html) methods.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  gains `subset`: the result holds only the selected columns, while
  every eligible column still takes part in the fit.

### Changes in results

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  with `post_imp = TRUE` sets a cell that K-NN leaves missing to the
  mean of its column’s observed values. 1.2.0 took a mean that included
  the K-NN imputed values. These cells change value, in
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  and
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  with K-NN too.

- The column means that `post_imp` and `on_infeasible = "mean"` set are
  now computed with base R’s
  [`.colMeans()`](https://rdrr.io/r/base/colSums.html), so they can
  differ from 1.2.0’s in the last bits.

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md):
  every group now reads `obj` as it was passed. 1.2.0’s sequential path
  wrote each group’s imputed values into `obj` before the next group
  ran, so a later group could use them as observed values in its
  auxiliary columns.

- [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  and
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  with PCA: a window or group whose missing target values all lie in
  columns that
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  leaves out is infeasible and goes through `on_infeasible`, as with
  K-NN. 1.2.0 fitted it and set the targets to their column means.

- [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  returns the observed cells of `obj` unchanged. Where three or more
  windows overlapped, 1.2.0 averaged every window’s copy of a column,
  which could move an observed value in its last bit.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  now calls R’s BLAS and LAPACK directly instead of going through
  Armadillo, so its results differ from 1.2.0’s in the last digits under
  the exact solver, and by up to about 1e-8 under LOBPCG.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  on a matrix with no missing values returns a `slideimp_results` object
  like any other call, where 1.2.0 returned the input unchanged.

### Minor improvements and fixes

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  no longer crashes the R session when `ncp` exceeds the rank of the
  completed data, as with replicate samples.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md):
  finite input whose squares overflow stops with “the cross-product of
  the data is not finite; values this large overflow when squared,
  rescale the input”, before either solver. A LAPACK failure reads “the
  symmetric eigensolver failed (LAPACK dsyevr, info = )”. The 1.2.0
  error “`eig_sym_sel` failed to converge” is gone.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  refuses eligible columns holding more than `.Machine$integer.max`
  missing values.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  scans for infinities once per call rather than once per `nb.init`
  initialization. When the input has more columns than rows, only the
  loadings of the columns being imputed are updated inside the EM loop.

- `pin_blas` in
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  and
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  now defaults to `TRUE`, pinning BLAS and OpenMP threads to 1 on mirai
  daemons; a sequential run pins nothing. RhpcBLASctl moves from
  `Suggests` to `Imports`.

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  under mirai sends `obj` to each daemon once, with
  [`mirai::everywhere()`](https://mirai.r-lib.org/reference/everywhere.html),
  so memory grows with the number of daemons. bigmemory is no longer
  imported, and mirai must be `>= 2.5.0`. A new error reports when `obj`
  did not reach every daemon.

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  with `.f = "knn_imp"` or `"slide_imp"` fails an iteration when the
  result is not the width of the repetition’s `subset`. That happens
  when a daemon runs an older slideimp; 1.2.0 read the wrong cells
  without an error.

- [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md),
  and so
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md),
  no longer draws a column whose variance is above zero but below
  `.Machine$double.eps`, the threshold
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  uses for eligibility.

## slideimp 1.2.0

CRAN release: 2026-06-16

### Breaking changes

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  no longer supports Ball Tree-based imputation due to lower accuracy
  compared to brute-force K-NN.

### New features

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  gained a `.progress` argument to display imputation progress.

### Minor improvements and fixes

- Fixed
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  tests that broke due to upstream changes in
  [`missMDA::imputePCA()`](https://rdrr.io/pkg/missMDA/man/imputePCA.html)
  and
  [`FactoMineR::svd.triplet()`](https://rdrr.io/pkg/FactoMineR/man/svd.triplet.html).

- Implemented an internal `.crate()` function, a minimal version of
  `carrier::crate()`.

- `DESCRIPTION` now describes the LOBPCG eigensolver and its warm-start
  behavior across EM iterations in
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).

- Minor wording improvements.

- Added citation information.

## slideimp 1.1.0

CRAN release: 2026-05-01

### Breaking changes

- Minimum R version is now 4.3.0 (was 4.1.0).
- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  no longer caches pair-wise distances. The `use_cache` argument has
  been removed.

### New features

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  gained a `solver` argument to select the LOBPCG eigensolver,
  configured via the new
  [`lobpcg_control()`](https://hhp94.github.io/slideimp/reference/lobpcg_control.md)
  helper.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  gained a `clamp` argument to bound imputed values, e.g. `c(0, 1)` for
  methylation beta values. Also threaded through
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  and
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).

### Minor improvements and fixes

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  now uses [RcppThread](https://github.com/tnagler/RcppThread) instead
  of OpenMP, enabling parallelism on macOS.

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  now infers the subset from `na_loc`, speeding up tuning for
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  and
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md).

- [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  is now an S3 generic, replacing the previous register-on-load pattern
  with `{slideimp.extra}`.

- New article describing how to speed up PCA imputation by choosing
  `solver` and `threshold`.

- Fixed numerical tolerance check failure on CRAN’s ATLAS configuration.

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
  and `mean_imp_col()` have been overhauled to use the faster
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
