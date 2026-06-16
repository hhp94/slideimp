# Grouped K-NN or PCA Imputation

Perform K-NN or PCA imputation independently within feature groups.

## Usage

``` r
group_imp(
  obj,
  group,
  subset = NULL,
  allow_unmapped = FALSE,
  k = NULL,
  ncp = NULL,
  method = NULL,
  cores = 1,
  .progress = TRUE,
  min_group_size = NULL,
  colmax = NULL,
  post_imp = NULL,
  dist_pow = NULL,
  scale = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  row.w = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL,
  solver = NULL,
  lobpcg_control = NULL,
  clamp = NULL,
  pin_blas = FALSE,
  na_check = TRUE,
  on_infeasible = c("error", "skip", "mean")
)
```

## Arguments

- obj:

  A numeric matrix with samples in rows and features in columns.

- group:

  Specification of how features should be grouped for imputation.
  Accepted formats are:

  - A character scalar naming a supported Illumina platform; see Note.

  - A long-format `data.frame` with columns `group` and `feature`.

  - A list-column `data.frame` with a `feature` list-column. Optional
    list-columns are `aux`, for auxiliary feature names, and
    `parameters`, for group-specific parameter lists.

- subset:

  Optional character vector of feature names to impute. If `NULL`, all
  grouped features are imputed. Features in a group but not in `subset`
  are demoted to auxiliary columns for that group. Groups left with zero
  features after demotion are dropped with a message.

- allow_unmapped:

  Logical. If `FALSE`, every column in `colnames(obj)` must appear in
  `group`. If `TRUE`, columns with no group assignment are left
  untouched and are not used as auxiliary columns.

- k:

  Integer or `NULL`. Number of nearest neighbors for K-NN imputation. If
  `NULL`, `k` may be supplied through `group$parameters`.

- ncp:

  Integer or `NULL`. Number of components for PCA imputation. If `NULL`,
  `ncp` may be supplied through `group$parameters`.

- method:

  Character or `NULL`. For K-NN imputation, one of `"euclidean"` or
  `"manhattan"`. For PCA imputation, one of `"regularized"` or `"EM"`.
  If `NULL`, the corresponding backend default is used unless overridden
  by `group$parameters`.

- cores:

  Integer. Number of cores for K-NN imputation only. For PCA imputation,
  use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  to parallelize across groups.

- .progress:

  Logical. If `TRUE`, show progress.

- min_group_size:

  Integer or `NULL`. Minimum total number of columns per group, counting
  both features and auxiliary columns. Groups smaller than this are
  padded with randomly sampled columns from `obj`.

- colmax:

  Numeric scalar between `0` and `1`. Columns with a missing-data
  proportion greater than `colmax` are excluded from the main imputation
  method. Excluded columns are left unchanged unless `post_imp = TRUE`,
  in which case remaining missing values are replaced by column means
  when possible.

- post_imp:

  Logical. If `TRUE`, replace missing values remaining after the main
  imputation method with column means when possible.

- dist_pow:

  Numeric. Power used to penalize more distant neighbors in the weighted
  average. `dist_pow = 0` gives an unweighted average of the nearest
  neighbors.

- scale:

  Logical. If `TRUE`, columns are scaled to unit variance.

- coeff.ridge:

  Numeric. Ridge regularization, used only when
  `method = "regularized"`. Values `< 1` move toward EM PCA; values
  `> 1` move toward mean imputation.

- threshold:

  Numeric. Convergence threshold.

- row.w:

  Row weights, normalized to sum to `1`. `NULL` (equal weights), a
  positive numeric vector of length `nrow(obj)`, or `"n_miss"`
  (down-weight rows with more missing values).

- seed:

  Integer, numeric, or `NULL`. Random seed for reproducibility.

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations.

- miniter:

  Integer. Minimum number of iterations.

- solver:

  Character. Eigensolver: `"auto"` (default), `"exact"`, or `"lobpcg"`.
  `"auto"` runs a short timed probe and picks `"lobpcg"` only when
  clearly faster. Consecutive EM calls warm-start LOBPCG with both the
  previous eigenblock and search direction. When `nb.init > 1`, the auto
  choice from the first init is reused. See Performance tips.

- lobpcg_control:

  A list of LOBPCG eigensolver control options, usually created by
  [`lobpcg_control()`](https://hhp94.github.io/slideimp/reference/lobpcg_control.md).
  A plain named list is also accepted. Ignored when `solver = "exact"`.

- clamp:

  Optional numeric vector `c(lower, upper)` bounding PCA-imputed values
  (use `-Inf`/`Inf` for one-sided, `NULL` for none). E.g., `c(0, 1)` for
  DNAm beta values. Observed values are not clamped.

- pin_blas:

  Logical. If `TRUE`, pin BLAS threads to 1 to reduce contention when
  using parallel PCA on systems linked with multithreaded BLAS.

- na_check:

  Logical. If `TRUE`, check whether the returned matrix still contains
  missing values.

- on_infeasible:

  Character. One of `"error"`, `"skip"`, or `"mean"`. Controls behavior
  when a group is infeasible for imputation, for example when `k` or
  `ncp` exceeds the number of usable columns after applying `colmax`.

## Value

A numeric matrix of the same dimensions as `obj`, with missing values
imputed. The returned object has class `slideimp_results`.

## Details

`group_imp()` performs K-NN or PCA imputation on feature groups
independently, which can substantially reduce runtime for large
matrices.

Specify `k` and related arguments to use K-NN imputation, or `ncp` and
related arguments to use PCA imputation. If both `k` and `ncp` are
`NULL`, `group$parameters` must supply either `k` or `ncp` for every
group.

Group-specific parameters in `group$parameters` take priority over
global arguments. Global arguments fill in any missing values. All
groups must use the same imputation method.

For method-specific arguments inherited from
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md) or
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md),
`NULL` means the backend default is used unless overridden by
`group$parameters`.

Per-group `k` is capped at the number of usable columns in the group
minus one. Per-group `ncp` is capped at the maximum feasible number of
PCA components for that group's submatrix. A warning is issued when
capping occurs.

## Note

A character scalar can be passed to `group` to name a supported Illumina
platform, such as `"EPICv2"` or `"EPICv2_deduped"`. This requires the
optional `slideimp.extra` package to be installed. Supported platforms
are listed in the `slideimp_arrays` object in the `slideimp.extra`
package.

## Parallelization

- K-NN: use the `cores` argument. If `mirai` daemons are active, `cores`
  is automatically set to `1` to avoid nested parallelism.

- PCA: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

When running PCA imputation in parallel with `mirai`, set
`pin_blas = TRUE` in
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
or `group_imp()` to prevent BLAS threads from oversubscribing CPU cores.
This relies on `RhpcBLASctl` and works with OpenBLAS and MKL (typical on
Linux, and on Windows after an OpenBLAS swap). `pin_blas = TRUE` may
have no effect on macOS.

## PCA Performance tips

Speed comes from three levers: `solver` (through LOBPCG with
warm-start), `threshold`, and `scale`. Tune these first, then accuracy
parameters (`ncp`, `coeff.ridge`) on a representative subset.

**Exact vs. LOBPCG with warm-start.** Whether `"lobpcg"` beats `"exact"`
depends on size and low-rankness: `"lobpcg"` is preferred for large,
approximately low-rank matrices with small `ncp`, and `"exact"` for
small matrices (including
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
windows), where it is faster and more robust. Separately, the warm-start
makes each successive solve cheap:
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
warm-starts LOBPCG with the previous eigenblock and search direction, so
once imputed values stabilize, later solves converge in a few
iterations. The payoff therefore grows with the number of EM iterations,
independent of low-rankness. `solver = "auto"` (default) probes both and
is a safe start.

**Threshold.** The default `1e-6` is conservative; `1e-5` is often
faster with very similar values.

**Scale.** For columns on a common scale (e.g., DNAm beta values in
`[0, 1]`), `scale = FALSE` can be faster and more accurate.

**Parallel and BLAS.** In parallel via
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
or `group_imp()` with a multithreaded BLAS, set `pin_blas = TRUE` to
avoid thread oversubscription. On Windows, the stock BLAS can be slow.
Advanced users can swap in
[OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).

See [Speeding up PCA
imputation](https://hhp94.github.io/slideimp/articles/speeding-up-pca-imputation.html)
for the full workflow.

## See also

[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md),
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)

## Examples

``` r
set.seed(1234)
to_test <- sim_mat(10, 20, perc_total_na = 0.05, perc_col_na = 1)
obj <- to_test$input
group <- to_test$col_group
head(group)
#>    feature  group
#> 1 feature1 group2
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group2
#> 5 feature5 group1
#> 6 feature6 group1

# Simple grouped K-NN imputation
results <- group_imp(obj, group = group, k = 2, .progress = FALSE)
#> Imputing 2 groups using KNN.
#> Running mode: sequential
results
#> Method: group_imp (KNN imputation)
#> Dimensions: 10 x 20
#> 
#>          feature1  feature2  feature3  feature4   feature5  feature6
#> sample1 0.1568098 0.3213953 0.2768746 1.0000000 0.07839046 0.4219375
#> sample2 0.4098416 0.8658918 0.8221066 0.6396885 0.68926345 1.0000000
#> sample3 0.6801685 1.0000000 1.0000000 0.9953450 0.75246030 0.6958550
#> sample4 0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.0000000
#> sample5 0.9639671 0.6409137 0.5111760 0.7184678 0.81815288 0.5883255
#> sample6 0.7031741 0.3731062 0.6782811 0.7542872 0.99910554 0.9070005
#> # Showing 6 of 10 rows and 6 of 20 columns

# Impute only a subset of features
subset_features <- sample(to_test$col_group$feature, size = 10)
knn_subset <- group_imp(
  obj,
  group = group,
  subset = subset_features,
  k = 2,
  .progress = FALSE
)
#> Imputing 2 groups using KNN.
#> Running mode: sequential

# Use prep_groups() to inspect and edit per-group parameters
prepped <- prep_groups(colnames(obj), group)
prepped$parameters <- lapply(seq_len(nrow(prepped)), function(i) list(k = 2))
prepped$parameters[[2]]$k <- 4
knn_grouped <- group_imp(obj, group = prepped, .progress = FALSE)
#> Imputing 2 groups using KNN.
#> Running mode: sequential

if (FALSE) { # interactive() && requireNamespace("mirai", quietly = TRUE)
# PCA imputation with mirai parallelism
mirai::daemons(2)
pca_grouped <- group_imp(obj, group = group, ncp = 2)
mirai::daemons(0)
pca_grouped
}
```
