# Grouped K-NN or PCA Imputation

Perform K-NN or PCA imputation independently on feature groups (e.g., by
chromosomes, flanking probes, or clustering-based groups).

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
  tree = NULL,
  max_cache = NULL,
  scale = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  row.w = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL,
  pin_blas = FALSE,
  na_check = TRUE,
  on_infeasible = c("error", "skip", "mean")
)
```

## Arguments

- obj:

  A numeric matrix with **samples in rows** and **features in columns**.

- group:

  Specification of how features should be grouped for imputation.
  Accepts three formats:

  - `character`: string naming a supported Illumina platform; see the
    Note section.

  - `data.frame` (Long format):

    - `group`: Column identifying the group for each feature.

    - `feature`: Character column of individual feature names.

  - `data.frame` (List-column format):

    - `feature`: List-column of character vectors to impute. A row is a
      group.

    - `aux`: (Optional) List-column of auxiliary names used for context.

    - `parameters`: (Optional) List-column of group-specific parameter
      lists.

- subset:

  Character vector of feature names to impute (default `NULL` means
  impute all features). Must be a subset of `obj_cn` (`colnames(obj)`)
  and must appear in at least one group's `feature`. Features in a group
  but not in `subset` are demoted to auxiliary columns for that group.
  Groups left with zero features after demotion are dropped with a
  message.

- allow_unmapped:

  Logical. If `FALSE`, every column in `colnames(obj)` must appear in
  `group`. If `TRUE`, columns with no group assignment are left
  untouched (neither imputed nor used as auxiliary columns) and a
  message is issued instead of an error.

- k:

  Integer. Number of nearest neighbors for imputation. 10 is a good
  starting point.

- ncp:

  Integer. Number of components used to predict the missing entries.

- method:

  For K-NN imputation: distance metric to use (`"euclidean"` or
  `"manhattan"`). For PCA imputation: regularization imputation
  algorithm (`"regularized"` or `"EM"`).

- cores:

  The number of OpenMP cores for K-NN imputation **only**. For PCA or
  mirai-based parallelism, use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead.

- .progress:

  Show imputation progress (default `TRUE`).

- min_group_size:

  Integer or `NULL`. Minimum column count (features + aux) per group.
  Groups smaller than this are padded with randomly sampled columns from
  `obj`. Passed to
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  internally.

- colmax:

  Numeric. A number from 0 to 1. Threshold of column-wise missing data
  rate above which imputation is skipped.

- post_imp:

  Boolean. Whether to impute remaining missing values (those that failed
  imputation) using column means.

- dist_pow:

  Numeric. The amount of penalization for further away nearest neighbors
  in the weighted average. `dist_pow = 0` (default) is the simple
  average of the nearest neighbors.

- tree:

  Logical. `FALSE` (default) uses brute-force K-NN. `TRUE` uses `mlpack`
  BallTree.

- max_cache:

  Numeric. Maximum allowed cache size in GB (default `4`). When greater
  than `0`, pairwise distances between columns with missing values are
  pre-computed and cached, which is faster for moderate-sized data but
  uses O(m^2) memory where m is the number of columns with missing
  values. Set to `0` to disable caching and trade speed for lower memory
  usage.

- scale:

  Logical. If `TRUE` (default), variables are scaled to have unit
  variance.

- coeff.ridge:

  Numeric. Ridge regularization coefficient (default is 1). Only used if
  `method = "regularized"`. Values \< 1 regularize less (closer to EM);
  values \> 1 regularize more (closer to mean imputation).

- threshold:

  Numeric. The threshold for assessing convergence.

- row.w:

  Row weights (internally normalized to sum to 1). Can be one of:

  - `NULL` (default): All rows weighted equally.

  - A numeric vector: Custom positive weights of length `nrow(obj)`.

  - `"n_miss"`: Rows with more missing values receive lower weight.

- seed:

  Numeric or `NULL`. Random seed for reproducibility.

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations for the algorithm.

- miniter:

  Integer. Minimum number of iterations for the algorithm.

- pin_blas:

  Logical. If `TRUE`, pin BLAS threads to 1 to reduce contention when
  using parallel PCA on systems linked with multi-threaded BLAS.

- na_check:

  Boolean. Check for leftover `NA` values in the results or not
  (internal use).

- on_infeasible:

  Character, one of `"error"` (default on `group_imp()`), `"skip"`, or
  `"mean"` (default on
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)).
  Controls behaviour when a group is infeasible for imputation, e.g.,
  `k`/`ncp` exceeds the number of usable columns after applying
  `colmax`, or all subset columns in the group exceed `colmax`.

## Value

A numeric matrix of the same dimensions as `obj` with missing values
imputed.

## Details

Performs K-NN or PCA imputation on groups of features independently,
which significantly reduces imputation time for large datasets.

Specify `k` and related arguments to use K-NN, or `ncp` and related
arguments for PCA imputation. If both `k` and `ncp` are `NULL`,
`group$parameters` must supply either `k` or `ncp` for every group.

### Parameter resolution

Group-wise parameters (in `group$parameters`) take priority; global
arguments (`k`, `ncp`, `method`, etc.) fill in any gaps. All groups must
use the same imputation method. Per-group `k` is capped at
`group_size - 1` and `ncp` at `min(nrow(group) - 2L, ncol(group) - 1L)`,
with a warning when capping occurs.

### Grouping strategies

- Chromosomal grouping to break down the search space.

- Flanking-probe groups for spatially local imputation.

- Column-clustering to form correlation-based groups.

## Note

A `character` string can be passed to `group` to name a supported
Illumina platform (e.g., `"EPICv2"`, `"EPICv2_deduped"`), which fetches
the manifest automatically. This requires the `slideimp.extra` package
(available on GitHub; see its README for installation instructions).
Supported platforms are listed in `slideimp.extra::slideimp_arrays`.

## Parallelization

- **K-NN**: use the `cores` argument (requires OpenMP). If
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  are active, `cores` is automatically set to 1 to avoid nested
  parallelism.

- **PCA**: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

On macOS, OpenMP is typically unavailable and `cores` falls back to

1.  Use
    [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
    for parallelization instead.

On Linux with OpenBLAS or MKL, set `pin_blas = TRUE` when running
parallel PCA to prevent BLAS threads and `mirai` workers competing for
cores.

## See also

[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)

## Examples

``` r
# Generate example data with missing values
set.seed(1234)
to_test <- sim_mat(10, 20, perc_total_na = 0.05, perc_col_na = 1)
obj <- to_test$input
group <- to_test$col_group # metadata that maps `colnames(obj)` to groups
head(group)
#>    feature  group
#> 1 feature1 group2
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group2
#> 5 feature5 group1
#> 6 feature6 group1

# Simple grouped K-NN imputation
results <- group_imp(obj, group = group, k = 2)
#> Imputing 2 group(s) using KNN.
#> Running Mode: sequential ...

# Impute only a subset of features
subset_features <- sample(to_test$col_group$feature, size = 10)
knn_subset <- group_imp(obj, group = group, subset = subset_features, k = 2)
#> Imputing 2 group(s) using KNN.
#> Running Mode: sequential ...

# Use prep_groups() to inspect and tweak per-group parameters
prepped <- prep_groups(colnames(obj), group)
prepped$parameters <- lapply(seq_len(nrow(prepped)), \(i) list(k = 2))
prepped$parameters[[2]]$k <- 4
knn_grouped <- group_imp(obj, group = prepped, cores = 2)
#> Imputing 2 group(s) using KNN.
#> Running Mode: parallel (OpenMP within groups)...
if (FALSE) { # interactive() && requireNamespace("mirai", quietly = TRUE)
# PCA imputation with mirai parallelism
mirai::daemons(2)
pca_grouped <- group_imp(obj, group = group, ncp = 2)
mirai::daemons(0)
pca_grouped
}
```
