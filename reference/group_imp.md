# Grouped K-NN or PCA Imputation

Perform K-NN or PCA imputation independently on feature groups (e.g. by
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
  pin_blas = FALSE
)
```

## Arguments

- obj:

  A numeric matrix with **samples in rows** and **features in columns**.

- group:

  Specification of how features should be grouped for imputation.
  Accepts three formats:

  - `character`: (Requires the `{slideimp.extra}` package, available on
    GitHub). A single string naming a supported Illumina platform (e.g.,
    `"EPICv2"`, `"EPICv2_deduped"`). The manifest is fetched
    automatically. For installation instructions, see the package
    README. Supported platforms can be viewed via
    [`slideimp.extra::slideimp_arrays`](https://rdrr.io/pkg/slideimp.extra/man/slideimp_arrays.html).

  - `data.frame` (Long format):

    - `group`: Column identifying the group for each feature.

    - `feature`: Character column of individual feature names.

  - `data.frame` (List-column format):

    - `feature`: List-column of character vectors to impute (e.g., from
      [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)).

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

  Logical. If `FALSE`, every column in `colnames(obj)` *must* appear in
  `group`. If `TRUE`, columns that have no group assignment are left
  untouched (neither imputed nor used as auxiliary columns) and a
  message is issued instead of an error.

- k:

  Number of nearest neighbors for imputation. 10 is a good starting
  point.

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

  Show imputation progress (default = `TRUE`).

- min_group_size:

  Integer or `NULL`. Minimum number of columns (features + aux) per
  group. Groups smaller than this are padded with randomly sampled
  columns from `obj`. Passed to
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  internally.

- colmax:

  A number from 0 to 1. Threshold of column-wise missing data rate above
  which K-NN imputation is skipped.

- post_imp:

  Whether to impute remaining missing values (those that failed K-NN
  imputation) using column means (default = `TRUE`).

- dist_pow:

  The amount of penalization for further away nearest neighbors in the
  weighted average. `dist_pow = 0` (default) is the simple average of
  the nearest neighbors.

- tree:

  Logical. `FALSE` (default) = brute-force K-NN. `TRUE` = use `{mlpack}`
  BallTree.

- max_cache:

  Maximum allowed cache size in GB (default `4`). When greater than `0`,
  pairwise distances between columns with missing values are
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

  Numeric or `NULL`. Random seed for reproducibility when padding for
  `min_group_size` and passed to
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations for the algorithm.

- miniter:

  Integer. Minimum number of iterations for the algorithm.

- pin_blas:

  Logical. If `TRUE`, pin BLAS thread to 1 to help with parallel
  performance on systems linked with multi-threaded BLAS.

## Value

A numeric matrix of the same dimensions as `obj` with missing values
imputed.

## Details

This function performs K-NN or PCA imputation on groups of features
independently, which significantly reduces imputation time for large
datasets.

Specify `k` and related arguments to use K-NN, `ncp` and related
arguments for PCA imputation. If `k` and `ncp` are both `NULL`, then
`group$parameters` must contain either `k` or `ncp` for every group.

### Parameter resolution

Group-wise parameters (in `group$parameters`) take priority. Global
arguments (`k`, `ncp`, `method`, etc.) fill in any gaps where a group
has no value set. All groups must agree on the imputation method (all
KNN or all PCA). Per-group `k` is capped at `group_size - 1` and `ncp`
at `min(nrow(group) - 2L, ncol(group) - 1L)`, with a warning when
capping occurs.

### Strategies for grouping

- Breaking down search space by chromosomes

- Grouping features with their flanking values/neighbors

- Using clusters identified by column clustering techniques

## Parallelization

Parallelization behavior depends on the imputation method:

- **KNN**: use the `cores` argument (if OpenMP is available). If
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  are also active, `cores` is automatically set to 1 to avoid nested
  parallelism.

- **PCA**: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

**Linux / OpenBLAS / MKL users:** If your machine uses a multi-threaded
BLAS (e.g., OpenBLAS or Intel MKL), set `pin_blas = TRUE` when tuning
PCA imputation in parallel. Without it, BLAS threads and `mirai` workers
compete for cores, which can cause slowdowns (CPU thrashing).

**macOS users:** OpenMP is typically unavailable on macOS unless
manually configured. `cores` will fall back to 1 automatically; use
[`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html) for
parallelization instead.

## See also

[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)

## Examples

``` r
# Generate example data with missing values
set.seed(1234)
to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
obj <- to_test$input
group <- to_test$col_group # metadata that maps `colnames(obj)` to groups
head(group)
#>    feature  group
#> 1 feature1 group2
#> 2 feature2 group2
#> 3 feature3 group2
#> 4 feature4 group2
#> 5 feature5 group2
#> 6 feature6 group2

# Simple grouped K-NN imputation
results <- group_imp(obj, group = group, k = 2)
#> Imputing 2 group(s) using KNN.
#> Running Mode: sequential ...

# Impute only a subset of features
subset_features <- sample(to_test$col_group$feature, size = 10)
knn_subset <- group_imp(obj, group = group, subset = subset_features, k = 5)
#> Imputing 2 group(s) using KNN.
#> Running Mode: sequential ...

# Use prep_groups() to inspect and tweak per-group parameters
prepped <- prep_groups(colnames(obj), group)
prepped$parameters <- lapply(seq_len(nrow(prepped)), \(i) list(k = 5))
prepped$parameters[[2]]$k <- 10
knn_grouped <- group_imp(obj, group = prepped, cores = 2)
#> Imputing 2 group(s) using KNN.
#> Running Mode: parallel (OpenMP within groups)...

# PCA imputation with mirai parallelism
mirai::daemons(2)
pca_grouped <- group_imp(obj, group = group, ncp = 2)
#> Imputing 2 group(s) using PCA.
#> Running Mode: parallel (mirai across groups)...
#> Tip: set `pin_blas = TRUE` may improve parallel performance.
mirai::daemons(0)
pca_grouped
#> slideimp_results (PCA)
#> Dimensions: 50 x 20
#> 
#>           feature1   feature2  feature3  feature4  feature5  feature6
#> sample1 0.05685662 0.39194101 0.2995808 0.5266989 0.3693625 0.3372606
#> sample2 0.49992298 0.54447261 0.5504182 0.6699885 0.7143775 0.7084667
#> sample3 0.54575968 0.78179817 0.8510282 0.7503179 0.7469653 0.8390176
#> sample4 0.00000000 0.07766851 0.2063256 0.2967793 0.2200145 0.2558501
#> sample5 0.59340652 0.50932703 0.5127284 0.6909102 0.5803437 0.7547226
#> sample6 0.62356150 0.64673748 0.5254893 0.7499320 0.6118071 0.6436882
#> 
#> # Showing [1:6, 1:6] of full matrix
```
