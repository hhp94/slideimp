
<!-- README.md is generated from README.Rmd. Please edit that file -->

# slideimp

<!-- badges: start -->

[![R-CMD-check](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`{slideimp}` is a lightweight R package for fast K-NN and PCA imputation
of missing values in high-dimensional numeric matrices.

**Core functions**

- `group_imp()`: Parallelizable group-wise (e.g., by chromosomes or
  column clusters) K-NN or PCA imputation with optional auxiliary
  features and group-wise parameters.
  - (Basic): `group_imp()` requires a `data.frame` for the `group`
    argument. This `data.frame` must contains a `feature` column
    (character strings of feature names) and a `group` column (defining
    the group).
  - (Advanced): `prep_groups()`: `group_imp()`’s helper function to
    create groups based on a mapping data.frame (i.e., Illumina
    manifests). See
    [`{slideimp.extra}`](https://github.com/hhp94/slideimp.extra) on
    GitHub for tools to process common Illumina manifests.
- `slide_imp()`: Sliding window K-NN or PCA imputation for extremely
  high-dimensional numeric matrices (WGBS/EM-seq) with ordered features
  (i.e., by genomic position).
- `tune_imp()`: Parallelizable hyperparameter tuning with repeated
  cross-validation; works with built-in or custom imputation functions.
- `knn_imp()`: Full-matrix K-NN imputation with multi-core
  parallelization, [`{mlpack}`](https://mlpack.org/) KD/Ball-Tree
  nearest neighbor implementation (for data with very low missing rates
  and extremely high dimensions), and optional subset imputation (ideal
  for epigenetic clock calculations).
- `pca_imp()`: Optimized version of
  [`missMDA::imputePCA()`](http://factominer.free.fr/missMDA/PCA.html)
  for high-dimensional numeric matrices.

## Installation

The stable version of `{slideimp}` can be installed from CRAN using:

``` r
install.packages("slideimp")
```

You can install the development version of `{slideimp}` with:

``` r
pak::pkg_install("hhp94/slideimp")
```

## Workflow

Let’s simulate some DNA methylation (DNAm) microarray data from 2
chromosomes. All `{slideimp}` functions expect the input to be a numeric
matrix where variables are stored in the columns.

``` r
library(slideimp)
# Simulate data from 2 chromosomes. 20 samples, 50 CpGs
set.seed(1234)
sim_obj <- sim_mat(n = 20, p = 50, perc_total_na = 0.3, perc_col_na = 1, n_col_groups = 2)

sim_obj
#> $col_group (2 column groups)
#>    feature  group
#> 1 feature1 group2
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group2
#> 5 feature5 group1
#> 6 feature6 group1
#> 
#> $row_group (1 row groups)
#>    sample  group
#> 1 sample1 group1
#> 2 sample2 group1
#> 3 sample3 group1
#> 4 sample4 group1
#> 5 sample5 group1
#> 6 sample6 group1
#> 
#> $input (20 x 50)
#>          feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.3486482 0.7385414 0.4077444        NA        NA 0.2434143
#> sample2        NA        NA        NA 0.3722729 0.5061132 0.3923603
#> sample3 0.7185848 0.7351035        NA        NA 0.7634236 1.0000000
#> sample4        NA 0.0000000 0.0000000 0.0000000        NA 0.2106358
#> sample5 0.5228314 0.5306182 0.5685354 0.5383513 0.4680080 0.8518388
#> sample6 0.3768380        NA        NA 0.5276245 0.6722794 0.5740639
#> # Showing [1:6, 1:6] of full matrix

# Extract the $input attribute
obj <- sim_obj$input
```

We can estimate the prediction accuracy of different methods and tune
hyperparameters prior to imputation with `tune_imp()`. For example tune
the results using 2 repeats (`rep = 2`) for illustration (increase in
actual analyses).

``` r
# Our parameters grid contains 2 values: k = 5 and 20.
knn_params <- data.frame(k = c(5, 20))
# Parallelization with OpenMP is controlled by `cores` only for knn or slideimp knn with OpenMP only.
tune_knn <- tune_imp(obj, parameters = knn_params, .f = "knn_imp", cores = 2, rep = 2)
#> Tuning knn_imp
#> Step 1/2: Injecting NA
#> Running Mode: parallel...
#> Step 2/2: Tuning
compute_metrics(tune_knn)
#>        k param_set rep error   n n_miss .metric .estimator .estimate
#> mae    5         1   1  <NA> 100      0     mae   standard 0.1670189
#> rmse   5         1   1  <NA> 100      0    rmse   standard 0.2160274
#> mae1  20         2   1  <NA> 100      0     mae   standard 0.1648807
#> rmse1 20         2   1  <NA> 100      0    rmse   standard 0.1998409
#> mae2   5         1   2  <NA> 100      0     mae   standard 0.1928833
#> rmse2  5         1   2  <NA> 100      0    rmse   standard 0.2538705
#> mae3  20         2   2  <NA> 100      0     mae   standard 0.1595128
#> rmse3 20         2   2  <NA> 100      0    rmse   standard 0.2125974
```

For K-NN without OpenMP, PCA, and custom functions (see vignette), setup
parallelization with `mirai::daemons()`.

``` r
mirai::daemons(2) # 2 Cores
# Note, for PCA and custom functions, cores is controlled by `mirai::daemons()`
# and the `cores` argument is ignored.

# PCA imputation.
pca_params <- data.frame(ncp = c(1, 5))
tune_pca <- tune_imp(obj, parameters = pca_params, .f = "pca_imp", rep = 2)

mirai::daemons(0) # Close daemons
```

Then, preferably perform imputation by group with `group_imp()` if the
variables can be meaningfully grouped (e.g., by chromosomes).

- `group_imp()` allows imputation to be performed separately within
  defined groups (e.g., by chromosome), which significantly reduces
  runtime and can increase accuracy for both K-NN and PCA imputation.
  For the basic API, use a data.frame that contain a `feature` column
  (character strings of feature names) and a `group` column (defining
  the group).
- PCA-based imputation with `group_imp()` can be parallelized using the
  `{mirai}` package, similar to how parallelization is done with
  `tune_imp()`.

``` r
# We picked k = 5 (use parameters with lowest errors in actual analyses)
knn_group_results <- group_imp(obj, group = sim_obj$col_group, k = 5, cores = 2)

# Similar to `tune_imp`, parallelization is controlled by `mirai::daemons()`
mirai::daemons(2)
knn_group_results <- group_imp(obj, group = sim_obj$col_group, ncp = 3)
mirai::daemons(0)
```

Alternatively, full matrix imputation can be performed using `knn_imp()`
or `pca_imp()`.

``` r
full_knn_results <- knn_imp(obj = obj, k = 5)
full_pca_results <- pca_imp(obj = obj, ncp = 5)
```

## Sliding Window Imputation

Sliding window imputation can be performed using `slide_imp()`.
**Note:** DNAm WGBS/EM-seq data should be grouped by chromosomes
(i.e. run `slide_imp()` separately on each chromosome) before
imputation. See the package vignette for more details. This function is
**not** for Illumina DNAm microarrays, see `group_imp()`.

``` r
# Simulate some data
chr1_beta <- sim_mat(n = 10, p = 2000, perc_total_na = 0.3, perc_col_na = 1, n_col_groups = 1)$input
dim(chr1_beta)
#> [1]   10 2000
chr1_beta[1:5, 1:5]
#>           feature1  feature2   feature3  feature4  feature5
#> sample1 0.85350514 0.9300089 0.53948768 0.1757192 0.2016609
#> sample2         NA        NA         NA 0.7618257 0.6474535
#> sample3         NA        NA 0.03834917        NA 0.3475937
#> sample4 0.03859717 0.2061508 0.50760195        NA        NA
#> sample5 0.49341163        NA 0.61117531 0.8327905 1.0000000
```

`slide_imp()` parameters:

- `location`: **required** - a sorted numeric vector of length
  `ncol(obj)` giving the position of each column (e.g. genomic
  coordinates in *bp*).
- `window_size`: **required** - width of each sliding window (same unit
  as `location`).
- `overlap_size`: **optional** - overlap width between consecutive
  windows (same units as `location`). Must be strictly less than
  `window_size`.
- `min_window_n`: **required** - minimum number of columns a window must
  contain to be imputed. Windows smaller than this are dropped. Must be
  greater than `k` (for KNN) or `ncp` (for PCA).
- `dry_run`: **optional** - return just the calculated windows to
  examine which windows are included.
- `k`: **required** - *(specifying KNN imputation)* Number of nearest
  neighbors to use inside each window.
- `ncp`: **required** - *(specifying PCA imputation)* Number of
  principal components to retain. Use this instead of `k` when
  performing sliding-window PCA imputation.
- `subset`: **optional** - impute just a subset of features (i.e., just
  clock CpGs).
- `flank`: **optional** - to build flanking windows of `window_size`
  around features provided in `subset`.

First, let’s perform a dry run to examine the windows that will be
imputed by `slide_imp`.

``` r
location <- seq_len(ncol(chr1_beta)) # 1, 2, ..., 2000 for this simulated chromosome

slide_imp(
  obj = chr1_beta,
  location = location,
  window_size = 50,
  overlap_size = 5,
  min_window_n = 11, # must be > k = 10
  dry_run = TRUE
)
#> # slideimp table: 45 x 4
#>  start end window_n  subset_local
#>      1  50       50 <double [50]>
#>     46  95       50 <double [50]>
#>     91 140       50 <double [50]>
#>    136 185       50 <double [50]>
#>    181 230       50 <double [50]>
#>    226 275       50 <double [50]>
#>    271 320       50 <double [50]>
#>    316 365       50 <double [50]>
#>    361 410       50 <double [50]>
#>    406 455       50 <double [50]>
#> # ... with 35 more rows
```

Then, we can perform the imputation using the given parameters.

``` r
slide_imp(
  obj = chr1_beta,
  location = location,
  window_size = 50,
  overlap_size = 5,
  min_window_n = 11, # must be > k = 10
  k = 10,
  cores = 2,
  .progress = FALSE
)
#> slideimp_results (KNN)
#> Dimensions: 10 x 2000
#> 
#>           feature1   feature2   feature3  feature4  feature5
#> sample1 0.85350514 0.93000887 0.53948768 0.1757192 0.2016609
#> sample2 0.60408183 0.77061812 0.73823976 0.7618257 0.6474535
#> sample3 0.34882700 0.02388497 0.03834917 0.2039694 0.3475937
#> sample4 0.03859717 0.20615085 0.50760195 0.4065345 0.3499150
#> sample5 0.49341163 0.66639278 0.61117531 0.8327905 1.0000000
#> 
#> # Showing [1:5, 1:5] of full matrix
```
