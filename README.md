
<!-- README.md is generated from README.Rmd. Please edit that file -->

# slideimp

<!-- badges: start -->

[![R-CMD-check](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`{slideimp}` is a lightweight R package for fast K-NN and PCA imputation
of missing values in high-dimensional numeric matrices.

## Installation

The stable version of `{slideimp}` can be installed from CRAN using:

``` r
install.packages("slideimp")
```

You can install the development version of `{slideimp}` with:

``` r
pak::pkg_install("hhp94/slideimp")
```

You can install the optional
[`{slideimp.extra}`](https://github.com/hhp94/slideimp.extra) package
(which provides lightweight Illumina manifests) with:

``` r
pak::pkg_install("hhp94/slideimp.extra")
```

## PCA or KNN imputation of Illumina DNAm microarrays

- This is a simulated beta matrix from the
  [MSA](https://www.illumina.com/products/by-type/microarray-kits/infinium-methylation-screening-array.html)
  microarray. CpGs are in columns.

``` r
dim(MSA_beta_matrix)
# [1]     20 281797
MSA_beta_matrix[1:4, 1:4]
#         cg06185909_TC11 cg18975462_BC11 cg20516119_TC11 cg10149399_BC11
# sample1              NA       0.5023435       0.3835431              NA
# sample2       0.4907466       0.5095459       0.9025816       0.4313347
# sample3       0.6885036              NA       0.7646753       0.4498772
# sample4       0.0000000       0.0000000              NA       0.0000000
```

- Chromosome-wise imputation of Illumina microarrays can be performed
  with a single function call (requires the Illumina manifests,
  conveniently provided by the
  [`slideimp.extra`](https://github.com/hhp94/slideimp.extra) package).

- This example demonstrates PCA imputation. To use KNN imputation
  instead, supply the `k` argument.

``` r
library(slideimp)
library(slideimp.extra)
library(mirai)

imputed <- group_imp(
  obj = MSA_beta_matrix,
  group = "MSA", # <- this feature requires the `slideimp.extra` package
  ncp = 10, # <- change to `k` for KNN imputation
  .progress = FALSE # <- turn on to monitor progress of longer running jobs
)
# Found cleaned manifest for 'MSA'
# Groups 24 dropped: no features remaining after matching obj columns.
# Imputing 25 group(s) using PCA.
# Running Mode: sequential ...

print(imputed, n = 4, p = 4)
# slideimp_results (PCA)
# Dimensions: 20 x 281797
#
#         cg06185909_TC11 cg18975462_BC11 cg20516119_TC11 cg10149399_BC11
# sample1       0.1517542       0.5023435      0.38354308       0.2067731
# sample2       0.4907466       0.5095459      0.90258164       0.4313347
# sample3       0.6885036       0.7339375      0.76467530       0.4498772
# sample4       0.0000000       0.0000000      0.05230101       0.0000000
#
# # Showing [1:4, 1:4] of full matrix
```

- Tips
  - Remove the “ctl” probes (control probes) before imputation with
    `obj <- obj[, !grepl("^ctl", colnames(obj))]` or set
    `allow_unmapped = TRUE` to ignore them.
  - Deduplicate the MSA and EPICv2 data after imputation with
    `slideimp.extra::dedup_matrix()`.
  - Speed up KNN imputation with OpenMP by using the `cores` argument
    instead of `{mirai}` (available by default on Windows and Linux). If
    you only need clock CpGs, provide the `subset` argument to ignore
    all other probes.

## Extended Workflow

- We simulate DNA methylation (DNAm) microarray data from 2 chromosomes.
  All `{slideimp}` functions expect the input to be a numeric matrix
  where variables are stored in the columns.

``` r
library(slideimp)
set.seed(1234)
sim_obj <- sim_mat(n = 20, p = 100, n_col_groups = 2)

# Extract the $input attribute
obj <- sim_obj$input

# Extract metadata to map the features to groups
group_df <- sim_obj$col_group
```

- Hyperparameters are tuned using `tune_imp()`. We evaluate the
  following options with grid search:

  - Number of nearest neighbors (`k`): 5 or 20
  - Inverse distance power (`dist_pow`) for weighted averaging: 1 or 2
    (higher values assign lower weights to more distant neighbors)

- Tuning is performed on a subset of the data. We use 10 repeats
  (`n_reps = 10`) of cross-validation for evaluation. We re-impute 50
  observed values (`num_na = 50`) to compute the rmse and mae. Increase
  both `n_reps` and `num_na` in real analyses to increase reliability of
  error estimates.

- **Note:** Parallelization via the `cores` argument is only available
  for KNN imputation with OpenMP.

``` r
knn_params <- expand.grid(k = c(5, 20), dist_pow = c(1, 2))
group2_columns <- subset(group_df, group == "group2")
group2_only <- obj[, group2_columns$feature]

tune_knn <- tune_imp(
  group2_only,
  parameters = knn_params,
  .f = "knn_imp",
  cores = 8,
  n_reps = 10,
  num_na = 50
)
#> Tuning knn_imp
#> Step 1/2: Resolving NA locations
#> Running Mode: parallel...
#> Step 2/2: Tuning
```

- Calculate errors using `compute_metrics()` or
  [`{yardstick}`](https://github.com/tidymodels/yardstick) functions.

``` r
metrics <- compute_metrics(tune_knn)

# equivalently: dplyr::summarize(metrics, ..., .by = c(k, dist_pow, .metric))
sum_metrics <- do.call(
  data.frame,
  aggregate(
    .estimate ~ k + dist_pow + .metric,
    data = metrics,
    FUN = function(x) {
      c(
        n = length(x),
        mean_error = mean(x),
        sd_error = sd(x)
      )
    }
  )
)

sum_metrics[order(sum_metrics$.estimate.mean_error), ]
#>    k dist_pow .metric .estimate.n .estimate.mean_error .estimate.sd_error
#> 2 20        1     mae          10            0.1501427         0.01599598
#> 4 20        2     mae          10            0.1519464         0.01638402
#> 1  5        1     mae          10            0.1656087         0.01683512
#> 3  5        2     mae          10            0.1669792         0.01666404
#> 6 20        1    rmse          10            0.1897722         0.01915653
#> 8 20        2    rmse          10            0.1918765         0.01936392
#> 5  5        1    rmse          10            0.2081833         0.01949624
#> 7  5        2    rmse          10            0.2101955         0.01884442
```

- For K-NN without OpenMP, PCA imputation, and custom functions (see
  vignette), set up parallelization with `mirai::daemons()`.
- **Note:** For machines with multi-threaded BLAS, turn on
  `pin_blas = TRUE` when tuning PCA imputation in parallel to avoid
  thrashing.

``` r
mirai::daemons(2) # 2 Cores

# PCA imputation.
pca_params <- data.frame(ncp = c(1, 5))
# For machines with multi-threaded BLAS, turn on `pin_blas = TRUE`
tune_pca <- tune_imp(obj, parameters = pca_params, .f = "pca_imp", n_reps = 10, num_na = 50)

mirai::daemons(0) # Close daemons
```

- Then perform imputation with `group_imp()` using the best parameters.

``` r
knn_group_results <- group_imp(obj, group = group_df, k = 20, dist_pow = 1, cores = 2)
#> Imputing 2 group(s) using KNN.
#> Running Mode: parallel (OpenMP within groups)...
knn_group_results
#> slideimp_results (KNN)
#> Dimensions: 20 x 100
#> 
#>          feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.3486482 0.7385414 0.4077444 0.1607935 0.3924661 0.2434143
#> sample2 0.5338935 0.4724364 0.9663621 0.4788070 0.5061132 0.3923603
#> sample3 0.7185848 0.7351035 0.6724479 0.3162537 0.7634236 1.0000000
#> sample4 0.1734418 0.0000000 0.0000000 0.0000000 0.0000000 0.2106358
#> sample5 0.5388440 0.5306182 0.5685354 0.5383513 0.4680080 0.8518388
#> sample6 0.3768380 0.5570723 0.8764909 0.5276245 0.6722794 0.5740639
#> 
#> # Showing [1:6, 1:6] of full matrix
```

- PCA imputation with `group_imp()` can be parallelized using the
  `{mirai}` package, similar to `tune_imp()`. `pin_blas = TRUE` can help
  PCA imputation on multi-threaded BLAS machines.

``` r
# Similar to `tune_imp`, parallelization is controlled by `mirai::daemons()`
mirai::daemons(2)
pca_group_results <- group_imp(obj, group = group_df, ncp = 10)
mirai::daemons(0)
```

- Alternatively, full matrix imputation can be performed using
  `knn_imp()` or `pca_imp()`.

``` r
full_knn_results <- knn_imp(obj = obj, k = 20)
full_pca_results <- pca_imp(obj = obj, ncp = 10)
```

## Sliding Window Imputation

- `slide_imp()` performs sliding window imputation.
- This function is **not** for Illumina DNAm microarrays; see
  `group_imp()`.
- **Note:**
  - DNAm WGBS/EM-seq data should be grouped by chromosomes (i.e., run
    `slide_imp()` separately on each chromosome) before imputation. See
    the package vignette for more details.
  - See vignette for details about selecting `window_size` and
    `overlap_size` with `tune_imp()`.

``` r
# Simulate some data
chr1_beta <- sim_mat(n = 10, p = 2000)$input
dim(chr1_beta)
#> [1]   10 2000
chr1_beta[1:5, 1:5]
#>          feature1  feature2  feature3   feature4  feature5
#> sample1 0.7500638 0.5323295 0.6095626 0.96762386 0.5149855
#> sample2 0.2809107 0.8695599        NA         NA 0.6040981
#> sample3 0.9409348 0.5445597 0.6432675 1.00000000 0.5613868
#> sample4 0.5946795 0.0000000        NA 0.07237333 0.4410413
#> sample5 0.8664253 0.6206139 0.3444691 0.52025046 0.5220036
```

- `slide_imp()` parameters:

  - `location`: **required** - a sorted numeric vector of length
    `ncol(obj)` giving the position of each column (e.g., genomic
    coordinates in *bp*).
  - `window_size`: **required** - width of each sliding window (same
    unit as `location`).
  - `overlap_size`: **optional** - overlap width between consecutive
    windows (same units as `location`). Must be strictly less than
    `window_size`.
  - `min_window_n`: **required** - minimum number of columns a window
    must contain to be imputed. Windows smaller than this are dropped.
    Must be greater than `k` (for KNN) or `ncp` (for PCA).
  - `dry_run`: **optional** - return just the calculated windows to
    examine which windows are included.
  - `k`: **required** - *(specifying KNN imputation)* number of nearest
    neighbors to use inside each window.
  - `ncp`: **required** - *(specifying PCA imputation)* number of
    principal components to retain. Use this instead of `k` when
    performing sliding-window PCA imputation.
  - `subset`: **optional** - impute just a subset of features (i.e.,
    just clock CpGs).
  - `flank`: **optional** - build flanking windows of `window_size`
    around features provided in `subset`.

- First, let’s perform a dry run to examine the windows that will be
  imputed by `slide_imp`.

``` r
location <- seq_len(ncol(chr1_beta)) # 1, 2, ..., 2000 for this simulated chromosome

slide_imp(
  obj = chr1_beta,
  location = location,
  window_size = 50, # select with `tune_imp()`
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

- Then, we can perform the imputation using the given parameters.

``` r
slide_imp(
  obj = chr1_beta,
  location = location,
  window_size = 50, # select with `tune_imp()`
  overlap_size = 5,
  min_window_n = 11, # must be > k = 10
  k = 10, # select with `tune_imp()`
  cores = 2,
  .progress = FALSE
)
#> slideimp_results (KNN)
#> Dimensions: 10 x 2000
#> 
#>          feature1  feature2  feature3   feature4  feature5  feature6
#> sample1 0.7500638 0.5323295 0.6095626 0.96762386 0.5149855 1.0000000
#> sample2 0.2809107 0.8695599 0.6324029 0.62469147 0.6040981 0.1583159
#> sample3 0.9409348 0.5445597 0.6432675 1.00000000 0.5613868 0.3054879
#> sample4 0.5946795 0.0000000 0.5837423 0.07237333 0.4410413 0.6101870
#> sample5 0.8664253 0.6206139 0.3444691 0.52025046 0.5220036 0.7464794
#> sample6 0.8157626 0.3053222 0.7227880 0.57711498 0.4576367 0.0000000
#> 
#> # Showing [1:6, 1:6] of full matrix
```

## Core functions

- `group_imp()`: Parallelizable group-wise (e.g., by chromosomes or
  column clusters) K-NN or PCA imputation with optional auxiliary
  features and group-wise parameters. The `group` argument supports:
  - Choices such as `"EPICv2"` or `"MSA"` (requires the
    [`{slideimp.extra}`](https://github.com/hhp94/slideimp.extra)
    package).
  - (Basic): a `data.frame` containing a `feature` column (character
    strings of feature names) and a `group` column (defining the group).
  - (Advanced): use `prep_groups()` to create groups based on a mapping
    data.frame (i.e., Illumina manifests) with custom group-wise
    parameters.
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
