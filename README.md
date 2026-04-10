
<!-- README.md is generated from README.Rmd. Please edit that file -->

# slideimp

<!-- badges: start -->

[![R-CMD-check](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`{slideimp}` is a lightweight R package for fast K-NN and PCA imputation
of missing values in high-dimensional numeric matrices.

**Core functions**

- `knn_imp()`: Full-matrix K-NN imputation with multi-core
  parallelization, [`{mlpack}`](https://mlpack.org/) KD/Ball-Tree
  nearest neighbor implementation (for data with very low missing rates
  and extremely high dimensions), and optional subset imputation (ideal
  for epigenetic clock calculations).
- `pca_imp()`: Optimized version of
  [`missMDA::imputePCA()`](http://factominer.free.fr/missMDA/PCA.html)
  for high-dimensional numeric matrices.
- `slide_imp()`: Sliding window K-NN or PCA imputation for extremely
  high-dimensional numeric matrices with ordered features (i.e., by
  genomic position).
- `group_imp()`: Parallelizable group-wise (e.g., by chromosomes or
  column clusters) K-NN or PCA imputation with optional auxiliary
  features and group-wise parameters.
  - `group_features()`: `group_imp()`’s helper function to create groups
    based on a mapping data.frame (i.e., Illumina manifests). See
    [`{slideimp.extra}`](https://github.com/hhp94/slideimp.extra) on
    GitHub for tools to process common Illumina manifests.
- `tune_imp()`: Parallelizable hyperparameter tuning with repeated
  cross-validation; works with built-in or custom imputation functions.

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
# Simulate data from 2 chromosomes
set.seed(1234)
sim_obj <- sim_mat(m = 20, n = 50, perc_total_na = 0.3, perc_col_na = 1, nchr = 2)
# Here we see that variables are stored in rows
sim_obj$input[1:5, 1:5]
#>              s1        s2        s3        s4        s5
#> feat1 0.2391314 0.0000000 0.5897476 0.4201222        NA
#> feat2        NA 0.2810446 0.3677927        NA 0.6387734
#> feat3 0.7203854 0.1600776 0.5027545        NA 0.5556735
#> feat4 0.0000000 0.1816453 0.3608640 0.3356484 0.6394179
#> feat5 0.5827582 0.3774313 0.2801131 0.5047049 0.5761809

# So we t() to put the variables in columns
obj <- t(sim_obj$input)
```

We can optionally estimate the prediction accuracy of different methods
and tune hyperparameters prior to imputation with `tune_imp()`.

For custom functions (`.f` argument), the `parameters` data.frame must
include the columns corresponding to the arguments passed to the custom
function. The custom function must accept `obj` as the first argument
and return a matrix with the same dimensions as `obj`.

We tune the results using 2 repeats (`rep = 2`) for illustration
(increase in actual analyses).

``` r
knn_params <- tibble::tibble(k = c(5, 20))
# Parallelization with OpenMP is controlled by `cores` only for knn or slideimp knn with OpenMP only.
tune_knn <- tune_imp(obj, parameters = knn_params, .f = "knn_imp", cores = 2, rep = 2)
#> Tuning knn_imp
#> Step 1/2: Injecting NA
#> Running Mode: parallel...
#> Step 2/2: Tuning
compute_metrics(tune_knn)
#> # A tibble: 8 × 9
#>       k cores param_set   rep     n n_miss .metric .estimator .estimate
#>   <dbl> <dbl>     <int> <int> <int>  <int> <chr>   <chr>          <dbl>
#> 1     5     2         1     1   100      0 mae     standard       0.178
#> 2     5     2         1     1   100      0 rmse    standard       0.225
#> 3    20     2         2     1   100      0 mae     standard       0.149
#> 4    20     2         2     1   100      0 rmse    standard       0.190
#> 5     5     2         1     2   100      0 mae     standard       0.202
#> 6     5     2         1     2   100      0 rmse    standard       0.259
#> 7    20     2         2     2   100      0 mae     standard       0.172
#> 8    20     2         2     2   100      0 rmse    standard       0.219
```

For K-NN without OpenMP, PCA, and custom functions, setup
parallelization with `mirai::daemons()`.

``` r
mirai::daemons(2) # 2 Cores
# Note, for PCA and custom functions, cores is controlled by the `mirai::daemons()`
# and the `cores` argument is ignored.

# PCA imputation.
pca_params <- tibble::tibble(ncp = c(1, 5))
tune_pca <- tune_imp(obj, parameters = pca_params, .f = "pca_imp", rep = 2)

# The parameters have `mean` and `sd` columns.
custom_params <- tibble::tibble(mean = 1, sd = 0)
# This function impute data with rnorm values of different `mean` and `sd`.
custom_function <- function(obj, mean, sd) {
  missing <- is.na(obj)
  obj[missing] <- rnorm(sum(missing), mean = mean, sd = sd)
  return(obj)
}
tune_custom <- tune_imp(obj, parameters = custom_params, .f = custom_function, rep = 2)

mirai::daemons(0) # Close daemons
```

Then, preferably perform imputation by group with `group_imp()` if the
variables can be meaningfully grouped (e.g., by chromosomes).

- `group_imp()` allows imputation to be performed separately within
  defined groups (e.g., by chromosome), which significantly reduces
  runtime and can increase accuracy for both K-NN and PCA imputation.
- `group_imp()` requires a `group` tibble, *preferably* created with
  `group_features()`, with three list-columns:
  - `features`: **required** - a list-column where each element is a
    character vector of variable names to be imputed together.
  - `aux`: **optional** - auxiliary variables to include in each group.
  - `parameters`: **optional** - group-specific imputation parameters.
- In this example, we have data from 2 chromosomes so the `group` tibble
  should have **two rows** (one per chromosome), with the corresponding
  variables listed in the `features` column for each row.

PCA-based imputation with `group_imp()` can be parallelized using the
`{mirai}` package, similar to how parallelization is done with
`tune_imp()`.

``` r
# Use the `group_features()` helper function
group_df <- group_features(obj, sim_obj$col_group)
group_df

# We choose K-NN imputation, k = 5, from the `tune_imp` results.
knn_group_results <- group_imp(obj, group = group_df, k = 5, cores = 2)

# Similar to `tune_imp`, parallelization is controlled by `mirai::daemons()`
mirai::daemons(2)
knn_group_results <- group_imp(obj, group = group_df, ncp = 3)
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
imputation. See the package vignette for more details.

``` r
chr1_beta <- t(sim_mat(m = 10, n = 2000, perc_total_na = 0.3, perc_col_na = 1, nchr = 1)$input)
dim(chr1_beta)
#> [1]   10 2000
chr1_beta[1:5, 1:5]
#>        feat1     feat2     feat3     feat4     feat5
#> s1        NA 0.7297743        NA        NA 0.3968039
#> s2 0.7346970        NA 0.5669140 0.3236858 0.3932419
#> s3        NA        NA        NA 0.3108793        NA
#> s4 0.5401526 0.5779956 0.4271064        NA 0.3309645
#> s5 0.6457875        NA 0.7308792 0.4803642 0.5929590
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
- `k`: **required** - *(specifying KNN imputation)* Number of nearest
  neighbors to use inside each window.
- `ncp`: **required** - *(specifying PCA imputation)* Number of
  principal components to retain. Use this instead of `k` when
  performing sliding-window PCA imputation.

Use the function `compute_windows()` to examine the windows calculated
by `slide_imp()` before imputation.

``` r
location <- seq_len(ncol(chr1_beta)) # 1, 2, ..., 2000 for this simulated chromosome

# Examine the windows calculated by `slide_imp()`
windows <- compute_windows(
  location,
  window_size = 50,
  overlap_size = 5,
  min_window_n = 11
)
windows
#> # A tibble: 45 × 4
#>    start   end window_n keep 
#>    <int> <int>    <int> <lgl>
#>  1     1    50       50 TRUE 
#>  2    46    95       50 TRUE 
#>  3    91   140       50 TRUE 
#>  4   136   185       50 TRUE 
#>  5   181   230       50 TRUE 
#>  6   226   275       50 TRUE 
#>  7   271   320       50 TRUE 
#>  8   316   365       50 TRUE 
#>  9   361   410       50 TRUE 
#> 10   406   455       50 TRUE 
#> # ℹ 35 more rows

# Then do the imputation
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
#> SlideImpImputedMatrix (KNN)
#> Dimensions: 10 x 2000
#> 
#>        feat1     feat2     feat3     feat4     feat5
#> s1 0.5067435 0.7297743 0.5884198 0.5063839 0.3968039
#> s2 0.7346970 0.4551576 0.5669140 0.3236858 0.3932419
#> s3 0.5625864 0.4790436 0.5316400 0.3108793 0.5234974
#> s4 0.5401526 0.5779956 0.4271064 0.5551127 0.3309645
#> s5 0.6457875 0.4006866 0.7308792 0.4803642 0.5929590
#> 
#> # Showing [1:5, 1:5] of full matrix
```
