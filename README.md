
<!-- README.md is generated from README.Rmd. Please edit that file -->

# slideimp

<!-- badges: start -->

[![R-CMD-check](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`{slideimp}` is a lightweight R package for fast k-NN and PCA imputation
of missing values in high-dimensional numeric matrices (e.g., intensive
longitudinal or epigenetic data).

**Core functions**

- `knn_imp()`: full-matrix k-NN with multi-core parallelization,
  `{mlpack}` tree acceleration (KD/Ball-Tree, best when \<20% missing
  per feature), and optional subset imputation (ideal for epigenetic
  clock calculations).
- `pca_imp()`: optimized version of
  [`missMDA::imputePCA()`](http://factominer.free.fr/missMDA/PCA.html)
  for high-dimensional numeric matrices.
- `slide_imp()`: sliding-window k-NN or PCA imputation of extremely high
  dimensional numeric matrices with ordered features (by time or genomic
  position).
- `group_imp()`: group-wise (e.g. by chromosomes, or column clusters)
  k-NN or PCA imputation with optional auxiliary features for small
  groups.
- `tune_imp()`: hyperparameter tuning with cross-validation; works with
  built-in or custom imputation functions.

## Installation

The stable version of `{slideimp}` can be installed from CRAN using:

``` r
install.packages("slideimp")
```

You can install the development version of `{slideimp}` from
[GitHub](https://github.com/hhp94/slideimp) with:

``` r
install.packages("remotes")
remotes::install_github("hhp94/slideimp")
```

## Full k-NN imputation `knn_imp()` and PCA imputation with `pca_imp()`

Use the workhorse function of the package, `knn_imp()` to perform full
k-NN imputation.

``` r
library(slideimp)

data(khanmiss1)

# Transpose for samples in rows, features in columns
imputed_full <- knn_imp(t(khanmiss1), k = 3, method = "euclidean", cores = 1)
imputed_full
#> ImputedMatrix (KNN)
#> Dimensions: 63 x 2308
#> 
#>           g1   g2   g3   g4   g5
#> sample1 1873 1251  314 1324  776
#> sample2   57 1350 1758 1428  476
#> sample3   53 1140  162 1468  679
#> sample4 2059 1385 1857 1250  772
#> sample5 1537 1261 1939 1666 1307
#> 
#> # Showing [1:5, 1:5] of full matrix
```

`knn_imp()` yields the same results as
[`impute::impute.knn()`](https://bioconductor.org/packages/release/bioc/html/impute.html)
without post imputation. `knn_imp()` is faster in larger data and is as
fast in smaller data given `cores = 1`. `knn_imp()` scales well with
multiple `cores`.

``` r
library(bench)

set.seed(1234)

mark(
  # Single Core
  knn_imp(t(khanmiss1), k = 3, method = "euclidean")[, ],
  # Multiple Cores
  knn_imp(t(khanmiss1), k = 3, cores = 4, method = "euclidean")[, ],
  # impute::impute.knn
  t(impute::impute.knn(khanmiss1, k = 3, maxp = nrow(khanmiss1))$data),
  # Ensure results are the same between function calls
  check = TRUE
) |>
  dplyr::mutate(expression = c("knn_1", "knn_4", "impute.knn_1")) |>
  dplyr::select(-dplyr::where(is.list))
#> # A tibble: 3 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <chr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 knn_1         13.62ms  14.34ms      70.1     9.5MB     38.2
#> 2 knn_4          4.54ms   4.95ms     202.      9.5MB     70.1
#> 3 impute.knn_1  14.66ms  14.81ms      66.8    12.4MB     36.7
```

Using tree-based k-NN with weighted imputation:

``` r
imputed_tree <- knn_imp(
  t(khanmiss1),
  k = 5,
  tree = "kd", # KD-Tree via mlpack
  dist_pow = 2
)
imputed_tree
#> ImputedMatrix (KNN)
#> Dimensions: 63 x 2308
#> 
#>           g1   g2   g3   g4   g5
#> sample1 1873 1251  314 1324  776
#> sample2   57 1350 1758 1428  476
#> sample3   53 1140  162 1468  679
#> sample4 2059 1385 1857 1250  772
#> sample5 1537 1261 1939 1666 1307
#> 
#> # Showing [1:5, 1:5] of full matrix
```

Similarly, `pca_imp()` is an optimized version of
[`missMDA::imputePCA()`](http://factominer.free.fr/missMDA/PCA.html) for
high-dimensional numeric matrix.

``` r
bench::mark(
  pca_imp(t(khanmiss1)[, 1:500], ncp = 2)[, ],
  missMDA::imputePCA(t(khanmiss1)[, 1:500], ncp = 2)$completeObs,
  check = TRUE
) |>
  dplyr::mutate(expression = c("pca_imp", "imputePCA")) |>
  dplyr::select(-dplyr::where(is.list))
#> Warning: Some expressions had a GC in every iteration; so filtering is
#> disabled.
#> # A tibble: 2 × 6
#>   expression      min   median `itr/sec` mem_alloc `gc/sec`
#>   <chr>      <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 pca_imp      13.8ms   15.4ms      40.9    36.6MB     35.1
#> 2 imputePCA    34.7ms   36.7ms      20.7   164.8MB     28.2
```

## Sliding window k-NN and PCA imputation with `slide_imp()`

Epigenetic datasets such as WGBS or EM-seq are very spatially
correlated, `slide_imp()` allows the imputation of the full epigenome
with good accuracy and speed by limiting the neighbor search space to
sliding windows.

Use k-NN by specifying `k`, PCA by specifying `ncp`. Only arguments
applicable to each method are applied.

``` r
# Simulating WGBS data with 1000 CpGs and 10 samples
set.seed(1234)
beta_matrix <- t(sim_mat(m = 10, n = 1000, perc_NA = 0.3, perc_col_NA = 1)$input)
beta_matrix[1:10, 1:5]
#>         feat1     feat2     feat3     feat4     feat5
#> s1  0.3320706 0.5572683 0.6796919 0.1593403 0.5802804
#> s2  0.3046741 0.5442470 0.2515999 0.5973359 0.6080809
#> s3  0.3467045        NA        NA 0.7007942 0.2352033
#> s4         NA 0.5771749 0.5437410 0.7921855 0.7523102
#> s5  0.3395358 0.3664415        NA        NA 0.6398042
#> s6  0.3144117        NA 0.2995698 0.3984796 0.3572444
#> s7  0.4963455 0.4612003 0.3311917        NA        NA
#> s8  0.5039341 0.5381084 0.3811143 0.6207009 0.6710451
#> s9         NA        NA 0.3455449        NA        NA
#> s10        NA 0.5498897        NA 0.3738964        NA

# Sliding Window k-NN imputation by specifying `k`
knn_imputed <- slide_imp(beta_matrix, n_feat = 500, n_overlap = 10, k = 10)
knn_imputed
#> ImputedMatrix (KNN)
#> Dimensions: 10 x 1000
#> 
#>        feat1     feat2     feat3     feat4     feat5
#> s1 0.3320706 0.5572683 0.6796919 0.1593403 0.5802804
#> s2 0.3046741 0.5442470 0.2515999 0.5973359 0.6080809
#> s3 0.3467045 0.5617059 0.4789529 0.7007942 0.2352033
#> s4 0.3858128 0.5771749 0.5437410 0.7921855 0.7523102
#> s5 0.3395358 0.3664415 0.5779475 0.4324594 0.6398042
#> 
#> # Showing [1:5, 1:5] of full matrix

# Sliding Window PCA imputation by specifying `ncp`
pca_imputed <- slide_imp(beta_matrix, n_feat = 500, n_overlap = 10, ncp = 2)
pca_imputed
#> ImputedMatrix (PCA)
#> Dimensions: 10 x 1000
#> 
#>        feat1     feat2     feat3     feat4     feat5
#> s1 0.3320706 0.5572683 0.6796919 0.1593403 0.5802804
#> s2 0.3046741 0.5442470 0.2515999 0.5973359 0.6080809
#> s3 0.3467045 0.4949746 0.3892889 0.7007942 0.2352033
#> s4 0.3781874 0.5771749 0.5437410 0.7921855 0.7523102
#> s5 0.3395358 0.3664415 0.3919350 0.5182620 0.6398042
#> 
#> # Showing [1:5, 1:5] of full matrix
```

## Grouped imputation with `group_imp()`

Features can be grouped by chromosomes (in epigenetics) or by cluster
memberships identified through column clustering algorithms, such as
k-means. This significantly decreases imputation time and potentially
increases imputation accuracy.

`group_imp()` requires the `group` data.frame (preferably a `tibble` for
the handling of list columns), which is defined by two list-columns:
`features` and `aux` (optional) where each row is a group. The
`features` column includes all features to be imputed, while the `aux`
columns contain features used to assist with imputation but excluded
from the final results. For example, when imputing chrM, which may have
only a few CpGs, random CpGs from other chromosomes can be included in
the `aux` column to improve the quality of imputation.

Similar to `slide_imp()`, use k-NN by specifying `k`, PCA by specifying
`ncp`.

``` r
# Simulate data from 2 chromosomes
set.seed(1234)
to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1, nchr = 2)

# `group_1` will be all the CpGs on Chr1. Same for `group_2`
group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

# Impute only first 3 columns of group 1, the rest are aux. Group 2 does 4 features.
# Also optionally vary the parameters by group
knn_df <- tibble::tibble(
  features = list(group_1[1:3], group_2[1:4]),
  aux = list(group_1, group_2),
  parameters = list(list(k = 3, dist_pow = 0), list(k = 4, method = "manhattan"))
)

knn_df
#> # A tibble: 2 × 3
#>   features  aux        parameters      
#>   <list>    <list>     <list>          
#> 1 <chr [3]> <chr [28]> <named list [2]>
#> 2 <chr [4]> <chr [22]> <named list [2]>

# Run grouped imputation. t() to put features on the columns
obj <- t(to_test$input)
knn_grouped <- group_imp(obj, group = knn_df, k = 5)
#> Running with group-wise parameters
#> Imputing group 1/2
#> Imputing group 2/2
knn_grouped
#> ImputedMatrix (KNN)
#> Dimensions: 20 x 50
#> 
#>        feat1     feat2     feat3     feat4     feat5
#> s1 0.2391314 0.4411043 0.7203854 0.0000000 0.5827582
#> s2 0.0000000 0.2810446 0.1600776 0.1816453 0.3774313
#> s3 0.5897476 0.3677927 0.5027545 0.3608640 0.2801131
#> s4 0.4201222 0.5643869 0.4419687 0.3356484 0.5047049
#> s5 0.4597799 0.6387734 0.5556735 0.6394179 0.5761809
#> 
#> # Showing [1:5, 1:5] of full matrix

# Each group will be imputed with `pca_imp` and ncp = 2
pca_df <- tibble::tibble(
  features = list(group_1[1:3], group_2[1:4])
)

pca_grouped <- group_imp(obj, group = pca_df, ncp = 2)
#> Running with the same parameters for all groups
#> Imputing group 1/2
#> Imputing group 2/2
pca_grouped
#> ImputedMatrix (PCA)
#> Dimensions: 20 x 50
#> 
#>        feat1     feat2     feat3     feat4     feat5
#> s1 0.2391314 0.5727852 0.7203854 0.0000000 0.5827582
#> s2 0.0000000 0.2810446 0.1600776 0.1816453 0.3774313
#> s3 0.5897476 0.3677927 0.5027545 0.3608640 0.2801131
#> s4 0.4201222 0.5628569 0.4786760 0.3356484 0.5047049
#> s5 0.4382530 0.6387734 0.5556735 0.6394179 0.5761809
#> 
#> # Showing [1:5, 1:5] of full matrix
```

## Parameter Tuning with `tune_imp()`

Use `tune_imp()` to tune hyperparameters by cross-validation. This
function repeatedly injects `NA` and evaluates imputation accuracy. This
function supports built-in methods or custom functions. Available
built-in methods are `slide_imp`, `knn_imp`, and `pca_imp`.

``` r
# This tibble defines the hyperparameters to tune for. `knn_imp` requires `k`, but other parameters like `dist_pow` and `method` can also be tuned.
parameters <- dplyr::tibble(
  k = c(5, 10),
  method = "euclidean",
  dist_pow = c(2, 5)
)
parameters
#> # A tibble: 2 × 3
#>       k method    dist_pow
#>   <dbl> <chr>        <dbl>
#> 1     5 euclidean        2
#> 2    10 euclidean        5

obj_t <- t(khanmiss1)

# 3 repeats, each time inject 100 NAs
results <- tune_imp(obj_t, parameters, rep = 3, .f = "knn_imp", num_na = 100)

# Compute metrics with {yardstick}
library(yardstick)
met_set <- metric_set(mae, rmse, rsq)
results$metrics <- lapply(
  results$result,
  function(x) met_set(x, truth = truth, estimate = estimate)
)

head(
  dplyr::select(
    tidyr::unnest(dplyr::select(results, -result), cols = "metrics"),
    all_of(names(parameters)), contains(".")
  )
)
#> # A tibble: 6 × 6
#>       k method    dist_pow .metric .estimator .estimate
#>   <dbl> <chr>        <dbl> <chr>   <chr>          <dbl>
#> 1     5 euclidean        2 mae     standard     434.   
#> 2     5 euclidean        2 rmse    standard     553.   
#> 3     5 euclidean        2 rsq     standard       0.251
#> 4    10 euclidean        5 mae     standard     420.   
#> 5    10 euclidean        5 rmse    standard     527.   
#> 6    10 euclidean        5 rsq     standard       0.315
```

To tune a custom imputation function, a function has to take a numeric
matrix `obj` as the first argument, return a matrix with the same
dimension as `obj`, and parameters as arguments.

``` r
# This custom function imputes NAs with rnorm values. First argument has to be `obj`, which is a numeric matrix. The returned object has to be the same dimension and data type as `obj`.
custom_imp <- function(obj, mean = 0, sd = 1) {
  na_pos <- is.na(obj)
  obj[na_pos] <- rnorm(sum(na_pos), mean = mean, sd = sd)
  return(obj)
}

# This tibble defines the hyperparameters grid similar to the built-in case.
parameters_custom <- dplyr::tibble(
  mean = c(0, 1),
  sd = c(1, 2)
)

# Similarly, use `{yardstick}` to get the prediction metrics
results_custom <- tune_imp(obj, parameters_custom, .f = custom_imp, rep = 2, num_na = 20)
results_custom$metrics <- lapply(
  results_custom$result,
  function(x) met_set(x, truth = truth, estimate = estimate)
)

head(
  dplyr::select(
    tidyr::unnest(dplyr::select(results_custom, -result), cols = "metrics"),
    dplyr::all_of(names(parameters_custom)), dplyr::contains(".")
  )
)
#> # A tibble: 6 × 5
#>    mean    sd .metric .estimator .estimate
#>   <dbl> <dbl> <chr>   <chr>          <dbl>
#> 1     0     1 mae     standard      0.999 
#> 2     0     1 rmse    standard      1.15  
#> 3     0     1 rsq     standard      0.0183
#> 4     1     2 mae     standard      1.99  
#> 5     1     2 rmse    standard      2.49  
#> 6     1     2 rsq     standard      0.0969
```

For more details, see the function documentation (e.g., `?slide_imp`,
`?group_imp`, `?pca_imp`, `?knn_imp`, `?tune_imp`).
