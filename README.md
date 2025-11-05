
<!-- README.md is generated from README.Rmd. Please edit that file -->

# slide_imp

<!-- badges: start -->

[![R-CMD-check](https://github.com/hhp94/slide_imp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hhp94/slide_imp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`{slide_imp}` is a lightweight R package for k-nearest neighbors (k-NN)
imputation of missing values in high-dimensional numeric matrices, such
as those from intensively sampled longitudinal data or epigenetics.
`knn_imp()` implements a full k-NN imputation. `slide_imp()` implements a
sliding window k-NN imputation for data where features are ordered
(e.g., by time or distance). `group_imp()` implements a group-wise
imputation for data that can be broken into meaningful groups, like
chromosomes or results of column clustering algorithms.

The package builds on the efficient k-NN imputation algorithm
(Bioconductor’s
[`{impute}`](https://www.bioconductor.org/packages/release/bioc/html/impute.html)
package) and adds enhancements: parallelization and tree-based methods
for speed, weighted imputation, multiple imputation strategies, and
built-in tuning tools. `{slide_imp}` expects matrices with samples in
rows and features in columns.

Key features include:

- **Sliding Window k-NN Imputation**: Break large data into overlapping
  windows for computationally feasible imputation while maintaining
  local structures (e.g., intensively sampled longitudinal data or
  epigenetics data).
- **Group k-NN Imputation**: Break large data into groups such as
  chromosomes or clusters identified by column clustering algorithms.
- **Full Matrix k-NN Imputation**: Standard k-NN for smaller data, with
  multi-core parallelization over columns with missing values.
- **Tree-Based k-NN**: Integration with
  [`{mlpack}`](https://www.mlpack.org/) for KD-Tree or Ball-Tree
  methods, accelerating imputation in high dimensions.
- **Subset Imputation**: Only impute a subset of columns to save time.
  Important for applications such as epigenetics clocks calculations.
- **Weighted Imputation**: Use inverse-distance weighting for more
  accurate averages, with tunable penalties.
- **Multiple Imputation**: Support for Predictive Mean Matching (PMM) or
  bootstrap resampling from nearest neighbors.
- **Fallback Imputation**: Optional post-k-NN mean imputation by column
  to handle remaining NAs.
- **Parameter Tuning**: Inject artificial NAs to evaluate and tune
  hyperparameters, with support for custom imputation functions.
- **Big Matrix Support**: Compatible with
  [`{bigmemory}`](https://doi.org/10.32614/CRAN.package.bigmemory) for
  file-backed matrices to handle massive data.

## Installation

The stable version of `{slide_imp}` can be installed from CRAN using:

``` r
install.packages("slide_imp")
```

You can install the development version of `{slide_imp}` from
[GitHub](https://github.com/hhp94/slide_imp) with:

``` r
install.packages("remotes")
remotes::install_github("hhp94/slide_imp")
```

## Full k-NN imputation `knn_imp()`

Use the workhorse function of the package, `knn_imp()` to perform full
k-NN imputation.

``` r
library(slide_imp)

data(khanmiss1)
# Transpose for samples in rows, features in columns
imputed_full <- knn_imp(t(khanmiss1), k = 3, method = "euclidean", cores = 1)
imputed_full
#> KnnImpList: List of 1 imputation of a 63 x 2308 matrix
#> 
#> Preview of imputation 1:
#>           g1   g2   g3   g4   g5   g6   g7   g8   g9  g10
#> sample1 1873 1251  314 1324  776 1901 2048 1513 1558 1796
#> sample2   57 1350 1758 1428  476 1521 2104   85 1784 1598
#> sample3   53 1140  162 1468  679   14 2048 1519 1631 1798
#> sample4 2059 1385 1857 1250  772 2052 2141 1969  243 2079
#> sample5 1537 1261 1939 1666 1307 1705 2137 1910 1499 1673
#> 
#> [ ... 58 more rows, 2298 more columns ]
```

`knn_imp()` yields the same results as `impute::impute.knn()` without
post imputation. `knn_imp()` is faster in larger data and is as fast in
smaller data given `cores = 1`. `knn_imp()` scales well with multiple
`cores`.

``` r
set.seed(1234)
obj_t <- t(khanmiss1)

bench::mark(
  # Single Core
  knn_imp(obj_t, k = 3, rowmax = 1, method = "euclidean")[[1]],
  # Multiple Cores
  knn_imp(obj_t, k = 3, rowmax = 1, cores = 4, method = "euclidean")[[1]],
  # impute::impute.knn
  t(impute::impute.knn(khanmiss1, k = 3, rowmax = 1, maxp = nrow(khanmiss1))$data),
  # Ensure results are the same between function calls
  check = TRUE
) |>
  dplyr::mutate(expression = c("knn_1", "knn_4", "impute.knn_1")) |>
  dplyr::select(-dplyr::where(is.list))
#> # A tibble: 3 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <chr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 knn_1         14.39ms   15.2ms      65.9    7.28MB     14.7
#> 2 knn_4          5.21ms    5.5ms     180.     7.28MB     43.7
#> 3 impute.knn_1  14.75ms   15.2ms      65.5   12.42MB     39.3
```

Using tree-based k-NN with weighted imputation and multiple imputation
(PMM):

``` r
imputed_tree <- knn_imp(
  t(khanmiss1),
  k = 5,
  tree = "kd", # KD-Tree via mlpack
  weighted = TRUE, # Inverse-distance weighting, further neighbors are less influential
  dist_pow = 2, # Harsher penalty for distant neighbors
  n_imp = 3, # 3 imputations
  n_pmm = 10 # PMM with 10 donors. Set n_pmm = 0 to use neighbor bootstrapping instead of PMM
)
imputed_tree
#> KnnImpList: List of 3 imputations of a 63 x 2308 matrix
#> 
#> Preview of imputation 1:
#>           g1   g2   g3   g4   g5   g6   g7   g8   g9  g10
#> sample1 1873 1251  314 1324  776 1901 2048 1513 1558 1796
#> sample2   57 1350 1758 1428  476 1521 2104   85 1784 1598
#> sample3   53 1140  162 1468  679   14 2048 1519 1631 1798
#> sample4 2059 1385 1857 1250  772 2052 2141 1969  243 2079
#> sample5 1537 1261 1939 1666 1307 1705 2137 1910 1499 1673
#> 
#> [ ... 58 more rows, 2298 more columns ]
#> 
#> [ ... and 2 more imputations ]
```

## Sliding window k-NN imputation with `slide_imp()`

Epigenetic datasets such as WGBS or EM-seq are very spatially
correlated, this method allows the imputation of the full epigenome with
good accuracy and speed by limiting the neighbor search space to sliding
windows.

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

imputed <- slide_imp(beta_matrix, n_feat = 500, n_overlap = 10, k = 10)
imputed
#> slide_impList: List of 1 imputation of a 10 x 1000 big.matrix
#> 
#> Preview of imputation 1:
#>     feat1  feat2  feat3  feat4  feat5  feat6  feat7  feat8  feat9 feat10
#> s1 0.3321 0.5573 0.6797 0.1593 0.5803 0.5920 0.4280 0.4323 0.4296 0.3802
#> s2 0.3047 0.5442 0.2516 0.5973 0.6081 0.1933 0.6456 0.4606 0.3892 0.4934
#> s3 0.3467 0.5617 0.4790 0.7008 0.2352 0.3350 0.2214 0.5486 0.4330 0.2052
#> s4 0.3858 0.5772 0.5437 0.7922 0.7523 0.4356 0.3132 0.5742 0.5731 0.5602
#> s5 0.3395 0.3664 0.5779 0.4325 0.6398 0.5558 0.5449 0.2616 0.5279 0.4761
#> 
#> [ ... 5 more rows, 990 more columns ]
```

Optional simple mean imputation as a fallback or baseline that enables
multi-step imputation strategies

``` r
# Inject extra NA into simulated data to make k-NN fail for first imputation
set.seed(1234)
obj <- t(sim_mat(n = 1000, m = 100, perc_NA = 0.8, perc_col_NA = 1)$input)
# Step 1: Column-wise imputation (impute the features using neighbor features
# of the same person/sample). Disable fall back imputation with `post_imp = FALSE`
imputed_by_col <- knn_imp(obj, cores = 4, k = 10, post_imp = FALSE)
# Step 2: Then if values are still missing, impute by rows. (impute the same features
# using values from OTHER people/samples).
imputed_by_row <- knn_imp(t(imputed_by_col[[1]]), cores = 4, k = 10, post_imp = FALSE)
# Step 3: Lastly, impute by column mean for any remaining missing.
imputed_mean <- mean_imp_col(t(imputed_by_row[[1]]))
sum(is.na(imputed_mean))
#> [1] 0
```

## Grouped imputation with `group_imp()`

Features can be grouped by chromosomes (in epigenetics) or by cluster
memberships identified through column clustering algorithms, such as
k-means. This reduces the nearest neighbor search space for each group,
significantly decreasing imputation time. The `group` data.frame is
defined by two list-columns: `features` and `aux` where each row is a
group. The `features` column includes all features to be imputed, while
the `aux` columns contains features used to assist with imputation but
excluded from the final results. For example, when imputing chrM, which
may have only a few CpGs, random CpGs from other chromosomes can be
included in the `aux` column to improve the quality of imputation.

``` r
# Simulate data from 2 chromosomes
set.seed(1234)
to_test <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1, nchr = 2)

# `group_1` will be all the CpGs on Chr1. Same for `group_2`
group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id

# Impute only first 3 values of group 1, the rest are aux. Group 2 does 4 features.
# Also optionally vary the parameters by group
group_df <- tibble::tibble(
  features = list(group_1[1:3], group_2[1:4]),
  aux = list(group_1, group_2),
  parameters = list(list(k = 3, weighted = TRUE), list(k = 4, method = "manhattan"))
)

group_df
#> # A tibble: 2 × 3
#>   features  aux        parameters      
#>   <list>    <list>     <list>          
#> 1 <chr [3]> <chr [28]> <named list [2]>
#> 2 <chr [4]> <chr [22]> <named list [2]>

# Run grouped imputation. t() to put features on the columns
obj <- t(to_test$input)
grouped_results <- group_imp(obj, group = group_df, k = 5)
#> Running with group-wise parameters
#> Imputing group 1/2
#> Imputing group 2/2
grouped_results
#> KnnImpList: List of 1 imputation of a 20 x 50 matrix
#> Imputed columns: 7 of 50 total columns
#> 
#> Preview of imputation 1:
#>     feat1  feat2  feat3  feat4  feat5  feat6  feat7  feat8  feat9 feat10
#> s1 0.2391 0.4129 0.7204 0.0000 0.5828 0.5989 0.3719 0.3778 0.5896 0.3057
#> s2 0.0000 0.2810 0.1601 0.1816 0.3774 0.5440 0.7931 0.2371 0.4718 0.1488
#> s3 0.5897 0.3678 0.5028 0.3609 0.2801 0.5280 0.2626 0.5283 0.5749     NA
#> s4 0.4201 0.5235 0.4435 0.3356 0.5047 0.7150 0.4443 0.3055 0.3366 0.4179
#> s5 0.4598 0.6388 0.5557 0.6394 0.5762 0.6491 0.8250     NA 0.6304 0.3444
#> 
#> [ ... 15 more rows, 40 more columns ]
```

## File-backed `big.matrix`

`slide_imp()` supports passing a `big.matrix` object or the path to its
descriptor file. Specify `output` for the result to change the backend
to using file-backed big.matrix as well to minimize memory at a cost of
performance.

**NOTE**: See `?restore_dimnames` if output’s dimnames are stripped.

``` r
library(bigmemory)

# IMPORTANT: Enable dimnames support for big.matrix objects
options(bigmemory.allow.dimnames = TRUE)

# Load example data
data(khanmiss1)
mat <- t(khanmiss1) # samples as rows, features as columns
temp_dir <- withr::local_tempdir()
```

### Create file-backed big.matrix

``` r
# Convert t(khanmiss1) to big.matrix with backing files for large data
big_mat <- bigmemory::as.big.matrix(
  mat,
  type = "double",
  backingfile = "khan.bin",
  descriptorfile = "khan.desc",
  backingpath = temp_dir
)
```

### Method 1: Impute using big.matrix object

Impute a `bigmemory::big.matrix` and optionally save results to
file-backed big.matrix

``` r
imputed_obj <- slide_imp(
  obj = big_mat,
  n_feat = 100,
  n_overlap = 10,
  k = 10,
  overwrite = TRUE, # Overwrite any existing results
  output = file.path(temp_dir, "imputed.bin")
)

# Access the imputed data (returns list of big.matrix objects)
imputed_obj
#> slide_impList: List of 1 imputation of a 63 x 2308 big.matrix
#> 
#> Preview of imputation 1:
#>           g1   g2   g3   g4   g5   g6   g7   g8   g9  g10
#> sample1 1873 1251  314 1324  776 1901 2048 1513 1558 1796
#> sample2   57 1350 1758 1428  476 1521 2104   85 1784 1598
#> sample3   53 1140  162 1468  679   14 2048 1519 1631 1798
#> sample4 2059 1385 1857 1250  772 2052 2141 1969  243 2079
#> sample5 1537 1261 1939 1666 1307 1705 2137 1910 1499 1673
#> 
#> [ ... 58 more rows, 2298 more columns ]
```

### ⚠️ Regarding objects pointing to the same bigmemory matrices

On Windows, file-backed objects hold file locks. You MUST remove any
variables that point to the same files before overwriting them. This
only removes the R variables, not the files on disk.

In general, always remove references to big.matrix objects before
attempting to overwrite or modify their backing files.

``` r
rm(imputed_obj)
invisible(gc(verbose = FALSE)) # Force garbage collection to release file handles
```

### Method 2: Impute using descriptor file path

Alternatively, pass descriptor file path instead of object. This is
useful for distributed computing or when the object isn’t in memory.

``` r
desc_path <- file.path(temp_dir, "khan.desc")

imputed_path <- slide_imp(
  obj = desc_path, # Using path instead of object
  n_feat = 100,
  n_overlap = 10,
  k = 10,
  overwrite = TRUE,
  output = file.path(temp_dir, "imputed.bin")
)

imputed_path
#> slide_impList: List of 1 imputation of a 63 x 2308 big.matrix
#> 
#> Preview of imputation 1:
#>           g1   g2   g3   g4   g5   g6   g7   g8   g9  g10
#> sample1 1873 1251  314 1324  776 1901 2048 1513 1558 1796
#> sample2   57 1350 1758 1428  476 1521 2104   85 1784 1598
#> sample3   53 1140  162 1468  679   14 2048 1519 1631 1798
#> sample4 2059 1385 1857 1250  772 2052 2141 1969  243 2079
#> sample5 1537 1261 1939 1666 1307 1705 2137 1910 1499 1673
#> 
#> [ ... 58 more rows, 2298 more columns ]
```

## Parameters Tuning with `tune_imp()`

Use `tune_imp()` to tune hyperparameters by injecting artificial NAs and
evaluating imputation accuracy. This function supports built-in methods
or custom functions.

``` r
# This tibble defines the hyperparameters to tune for. `knn_imp` requires `k`, but other parameters like `dist_pow` and `method` can also be tuned.
parameters <- dplyr::tibble(
  k = c(5, 10),
  method = "euclidean",
  weighted = TRUE,
  dist_pow = c(2, 5),
  post_imp = TRUE
)
parameters
#> # A tibble: 2 × 5
#>       k method    weighted dist_pow post_imp
#>   <dbl> <chr>     <lgl>       <dbl> <lgl>   
#> 1     5 euclidean TRUE            2 TRUE    
#> 2    10 euclidean TRUE            5 TRUE

obj_t <- t(khanmiss1)

# 3 repeats, each time inject 100 NAs
results <- tune_imp(obj_t, parameters, rep = 3, .f = "knn_imp", num_na = 100)

# Compute metrics with {yardstick}
library(yardstick)
met_set <- metric_set(mae, rmse, rsq)
results$metrics <- lapply(results$result, function(x) met_set(x, truth = truth, estimate = estimate))
head(
  dplyr::select(
    tidyr::unnest(dplyr::select(results, -result), cols = "metrics"),
    all_of(names(parameters)), contains(".")
  )
)
#> # A tibble: 6 × 8
#>       k method    weighted dist_pow post_imp .metric .estimator .estimate
#>   <dbl> <chr>     <lgl>       <dbl> <lgl>    <chr>   <chr>          <dbl>
#> 1     5 euclidean TRUE            2 TRUE     mae     standard     344.   
#> 2     5 euclidean TRUE            2 TRUE     rmse    standard     492.   
#> 3     5 euclidean TRUE            2 TRUE     rsq     standard       0.390
#> 4    10 euclidean TRUE            5 TRUE     mae     standard     354.   
#> 5    10 euclidean TRUE            5 TRUE     rmse    standard     495.   
#> 6    10 euclidean TRUE            5 TRUE     rsq     standard       0.380
```

To tune a custom imputation function, a function has to take an numeric
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
#> 1     0     1 mae     standard      0.869 
#> 2     0     1 rmse    standard      0.981 
#> 3     0     1 rsq     standard      0.0480
#> 4     1     2 mae     standard      1.71  
#> 5     1     2 rmse    standard      2.10  
#> 6     1     2 rsq     standard      0.130
```

For more details, see the function documentation (e.g., `?slide_imp`,
`?group_imp`, `?knn_imp`, `?tune_imp`).

# Developer notes for [`{mlpack}`](https://www.mlpack.org/)

`{mlpack}` is an awesome library and works seamlessly with integration
in R. For other developers who are looking to integrate `{mlpack}` into
their own work, make sure to structure the .cpp file as follows. Most
importantly is the inclusion of `<mlpack.h>` first before including any
other `mlpack` headers. This appropriately directs the output of
`std::cout` and `std::cerr` to R and won’t trip `R CMD check` (i.e.,
`'_ZSt4cerr', possibly from 'std::cerr' (C++)` or
`'_ZSt4cout', possibly from 'std::out' (C++)`). Full implementation is
found under `src/impute_knn_mlpack.cpp`. I would like to thank the
developers and the community for the maintenance and development of
`mlpack`. I would like to especially thank Dirk Eddelbuettel for his
development of the R + C++ environment and guidance throughout the
development of this package.

``` cpp
// [[Rcpp::depends(mlpack, RcppArmadillo)]]
#include <mlpack.h>
#include "imputed_value.h"
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>
#include <RcppArmadillo.h>
```
