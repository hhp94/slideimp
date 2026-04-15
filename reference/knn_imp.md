# K-Nearest Neighbor Imputation for Numeric Matrices

Imputes missing values in a numeric matrix using k-nearest neighbors
(KNN).

## Usage

``` r
knn_imp(
  obj,
  k,
  colmax = 0.9,
  method = c("euclidean", "manhattan"),
  cores = 1,
  post_imp = TRUE,
  subset = NULL,
  dist_pow = 0,
  tree = FALSE,
  max_cache = 4
)
```

## Arguments

- obj:

  A numeric matrix with **samples in rows** and **features in columns**.

- k:

  Number of nearest neighbors for imputation. 10 is a good starting
  point.

- colmax:

  A number from 0 to 1. Threshold of column-wise missing data rate above
  which K-NN imputation is skipped.

- method:

  Either "euclidean" (default) or "manhattan". Distance metric for
  nearest neighbor calculation.

- cores:

  Number of cores for KNN parallelization (OpenMP). On macOS, OpenMP may
  need additional compiler configuration.

- post_imp:

  Whether to impute remaining missing values (those that failed K-NN
  imputation) using column means (default = `TRUE`).

- subset:

  Character vector of column names or integer vector of column indices
  specifying which columns to impute.

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

## Value

A numeric matrix of the same dimensions as `obj` with missing values
imputed.

## Details

This function performs imputation **column-wise** (using rows as
observations).

When `dist_pow > 0`, imputed values are computed as distance-weighted
averages where weights are inverse distances raised to the power of
`dist_pow`.

The `tree` parameter (when `TRUE`) uses a BallTree for faster neighbor
search via `{mlpack}` but **requires pre-filling** missing values with
column means. This can introduce a small bias when missingness is high.

## Performance Optimization

- **`tree = FALSE`** (default, brute-force KNN): Always safe and usually
  faster for small-to-moderate data or high-dimensional cases.

- **`tree = TRUE`** (BallTree KNN): Only use when imputation runtime
  becomes prohibitive and missingness is low (\<5% missing). Tree
  construction has overhead.

- **Subset imputation**: Use the `subset` parameter for efficiency when
  only specific columns need imputation (e.g., epigenetic clocks CpGs).

## References

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
Botstein D, Altman RB (2001). Missing value estimation methods for DNA
microarrays. Bioinformatics 17(6): 520-525.

## Examples

``` r
data(khanmiss1)
sum(is.na(khanmiss1))
#> [1] 1282

# Basic K-NN imputation (khanmiss1 has genes in rows, so transpose)
t_khanmiss1 <- t(khanmiss1)
result <- knn_imp(t_khanmiss1, k = 10)
result
#> slideimp_results (KNN)
#> Dimensions: 63 x 2308
#> 
#>           g1   g2   g3   g4   g5   g6
#> sample1 1873 1251  314 1324  776 1901
#> sample2   57 1350 1758 1428  476 1521
#> sample3   53 1140  162 1468  679   14
#> sample4 2059 1385 1857 1250  772 2052
#> sample5 1537 1261 1939 1666 1307 1705
#> sample6 1819 1526 1640 1795 1454 1643
#> 
#> # Showing [1:6, 1:6] of full matrix
```
