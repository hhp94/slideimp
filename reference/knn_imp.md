# K-Nearest Neighbor Imputation for Numeric Matrices

Impute missing values in a numeric matrix using k-nearest neighbors
(K-NN).

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
  max_cache = 4,
  na_check = TRUE
)
```

## Arguments

- obj:

  A numeric matrix with **samples in rows** and **features in columns**.

- k:

  Integer. Number of nearest neighbors for imputation. 10 is a good
  starting point.

- colmax:

  Numeric. A number from 0 to 1. Threshold of column-wise missing data
  rate above which imputation is skipped.

- method:

  Character. Either "euclidean" (default) or "manhattan". Distance
  metric for nearest neighbor calculation.

- cores:

  Integer. Number of cores for K-NN parallelization (OpenMP). On macOS,
  OpenMP may need additional compiler configuration.

- post_imp:

  Boolean. Whether to impute remaining missing values (those that failed
  imputation) using column means.

- subset:

  Character. Vector of column names or integer vector of column indices
  specifying which columns to impute.

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

- na_check:

  Boolean. Check for leftover `NA` values in the results or not
  (internal use).

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

- **`tree = FALSE`** (default, brute-force K-NN): Always safe and
  usually faster for small to moderate data or high-dimensional cases.

- **`tree = TRUE`** (BallTree K-NN): Only use when imputation run time
  becomes prohibitive and missingness is low (\<5% missing).

- **Subset imputation**: Use the `subset` parameter for efficiency when
  only specific columns need imputation (e.g., epigenetic clock CpGs).

## References

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
Botstein D, Altman RB (2001). Missing value estimation methods for DNA
microarrays. Bioinformatics 17(6): 520-525.

## Examples

``` r
# Basic K-NN imputation
obj <- sim_mat(20, 20, perc_col_na = 1)$input
sum(is.na(obj))
#> [1] 40
result <- knn_imp(obj, k = 10)
result
#> Method: KNN imputation
#> Dimensions: 20 x 20
#> 
#>          feature1  feature2  feature3  feature4  feature5   feature6
#> sample1 0.5735747 0.7626304 0.9963202 0.8746018 0.3284892 0.65258962
#> sample2 0.7822131 0.5781703 0.5488623 0.7839500 0.9539087 0.56076182
#> sample3 0.1483321 0.0000000 0.3522334 0.0000000 0.1723249 0.04672928
#> sample4 0.4022741 0.2915725 0.4836761 0.4337060 0.7217295 0.11032665
#> sample5 0.6369624 0.6085569 0.5528170 0.7578386 0.2057219 0.60409777
#> sample6 0.6537480 0.3639402 0.5886658 0.8679417 0.6366236 0.22544169
#> 
#> # Showing [1:6, 1:6] of full matrix
```
