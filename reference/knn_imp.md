# K-Nearest Neighbors Imputation for Numeric Matrices

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
  na_check = TRUE,
  .progress = FALSE
)
```

## Arguments

- obj:

  A numeric matrix with samples in rows and features in columns.

- k:

  Integer. Number of nearest neighbors to use for imputation.

- colmax:

  Numeric scalar between `0` and `1`. Columns with a missing-data
  proportion greater than `colmax` are not imputed.

- method:

  Character. Distance metric for nearest-neighbor calculation: either
  `"euclidean"` or `"manhattan"`.

- cores:

  Integer. Number of cores to use for parallel computation. Defaults to
  `1`.

- post_imp:

  Logical. If `TRUE`, replace any remaining missing values with column
  means after K-NN imputation.

- subset:

  Optional character or integer vector specifying columns to impute. If
  `NULL`, all eligible columns are imputed.

- dist_pow:

  Numeric. Power used to penalize more distant neighbors in the weighted
  average. `dist_pow = 0` gives an unweighted average of the nearest
  neighbors.

- tree:

  Logical. If `FALSE`, use brute-force K-NN. If `TRUE`, use ball-tree
  K-NN via `mlpack`.

- na_check:

  Logical. If `TRUE`, check whether the result still contains missing
  values.

- .progress:

  Logical. If `TRUE`, show imputation progress.

## Value

A numeric matrix of the same dimensions as `obj`, with missing values
imputed. The returned object has class `slideimp_results`.

## Details

`knn_imp()` performs imputation column-wise, treating rows as
observations and columns as features.

When `dist_pow > 0`, imputed values are computed as distance-weighted
averages. Weights are inverse distances raised to the power of
`dist_pow`.

If `tree = TRUE`, nearest neighbors are found with a ball tree via the
`mlpack` package. This can be faster for some large, low-missingness
data sets, but it requires initially filling missing values with column
means, which can introduce bias when missingness is high.

## Performance optimization

- `tree = FALSE` uses brute-force K-NN. This is always safe and is often
  faster for small to moderate data sets or high-dimensional data.

- `tree = TRUE` uses ball-tree K-NN. Consider this only when run time is
  prohibitive and missingness is low, for example less than 5%.

- Use `subset` when only specific columns need imputation.

## References

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
Botstein D, Altman RB (2001). Missing value estimation methods for DNA
microarrays. *Bioinformatics*, 17(6), 520-525.
[doi:10.1093/bioinformatics/17.6.520](https://doi.org/10.1093/bioinformatics/17.6.520)

## Examples

``` r
set.seed(123)
obj <- sim_mat(20, 20, perc_col_na = 1)$input
sum(is.na(obj))
#> [1] 40

result <- knn_imp(obj, k = 10, .progress = FALSE)
result
#> Method: KNN imputation
#> Dimensions: 20 x 20
#> 
#>           feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.08885928 0.0946424 0.5010982 0.2257198 0.3310293 0.3919172
#> sample2 0.35087647 0.2569208 0.4441198 0.3953534 0.5894579 0.2589739
#> sample3 0.56864697 0.4021824 0.7354948 0.6422099 0.8454413 0.6652821
#> sample4 0.30420093 0.7886995 0.3732225 0.5291319 0.5289648 0.4384578
#> sample5 0.34030847 0.6095144 0.3741498 0.3364915 0.4772101 0.8290134
#> sample6 0.45667473 0.4614949 0.8676809 0.7274044 0.9167450 0.6643710
#> # Showing 6 of 20 rows and 6 of 20 columns
```
