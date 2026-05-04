# Impute Numeric Matrix with PCA Imputation

Impute missing values in a numeric matrix using (regularized) iterative
PCA.

## Usage

``` r
pca_imp(
  obj,
  ncp = 2,
  scale = TRUE,
  method = c("regularized", "EM"),
  coeff.ridge = 1,
  row.w = NULL,
  threshold = 1e-06,
  seed = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  colmax = 0.9,
  post_imp = TRUE,
  na_check = TRUE
)
```

## Arguments

- obj:

  A numeric matrix with **samples in rows** and **features in columns**.

- ncp:

  Integer. Number of components used to predict the missing entries.

- scale:

  Logical. If `TRUE` (default), variables are scaled to have unit
  variance.

- method:

  Character. Either `"regularized"` (default) or `"EM"`.

- coeff.ridge:

  Numeric. Ridge regularization coefficient (default is 1). Only used if
  `method = "regularized"`. Values \< 1 regularize less (closer to EM);
  values \> 1 regularize more (closer to mean imputation).

- row.w:

  Row weights (internally normalized to sum to 1). Can be one of:

  - `NULL` (default): All rows weighted equally.

  - A numeric vector: Custom positive weights of length `nrow(obj)`.

  - `"n_miss"`: Rows with more missing values receive lower weight.

- threshold:

  Numeric. The threshold for assessing convergence.

- seed:

  Integer. Random number generator seed.

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations for the algorithm.

- miniter:

  Integer. Minimum number of iterations for the algorithm.

- colmax:

  Numeric. A number from 0 to 1. Threshold of column-wise missing data
  rate above which imputation is skipped.

- post_imp:

  Boolean. Whether to impute remaining missing values (those that failed
  imputation) using column means.

- na_check:

  Boolean. Check for leftover `NA` values in the results or not
  (internal use).

## Value

A numeric matrix of the same dimensions as `obj` with missing values
imputed.

## Details

This algorithm is based on the original
[`missMDA::imputePCA`](https://rdrr.io/pkg/missMDA/man/imputePCA.html)
function and is optimized for tall or wide numeric matrices.

## References

Josse, J. & Husson, F. (2013). Handling missing values in exploratory
multivariate data analysis methods. Journal de la SFdS. 153 (2), pp.
79-99.

Josse, J. and Husson, F. (2016). missMDA: A Package for Handling Missing
Values in Multivariate Data Analysis. Journal of Statistical Software,
70 (1), pp 1-31.
[doi:10.18637/jss.v070.i01](https://doi.org/10.18637/jss.v070.i01)

## Author

Francois Husson and Julie Josse (original `missMDA` implementation).

## Examples

``` r
obj <- sim_mat(10, 10)$input
sum(is.na(obj))
#> [1] 10
obj[1:4, 1:4]
#>          feature1  feature2  feature3   feature4
#> sample1 0.1993174 0.0000000 0.5194403 0.19248733
#> sample2 0.5861872 0.7793109 1.0000000 0.05740143
#> sample3 0.2991663 0.4695521 0.8906237 0.53696051
#> sample4 0.3010367 0.3079201 0.0000000 0.00000000
# Randomly initialize missing values 5 times (1st time is mean).
pca_imp(obj, ncp = 2, nb.init = 5)
#> Method: PCA imputation
#> Dimensions: 10 x 10
#> 
#>          feature1  feature2  feature3   feature4  feature5  feature6
#> sample1 0.1993174 0.0000000 0.5194403 0.19248733 0.0000000 0.3136988
#> sample2 0.5861872 0.7793109 1.0000000 0.05740143 0.1396488 0.3300150
#> sample3 0.2991663 0.4695521 0.8906237 0.53696051 0.6138600 0.3213725
#> sample4 0.3010367 0.3079201 0.0000000 0.00000000 0.4679464 0.4484497
#> sample5 0.3898840 0.8718041 0.8782409 0.81735982 1.0000000 0.7495673
#> sample6 0.3191516 0.6561425 0.2994724 0.89020793 0.4422327 0.6983652
#> 
#> # Showing [1:6, 1:6] of full matrix
```
