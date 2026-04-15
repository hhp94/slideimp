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
  miniter = 5
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

## Value

A numeric matrix of the same dimensions as `obj` with missing values
imputed.

## Details

This algorithm is based on the original
[`missMDA::imputePCA`](https://rdrr.io/pkg/missMDA/man/imputePCA.html)
function and is optimized for tall/wide numeric matrices.

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
data("khanmiss1")

# Transpose to put genes on columns.
# Randomly initialize missing values 5 times (1st time is mean).
pca_imp(t(khanmiss1), ncp = 2, nb.init = 5)
#> slideimp_results (PCA)
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
