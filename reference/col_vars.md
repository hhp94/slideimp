# Calculate Matrix Column Variances

Computes the sample variance for each column of a numeric matrix.

## Usage

``` r
col_vars(mat, cores = 1)
```

## Arguments

- mat:

  A numeric matrix.

- cores:

  Integer. Number of cores to use for parallel computation. Defaults to
  `1`.

## Value

A numeric vector of column variances, named when `mat` has column names.

## Details

Columns with fewer than two distinct non-missing values are assigned
`NA`.

## Examples

``` r
set.seed(123)
mat <- matrix(rnorm(4 * 10), ncol = 4)
mat[1, 1] <- NA
mat[1:8, 2] <- NA
mat[1:9, 3] <- NA
mat[, 4] <- NA
mat
#>              [,1]       [,2]     [,3] [,4]
#>  [1,]          NA         NA       NA   NA
#>  [2,] -0.23017749         NA       NA   NA
#>  [3,]  1.55870831         NA       NA   NA
#>  [4,]  0.07050839         NA       NA   NA
#>  [5,]  0.12928774         NA       NA   NA
#>  [6,]  1.71506499         NA       NA   NA
#>  [7,]  0.46091621         NA       NA   NA
#>  [8,] -1.26506123         NA       NA   NA
#>  [9,] -0.68685285  0.7013559       NA   NA
#> [10,] -0.44566197 -0.4727914 1.253815   NA

col_vars(mat)
#> [1] 0.9673957 0.6893110        NA        NA
apply(mat, 2, var, na.rm = TRUE)
#> [1] 0.9673957 0.6893110        NA        NA
```
