# Calculate Matrix Column Variances

Compute the sample variance for each column of a numeric matrix.

## Usage

``` r
col_vars(obj, cores = 1)
```

## Arguments

- obj:

  A numeric matrix.

- cores:

  Integer. Number of cores to use for parallel computation. Defaults to
  `1`.

## Value

A numeric vector of column variances, named if `obj` has column names.

## Details

Columns with fewer than two non-missing values are assigned `NA`.

## Examples

``` r
set.seed(123)
obj <- matrix(rnorm(5 * 10), ncol = 5)
obj[1, 1] <- NA
obj[1:8, 2] <- NA
obj[1:8, 3] <- NA
obj[9, 3] <- obj[10, 3]
obj[1:9, 4] <- NA
obj[, 5] <- NA
obj
#>              [,1]       [,2]     [,3]      [,4] [,5]
#>  [1,]          NA         NA       NA        NA   NA
#>  [2,] -0.23017749         NA       NA        NA   NA
#>  [3,]  1.55870831         NA       NA        NA   NA
#>  [4,]  0.07050839         NA       NA        NA   NA
#>  [5,]  0.12928774         NA       NA        NA   NA
#>  [6,]  1.71506499         NA       NA        NA   NA
#>  [7,]  0.46091621         NA       NA        NA   NA
#>  [8,] -1.26506123         NA       NA        NA   NA
#>  [9,] -0.68685285  0.7013559 1.253815        NA   NA
#> [10,] -0.44566197 -0.4727914 1.253815 -0.380471   NA

col_vars(obj)
#> [1] 0.9673957 0.6893110 0.0000000        NA        NA
apply(obj, 2, var, na.rm = TRUE)
#> [1] 0.9673957 0.6893110 0.0000000        NA        NA
```
