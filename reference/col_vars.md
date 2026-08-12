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

`NA`/`NaN` are treated as missing and dropped (equivalent to
`var(x, na.rm = TRUE)`). `Inf`/`-Inf` are **not** missing. They enter
the arithmetic, so a column's variance can be `NaN`, matching base R.

## Examples

``` r
set.seed(123)
obj <- matrix(rnorm(7 * 10), ncol = 7)
obj[1, 1] <- Inf
obj[1, 2] <- NA
obj[1:8, 3] <- NA
obj[8, 3] <- Inf
obj[1:8, 4] <- NA
obj[1:8, 5] <- NA
obj[9, 5] <- obj[10, 5]
obj[1:9, 6] <- NA
obj[, 7] <- NA
obj
#>              [,1]       [,2]      [,3]       [,4]        [,5]      [,6] [,7]
#>  [1,]         Inf         NA        NA         NA          NA        NA   NA
#>  [2,] -0.23017749  0.3598138        NA         NA          NA        NA   NA
#>  [3,]  1.55870831  0.4007715        NA         NA          NA        NA   NA
#>  [4,]  0.07050839  0.1106827        NA         NA          NA        NA   NA
#>  [5,]  0.12928774 -0.5558411        NA         NA          NA        NA   NA
#>  [6,]  1.71506499  1.7869131        NA         NA          NA        NA   NA
#>  [7,]  0.46091621  0.4978505        NA         NA          NA        NA   NA
#>  [8,] -1.26506123 -1.9666172       Inf         NA          NA        NA   NA
#>  [9,] -0.68685285  0.7013559 -1.138137 -0.3059627 -0.08336907        NA   NA
#> [10,] -0.44566197 -0.4727914  1.253815 -0.3804710 -0.08336907 0.2159416   NA

col_vars(obj)
#> [1]         NaN 1.069079431         NaN 0.002775746 0.000000000          NA
#> [7]          NA
apply(obj, 2, var, na.rm = TRUE)
#> [1]         NaN 1.069079431         NaN 0.002775746 0.000000000          NA
#> [7]          NA
```
