# Calculate Matrix Column Variance

`col_vars` computes the sample variance for each column of a numeric
matrix.

## Usage

``` r
col_vars(mat, cores = 1)
```

## Arguments

- mat:

  A numeric matrix.

- cores:

  Number of cores to use for parallel computation. Defaults to 1.

## Value

`col_vars` returns a named numeric vector of column variances.

## Details

Variances for columns with one unique value after dropping `NA` are set
to `NA`.

## Examples

``` r
mat <- matrix(rnorm(4 * 10), ncol = 4)
mat[1, 1] <- NA
mat[1:8, 2] <- NA
mat[1:9, 3] <- NA
mat[, 4] <- NA
mat
#>               [,1]       [,2]    [,3] [,4]
#>  [1,]           NA         NA      NA   NA
#>  [2,]  0.255317055         NA      NA   NA
#>  [3,] -2.437263611         NA      NA   NA
#>  [4,] -0.005571287         NA      NA   NA
#>  [5,]  0.621552721         NA      NA   NA
#>  [6,]  1.148411606         NA      NA   NA
#>  [7,] -1.821817661         NA      NA   NA
#>  [8,] -0.247325302         NA      NA   NA
#>  [9,] -0.244199607  0.5429963      NA   NA
#> [10,] -0.282705449 -0.9140748 -1.5124   NA
col_vars(mat)
#> [1] 1.277663 1.061528       NA       NA
apply(mat, 2, var, na.rm = TRUE)
#> [1] 1.277663 1.061528       NA       NA
```
