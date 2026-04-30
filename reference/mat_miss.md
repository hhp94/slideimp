# Column or Row Missing Counts and Proportions

Calculate the number or proportion of missing values per column or per
row of a numeric matrix without allocating a full logical mask matrix.

## Usage

``` r
mat_miss(obj, col = TRUE, prop = FALSE)
```

## Arguments

- obj:

  A numeric matrix.

- col:

  Logical. If `TRUE`, compute column-wise. If `FALSE`, compute row-wise.

- prop:

  Logical. If `FALSE`, return missing-value counts. If `TRUE`, return
  missing-value proportions.

## Value

A numeric vector containing missing-value counts or proportions for
columns or rows, named when the corresponding dimension names are
present.

## Examples

``` r
obj <- matrix(c(1, NA, 3, 4, NA, 6, NA, 8, 9), nrow = 3)
obj
#>      [,1] [,2] [,3]
#> [1,]    1    4   NA
#> [2,]   NA   NA    8
#> [3,]    3    6    9

# column missing counts
mat_miss(obj)
#> [1] 1 1 1

# row missing counts
mat_miss(obj, col = FALSE)
#> [1] 1 2 0

# column missing proportions
mat_miss(obj, prop = TRUE)
#> [1] 0.3333333 0.3333333 0.3333333
```
