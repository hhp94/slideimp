# Column or Row Missing Count (or Proportion)

Calculate the number (or proportion) of missing values per column or per
row of a numeric matrix without realizing the full mask matrix.

## Usage

``` r
mat_miss(obj, col = TRUE, prop = FALSE)
```

## Arguments

- obj:

  A numeric matrix.

- col:

  Logical. If `TRUE` (default), compute per-column statistics; if
  `FALSE`, compute per-row statistics.

- prop:

  Logical. If `FALSE` (default), return raw missing counts; if `TRUE`,
  return missing proportion.

## Value

A named numeric vector of either the column or row missing count or
proportion.

## Examples

``` r
mat <- matrix(c(1, NA, 3, 4, NA, 6, NA, 8, 9), nrow = 3)
mat
#>      [,1] [,2] [,3]
#> [1,]    1    4   NA
#> [2,]   NA   NA    8
#> [3,]    3    6    9

# column missing counts (default)
mat_miss(mat)
#> [1] 1 1 1

# row missing counts
mat_miss(mat, col = FALSE)
#> [1] 1 2 0

# column missing proportions
mat_miss(mat, prop = TRUE)
#> [1] 0.3333333 0.3333333 0.3333333
```
