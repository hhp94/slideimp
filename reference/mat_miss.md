# Column or Row Missing Counts and Proportions

Calculate the number or proportion of missing values per column or per
row of a numeric matrix, without allocating a logical matrix the size of
`obj`.

## Usage

``` r
mat_miss(obj, col = TRUE, prop = FALSE)
```

## Arguments

- obj:

  A double matrix.

- col:

  Logical. Whether to compute per column (`TRUE`) or per row (`FALSE`).

- prop:

  Logical. Whether to return missing-value proportions (`TRUE`) or
  counts (`FALSE`).

## Value

When `prop = FALSE`, an integer vector of missing-value counts. When
`prop = TRUE`, a double vector of missing-value proportions. The vector
has one element per column when `col = TRUE`, or per row otherwise,
named from `colnames(obj)` or `rownames(obj)` when set.

## Details

`NA` and `NaN` both count as missing. `Inf` and `-Inf` do not.

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
