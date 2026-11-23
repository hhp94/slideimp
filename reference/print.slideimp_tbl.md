# Print Method for `slideimp_tbl` Objects

Print a `slideimp_tbl` data frame. Each cell of a list-column shows as
its type and size.

## Usage

``` r
# S3 method for class 'slideimp_tbl'
print(x, n = NULL, ...)
```

## Arguments

- x:

  A `slideimp_tbl` object.

- n:

  Integer or `NULL`. Number of rows to show, 10 if `NULL`.

- ...:

  Further arguments passed to
  [`print()`](https://rdrr.io/r/base/print.html) for the table.

## Value

`x`, invisibly.

## Examples

``` r
sim <- sim_mat(n = 10, p = 20)
tbl <- prep_groups(colnames(sim$input), sim$col_group)
class(tbl)
#> [1] "slideimp_tbl" "data.frame"  
print(tbl)
#> # slideimp table: 2 x 4
#>   group          feature             aux parameters
#>  group1 <character [14]> <character [0]> <list [0]>
#>  group2  <character [6]> <character [0]> <list [0]>
```
