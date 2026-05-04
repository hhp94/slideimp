# Print a `slideimp_tbl` Object

Print `slideimp_tbl` objects, which inherit from `data.frame`, with
compact display of list-columns.

## Usage

``` r
# S3 method for class 'slideimp_tbl'
print(x, n = NULL, ...)
```

## Arguments

- x:

  A `slideimp_tbl` object.

- n:

  Number of rows to show. If `NULL`, a default is used.

- ...:

  Not used.

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
