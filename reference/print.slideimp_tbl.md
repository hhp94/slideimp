# Print a `slideimp_tbl` Object

Print `slideimp_tbl` objects (which inherit `data.frame`) with nicer
looking list-columns (similar to `tibble`).

## Usage

``` r
# S3 method for class 'slideimp_tbl'
print(x, n = NULL, ...)
```

## Arguments

- x:

  A `slideimp_tbl` object.

- n:

  Number of rows to show. Defaults to 10.

- ...:

  Not used.

## Value

Invisible `x`.

## Examples

``` r
mat <- sim_mat(n = 10, p = 500)
set.seed(1234)
results <- tune_imp(mat$input, parameters = data.frame(k = 5), .f = "knn_imp")
#> Tuning knn_imp
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 250 (~5% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running Mode: sequential...
#> Step 2/2: Tuning
class(results)
#> [1] "slideimp_tune" "slideimp_tbl"  "data.frame"   
print(results)
#> # slideimp table: 1 x 5
#>  k param_set rep_id         result error
#>  5         1      1 <df [250 x 2]>  <NA>
```
