# Print a `slideimp_results` Object

Print a `slideimp_results` Object

## Usage

``` r
# S3 method for class 'slideimp_results'
print(x, n = 6L, p = 6L, ...)
```

## Arguments

- x:

  An `slideimp_results` object

- n:

  Number of rows to print

- p:

  Number of cols to print

- ...:

  Not used

## Value

Invisible `x`

## Examples

``` r
data(khanmiss1)
t_khanmiss1 <- t(khanmiss1)
result <- knn_imp(t_khanmiss1, k = 5)
class(result)
#> [1] "slideimp_results" "matrix"           "array"           
print(result, n = 6, p = 6)
#> slideimp_results (KNN)
#> Dimensions: 63 x 2308
#> 
#>           g1   g2   g3   g4   g5   g6
#> sample1 1873 1251  314 1324  776 1901
#> sample2   57 1350 1758 1428  476 1521
#> sample3   53 1140  162 1468  679   14
#> sample4 2059 1385 1857 1250  772 2052
#> sample5 1537 1261 1939 1666 1307 1705
#> sample6 1819 1526 1640 1795 1454 1643
#> 
#> # Showing [1:6, 1:6] of full matrix
```
