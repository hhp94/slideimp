# Print a `slideimp_results` Object

Print the output of
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md),
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md),
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md).

## Usage

``` r
# S3 method for class 'slideimp_results'
print(x, n = 6L, p = 6L, ...)
```

## Arguments

- x:

  A `slideimp_results` object.

- n:

  Number of rows to print.

- p:

  Number of cols to print.

- ...:

  Not used.

## Value

Invisible `x`.

## Examples

``` r
set.seed(1234)
mat <- sim_mat(n = 10, p = 10)
result <- knn_imp(mat$input, k = 5)
class(result)
#> [1] "slideimp_results" "matrix"           "array"           
print(result, n = 6, p = 6)
#> Method: KNN imputation
#> Dimensions: 10 x 10
#> 
#>          feature1   feature2  feature3  feature4   feature5  feature6
#> sample1 0.1568098 0.47991216 0.8510938 1.0000000 0.07839046 0.5527338
#> sample2 0.4098416 0.66120576 0.8221066 0.6396885 0.68926345 1.0000000
#> sample3 0.6801685 1.00000000 1.0000000 0.9953450 0.75246030 0.6958550
#> sample4 0.0000000 0.01200625 0.0000000 0.0000000 0.00000000 0.0000000
#> sample5 0.9639671 0.64091372 0.5111760 0.7184678 0.81815288 0.5883255
#> sample6 0.7031741 0.37310619 0.6782811 0.7542872 0.99910554 0.9070005
#> 
#> # Showing [1:6, 1:6] of full matrix
```
