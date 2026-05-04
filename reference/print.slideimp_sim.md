# Print a `slideimp_sim` Object

Print the output of
[`sim_mat()`](https://hhp94.github.io/slideimp/reference/sim_mat.md).

## Usage

``` r
# S3 method for class 'slideimp_sim'
print(x, n = 6L, p = 6L, ...)
```

## Arguments

- x:

  A `slideimp_sim` object.

- n:

  Number of rows of each component to show.

- p:

  Number of columns of `input` to show.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(123)
sim_data <- sim_mat(n = 50, p = 10, rho = 0.5)
class(sim_data)
#> [1] "slideimp_sim"
print(sim_data)
#> $col_group (2 column groups)
#>    feature  group
#> 1 feature1 group1
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group1
#> 5 feature5 group2
#> 6 feature6 group2
#> # Showing 6 of 10 rows
#> 
#> $row_group (1 row groups)
#>    sample  group
#> 1 sample1 group1
#> 2 sample2 group1
#> 3 sample3 group1
#> 4 sample4 group1
#> 5 sample5 group1
#> 6 sample6 group1
#> # Showing 6 of 50 rows
#> 
#> $input (50 x 10)
#>          feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.3855027 0.2859279        NA 0.7034804 0.3949006 0.3390209
#> sample2 0.3939132 0.5091468 0.5414625 0.6130915 0.4153076 0.3899470
#> sample3 0.7020675 0.7302534 0.7618661 0.6474409 0.6996583 0.6687105
#> sample4 0.6887438 0.4568955 0.3007332 0.5369357 0.5503469 0.3900968
#> sample5 0.4220865 0.3630904 0.4552229        NA 0.7723463 0.5073267
#> sample6 1.0000000 0.7918423 0.6874920 0.6385426 0.7579935 0.9167009
#> # Showing 6 of 50 rows and 6 of 10 columns
```
