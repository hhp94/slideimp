# Column Mean Imputation

Impute missing values in a matrix by replacing them with the mean of
their respective columns.

## Usage

``` r
mean_imp_col(obj, subset = NULL, cores = 1)
```

## Arguments

- obj:

  A numeric matrix.

- subset:

  Optional character or integer vector specifying columns to impute. If
  `NULL`, all columns are imputed.

- cores:

  Integer. Number of cores to use for parallel computation. Defaults to
  `1`.

## Value

A numeric matrix of the same dimensions as `obj`, with missing values in
the selected columns replaced by column means.

## Details

Columns with no observed values cannot be imputed by their column mean
and are left unchanged.

## Examples

``` r
obj <- matrix(c(1, 2, NA, 4, NA, 6, NA, 8, 9, NA, NA, NA), nrow = 3)
colnames(obj) <- c("A", "B", "C", "D")
obj
#>       A  B  C  D
#> [1,]  1  4 NA NA
#> [2,]  2 NA  8 NA
#> [3,] NA  6  9 NA

# impute missing values with column means
mean_imp_col(obj)
#>        A B   C  D
#> [1,] 1.0 4 8.5 NA
#> [2,] 2.0 5 8.0 NA
#> [3,] 1.5 6 9.0 NA

# impute only specific columns by name
mean_imp_col(obj, subset = c("A", "C"))
#>        A  B   C  D
#> [1,] 1.0  4 8.5 NA
#> [2,] 2.0 NA 8.0 NA
#> [3,] 1.5  6 9.0 NA

# impute only specific columns by index
mean_imp_col(obj, subset = c(1, 3))
#>        A  B   C  D
#> [1,] 1.0  4 8.5 NA
#> [2,] 2.0 NA 8.0 NA
#> [3,] 1.5  6 9.0 NA
```
