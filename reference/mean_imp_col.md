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

## Examples

``` r
mat <- matrix(c(1, 2, NA, 4, NA, 6, NA, 8, 9), nrow = 3)
colnames(mat) <- c("A", "B", "C")
mat
#>       A  B  C
#> [1,]  1  4 NA
#> [2,]  2 NA  8
#> [3,] NA  6  9

# Impute missing values with column means
imputed_mat <- mean_imp_col(mat)
imputed_mat
#>        A B   C
#> [1,] 1.0 4 8.5
#> [2,] 2.0 5 8.0
#> [3,] 1.5 6 9.0

# Impute only specific columns by name
imputed_subset <- mean_imp_col(mat, subset = c("A", "C"))
imputed_subset
#>        A  B   C
#> [1,] 1.0  4 8.5
#> [2,] 2.0 NA 8.0
#> [3,] 1.5  6 9.0

# Impute only specific columns by index
imputed_idx <- mean_imp_col(mat, subset = c(1, 3))
imputed_idx
#>        A  B   C
#> [1,] 1.0  4 8.5
#> [2,] 2.0 NA 8.0
#> [3,] 1.5  6 9.0
```
