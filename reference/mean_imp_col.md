# Column Mean Imputation

Impute missing values in a matrix by replacing them with the mean of
their respective columns.

## Usage

``` r
mean_imp_col(obj, subset = NULL, cores = 1)
```

## Arguments

- obj:

  A numeric matrix with **samples in rows** and **features in columns**.

- subset:

  Character. Vector of column names or integer vector of column indices
  specifying which columns to impute.

- cores:

  Integer. Number of cores for K-NN parallelization (OpenMP). On macOS,
  OpenMP may need additional compiler configuration.

## Value

A numeric matrix of the same dimensions as `obj` with missing values in
the specified columns replaced by column means.

## Examples

``` r
# Create example matrix with missing values
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
