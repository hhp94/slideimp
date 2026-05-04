# Sample Missing-Value Locations with Constraints

Sample matrix indices for `NA` injection while respecting row and column
missingness limits and avoiding zero-variance columns.

## Usage

``` r
sample_na_loc(
  obj,
  n_cols = NULL,
  n_rows = 2L,
  num_na = NULL,
  n_reps = 1L,
  rowmax = 0.9,
  colmax = 0.9,
  na_col_subset = NULL,
  max_attempts = 100
)
```

## Arguments

- obj:

  A numeric matrix.

- n_cols:

  Integer or `NULL`. Number of columns to receive injected missing
  values. Must be supplied when `num_na = NULL`.

- n_rows:

  Integer. Target number of missing values to inject per selected
  column.

- num_na:

  Integer or `NULL`. Total number of missing values to inject per
  repetition. If supplied, `n_cols` is derived from `num_na` and
  `n_rows`, and missing values are distributed as evenly as possible
  across columns.

- n_reps:

  Integer. Number of independent repetitions.

- rowmax:

  Numeric scalar between `0` and `1`. Maximum allowed missing-data
  proportion per row after injection.

- colmax:

  Numeric scalar between `0` and `1`. Maximum allowed missing-data
  proportion per column after injection.

- na_col_subset:

  Optional integer or character vector restricting which columns are
  eligible for missing-value injection.

- max_attempts:

  Integer. Maximum number of resampling attempts per repetition before
  giving up.

## Value

A list of length `n_reps`. Each element is a two-column integer matrix
with row and column indices for sampled `NA` locations.

## Details

The function uses a greedy stochastic search for valid `NA` locations.
It ensures that:

- Total missingness per row and column does not exceed `rowmax` and
  `colmax`.

- At least two distinct observed values are preserved in every affected
  column.

## Examples

``` r
set.seed(123)
mat <- matrix(runif(100), nrow = 10)

# Sample 5 missing values across 5 columns
locs <- sample_na_loc(mat, n_cols = 5, n_rows = 1)
locs
#> [[1]]
#>      row col
#> [1,]   1  10
#> [2,]   2   2
#> [3,]   5   9
#> [4,]   3   5
#> [5,]   5   3
#> 

# Inject the missing values from the first repetition
mat[locs[[1]]] <- NA
mat
#>            [,1]       [,2]      [,3]       [,4]      [,5]       [,6]       [,7]
#>  [1,] 0.2875775 0.95683335 0.8895393 0.96302423 0.1428000 0.04583117 0.66511519
#>  [2,] 0.7883051         NA 0.6928034 0.90229905 0.4145463 0.44220007 0.09484066
#>  [3,] 0.4089769 0.67757064 0.6405068 0.69070528        NA 0.79892485 0.38396964
#>  [4,] 0.8830174 0.57263340 0.9942698 0.79546742 0.3688455 0.12189926 0.27438364
#>  [5,] 0.9404673 0.10292468        NA 0.02461368 0.1524447 0.56094798 0.81464004
#>  [6,] 0.0455565 0.89982497 0.7085305 0.47779597 0.1388061 0.20653139 0.44851634
#>  [7,] 0.5281055 0.24608773 0.5440660 0.75845954 0.2330341 0.12753165 0.81006435
#>  [8,] 0.8924190 0.04205953 0.5941420 0.21640794 0.4659625 0.75330786 0.81238951
#>  [9,] 0.5514350 0.32792072 0.2891597 0.31818101 0.2659726 0.89504536 0.79434232
#> [10,] 0.4566147 0.95450365 0.1471136 0.23162579 0.8578277 0.37446278 0.43983169
#>               [,8]      [,9]      [,10]
#>  [1,] 0.7544751586 0.2436195         NA
#>  [2,] 0.6292211316 0.6680556 0.65310193
#>  [3,] 0.7101824014 0.4176468 0.34351647
#>  [4,] 0.0006247733 0.7881958 0.65675813
#>  [5,] 0.4753165741        NA 0.32037324
#>  [6,] 0.2201188852 0.4348927 0.18769112
#>  [7,] 0.3798165377 0.9849570 0.78229430
#>  [8,] 0.6127710033 0.8930511 0.09359499
#>  [9,] 0.3517979092 0.8864691 0.46677904
#> [10,] 0.1111354243 0.1750527 0.51150546
```
