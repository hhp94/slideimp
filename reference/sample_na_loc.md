# Sample Missing-Value Locations with Constraints

Draws matrix positions for `NA` injection. The positions respect row and
column limits on missing values and skip constant columns.

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

  A double matrix.

- n_cols:

  Integer or `NULL`. Number of columns that receive missing values.

- n_rows:

  Integer. Minimum number of missing values per selected column.

- num_na:

  Integer or `NULL`. Total number of missing values per repetition. When
  given, `n_cols` is ignored and the values are spread as evenly as
  possible over `num_na %/% n_rows` columns.

- n_reps:

  Integer. Number of independent repetitions.

- rowmax:

  Numeric scalar between `0` and `1`. Maximum proportion of missing
  values per row after injection.

- colmax:

  Numeric scalar between `0` and `1`. Maximum proportion of missing
  values per column after injection.

- na_col_subset:

  Optional character or integer vector. The columns that can receive
  missing values, by name or by number.

- max_attempts:

  Integer. Maximum number of sampling attempts per repetition before an
  error.

## Value

A list of length `n_reps`. Each element is a two-column integer matrix,
with columns `row` and `col`, of the sampled `NA` locations.

## Details

`sample_na_loc()` places the `NA` locations by a greedy random search.
Every result meets these limits:

- The proportion of missing values in each row and column, existing ones
  included, stays at or under `rowmax` and `colmax`.

- Every selected column keeps at least two distinct observed values.

- A column with a variance below `.Machine$double.eps` counts as
  constant and is never selected.
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  uses the same threshold to leave a column out of its fit.

A column that cannot take `n_rows` missing values under these limits is
skipped without an error.

When both `n_cols` and `num_na` are `NULL`, `num_na` defaults to 5% of
the cells of `obj`, with a maximum of 500. The default never asks for
more columns than the number that can take `n_rows` missing values.

## See also

- [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md),
  which takes the result as `na_loc`.

- [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md),
  which takes a result with `n_reps = 1` as `na_loc`.

## Examples

``` r
set.seed(123)
mat <- matrix(runif(100), nrow = 10)

# Sample 5 missing values across 5 columns
locs <- sample_na_loc(mat, n_cols = 5, n_rows = 1)
locs
#> [[1]]
#>      row col
#> [1,]   8   6
#> [2,]   5   1
#> [3,]   3   2
#> [4,]   2   4
#> [5,]   4   8
#> 

# Inject the missing values from the first repetition
mat[locs[[1]]] <- NA
mat
#>            [,1]       [,2]      [,3]       [,4]      [,5]       [,6]       [,7]
#>  [1,] 0.2875775 0.95683335 0.8895393 0.96302423 0.1428000 0.04583117 0.66511519
#>  [2,] 0.7883051 0.45333416 0.6928034         NA 0.4145463 0.44220007 0.09484066
#>  [3,] 0.4089769         NA 0.6405068 0.69070528 0.4137243 0.79892485 0.38396964
#>  [4,] 0.8830174 0.57263340 0.9942698 0.79546742 0.3688455 0.12189926 0.27438364
#>  [5,]        NA 0.10292468 0.6557058 0.02461368 0.1524447 0.56094798 0.81464004
#>  [6,] 0.0455565 0.89982497 0.7085305 0.47779597 0.1388061 0.20653139 0.44851634
#>  [7,] 0.5281055 0.24608773 0.5440660 0.75845954 0.2330341 0.12753165 0.81006435
#>  [8,] 0.8924190 0.04205953 0.5941420 0.21640794 0.4659625         NA 0.81238951
#>  [9,] 0.5514350 0.32792072 0.2891597 0.31818101 0.2659726 0.89504536 0.79434232
#> [10,] 0.4566147 0.95450365 0.1471136 0.23162579 0.8578277 0.37446278 0.43983169
#>            [,8]      [,9]      [,10]
#>  [1,] 0.7544752 0.2436195 0.13069569
#>  [2,] 0.6292211 0.6680556 0.65310193
#>  [3,] 0.7101824 0.4176468 0.34351647
#>  [4,]        NA 0.7881958 0.65675813
#>  [5,] 0.4753166 0.1028646 0.32037324
#>  [6,] 0.2201189 0.4348927 0.18769112
#>  [7,] 0.3798165 0.9849570 0.78229430
#>  [8,] 0.6127710 0.8930511 0.09359499
#>  [9,] 0.3517979 0.8864691 0.46677904
#> [10,] 0.1111354 0.1750527 0.51150546
```
