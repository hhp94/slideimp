# Sample Missing Value Locations with Constraints

Samples indices for `NA` injection into a matrix while maintaining
row/column missing value budgets and avoiding zero-variance columns.

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

  A numeric matrix with **samples in rows** and **features in columns**.

- n_cols:

  Integer. The number of columns to receive injected `NA` per
  repetition. Ignored when `num_na` is supplied (in which case `n_cols`
  is derived as `num_na %/% n_rows`). Must be provided if `num_na` is
  `NULL`. Ignored in
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  when `na_loc` is supplied.

- n_rows:

  Integer. The target number of `NA` values per column (default `2L`).

  - When `num_na` is supplied: used as the base size. Most columns
    receive exactly `n_rows` missing values; `num_na %% n_rows` columns
    receive one extra. If there's only one column, it receives all the
    remainder.

  - When `num_na` is `NULL`: every selected column receives exactly
    `n_rows` `NA`. Ignored in
    [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
    when `na_loc` is supplied.

- num_na:

  Integer. Total number of missing values to inject per repetition. If
  supplied, `n_cols` is computed automatically and missing values are
  distributed as evenly as possible, using `n_rows` as the base size
  (`num_na` must be at least `n_rows`). If omitted but `n_cols` is
  supplied, the total injected is `n_cols * n_rows`. If `num_na`,
  `n_cols`, and `na_loc` are all `NULL`,
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  defaults to roughly 5% of total cells, capped at 500.
  `sample_na_loc()` has no default. Ignored in
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  when `na_loc` is supplied.

- n_reps:

  Integer. Number of repetitions for random NA injection (default `1`).

- rowmax, colmax:

  Numbers between 0 and 1. NA injection cannot create rows/columns with
  a higher proportion of missing values than these thresholds.

- na_col_subset:

  Optional integer or character vector restricting which columns of
  `obj` are eligible for NA injection.

  - If `NULL` (default): all columns are eligible.

  - If character: values must exist in `colnames(obj)`.

  - If integer/numeric: values must be valid 1-based column indices. The
    vector must be unique and must contain at least `n_cols` columns (or
    the number derived from `num_na`). Ignored in
    [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
    when `na_loc` is supplied.

- max_attempts:

  Integer. Maximum number of resampling attempts per repetition before
  giving up due to row-budget exhaustion (default `100`).

## Value

A list of length `n_reps`. Each element is a two-column integer matrix
(`row`, `col`) representing the coordinates of the sampled `NA`
locations.

## Details

The function uses a greedy stochastic search for valid `NA` locations.
It ensures that:

- Total missingness per row and column does not exceed `rowmax` and
  `colmax`.

- At least two distinct observed values are preserved in every column to
  ensure the column maintains non-zero variance.

## Examples

``` r
mat <- matrix(runif(100), nrow = 10)

# Sample 5 `NA` across 5 columns (1 per column)
locs <- sample_na_loc(mat, n_cols = 5, n_rows = 1)
locs
#> [[1]]
#>      row col
#> [1,]   2   4
#> [2,]   9   1
#> [3,]   7   5
#> [4,]   8   3
#> [5,]   5   7
#> 

# Inject the `NA` from the first repetition
mat[locs[[1]]] <- NA
mat
#>            [,1]      [,2]      [,3]        [,4]       [,5]      [,6]       [,7]
#>  [1,] 0.8376339 0.4038305 0.1453584 0.191754386 0.92877424 0.1981642 0.57884620
#>  [2,] 0.4874663 0.8675411 0.4213214          NA 0.26107766 0.8810928 0.54187464
#>  [3,] 0.1103370 0.0633672 0.2007586 0.003423417 0.44312295 0.3643178 0.76257168
#>  [4,] 0.3513630 0.6840599 0.9788740 0.931021894 0.70764482 0.2699500 0.14225535
#>  [5,] 0.7610630 0.6885818 0.4949632 0.689435125 0.90294457 0.5899745         NA
#>  [6,] 0.3896648 0.3139780 0.3744674 0.252752172 0.51312311 0.1654506 0.80370827
#>  [7,] 0.4565140 0.6385193 0.8503906 0.665434266         NA 0.4890698 0.94901189
#>  [8,] 0.1073443 0.6412363        NA 0.476809311 0.04979171 0.9633836 0.03845219
#>  [9,]        NA 0.4929685 0.3577581 0.904765149 0.24939016 0.5796029 0.66919590
#> [10,] 0.9880969 0.8560004 0.3736045 0.305509669 0.78848105 0.7907979 0.16308326
#>            [,8]       [,9]      [,10]
#>  [1,] 0.1163298 0.66799949 0.12984316
#>  [2,] 0.1324231 0.99507695 0.11584997
#>  [3,] 0.3910245 0.90379414 0.03223105
#>  [4,] 0.2157433 0.34029998 0.28859893
#>  [5,] 0.6970320 0.25729868 0.57751060
#>  [6,] 0.1571268 0.32218851 0.03855550
#>  [7,] 0.3405160 0.49326583 0.72626821
#>  [8,] 0.3331250 0.01892372 0.89637255
#>  [9,] 0.6349510 0.28512190 0.73863998
#> [10,] 0.4057082 0.65870852 0.40522444
```
