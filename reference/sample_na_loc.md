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
  is derived as `as.integer(num_na %/% n_rows)`). Must be provided if
  `num_na` is `NULL`. Ignored in
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  when `na_loc` is supplied.

- n_rows:

  Integer. The target number of `NA` per column (default `2L`).

  - When `num_na` is supplied: used as the base size. Most columns
    receive exactly `n_rows` `NA`; `num_na %% n_rows` columns receive
    `n_rows + 1`. If there's only one column then that column receives
    all the remainder.

  - When `num_na` is `NULL`: every selected column receives exactly
    `n_rows` `NA`. Ignored in
    [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
    when `na_loc` is supplied.

- num_na:

  Integer. Total number of missing values to inject per repetition. If
  supplied, `n_cols` is computed automatically and the `NA` are
  distributed as evenly as possible using `n_rows` as the base (`num_na`
  must be `>= n_rows`). If omitted but `n_cols` is supplied, exactly
  `n_cols * n_rows` `NA` are injected. If `num_na`, `n_cols`, and
  `na_loc` are all `NULL`, defaults to ~5% of total cells (capped
  at 500) in
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).
  `sample_na_loc()` has no default `num_na`. Ignored in
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
    the number derived from `num_na`).

  Ignored in
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

The function uses a greedy stochastic search valid `NA` locations. It
ensures that:

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
#> [1,]   3   6
#> [2,]   9   3
#> [3,]   3   2
#> [4,]   5   8
#> [5,]   3   4
#> 

# Inject the `NA` from the first repetition
mat[locs[[1]]] <- NA
mat
#>            [,1]         [,2]       [,3]       [,4]      [,5]      [,6]
#>  [1,] 0.5044793 0.0008993864 0.43429326 0.22620728 0.7797749 0.9254811
#>  [2,] 0.4542898 0.7648822351 0.34157126 0.02528094 0.3495335 0.6082241
#>  [3,] 0.5930295           NA 0.04843577         NA 0.3519479        NA
#>  [4,] 0.6524231 0.8888526005 0.74858641 0.64403813 0.5107448 0.8905952
#>  [5,] 0.1709622 0.7550757593 0.18985401 0.89009857 0.1295436 0.3780364
#>  [6,] 0.1292282 0.0657838900 0.79292257 0.44571449 0.7847482 0.5121671
#>  [7,] 0.0394463 0.1643661999 0.83293183 0.42431843 0.5675350 0.7526136
#>  [8,] 0.1044847 0.4768695564 0.69533197 0.03300803 0.8175090 0.7031305
#>  [9,] 0.3594383 0.5060374299         NA 0.71031690 0.3599687 0.2534801
#> [10,] 0.2879419 0.4133408088 0.67032239 0.66022914 0.3902879 0.2829546
#>            [,7]       [,8]       [,9]      [,10]
#>  [1,] 0.4928392 0.38979175 0.29948597 0.16996773
#>  [2,] 0.2184650 0.17013852 0.84989299 0.92530165
#>  [3,] 0.6902613 0.02169921 0.43049085 0.32640447
#>  [4,] 0.4112966 0.41161060 0.51858746 0.15314706
#>  [5,] 0.3430054         NA 0.28888987 0.92561734
#>  [6,] 0.6131221 0.67631006 0.36566163 0.21208089
#>  [7,] 0.1128234 0.55242846 0.70947389 0.43831051
#>  [8,] 0.2256672 0.84392291 0.76371036 0.03751617
#>  [9,] 0.6227623 0.04659831 0.05750697 0.15345698
#> [10,] 0.3846486 0.95434221 0.31801622 0.15895596
```
