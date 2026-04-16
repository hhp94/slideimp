# Tune Parameters for Imputation Methods

Tunes hyperparameters for imputation methods such as
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md),
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md), or
user-supplied custom functions by repeated cross-validation. For
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md),
tune
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md) or
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md) on
a single group.

## Usage

``` r
tune_imp(
  obj,
  parameters = NULL,
  .f,
  na_loc = NULL,
  num_na = NULL,
  n_reps = 1,
  n_cols = NULL,
  n_rows = 2,
  rowmax = 0.9,
  colmax = 0.9,
  na_col_subset = NULL,
  max_attempts = 100,
  .progress = TRUE,
  cores = 1,
  location = NULL,
  pin_blas = FALSE
)
```

## Arguments

- obj:

  A numeric matrix with **samples in rows** and **features in columns**.

- parameters:

  A data.frame specifying parameter combinations to tune, where each
  column represents a parameter accepted by `.f` (excluding `obj`). List
  columns are supported for complex parameters. Duplicate rows are
  automatically removed. `NULL` is treated as tuning the function with
  its default parameters.

- .f:

  Either `"knn_imp"`, `"pca_imp"`, `"slide_imp"`, or a custom function
  specifying the imputation method to tune.

- na_loc:

  Optional. Pre-defined missing value locations to bypass random NA
  injection with
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md).
  Accepted formats include:

  - A two-column integer matrix (row, column indices).

  - A numeric vector of linear locations.

  - A list where each element is one of the above formats (one per
    repetition).

- num_na:

  Integer. Total number of missing values to inject per repetition. If
  supplied, `n_cols` is computed automatically and missing values are
  distributed as evenly as possible, using `n_rows` as the base size
  (`num_na` must be at least `n_rows`). If omitted but `n_cols` is
  supplied, the total injected is `n_cols * n_rows`. If `num_na`,
  `n_cols`, and `na_loc` are all `NULL`, `tune_imp()` defaults to
  roughly 5% of total cells, capped at 500.
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  has no default. Ignored in `tune_imp()` when `na_loc` is supplied.

- n_reps:

  Integer. Number of repetitions for random NA injection (default `1`).

- n_cols:

  Integer. The number of columns to receive injected `NA` per
  repetition. Ignored when `num_na` is supplied (in which case `n_cols`
  is derived as `num_na %/% n_rows`). Must be provided if `num_na` is
  `NULL`. Ignored in `tune_imp()` when `na_loc` is supplied.

- n_rows:

  Integer. The target number of `NA` values per column (default `2L`).

  - When `num_na` is supplied: used as the base size. Most columns
    receive exactly `n_rows` missing values; `num_na %% n_rows` columns
    receive one extra. If there's only one column, it receives all the
    remainder.

  - When `num_na` is `NULL`: every selected column receives exactly
    `n_rows` `NA`. Ignored in `tune_imp()` when `na_loc` is supplied.

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
    the number derived from `num_na`). Ignored in `tune_imp()` when
    `na_loc` is supplied.

- max_attempts:

  Integer. Maximum number of resampling attempts per repetition before
  giving up due to row-budget exhaustion (default `100`).

- .progress:

  Logical. Show a progress bar during tuning (default `TRUE`).

- cores:

  Controls the number of cores to parallelize over for K-NN and
  sliding-window K-NN imputation with OpenMP. For other methods, use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead.

- location:

  Required only for `slide_imp`. Numeric vector of column locations.

- pin_blas:

  Logical. Pin BLAS threads to 1 during parallel tuning (default
  `FALSE`).

## Value

A `data.frame` of class `slideimp_tune` containing:

- `...`: All columns originally provided in `parameters`.

- `param_set`: An integer ID for the unique parameter combination.

- `rep_id`: An integer indicating the repetition index.

- `result`: A nested list-column where each element is a data.frame
  containing `truth` (original values) and `estimate` (imputed values).

- `error`: A character column containing the error message if the
  iteration failed, otherwise `NA`.

## Details

The function supports tuning for built-in methods (`"slide_imp"`,
`"knn_imp"`, `"pca_imp"`) or custom functions provided via `.f`.

When `.f` is a character string, the columns in `parameters` are
validated against the chosen method's requirements:

- `"knn_imp"`: requires `k` in `parameters`

- `"pca_imp"`: requires `ncp` in `parameters`

- `"slide_imp"`: requires `window_size`, `overlap_size`, and
  `min_window_n`, plus exactly one of `k` or `ncp`

When `.f` is a custom function, the columns in `parameters` must
correspond to the arguments of `.f` (excluding the `obj` argument). The
custom function must accept `obj` (a numeric matrix) as its first
argument and return a numeric matrix of identical dimensions.

Tuning results can be evaluated using the `yardstick` package or
[`compute_metrics()`](https://hhp94.github.io/slideimp/reference/compute_metrics.md).

## Parallelization

- **K-NN**: use the `cores` argument (requires OpenMP). If
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  are active, `cores` is automatically set to 1 to avoid nested
  parallelism.

- **PCA**: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

On macOS, OpenMP is typically unavailable and `cores` falls back to

1.  Use
    [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
    for parallelization instead.

On Linux with OpenBLAS or MKL, set `pin_blas = TRUE` when running
parallel PCA to prevent BLAS threads and `mirai` workers competing for
cores.

## Examples

``` r
# Setup example data. Increase `num_na` (500) and `n_reps` (10-30) in real
# analyses
obj <- sim_mat(10, 50)$input

# 1. Tune K-NN imputation with random NA injection
params_knn <- data.frame(k = c(2, 4))
results <- tune_imp(obj, params_knn, .f = "knn_imp", n_reps = 1, num_na = 10)
#> Tuning knn_imp
#> Step 1/2: Resolving NA locations
#> Running Mode: sequential...
#> Step 2/2: Tuning
compute_metrics(results)
#>   k param_set rep_id error  n n_miss .metric .estimator .estimate
#> 1 2         1      1  <NA> 10      0     mae   standard 0.1744917
#> 2 2         1      1  <NA> 10      0    rmse   standard 0.2441705
#> 3 4         2      1  <NA> 10      0     mae   standard 0.1528578
#> 4 4         2      1  <NA> 10      0    rmse   standard 0.2000098

# 2. Tune with fixed NA positions
na_positions <- list(
  matrix(c(1, 2, 3, 1, 1, 1), ncol = 2),
  matrix(c(2, 3, 4, 2, 2, 2), ncol = 2)
)

results_fixed <- tune_imp(
  obj,
  data.frame(k = 2),
  .f = "knn_imp",
  na_loc = na_positions
)
#> Tuning knn_imp
#> Step 1/2: Resolving NA locations
#> Running Mode: sequential...
#> Step 2/2: Tuning

# 3. Custom imputation function
custom_fill <- function(obj, val = 0) {
  obj[is.na(obj)] <- val
  obj
}
tune_imp(obj, data.frame(val = c(0, 1)), .f = custom_fill, num_na = 10)
#> Tuning custom function
#> Step 1/2: Resolving NA locations
#> Running Mode: sequential...
#> Step 2/2: Tuning
#> # slideimp table: 2 x 5
#>  val param_set rep_id        result error
#>    0         1      1 <df [10 x 2]>  <NA>
#>    1         2      1 <df [10 x 2]>  <NA>
if (FALSE) { # interactive() && requireNamespace("mirai", quietly = TRUE)
# 4. Parallel tuning (requires mirai package)
mirai::daemons(2)
parameters_custom <- data.frame(mean = c(0, 1), sd = c(1, 1))

# Define a simple custom function for illustration
custom_imp <- function(obj, mean, sd) {
  na_pos <- is.na(obj)
  obj[na_pos] <- stats::rnorm(sum(na_pos), mean = mean, sd = sd)
  obj
}

results_p <- tune_imp(
  obj, parameters_custom, .f = custom_imp, n_reps = 1, num_na = 10
)
mirai::daemons(0) # Close workers
}
```
