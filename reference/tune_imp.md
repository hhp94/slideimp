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

  The imputation method to tune. Either a character string (`"knn_imp"`,
  `"pca_imp"`, or `"slide_imp"`) or a custom function.

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
  supplied, `n_cols` is computed automatically and the `NA` are
  distributed as evenly as possible using `n_rows` as the base (`num_na`
  must be `>= n_rows`). If omitted but `n_cols` is supplied, exactly
  `n_cols * n_rows` `NA` are injected. If `num_na`, `n_cols`, and
  `na_loc` are all `NULL`, defaults to ~5% of total cells (capped
  at 500) in `tune_imp()`.
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  has no default `num_na`. Ignored in `tune_imp()` when `na_loc` is
  supplied.

- n_reps:

  Integer. Number of repetitions for random NA injection (default `1`).

- n_cols:

  Integer. The number of columns to receive injected `NA` per
  repetition. Ignored when `num_na` is supplied (in which case `n_cols`
  is derived as `as.integer(num_na %/% n_rows)`). Must be provided if
  `num_na` is `NULL`. Ignored in `tune_imp()` when `na_loc` is supplied.

- n_rows:

  Integer. The target number of `NA` per column (default `2L`).

  - When `num_na` is supplied: used as the base size. Most columns
    receive exactly `n_rows` `NA`; `num_na %% n_rows` columns receive
    `n_rows + 1`. If there's only one column then that column receives
    all the remainder.

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
    the number derived from `num_na`).

  Ignored in `tune_imp()` when `na_loc` is supplied.

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

A `data.frame` of class
`c("slideimp_tune", "slideimp_tbl", "data.frame")` containing the
following columns:

- `...`: All columns originally provided in `parameters`.

- `param_set`: An integer ID representing the unique parameter
  combination.

- `rep_id`: An integer indicating the repetition index.

- `result`: A nested list-column where each element is a data.frame
  containing `truth` (original values) and `estimate` (imputed values).

- `error`: A character column containing the error message if the
  iteration failed, otherwise `NA`.

## Details

The function supports tuning for built-in imputation methods
(`"slide_imp"`, `"knn_imp"`, `"pca_imp"`) or custom functions provided
via `.f`.

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

Parallelization behavior depends on the imputation method:

- **KNN**: use the `cores` argument (if OpenMP is available). If
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  are also active, `cores` is automatically set to 1 to avoid nested
  parallelism.

- **PCA**: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

**Linux / OpenBLAS / MKL users:** If your machine uses a multi-threaded
BLAS (e.g., OpenBLAS or Intel MKL), set `pin_blas = TRUE` when tuning
PCA imputation in parallel. Without it, BLAS threads and `mirai` workers
compete for cores, which can cause slowdowns (CPU thrashing).

**macOS users:** OpenMP is typically unavailable on macOS unless
manually configured. `cores` will fall back to 1 automatically; use
[`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html) for
parallelization instead.

## Examples

``` r
# Setup example data
data(khanmiss1)
obj <- t(khanmiss1)[1:20, sample.int(nrow(khanmiss1), size = 200)]

# 1. Tune K-NN imputation with random NA injection
params_knn <- data.frame(k = c(5, 10))
results <- tune_imp(obj, params_knn, .f = "knn_imp", n_reps = 1, num_na = 20)
#> Tuning knn_imp
#> Step 1/2: Resolving NA locations
#> Running Mode: sequential...
#> Step 2/2: Tuning
compute_metrics(results)
#>    k param_set rep_id error  n n_miss .metric .estimator .estimate
#> 1  5         1      1  <NA> 20      0     mae   standard  463.8675
#> 2  5         1      1  <NA> 20      0    rmse   standard  634.3499
#> 3 10         2      1  <NA> 20      0     mae   standard  485.9833
#> 4 10         2      1  <NA> 20      0    rmse   standard  625.0286

# 2. Tune with fixed NA positions
na_positions <- list(
  matrix(c(1, 2, 3, 1, 1, 1), ncol = 2),
  matrix(c(2, 3, 4, 2, 2, 2), ncol = 2)
)
results_fixed <- tune_imp(
  obj,
  data.frame(k = 10),
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

# 4. Parallel tuning (requires mirai package)
mirai::daemons(2)
parameters_custom <- data.frame(mean = c(0, 1), sd = c(1, 1))

# Define a simple custom function for illustration
custom_imp <- function(obj, mean, sd) {
  na_pos <- is.na(obj)
  obj[na_pos] <- stats::rnorm(sum(na_pos), mean = mean, sd = sd)
  obj
}

results_p <- tune_imp(obj, parameters_custom, .f = custom_imp, n_reps = 1, num_na = 10)
#> Tuning custom function
#> Step 1/2: Resolving NA locations
#> Running Mode: parallel...
#> Step 2/2: Tuning
#> Tip: set `pin_blas = TRUE` may improve parallel performance.
mirai::daemons(0) # Close workers
```
