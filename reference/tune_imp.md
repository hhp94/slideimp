# Tune Imputation Method Parameters

Tune method specific hyperparameters by repeatedly masking observed
values, imputing them, and comparing the imputed values with the
original values.

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

  A numeric matrix.

- parameters:

  A `data.frame` specifying parameter combinations to tune. Each column
  should be a parameter accepted by `.f`, excluding `obj`. List-columns
  are supported for complex parameters. Duplicate rows are removed.
  `NULL` is treated as a single parameter set with no additional
  arguments, which is useful for functions whose required arguments all
  have defaults.

- .f:

  One of `"knn_imp"`, `"pca_imp"`, or `"slide_imp"`, or a custom
  imputation function.

- na_loc:

  Optional predefined missing-value locations. Accepted formats are a
  two-column integer matrix of row and column indices, a numeric vector
  of linear positions, or a list whose elements are either of those
  formats.

- num_na:

  Integer or `NULL`. Total number of missing values to inject per
  repetition. If supplied, `n_cols` is derived from `num_na` and
  `n_rows`, and missing values are distributed as evenly as possible
  across columns.

- n_reps:

  Integer. Number of independent repetitions.

- n_cols:

  Integer or `NULL`. Must be supplied when both `num_na` and `na_loc`
  are `NULL`, unless the automatic default applies.

- n_rows:

  Integer. Target number of missing values to inject per selected
  column.

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

- .progress:

  Logical. If `TRUE`, show progress during tuning.

- cores:

  Integer. Number of cores to use for K-NN and sliding-window K-NN
  imputation. For other methods, use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html).

- location:

  Numeric vector of column locations. Required when `.f = "slide_imp"`.

- pin_blas:

  Logical. If `TRUE`, pin BLAS threads to 1 during parallel tuning to
  reduce thread contention.

## Value

A data frame of class `slideimp_tune` containing:

- columns originally provided in `parameters`;

- `param_set`, an integer ID for each unique parameter combination;

- `rep_id`, an integer repetition index;

- `result`, a list-column where each element is a data frame containing
  `truth` and `estimate` columns;

- `error`, a character column containing the error message if the
  iteration failed, otherwise `NA`.

## Details

Built-in methods can be selected by passing `.f = "knn_imp"`,
`.f = "pca_imp"`, or `.f = "slide_imp"`. A custom function can also be
supplied. Custom functions must accept `obj` as their first argument and
return a numeric matrix with the same dimensions as `obj`.

When `na_loc` is supplied, `num_na`, `n_cols`, `n_rows`, and
`na_col_subset` are ignored.

When `.f` is a character string, columns in `parameters` are validated
against the selected method:

- `"knn_imp"` requires `k`.

- `"pca_imp"` requires `ncp`.

- `"slide_imp"` requires `window_size`, `overlap_size`, and
  `min_window_n`, plus exactly one of `k` or `ncp`.

To tune parameters for grouped imputation, tune
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md) or
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md) on
representative groups, then pass the selected parameters to
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).

The top-level `rowmax` and `colmax` arguments control random
missing-value injection performed by
[`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md).
To tune or pass an imputation method's own `colmax` argument, include a
`colmax` column in `parameters`.

Tuning results can be summarized with
[`compute_metrics()`](https://hhp94.github.io/slideimp/reference/compute_metrics.md)
or evaluated with external packages such as `yardstick`.

## Parallelization

- K-NN: use the `cores` argument. If `mirai` daemons are active, `cores`
  is automatically set to `1` to avoid nested parallelism.

- PCA: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

When running PCA imputation in parallel with `mirai`, set
`pin_blas = TRUE` in `tune_imp()` or
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
to prevent BLAS threads from oversubscribing CPU cores. This relies on
`RhpcBLASctl` and works with OpenBLAS and MKL (typical on Linux, and on
Windows after an OpenBLAS swap). `pin_blas = TRUE` may have no effect on
macOS.

## Performance tips

[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
relies heavily on linear algebra. On Windows, the default BLAS shipped
with R may be slow for large matrices. Advanced users can replace it
with [OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).

PCA imputation speed depends on the eigensolver selected by `solver` and
the convergence threshold `threshold`. The exact solver is selected with
`solver = "exact"`. The iterative LOBPCG solver is selected with
`solver = "lobpcg"`. The default, `solver = "auto"`, performs a short
timed probe and chooses LOBPCG only when it is clearly faster.

For large or approximately low-rank genomic matrices, it can be useful
to benchmark `solver = "exact"` against `solver = "lobpcg"` on a
representative subset, such as chromosome 22, before tuning
accuracy-related parameters. For
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md),
this may include `window_size` and `overlap_size`.

The default `threshold = 1e-6` is conservative. In many genomic
datasets, `threshold = 1e-5` can be faster while giving very similar
imputed values. Check this on a representative subset before using the
relaxed threshold in a full analysis.

See the pkgdown article [Speeding up PCA
imputation](https://hhp94.github.io/slideimp/articles/speeding-up-pca-imputation.html)
for a full workflow.

## Examples

``` r
set.seed(123)

# Simulate some data
obj <- sim_mat(10, 50)$input

# Tune K-NN imputation with random missing-value injection.
# Use larger `num_na` and `n_reps` values for real analyses.
params_knn <- data.frame(k = c(2, 4))
results <- tune_imp(
  obj,
  params_knn,
  .f = "knn_imp",
  n_reps = 1,
  num_na = 10,
  .progress = FALSE
)
#> Tuning `knn_imp()`
#> Step 1/2: Resolving NA locations
#> Running mode: sequential
#> Step 2/2: Tuning
compute_metrics(results)
#>   k .progress param_set rep_id error  n n_miss .metric .estimator .estimate
#> 1 2     FALSE         1      1  <NA> 10      0     mae   standard 0.2075860
#> 2 2     FALSE         1      1  <NA> 10      0    rmse   standard 0.2841912
#> 3 4     FALSE         2      1  <NA> 10      0     mae   standard 0.1829599
#> 4 4     FALSE         2      1  <NA> 10      0    rmse   standard 0.2202581

# Tune with fixed missing-value positions
na_positions <- list(
  matrix(c(1, 2, 3, 1, 1, 1), ncol = 2),
  matrix(c(2, 3, 4, 2, 2, 2), ncol = 2)
)

results_fixed <- tune_imp(
  obj,
  data.frame(k = 2),
  .f = "knn_imp",
  na_loc = na_positions,
  .progress = FALSE
)
#> Tuning `knn_imp()`
#> Step 1/2: Resolving NA locations
#> Running mode: sequential
#> Step 2/2: Tuning

# Custom imputation function
custom_fill <- function(obj, val = 0) {
  obj[is.na(obj)] <- val
  obj
}

tune_imp(
  obj,
  data.frame(val = c(0, 1)),
  .f = custom_fill,
  num_na = 10,
  .progress = FALSE
)
#> Tuning custom function
#> Step 1/2: Resolving NA locations
#> Running mode: sequential
#> Step 2/2: Tuning
#> # slideimp table: 2 x 5
#>  val param_set rep_id        result error
#>    0         1      1 <df [10 x 2]>  <NA>
#>    1         2      1 <df [10 x 2]>  <NA>

if (FALSE) { # interactive() && requireNamespace("mirai", quietly = TRUE)
# Parallel tuning with mirai
mirai::daemons(2)

parameters_custom <- data.frame(mean = c(0, 1), sd = c(1, 1))

custom_imp <- function(obj, mean, sd) {
  na_pos <- is.na(obj)
  obj[na_pos] <- stats::rnorm(sum(na_pos), mean = mean, sd = sd)
  obj
}

results_p <- tune_imp(
  obj,
  parameters_custom,
  .f = custom_imp,
  n_reps = 1,
  num_na = 10,
  .progress = FALSE
)

mirai::daemons(0)
}
```
