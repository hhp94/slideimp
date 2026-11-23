# Hyperparameter Tuning for Imputation Methods

Tune method-specific hyperparameters by repeatedly masking observed
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
  pin_blas = TRUE
)
```

## Arguments

- obj:

  A double matrix.

- parameters:

  A data frame or `NULL`. One row per parameter set, one column per
  argument of `.f` other than `obj`. List-columns hold arguments that
  are not scalars. Duplicate rows are removed. `NULL` is allowed only
  for a custom `.f`, which is then called with its defaults.

- .f:

  Character or function. One of `"knn_imp"`, `"pca_imp"`, or
  `"slide_imp"`, or a custom imputation function, as described in
  Details.

- na_loc:

  A matrix, numeric vector, list, or `NULL`. Positions to mask in place
  of sampled ones. A two-column matrix holds row and column indices, a
  numeric vector holds linear positions, and a list holds one of those
  per repetition.

- num_na:

  Integer or `NULL`. Total number of missing values per repetition. When
  given, `n_cols` is ignored and the values are spread as evenly as
  possible over `num_na %/% n_rows` columns.

- n_reps:

  Integer. Number of independent repetitions.

- n_cols:

  Integer or `NULL`. Number of columns that receive missing values.

- n_rows:

  Integer. Minimum number of missing values per selected column.

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

- .progress:

  Logical. If `TRUE`, show progress during tuning.

- cores:

  Integer. Number of cores for K-NN imputation, including
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  with `k`. For other methods, use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html).

- location:

  Numeric or `NULL`. The position of each column of `obj`, sorted in
  ascending order. `.f = "slide_imp"` needs it, and the other methods
  ignore it.

- pin_blas:

  Logical. If `TRUE`, pin BLAS threads to 1 on the `mirai` daemons in a
  parallel run. Defaults to `TRUE`.

## Value

A data frame of class `slideimp_tune` with these columns:

- the columns of `parameters`.

- `param_set`, an integer ID for each unique parameter combination.

- `rep_id`, an integer repetition index.

- `result`, a list-column of data frames, each with `truth` and
  `estimate` columns.

- `error`, a character column with the error message of a failed
  iteration, or `NA`.

## Details

A custom `.f` must accept `obj` as its first argument and return a
numeric matrix with the same dimensions as `obj` and no `Inf` values.
`tune_imp()` calls it once with the first parameter set before tuning
starts, and stops if that call fails. After that, an iteration that
fails or breaks this rule is recorded in the `error` column.

The columns of `parameters` for a custom function must not use the names
of the columns that `tune_imp()` and
[`compute_metrics()`](https://hhp94.github.io/slideimp/reference/compute_metrics.md)
add: `param_set`, `rep_id`, `result`, `error`, `n`, `n_miss`, `.metric`,
`.estimator`, and `.estimate`.

When `mirai` daemons are set, a custom function runs on the daemons. It
finds the packages attached on the daemons, but not the variables of the
session that defined it. Pass such values through `parameters`.

When `mirai` daemons are set, each daemon receives one copy of `obj`
before tuning starts and keeps it until tuning ends. Each iteration
masks cells in a second copy. Memory use therefore grows with the number
of daemons.

When `na_loc` is supplied, `num_na`, `n_reps`, `n_cols`, `n_rows`,
`rowmax`, `colmax`, `na_col_subset`, and `max_attempts` are ignored. The
positions are not checked against `rowmax` or `colmax`.

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

## Missing values and non-finite input

`obj` must be a double matrix. An integer or logical matrix is rejected,
not converted. Convert one with `storage.mode(obj) <- "double"`.

`NA` and `NaN` are both treated as missing. The package does not
distinguish between them.

`Inf` and `-Inf` are rejected on entry, and so is a column with no
observed values (every cell `NA` or `NaN`). Locate them with
`which(is.infinite(obj), arr.ind = TRUE)` and
`which(colSums(!is.na(obj)) == 0)`.

## Parallelization

- K-NN: use the `cores` argument. If `mirai` daemons are active, `cores`
  is set to `1` to avoid nested parallelism.

- PCA: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

- With `mirai` daemons, each daemon runs BLAS on one thread unless
  `pin_blas = FALSE`.

## PCA Performance tips

Three arguments set the speed: `solver`, `threshold`, and `scale`. Tune
these first, then the accuracy arguments `ncp` and `coeff.ridge`, on a
representative subset.

**Exact vs. LOBPCG.** `"lobpcg"` suits large, approximately low-rank
matrices with small `ncp`. `"exact"` suits small matrices, such as
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
windows, where it is faster. It is also more reliable: it always
converges, to eigenvectors exact to machine precision, while LOBPCG
stops at the tolerance set in
[`lobpcg_control()`](https://hhp94.github.io/slideimp/reference/lobpcg_control.md)
and can take many iterations, or reach `maxiter`, when the eigenvalues
around the `ncp`-th one are close together. The warm start makes each
later LOBPCG solve cheap:
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
starts LOBPCG from the previous eigenblock and search direction, so once
the imputed values stabilize, each solve takes a few iterations. The
gain grows with the number of EM iterations, whatever the rank.
`solver = "auto"`, the default, chooses between the two and is a safe
start.

**Threshold.** The default `1e-6` is conservative. `1e-5` is often
faster and gives very similar values.

**Scale.** For columns on a common scale, such as DNA methylation beta
values in `[0, 1]`, `scale = FALSE` can be faster and more accurate.

**BLAS.** On Windows, the BLAS that ships with R can be slow, and
[OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows) can
replace it.

See [Speeding up PCA
imputation](https://hhp94.github.io/slideimp/articles/speeding-up-pca-imputation.html)
for the full workflow.

## See also

[`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
to draw the `na_loc` locations in advance.

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
#>   k param_set rep_id error  n n_miss .metric .estimator .estimate
#> 1 2         1      1  <NA> 10      0     mae   standard 0.1735950
#> 2 2         1      1  <NA> 10      0    rmse   standard 0.2306389
#> 3 4         2      1  <NA> 10      0     mae   standard 0.2173526
#> 4 4         2      1  <NA> 10      0    rmse   standard 0.2634978

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
#> Tuning a custom function
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
