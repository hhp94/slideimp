# Compute Prediction Accuracy Metrics

Computes accuracy metrics for each parameter set and repetition in the
output of
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).

## Usage

``` r
compute_metrics(results, metrics = c("mae", "rmse"))

# S3 method for class 'data.frame'
compute_metrics(results, metrics = c("mae", "rmse"))

# S3 method for class 'slideimp_tune'
compute_metrics(results, metrics = c("mae", "rmse"))
```

## Arguments

- results:

  A `slideimp_tune` object or a data frame. The output of
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md),
  or a data frame whose `result` list-column holds data frames with
  `truth` and `estimate` columns.

- metrics:

  Character. The metrics to compute, any of `"mae"`, `"rmse"`, `"mape"`,
  `"bias"`, `"rsq"`, and `"rsq_trad"`. Defaults to `c("mae", "rmse")`.

## Value

A data frame with one row per row of `results` and metric. It keeps
every column of `results` except `result`, and adds `n`, the number of
masked values, `n_miss`, the number of those with no estimate,
`.metric`, `.estimator`, and `.estimate`.

## Details

Each metric is taken over the masked values of one row of `results`. A
masked value with no estimate is left out, and `n_miss` counts it.

- `"mae"`: mean absolute error.

- `"rmse"`: root mean squared error.

- `"mape"`: mean absolute percentage error, in percent. It divides by
  the true value, so it is not finite (`Inf` or `NaN`) for a repetition
  in which a masked value is exactly `0`. Use another metric for data
  that contain zeros.

- `"bias"`: mean of the estimate minus the true value.

- `"rsq"`: squared Pearson correlation of the estimates and the true
  values. It is `NA` when either is constant or fewer than two values
  have an estimate.

- `"rsq_trad"`: one minus the residual sum of squares over the total sum
  of squares. It is `NA` when the true values are constant or fewer than
  two values have an estimate.

## Examples

``` r
set.seed(1234)
obj <- sim_mat(20, 30)$input

results <- tune_imp(
  obj = obj,
  parameters = data.frame(k = 5),
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
#> 1 5         1      1  <NA> 10      0     mae   standard 0.1074189
#> 2 5         1      1  <NA> 10      0    rmse   standard 0.1302340
```
