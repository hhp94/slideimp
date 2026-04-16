# Compute Prediction Accuracy Metrics

Computes prediction accuracy metrics for results from
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

  A `slideimp_tune` data.frame from
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).
  Must contain a `result` list-column with data.frames that have `truth`
  and `estimate` columns.

- metrics:

  A character vector of metric names to compute. Defaults to
  `c("mae", "rmse")`. Also available: `"mape"`, `"bias"`, `"rsq"`, and
  `"rsq_trad"`.

## Value

A data.frame with the original parameters along with unnested metrics:
`.metric`, `.estimator`, and `.estimate`.

## Details

For alternative or faster metrics, see the `{yardstick}` package.

## Examples

``` r
obj <- sim_mat(100, 100)$input

set.seed(1234)
results <- tune_imp(
  obj = obj,
  parameters = data.frame(k = 10),
  .f = "knn_imp",
  n_reps = 1,
  num_na = 20
)
#> Tuning knn_imp
#> Step 1/2: Resolving NA locations
#> Running Mode: sequential...
#> Step 2/2: Tuning

compute_metrics(results)
#>    k param_set rep_id error  n n_miss .metric .estimator .estimate
#> 1 10         1      1  <NA> 20      0     mae   standard 0.1046064
#> 2 10         1      1  <NA> 20      0    rmse   standard 0.1393337
```
