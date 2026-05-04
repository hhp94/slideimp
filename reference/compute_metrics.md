# Compute Prediction Accuracy Metrics

Compute prediction accuracy metrics for results from
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

  A `slideimp_tune` data frame from
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).
  Must contain a `result` list-column whose elements are data frames
  with `truth` and `estimate` columns.

- metrics:

  Character vector of metric names to compute. Defaults to
  `c("mae", "rmse")`. Available metrics are `"mae"`, `"rmse"`, `"mape"`,
  `"bias"`, `"rsq"`, and `"rsq_trad"`.

## Value

A data frame containing the original parameter columns along with
unnested metric columns: `.metric`, `.estimator`, and `.estimate`.

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
#>   k .progress param_set rep_id error  n n_miss .metric .estimator .estimate
#> 1 5     FALSE         1      1  <NA> 10      0     mae   standard 0.1215382
#> 2 5     FALSE         1      1  <NA> 10      0    rmse   standard 0.1676796
```
