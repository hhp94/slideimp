# Speeding up PCA imputation

``` r

library(slideimp)
```

## Tips to speed up PCA imputation

[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md) is
the PCA imputation workhorse used by
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md),
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md),
and
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
whenever PCA imputation is requested. This article describes how to make
PCA imputation faster by choosing `threshold`, `scale`, `solver`, and a
parallel backend with [mirai](https://mirai.r-lib.org).

Before any further tuning, try running
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
with `.progress = TRUE` on a meaningful subset of variables
(e.g. chromosome 22) to gauge how fast it runs (interruptible with
Ctrl+C mid-run).

### Overview

For large numeric matrices, runtime can be significantly reduced by
tuning three arguments: `threshold`, `scale`, and `solver`.

``` r

set.seed(1234) # simulate some data
sim <- sim_mat(
  n = 80,
  p = 1500,
  rho = 0.7,
  perc_total_na = 0.10,
  perc_col_na = 0.8
)
sim_df <- sim$col_group
obj <- sim$input
```

#### `threshold`

The convergence threshold for the PCA imputation algorithm. The default
`threshold = 1e-6` matches the conservative behavior of
[missMDA](http://factominer.free.fr/missMDA/index.md). For many genomic
matrices, `threshold = 1e-5` is much faster while producing very similar
imputed values.

#### `scale`

Whether columns are centered and scaled before PCA. The default
`scale = TRUE` adds per-iteration overhead. For data where columns are
on the same scale (e.g., DNAm beta values in `[0, 1]`), `scale = TRUE`
may amplify low-variance and noise-dominated columns. For such data,
`scale = FALSE` can be both faster and more accurate.

#### `solver`

[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
supports two PCA eigensolvers:

- `solver = "exact"` uses an exact solver. Robust and typically faster
  for small to moderate matrices.
- `solver = "lobpcg"` uses an iterative eigensolver with warm starts. It
  can be faster for large, approximately low-rank matrices with small
  ncp, especially when many EM iterations are needed. For small or
  moderate matrices, “exact” is often faster.
- `solver = "auto"` (default): runs a short timed probe of both solvers
  and picks LOBPCG only when it is clearly faster.

For most users, `solver = "auto"` is a good default. However, consider
setting the solver manually when:

- You are running many independent
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  calls (i.e., in
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  or
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)),
  where repeated probing adds overhead.
- You know one solver should dominate but `"auto"` is selecting the
  wrong one.

If neither solver is clearly faster on a representative subset, pick
`solver = "exact"` for robustness.

### Tuning workflow

Before an expensive analysis, a practical
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
tuning workflow is:

1.  Choose a representative subset of the data, for example chromosome
    22.
2.  Time the `exact` and `lobpcg` solvers on the subset at a fixed
    `threshold`.
3.  Use
    [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
    to check that a relaxed `threshold` and `scale = FALSE` do not
    materially change cross-validation error.
4.  Fix `solver`, `threshold`, and `scale`, then tune the other
    parameters (`ncp` and `coeff.ridge`, and `window_size` and
    `overlap_size` for
    [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)).

#### Create a representative test matrix

In a real DNA methylation analysis, use a representative chromosome or
subset. Here, we treat `group1` of the simulated data as the
representative chromosome.

``` r

set.seed(12345)

obj_probe <- obj[, subset(sim_df, group == "group1")$feature]

dim(obj_probe)
#> [1]  80 731
mean(is.na(obj_probe))
#> [1] 0.09849521
```

#### Time the solvers

Use a small `ncp` and a relaxed `threshold = 1e-5` for the probe. This
`ncp` is not necessarily the final `ncp`. Both solvers use the same
convergence criterion, so their relative timing is representative.

``` r

ncp_probe <- 5

exact <- system.time(invisible(
  pca_imp(
    obj = obj_probe,
    ncp = ncp_probe,
    solver = "exact",
    threshold = 1e-5,
    na_check = FALSE
  )
))

iterative <- system.time(invisible(
  pca_imp(
    obj = obj_probe,
    ncp = ncp_probe,
    solver = "lobpcg",
    threshold = 1e-5,
    na_check = FALSE
  )
))

# diff
exact[["elapsed"]] - iterative[["elapsed"]]
#> [1] 0.001
```

#### Check `threshold`

Use
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
to compare cross-validation error across `threshold` values. In this
example, relaxing the threshold to `1e-5` does not meaningfully change
the error.

``` r

thresh_df <- data.frame(threshold = c(1e-5, 1e-6), ncp = ncp_probe)
res <- tune_imp(obj = obj_probe, parameters = thresh_df, .f = "pca_imp")
#> Tuning `pca_imp()`
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 500 (~0.9% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running mode: sequential
#> Step 2/2: Tuning

compute_metrics(res)
#>   threshold ncp .progress param_set rep_id error   n n_miss .metric .estimator
#> 1     1e-05   5     FALSE         1      1  <NA> 500      0     mae   standard
#> 2     1e-05   5     FALSE         1      1  <NA> 500      0    rmse   standard
#> 3     1e-06   5     FALSE         2      1  <NA> 500      0     mae   standard
#> 4     1e-06   5     FALSE         2      1  <NA> 500      0    rmse   standard
#>    .estimate
#> 1 0.09327532
#> 2 0.11584434
#> 3 0.09326740
#> 4 0.11585427
```

#### Check `scale`

`scale = FALSE` can be both faster and more accurate for certain DNAm
data. We verify this with
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).

``` r

scale_df <- data.frame(
  scale = c(TRUE, FALSE),
  ncp = ncp_probe,
  threshold = 1e-5
)
res_scale <- tune_imp(obj = obj_probe, parameters = scale_df, .f = "pca_imp")
#> Tuning `pca_imp()`
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 500 (~0.9% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running mode: sequential
#> Step 2/2: Tuning

compute_metrics(res_scale)
#>   scale ncp threshold .progress param_set rep_id error   n n_miss .metric
#> 1  TRUE   5     1e-05     FALSE         1      1  <NA> 500      0     mae
#> 2  TRUE   5     1e-05     FALSE         1      1  <NA> 500      0    rmse
#> 3 FALSE   5     1e-05     FALSE         2      1  <NA> 500      0     mae
#> 4 FALSE   5     1e-05     FALSE         2      1  <NA> 500      0    rmse
#>   .estimator  .estimate
#> 1   standard 0.09281453
#> 2   standard 0.11612999
#> 3   standard 0.09290920
#> 4   standard 0.11589055
```

#### Tune accuracy parameters with the chosen settings

After choosing `solver`, `threshold`, and `scale`, keep them fixed while
tuning accuracy-related PCA parameters (`ncp` and `coeff.ridge`). Here,
`solver = "lobpcg"` is not clearly faster, and `scale = FALSE` has
little effect on RMSE. On this small, well-conditioned subset `exact`
wins. On a large, approximately low-rank matrix that needs many EM
iterations, `lobpcg` with warm start is typically faster.

``` r

chosen_solver <- "exact"
chosen_threshold <- 1e-5
chosen_scale <- FALSE

# tune ncp and coeff.ridge
params <- expand.grid(
  ncp = c(2, 4, 6),
  coeff.ridge = c(0.8, 1, 1.2)
)

params$solver <- chosen_solver
params$threshold <- chosen_threshold
params$scale <- chosen_scale

tune_pca <- tune_imp(
  obj = obj_probe,
  .f = "pca_imp",
  parameters = params,
  n_reps = 1 # <- increase to a higher number in real analyses
)
#> Tuning `pca_imp()`
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 500 (~0.9% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running mode: sequential
#> Step 2/2: Tuning

metrics <- compute_metrics(tune_pca, metrics = "rmse")

aggregate(
  .estimate ~ ncp + coeff.ridge,
  data = metrics,
  FUN = mean
)
#>   ncp coeff.ridge .estimate
#> 1   2         0.8 0.1158408
#> 2   4         0.8 0.1166093
#> 3   6         0.8 0.1177499
#> 4   2         1.0 0.1157635
#> 5   4         1.0 0.1162633
#> 6   6         1.0 0.1169795
#> 7   2         1.2 0.1157307
#> 8   4         1.2 0.1160604
#> 9   6         1.2 0.1164008
```

#### Apply to `slide_imp()`

The same idea applies when tuning sliding-window PCA imputation.
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
typically operates on small windows with few samples, so
`solver = "exact"` is usually a good default. With `solver = "auto"`,
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
probes only the first window and reuses the chosen solver for the rest.

``` r

beta_matrix <- obj_probe # Treat `obj_probe` as a `slide_imp()` input
locations <- seq_len(ncol(beta_matrix))
slide_params <- expand.grid(
  ncp = c(2, 4),
  coeff.ridge = c(0.8, 1, 1.2),
  window_size = c(50, 500) # Adjust based on units of `locations`
)

slide_params$overlap_size <- 5
slide_params$min_window_n <- 20
slide_params$solver <- "exact" # Usually better for `slide_imp()`
slide_params$threshold <- chosen_threshold
slide_params$scale <- chosen_scale

tune_slide_pca <- tune_imp(
  obj = beta_matrix,
  .f = "slide_imp",
  parameters = slide_params,
  n_reps = 1, # <- increase to a higher number in real analyses
  location = locations
)
#> Tuning `slide_imp()`
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 500 (~0.9% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running mode: sequential
#> Step 2/2: Tuning

slide_metrics <- compute_metrics(tune_slide_pca, metrics = "rmse")

aggregate(
  .estimate ~ ncp + coeff.ridge + window_size,
  data = slide_metrics,
  FUN = mean
)
#>    ncp coeff.ridge window_size .estimate
#> 1    2         0.8          50 0.1214082
#> 2    4         0.8          50 0.1258193
#> 3    2         1.0          50 0.1208723
#> 4    4         1.0          50 0.1242678
#> 5    2         1.2          50 0.1203357
#> 6    4         1.2          50 0.1229213
#> 7    2         0.8         500 0.1171656
#> 8    4         0.8         500 0.1180867
#> 9    2         1.0         500 0.1169246
#> 10   4         1.0         500 0.1175901
#> 11   2         1.2         500 0.1167663
#> 12   4         1.2         500 0.1172382
```

Finally, use the selected parameters when imputing the full dataset:

``` r

slide_results <- slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 500,
  overlap_size = 5,
  min_window_n = 20,
  ncp = 2,
  coeff.ridge = 1.2,
  solver = "exact",
  threshold = chosen_threshold,
  scale = chosen_scale,
  .progress = FALSE
)
```

### Parallel PCA imputation and BLAS threading

PCA imputation can also be sped up by running independent imputation
tasks in parallel with
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
or
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).
For PCA imputation, we use [mirai](https://mirai.r-lib.org)
parallelization rather than the `cores` argument.

When PCA imputation is run in parallel, each worker may call
multithreaded linear algebra routines through BLAS. If each worker uses
multiple BLAS threads, the total number of active CPU threads can far
exceed the number of physical cores, for example, 8
[mirai](https://mirai.r-lib.org) workers each using 8 BLAS threads will
request 64 threads. This **thread oversubscription** can make code
slower rather than faster.

Set `pin_blas = TRUE` to avoid this. It uses
[RhpcBLASctl](https://prs.ism.ac.jp/~nakama/Rhpc/) to limit BLAS threads
inside workers, and works with common BLAS libraries such as OpenBLAS
and MKL.

``` r

library(mirai)

mirai::daemons(4)

pca_group_results <- group_imp(
  obj = obj,
  group = sim_df,
  ncp = 2,
  solver = chosen_solver,
  threshold = chosen_threshold,
  scale = chosen_scale,
  pin_blas = TRUE # <- Set to `TRUE` if R is linked to multi-threaded BLAS/LAPACK
)
#> Imputing 2 groups using PCA.
#> Running mode: mirai

mirai::daemons(0)
```

The same `pin_blas = TRUE` argument is available in
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).

### Practical recommendations

- Start with `solver = "auto"` unless you are about to run many calls
  (e.g., with
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)),
  where per-call probe overhead adds up.
- If PCA imputation is slow, benchmark `"exact"` vs `"lobpcg"` on a
  representative subset and force the faster solver.
- Try `threshold = 1e-5` if it produces similar imputed values or
  cross-validation error.
- For data where columns already share a common scale (e.g., DNAm beta
  values), try `scale = FALSE` and verify with cross-validation.
- For
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md),
  `solver = "exact"` is usually the right choice.
- Choose speed settings first, then tune accuracy parameters such as
  `ncp`, `coeff.ridge`, and `window_size`.
- When running PCA imputation in parallel with
  [mirai](https://mirai.r-lib.org) and a multithreaded BLAS, use
  `pin_blas = TRUE` to avoid thread oversubscription.
- Benchmark on a representative subset. The best setup depends on the
  hardware, BLAS library, matrix size, and missing-data pattern.
