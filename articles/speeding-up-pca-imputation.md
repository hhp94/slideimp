# Speeding up PCA imputation

``` r

library(slideimp)
```

[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
does the PCA imputation in
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md),
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md),
and
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
when they impute with PCA. This article shows how to make PCA imputation
faster with the `threshold`, `scale`, and `solver` arguments, and with a
parallel backend from [mirai](https://mirai.r-lib.org).

Before you tune anything, run
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
with `.progress = TRUE` on a meaningful subset of the variables, for
example chromosome 22, to see how fast it runs. Ctrl+C stops a run
before it ends.

## Arguments that set the speed

For large numeric matrices, three arguments can reduce the runtime by a
large amount: `threshold`, `scale`, and `solver`.

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

### `threshold`

`threshold` is the convergence threshold of the PCA imputation
algorithm. The default, `threshold = 1e-6`, matches the conservative
behavior of [missMDA](http://factominer.free.fr/missMDA/index.md). For
many genomic matrices, `threshold = 1e-5` is much faster and gives very
similar imputed values.

### `scale`

The columns are always centered before PCA, and `scale` sets whether
they are also scaled to unit variance. The default `scale = TRUE` adds
work to every iteration. For data with columns on the same scale, such
as DNAm beta values in `[0, 1]`, `scale = TRUE` can amplify columns with
low variance and columns that are mostly noise. For such data,
`scale = FALSE` can be both faster and more accurate.

### `solver`

`solver` takes three values: two PCA eigensolvers, and `"auto"`, which
chooses between them:

- `solver = "exact"` uses an exact solver. It always converges, and it
  is typically faster for small to moderate matrices.
- `solver = "lobpcg"` uses an iterative eigensolver with warm starts. It
  can be faster for large, approximately low-rank matrices with a small
  `ncp`, especially when the imputation needs many EM iterations.
- `solver = "auto"`, the default, uses `"exact"` for small problems. For
  other problems, it runs a short timed probe of both solvers and picks
  `"lobpcg"` only when it is clearly faster.

`solver = "auto"` is a good default for most users. Set `solver` to
`"exact"` or `"lobpcg"` when:

- You run many independent
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  calls, for example in
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  or
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).
  Each call can run its own probe, and the probes add time.
- You know which solver is faster on your data, and `"auto"` picks the
  other one.

If neither solver is clearly faster on a representative subset, pick
`solver = "exact"`, because it always converges.

## Tuning workflow

The best settings depend on the hardware, the BLAS library, the matrix
size, and the missing-data pattern. Before an expensive analysis, tune
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md) in
these steps:

1.  Choose a representative subset of the data, for example chromosome
    22.
2.  Time the `exact` and `lobpcg` solvers on the subset at a fixed
    `threshold`.
3.  Use
    [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
    to check that a relaxed `threshold` and `scale = FALSE` do not
    materially change the cross-validation error.
4.  Fix `solver`, `threshold`, and `scale`, then tune the other
    parameters (`ncp` and `coeff.ridge`, and `window_size` and
    `overlap_size` for
    [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)).

### Create a representative test matrix

In a real DNA methylation analysis, use a representative chromosome or
subset. Here, `group1` of the simulated data stands in for that
chromosome.

``` r

set.seed(12345)

obj_probe <- obj[, subset(sim_df, group == "group1")$feature]

dim(obj_probe)
#> [1]  80 731
mean(is.na(obj_probe))
#> [1] 0.09849521
```

### Time the solvers

Use a small `ncp` and a relaxed `threshold = 1e-5` for the timing. This
`ncp` does not have to be the final `ncp`. Both solvers use the same
convergence criterion, so their relative timing is representative.

Run this timing on your own representative subset. On the simulated data
of this article, the two solvers take a similar time, so the chunk below
is not run here. A positive difference means that `"lobpcg"` was faster.

``` r

ncp_probe <- 5
```

``` r

exact <- system.time(invisible(
  pca_imp(
    obj = obj_probe,
    ncp = ncp_probe,
    solver = "exact",
    threshold = 1e-5
  )
))

iterative <- system.time(invisible(
  pca_imp(
    obj = obj_probe,
    ncp = ncp_probe,
    solver = "lobpcg",
    threshold = 1e-5
  )
))

# diff
exact[["elapsed"]] - iterative[["elapsed"]]
```

### Check `threshold`

Use
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
to compare the cross-validation error across `threshold` values. In this
example, the error at `1e-5` is very close to the error at `1e-6`.

``` r

thresh_df <- data.frame(threshold = c(1e-5, 1e-6), ncp = ncp_probe)
res <- tune_imp(obj = obj_probe, parameters = thresh_df, .f = "pca_imp")
#> Tuning `pca_imp()`
#> Step 1/2: Resolving NA locations
#> Using default `num_na` = 500 (~0.85% of cells).
#> ℹ Increase `num_na` for more stable metrics, or decrease it when `obj` already
#>   has many missing values.
#> Running mode: sequential
#> Step 2/2: Tuning

compute_metrics(res)
#>   threshold ncp param_set rep_id error   n n_miss .metric .estimator  .estimate
#> 1     1e-05   5         1      1  <NA> 500      0     mae   standard 0.09335752
#> 2     1e-05   5         1      1  <NA> 500      0    rmse   standard 0.11652556
#> 3     1e-06   5         2      1  <NA> 500      0     mae   standard 0.09320836
#> 4     1e-06   5         2      1  <NA> 500      0    rmse   standard 0.11641720
```

### Check `scale`

`scale = FALSE` can be faster and more accurate for some DNAm data, but
it is not always more accurate. Compare its error against the default
with
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
does not report the time, so time the two settings as in the solver
timing above.

``` r

scale_df <- data.frame(
  scale = c(TRUE, FALSE),
  ncp = ncp_probe,
  threshold = 1e-5
)
res_scale <- tune_imp(obj = obj_probe, parameters = scale_df, .f = "pca_imp")
#> Tuning `pca_imp()`
#> Step 1/2: Resolving NA locations
#> Using default `num_na` = 500 (~0.85% of cells).
#> ℹ Increase `num_na` for more stable metrics, or decrease it when `obj` already
#>   has many missing values.
#> Running mode: sequential
#> Step 2/2: Tuning

compute_metrics(res_scale)
#>   scale ncp threshold param_set rep_id error   n n_miss .metric .estimator
#> 1  TRUE   5     1e-05         1      1  <NA> 500      0     mae   standard
#> 2  TRUE   5     1e-05         1      1  <NA> 500      0    rmse   standard
#> 3 FALSE   5     1e-05         2      1  <NA> 500      0     mae   standard
#> 4 FALSE   5     1e-05         2      1  <NA> 500      0    rmse   standard
#>    .estimate
#> 1 0.08928488
#> 2 0.11088626
#> 3 0.08957511
#> 4 0.11124100
```

In this example, the error with `scale = FALSE` is very close to the
error with the default, and slightly higher.

### Tune accuracy parameters with the chosen settings

After you choose `solver`, `threshold`, and `scale`, keep them fixed
while you tune the accuracy arguments `ncp` and `coeff.ridge`. On the
simulated data here, the two solvers take a similar time, so this
example keeps `"exact"`. On a large, approximately low-rank matrix that
needs many EM iterations, `"lobpcg"` with warm starts is typically
faster. The example also uses `scale = FALSE`, because its error is very
close to the default and `scale = FALSE` can be faster.

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
#> Using default `num_na` = 500 (~0.85% of cells).
#> ℹ Increase `num_na` for more stable metrics, or decrease it when `obj` already
#>   has many missing values.
#> Running mode: sequential
#> Step 2/2: Tuning

metrics <- compute_metrics(tune_pca, metrics = "rmse")

aggregate(
  .estimate ~ ncp + coeff.ridge,
  data = metrics,
  FUN = mean
)
#>   ncp coeff.ridge .estimate
#> 1   2         0.8 0.1124445
#> 2   4         0.8 0.1128692
#> 3   6         0.8 0.1133943
#> 4   2         1.0 0.1122923
#> 5   4         1.0 0.1125142
#> 6   6         1.0 0.1127536
#> 7   2         1.2 0.1121828
#> 8   4         1.2 0.1122631
#> 9   6         1.2 0.1123292
```

### Apply to `slide_imp()`

The same steps apply to sliding-window PCA imputation.
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
typically works on small windows, each with few columns, so
`solver = "exact"` is usually a good default. With `solver = "auto"`,
the first window that completes the solver probe picks the solver, and
every later window uses it.

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
#> Using default `num_na` = 500 (~0.85% of cells).
#> ℹ Increase `num_na` for more stable metrics, or decrease it when `obj` already
#>   has many missing values.
#> Running mode: sequential
#> Step 2/2: Tuning

slide_metrics <- compute_metrics(tune_slide_pca, metrics = "rmse")

aggregate(
  .estimate ~ ncp + coeff.ridge + window_size,
  data = slide_metrics,
  FUN = mean
)
#>    ncp coeff.ridge window_size .estimate
#> 1    2         0.8          50 0.1111740
#> 2    4         0.8          50 0.1162099
#> 3    2         1.0          50 0.1106076
#> 4    4         1.0          50 0.1147531
#> 5    2         1.2          50 0.1104182
#> 6    4         1.2          50 0.1134466
#> 7    2         0.8         500 0.1107731
#> 8    4         0.8         500 0.1134791
#> 9    2         1.0         500 0.1104863
#> 10   4         1.0         500 0.1126640
#> 11   2         1.2         500 0.1102355
#> 12   4         1.2         500 0.1119397
```

Finally, impute with the selected parameters. Here the call runs on
`beta_matrix`, the `group1` subset. In a real analysis, pass the full
matrix:

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

## Parallel PCA imputation and BLAS threading

[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
and
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
can also make PCA imputation faster by running independent imputation
tasks in parallel on [mirai](https://mirai.r-lib.org) daemons. Their
`cores` argument applies to K-NN imputation only.

When PCA imputation runs in parallel, each daemon can call multithreaded
linear algebra routines through BLAS. If each daemon uses several BLAS
threads, the total number of active CPU threads can far exceed the
number of physical cores. For example, 8
[mirai](https://mirai.r-lib.org) daemons that each use 8 BLAS threads
request 64 threads. This **thread oversubscription** can make code
slower rather than faster.

[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
and
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
prevent this by default: `pin_blas = TRUE` uses
[RhpcBLASctl](https://prs.ism.ac.jp/~nakama/Rhpc/) to limit BLAS to one
thread on each daemon. On Linux, this works with common BLAS libraries
such as OpenBLAS, FlexiBLAS, and MKL. Set `pin_blas = FALSE` to let each
daemon use a threaded BLAS, for example when a few daemons run on many
cores.

On Windows, with OpenBLAS in place of R’s own BLAS, `pin_blas = TRUE`
does not limit every BLAS and LAPACK call that PCA imputation makes. To
keep each daemon on one thread, set `OPENBLAS_NUM_THREADS=1` before R
starts, for example in `.Renviron`. Or do not start
[mirai](https://mirai.r-lib.org) daemons, so the tasks run in sequence,
and let the multithreaded OpenBLAS use the cores.

``` r

library(mirai)

mirai::daemons(4)

pca_group_results <- group_imp(
  obj = obj,
  group = sim_df,
  ncp = 2,
  solver = chosen_solver,
  threshold = chosen_threshold,
  scale = chosen_scale
)
#> Imputing 2 groups using PCA.
#> Running mode: mirai

mirai::daemons(0)
```
