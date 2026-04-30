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
PCA imputation faster by choosing the convergence threshold, the PCA
eigensolver, and a parallel backend with
[mirai](https://mirai.r-lib.org).

### Overview

Most users can keep the defaults:

``` r
set.seed(1234)
sim <- sim_mat(
  n = 80,
  p = 1500,
  rho = 0.7,
  perc_total_na = 0.10,
  perc_col_na = 0.8
)
sim_df <- sim$col_group
obj <- sim$input
pca_imp(obj, ncp = 2)
#> Method: PCA imputation
#> Dimensions: 80 x 1500
#> 
#>           feature1  feature2   feature3  feature4  feature5  feature6
#> sample1 0.25105224 0.2690159 0.22175022 0.3808861 0.1697869 0.3006772
#> sample2 0.46488252 0.5570005 0.56877205 0.6171557 0.5943715 0.6692621
#> sample3 0.46781185 0.6658607 0.53191864 0.6945785 0.6875948 0.5242096
#> sample4 0.07870941 0.1605554 0.07784548 0.2108805 0.0000000 0.2506738
#> sample5 0.52579469 0.6477465 0.63987522 0.6117717 0.6680913 0.6045150
#> sample6 0.58204683 0.6815553 0.59358698 0.5741023 0.6504755 0.8734424
#> # Showing 6 of 80 rows and 6 of 1500 columns
```

However, for large genomic matrices, runtime can depend strongly on two
arguments:

- `threshold`: the convergence threshold for the PCA imputation
  algorithm.
  - The default is `threshold = 1e-6`, matching the conservative
    behavior of [missMDA](http://factominer.free.fr/missMDA/index.md).
  - For many genomic matrices, `threshold = 1e-5` can be much faster
    while producing very similar imputed values.
- `solver`: the eigensolver used inside
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).
  - `solver = "auto"` runs a short timed probe of both solvers and picks
    the faster one.
  - `solver = "exact"` forces the exact solver.
  - `solver = "lobpcg"` forces the iterative LOBPCG solver.

### Exact versus iterative PCA solvers

[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
supports two PCA eigensolvers:

- The exact solver, selected by `solver = "exact"`.
- The iterative LOBPCG solver, selected by `solver = "lobpcg"`.

The exact solver is robust and often fast for small to moderate
matrices. The iterative LOBPCG solver can be faster for large or
approximately low-rank genomic matrices, especially when only a small
number of principal components is needed.

The default `solver = "auto"` runs a small number of probe iterations of
each solver at the start of the run, measures wall-clock time, and picks
LOBPCG only if it is clearly faster than the exact solver. This decision
is data-dependent and is made per call.

For most users, `solver = "auto"` is a good default. Manual benchmarking
is still useful in two situations:

- Before a long
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  run, where amortizing the probe overhead across many calls matters
  less than knowing the right solver up front.
- When you suspect the auto choice is wrong (for example, very large or
  very low-rank matrices where LOBPCG should clearly win, or very small
  windows where the exact solver should clearly win).

### Recommended workflow

A practical workflow is:

1.  Choose a representative subset of the data, for example chromosome
    22.
2.  Time the `exact` and `lobpcg` solvers on the subset at a fixed
    `threshold`.
3.  Use
    [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
    to check that a relaxed `threshold` does not materially change
    cross-validation error.
4.  Use the chosen `solver` and `threshold` when tuning accuracy-related
    parameters such as `ncp`, `coeff.ridge`, `window_size`, or
    `overlap_size`.

### Create a representative test matrix

In a real methylation analysis, use a representative chromosome or
subset. Here, we designate `group1` of the simulated data as the
representative chromosome.

``` r
set.seed(1234)

obj_probe <- obj[, subset(sim_df, group == "group1")$feature]

dim(obj_probe)
#> [1]  80 731
mean(is.na(obj_probe))
#> [1] 0.09849521
```

### Time solver and threshold settings

First, choose a small `ncp` for the speed probe. This does not need to
be the final `ncp`. The goal is to identify whether the exact or
iterative solver is faster on this kind of data. We use a relaxed
`threshold = 1e-5` to keep the probe quick. Both solvers stop at the
same convergence criterion, so their relative timing is representative.

``` r
ncp_probe <- 5
```

``` r
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

exact - iterative
#>          user        system       elapsed 
#> -2.000000e-03 -6.000000e-03 -2.220446e-16
```

Second, use
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
to compare cross-validation error across `threshold` values. In this
example, relaxing the threshold to `1e-5` did not meaningfully change
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
#>   threshold ncp param_set rep_id error   n n_miss .metric .estimator  .estimate
#> 1     1e-05   5         1      1  <NA> 500      0     mae   standard 0.09174255
#> 2     1e-05   5         1      1  <NA> 500      0    rmse   standard 0.11499555
#> 3     1e-06   5         2      1  <NA> 500      0     mae   standard 0.09174508
#> 4     1e-06   5         2      1  <NA> 500      0    rmse   standard 0.11500090
```

If `threshold = 1e-5` changes cross-validation error or imputed values
more than expected, use the more conservative `threshold = 1e-6`.

If neither solver is clearly faster on the representative subset,
`solver = "auto"` is fine. If one solver is clearly faster, force it
explicitly to skip the auto probe overhead.

### Use the chosen settings with `tune_imp()`

After choosing `solver` and `threshold`, keep them fixed while tuning
accuracy-related PCA parameters (i.e., `ncp` and `coeff.ridge`):

``` r
chosen_solver <- "exact"
chosen_threshold <- 1e-5

params <- expand.grid(
  ncp = c(2, 4, 6),
  coeff.ridge = c(0.8, 1, 1.2)
)

params$solver <- chosen_solver
params$threshold <- chosen_threshold

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
#> 1   2         0.8 0.1125669
#> 2   4         0.8 0.1132658
#> 3   6         0.8 0.1148464
#> 4   2         1.0 0.1123355
#> 5   4         1.0 0.1128064
#> 6   6         1.0 0.1138270
#> 7   2         1.2 0.1121501
#> 8   4         1.2 0.1124116
#> 9   6         1.2 0.1130121
```

### Use the chosen settings with `slide_imp()`

The same idea applies when tuning sliding-window PCA imputation. First
choose the speed settings on a representative subset, then include those
fixed settings in the `parameters` data frame while tuning window and
PCA parameters.

[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
typically operates on small windows with few samples, so the `"exact"`
solver is usually the better default. Forcing `solver = "exact"` also
avoids paying the auto probe cost on every window.

``` r
beta_matrix <- obj_probe # We treat `obj_probe` as a slide_imp() input instead
locations <- seq_len(ncol(beta_matrix))
slide_params <- expand.grid(
  ncp = c(2, 4),
  coeff.ridge = c(0.8, 1, 1.2),
  window_size = c(50, 500) # Change based on units of `locations`
)

slide_params$overlap_size <- 5
slide_params$min_window_n <- 20
slide_params$solver <- "exact" # Usually better for `slide_imp()`
slide_params$threshold <- chosen_threshold

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
#> 1    2         0.8          50 0.1218555
#> 2    4         0.8          50 0.1257166
#> 3    2         1.0          50 0.1213016
#> 4    4         1.0          50 0.1242115
#> 5    2         1.2          50 0.1208738
#> 6    4         1.2          50 0.1228126
#> 7    2         0.8         500 0.1170125
#> 8    4         0.8         500 0.1185850
#> 9    2         1.0         500 0.1168271
#> 10   4         1.0         500 0.1179256
#> 11   2         1.2         500 0.1166733
#> 12   4         1.2         500 0.1173506
```

Finally, use the selected speed and accuracy parameters when imputing
the full dataset:

``` r
slide_results <- slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 500,
  overlap_size = 5,
  min_window_n = 20,
  ncp = 2,
  coeff.ridge = 1.2,
  solver = "exact", # Usually better for `slide_imp()`
  threshold = chosen_threshold,
  .progress = FALSE
)
```

### Parallel PCA imputation and BLAS threading

PCA imputation can also be sped up by running independent imputation
tasks in parallel within
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).
For PCA imputation, use [mirai](https://mirai.r-lib.org) parallelization
rather than the `cores` argument.

The `cores` argument is for K-NN imputation. For PCA imputation, start
[mirai](https://mirai.r-lib.org) daemons before calling functions such
as
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
or
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).

``` r
library(mirai)

# Start 4 mirai workers.
mirai::daemons(4)

# Run PCA tuning in parallel.
tune_pca <- tune_imp(
  obj = obj_probe,
  .f = "pca_imp",
  parameters = params,
  n_reps = 1,
  .progress = FALSE # <- `TRUE` to monitor longer running jobs
)
#> Tuning `pca_imp()`
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 500 (~0.9% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running mode: mirai
#> Step 2/2: Tuning
#> Tip: set `pin_blas = TRUE` may improve parallel performance.

# Stop mirai workers when finished.
mirai::daemons(0)
```

When PCA imputation is run in parallel, each worker may call
multi-threaded linear algebra routines through BLAS. If those workers
each use multiple BLAS threads, the total number of active CPU threads
can become much larger than the number of physical cores. This is called
**thread oversubscription**, and it can make code slower rather than
faster.

For example, 8 [mirai](https://mirai.r-lib.org) workers each using 8
BLAS threads can try to use 64 CPU threads.

To avoid this, set `pin_blas = TRUE` in
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
or
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
when running PCA imputation with [mirai](https://mirai.r-lib.org). This
requires [RhpcBLASctl](https://prs.ism.ac.jp/~nakama/Rhpc/):

``` r
mirai::daemons(4)

tune_slide_pca <- tune_imp(
  obj = beta_matrix,
  .f = "slide_imp",
  parameters = slide_params,
  n_reps = 1,
  location = locations,
  pin_blas = TRUE # <- Set to `TRUE` if R is linked to multi-threaded BLAS/LAPACK
)
#> Tuning `slide_imp()`
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 500 (~0.9% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running mode: mirai
#> Step 2/2: Tuning

mirai::daemons(0)
```

For group-wise PCA imputation:

``` r
mirai::daemons(4)

pca_group_results <- group_imp(
  obj = obj,
  group = sim_df,
  ncp = 2,
  solver = chosen_solver,
  threshold = chosen_threshold,
  pin_blas = TRUE # <- Set to `TRUE` if R is linked to multi-threaded BLAS/LAPACK
)
#> Imputing 2 groups using PCA.
#> Running mode: mirai

mirai::daemons(0)
```

`pin_blas = TRUE` uses
[RhpcBLASctl](https://prs.ism.ac.jp/~nakama/Rhpc/) to limit BLAS threads
inside parallel workers. This works with common BLAS libraries such as
OpenBLAS and MKL, which are typical on Linux and also available on
Windows after replacing R’s default BLAS. It may have no effect on
macOS, depending on the BLAS implementation.

### Practical recommendations

- Start with `solver = "auto"` unless PCA imputation is slow or you are
  about to run many calls (e.g., with
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md))
  where the per-call probe overhead adds up.
- If PCA imputation is slow, benchmark `solver = "exact"` against
  `solver = "lobpcg"` on a representative subset and force the winner.
- Try `threshold = 1e-5` or higher, depending on your error tolerance.
- Keep `threshold = 1e-5` only if it produces similar imputed values or
  similar cross-validation error.
- LOBPCG is most likely to help when the matrix is large, approximately
  low-rank, and `ncp` is small.
- For
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md),
  force `solver = "exact"`. Windows are usually small enough that exact
  wins, and forcing it avoids the per-window probe.
- Choose speed settings first, then tune accuracy-related parameters
  such as `ncp`, `coeff.ridge`, and `window_size`.
- If running PCA imputation sequentially, BLAS multi-threading may help,
  and `pin_blas` is unnecessary.
- If running many PCA imputations in parallel with
  [mirai](https://mirai.r-lib.org), use `pin_blas = TRUE`.
- Avoid using many [mirai](https://mirai.r-lib.org) workers and many
  BLAS threads at the same time.
- Benchmark on a representative subset, because the best setup depends
  on hardware, BLAS library, matrix size, and missing-data pattern.
