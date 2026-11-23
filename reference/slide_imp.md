# Sliding-Window K-NN or PCA Imputation

Impute missing values in a numeric matrix whose columns are ordered by
position, using K-NN or PCA imputation within sliding windows. For
Illumina DNA methylation microarrays, use
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
instead.

## Usage

``` r
slide_imp(
  obj,
  location,
  window_size,
  overlap_size = 0,
  flank = FALSE,
  min_window_n,
  subset = NULL,
  dry_run = FALSE,
  k = NULL,
  cores = 1,
  dist_pow = 0,
  ncp = NULL,
  scale = TRUE,
  coeff.ridge = 1,
  threshold = 1e-06,
  seed = NULL,
  row.w = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  solver = c("auto", "exact", "lobpcg"),
  lobpcg_control = NULL,
  clamp = NULL,
  method = NULL,
  .progress = TRUE,
  colmax = 0.9,
  post_imp = TRUE,
  on_infeasible = c("skip", "error", "mean")
)
```

## Arguments

- obj:

  A double matrix. Samples are in rows and features in columns.

- location:

  Numeric. The position of each column of `obj`, such as a genomic
  coordinate, sorted in ascending order.

- window_size:

  Numeric. The width of each window, in the units of `location`. With
  `flank = TRUE`, the distance each window reaches on each side of its
  feature.

- overlap_size:

  Numeric. The overlap between consecutive windows, in the units of
  `location`, less than `window_size`. Ignored when `flank = TRUE`.

- flank:

  Logical. If `FALSE`, slide windows across all columns. If `TRUE`,
  build one window around each column in `subset`, as described under
  Details.

- min_window_n:

  Integer. The smallest number of columns, at least `2`, that a window
  needs to be kept. `k` and `ncp` must be less than `min_window_n`.

- subset:

  Optional character or integer vector. The columns to impute and
  return, by name or by number, or `NULL` for every column. Must be
  supplied when `flank = TRUE`. Defaults to `NULL`.

- dry_run:

  Logical. If `TRUE`, return the windows that `min_window_n` and
  `subset` leave, as a `slideimp_tbl`, without imputing. `k` and `ncp`
  are not needed then.

- k:

  Integer or `NULL`. Number of nearest neighbors for K-NN imputation.

- cores:

  Integer. Number of cores for K-NN imputation. Defaults to `1`.

- dist_pow:

  Numeric. Power used to penalize more distant neighbors in the weighted
  average, where `0` gives an unweighted average. The weighting rule is
  given under Details in
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md).

- ncp:

  Integer or `NULL`. Number of components for PCA imputation.

- scale:

  Logical. If `TRUE`, columns are scaled to unit variance.

- coeff.ridge:

  Numeric. The ridge regularization for `method = "regularized"`,
  ignored otherwise. Values below `1` move toward EM PCA, and values
  above `1` move toward mean imputation.

- threshold:

  Numeric. Convergence threshold.

- seed:

  Integer, numeric, or `NULL`. Random seed for reproducibility. As in
  [`missMDA::imputePCA()`](https://rdrr.io/pkg/missMDA/man/imputePCA.html),
  initialization `i` is seeded with `seed * (i - 1)`. So `seed = 0`
  gives every random initialization the same start, and `nb.init` then
  yields two distinct solves: the mean initialization and one random
  one. With `nb.init > 1`, `seed` can be at most
  `.Machine$integer.max %/% (nb.init - 1)`.

- row.w:

  Numeric vector, `"n_miss"`, or `NULL`. The row weights, normalized to
  sum to `1`: positive values of length `nrow(obj)`, `"n_miss"` to
  down-weight rows with more missing values, or `NULL` for equal
  weights.

- nb.init:

  Integer. Number of initializations. The first is always mean
  imputation, and the others are random.

- maxiter:

  Integer. Maximum number of iterations.

- miniter:

  Integer. Minimum number of iterations, at most `maxiter`.

- solver:

  Character. One of `"auto"`, `"exact"`, or `"lobpcg"`. The eigensolver.
  `"auto"` uses `"exact"` for small problems and otherwise runs a short
  timed probe, picking `"lobpcg"` only when it is clearly faster. With
  `nb.init > 1`, the choice from the first initialization is reused. See
  the PCA Performance tips section.

- lobpcg_control:

  A list or `NULL`. LOBPCG control options, from
  [`lobpcg_control()`](https://hhp94.github.io/slideimp/reference/lobpcg_control.md)
  or as a plain named list, or `NULL` for the defaults. Ignored when
  `solver = "exact"`.

- clamp:

  Numeric vector of length 2 or `NULL`. The bounds `c(lower, upper)` for
  PCA-imputed values, such as `c(0, 1)` for DNA methylation beta values,
  or `NULL` for none. Use `-Inf` or `Inf` for a one-sided bound.
  Observed values are not clamped.

- method:

  Character or `NULL`. One of `"euclidean"` or `"manhattan"` for K-NN
  imputation, or `"regularized"` or `"EM"` for PCA imputation. `NULL`
  uses `"euclidean"` for K-NN and `"regularized"` for PCA.

- .progress:

  Logical. If `TRUE`, show imputation progress.

- colmax:

  Numeric scalar between `0` and `1`. Columns with a proportion of
  missing values greater than `colmax` take no part in the imputation,
  as described under Column mean fill.

- post_imp:

  Logical. If `TRUE`, column means replace the missing values that
  remain after imputation, as described under Column mean fill.

- on_infeasible:

  Character. One of `"skip"`, `"error"`, or `"mean"`. The action for a
  window that is infeasible for imputation, for example because `k` or
  `ncp` is too large for the columns that `colmax` leaves. Defaults to
  `"skip"`.

## Value

If `dry_run = FALSE`, a numeric matrix of class `slideimp_results` with
the rows of `obj` and the columns that `subset` selects, with missing
values imputed.

If `dry_run = TRUE`, a data frame of class `slideimp_tbl` with one row
per window and the columns `start` and `end` (the first and last column
of `obj` in the window), `window_n`, `subset_local` and, when
`flank = TRUE`, `target`.

## Details

`slide_imp()` splits the columns of `obj` into windows by `location` and
imputes each window on its own. There are two windowing modes:

- `flank = FALSE`: greedily partition `location` into windows of width
  `window_size` with the requested `overlap_size` between consecutive
  windows. A window that would add no column to the window before it is
  not created. A column in more than one window takes the mean of their
  imputed values.

- `flank = TRUE`: create one window per feature in `subset`. Each window
  covers the `location` values within `window_size` of that feature, on
  both sides, so it spans twice `window_size`.

Supply exactly one of `k` and `ncp`. `k` imputes each window with
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
and `ncp` imputes it with
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).
The arguments of the other method are ignored.

The result holds the `subset` columns in the order they have in `obj`.
In each window, every column the method can use takes part in the
imputation, in `subset` or not: for K-NN every column within `colmax`,
and for PCA every column within `colmax` whose variance is not zero. A
name that is not in `colnames(obj)`, or a `subset` that selects no
column, is an error.

`on_infeasible = "skip"` does not impute an infeasible window, `"mean"`
fills it with column means, and `"error"` stops at the first such
window. With `"skip"` and `"mean"`, `slide_imp()` raises one warning of
class `slideimp_fallback` that lists the affected windows.

An infeasible window adds nothing to the `mean_filled_cols` and
`mean_filled_cells` attributes. Where windows overlap,
`mean_filled_cells` lists each cell that at least one window filled with
its column mean.

With `solver = "auto"`, the first window that completes the solver probe
picks the solver, and every later window uses it.

`slide_imp()` rejects a column with no observed values only when a
window covers it. A column that no window covers is returned unchanged.

## Column mean fill

A column with a proportion of missing values greater than `colmax` takes
no part in the imputation. The method neither uses the column to impute
other columns nor imputes it.

With `post_imp = TRUE`, column means replace the missing values that
remain after the method has run. Two kinds of missing value can remain:
the missing values of the columns that took no part, and the cells for
which the method found no value. Each column mean is the mean of the
column's observed values. With `post_imp = FALSE`, both kinds stay
missing.

The result records what the column means replaced, in two attributes:

- `mean_filled_cols`: an integer vector with the index, in the result,
  of each column that took no part and was filled.

- `mean_filled_cells`: an integer matrix with the columns `row` and
  `col`, one row for each cell that the method could not impute and that
  was filled.

Both attributes are empty when no column mean was used.
[`print()`](https://rdrr.io/r/base/print.html) reports them when they
are not empty.

## Missing values and non-finite input

`obj` must be a double matrix. An integer or logical matrix is rejected,
not converted. Convert one with `storage.mode(obj) <- "double"`.

`NA` and `NaN` are both treated as missing. The package does not
distinguish between them.

`Inf` and `-Inf` are rejected on entry, and so is a column with no
observed values (every cell `NA` or `NaN`). Locate them with
`which(is.infinite(obj), arr.ind = TRUE)` and
`which(colSums(!is.na(obj)) == 0)`.

## PCA Performance tips

Three arguments set the speed: `solver`, `threshold`, and `scale`. Tune
these first, then the accuracy arguments `ncp` and `coeff.ridge`, on a
representative subset.

**Exact vs. LOBPCG.** `"lobpcg"` suits large, approximately low-rank
matrices with small `ncp`. `"exact"` suits small matrices, such as
`slide_imp()` windows, where it is faster. It is also more reliable: it
always converges, to eigenvectors exact to machine precision, while
LOBPCG stops at the tolerance set in
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

## Examples

``` r
set.seed(1234)

# Example data with 20 samples and 100 ordered columns
beta_matrix <- sim_mat(20, 100)$input
location <- 1:100

# First perform a dry run to inspect the calculated windows
window_statistics <- slide_imp(
  beta_matrix,
  location = location,
  window_size = 50,
  overlap_size = 10,
  min_window_n = 10,
  dry_run = TRUE,
  .progress = FALSE
)
window_statistics
#> # slideimp table: 3 x 4
#>  start end window_n  subset_local
#>      1  50       50 <double [50]>
#>     41  90       50 <double [50]>
#>     81 100       20 <double [20]>

# Sliding-window K-NN imputation
imputed_knn <- slide_imp(
  beta_matrix,
  location = location,
  k = 5,
  window_size = 50,
  overlap_size = 10,
  min_window_n = 10,
  .progress = FALSE
)
imputed_knn
#> Method: slide_imp (K-NN imputation)
#> Dimensions: 20 x 100
#> 
#>          feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.3486482 0.7385414 0.4077444 0.1607935 0.3924661 0.2434143
#> sample2 0.5338935 0.4724364 0.9663621 0.4922511 0.5061132 0.3923603
#> sample3 0.7185848 0.7351035 0.6724479 0.3162537 0.7634236 1.0000000
#> sample4 0.1734418 0.0000000 0.0000000 0.0000000 0.0000000 0.2106358
#> sample5 0.4499620 0.5306182 0.5685354 0.5383513 0.4680080 0.8518388
#> sample6 0.3768380 0.5570723 0.8764909 0.5276245 0.6722794 0.5740639
#> # Showing 6 of 20 rows and 6 of 100 columns

# Sliding-window PCA imputation
imputed_pca <- slide_imp(
  beta_matrix,
  location = location,
  ncp = 2,
  window_size = 50,
  overlap_size = 10,
  min_window_n = 10,
  .progress = FALSE
)
imputed_pca
#> Method: slide_imp (PCA imputation)
#> Dimensions: 20 x 100
#> 
#>          feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.3486482 0.7385414 0.4077444 0.1607935 0.3924661 0.2434143
#> sample2 0.5338935 0.4724364 0.9663621 0.4242777 0.5061132 0.3923603
#> sample3 0.7185848 0.7351035 0.6724479 0.3162537 0.7634236 1.0000000
#> sample4 0.1734418 0.0000000 0.0000000 0.0000000 0.0000000 0.2106358
#> sample5 0.4971181 0.5306182 0.5685354 0.5383513 0.4680080 0.8518388
#> sample6 0.3768380 0.5570723 0.8764909 0.5276245 0.6722794 0.5740639
#> # Showing 6 of 20 rows and 6 of 100 columns

# K-NN imputation with flanking windows
imputed_flank <- slide_imp(
  beta_matrix,
  location = location,
  k = 2,
  window_size = 30,
  flank = TRUE,
  subset = c(10, 30, 70),
  min_window_n = 5,
  .progress = FALSE
)
imputed_flank
#> Method: slide_imp (K-NN imputation)
#> Dimensions: 20 x 3
#> Note: 1 cell could not be imputed. Column means replaced it.
#>   See attribute "mean_filled_cells".
#> 
#>         feature10  feature30 feature70
#> sample1 0.2539968 0.23569167 0.4917902
#> sample2 0.5898552 0.74940536 0.3084749
#> sample3 0.6484226 0.56464053 1.0000000
#> sample4 0.0711946 0.04345129 0.1101074
#> sample5 0.5436376 0.62011849 0.7986741
#> sample6 0.6477420 0.58663423 0.7320680
#> # Showing 6 of 20 rows and 3 of 3 columns
```
