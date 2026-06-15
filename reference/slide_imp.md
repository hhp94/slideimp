# Sliding-Window K-NN or PCA Imputation

Perform sliding-window K-NN or PCA imputation on a numeric matrix whose
columns are meaningfully ordered. Not intended for Illumina DNA
methylation microarrays, use
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
  na_check = TRUE,
  on_infeasible = c("skip", "error", "mean")
)
```

## Arguments

- obj:

  A numeric matrix with samples in rows and features in columns.

- location:

  A sorted numeric vector of length `ncol(obj)` giving the position of
  each column, such as genomic coordinates.

- window_size:

  Numeric. Window width in the same units as `location`.

- overlap_size:

  Numeric. Overlap between consecutive windows in the same units as
  `location`. Must be less than `window_size`. Ignored when
  `flank = TRUE`.

- flank:

  Logical. If `FALSE`, imputation uses sliding windows across the full
  matrix. If `TRUE`, one window of width `window_size` is created for
  each feature listed in `subset`; `overlap_size` is ignored. `subset`
  must be supplied when `flank = TRUE`.

- min_window_n:

  Integer. Minimum number of columns a window must contain to be
  considered for imputation. For non-dry runs, the selected `k` or `ncp`
  value must be feasible for the usable columns in each retained window
  after applying `colmax`.

- subset:

  Optional character or integer vector specifying columns to impute. If
  `NULL`, all eligible columns are imputed. Required when
  `flank = TRUE`.

- dry_run:

  Logical. If `TRUE`, skip imputation and return a `slideimp_tbl`
  describing the windows that would be used after all filtering rules
  are applied. In this mode, `k` and `ncp` are not required.

- k:

  Integer. Number of nearest neighbors to use for K-NN imputation.

- cores:

  Integer. Number of cores to use for K-NN imputation.

- dist_pow:

  Numeric. Power used to penalize more distant neighbors in the weighted
  average. `dist_pow = 0` gives an unweighted average of the nearest
  neighbors.

- ncp:

  Integer. Number of principal components used to predict missing
  entries.

- scale:

  Logical. If `TRUE`, columns are scaled to unit variance.

- coeff.ridge:

  Numeric. Ridge regularization, used only when
  `method = "regularized"`. Values `< 1` move toward EM PCA; values
  `> 1` move toward mean imputation.

- threshold:

  Numeric. Convergence threshold.

- seed:

  Integer, numeric, or `NULL`. Random seed for reproducibility.

- row.w:

  Row weights, normalized to sum to `1`. `NULL` (equal weights), a
  positive numeric vector of length `nrow(obj)`, or `"n_miss"`
  (down-weight rows with more missing values).

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations.

- miniter:

  Integer. Minimum number of iterations.

- solver:

  Character. Eigensolver: `"auto"` (default), `"exact"`, or `"lobpcg"`.
  `"auto"` runs a short timed probe and picks `"lobpcg"` only when
  clearly faster. Consecutive EM calls warm-start LOBPCG with both the
  previous eigenblock and search direction. When `nb.init > 1`, the auto
  choice from the first init is reused. See Performance tips.

- lobpcg_control:

  A list of LOBPCG eigensolver control options, usually created by
  [`lobpcg_control()`](https://hhp94.github.io/slideimp/reference/lobpcg_control.md).
  A plain named list is also accepted. Ignored when `solver = "exact"`.

- clamp:

  Optional numeric vector `c(lower, upper)` bounding PCA-imputed values
  (use `-Inf`/`Inf` for one-sided, `NULL` for none). E.g., `c(0, 1)` for
  DNAm beta values. Observed values are not clamped.

- method:

  Character or `NULL`. For K-NN imputation, one of `"euclidean"` or
  `"manhattan"`. For PCA imputation, one of `"regularized"` or `"EM"`.
  If `NULL`, the corresponding backend default is used.

- .progress:

  Logical. If `TRUE`, show imputation progress.

- colmax:

  Numeric scalar between `0` and `1`. Columns with a missing-data
  proportion greater than `colmax` are excluded from the main imputation
  method. Excluded columns are left unchanged unless `post_imp = TRUE`,
  in which case remaining missing values are replaced by column means
  when possible.

- post_imp:

  Logical. If `TRUE`, replace missing values remaining after the main
  imputation method with column means when possible.

- na_check:

  Logical. If `TRUE`, check whether the returned matrix still contains
  missing values.

- on_infeasible:

  Character. One of `"skip"`, `"error"`, or `"mean"`. Controls behavior
  when a window is infeasible for imputation, for example when `k` or
  `ncp` exceeds the number of usable columns after applying `colmax`.

## Value

If `dry_run = FALSE`, a numeric matrix of the same dimensions as `obj`,
with missing values imputed. The returned object has class
`slideimp_results`.

If `dry_run = TRUE`, a data frame of class `slideimp_tbl` with columns
`start`, `end`, and `window_n`, plus `subset_local` and, when
`flank = TRUE`, `target`.

## Details

The sliding-window approach divides the input matrix into smaller
segments based on `location` values and applies imputation to each
window independently. Values in overlapping regions are averaged across
windows to produce the final imputed result.

Two windowing modes are supported:

- `flank = FALSE`: greedily partition `location` into windows of width
  `window_size` with the requested `overlap_size` between consecutive
  windows.

- `flank = TRUE`: create one window per feature in `subset`, centered on
  that feature using the supplied `window_size`.

Specify `k` and related arguments to use
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md), or
`ncp` and related arguments to use
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).

## PCA Performance tips

Speed comes from three levers: `solver` (through LOBPCG with
warm-start), `threshold`, and `scale`. Tune these first, then accuracy
parameters (`ncp`, `coeff.ridge`) on a representative subset.

**Exact vs. LOBPCG with warm-start.** Whether `"lobpcg"` beats `"exact"`
depends on size and low-rankness: prefer `"lobpcg"` for large,
approximately low-rank matrices with small `ncp`, and `"exact"` for
small matrices (including `slide_imp()` windows), where it is faster and
more robust. Separately, the warm-start makes each successive solve
cheap:
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
warm-starts LOBPCG with the previous eigenblock and search direction, so
once imputed values stabilize, later solves converge in a few
iterations. The payoff therefore grows with the number of EM iterations,
independent of low-rankness. `solver = "auto"` (default) probes both and
is a safe start.

**Threshold.** The default `1e-6` is conservative; `1e-5` is often
faster with very similar values.

**Scale.** For columns on a common scale (e.g., DNAm beta values in
`[0, 1]`), `scale = FALSE` can be faster and more accurate.

**Parallel and BLAS.** In parallel via
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
or
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
with a multithreaded BLAS, set `pin_blas = TRUE` to avoid thread
oversubscription. On Windows, the stock BLAS can be slow. Advanced users
can swap in
[OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).

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
#> Method: slide_imp (KNN imputation)
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
#> Method: slide_imp (KNN imputation)
#> Dimensions: 20 x 100
#> 
#>          feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.3486482 0.7385414 0.4077444 0.1607935 0.3924661 0.2434143
#> sample2 0.5338935 0.4724364 0.9663621        NA 0.5061132 0.3923603
#> sample3 0.7185848 0.7351035 0.6724479 0.3162537 0.7634236 1.0000000
#> sample4 0.1734418 0.0000000 0.0000000 0.0000000 0.0000000 0.2106358
#> sample5        NA 0.5306182 0.5685354 0.5383513 0.4680080 0.8518388
#> sample6 0.3768380 0.5570723 0.8764909 0.5276245 0.6722794 0.5740639
#> # Showing 6 of 20 rows and 6 of 100 columns
```
