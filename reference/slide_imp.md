# Sliding-Window K-NN or PCA Imputation

Perform sliding-window K-NN or PCA imputation on a numeric matrix whose
columns are meaningfully ordered.

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
  solver = c("auto", "dsyevr", "lobpcg"),
  lobpcg_control = NULL,
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
  imputed. `k` and `ncp` must be smaller than `min_window_n`.

- subset:

  Optional character or integer vector specifying columns to impute. If
  `NULL`, all eligible columns are imputed. Required when
  `flank = TRUE`.

- dry_run:

  Logical. If `TRUE`, skip imputation and return a `slideimp_tbl`
  describing the windows that would be used after all filtering rules
  are applied. In this mode, `k` and `ncp` are not required.

- k:

  Integer or `NULL`. Number of nearest neighbors for K-NN imputation.
  Supply `k` to use K-NN imputation.

- cores:

  Integer. Number of cores to use for K-NN imputation.

- dist_pow:

  Numeric. Power used to penalize more distant neighbors in the weighted
  average. `dist_pow = 0` gives an unweighted average of the nearest
  neighbors.

- ncp:

  Integer or `NULL`. Number of components for PCA imputation. Supply
  `ncp` to use PCA imputation.

- scale:

  Logical. If `TRUE`, columns are scaled to unit variance.

- coeff.ridge:

  Numeric. Ridge regularization coefficient. Only used when
  `method = "regularized"`. Values less than `1` regularize less, moving
  closer to EM PCA. Values greater than `1` regularize more, moving
  closer to mean imputation.

- threshold:

  Numeric. Convergence threshold.

- seed:

  Integer, numeric, or `NULL`. Random seed for reproducibility.

- row.w:

  Row weights, internally normalized to sum to `1`. Can be:

  - `NULL`: all rows are weighted equally.

  - A numeric vector of positive weights with length `nrow(obj)`.

  - `"n_miss"`: rows with more missing values receive lower weight.

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations.

- miniter:

  Integer. Minimum number of iterations.

- solver:

  Character. Eigensolver selection. One of `"auto"`, `"dsyevr"`, or
  `"lobpcg"`. `"dsyevr"` uses the exact solver. `"lobpcg"` uses the
  iterative LOBPCG solver. If `"auto"`, LOBPCG is used when the smaller
  input dimension is at least 500 and `ncp <= 50`; otherwise the exact
  solver is used.

- lobpcg_control:

  A list of LOBPCG eigensolver control options, usually created by
  [`lobpcg_control()`](https://hhp94.github.io/slideimp/reference/lobpcg_control.md).
  A plain named list is also accepted. If `NULL`, defaults are chosen
  according to `solver`.

- method:

  Character or `NULL`. For K-NN imputation, one of `"euclidean"` or
  `"manhattan"`. For PCA imputation, one of `"regularized"` or `"EM"`.
  If `NULL`, the corresponding backend default is used.

- .progress:

  Logical. If `TRUE`, show progress.

- colmax:

  Numeric scalar between `0` and `1`. Columns with a missing-data
  proportion greater than `colmax` are not imputed.

- post_imp:

  Logical. If `TRUE`, replace any remaining missing values with column
  means after imputation.

- na_check:

  Logical. If `TRUE`, check whether the result still contains missing
  values.

- on_infeasible:

  Character. One of `"skip"`, `"error"`, or `"mean"`. Controls behavior
  when a window is infeasible for imputation, for example when `k` or
  `ncp` exceeds the number of usable columns after applying `colmax`.

## Value

A numeric matrix of the same dimensions as `obj`, with missing values
imputed. The returned object has class `slideimp_results`.

When `dry_run = TRUE`, returns a `data.frame` of class `slideimp_tbl`
with columns `start`, `end`, and `window_n`, plus `subset_local` and,
when `flank = TRUE`, `target`.

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

## Performance tips

[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
relies heavily on linear algebra. On Windows, the default BLAS shipped
with R may be slow for large matrices. Advanced users can replace it
with [OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).

PCA imputation speed depends on the eigensolver selected by `solver` and
the convergence threshold `threshold`. The exact solver is selected with
`solver = "dsyevr"`. The iterative LOBPCG solver is selected with
`solver = "lobpcg"`. The default, `solver = "auto"`, uses a conservative
internal rule.

For large or approximately low-rank genomic matrices, it can be useful
to benchmark `solver = "dsyevr"` against `solver = "lobpcg"` on a
representative subset, such as chromosome 22, before tuning
accuracy-related parameters such as `ncp`, `coeff.ridge`, `window_size`,
or `overlap_size`.

The default `threshold = 1e-6` is conservative. In some genomic
datasets, `threshold = 1e-5` can be faster while giving very similar
imputed values. Check this on a representative subset before using the
relaxed threshold in a full analysis.

See the pkgdown article "Speeding up PCA imputation" for a full
workflow.

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
