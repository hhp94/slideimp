# Sliding Window K-NN or PCA Imputation

Performs sliding window K-NN or PCA imputation of large numeric matrices
column-wise. This method assumes that columns are meaningfully sorted by
`location`.

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
  max_cache = 4,
  ncp = NULL,
  scale = TRUE,
  coeff.ridge = 1,
  seed = NULL,
  row.w = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
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

  A numeric matrix with **samples in rows** and **features in columns**.

- location:

  A sorted numeric vector of length `ncol(obj)` giving the position of
  each column (e.g., genomic coordinates). Used to define sliding
  windows.

- window_size:

  Window width in the same units as `location`.

- overlap_size:

  Overlap between consecutive windows in the same units as `location`.
  Must be less than `window_size`. Default is `0`. Ignored when
  `flank = TRUE`.

- flank:

  Logical. If `TRUE`, instead of sliding windows across the whole
  matrix, one window of width `window_size` is created flanking each
  feature listed in `subset`. In this mode `overlap_size` is ignored.
  Requires `subset` to be provided. Default = `FALSE`.

- min_window_n:

  Minimum number of columns a window must contain to be imputed. Windows
  smaller than this are not imputed. `k` and `ncp` must also be smaller
  than `min_window_n`.

- subset:

  Character. Vector of column names or integer vector of column indices
  specifying which columns to impute.

- dry_run:

  Logical. If `TRUE`, skip imputation and return a `slideimp_tbl` object
  of the windows that *would* be used (after all dropping rules).
  `k`/`ncp` are not required in this mode. Columns: `start`, `end`,
  `window_n`, plus `subset_local` (list-column of local subset indices)
  when `flank = FALSE`, or `target` and `subset_local` when
  `flank = TRUE`. Default = `FALSE`.

- k:

  Integer. Number of nearest neighbors for imputation. 10 is a good
  starting point.

- cores:

  Integer. Number of cores for K-NN parallelization (OpenMP). On macOS,
  OpenMP may need additional compiler configuration.

- dist_pow:

  Numeric. The amount of penalization for further away nearest neighbors
  in the weighted average. `dist_pow = 0` (default) is the simple
  average of the nearest neighbors.

- max_cache:

  Numeric. Maximum allowed cache size in GB (default `4`). When greater
  than `0`, pairwise distances between columns with missing values are
  pre-computed and cached, which is faster for moderate-sized data but
  uses O(m^2) memory where m is the number of columns with missing
  values. Set to `0` to disable caching and trade speed for lower memory
  usage.

- ncp:

  Integer. Number of components used to predict the missing entries.

- scale:

  Logical. If `TRUE` (default), variables are scaled to have unit
  variance.

- coeff.ridge:

  Numeric. Ridge regularization coefficient (default is 1). Only used if
  `method = "regularized"`. Values \< 1 regularize less (closer to EM);
  values \> 1 regularize more (closer to mean imputation).

- seed:

  Integer. Random number generator seed.

- row.w:

  Row weights (internally normalized to sum to 1). Can be one of:

  - `NULL` (default): All rows weighted equally.

  - A numeric vector: Custom positive weights of length `nrow(obj)`.

  - `"n_miss"`: Rows with more missing values receive lower weight.

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations for the algorithm.

- miniter:

  Integer. Minimum number of iterations for the algorithm.

- method:

  For K-NN imputation: distance metric to use (`"euclidean"` or
  `"manhattan"`). For PCA imputation: regularization imputation
  algorithm (`"regularized"` or `"EM"`).

- .progress:

  Show progress bar (default = `TRUE`).

- colmax:

  Numeric. A number from 0 to 1. Threshold of column-wise missing data
  rate above which imputation is skipped.

- post_imp:

  Boolean. Whether to impute remaining missing values (those that failed
  imputation) using column means.

- na_check:

  Boolean. Check for leftover `NA` values in the results or not
  (internal use).

- on_infeasible:

  Character, one of `"error"` (default on
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)),
  `"skip"`, or `"mean"` (default on `slide_imp()`). Controls behaviour
  when a group is infeasible for imputation, e.g., `k`/`ncp` exceeds the
  number of usable columns after applying `colmax`, or all subset
  columns in the group exceed `colmax`.

## Value

A numeric matrix of the same dimensions as `obj` with missing values
imputed. When `dry_run = TRUE`, returns a `data.frame` of class
`slideimp_tbl` with columns `start`, `end`, `window_n`, plus
`subset_local` (and `target` when `flank = TRUE`).

## Details

The sliding window approach divides the input matrix into smaller
segments based on `location` values and applies imputation to each
window independently. Values in overlapping areas are averaged across
windows to produce the final imputed result.

Two windowing modes are supported:

- `flank = FALSE` (default): Greedily partitions the `location` vector
  into windows of width `window_size` with the requested `overlap_size`
  between consecutive windows.

- `flank = TRUE`: Creates one window per feature in `subset` that
  exactly flanks that specific feature using the supplied `window_size`.

Specify `k` and related arguments to use
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md),
`ncp` and related arguments for
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).

## Examples

``` r
# Generate sample data with missing values with 20 samples and 100 columns
# where the column order is sorted (i.e., by genomic position)
set.seed(1234)
beta_matrix <- sim_mat(20, 100)$input
location <- 1:100

# It's very useful to first perform a dry run to examine the calculated windows
windows_statistics <- slide_imp(
  beta_matrix,
  location = location,
  window_size = 50,
  overlap_size = 10,
  min_window_n = 10,
  dry_run = TRUE
)
windows_statistics
#> # slideimp table: 3 x 4
#>  start end window_n  subset_local
#>      1  50       50 <double [50]>
#>     41  90       50 <double [50]>
#>     81 100       20 <double [20]>

# Sliding Window K-NN imputation by specifying `k` (sliding windows)
imputed_knn <- slide_imp(
  beta_matrix,
  location = location,
  k = 5,
  window_size = 50,
  overlap_size = 10,
  min_window_n = 10,
  scale = FALSE # This argument belongs to PCA imputation and will be ignored
)
#> Step 1/2: Imputing
#>  Processing window 1 of 3
#>  Processing window 2 of 3
#>  Processing window 3 of 3
#> Step 2/2: Averaging overlapping regions
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
#> 
#> # Showing [1:6, 1:6] of full matrix

# Sliding Window PCA imputation by specifying `ncp` (sliding windows)
pca_knn <- slide_imp(
  beta_matrix,
  location = location,
  ncp = 2,
  window_size = 50,
  overlap_size = 10,
  min_window_n = 10
)
#> Step 1/2: Imputing
#>  Processing window 1 of 3
#>  Processing window 2 of 3
#>  Processing window 3 of 3
#> Step 2/2: Averaging overlapping regions
pca_knn
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
#> 
#> # Showing [1:6, 1:6] of full matrix

# Sliding Window K-NN imputation with flanking windows (flank = TRUE)
# Only the columns listed in `subset` are imputed; each uses its own
# centered window of width `window_size`.
imputed_flank <- slide_imp(
  beta_matrix,
  location = location,
  k = 2,
  window_size = 30,
  flank = TRUE,
  subset = c(10, 30, 70),
  min_window_n = 5
)
#> Step 1/2: Imputing
#>  Processing window 1 of 3
#>  Processing window 2 of 3
#>  Processing window 3 of 3
#> Step 2/2: Averaging overlapping regions
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
#> 
#> # Showing [1:6, 1:6] of full matrix
```
