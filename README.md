
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SlideKnn

<!-- badges: start -->

<!-- badges: end -->

`{SlideKnn}` is an efficient R package for k-nearest neighbors (k-NN)
imputation of missing values in high-dimensional numeric matrices, such
as those from intensive longitudinal data or epi-genomics. It introduces
a sliding window approach to handle very large feature sets while
preserving local structure, making it suitable for data where features
are ordered (e.g., by time or distance).

The package builds on proven k-NN concepts (Bioconductor’s `impute`
package) but adds enhancements: parallelization for speed, tree-based
methods for efficiency, weighted imputation, multiple imputation
strategies, and built-in tuning tools. It’s designed for matrices with
samples in rows and features in columns.

Key features include:

- **Sliding Window k-NN Imputation**: Break large data into overlapping
  windows for computationally feasible imputation while maintaining
  local structures (e.g., intensively sampled longitudinal data or
  epi-genomic data).
- **Full Matrix k-NN Imputation**: Standard k-NN for smaller data, with
  multi-core parallelization over columns with missing values.
- **Tree-Based k-NN**: Integration with `{mlpack}` for KD-Tree or
  Ball-Tree methods, accelerating imputation in high dimensions.
- **Subset Imputation**: Only impute a subset of columns to save time.
  Important for applications such as epi-genetic clocks calculations.
- **Weighted Imputation**: Use inverse-distance weighting for more
  accurate averages, with tunable penalties.
- **Multiple Imputation**: Support for Predictive Mean Matching (PMM) or
  bootstrap resampling from nearest neighbors.
- **Fallback Imputation**: Optional post-k-NN mean imputation by column
  to handle remaining NAs.
- **Parameter Tuning**: Inject artificial NAs to evaluate and tune
  hyperparameters, with support for custom imputation functions.
- **Big Matrix Support**: Compatible with `{bigmemory}` for file-backed
  matrices to handle massive data without loading everything into RAM.

## Installation

You can install the development version of `{SlideKnn}` from
[GitHub](https://github.com/hhp94/SlideKnn) with:

``` r
# install.packages("remotes")
remotes::install_github("hhp94/SlideKnn")
```

## Example

Load the package and use the built-in `khanmiss1` data (see `?khanmiss1`
for details). For full matrix k-NN imputation without sliding windows:

``` r
library(SlideKnn)

data(khanmiss1)
# Transpose for samples in rows, features in columns
imputed_full <- knn_imp(t(khanmiss1), k = 3, method = "euclidean", cores = 1)

sum(is.na(imputed_full[[1]]))
#> [1] 0
```

Yield same results as `{impute::impute.knn()}`. Is faster in larger data
and almost as fast in smaller data. Scaling well with multiple `cores`.

``` r
set.seed(1234)
obj_t <- khanmiss1

bench::mark(
  # Single Core
  knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "euclidean")[[1]],
  # Multiple Cores
  knn_imp(t(khanmiss1), k = 3, rowmax = 1, cores = 4, method = "euclidean")[[1]],
  t(impute::impute.knn(khanmiss1, k = 3, rowmax = 1, maxp = nrow(khanmiss1))$data)
) |>
  dplyr::mutate(expression = c("knn_1", "knn_4", "impute.knn_1")) |>
  dplyr::select(-dplyr::where(is.list))
#> # A tibble: 3 × 6
#>   expression        min   median `itr/sec` mem_alloc `gc/sec`
#>   <chr>        <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#> 1 knn_1         15.13ms  15.94ms      63.0    9.47MB     30.0
#> 2 knn_4          5.76ms   6.11ms     162.     9.47MB     75.2
#> 3 impute.knn_1  14.62ms  14.77ms      67.1   12.42MB     28.8
```

Sliding window k-NN imputation for epi-genetics data with 1000 CpGs and
10 Samples

``` r
set.seed(1234)
beta_matrix <- t(sim_mat(m = 10, n = 1000, perc_NA = 0.3, perc_col_NA = 1)$input)
beta_matrix[1:10, 1:5]
#>         feat1     feat2     feat3     feat4     feat5
#> s1  0.3320706 0.5572683 0.6796919 0.1593403 0.5802804
#> s2  0.3046741 0.5442470 0.2515999 0.5973359 0.6080809
#> s3  0.3467045        NA        NA 0.7007942 0.2352033
#> s4         NA 0.5771749 0.5437410 0.7921855 0.7523102
#> s5  0.3395358 0.3664415        NA        NA 0.6398042
#> s6  0.3144117        NA 0.2995698 0.3984796 0.3572444
#> s7  0.4963455 0.4612003 0.3311917        NA        NA
#> s8  0.5039341 0.5381084 0.3811143 0.6207009 0.6710451
#> s9         NA        NA 0.3455449        NA        NA
#> s10        NA 0.5498897        NA 0.3738964        NA

imputed <- SlideKnn(beta_matrix, n_feat = 500, n_overlap = 10, k = 10)
sum(is.na(imputed[[1]]))
#> [1] 0
```

Using tree-based k-NN with weighting and multiple imputation (PMM):

``` r
imputed_tree <- knn_imp(
  t(khanmiss1),
  k = 5,
  tree = "kd", # KD-Tree via mlpack
  weighted = TRUE, # Inverse-distance weighting
  dist_pow = 2, # Harsher penalty for distant neighbors
  n_imp = 3, # 3 imputations
  n_pmm = 10 # PMM with 10 donors
)
#> [INFO ] 55839 node combinations were scored.
#> [INFO ] 480166 base cases were calculated.
length(imputed_tree) # 3 imputed matrices
#> [1] 3
```

Optional simple mean imputation as a fallback or baseline that enables
multi-step imputation strategies

``` r
# Inject extra NA into simulated data to make k-NN fails for first imputation
set.seed(1234)
obj <- t(sim_mat(n = 1000, m = 100, perc_NA = 0.8, perc_col_NA = 1)$input)
# Step 1: Column-wise imputation. Disable fall back imputation with `post_imp == FALSE`
imputed_by_col <- knn_imp(obj, cores = 4, k = 10, post_imp = FALSE)
# Step 2: Row-wise imputation. Then if values are still missing, impute by rows
imputed_by_row <- knn_imp(t(imputed_by_col[[1]]), cores = 4, k = 10, post_imp = FALSE)
# Step 3: Lastly, impute by column mean
imputed_mean <- mean_impute_col(t(imputed_by_row[[1]]))
sum(is.na(imputed_mean))
#> [1] 0
```

For very large matrices that don’t fit in memory, use `{bigmemory}` to
create file-backed objects. `SlideKnn` supports passing a `big.matrix`
object or the path to its descriptor file. Always specify `output` for
the result (also file-backed).

``` r
library(bigmemory)

data(khanmiss1)
mat <- t(khanmiss1) # samples rows, features cols

temp_dir <- withr::local_tempdir()

# Create big.matrix with backing and descriptor files
big_mat <- bigmemory::as.big.matrix(
  mat,
  type = "double",
  backingfile = "khan.bin",
  descriptorfile = "khan.desc",
  backingpath = temp_dir
)

# Impute using the big.matrix object (returns list of big.matrix)
imputed_obj <- SlideKnn(
  obj = big_mat,
  n_feat = 100,
  n_overlap = 10,
  k = 10,
  output = file.path(temp_dir, "imputed.bin")
)
sum(is.na(imputed_obj[[1]][, ])) # Access the big.matrix result
#> [1] 0

# Alternatively, impute using the descriptor path
desc_path <- file.path(temp_dir, "khan.desc")
imputed_path <- SlideKnn(
  obj = desc_path,
  n_feat = 100,
  n_overlap = 10,
  k = 10,
  output = file.path(temp_dir, "imputed.bin")
)
sum(is.na(imputed_path[[1]][, ]))
#> [1] 0
```

Note: Results are lists of `big.matrix` objects (one per `n_imp`). Use
`bigmemory::attach.big.matrix()` if needed to reload from descriptor
later. Set `strip_dimnames = TRUE` for multi-core efficiency, and
restore dimnames post-imputation if required (same dimnames as `obj`).

## Parameter Tuning

Use `tune_imp()` to tune hyperparameters by injecting artificial NAs and
evaluating imputation accuracy. It supports built-in methods or custom
functions.

``` r
parameters <- dplyr::tibble(
  k = c(5, 10),
  method = "euclidean",
  weighted = TRUE,
  dist_pow = c(2, 5),
  post_imp = TRUE
)
parameters
#> # A tibble: 2 × 5
#>       k method    weighted dist_pow post_imp
#>   <dbl> <chr>     <lgl>       <dbl> <lgl>   
#> 1     5 euclidean TRUE            2 TRUE    
#> 2    10 euclidean TRUE            5 TRUE

# Subset for faster tuning
obj <- t(khanmiss1)

# 3 repeats, each time inject 50 NAs
results <- tune_imp(obj, parameters, rep = 3, .f = "knn_imp", num_na = 100)

# Compute metrics with {yardstick}
library(yardstick)
met_set <- metric_set(mae, rmse, rsq)
results$metrics <- lapply(results$result, function(x) met_set(x, truth = truth, estimate = estimate))
head(
  dplyr::select(
    tidyr::unnest(dplyr::select(results, -result), cols = "metrics"),
    all_of(names(parameters)), contains(".")
  )
)
#> # A tibble: 6 × 8
#>       k method    weighted dist_pow post_imp .metric .estimator .estimate
#>   <dbl> <chr>     <lgl>       <dbl> <lgl>    <chr>   <chr>          <dbl>
#> 1     5 euclidean TRUE            2 TRUE     mae     standard     344.   
#> 2     5 euclidean TRUE            2 TRUE     rmse    standard     492.   
#> 3     5 euclidean TRUE            2 TRUE     rsq     standard       0.390
#> 4    10 euclidean TRUE            5 TRUE     mae     standard     354.   
#> 5    10 euclidean TRUE            5 TRUE     rmse    standard     495.   
#> 6    10 euclidean TRUE            5 TRUE     rsq     standard       0.380
```

Tuning a custom imputation function:

``` r
# This custom function impute NAs with rnorm values
custom_imp <- function(obj, mean = 0, sd = 1) {
  na_pos <- is.na(obj)
  obj[na_pos] <- rnorm(sum(na_pos), mean = mean, sd = sd)
  return(obj)
}

parameters_custom <- dplyr::tibble(
  mean = c(0, 1),
  sd = c(1, 2)
)

results_custom <- tune_imp(obj, parameters_custom, .f = custom_imp, rep = 2, num_na = 20)
results_custom$metrics <- lapply(results_custom$result, function(x) met_set(x, truth = truth, estimate = estimate))
head(
  dplyr::select(
    tidyr::unnest(dplyr::select(results_custom, -result), cols = "metrics"),
    dplyr::all_of(names(parameters_custom)), dplyr::contains(".")
  )
)
#> # A tibble: 6 × 5
#>    mean    sd .metric .estimator  .estimate
#>   <dbl> <dbl> <chr>   <chr>           <dbl>
#> 1     0     1 mae     standard   1038.     
#> 2     0     1 rmse    standard   1238.     
#> 3     0     1 rsq     standard      0.00312
#> 4     1     2 mae     standard   1037.     
#> 5     1     2 rmse    standard   1238.     
#> 6     1     2 rsq     standard      0.0810
```

For more details, see the function documentation (e.g., `?SlideKnn`,
`?knn_imp`, `?tune_imp`).
