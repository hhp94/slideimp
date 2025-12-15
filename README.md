
<!-- README.md is generated from README.Rmd. Please edit that file -->

# slideimp

<!-- badges: start -->

[![R-CMD-check](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hhp94/slideimp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`{slideimp}` is a lightweight R package for fast K-NN and PCA imputation
of missing values in high-dimensional numeric matrices (such as
intensive longitudinal or epigenetic data).

**Core functions**

- `knn_imp()`: Full-matrix K-NN imputation with multi-core
  parallelization, [`{mlpack}`](https://mlpack.org/) KD/Ball-Tree
  nearest neighbor implementation (for data with very low missing rates
  and extremely high dimensions), and optional subset imputation (ideal
  for epigenetic clock calculations).
- `pca_imp()`: Optimized version of
  [`missMDA::imputePCA()`](http://factominer.free.fr/missMDA/PCA.html)
  for high-dimensional numeric matrices.
- `slide_imp()`: Sliding window K-NN or PCA imputation for extremely
  high-dimensional numeric matrices with ordered features (by time or
  genomic position).
- `group_imp()`: Parallelizable group-wise (e.g., by chromosomes or
  column clusters) K-NN or PCA imputation with optional auxiliary
  features and group-wise parameters.
- `tune_imp()`: Parallelizable hyperparameter tuning with
  cross-validation; works with built-in or custom imputation functions.

## Installation

The stable version of `{slideimp}` can be installed from CRAN using:

``` r
install.packages("slideimp")
```

You can install the development version of `{slideimp}` from
[GitHub](https://github.com/hhp94/slideimp) with:

``` r
install.packages("remotes")
remotes::install_github("hhp94/slideimp")
```

## Workflow

Let’s simulate some DNA methylation (DNAm) microarray data from 2
chromosomes. All `{slideimp}` functions expect the input to be numeric
matrices where variables are stored in the columns.

``` r
library(slideimp)
# Simulate data from 2 chromosomes
set.seed(1234)
sim_obj <- sim_mat(m = 20, n = 50, perc_NA = 0.3, perc_col_NA = 1, nchr = 2)
# Here we see that variables are stored in rows
sim_obj$input[1:5, 1:5]
#>              s1        s2        s3        s4        s5
#> feat1 0.2391314 0.0000000 0.5897476 0.4201222        NA
#> feat2        NA 0.2810446 0.3677927        NA 0.6387734
#> feat3 0.7203854 0.1600776 0.5027545        NA 0.5556735
#> feat4 0.0000000 0.1816453 0.3608640 0.3356484 0.6394179
#> feat5 0.5827582 0.3774313 0.2801131 0.5047049 0.5761809

# So we t() to put the variables in columns
obj <- t(sim_obj$input)
```

Estimate the prediction accuracy of different methods and tune
hyperparameters. Use the .f argument for custom function.

For custom functions, the `parameters` data.frame must include the
columns corresponding to the arguments passed to the custom function.
The custom function must accept `obj` as the first argument and return
an object with the same dimensions.

We tune the results using 2 repeats (`rep = 2`) for illustration
(increase in actual analyses).

``` r
knn_params <- tibble::tibble(k = c(5, 20))
tune_knn <- tune_imp(obj, parameters = knn_params, cores = 2, rep = 2)
#> Tuning knn_imp
#> Step 1/2: Injecting NA
#> Step 2/2: Tuning
compute_metrics(tune_knn)
#> # A tibble: 12 × 7
#>        k cores param_set   rep .metric .estimator .estimate
#>    <dbl> <dbl>     <int> <int> <chr>   <chr>          <dbl>
#>  1     5     2         1     1 mae     standard     0.178  
#>  2     5     2         1     1 rmse    standard     0.225  
#>  3     5     2         1     1 rsq     standard     0.00454
#>  4    20     2         2     1 mae     standard     0.149  
#>  5    20     2         2     1 rmse    standard     0.190  
#>  6    20     2         2     1 rsq     standard     0.0172 
#>  7     5     2         1     2 mae     standard     0.202  
#>  8     5     2         1     2 rmse    standard     0.259  
#>  9     5     2         1     2 rsq     standard     0.00960
#> 10    20     2         2     2 mae     standard     0.172  
#> 11    20     2         2     2 rmse    standard     0.219  
#> 12    20     2         2     2 rsq     standard     0.0850
```

For PCA and custom functions, setup parallelization with
`mirai::daemons(n_cores)`.

``` r
mirai::daemons(2)
pca_params <- tibble::tibble(ncp = c(1, 5))
tune_pca <- tune_imp(obj, parameters = pca_params, cores = 2, rep = 2)

# The parameters have `mean` and `sd` columns.
custom_params <- tibble::tibble(mean = 1, sd = 0)
# This function impute data with rnorm value of different `mean` and `sd`.
custom_function <- function(obj, mean, sd) {
  missing <- is.na(obj)
  obj[missing] <- rnorm(sum(missing), mean = mean, sd = sd)
  return(obj)
}
tune_custom <- tune_imp(obj, parameters = custom_params, .f = custom_function, cores = 2, rep = 2)

mirai::daemons(0) # Close daemons
```

Then, preferably perform imputation by group with `group_imp` if the
variables can be meaningfully grouped (e.g., by chromosome). `group_imp`
requires the `group` data frame, which contains 3 list columns: 1)
`features`, 2) (optional) `aux`, and 3) (optional) `parameters`. Each
element of the list column `features` is a character vector of variables
to be imputed. Here, we have 2 chromosomes, so the `group` tibble has
two groups (i.e., two rows). PCA imputation can be parallelized with
`{mirai}` similar to `tune_imp`.

``` r
group_df <- tibble::tibble(
  features = lapply(c("chr1", "chr2"), \(x) subset(sim_obj$group_feature, group == x)$feature_id)
)
group_df
# We choose K-NN imputation, k = 5, from the `tune_imp` results.
group_imp(obj, group = group_df, k = 5)

mirai::daemons(2)
pca_results <- group_imp(obj, group = group_df, ncp = 3, cores = 2)
mirai::daemons(0)
```

Full matrix imputation can be performed using `knn_imp` or `pca_imp`.

``` r
full_knn_results <- knn_imp(obj = obj, k = 5)
full_pca_results <- pca_imp(obj = obj, ncp = 5)
```

Sliding window imputation can be performed using `slide_imp`. DNAm
WGBS/EM-seq data should be grouped by chromosome and converted into
either beta or M values before sliding window imputation.

``` r
chr1_beta <- t(sim_mat(m = 10, n = 2000, perc_NA = 0.3, perc_col_NA = 1, nchr = 1)$input)
dim(chr1_beta)
#> [1]   10 2000
chr1_beta[1:5, 1:5]
#>        feat1     feat2     feat3     feat4     feat5
#> s1        NA 0.5207813        NA        NA 0.5614669
#> s2        NA 0.3708712 0.4521467 0.3680747 0.6094274
#> s3        NA 0.3900858        NA        NA 0.3807988
#> s4 0.5407613 0.5204988        NA 0.7050222 0.3691134
#> s5 0.3733871        NA 0.3972191 0.3864667 0.5326897

# Tune the results using the first 50 variables
slide_knn_params <- tibble::tibble(n_feat = c(20, 50), n_overlap = 5, k = 10)
slide_knn_tune <- tune_imp(
  chr1_beta[, 1:50],
  parameters = slide_knn_params,
  .f = "slide_imp",
  cores = 2,
  rep = 2
)
#> Tuning slide_imp
#> Step 1/2: Injecting NA
#> Step 2/2: Tuning
# compute_metrics(slide_knn_tune)

# From the tune results, choose window size of 50, overlap of size 5 between windows,
# K-NN imputation using k = 10. Specify `ncp` for sliding window PCA imputation.
slide_imp(obj = chr1_beta, n_feat = 50, n_overlap = 5, k = 10, cores = 2, .progress = FALSE)
#> ImputedMatrix (KNN)
#> Dimensions: 10 x 2000
#> 
#>        feat1     feat2     feat3     feat4     feat5
#> s1 0.5391722 0.5207813 0.5351055 0.5716223 0.5614669
#> s2 0.5124904 0.3708712 0.4521467 0.3680747 0.6094274
#> s3 0.4347679 0.3900858 0.4207730 0.4928564 0.3807988
#> s4 0.5407613 0.5204988 0.4919336 0.7050222 0.3691134
#> s5 0.3733871 0.4298032 0.3972191 0.3864667 0.5326897
#> 
#> # Showing [1:5, 1:5] of full matrix
```
