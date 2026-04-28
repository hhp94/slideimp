# slideimp

``` r
library(slideimp)
set.seed(1234)
```

## Comparing Methods Using a Shared Cross-Validation Split with `tune_imp()`

- Simulate some data

``` r
# 20 rows, 1000 columns, all columns have at least some NA
sim_obj <- sim_mat(n = 20, p = 1000, perc_col_na = 1)
obj <- sim_obj$input
obj[1:4, 1:4]
#>          feature1  feature2  feature3  feature4
#> sample1 0.3486482 0.7385414 0.4077444 0.1607935
#> sample2 0.5338935 0.4724364 0.9663621 0.3722729
#> sample3 0.7185848 0.7351035 0.6724479 0.3162537
#> sample4 0.1734418        NA 0.0000000 0.0000000
```

- Instead of random NA sampling inside
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md),
  which internally calls
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md),
  we can generate the NA locations up front and pass them to
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md).
- Here, we randomly select 5 rows from each of 200 random columns
  (features) for 5 repetitions using
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md).
  To specify just certain columns (i.e., clock CpGs), provide the
  `na_col_subset` argument.
  - `na_loc` has 5 elements, one for each repeat. Each element is a
    matrix whose rows store the row and column index of a missing value.

``` r
na_loc <- sample_na_loc(obj, n_cols = 200, n_rows = 5, n_reps = 5)
length(na_loc)
#> [1] 5
na_loc[[1]][1:6, ]
#>      row col
#> [1,]  18 661
#> [2,]  17 661
#> [3,]  20 661
#> [4,]  15 661
#> [5,]   6 661
#> [6,]  10 293
```

- Then we can compare 1) PCA, 2) K-NN imputation, and 3) a custom
  method, since the cross-validation missing values are the same for all
  methods.
- **Note**: The custom function requires the first argument to be `obj`,
  must return an object with the same dimensions as `obj`, and all
  subsequent arguments must match the column names of the `parameters`
  data frame.

``` r
# This custom function imputes missing values with random normal values and takes
# `mean` and `sd` as params
rnorm_imp <- function(obj, mean, sd) {
  na <- is.na(obj)
  obj[na] <- rnorm(sum(na), mean = mean, sd = sd) # <- impute values with rnorm
  return(obj) # <- return an imputed object with the same dim as obj
}

pca_tune <- tune_imp(
  obj,
  .f = "pca_imp",
  na_loc = na_loc,
  parameters = data.frame(ncp = 10)
)
#> Tuning `pca_imp()`
#> Step 1/2: Resolving NA locations
#> Running mode: sequential
#> Step 2/2: Tuning

knn_tune <- tune_imp(
  obj,
  .f = "knn_imp",
  na_loc = na_loc,
  parameters = data.frame(k = 10)
)
#> Tuning `knn_imp()`
#> Step 1/2: Resolving NA locations
#> Running mode: sequential
#> Step 2/2: Tuning

rnorm_tune <- tune_imp(
  obj,
  .f = rnorm_imp,
  na_loc = na_loc,
  parameters = data.frame(mean = 0, sd = 1) # must match with arguments of `rnorm_imp`
)
#> Tuning custom function
#> Step 1/2: Resolving NA locations
#> Running mode: sequential
#> Step 2/2: Tuning
```

- Pick the method with the lowest cross-validation error:

``` r
mean(compute_metrics(pca_tune, metrics = "rmse")$.estimate)
#> [1] 0.2048011
mean(compute_metrics(knn_tune, metrics = "rmse")$.estimate)
#> [1] 0.1962497
mean(compute_metrics(rnorm_tune, metrics = "rmse")$.estimate)
#> [1] 1.110258
```

## Group-Wise Imputation with Small-Group Padding and Group-Wise Parameters with `group_imp()`

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  allows imputation to be performed separately within defined groups
  (e.g., by chromosome), which significantly reduces run time and can
  increase accuracy for both K-NN and PCA imputation.
- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  requires the `group` argument, which maps `colnames(obj)` to groups.
  This can be created up front with
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  for advanced features such as group-wise parameters and padding of
  small groups with random features from other groups.
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  returns a list-column data frame with:
  - `features`: **required** - a list-column where each element is a
    character vector of variable names to be imputed together.
  - `aux`: **optional** - auxiliary variables to include in each group.
    These are only used to augment the imputation quality of features
    and are not imputed themselves. If one group is too small,
    e.g. `chrM`, `aux` is used to pad the group by randomly drawing
    features from other groups to meet `min_group_size`.
  - `parameters`: **optional** - group-specific imputation parameters.
- First we simulate data from 2 groups. We then create `group3` with
  only 1 feature to show how `min_group_size` pads it using the `aux`
  list column.

``` r
sim_obj <- sim_mat(n = 20, p = 50, n_col_groups = 2)

# Matrix to be imputed
obj <- sim_obj$input
obj[1:5, 1:4]
#>           feature1  feature2  feature3  feature4
#> sample1 0.03031619 0.1352710 0.4653579 0.1143298
#> sample2 0.50962966 0.2313022 0.9737881 0.5198844
#> sample3 0.69068459 0.5020337 0.7863144 0.9997625
#> sample4 0.40710170 0.4094603 0.3228294 0.2529592
#> sample5 0.73113251 0.4108405 0.6385945        NA

# Metadata, i.e., which features belong to which group
meta <- sim_obj$col_group
meta[1:5, ]
#>    feature  group
#> 1 feature1 group1
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group2
#> 5 feature5 group2

# We put feature 1 in `group3`
meta[1, 2] <- "group3"
meta[1:5, ]
#>    feature  group
#> 1 feature1 group3
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group2
#> 5 feature5 group2
```

- Then we generate the `group` parameter for
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  up front. We can see that `group3` has been padded to have 10 columns.
- For each group, we also assign a different number of nearest neighbors
  as a demonstration.

``` r
set.seed(1234)
group_imp_df <- prep_groups(colnames(obj), group = meta, min_group_size = 10)
group_imp_df$parameters <- list(list(k = 3), list(k = 4), list(k = 5))
group_imp_df
#> # slideimp table: 3 x 4
#>   group          feature             aux parameters
#>  group1 <character [26]> <character [0]> <list [1]>
#>  group2 <character [23]> <character [0]> <list [1]>
#>  group3  <character [1]> <character [9]> <list [1]>
```

- Finally, we impute `obj` using the modified `group_imp_df`. The
  `k = 10` passed to
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  is ignored since all groups have group-wise `k` specified.

``` r
knn_results <- group_imp(obj, group = group_imp_df, cores = 4, k = 10)
#> Imputing 3 groups using KNN.
#> Running mode: threaded (4 cores)
print(knn_results, p = 4)
#> Method: group_imp (KNN imputation)
#> Dimensions: 20 x 50
#> 
#>           feature1  feature2  feature3  feature4
#> sample1 0.03031619 0.1352710 0.4653579 0.1143298
#> sample2 0.50962966 0.2313022 0.9737881 0.5198844
#> sample3 0.69068459 0.5020337 0.7863144 0.9997625
#> sample4 0.40710170 0.4094603 0.3228294 0.2529592
#> sample5 0.73113251 0.4108405 0.6385945 0.4472338
#> sample6 0.48993078 0.4950764 0.6553896 0.6737998
#> # Showing 6 of 20 rows and 4 of 50 columns
```

## Sliding Window Imputation for WGBS/EM-seq Data with `slide_imp()`

### Select `window_size`, `overlap_size`, and PCA/K-NN Parameters

- We simulate the output of the `{methylKit}` package.
  - **Note:** WGBS/EM-seq data should be grouped by chromosome before
    performing sliding window imputation.
  - Here, we simulate 1000 sites.
- Importantly, we simulate the location vector of each feature (genomic
  position) with varying spacing of 50 to 500 bp apart from each other.
- The `locations` vector contains the genomic position of each feature
  (column). It is used to determine which columns are grouped together
  given a window size.

``` r
set.seed(1234)
sample_names <- paste0("S", 1:10)
n_sites <- 1000

# Simulate positions with 50–500 bp between each site
distances_between <- sample(50:500, size = n_sites, replace = TRUE)
locations <- cumsum(distances_between) # <- important, location vector

methyl <- data.frame(
  chr = "chr1",
  start = locations,
  end = locations,
  strand = "+"
)

for (i in seq_along(sample_names)) {
  methyl[[paste0("numCs", i)]] <- sample.int(100, size = n_sites, replace = TRUE)
  methyl[[paste0("numTs", i)]] <- sample.int(100, size = n_sites, replace = TRUE)
  methyl[[paste0("coverage", i)]] <- methyl[[paste0("numCs", i)]] + methyl[[paste0("numTs", i)]]
}

methyl[1:5, 1:10]
#>    chr start  end strand numCs1 numTs1 coverage1 numCs2 numTs2 coverage2
#> 1 chr1   333  333      +     83     68       151     84     36       120
#> 2 chr1   718  718      +     82      2        84     54     99       153
#> 3 chr1  1173 1173      +     48      2        50     44     16        60
#> 4 chr1  1323 1323      +     70     48       118     48     56       104
#> 5 chr1  1483 1483      +     88     15       103     83     19       102
```

- First, we convert the data into a beta matrix with features in the
  columns.

``` r
numCs_matrix <- as.matrix(methyl[, paste0("numCs", seq_along(sample_names))])
cov_matrix <- as.matrix(methyl[, paste0("coverage", seq_along(sample_names))])
beta_matrix <- numCs_matrix / cov_matrix

colnames(beta_matrix) <- sample_names
rownames(beta_matrix) <- methyl$start

beta_matrix <- t(beta_matrix)
# Set 10% of the data to missing
set.seed(1234)
beta_matrix[sample.int(length(beta_matrix), floor(length(beta_matrix) * 0.1))] <- NA
beta_matrix[1:4, 1:4]
#>           333       718       1173      1323
#> S1 0.54966887 0.9761905 0.96000000 0.5932203
#> S2 0.70000000 0.3529412 0.73333333 0.4615385
#> S3 0.06521739 0.7053571         NA 0.6507937
#> S4 0.69444444 0.2846154 0.04347826 0.7636364
```

- Then, in a real dataset, we would tune hyperparameters using `chr22`.
  Here, as a demonstration, we use the whole data since the size is
  small.
- We use 2 repetitions of cross-validation (increase to 10–30 in real
  analyses). We are selecting between:
  - `ncp` (number of principal components) of 2 or 4, indicating that we
    are performing sliding PCA imputation. Pass `k` for sliding K-NN
    imputation.
  - `window_size` of 5,000 or 10,000 bp.
  - `overlap_size` fixed at 1,000 bp (does not affect results much in
    real analyses).

``` r
params <- expand.grid(ncp = c(2, 4), window_size = c(5000, 10000))
params$overlap_size <- 1000
params$min_window_n <- 20 # windows with less than 20 columns are dropped

# Increase n_reps from 2 in actual analyses and use another chromosome (i.e., chr22)
tune_slide_pca <- tune_imp(
  obj = beta_matrix,
  parameters = params,
  .f = "slide_imp",
  n_reps = 2,
  location = locations
)
#> Tuning `slide_imp()`
#> Step 1/2: Resolving NA locations
#> ℹ Using default `num_na` = 500 (~5% of cells).
#>   Increase for more reliability or decrease if missing is dense.
#> Running mode: sequential
#> Step 2/2: Tuning

metrics <- compute_metrics(tune_slide_pca)

aggregate(.estimate ~ .metric + ncp + window_size, data = metrics, FUN = mean)
#>   .metric ncp window_size .estimate
#> 1     mae   2        5000 0.2257168
#> 2    rmse   2        5000 0.2838435
#> 3     mae   4        5000 0.2395280
#> 4    rmse   4        5000 0.2984013
#> 5     mae   2       10000 0.2123151
#> 6    rmse   2       10000 0.2642956
#> 7     mae   4       10000 0.2223801
#> 8    rmse   4       10000 0.2756162
```

- Finally, we can use
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  to impute the full `beta_matrix`. Use the best parameter combination
  from the cross-validation metrics.
- First, we use `dry_run = TRUE` to examine the columns to be imputed.
  - `start` and `end` are window location vectors.
  - `window_n` is the number of features included in the window.

``` r
slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 5000,
  overlap_size = 1000,
  ncp = 2,
  min_window_n = 20,
  dry_run = TRUE # <- dry_run to inspect the windows
)
#> Dropping 43 window(s) with fewer than 20 columns.
#> # slideimp table: 24 x 4
#>  start end window_n  subset_local
#>     15  34       20 <double [20]>
#>    104 124       21 <double [21]>
#>    188 209       22 <double [22]>
#>    205 225       21 <double [21]>
#>    222 241       20 <double [20]>
#>    239 258       20 <double [20]>
#>    254 273       20 <double [20]>
#>    369 388       20 <double [20]>
#>    385 404       20 <double [20]>
#>    402 424       23 <double [23]>
#> # ... with 14 more rows
```

- Turn off `dry_run` to impute the data

``` r
slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 5000,
  overlap_size = 1000,
  ncp = 2,
  min_window_n = 20,
  dry_run = FALSE,
  .progress = FALSE
)
#> Method: slide_imp (PCA imputation)
#> Dimensions: 10 x 1000
#> Note: requested columns still contain NA values
#> 
#>           333       718       1173      1323      1483       1925
#> S1 0.54966887 0.9761905 0.96000000 0.5932203        NA 0.08695652
#> S2 0.70000000 0.3529412 0.73333333 0.4615385 0.8137255 0.45454545
#> S3 0.06521739 0.7053571         NA 0.6507937 0.4393064 0.37735849
#> S4 0.69444444 0.2846154 0.04347826 0.7636364 0.9595960 0.54621849
#> S5 0.83505155 0.5777778 0.47517730        NA 0.7368421 0.21666667
#> S6 0.26612903 0.3451327 0.53072626 0.5000000 0.6465517 0.37288136
#> # Showing 6 of 10 rows and 6 of 1000 columns
```

### Subset and Flanking Mode

- Use this mode when you only need to impute specific target features
  (e.g., clock CpGs) rather than the entire dataset.
  - Pass the desired feature names to the `subset` argument. Only
    windows containing these features will be imputed.
  - Set `flank = TRUE` to build windows *centered* on each feature in
    the `subset`. Each window will extend `window_size` bp on either
    side of the target feature (flanking mode).
  - In this mode, the `overlap_size` argument is ignored.
- In this example, we only want to impute the features `"1323"` and
  `"33810"` by creating 5,000 bp flanking windows around each feature:

``` r
slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 5000,
  ncp = 2,
  min_window_n = 20,
  subset = c("1323", "33810"),
  flank = TRUE,
  dry_run = TRUE
)
#> # slideimp table: 2 x 5
#>  start end window_n target  subset_local
#>      1  22       22      4 <integer [1]>
#>    101 138       38    122 <integer [1]>
```
