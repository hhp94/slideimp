# slideimp

``` r

library(slideimp)
set.seed(1234)
```

## Comparing Methods on One Cross-Validation Split with `tune_imp()`

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

- By default,
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  samples the NA locations itself with
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md).
  To compare methods on one split, generate the locations up front with
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  and pass them to each
  [`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
  call.
- Here,
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  selects 5 rows in each of 200 random columns (features), for 5
  repetitions. To use only certain columns (e.g., clock CpGs), pass the
  `na_col_subset` argument.
  - `na_loc` has 5 elements, one for each repetition. Each element is a
    matrix with one row per masked value, which holds the row and column
    index of that value.

``` r

na_loc <- sample_na_loc(obj, n_cols = 200, n_rows = 5, n_reps = 5)
length(na_loc)
#> [1] 5
na_loc[[1]][1:6, ]
#>      row col
#> [1,]   1 154
#> [2,]   5 154
#> [3,]  12 154
#> [4,]  11 154
#> [5,]   2 154
#> [6,]  17 704
```

- Every method gets the same masked values, so PCA, K-NN and a custom
  method can be compared directly.
- A custom function must take `obj` as its first argument and return an
  object with the same dimensions as `obj`. Its other arguments must
  match the column names of the `parameters` data frame.

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
#> Tuning a custom function
#> Step 1/2: Resolving NA locations
#> Running mode: sequential
#> Step 2/2: Tuning
```

- Pick the method with the lowest cross-validation error:

``` r

mean(compute_metrics(pca_tune, metrics = "rmse")$.estimate)
#> [1] 0.2004882
mean(compute_metrics(knn_tune, metrics = "rmse")$.estimate)
#> [1] 0.1935046
mean(compute_metrics(rnorm_tune, metrics = "rmse")$.estimate)
#> [1] 1.131702
```

## Group-Wise Imputation with `group_imp()`

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  imputes each group of columns (e.g., each chromosome) separately. This
  can significantly increase accuracy for both K-NN and PCA imputation
  ([Pham et al. 2026](https://doi.org/10.1093/bioinformatics/btag318)),
  and it reduces run time:
  - K-NN searches for neighbors within the group only, so it is much
    faster.
  - PCA is usually faster: each group forms a smaller cross-product, has
    less memory pressure, tends to converge in fewer iterations, and
    groups can run in parallel.
- The `group` argument maps `colnames(obj)` to groups. Build it with
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  to set parameters per group or to pad small groups.
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  returns a data frame with these columns:
  - `group`: the group labels.
  - `feature`: a list-column. Each element holds the column names
    imputed together.
  - `aux`: a list-column of auxiliary columns. They help impute the
    group’s features and are not imputed themselves. When a group such
    as `chrM` has fewer than `min_group_size` columns,
    [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
    adds columns drawn at random from other groups to its `aux` until it
    reaches that size.
  - `parameters`: optional group-specific imputation parameters.
- First, simulate data with 2 groups. Then move feature 1 into a new
  group, `group3`, to show how `min_group_size` pads a one-feature group
  through the `aux` list-column.

``` r

sim_obj <- sim_mat(n = 20, p = 50, n_col_groups = 2)

# Matrix to be imputed
obj <- sim_obj$input
obj[1:5, 1:4]
#>          feature1  feature2  feature3  feature4
#> sample1 0.0000000        NA 0.2437860 0.0000000
#> sample2 0.9022484 0.6724700 0.5377050 0.8242256
#> sample3 1.0000000 0.6652907 0.3189869 1.0000000
#> sample4 0.8367542 0.1726599 0.2097971 0.4225234
#> sample5 0.4962462 0.6250307 0.4247322 0.6390391

# Metadata, i.e., which features belong to which group
meta <- sim_obj$col_group
meta[1:5, ]
#>    feature  group
#> 1 feature1 group1
#> 2 feature2 group2
#> 3 feature3 group1
#> 4 feature4 group1
#> 5 feature5 group1

# Put feature 1 in `group3`
meta[1, 2] <- "group3"
meta[1:5, ]
#>    feature  group
#> 1 feature1 group3
#> 2 feature2 group2
#> 3 feature3 group1
#> 4 feature4 group1
#> 5 feature5 group1
```

- Then build the `group` argument with
  [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md).
  The output shows `group3` padded to 10 columns: 1 feature and 9 `aux`
  columns.
- As a demonstration, each group also gets a different number of nearest
  neighbors, `k`.

``` r

set.seed(1234)
group_imp_df <- prep_groups(colnames(obj), group = meta, min_group_size = 10)
group_imp_df$parameters <- list(list(k = 3), list(k = 4), list(k = 5))
group_imp_df
#> # slideimp table: 3 x 4
#>   group          feature             aux parameters
#>  group1 <character [25]> <character [0]> <list [1]>
#>  group2 <character [24]> <character [0]> <list [1]>
#>  group3  <character [1]> <character [9]> <list [1]>
```

- Finally, impute `obj` with the modified `group_imp_df`. Every group
  sets its own `k`, so the `k = 10` passed to
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  is ignored.

``` r

knn_results <- group_imp(obj, group = group_imp_df, cores = 2, k = 10)
#> Imputing 3 groups using K-NN.
#> Running mode: threaded (2 cores)
print(knn_results, p = 4)
#> Method: group_imp (K-NN imputation)
#> Dimensions: 20 x 50
#> 
#>          feature1  feature2  feature3  feature4
#> sample1 0.0000000 0.2177326 0.2437860 0.0000000
#> sample2 0.9022484 0.6724700 0.5377050 0.8242256
#> sample3 1.0000000 0.6652907 0.3189869 1.0000000
#> sample4 0.8367542 0.1726599 0.2097971 0.4225234
#> sample5 0.4962462 0.6250307 0.4247322 0.6390391
#> sample6 0.5124205 0.7016570 0.1918216 0.5038718
#> # Showing 6 of 20 rows and 4 of 50 columns
```

## Multiple Imputation of WGBS/EM-seq Data with `mipca_flank()`

- For whole-genome methylation data, the Bayesian PCA sampler is
  recommended because such data are sparsely observed.
  [`mipca_flank()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
  draws multiple imputations of target CpGs, each fitted on the CpGs
  within `window_size` bp of it.
- Split the data by chromosome first. This simulated beta matrix stands
  in for one chromosome, with 50 to 500 bp between neighboring CpGs.
  `location` holds the position of each column.

``` r

set.seed(1234)
chr1_beta <- sim_mat(n = 30, p = 500)$input
location <- cumsum(sample(50:500, ncol(chr1_beta), replace = TRUE))
chr1_beta[1:4, 1:4]
#>          feature1  feature2   feature3  feature4
#> sample1 0.5409157 0.4077444 0.09981068 0.2434143
#> sample2 0.5224927 0.9663621 0.49824360 0.3923603
#> sample3 0.6355535        NA 1.00000000 1.0000000
#> sample4 0.0000000 0.0000000 0.08115191 0.2106358
```

- [`mipca_flank_plan()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
  finds the window of each target. A target whose window has fewer than
  `min_window_n` columns is not imputed, and a target with no missing
  values is marked `clean`.

``` r

targets <- colnames(chr1_beta)[c(100, 250, 400)]
plan <- mipca_flank_plan(
  chr1_beta,
  location,
  window_size = 2000, # bp on each side of the target
  subset = targets,
  min_window_n = 10
)
plan
#> si_flank_plan: 3 targets on a 30 x 500 matrix, nothing fitted
#>   status: impute 2, clean 1
#>  target       name n_miss start end window_n status
#>     100 feature100      6    94 105       12 impute
#>     250 feature250      0   244 256       13  clean
#>     400 feature400      6   393 406       14 impute
#>   next: mipca_flank(x, rank, seed, transform) imputes.
#>         mipca_flank_window(x, subset) gives one target's window to tune on.
```

- Choose `rank` on one target’s window with
  [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md).
  A window has few columns, so lower `n_cols`, the number of columns
  that get held-out values.

``` r

win <- mipca_flank_window(plan, targets[1])
tuned <- mipca_tune(win, rank = 1:3, transform = "logit_biwhiten", n_cols = 5, seed = 1)
tuned
#> si_tune: 3 ranks on a 30 x 12 matrix, holdout 20 cells over 5 columns
#>  rank    crps ess_p10  ess_ok rhat_max implied_nu implied_lo implied_hi  time
#>     1 0.08725  224.08 0.03571    1.050      2.231      2.143      3.049 0.014
#>     2 0.10157   97.22 0.00000    1.132      2.435      2.281      4.789 0.010
#>     3 0.09657   87.45 0.00000    1.160      3.477      2.698      5.068 0.010
#>   next: mipca_fit(x, obj) fits the rank with the lowest crps.
#>         mipca_fit(x, obj, rank) fits another rank.
```

- This is simulated data, so its convergence figures (`ess_ok`,
  `rhat_max`) say nothing about a real analysis. Check them on your own
  data, as the article below describes.
- Then fit every window with the chosen `rank` and keep the targets’
  draws. Each draw is one completed matrix of the targets.

``` r

fit <- mipca_flank(
  plan,
  rank = tuned$rank[which.min(tuned$crps)],
  seed = 1,
  transform = "logit_biwhiten",
  .progress = FALSE
)
fit
#> si_mi fit: 30 x 3, 12 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
dim(as.matrix(fit, draw = 1))
#> [1] 30  3
```

- See the [multiple imputation
  article](https://hhp94.github.io/slideimp/articles/multiple-imputation-beta.html)
  for `nu`, convergence checks, and pooling the draws.
