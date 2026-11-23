# Grouped K-NN or PCA Imputation

Perform K-NN or PCA imputation independently within feature groups.

## Usage

``` r
group_imp(
  obj,
  group,
  subset = NULL,
  allow_unmapped = FALSE,
  k = NULL,
  ncp = NULL,
  method = NULL,
  cores = 1,
  .progress = TRUE,
  min_group_size = NULL,
  colmax = NULL,
  post_imp = NULL,
  dist_pow = NULL,
  scale = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  row.w = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL,
  solver = NULL,
  lobpcg_control = NULL,
  clamp = NULL,
  pin_blas = TRUE,
  on_infeasible = c("error", "skip", "mean")
)
```

## Arguments

- obj:

  A double matrix. Samples are in rows and features in columns.

- group:

  A data frame or a character scalar. The feature groups, in one of
  these forms:

  - A character scalar naming a supported Illumina platform, such as
    `"EPICv2"`. This needs the `slideimp.extra` package.

  - A long-format data frame with columns `group` and `feature`.

  - A list-column data frame with a `feature` list-column. Optional
    list-columns are `aux`, for auxiliary feature names, and
    `parameters`, for group-specific parameter lists.

- subset:

  Optional character vector. The features to impute and return, by name,
  or `NULL` to impute every grouped feature and return every column. A
  group's other features become its auxiliary columns, and a group left
  with no feature is dropped with a message.

- allow_unmapped:

  Logical. If `FALSE`, every column in `colnames(obj)` must appear in
  `group`. If `TRUE`, columns with no group assignment are left
  untouched and are not used as auxiliary columns.

- k:

  Integer or `NULL`. Number of nearest neighbors for K-NN imputation. If
  `NULL`, `k` may be supplied through `group$parameters`.

- ncp:

  Integer or `NULL`. Number of components for PCA imputation. If `NULL`,
  `ncp` may be supplied through `group$parameters`.

- method:

  Character or `NULL`. One of `"euclidean"` or `"manhattan"` for K-NN
  imputation, or `"regularized"` or `"EM"` for PCA imputation. `NULL`
  uses `"euclidean"` for K-NN and `"regularized"` for PCA.

- cores:

  Integer. Number of cores for K-NN imputation only, as described under
  Parallelization.

- .progress:

  Logical. If `TRUE`, show progress.

- min_group_size:

  Integer or `NULL`. Minimum total number of columns per group, counting
  both features and auxiliary columns. Groups smaller than this are
  padded with columns sampled at random from the other grouped columns
  of `obj`.

- colmax:

  Numeric scalar between `0` and `1`. Columns with a proportion of
  missing values greater than `colmax` take no part in the imputation,
  as described under Column mean fill.

- post_imp:

  Logical. If `TRUE`, column means replace the missing values that
  remain after imputation, as described under Column mean fill.

- dist_pow:

  Numeric. Power used to penalize more distant neighbors in the weighted
  average, where `0` gives an unweighted average. The weighting rule is
  given under Details in
  [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md).

- scale:

  Logical. If `TRUE`, columns are scaled to unit variance.

- coeff.ridge:

  Numeric. The ridge regularization for `method = "regularized"`,
  ignored otherwise. Values below `1` move toward EM PCA, and values
  above `1` move toward mean imputation.

- threshold:

  Numeric. Convergence threshold.

- row.w:

  Numeric vector, `"n_miss"`, or `NULL`. The row weights, normalized to
  sum to `1`: positive values of length `nrow(obj)`, `"n_miss"` to
  down-weight rows with more missing values, or `NULL` for equal
  weights.

- seed:

  Integer, numeric, or `NULL`. Random seed used when padding small
  groups and, for PCA imputation, passed to
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).

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

- pin_blas:

  Logical. If `TRUE`, pin BLAS threads to 1 on the `mirai` daemons in a
  parallel run. Defaults to `TRUE`.

- on_infeasible:

  Character. One of `"error"`, `"skip"`, or `"mean"`. The action for a
  group that is infeasible for imputation, for example because `k` or
  `ncp` is too large for the columns that `colmax` leaves. Defaults to
  `"error"`.

## Value

A `slideimp_results` object: a numeric matrix with the rows of `obj` and
the columns that `subset` selects, with the missing values of the
grouped features imputed.

## Details

Imputing each group on its own can reduce the runtime for large
matrices.

Specify `k` and related arguments to use K-NN imputation, or `ncp` and
related arguments to use PCA imputation. If both `k` and `ncp` are
`NULL`, `group$parameters` must supply either `k` or `ncp` for every
group. Passing an argument of the other method, such as `scale` with
`k`, is an error.

Group-specific parameters in `group$parameters` take priority over
global arguments, and global arguments fill in the parameters a group
does not set. All groups must use the same imputation method. A global
argument shared with
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md) or
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
that is `NULL` takes that function's default.

Per-group `k` is capped at the group's column count minus one, counting
features and auxiliary columns. Per-group `ncp` is capped at the smaller
of that cap and `nrow(obj) - 2`. `group_imp()` warns when it caps a
value.

A group is infeasible when it is too small for any valid `k` or `ncp`,
or when
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md) or
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
cannot impute it, for example because `colmax` leaves too few columns.
With `on_infeasible = "skip"`, an infeasible group keeps its missing
values. With `"mean"`, column means fill them, and with `"error"`,
`group_imp()` stops with an error. A group whose features have no
missing values is not imputed and is never infeasible. An infeasible
group adds nothing to the `mean_filled_cols` and `mean_filled_cells`
attributes of the result.

With `subset`, the result holds only the `subset` columns, in the order
they have in `obj`. The other features still take part in the
imputation, as auxiliary columns of their groups.

The `mirai` daemons run R's `"L'Ecuyer-CMRG"` generator. So with
`nb.init > 1`, a parallel PCA run draws different random initializations
from a sequential run under R's default generator. The two agree when
the session also runs `RNGkind("L'Ecuyer-CMRG")`.

## Note

A platform name for `group`, such as `"EPICv2"` or `"EPICv2_deduped"`,
needs the optional `slideimp.extra` package. Its `slideimp_arrays`
object lists the supported platforms.

## Parallelization

- K-NN: use the `cores` argument. If `mirai` daemons are active, `cores`
  is set to `1` to avoid nested parallelism.

- PCA: use
  [`mirai::daemons()`](https://mirai.r-lib.org/reference/daemons.html)
  instead of `cores`.

- With `mirai` daemons, each daemon runs BLAS on one thread unless
  `pin_blas = FALSE`.

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
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
windows, where it is faster. It is also more reliable: it always
converges, to eigenvectors exact to machine precision, while LOBPCG
stops at the tolerance set in
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

## See also

- [`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
  to inspect and edit the groups and their parameters.

- [`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
  and
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  to impute a single matrix.

## Examples

``` r
set.seed(1234)
to_test <- sim_mat(10, 20, perc_total_na = 0.05, perc_col_na = 1)
obj <- to_test$input
group <- to_test$col_group
head(group)
#>    feature  group
#> 1 feature1 group2
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group2
#> 5 feature5 group1
#> 6 feature6 group1

# Simple grouped K-NN imputation
results <- group_imp(obj, group = group, k = 2, .progress = FALSE)
#> Imputing 2 groups using K-NN.
#> Running mode: sequential
results
#> Method: group_imp (K-NN imputation)
#> Dimensions: 10 x 20
#> 
#>          feature1  feature2  feature3  feature4   feature5  feature6
#> sample1 0.1568098 0.3213953 0.2768746 1.0000000 0.07839046 0.4219375
#> sample2 0.4098416 0.8658918 0.8221066 0.6396885 0.68926345 1.0000000
#> sample3 0.6801685 1.0000000 1.0000000 0.9953450 0.75246030 0.6958550
#> sample4 0.0000000 0.0000000 0.0000000 0.0000000 0.00000000 0.0000000
#> sample5 0.9639671 0.6409137 0.5111760 0.7184678 0.81815288 0.5883255
#> sample6 0.7031741 0.3731062 0.6782811 0.7542872 0.99910554 0.9070005
#> # Showing 6 of 10 rows and 6 of 20 columns

# Impute only a subset of features
subset_features <- sample(to_test$col_group$feature, size = 10)
knn_subset <- group_imp(
  obj,
  group = group,
  subset = subset_features,
  k = 2,
  .progress = FALSE
)
#> Imputing 2 groups using K-NN.
#> Running mode: sequential

# Use prep_groups() to inspect and edit per-group parameters
prepped <- prep_groups(colnames(obj), group)
prepped$parameters <- lapply(seq_len(nrow(prepped)), function(i) list(k = 2))
prepped$parameters[[2]]$k <- 4
knn_grouped <- group_imp(obj, group = prepped, .progress = FALSE)
#> Imputing 2 groups using K-NN.
#> Running mode: sequential

if (FALSE) { # interactive() && requireNamespace("mirai", quietly = TRUE)
# PCA imputation with mirai parallelism
mirai::daemons(2)
pca_grouped <- group_imp(obj, group = group, ncp = 2)
mirai::daemons(0)
pca_grouped
}
```
