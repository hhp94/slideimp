# K-Nearest Neighbor Imputation for Numeric Matrices

Impute missing values in a numeric matrix with K-nearest neighbors
(K-NN), searching the whole matrix for neighbors.

## Usage

``` r
knn_imp(
  obj,
  k,
  colmax = 0.9,
  method = c("euclidean", "manhattan"),
  cores = 1,
  post_imp = TRUE,
  subset = NULL,
  dist_pow = 0,
  .progress = FALSE
)
```

## Arguments

- obj:

  A double matrix. Samples are in rows and features in columns.

- k:

  Integer. Number of nearest neighbors for K-NN imputation.

- colmax:

  Numeric scalar between `0` and `1`. Columns with a proportion of
  missing values greater than `colmax` take no part in the imputation,
  as described under Column mean fill.

- method:

  Character. One of `"euclidean"` or `"manhattan"`. The distance used to
  find neighbors, as defined under Details.

- cores:

  Integer. Number of cores for K-NN imputation. Defaults to `1`.

- post_imp:

  Logical. If `TRUE`, column means replace the missing values that
  remain after imputation, as described under Column mean fill.

- subset:

  Optional character or integer vector. The columns to impute and
  return, by name or by number, or `NULL` for every column. Defaults to
  `NULL`.

- dist_pow:

  Numeric. Power used to penalize more distant neighbors in the weighted
  average, where `0` gives an unweighted average. The weighting rule is
  given under Details in `knn_imp()`.

- .progress:

  Logical. If `TRUE`, show imputation progress.

## Value

A `slideimp_results` object: a numeric matrix with the rows of `obj` and
the columns that `subset` selects, with missing values imputed.

## Details

`knn_imp()` imputes each column from the `k` columns nearest to it,
found by an exact brute-force search. The distance between the target
column and a candidate column is a mean over the rows where both are
observed. Candidates that share different numbers of rows with the
target therefore stay comparable. `method = "euclidean"` averages
squared differences and takes no square root. `method = "manhattan"`
averages absolute differences.

A candidate that shares no observed row with the target is never a
neighbor. A column can therefore have fewer than `k` neighbors, or none,
and K-NN does not impute a column with none.

A column with few observed values can be selected as a neighbor because
it agrees with the target on one or two shared rows. Such a neighbor has
no value at most of the rows to impute. A lower `colmax` keeps these
columns out of the neighbors.

When `dist_pow > 0`, imputed values are distance-weighted averages of
the selected neighbors: a neighbor at distance `d` is weighted
`(d_min / d)^dist_pow`, where `d_min` is the smallest positive distance
among the neighbors chosen for that column. `d_min` cancels out of the
average and only keeps the largest weight at `1`.

Under `"euclidean"`, `d` is a mean squared difference, so weighting by
`1 / d^dist_pow` equals weighting by `1 / d^(2 * dist_pow)` on the
usual, square-rooted euclidean distance. A `dist_pow` tuned for one
`method` therefore does not carry over to the other.

Neighbors are selected once per column. A neighbor that is missing at a
row drops out of the average for that row, which is taken over the
neighbors observed there. K-NN does not impute a row where none of the
neighbors is observed.

With `dist_pow > 0`, a neighbor at distance zero agrees with the target
on every shared row and would take an infinite weight. It takes weight
`1` instead, and the other neighbors take a small positive weight. Exact
matches therefore dominate the average, but the other neighbors stay in
it. When every neighbor is at distance zero, all weights are `1`.

The result holds the `subset` columns in the order they have in `obj`.
Every column within `colmax` can still be a neighbor, in `subset` or
not. A name that is not in `colnames(obj)`, or a `subset` that selects
no column, is an error.

## K-NN performance optimization

- Use `subset` when only specific columns need imputation.

- Use
  [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  or
  [`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
  for very large matrices.

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

## References

Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
Botstein D, Altman RB (2001). Missing value estimation methods for DNA
microarrays. *Bioinformatics*, 17(6), 520-525.
[doi:10.1093/bioinformatics/17.6.520](https://doi.org/10.1093/bioinformatics/17.6.520)

## See also

[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
to impute within feature groups.

## Examples

``` r
set.seed(123)
obj <- sim_mat(20, 20, perc_col_na = 1)$input
sum(is.na(obj))
#> [1] 40

# Select `k` with `tune_imp()`.
result <- knn_imp(obj, k = 10)
result
#> Method: K-NN imputation
#> Dimensions: 20 x 20
#> 
#>           feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.08885928 0.0946424 0.5010982 0.2257198 0.3310293 0.3919172
#> sample2 0.35087647 0.2569208 0.4441198 0.3953534 0.5894579 0.2589739
#> sample3 0.56864697 0.4021824 0.7354948 0.6422099 0.8454413 0.6652821
#> sample4 0.30420093 0.7886995 0.3732225 0.5291319 0.5289648 0.4384578
#> sample5 0.34030847 0.6095144 0.3741498 0.3364915 0.4772101 0.8290134
#> sample6 0.45667473 0.4614949 0.8676809 0.7274044 0.9167450 0.6643710
#> # Showing 6 of 20 rows and 6 of 20 columns

# Impute and return only the first three columns
knn_imp(obj, k = 10, subset = 1:3)
#> Method: K-NN imputation
#> Dimensions: 20 x 3
#> 
#>           feature1  feature2  feature3
#> sample1 0.08885928 0.0946424 0.5010982
#> sample2 0.35087647 0.2569208 0.4441198
#> sample3 0.56864697 0.4021824 0.7354948
#> sample4 0.30420093 0.7886995 0.3732225
#> sample5 0.34030847 0.6095144 0.3741498
#> sample6 0.45667473 0.4614949 0.8676809
#> # Showing 6 of 20 rows and 3 of 3 columns
```
