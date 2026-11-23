# PCA Imputation for Numeric Matrices

Impute missing values in a numeric matrix with regularized or
expectation-maximization (EM) PCA imputation. The eigensolver can be
LOBPCG, warm-started from the previous eigenblock and search direction.

## Usage

``` r
pca_imp(
  obj,
  ncp = 2,
  scale = TRUE,
  method = c("regularized", "EM"),
  coeff.ridge = 1,
  row.w = NULL,
  threshold = 1e-06,
  seed = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  solver = c("auto", "exact", "lobpcg"),
  lobpcg_control = NULL,
  colmax = 0.9,
  post_imp = TRUE,
  clamp = NULL,
  subset = NULL,
  .progress = FALSE
)
```

## Arguments

- obj:

  A double matrix. Samples are in rows and features in columns.

- ncp:

  Integer. Number of principal components used to impute missing values.

- scale:

  Logical. If `TRUE`, columns are scaled to unit variance.

- method:

  Character. One of `"regularized"` or `"EM"`. The PCA imputation
  method.

- coeff.ridge:

  Numeric. The ridge regularization for `method = "regularized"`,
  ignored otherwise. Values below `1` move toward EM PCA, and values
  above `1` move toward mean imputation.

- row.w:

  Numeric vector, `"n_miss"`, or `NULL`. The row weights, normalized to
  sum to `1`: positive values of length `nrow(obj)`, `"n_miss"` to
  down-weight rows with more missing values, or `NULL` for equal
  weights.

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

- colmax:

  Numeric scalar between `0` and `1`. Columns with a proportion of
  missing values greater than `colmax` take no part in the imputation,
  as described under Column mean fill.

- post_imp:

  Logical. If `TRUE`, column means replace the missing values that
  remain after imputation, as described under Column mean fill.

- clamp:

  Numeric vector of length 2 or `NULL`. The bounds `c(lower, upper)` for
  PCA-imputed values, such as `c(0, 1)` for DNA methylation beta values,
  or `NULL` for none. Use `-Inf` or `Inf` for a one-sided bound.
  Observed values are not clamped.

- subset:

  Optional character or integer vector. The columns to impute and
  return, by name or by number, or `NULL` for every column. Defaults to
  `NULL`.

- .progress:

  Logical. If `TRUE`, show imputation progress.

## Value

A `slideimp_results` object: a numeric matrix with the rows of `obj` and
the columns that `subset` selects, with missing values imputed.

## Details

`pca_imp()` reimplements the PCA imputation method from the `missMDA`
package by Francois Husson and Julie Josse, based on Josse and Husson
(2016).

A column whose observed values do not vary takes no part in PCA
imputation. Column mean fill treats it like a column over `colmax`, in
`post_imp` and in the `mean_filled_cols` attribute.

With `subset`, the result holds only the selected columns, in the order
they have in `obj`. The other columns take part in the imputation as
they do without `subset`. So `subset` makes the result smaller, not the
imputation faster. A name that is not in `colnames(obj)`, or a `subset`
that selects no column, is an error.

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
later LOBPCG solve cheap: `pca_imp()` starts LOBPCG from the previous
eigenblock and search direction, so once the imputed values stabilize,
each solve takes a few iterations. The gain grows with the number of EM
iterations, whatever the rank. `solver = "auto"`, the default, chooses
between the two and is a safe start.

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

Josse J, Husson F (2013). Handling missing values in exploratory
multivariate data analysis methods. *Journal de la SFdS*, 153(2), 79-99.

Josse J, Husson F (2016). missMDA: A Package for Handling Missing Values
in Multivariate Data Analysis. *Journal of Statistical Software*, 70(1),
1-31. [doi:10.18637/jss.v070.i01](https://doi.org/10.18637/jss.v070.i01)

## See also

- [`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
  to impute within feature groups.

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) for
  multiple imputation with PCA.

## Examples

``` r
set.seed(123)
obj <- sim_mat(10, 10)$input
sum(is.na(obj))
#> [1] 10
obj[1:4, 1:4]
#>          feature1   feature2  feature3  feature4
#> sample1 0.5784798 0.06296727 0.3155309 0.1199980
#> sample2 0.4991812 0.44077231 0.2120510 0.3257524
#> sample3 0.7709271 0.75477764 1.0000000 0.5099311
#> sample4 0.5068375 0.37347042 0.6018860 1.0000000

# Five initializations: the first is mean imputation, the other four are
# random. Select `ncp` with `tune_imp()`.
pca_imp(obj, ncp = 2, nb.init = 5, seed = 123)
#> Method: PCA imputation
#> Dimensions: 10 x 10
#> 
#>          feature1   feature2  feature3  feature4  feature5   feature6
#> sample1 0.5784798 0.06296727 0.3155309 0.1199980 0.1807391 0.31919912
#> sample2 0.4991812 0.44077231 0.2120510 0.3257524 0.1919521 0.14843947
#> sample3 0.7709271 0.75477764 1.0000000 0.5099311 0.6027900 0.75450992
#> sample4 0.5068375 0.37347042 0.6018860 1.0000000 0.5850265 0.08171420
#> sample5 0.4165827 0.42553421 0.6024750 0.7728095 0.2295133 0.08343629
#> sample6 1.0000000 0.93942257 0.9867413 0.5851340 1.0000000 1.00000000
#> # Showing 6 of 10 rows and 6 of 10 columns

# Impute every column, but return only the first three
pca_imp(obj, ncp = 2, subset = 1:3)
#> Method: PCA imputation
#> Dimensions: 10 x 3
#> 
#>          feature1   feature2  feature3
#> sample1 0.5784798 0.06296727 0.3155309
#> sample2 0.4991812 0.44077231 0.2120510
#> sample3 0.7709271 0.75477764 1.0000000
#> sample4 0.5068375 0.37347042 0.6018860
#> sample5 0.4165827 0.42553421 0.6024750
#> sample6 1.0000000 0.93942190 0.9867413
#> # Showing 6 of 10 rows and 3 of 3 columns
```
