# PCA Imputation for Numeric Matrices

Impute missing values in a numeric matrix using regularized or
expectation-maximization PCA imputation.

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
  na_check = TRUE,
  clamp = NULL
)
```

## Arguments

- obj:

  A numeric matrix with samples in rows and features in columns.

- ncp:

  Integer. Number of principal components used to predict missing
  entries.

- scale:

  Logical. If `TRUE`, columns are scaled to unit variance.

- method:

  Character. PCA imputation method: either `"regularized"` or `"EM"`.

- coeff.ridge:

  Numeric. Ridge regularization coefficient. Only used when
  `method = "regularized"`. Values less than `1` regularize less, moving
  closer to EM PCA. Values greater than `1` regularize more, moving
  closer to mean imputation.

- row.w:

  Row weights, internally normalized to sum to `1`. Can be:

  - `NULL`: all rows are weighted equally.

  - A numeric vector of positive weights with length `nrow(obj)`.

  - `"n_miss"`: rows with more missing values receive lower weight.

- threshold:

  Numeric. Convergence threshold.

- seed:

  Integer, numeric, or `NULL`. Random seed for reproducibility.

- nb.init:

  Integer. Number of random initializations. The first initialization is
  always mean imputation.

- maxiter:

  Integer. Maximum number of iterations.

- miniter:

  Integer. Minimum number of iterations.

- solver:

  Character. Eigensolver selection. One of `"auto"`, `"exact"`, or
  `"lobpcg"`. `"exact"` uses the exact solver. `"lobpcg"` uses the
  iterative LOBPCG solver with exact fallback. `"auto"` performs a short
  timed probe and chooses LOBPCG only if it is clearly faster than the
  exact solver. When `nb.init > 1`, the auto choice from the first PCA
  initialization is reused for subsequent PCA initializations.

- lobpcg_control:

  A list of LOBPCG eigensolver control options, usually created by
  [`lobpcg_control()`](https://hhp94.github.io/slideimp/reference/lobpcg_control.md).
  A plain named list is also accepted. Ignored when `solver = "exact"`.

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

- clamp:

  Optional numeric vector of length 2 giving lower and upper bounds for
  PCA-imputed values. Use `NULL` for no clamping. Use `c(0, 1)` for DNA
  methylation beta values. Use `c(lb, Inf)` for only lower bound
  clamping, or `c(-Inf, ub)` for only upper bound clamping. Clamping is
  applied only to values imputed by the PCA step, not to observed
  values.

## Value

A numeric matrix of the same dimensions as `obj`, with missing values
imputed. The returned object has class `slideimp_results`.

## Details

This algorithm is based on
[`missMDA::imputePCA()`](https://rdrr.io/pkg/missMDA/man/imputePCA.html)
and is optimized for tall or wide numeric matrices.

## Performance tips

`pca_imp()` relies heavily on linear algebra. On Windows, the default
BLAS shipped with R may be slow for large matrices. Advanced users can
replace it with
[OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).

PCA imputation speed depends on the eigensolver selected by `solver` and
the convergence threshold `threshold`. The exact solver is selected with
`solver = "exact"`. The iterative LOBPCG solver is selected with
`solver = "lobpcg"`. The default, `solver = "auto"`, performs a short
timed probe and chooses LOBPCG only when it is clearly faster.

For large or approximately low-rank genomic matrices, it can be useful
to benchmark `solver = "exact"` against `solver = "lobpcg"` on a
representative subset, such as chromosome 22, before tuning
accuracy-related parameters. For
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md),
this may include `window_size` and `overlap_size`.

The default `threshold = 1e-6` is conservative. In many genomic
datasets, `threshold = 1e-5` can be faster while giving very similar
imputed values. Check this on a representative subset before using the
relaxed threshold in a full analysis.

See the pkgdown article [Speeding up PCA
imputation](https://hhp94.github.io/slideimp/articles/speeding-up-pca-imputation.html)
for a full workflow.

## References

Josse J, Husson F (2013). Handling missing values in exploratory
multivariate data analysis methods. *Journal de la SFdS*, 153(2), 79-99.

Josse J, Husson F (2016). missMDA: A Package for Handling Missing Values
in Multivariate Data Analysis. *Journal of Statistical Software*, 70(1),
1-31. [doi:10.18637/jss.v070.i01](https://doi.org/10.18637/jss.v070.i01)

The PCA imputation algorithm is based on the original `missMDA`
implementation by Francois Husson and Julie Josse.

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

# Randomly initialize missing values 5 times. The first initialization is
# mean imputation.
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
```
