# PCA Imputation for Numeric Matrices

Impute missing values in a numeric matrix using regularized or
expectation-maximization (EM) PCA imputation. Supports warm-start LOBPCG
with both the previous eigenblock and search direction.

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
  clamp = NULL,
  .progress = FALSE
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

  Numeric. Ridge regularization, used only when
  `method = "regularized"`. Values `< 1` move toward EM PCA; values
  `> 1` move toward mean imputation.

- row.w:

  Row weights, normalized to sum to `1`. `NULL` (equal weights), a
  positive numeric vector of length `nrow(obj)`, or `"n_miss"`
  (down-weight rows with more missing values).

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

  Character. Eigensolver: `"auto"` (default), `"exact"`, or `"lobpcg"`.
  `"auto"` runs a short timed probe and picks `"lobpcg"` only when
  clearly faster. Consecutive EM calls warm-start LOBPCG with both the
  previous eigenblock and search direction. When `nb.init > 1`, the auto
  choice from the first init is reused. See Performance tips.

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

  Optional numeric vector `c(lower, upper)` bounding PCA-imputed values
  (use `-Inf`/`Inf` for one-sided, `NULL` for none). E.g., `c(0, 1)` for
  DNAm beta values. Observed values are not clamped.

- .progress:

  Logical. If `TRUE`, show imputation progress.

## Value

A numeric matrix of the same dimensions as `obj`, with missing values
imputed. The returned object has class `slideimp_results`.

## Details

This function reimplements the PCA imputation method from the `missMDA`
package by Francois Husson and Julie Josse, based on Josse and Husson
(2016).

## PCA Performance tips

Speed comes from three levers: `solver` (through LOBPCG with
warm-start), `threshold`, and `scale`. Tune these first, then accuracy
parameters (`ncp`, `coeff.ridge`) on a representative subset.

**Exact vs. LOBPCG with warm-start.** Whether `"lobpcg"` beats `"exact"`
depends on size and low-rankness: `"lobpcg"` is preferred for large,
approximately low-rank matrices with small `ncp`, and `"exact"` for
small matrices (including
[`slide_imp()`](https://hhp94.github.io/slideimp/reference/slide_imp.md)
windows), where it is faster and more robust. Separately, the warm-start
makes each successive solve cheap: `pca_imp()` warm-starts LOBPCG with
the previous eigenblock and search direction, so once imputed values
stabilize, later solves converge in a few iterations. The payoff
therefore grows with the number of EM iterations, independent of
low-rankness. `solver = "auto"` (default) probes both and is a safe
start.

**Threshold.** The default `1e-6` is conservative; `1e-5` is often
faster with very similar values.

**Scale.** For columns on a common scale (e.g., DNAm beta values in
`[0, 1]`), `scale = FALSE` can be faster and more accurate.

**Parallel and BLAS.** In parallel via
[`tune_imp()`](https://hhp94.github.io/slideimp/reference/tune_imp.md)
or
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
with a multithreaded BLAS, set `pin_blas = TRUE` to avoid thread
oversubscription. On Windows, the stock BLAS can be slow. Advanced users
can swap in
[OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).

See [Speeding up PCA
imputation](https://hhp94.github.io/slideimp/articles/speeding-up-pca-imputation.html)
for the full workflow.

## References

Josse J, Husson F (2013). Handling missing values in exploratory
multivariate data analysis methods. *Journal de la SFdS*, 153(2), 79-99.

Josse J, Husson F (2016). missMDA: A Package for Handling Missing Values
in Multivariate Data Analysis. *Journal of Statistical Software*, 70(1),
1-31. [doi:10.18637/jss.v070.i01](https://doi.org/10.18637/jss.v070.i01)

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
# mean imputation. Select `ncp` with `tune_imp()`.
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
