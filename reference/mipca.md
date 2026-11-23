# Multiple Imputation by Bayesian PCA

Draw multiple imputations for the missing values of a numeric matrix
under a low-rank model using the Audigier, Husson and Josse (BayesMIPCA)
data-augmentation sampler.

## Usage

``` r
mipca(
  obj,
  rank,
  transform,
  subset = NULL,
  keep_state = FALSE,
  prop_diag = 0.5,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  refresh = 0L,
  spectra = spectra_opts(),
  init = 1,
  seed = NULL,
  abb = NULL
)
```

## Arguments

- obj:

  A double matrix with samples in rows and features in columns. Row and
  column names are required and must be unique. `NA` and `NaN` values
  are imputed.

- rank:

  Integer. Number of principal components in the model. Use
  [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  to choose the rank on held-out values.

- transform:

  Character or a transform object. The scale the model is fitted on:
  `"identity"`, `"logit_biwhiten"`, or an object from
  [`tf_identity()`](https://hhp94.github.io/slideimp/reference/tf_identity.md),
  [`tf_sv()`](https://hhp94.github.io/slideimp/reference/tf_identity.md)
  or
  [`tf_biwhiten()`](https://hhp94.github.io/slideimp/reference/tf_identity.md).
  See Details in `mipca()`.

- subset:

  Optional character or integer vector. The columns whose draws are kept
  and returned, by name or by number, or `NULL` for every column.
  Defaults to `NULL`.

- keep_state:

  Logical. If `TRUE`, return a fit that
  [`update()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
  can extend at the cost of memory. Defaults to `FALSE`.

- prop_diag:

  Numeric scalar between `0` and `1`. Proportion of each column's
  missing values on which
  [`summary()`](https://hhp94.github.io/slideimp/reference/summary.si_mi.md)
  computes Rhat and ESS. `0` is refused. Defaults to `0.5`.

- iter_sampling:

  Integer. Number of draws kept per chain. Defaults to `100`.

- iter_warmup:

  Integer. Number of iterations discarded at the start of each chain,
  `0` or more. Defaults to `100`.

- thin:

  Integer. After warmup, each chain keeps one iteration in every `thin`,
  so it runs `iter_sampling * thin` iterations to keep `iter_sampling`
  draws. Defaults to `1`.

- nu:

  Numeric. Degrees of freedom of the noise drawn for a missing value:
  `Inf` for Gaussian noise, or a number above `2` for Student-t noise of
  the same variance. Defaults to `Inf`. Can be chosen by
  [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md).

- chains:

  Integer. Number of independent chains. Defaults to `4`.

- refresh:

  Integer. Print a progress line every `refresh` iterations of each
  chain, or `0` for none. Defaults to `0`.

- spectra:

  A `spectra_opts` object. Settings for the eigensolver, from
  [`spectra_opts()`](https://hhp94.github.io/slideimp/reference/spectra_opts.md).

- init:

  Numeric. How widely each chain's starting values are spread: each
  missing value starts at a random draw whose standard deviation is
  `init` times its column's observed standard deviation. Must be
  positive. Defaults to `1`.

- seed:

  Integer, numeric, or `NULL`. Seed for the sampler. With `NULL`, the
  seed is drawn from R's random number stream. Defaults to `NULL`.

- abb:

  Optional character or integer vector. Columns to draw by approximate
  Bayesian bootstrap instead of by the model.
  [`abb_cols()`](https://hhp94.github.io/slideimp/reference/abb_cols.md)
  proposes them. Defaults to `NULL`.

## Value

With `keep_state = FALSE`, a finalized fit of class `si_mi`. It holds
`chains * iter_sampling` draws of every missing value in the selected
columns. Use:

- [`as.matrix(fit, draw = b)`](https://hhp94.github.io/slideimp/reference/as.matrix.si_mi.md)
  for the completed matrix at draw `b`.

- [`summary()`](https://hhp94.github.io/slideimp/reference/summary.si_mi.md)
  for Rhat and ESS.

- [`cbind()`](https://hhp94.github.io/slideimp/reference/cbind.si_mi.md)
  to join fits made over different columns.

With `keep_state = TRUE`, a resumable fit of class `si_mi_pca`. Extend
it with
[`update()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md),
then finish it with
[`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md).

## Details

**The transform.** `"identity"` fits the values as they are.
`"logit_biwhiten"` is for proportions such as DNAm beta values, and
every observed value must lie in `[0, 1]`; see
[`tf_sv()`](https://hhp94.github.io/slideimp/reference/tf_identity.md)
and
[`tf_biwhiten()`](https://hhp94.github.io/slideimp/reference/tf_identity.md).
Draws are returned on the scale of `obj`.

**Which columns take part.** Every column of `obj` takes part in the
model, and every missing value is drawn at every iteration. `subset`
selects only which columns' draws are kept. The result holds the
selected columns in the order they have in `obj`.

**Size.** The draws hold one number per missing value kept, per chain,
per draw. Keep fewer columns with `subset`, or ask for fewer `chains` or
`iter_sampling`. If the selected columns have no missing values, the fit
is returned without sampling.

**Reproducibility.** With the same `seed`, or the same
[`set.seed()`](https://rdrr.io/r/base/Random.html) before a call with
`seed = NULL`, the draws are the same.

## Performance

Linking R to a faster BLAS can make the sampler significantly faster:
[OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows) on
Windows, or Apple's Accelerate on macOS.

## References

Audigier V, Husson F, Josse J (2016). Multiple imputation for continuous
variables using a Bayesian principal component analysis. *Journal of
Statistical Computation and Simulation*, 86(11), 2140-2156.
[doi:10.1080/00949655.2015.1104683](https://doi.org/10.1080/00949655.2015.1104683)

## See also

- [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  and
  [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
  to choose `rank` and `nu` on held-out values, then fit.

- [`mipca_flank_plan()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
  to fit each target column on a window of its neighbors by location.

- [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md)
  for a single imputation.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)
fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", seed = 1)
fit
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
completed <- lapply(c(1, 200, 400), function(b) as.matrix(fit, draw = b))
summary(fit)
#> si_mi_summary: Rhat and ESS for 50 cells in 10 columns
#>   prop_diag 0.5, 4 chains x 100 draws
#>           0% 25%    50%    75%   100%
#> rhat   0.996   1   1.01   1.01   1.03
#> ess  153.000 291 330.00 354.00 409.00
#> worst 6 cells by rhat:
#>  chunk col      name row     rhat      ess
#>      1  10 feature10  11 1.029468 230.5473
#>      1   7  feature7   4 1.017100 333.5527
#>      1  14 feature14   2 1.016376 348.8917
#>      1  10 feature10  28 1.016294 219.3950
#>      1  20 feature20   5 1.016222 266.4496
#>      1  15 feature15  35 1.015993 300.7472
```
