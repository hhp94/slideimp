# Extend or Finalize a Resumable Bayesian PCA Fit

A fit made with `keep_state = TRUE` is resumable.
[`update()`](https://rdrr.io/r/stats/update.html) adds draws to each of
its chains. `mipca_finalize()` turns it into a finalized fit.

## Usage

``` r
# S3 method for class 'si_mi_pca'
update(object, iter_sampling, ...)

# S3 method for class 'si_mi'
update(object, ...)

mipca_finalize(x, thin = 1L, prop_diag = 0.5)
```

## Arguments

- object:

  A `si_mi_pca` object. A resumable fit, from
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) or
  [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
  with `keep_state = TRUE`.

- iter_sampling:

  Integer. Number of draws to add to each chain.

- ...:

  Not used.

- x:

  A `si_mi_pca` object. A resumable fit, from
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) or
  [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
  with `keep_state = TRUE`.

- thin:

  Integer. Keep every `thin`-th draw of each chain, starting at the
  first. Defaults to `1`.

- prop_diag:

  Numeric scalar between `0` and `1`. Proportion of each column's
  missing values on which
  [`summary()`](https://hhp94.github.io/slideimp/reference/summary.si_mi.md)
  computes Rhat and ESS. `0` is refused. Defaults to `0.5`.

## Value

- [`update()`](https://rdrr.io/r/stats/update.html): `object`,
  invisibly, extended in place.

- `mipca_finalize()`: a finalized fit of class `si_mi`, as
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md)
  returns with `keep_state = FALSE`.

## Details

**Extending.** [`update()`](https://rdrr.io/r/stats/update.html)
continues each chain from where it stopped, without a new warmup. The
rank, `nu`, transform, thinning, eigensolver settings and `subset` of
the fit stay as they were. A fit of `a` draws extended by `b` holds the
same draws, to rounding, as a fit of `a + b` draws made in one run with
the same `seed`.

**In place.** [`update()`](https://rdrr.io/r/stats/update.html) changes
`object` itself and returns it invisibly, so
`update(fit, iter_sampling = 100)` needs no assignment. `b <- fit` does
not copy a resumable fit: `b` and `fit` are the same fit, and extending
one extends both.

**Finalizing.** `mipca_finalize()` keeps every `thin`-th draw and
computes Rhat and ESS on `prop_diag` of each column's missing values.
Use the `si_mi` fit it returns: the resumable fit can no longer be used.

**Failures.** If [`update()`](https://rdrr.io/r/stats/update.html) fails
partway, the fit may no longer be usable, and
[`print()`](https://rdrr.io/r/base/print.html) says so. Fit again with
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).

**Finalized fits.** A `si_mi` fit holds no saved state.
[`update()`](https://rdrr.io/r/stats/update.html) and `mipca_finalize()`
refuse it.

## See also

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  make a resumable fit with `keep_state = TRUE`.

- [`as.matrix()`](https://hhp94.github.io/slideimp/reference/as.matrix.si_mi.md)
  for the completed matrices of the finalized fit.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)
fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", keep_state = TRUE, seed = 1)
fit
#> si_mi_pca fit, resumable: 50 x 20, rank 2, 4 chains x 100 draws
#>   next: update() to extend, or mipca_finalize() to finish.
#>         summary() shows Rhat and ESS per cell.

# 100 more draws per chain, in place
update(fit, iter_sampling = 100)
fit
#> si_mi_pca fit, resumable: 50 x 20, rank 2, 4 chains x 200 draws
#>   next: update() to extend, or mipca_finalize() to finish.
#>         summary() shows Rhat and ESS per cell.

# Keep every second draw
final <- mipca_finalize(fit, thin = 2)
final
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
```
