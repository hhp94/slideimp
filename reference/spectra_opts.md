# Eigensolver Settings for Bayesian PCA Imputation

Settings for the eigensolver that
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md),
[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md),
[`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
and
[`mipca_flank()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
run at every iteration of the sampler. The defaults rarely need
changing.

## Usage

``` r
spectra_opts(tol = 1e-10, ncv = NULL, maxitr = 1000)
```

## Arguments

- tol:

  Numeric. Convergence tolerance of the eigenvalues, relative to their
  size. Defaults to `1e-10`.

- ncv:

  Integer or `NULL`. Number of Lanczos vectors the eigensolver keeps;
  more vectors converge in fewer iterations and use more memory.
  Defaults to `NULL`.

- maxitr:

  Integer. Maximum number of eigensolver iterations at each iteration of
  the sampler. Defaults to `1000`.

## Value

A `spectra_opts` object, to pass as `spectra`.

## Details

The eigensolver is Spectra's. These settings change how it converges.

With `ncv = NULL`, `ncv` is `2 * rank + 3`, at least `20` and at most
the smaller dimension of `obj` once the `abb` columns are left out. A
value given by hand must be more than `rank + 1` and at most that
dimension.

If a fit stops with "the PCA step did not converge", raise `maxitr` or
`ncv`.

## See also

[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md), whose
`spectra` argument this builds.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)

# More Lanczos vectors and iterations than the defaults
opts <- spectra_opts(ncv = 15, maxitr = 5000)
fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", spectra = opts, seed = 1)
fit
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
```
