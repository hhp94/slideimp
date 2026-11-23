# Summarize and Print a Bayesian PCA Imputation

[`summary()`](https://rdrr.io/r/base/summary.html) returns Rhat and ESS
for a sample of the missing values of a fit.
[`print()`](https://rdrr.io/r/base/print.html) shows a fit or a summary.

## Usage

``` r
# S3 method for class 'si_mi_pca'
print(x, ...)

# S3 method for class 'si_mi'
print(x, ...)

# S3 method for class 'si_mi_pca'
summary(object, prop_diag = 0.5, ...)

# S3 method for class 'si_mi'
summary(object, ...)

# S3 method for class 'si_mi_summary'
print(x, n = 6L, ...)
```

## Arguments

- x:

  A `si_mi`, `si_mi_pca` or `si_mi_summary` object.

- ...:

  Not used.

- object:

  A `si_mi` or `si_mi_pca` object. A fit from
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) or
  [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md).

- prop_diag:

  Numeric scalar between `0` and `1`. Proportion of each column's
  missing values on which Rhat and ESS are computed. `0` is refused.
  Only a resumable fit takes it. Defaults to `0.5`.

- n:

  Integer. Number of values with the highest Rhat to list. `0` lists
  none. Defaults to `6`.

## Value

- [`summary()`](https://rdrr.io/r/base/summary.html): an object of class
  `si_mi_summary`. Its attributes `prop_diag`, `nchains` and `ndraws`
  hold the proportion, the number of chains and the draws per chain, or
  `NA` where joined fits differ. Its element `cells` is a data frame
  with one row per value, in column order:

  - `chunk`: the fit the value comes from, `1` unless fits were joined
    with
    [`cbind()`](https://hhp94.github.io/slideimp/reference/cbind.si_mi.md).
    Finalized fits only.

  - `col`, `name`: the column number and name in the completed matrix.

  - `row`: the row number.

  - `rhat`, `ess`: Rhat and ESS.

- [`print()`](https://rdrr.io/r/base/print.html): `x`, invisibly.

## Details

**Rhat and ESS.** Rhat is the folded split Rhat of Vehtari et al.
(2021). ESS is computed over all chains together.
[`update()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
adds draws to a resumable fit, which raises ESS.

**Which values.** A column with `h` missing values gives
`ceiling(prop_diag * h)` of them, spread evenly along the column. The
sample is the same on every call. Columns listed in `abb` are left out.

**Finalized and resumable fits.** A finalized fit returns the values
computed when it was finalized, with the `prop_diag` given to
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) or
[`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md).
A resumable fit computes them on each call, with the `prop_diag` given
here.

**Short runs.** Rhat is `NA` with fewer than `4` draws per chain, and
ESS is `NA` with fewer than `12`.

**Joined fits.** For fits joined with
[`cbind()`](https://hhp94.github.io/slideimp/reference/cbind.si_mi.md),
`chunk` gives the fit each value comes from.
[`print()`](https://rdrr.io/r/base/print.html) shows `mixed` where the
fits differ in `prop_diag`, chains or draws.

**Printing.** [`print()`](https://rdrr.io/r/base/print.html) of a fit
shows its size and the calls to make next.
[`print()`](https://rdrr.io/r/base/print.html) of a summary shows the
number of values, the quantiles of `rhat` and `ess`, and the `n` values
with the highest Rhat.

## References

Vehtari A, Gelman A, Simpson D, Carpenter B, Burkner P-C (2021).
Rank-normalization, folding, and localization: an improved R-hat for
assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)

## See also

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  fit, and choose `prop_diag`.

- [`update()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
  to add draws to a resumable fit.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)
fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", seed = 1)
s <- summary(fit)
s
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
print(s, n = 3)
#> si_mi_summary: Rhat and ESS for 50 cells in 10 columns
#>   prop_diag 0.5, 4 chains x 100 draws
#>           0% 25%    50%    75%   100%
#> rhat   0.996   1   1.01   1.01   1.03
#> ess  153.000 291 330.00 354.00 409.00
#> worst 3 cells by rhat:
#>  chunk col      name row     rhat      ess
#>      1  10 feature10  11 1.029468 230.5473
#>      1   7  feature7   4 1.017100 333.5527
#>      1  14 feature14   2 1.016376 348.8917
head(s$cells)
#>   chunk col     name row     rhat      ess
#> 1     1   3 feature3   3 1.001578 377.5531
#> 2     1   3 feature3  13 1.007004 378.1096
#> 3     1   3 feature3  19 1.014247 291.7192
#> 4     1   3 feature3  38 1.000465 328.7449
#> 5     1   3 feature3  49 1.005982 322.9542
#> 6     1   4 feature4   8 1.011456 234.8404

# A resumable fit takes prop_diag at each call
fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", keep_state = TRUE, seed = 1)
summary(fit, prop_diag = 1)
#> si_mi_summary: Rhat and ESS for 100 cells in 10 columns
#>   prop_diag 1, 4 chains x 100 draws
#>           0% 25%    50%    75%   100%
#> rhat   0.994   1   1.01   1.01   1.03
#> ess  153.000 290 323.00 354.00 442.00
#> worst 6 cells by rhat:
#>  col      name row     rhat      ess
#>   10 feature10  11 1.029468 230.5473
#>    7  feature7  46 1.023850 168.3942
#>   17 feature17  26 1.018449 278.6523
#>    7  feature7   4 1.017100 333.5527
#>   14 feature14   2 1.016376 348.8917
#>   10 feature10  28 1.016294 219.3950
```
