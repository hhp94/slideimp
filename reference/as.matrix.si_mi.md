# Completed Matrix from a Bayesian PCA Imputation

Fill the missing values of `obj` with one draw from a finalized fit.

## Usage

``` r
# S3 method for class 'si_mi'
x[...]

# S3 method for class 'si_mi'
as.matrix(x, draw = 1L, ...)

# S3 method for class 'si_mi_pca'
as.matrix(x, ...)
```

## Arguments

- x:

  A `si_mi` object. A finalized fit, from
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md),
  [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
  or
  [`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md).

- ...:

  Not used.

- draw:

  Integer. The kept draw that fills the missing values, from `1` to the
  number of kept draws. Defaults to `1`.

## Value

A numeric matrix.

## Details

**The completed matrix.** A completed matrix is `obj` with its missing
values filled by one draw, on the scale of `obj`. The observed values
are unchanged, and the row and column names are kept. A fit made with
`subset` gives a matrix of the selected columns only.

**Which draws.** [`print()`](https://rdrr.io/r/base/print.html) shows
the number of kept draws. They are numbered chain by chain: with `k`
draws per chain, draws `1` to `k` are the first chain's, `k + 1` to
`2 * k` the second's, and so on.

**Multiple imputation.** Call
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) once per draw. Draws
evenly spaced from the first to the last are spread across the chains.
Each completed matrix is a full copy of `obj`, so `m` of them take `m`
times its memory. On a large `obj`, use each matrix in turn instead of
keeping all of them.

**Refusals.** A resumable fit, made with `keep_state = TRUE`, is
refused. Finalize it with
[`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
first. `[` does not work on a finalized fit: choose the columns with
`subset` when fitting, or subset a completed matrix.

## See also

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  make a finalized fit.

- [`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
  to finalize a resumable fit.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)
fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", seed = 1)
fit
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS

# The matrix completed by draw 10
x10 <- as.matrix(fit, draw = 10)

# Five completed matrices, at draws spread over the 400 kept
b <- round(seq(1, 400, length.out = 5))
completed <- lapply(b, function(d) as.matrix(fit, draw = d))
```
