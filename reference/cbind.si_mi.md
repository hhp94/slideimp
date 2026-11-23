# Join Bayesian PCA Imputations of Different Columns

Join finalized fits made on the same rows over different columns into
one fit, such as a large matrix fitted one chromosome at a time.

## Usage

``` r
# S3 method for class 'si_mi'
cbind(..., deparse.level = 1)
```

## Arguments

- ...:

  One or more `si_mi` objects. Finalized fits, from
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md),
  [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
  or
  [`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md).

- deparse.level:

  Not used.

## Value

A finalized fit of class `si_mi`, over the columns of every fit. A
single fit is returned as it is.

## Details

**The joined fit.** The columns come in the order the fits are given.
Draw `b` of the joined fit is draw `b` of every fit, so its completed
matrix is the completed matrices of the fits at draw `b`, side by side.

**What must match.** The fits must have the same row names in the same
order and the same number of kept draws, so fit and finalize each with
the same `chains`, `iter_sampling` and `thin`. A column can be in one
fit only. The rank, transform and other settings can differ.

**Seeds.** Give each fit its own `seed`, to reduce the correlation
between the draws of different fits.

**Refusals.** A resumable fit, made with `keep_state = TRUE`, is
refused. Finalize it with
[`mipca_finalize()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
first.

## See also

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  fit each part.

- [`as.matrix()`](https://hhp94.github.io/slideimp/reference/as.matrix.si_mi.md)
  and
  [`summary()`](https://hhp94.github.io/slideimp/reference/summary.si_mi.md),
  which take the joined fit as any other.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 40, beta = TRUE)
obj <- sim$input

# Fit two blocks of columns, each with its own seed, and join them
left <- mipca(obj[, 1:20], rank = 2, transform = "logit_biwhiten", seed = 1)
right <- mipca(obj[, 21:40], rank = 2, transform = "logit_biwhiten", seed = 2)
fit <- cbind(left, right)
fit
#> si_mi fit: 50 x 40, 200 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
dim(as.matrix(fit, draw = 1))
#> [1] 50 40
```
