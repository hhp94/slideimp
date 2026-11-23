# Fit Bayesian PCA Imputation at a Tuned Rank

Fit [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) at
a rank from a
[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
table, with the `nu`, transform and sampler settings of that tuning run.

## Usage

``` r
mipca_fit(
  object,
  obj,
  rank = NULL,
  ...,
  chains,
  iter_warmup,
  iter_sampling,
  thin,
  spectra,
  init,
  refresh = 100L,
  keep_state = FALSE,
  prop_diag = 0.5,
  seed = NULL,
  nu = NULL,
  transform,
  abb
)
```

## Arguments

- object:

  A `si_tune` object. The table
  [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  returns, or a subset of its rows.

- obj:

  A double matrix. The matrix `object` was tuned on, passed again. The
  table does not carry it.

- rank:

  Integer or `NULL`. A rank in `object`. With `NULL`, the rank with the
  lowest `crps`. Defaults to `NULL`.

- ...:

  Not used. Name every argument after `rank`.

- chains:

  Integer. Number of independent chains. Taken from `object` when left
  out.

- iter_warmup:

  Integer. Number of iterations discarded at the start of each chain,
  `0` or more. Taken from `object` when left out.

- iter_sampling:

  Integer. Number of draws kept per chain. Taken from `object` when left
  out.

- thin:

  Integer. After warmup, each chain keeps one iteration in every `thin`,
  as in
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).
  Taken from `object` when left out.

- spectra:

  A `spectra_opts` object. Settings for the eigensolver, from
  [`spectra_opts()`](https://hhp94.github.io/slideimp/reference/spectra_opts.md).
  Taken from `object` when left out.

- init:

  Numeric. How widely each chain's starting values are spread, as in
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).
  Taken from `object` when left out.

- refresh:

  Integer. Print a progress line every `refresh` iterations of each
  chain, or `0` for none. Defaults to `100`.

- keep_state:

  Logical. If `TRUE`, return a fit that
  [`update()`](https://hhp94.github.io/slideimp/reference/mipca_finalize.md)
  can extend at the cost of memory. Defaults to `FALSE`.

- prop_diag:

  Numeric scalar between `0` and `1`. Proportion of each column's
  missing values on which
  [`summary()`](https://hhp94.github.io/slideimp/reference/summary.si_mi.md)
  computes Rhat and ESS. `0` is refused. Defaults to `0.5`.

- seed:

  Integer, numeric, or `NULL`. Seed for the sampler. With `NULL`, the
  seed is drawn from R's random number stream. Defaults to `NULL`.

- nu:

  Numeric or `NULL`. Degrees of freedom of the noise drawn for a missing
  value, as in
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md). With
  `NULL`, the `implied_nu` of `rank` in `object`, or `5` if it is lower.
  Defaults to `NULL`.

- transform, abb:

  Not used. The fit uses the transform and `abb` stored in `object`. To
  change them, run
  [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  again, or call
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).

## Value

As [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md):
with `keep_state = FALSE`, a finalized fit of class `si_mi`. With
`keep_state = TRUE`, a resumable fit of class `si_mi_pca`.

## Details

**The matrix.** `obj` must be the same matrix `object` was tuned on.
Another matrix is refused. The held-out values are not held out in the
fit: every observed value of `obj` is used.

**The rank and `nu`.** With `rank = NULL`, the rank is the one with the
lowest `crps` among the ranks whose fit did not fail. A given `rank`
must be in `object`. When the `implied_nu` is below `5`, the fit uses
`5` and a message says so. A `nu` that is given is used as it is, with
the warning
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) gives
below `5`.

**The same fit as
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).** The
fit is the one
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) returns
for the same `obj`, `rank`, `nu`, transform, `abb`, sampler settings and
`seed`.

## See also

- [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  to make the table.

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  fit with a rank and `nu` of your own.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)
tuned <- mipca_tune(
  sim$input,
  rank = 1:3,
  transform = "logit_biwhiten",
  n_cols = 10,
  seed = 1
)

# The rank with the lowest crps
fit <- mipca_fit(tuned, sim$input, refresh = 0, seed = 1)
#> the implied nu at rank 1 is 2.6, under the floor of 5, so the fit uses nu = 5. Pass 'nu' to use another value
fit
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS

# Another tuned rank, with more draws per chain
fit3 <- mipca_fit(tuned, sim$input, rank = 3, iter_sampling = 200, refresh = 0, seed = 1)
#> the implied nu at rank 3 is 3.41, under the floor of 5, so the fit uses nu = 5. Pass 'nu' to use another value
```
