# Choose the Rank of Bayesian PCA Imputation on Held-Out Values

Hide a set of observed values, fit
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) at each
candidate rank, and score how well each rank's draws predict the hidden
values. Pass the result to
[`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
to fit the best rank.

## Usage

``` r
mipca_tune(
  obj,
  rank,
  na_loc = NULL,
  n_cols = NULL,
  n_rows = 4L,
  transform,
  n_diag = 10000L,
  iter_warmup = 100L,
  iter_sampling = 100L,
  chains = 4L,
  thin = 1,
  init = 1,
  refresh = 0L,
  spectra = spectra_opts(),
  seed = NULL,
  abb = NULL
)

# S3 method for class 'si_tune'
print(x, ...)
```

## Arguments

- obj:

  A double matrix with samples in rows and features in columns. Row and
  column names are required and must be unique. `NA` and `NaN` values
  are imputed.

- rank:

  Integer vector. The candidate ranks, each a whole number of at least
  `1`, given once.

- na_loc:

  A list or `NULL`. The values to hold out, as one two-column matrix of
  row and column positions in a list, the form
  [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  returns with `n_reps = 1`. Each must be an observed value of `obj`.
  With `NULL`, the positions are drawn from `n_cols` and `n_rows`.
  Defaults to `NULL`.

- n_cols:

  Integer or `NULL`. Number of columns that receive held-out values.
  With `NULL`, 10% of the columns, at least `50` and at most `500`.
  Defaults to `NULL`.

- n_rows:

  Integer. Number of held-out values in each of the `n_cols` columns.
  Leave it out when `na_loc` is given. Defaults to `4`.

- transform:

  Character or a transform object. The scale the model is fitted on:
  `"identity"`, `"logit_biwhiten"`, or an object from
  [`tf_identity()`](https://hhp94.github.io/slideimp/reference/tf_identity.md),
  [`tf_sv()`](https://hhp94.github.io/slideimp/reference/tf_identity.md)
  or
  [`tf_biwhiten()`](https://hhp94.github.io/slideimp/reference/tf_identity.md).
  See Details in
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).

- n_diag:

  Integer. Number of missing values, held-out values included, on which
  Rhat and ESS are computed for each rank. They are sampled across
  columns with few and with many missing values. Defaults to `10000`.

- iter_warmup:

  Integer. Number of iterations discarded at the start of each chain,
  `0` or more. Defaults to `100`.

- iter_sampling:

  Integer. Number of draws kept per chain. Defaults to `100`.

- chains:

  Integer. Number of independent chains. Defaults to `4`.

- thin:

  Integer. After warmup, each chain keeps one iteration in every `thin`,
  so it runs `iter_sampling * thin` iterations to keep `iter_sampling`
  draws. Defaults to `1`.

- init:

  Numeric. How widely each chain's starting values are spread: each
  missing value starts at a random draw whose standard deviation is
  `init` times its column's observed standard deviation. Must be
  positive. Defaults to `1`.

- refresh:

  Integer. Print a progress line every `refresh` iterations of each
  chain, or `0` for none. Defaults to `0`.

- spectra:

  A `spectra_opts` object. Settings for the eigensolver, from
  [`spectra_opts()`](https://hhp94.github.io/slideimp/reference/spectra_opts.md).

- seed:

  Integer, numeric, or `NULL`. Seed for the held-out values and the
  sampler. With `NULL`, both are drawn from R's random number stream.
  Defaults to `NULL`.

- abb:

  Optional character or integer vector. Columns to draw by approximate
  Bayesian bootstrap instead of by the model.
  [`abb_cols()`](https://hhp94.github.io/slideimp/reference/abb_cols.md)
  proposes them. Defaults to `NULL`.

- x:

  A `si_tune` object.

- ...:

  Not used.

## Value

A data frame of class `si_tune`, one row per rank, with columns:

- `rank`: the rank.

- `crps`: the mean score of the held-out values. Lower is better.

- `implied_nu`, `implied_lo`, `implied_hi`: the degrees of freedom the
  held-out values imply, and its 68% bootstrap interval.

- `robust_z`: the interquartile range of the standardized held-out
  values, divided by `1.349`. Near `1` for Gaussian noise.

- `ess_min`, `ess_p10`, `ess_med`: the minimum, 10% quantile and median
  ESS over the `n_diag` values.

- `ess_ok`: the proportion of them with an ESS of at least
  `100 * chains`.

- `rhat_med`, `rhat_p90`, `rhat_max`: the median, 90% quantile and
  maximum Rhat.

- `rhat_ok`: the proportion with an Rhat of at most `1.05`.

- `ndraws`: the number of draws, over all chains.

- `n_diag`: the number of values Rhat and ESS were computed on.

- `sigma`: a list column, the noise standard deviation at the end of
  each chain.

- `time`: the seconds the rank took.

- `error`: the message of a failed fit, or `NA`.

- `result`: a list column, a data frame with one row per held-out value:
  its `row` and `col`, the `truth`, the `mean` and `sd` of its draws,
  its `crps` and `pit`, and `truth_lat`, `mean_lat` and `sd_lat` on the
  scale the model is fitted on.

- `diag`: a list column, a data frame of the `row`, `col`, `rhat` and
  `ess` of each of the `n_diag` values.

The table also carries the transform, the sampler settings, `abb` and
the held-out positions, which
[`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
reads.

## Details

**The held-out values.** With `na_loc = NULL`,
[`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
picks `n_cols` columns and hides `n_rows` observed values in each. A
column is eligible if it is not constant and has at least `n_rows + 2`
observed values. Every rank is fitted on the same matrix, with the same
values hidden, and is scored on those values only. The holdout needs
more than two values. The `abb` columns take no part in tuning, and no
value in them is held out.

**The rank.** `crps` is the continuous ranked probability score of each
held-out value against its draws, on the scale of `obj`, averaged over
the held-out values. A lower score is better.
[`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
fits the rank with the lowest `crps` unless it is given another.

**The noise.** `implied_nu` is the Student-t degrees of freedom that the
spread of the held-out values implies, `Inf` when they look Gaussian.
`implied_lo` and `implied_hi` are its bootstrap interval. `implied_nu`
is the raw reading.
[`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
fits with `nu` set to `implied_nu`, or to `5` if the reading is lower.

**Failures.** A rank whose fit fails does not stop the call. Its row
holds the message in `error`, and its scores are `NA`.

**Reproducibility.** With the same `seed`, the held-out values and the
table are the same. A given `seed` leaves R's random number stream as it
was. With `seed = NULL`, the same
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call gives
the same table.

## See also

- [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
  to fit a rank from this table.

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  fit with a rank and `nu` of your own.

- [`sample_na_loc()`](https://hhp94.github.io/slideimp/reference/sample_na_loc.md)
  to draw held-out positions to pass as `na_loc`.

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
tuned
#> si_tune: 3 ranks on a 50 x 20 matrix, holdout 40 cells over 10 columns
#>  rank    crps ess_p10  ess_ok rhat_max implied_nu implied_lo implied_hi  time
#>     1 0.09652   267.0 0.07143    1.032      2.600      2.324      4.215 0.035
#>     2 0.09713   247.1 0.02857    1.040      2.842      2.366      4.286 0.036
#>     3 0.10109   172.0 0.01429    1.052      3.414      2.726      6.308 0.032
#>   next: mipca_fit(x, obj) fits the rank with the lowest crps.
#>         mipca_fit(x, obj, rank) fits another rank.
fit <- mipca_fit(tuned, sim$input, seed = 1)
#> the implied nu at rank 1 is 2.6, under the floor of 5, so the fit uses nu = 5. Pass 'nu' to use another value
#> Chain 1 Iteration:   1 / 200 [  0%] (Warmup)
#> Chain 1 Iteration: 100 / 200 [ 50%] (Warmup)
#> Chain 1 Iteration: 200 / 200 [100%] (Sampling)
#> Chain 1 finished in 0.0 seconds.
#> Chain 2 Iteration:   1 / 200 [  0%] (Warmup)
#> Chain 2 Iteration: 100 / 200 [ 50%] (Warmup)
#> Chain 2 Iteration: 200 / 200 [100%] (Sampling)
#> Chain 2 finished in 0.0 seconds.
#> Chain 3 Iteration:   1 / 200 [  0%] (Warmup)
#> Chain 3 Iteration: 100 / 200 [ 50%] (Warmup)
#> Chain 3 Iteration: 200 / 200 [100%] (Sampling)
#> Chain 3 finished in 0.0 seconds.
#> Chain 4 Iteration:   1 / 200 [  0%] (Warmup)
#> Chain 4 Iteration: 100 / 200 [ 50%] (Warmup)
#> Chain 4 Iteration: 200 / 200 [100%] (Sampling)
#> Chain 4 finished in 0.0 seconds.
fit
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
```
