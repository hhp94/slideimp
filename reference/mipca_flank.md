# Bayesian PCA Imputation on a Flanking Window per Column

Impute each target column from the columns near it by `location`, such
as the CpGs within some distance of a CpG. `mipca_flank_plan()` finds
each target's window, and `mipca_flank()` fits
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) on each
window and keeps the target's draws.

## Usage

``` r
mipca_flank_plan(
  obj,
  location,
  window_size,
  subset,
  min_window_n,
  colmax = 0.95,
  abb = NULL
)

# S3 method for class 'si_flank_plan'
print(x, n = 10L, ...)

mipca_flank_window(plan, subset)

mipca_flank(
  plan,
  rank,
  seed,
  transform,
  prop_diag = 0.5,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  spectra = spectra_opts(),
  init = 1,
  .progress = TRUE
)
```

## Arguments

- obj:

  A double matrix with samples in rows and features in columns. Row and
  column names are required and must be unique. `NA` and `NaN` values
  are imputed.

- location:

  Numeric. The position of each column of `obj`, such as a genomic
  coordinate, sorted in ascending order.

- window_size:

  Numeric. Distance from the target, in the units of `location`, that
  its window covers on each side.

- subset:

  Character or integer vector. In `mipca_flank_plan()`, the target
  columns. In `mipca_flank_window()`, one target of the plan. By name or
  by number.

- min_window_n:

  Integer. Minimum number of columns in a window, the target included,
  at least `2`. `mipca_flank()` needs at least `rank + 2`.

- colmax:

  Numeric scalar between `0` and `1`. Columns missing a larger
  proportion of their values are left out of every window, and are not
  imputed. Defaults to `0.95`.

- abb:

  Optional character or integer vector. Columns to draw by approximate
  Bayesian bootstrap instead of by the model. They are left out of every
  window.
  [`abb_cols()`](https://hhp94.github.io/slideimp/reference/abb_cols.md)
  proposes them. Defaults to `NULL`.

- x:

  A `si_flank_plan` object.

- n:

  Integer. Number of rows of the table to print. Defaults to `10`.

- ...:

  Not used.

- plan:

  A `si_flank_plan` object. The plan `mipca_flank_plan()` returns.

- rank:

  Integer. Number of principal components in the model. Use
  [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  to choose the rank on held-out values.

- seed:

  Integer. Seed for the sampler.

- transform:

  Character or a transform object. The scale the model is fitted on:
  `"identity"`, `"logit_biwhiten"`, or an object from
  [`tf_identity()`](https://hhp94.github.io/slideimp/reference/tf_identity.md),
  [`tf_sv()`](https://hhp94.github.io/slideimp/reference/tf_identity.md)
  or
  [`tf_biwhiten()`](https://hhp94.github.io/slideimp/reference/tf_identity.md).
  See Details in
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).

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

- spectra:

  A `spectra_opts` object. Settings for the eigensolver, from
  [`spectra_opts()`](https://hhp94.github.io/slideimp/reference/spectra_opts.md).

- init:

  Numeric. How widely each chain's starting values are spread: each
  missing value starts at a random draw whose standard deviation is
  `init` times its column's observed standard deviation. Must be
  positive. Defaults to `1`.

- .progress:

  Logical. If `TRUE`, report progress through the windows. Defaults to
  `TRUE`.

## Value

- `mipca_flank_plan()`: a `si_flank_plan` object, holding `obj`, the
  settings and the table `targets`.

- `mipca_flank()`: a finalized fit of class `si_mi` over the targets, as
  [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md)
  returns.

- `mipca_flank_window()`: a numeric matrix, the columns of the target's
  window.

- [`print()`](https://rdrr.io/r/base/print.html): `x`, invisibly.

## Details

**The window.** A target's window is the columns within `window_size` of
it on `location`, on both sides. Columns that are constant, have fewer
than two observed values, miss more than `colmax` of their values, or
are listed in `abb` are left out.

**The columns read.** `mipca_flank_plan()` reads a column only if a
window can hold it, or if it is a target listed in `abb`, and refuses an
`Inf` in any column it reads. It does not read the other columns, so an
`Inf` in one of them is not refused. `mipca_flank()` reads only the
targets listed in `abb` and the columns of the windows it fits.

**The plan.** `plan$targets` has one row per target, in column order:
its column number `target` and `name`, its number of missing values
`n_miss`, the first and last column of its window `start` and `end`, the
number of columns in the window `window_n`, and its `status`:

- `"impute"`: fitted on its window.

- `"clean"`: no missing values.

- `"abb"`: listed in `abb`, and drawn from its own observed values.

- `"inadmissible"`: constant, with fewer than two observed values, or
  missing more than `colmax`. Its missing values stay `NA`.

- `"window_small"`: fewer than `min_window_n` columns in its window. Its
  missing values stay `NA`.

A plan in which no target with missing values can be imputed is refused.
[`print()`](https://rdrr.io/r/base/print.html) shows the count per
status and the first `n` rows.

**The rank and `nu`.** Choose them on one target's window: pass
`mipca_flank_window(plan, target)` to
[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md),
then give the rank and `nu` to `mipca_flank()`. The default `n_cols` of
[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
holds out values in at least 50 columns, more than a window may have, so
pass a smaller `n_cols` there.

**The fit.** The result holds the targets in the order of
`plan$targets`. Its completed matrices hold `NA` where a target was not
imputed, and [`print()`](https://rdrr.io/r/base/print.html) counts those
columns. A target's draws depend on `seed` and its column number in
`obj` only, so dropping another target from the plan leaves them
unchanged.

**Failures.** A window whose fit fails stops the run, with the target's
name in the message.

## See also

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  fit the whole matrix at once.

- [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  to choose the rank on a window.

- [`cbind()`](https://hhp94.github.io/slideimp/reference/cbind.si_mi.md)
  to join the fits of several matrices, such as one per chromosome.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 200, beta = TRUE)
obj <- sim$input
location <- seq_len(ncol(obj))

# Windows of 20 columns on each side of three targets
plan <- mipca_flank_plan(
  obj,
  location,
  window_size = 20,
  subset = c(49, 102, 150),
  min_window_n = 10
)
plan
#> si_flank_plan: 3 targets on a 50 x 200 matrix, nothing fitted
#>   status: impute 3
#>  target       name n_miss start end window_n status
#>      49  feature49     10    29  69       41 impute
#>     102 feature102     10    82 122       41 impute
#>     150 feature150     10   130 170       41 impute
#>   next: mipca_flank(x, rank, seed, transform) imputes.
#>         mipca_flank_window(x, subset) gives one target's window to tune on.

# One target's window, to choose the rank on
win <- mipca_flank_window(plan, 150)
dim(win)
#> [1] 50 41

fit <- mipca_flank(plan, rank = 2, seed = 1, transform = "logit_biwhiten")
#> processing window 1 of 3
#> processing window 2 of 3
#> processing window 3 of 3
fit
#> si_mi fit: 50 x 3, 30 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
```
