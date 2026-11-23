# Multiple imputation of a beta matrix

``` r

library(slideimp)
```

> **Experimental**
>
> The Bayesian PCA sampler
> ([`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md),
> [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md),
> [`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
> and
> [`mipca_flank()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md))
> is experimental. Its arguments, defaults and output can change in a
> later release.

[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md) and
[`knn_imp()`](https://hhp94.github.io/slideimp/reference/knn_imp.md)
fill each missing value with one best guess. That is enough for many
uses, but an analysis run on a singly imputed matrix treats the guesses
as if they had been measured, so its standard errors come out too small.
Multiple imputation instead draws several plausible values for every
missing value, so the uncertainty about the missing cells carries
through to the final estimate.

slideimp does this with a Bayesian low-rank (PCA) model, the method of
[Audigier, Husson and Josse
(2016)](https://doi.org/10.1080/00949655.2015.1104683). This article
goes through each step for a beta-value matrix (DNA methylation, values
in `[0, 1]`), from tuning the model to a pooled estimate.

## Simulate a beta matrix

The matrix below is built on the logit scale - a mean per CpG with the
two modes typical of methylation, a rank-2 signal, and heavy-tailed
noise - and then mapped into `[0, 1]` with
[`plogis()`](https://rdrr.io/r/stats/Logistic.html). Half of the columns
have 16 of their 80 samples set to `NA`.

This data follows the model the sampler assumes, so the fit will look
good. It is here to show the steps, not to show how well the method
works on real data.

``` r

set.seed(2)
n <- 80
p <- 300
mu <- ifelse(runif(p) < 0.5, rnorm(p, -2, 0.7), rnorm(p, 2, 0.7))
U <- matrix(rnorm(n * 2), n, 2)
V <- matrix(rnorm(p * 2, 0, 0.5), p, 2)
noise <- 0.4 * sqrt(4 / 6) * matrix(rt(n * p, df = 6), n, p)
truth <- plogis(matrix(mu, n, p, byrow = TRUE) + U %*% t(V) + noise)
dimnames(truth) <- list(sprintf("s%02d", seq_len(n)), sprintf("cpg%03d", seq_len(p)))

obj <- truth
na_cols <- sample.int(p, p / 2)
for (j in na_cols) obj[sample.int(n, 16), j] <- NA
obj[1:4, 1:4]
#>         cpg001    cpg002    cpg003     cpg004
#> s01 0.21050697 0.9090441 0.9077718 0.06148590
#> s02 0.13928224 0.7599630 0.7580556 0.07896476
#> s03 0.21884599 0.8688346 0.9009699 0.05805358
#> s04 0.04057579 0.8940313 0.7209590 0.09907192
mean(is.na(obj))
#> [1] 0.1
```

Row and column names are required: the fit uses them to keep track of
which cell is which.

## Tune the model with `mipca_tune()`

Two choices shape the model: the `rank` (how many components) and `nu`,
the degrees of freedom of the noise. `Inf` gives Gaussian noise, and
smaller values allow heavier tails.

[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
hides some observed cells, imputes them at each candidate rank, and
scores the draws against the hidden truth with the continuous ranked
probability score (CRPS, lower is better). From the same draws it also
estimates the `nu` that each rank implies, so one call informs both
choices.

`transform = "logit_biwhiten"` is the transform for beta values. It maps
the data onto the logit scale (after squeezing exact 0 and 1 inwards so
the logit stays finite) and then rescales rows and columns. The sampler
works on that scale, and every value it hands back is mapped back into
`[0, 1]`. The argument has no default: pass `transform = "identity"` for
data that is already roughly Gaussian.

``` r

tuned <- mipca_tune(
  obj,
  rank = 1:4,
  transform = "logit_biwhiten",
  seed = 1
)
tuned
#> si_tune: 4 ranks on a 80 x 300 matrix, holdout 200 cells over 50 columns
#>  rank    crps ess_p10 ess_ok rhat_max implied_nu implied_lo implied_hi  time
#>     1 0.04873   323.3 0.3365    1.032        Inf        Inf        Inf 0.252
#>     2 0.02606   305.4 0.2554    1.036      9.587      6.334        Inf 0.256
#>     3 0.02598   290.6 0.1992    1.039     13.679      4.561        Inf 0.377
#>     4 0.02618   277.6 0.1369    1.039      9.954      5.475        Inf 0.422
#>   next: mipca_fit(x, obj) fits the rank with the lowest crps.
#>         mipca_fit(x, obj, rank) fits another rank.
```

How to read the table:

- `crps`: the mean score on the hidden cells. The best rank has the
  lowest `crps`.
- `ess_p10` and `ess_ok`: the 10th percentile of the effective sample
  size (ESS) over the values that Rhat and ESS are computed on, and the
  share of them with an ESS of at least 100 x `chains` (400 here, with
  the default 4 chains). Low values mean the chains are short, not that
  the rank is wrong. Run
  [`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
  again with a larger `iter_sampling` for longer chains.
- `rhat_max`: the worst Rhat over the same values. Values well above
  1.05 mean the chains disagree.
- `implied_nu`, `implied_lo` and `implied_hi`: the noise degrees of
  freedom that the hidden cells imply at that rank, and its 68%
  bootstrap interval.
- `time`: the seconds each rank took.

The table holds more columns than it prints.
[`?mipca_tune`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
lists them all. The `result` column holds the score of each hidden cell,
per rank:

``` r

head(tuned$result[[1]][, c("row", "col", "truth", "mean", "sd", "crps")])
#>   row col      truth       mean         sd       crps
#> 1   1   2 0.90904408 0.85357938 0.05664707 0.02975744
#> 2  15   2 0.89156589 0.92327401 0.05471297 0.02464390
#> 3  38   2 0.90441911 0.83230465 0.05491061 0.04333744
#> 4  60   2 0.90403138 0.89717857 0.05166960 0.01188086
#> 5  38  30 0.15984324 0.11436739 0.06370078 0.03423272
#> 6  40  30 0.08271408 0.04613398 0.04426940 0.02868168
```

## Fit the model with `mipca_fit()`

`mipca_fit(tuned, obj)` runs the final fit at the rank with the lowest
CRPS, with `nu` taken from that rank’s estimate and the chain settings
carried over from the tuner. Pass the matrix again: the tuning table
does not hold it. An implied `nu` below 5 is raised to 5, with a
message, because the sampler can become unstable below 5. The transform
is refitted on the full matrix.

``` r

fit <- mipca_fit(tuned, obj, seed = 2, refresh = 0)
fit
#> si_mi fit: 80 x 300, 2400 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
```

With the tuner’s default settings, this runs 4 chains with 100 warm-up
and 100 sampling iterations each, so 400 draws are kept for every
missing value.

To fit a rank and `nu` that do not come from a tuning table, call
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md)
directly.
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) also
takes `subset`, which limits which columns have their draws stored.
Every missing value is still drawn at every iteration, because the
low-rank structure ties the columns together, so `subset` saves memory,
not time. The stored draws take cells x draws of memory, so on a
genome-scale matrix the saving is large.

## Check convergence

[`summary()`](https://rdrr.io/r/base/summary.html) reports Rhat and ESS
for a sample of the imputed cells. Rhat compares the chains with each
other ([Gelman and Rubin 1992](https://doi.org/10.1214/ss/1177011136),
[Vehtari et al. 2021](https://doi.org/10.1214/20-BA1221)). ESS, the
effective sample size, estimates how many independent draws the chains
are worth ([Geyer 1992](https://doi.org/10.1214/ss/1177011137), [Vehtari
et al. 2021](https://doi.org/10.1214/20-BA1221)).

``` r

summary(fit)
#> si_mi_summary: Rhat and ESS for 1200 cells in 150 columns
#>   prop_diag 0.5, 4 chains x 100 draws
#>           0% 25% 50%    75%   100%
#> rhat   0.992   1   1   1.01   1.04
#> ess  119.000 325 360 393.00 536.00
#> worst 6 cells by rhat:
#>  chunk col   name row     rhat      ess
#>      1 137 cpg137  11 1.040830 118.6063
#>      1 120 cpg120  58 1.029801 220.7080
#>      1  59 cpg059  38 1.027809 200.7972
#>      1 118 cpg118  41 1.025974 213.3689
#>      1 113 cpg113  12 1.024725 156.2937
#>      1 236 cpg236  34 1.024028 154.6561
```

The five numbers for each row are the minimum, first quartile, median,
third quartile and maximum across cells. Most cells here have Rhat near
1 and an ESS in the hundreds, so the chains agree and the 400 draws are
worth a few hundred independent ones. The cells with the highest Rhat
are listed below the summary. If many of them sit above an Rhat of 1.05,
or the ESS is low, run longer chains, for example with
`mipca_fit(tuned, obj, iter_sampling = 400)`.

## Extract completed matrices

`as.matrix(fit, draw = b)` returns one completed matrix: the observed
values, with draw `b` filled into the missing values. The rest of this
article follows two CpGs with missing values.

``` r

two <- names(sort(colSums(is.na(obj)), decreasing = TRUE))[1:2]
completed_1 <- as.matrix(fit, draw = 1)
dim(completed_1)
#> [1]  80 300
head(completed_1[, two])
#>         cpg001    cpg005
#> s01 0.21050697 0.9605862
#> s02 0.13928224 0.9314553
#> s03 0.21884599 0.9629385
#> s04 0.04057579 0.8823267
#> s05 0.18268468 0.9643872
#> s06 0.21027450 0.9652667
```

Completed matrices are built one at a time, so the draws are never all
expanded into memory at once. Each one is a full copy of `obj`, so keep
only what the analysis needs. A typical analysis takes a few dozen draws
spread across the chains. The draws are numbered chain by chain, so
evenly spaced numbers take some from every chain. This example takes 20:

``` r

draws <- round(seq(1, 400, length.out = 20)) # print(fit) shows 400 draws
completed <- lapply(draws, function(b) as.matrix(fit, draw = b)[, two])
length(completed)
#> [1] 20
```

Each missing value now has a distribution rather than a single value.
For one missing cell, here are its 20 draws next to its true value in
`truth`:

``` r

miss_row <- which(is.na(obj[, two[1]]))[1]
cell_draws <- vapply(completed, function(m) m[miss_row, two[1]], numeric(1))
quantile(cell_draws, c(0.05, 0.5, 0.95))
#>        5%       50%       95% 
#> 0.1500064 0.2038251 0.3122086
truth[miss_row, two[1]]
#> [1] 0.2257819
```

## Pool an estimate

Run the analysis on every completed matrix and combine the results with
Rubin’s rules. As an example, estimate the mean beta value of the first
of the two CpGs:

``` r

est <- vapply(completed, function(m) mean(m[, two[1]]), numeric(1))
within <- vapply(completed, function(m) var(m[, two[1]]) / nrow(m), numeric(1))

M <- length(completed)
q_bar <- mean(est) # pooled estimate
w <- mean(within) # average within-imputation variance
b <- var(est) # between-imputation variance
total <- w + (1 + 1 / M) * b # total variance

c(estimate = q_bar, se = sqrt(total), se_single = sqrt(w))
#>   estimate         se  se_single 
#> 0.23271605 0.02261688 0.02231940
mean(truth[, two[1]])
#> [1] 0.2357462
```

`se_single` is the standard error you would report from one completed
matrix. The pooled `se` adds the between-imputation variance `b`, which
is the uncertainty about the missing cells that single imputation leaves
out. Here only 16 of 80 values are missing and the model fits well, so
the difference is small. It grows with the share of missing values.

## Columns near 0 or 1: `abb`

Some CpGs are almost always unmethylated or almost always methylated, so
nearly all of their observed values sit at or next to 0 or 1.
[`abb_cols()`](https://hhp94.github.io/slideimp/reference/abb_cols.md)
returns the column numbers of such columns. Passed as `abb`, they are
left out of the model, and their missing values are drawn from their own
observed values by approximate Bayesian bootstrap. This keeps the exact
0s and 1s such columns often hold, but it does not keep their relation
to the other columns.

The simulated matrix has no such column:

``` r

abb_cols(obj)
#> integer(0)
```

Move the first of the two CpGs towards 0 and the second towards 1, and
[`abb_cols()`](https://hhp94.github.io/slideimp/reference/abb_cols.md)
returns both:

``` r

obj_abb <- obj
obj_abb[, two[1]] <- obj_abb[, two[1]] / 100
obj_abb[, two[2]] <- 1 - (1 - obj_abb[, two[2]]) / 100
abb <- abb_cols(obj_abb)
colnames(obj_abb)[abb]
#> [1] "cpg001" "cpg005"
```

Pass them to
[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
as `abb`. It leaves the `abb` columns out of tuning, and
[`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
reads `abb` from the table. Every draw for an `abb` column is one of
that column’s observed values:

``` r

tuned_abb <- mipca_tune(
  obj_abb,
  rank = 1:4,
  transform = "logit_biwhiten",
  seed = 1,
  abb = abb
)
fit_abb <- mipca_fit(tuned_abb, obj_abb, seed = 2, refresh = 0)

miss <- is.na(obj_abb[, two[1]])
drawn <- as.matrix(fit_abb, draw = 1)[miss, two[1]]
all(drawn %in% obj_abb[!miss, two[1]])
#> [1] TRUE
```

## Impute target CpGs on flanking windows

When only some CpGs need draws, such as the CpGs of an epigenetic clock,
each one can be fitted on the CpGs near it instead of on the whole
matrix.
[`mipca_flank_plan()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
finds the window of each target: the columns within `window_size` of it
on `location`, on both sides. A target whose window has fewer than
`min_window_n` columns is not imputed.
[`mipca_flank()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md)
then fits
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) on each
window and keeps the target’s draws.

The simulated CpGs have no genomic positions, so the column order stands
in for them here. Use real coordinates, one chromosome at a time, on
real data.

``` r

location <- seq_len(ncol(obj))
plan <- mipca_flank_plan(
  obj,
  location,
  window_size = 20,
  subset = two,
  min_window_n = 10
)
plan
#> si_flank_plan: 2 targets on a 80 x 300 matrix, nothing fitted
#>   status: impute 2
#>  target   name n_miss start end window_n status
#>       1 cpg001     16     1  21       21 impute
#>       5 cpg005     16     1  25       25 impute
#>   next: mipca_flank(x, rank, seed, transform) imputes.
#>         mipca_flank_window(x, subset) gives one target's window to tune on.
```

To choose the rank and `nu` for the windows, tune on one target’s
window. A window has few columns, so lower `n_cols`, the number of
columns that get hidden cells. Its default is 10% of the columns, but at
least 50:

``` r

win <- mipca_flank_window(plan, two[1])
tuned_win <- mipca_tune(
  win,
  rank = 1:3,
  transform = "logit_biwhiten",
  n_cols = 10,
  seed = 1
)
tuned_win
#> si_tune: 3 ranks on a 80 x 21 matrix, holdout 40 cells over 10 columns
#>  rank    crps ess_p10 ess_ok rhat_max implied_nu implied_lo implied_hi  time
#>     1 0.04828   278.1  0.095    1.030      3.045      2.620      10.98 0.047
#>     2 0.03778   233.7  0.040    1.045        Inf      8.555        Inf 0.054
#>     3 0.03835   193.1  0.015    1.054        Inf        Inf        Inf 0.061
#>   next: mipca_fit(x, obj) fits the rank with the lowest crps.
#>         mipca_fit(x, obj, rank) fits another rank.
```

Then give the rank and `nu` to
[`mipca_flank()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md).
It warns about a `nu` below 5 but does not raise it, so raise it to 5
here, as
[`mipca_fit()`](https://hhp94.github.io/slideimp/reference/mipca_fit.md)
would:

``` r

best <- which.min(tuned_win$crps)
fit_flank <- mipca_flank(
  plan,
  rank = tuned_win$rank[best],
  nu = max(tuned_win$implied_nu[best], 5),
  seed = 1,
  transform = "logit_biwhiten",
  .progress = FALSE
)
fit_flank
#> si_mi fit: 80 x 2, 32 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
dim(as.matrix(fit_flank, draw = 1))
#> [1] 80  2
```

The fit holds the targets only, and
[`as.matrix()`](https://rdrr.io/r/base/matrix.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and the pooling above
work on it as on a fit from
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md).
[`cbind()`](https://rdrr.io/r/base/cbind.html) joins the fits of several
matrices, such as one per chromosome.

## References

Audigier, V., Husson, F. and Josse, J. (2016). Multiple imputation for
continuous variables using a Bayesian principal component analysis.
*Journal of Statistical Computation and Simulation*, 86(11), 2140-2156.
<https://doi.org/10.1080/00949655.2015.1104683>

Gelman, A. and Rubin, D. B. (1992). Inference from iterative simulation
using multiple sequences. *Statistical Science*, 7(4), 457-472.
<https://doi.org/10.1214/ss/1177011136>

Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical
Science*, 7(4), 473-483. <https://doi.org/10.1214/ss/1177011137>

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Burkner, P.-C.
(2021). Rank-normalization, folding, and localization: an improved R-hat
for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.
<https://doi.org/10.1214/20-BA1221>
