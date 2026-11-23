# Propose Columns to Draw by Approximate Bayesian Bootstrap

Find the columns of a matrix of proportions, such as DNAm beta values,
whose observed values lie almost all near `0` or almost all near `1`.
Pass the result as `abb` to
[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md),
[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
or
[`mipca_flank_plan()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md).

## Usage

``` r
abb_cols(obj, bound = 0.02, q = 0.95)
```

## Arguments

- obj:

  A double matrix with samples in rows and features in columns. Row and
  column names are required and must be unique. Every observed value
  must lie in `[0, 1]`.

- bound:

  Numeric. A value below `bound` counts as near `0`, and a value above
  `1 - bound` as near `1`. Must be above `0` and below `0.5`. Defaults
  to `0.02`.

- q:

  Numeric. Proportion of a column's observed values that must be near
  `0`, or near `1`, for the column to be listed. Must be above `0.5` and
  at most `1`. Defaults to `0.95`.

## Value

An integer vector of column numbers in `obj`, in increasing order. It is
empty when no column is listed, and `abb` takes it as no list.

## Details

**The listed columns.** A column listed in `abb` is left out of the
model, and its missing values are drawn from its own observed values by
approximate Bayesian bootstrap (Rubin and Schenker 1986). Unlike the
model, this keeps the exact `0`s and `1`s such columns often hold. Each
column is drawn on its own, so its relation to the other columns is not
kept.
[`mipca_tune()`](https://hhp94.github.io/slideimp/reference/mipca_tune.md)
tunes on the other columns, and
[`summary()`](https://hhp94.github.io/slideimp/reference/summary.si_mi.md)
does not report the listed ones.

**Choosing the list.** The defaults list only the most extreme columns.
Raise `bound` or lower `q` to list more, or pass any columns to `abb`. A
column with no observed values is never listed.

## References

Rubin DB, Schenker N (1986). Multiple imputation for interval estimation
from simple random samples with ignorable nonresponse. *Journal of the
American Statistical Association*, 81(394), 366-374.
[doi:10.1080/01621459.1986.10478280](https://doi.org/10.1080/01621459.1986.10478280)

## See also

- [`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md) to
  fit with the list.

- [`mipca_flank_plan()`](https://hhp94.github.io/slideimp/reference/mipca_flank.md),
  which leaves the listed columns out of every window.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)
obj <- sim$input

# Push two columns against 0 and 1
obj[, 3] <- obj[, 3] / 100
obj[, 4] <- 1 - (1 - obj[, 4]) / 100

abb <- abb_cols(obj)
abb
#> [1] 3 4

fit <- mipca(obj, rank = 2, transform = "logit_biwhiten", seed = 1, abb = abb)
fit
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
```
