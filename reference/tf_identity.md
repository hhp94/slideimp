# Transforms for Bayesian PCA Imputation

Transform `obj` prior to the Bayesian PCA imputation. Each draw is
mapped back to the scale of `obj` before it is returned.

## Usage

``` r
tf_identity()

tf_sv(n)

tf_biwhiten(base = tf_identity(), rows = TRUE)
```

## Arguments

- n:

  Integer. The sample size of the squeeze, usually `nrow(obj)`. At least
  `2`.

- base:

  A transform object. The scale the rescaling is applied on:
  `tf_identity()`, `tf_sv()` or a transform of your own, not another
  `tf_biwhiten()`. Defaults to `tf_identity()`.

- rows:

  Logical. If `TRUE`, rescale the rows as well as the columns. Defaults
  to `TRUE`.

## Value

A transform object, to pass as `transform`.

## Details

**The two names.** `transform = "identity"` is `tf_identity()`, and
`transform = "logit_biwhiten"` is
`tf_biwhiten(tf_sv(nrow(obj)), rows = TRUE)`. Build a transform with
these functions to change one part of either, such as the logit without
the rescaling, or the columns rescaled and the rows not.

**`tf_identity()`** fits the values as is.

**`tf_sv()`** is for proportions such as DNAm beta values. Each value
`y` is squeezed to `(y * (n - 1) + 0.5) / n`, which moves `0` and `1`
inside the interval (Smithson and Verkuilen, 2006), and then taken to
the logit scale. Draws are mapped back and clamped to `[0, 1]`. Every
observed value must lie in `[0, 1]`.

**`tf_biwhiten()`** rescales the columns, and with `rows = TRUE` the
rows, of the matrix on the scale of `base` to a common noise level
before the fit, and takes the scales back out of each draw. The scales
are estimated once from the observed values, with the column variances
shrunk as in Smyth (2004). The rescaling follows the idea of Landa,
Zhang and Kluger (2022), with the scales estimated differently.

**A transform of your own** is a list with `name`, a string, and two
functions, `forward(y, ...)` and `backward(z, ...)`, where `backward()`
undoes `forward()` value by value. Both must take `...`, since the fit
passes them column and row positions. It can be the `base` of
`tf_biwhiten()`.

## References

Smithson M, Verkuilen J (2006). A better lemon squeezer?
Maximum-likelihood regression with beta-distributed dependent variables.
*Psychological Methods*, 11(1), 54-71.
[doi:10.1037/1082-989X.11.1.54](https://doi.org/10.1037/1082-989X.11.1.54)

Smyth GK (2004). Linear models and empirical Bayes methods for assessing
differential expression in microarray experiments. *Statistical
Applications in Genetics and Molecular Biology*, 3(1), 1-25.
[doi:10.2202/1544-6115.1027](https://doi.org/10.2202/1544-6115.1027)

Landa B, Zhang TTCK, Kluger Y (2022). Biwhitening reveals the rank of a
count matrix. *SIAM Journal on Mathematics of Data Science*, 4(4),
1420-1446. [doi:10.1137/21M1456807](https://doi.org/10.1137/21M1456807)

## See also

[`mipca()`](https://hhp94.github.io/slideimp/reference/mipca.md), whose
`transform` argument these build.

## Examples

``` r
set.seed(1)
sim <- sim_mat(n = 50, p = 20, beta = TRUE)
obj <- sim$input

# The same model as transform = "logit_biwhiten"
tf_logit_biwhiten <- tf_biwhiten(tf_sv(nrow(obj)), rows = TRUE)

# The logit without the rescaling
tf_logit <- tf_sv(nrow(obj))

# The columns rescaled, the rows not
tf_cols <- tf_biwhiten(tf_sv(nrow(obj)), rows = FALSE)
fit <- mipca(obj, rank = 2, transform = tf_cols, seed = 1)
fit
#> si_mi fit: 50 x 20, 100 imputed cells, 400 draws
#>   next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS
```
