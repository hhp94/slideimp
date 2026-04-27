# LOBPCG Eigensolver Control Options

Construct a validated list of control options for the LOBPCG eigensolver
used by
[`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).
Most users do not need to call this directly.

## Usage

``` r
lobpcg_control(warmup_iters = 10L, tol = 1e-09, maxiter = 20)
```

## Arguments

- warmup_iters:

  Integer. Number of warm-up iterations before the main LOBPCG solve.
  Must be non-negative.

- tol:

  Numeric. Convergence tolerance for the LOBPCG eigensolver. Must be
  non-negative and finite.

- maxiter:

  Integer. Maximum number of LOBPCG iterations. Must be non-negative.
  Setting `maxiter = 0` disables LOBPCG and uses the exact solver
  instead.

## Value

A named list of class `"slideimp_lobpcg_control"` containing
`warmup_iters`, `tol`, and `maxiter`.

## Examples

``` r
# Use all defaults
lobpcg_control()
#> $warmup_iters
#> [1] 10
#> 
#> $tol
#> [1] 1e-09
#> 
#> $maxiter
#> [1] 20
#> 
#> attr(,"class")
#> [1] "slideimp_lobpcg_control"

# Override a single option
lobpcg_control(maxiter = 50)
#> $warmup_iters
#> [1] 10
#> 
#> $tol
#> [1] 1e-09
#> 
#> $maxiter
#> [1] 50
#> 
#> attr(,"class")
#> [1] "slideimp_lobpcg_control"

# Disable LOBPCG and use the exact solver
lobpcg_control(maxiter = 0)
#> $warmup_iters
#> [1] 10
#> 
#> $tol
#> [1] 1e-09
#> 
#> $maxiter
#> [1] 0
#> 
#> attr(,"class")
#> [1] "slideimp_lobpcg_control"

# Pass directly to pca_imp()
set.seed(123)
obj <- sim_mat(10, 10)$input
pca_imp(obj, ncp = 2, lobpcg_control = lobpcg_control(tol = 1e-9))
#> Method: PCA imputation
#> Dimensions: 10 x 10
#> 
#>          feature1   feature2  feature3  feature4  feature5   feature6
#> sample1 0.5784798 0.06296727 0.3155309 0.1199980 0.1807391 0.31919912
#> sample2 0.4991812 0.44077231 0.2120510 0.3257524 0.1919521 0.14843947
#> sample3 0.7709271 0.75477764 1.0000000 0.5099311 0.6027900 0.75450992
#> sample4 0.5068375 0.37347042 0.6018860 1.0000000 0.5850265 0.08171420
#> sample5 0.4165827 0.42553421 0.6024750 0.7728095 0.2295133 0.08343629
#> sample6 1.0000000 0.93942190 0.9867413 0.5851340 1.0000000 1.00000000
#> # Showing 6 of 10 rows and 6 of 10 columns

# Or use a named list
pca_imp(obj, ncp = 2, lobpcg_control = list(maxiter = 0))
#> Method: PCA imputation
#> Dimensions: 10 x 10
#> 
#>          feature1   feature2  feature3  feature4  feature5   feature6
#> sample1 0.5784798 0.06296727 0.3155309 0.1199980 0.1807391 0.31919912
#> sample2 0.4991812 0.44077231 0.2120510 0.3257524 0.1919521 0.14843947
#> sample3 0.7709271 0.75477764 1.0000000 0.5099311 0.6027900 0.75450992
#> sample4 0.5068375 0.37347042 0.6018860 1.0000000 0.5850265 0.08171420
#> sample5 0.4165827 0.42553421 0.6024750 0.7728095 0.2295133 0.08343629
#> sample6 1.0000000 0.93942190 0.9867413 0.5851340 1.0000000 1.00000000
#> # Showing 6 of 10 rows and 6 of 10 columns
```
