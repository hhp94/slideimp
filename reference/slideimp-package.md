# slideimp: Numeric Matrices K-NN and PCA Imputation

Provides K-nearest-neighbor, PCA, grouped, and sliding-window imputation
methods for numeric matrices, along with utilities for simulation and
parameter tuning.

## Column mean fill

A column with a proportion of missing values greater than `colmax` takes
no part in the imputation. The method neither uses the column to impute
other columns nor imputes it.

With `post_imp = TRUE`, column means replace the missing values that
remain after the method has run. Two kinds of missing value can remain:
the missing values of the columns that took no part, and the cells for
which the method found no value. Each column mean is the mean of the
column's observed values. With `post_imp = FALSE`, both kinds stay
missing.

The result records what the column means replaced, in two attributes:

- `mean_filled_cols`: an integer vector with the index, in the result,
  of each column that took no part and was filled.

- `mean_filled_cells`: an integer matrix with the columns `row` and
  `col`, one row for each cell that the method could not impute and that
  was filled.

Both attributes are empty when no column mean was used.
[`print()`](https://rdrr.io/r/base/print.html) reports them when they
are not empty.

## Missing values and non-finite input

`obj` must be a double matrix. An integer or logical matrix is rejected,
not converted. Convert one with `storage.mode(obj) <- "double"`.

`NA` and `NaN` are both treated as missing. The package does not
distinguish between them.

`Inf` and `-Inf` are rejected on entry, and so is a column with no
observed values (every cell `NA` or `NaN`). Locate them with
`which(is.infinite(obj), arr.ind = TRUE)` and
`which(colSums(!is.na(obj)) == 0)`.

## See also

Useful links:

- <https://github.com/hhp94/slideimp>

- <https://hhp94.github.io/slideimp/>

- Report bugs at <https://github.com/hhp94/slideimp/issues>

## Author

**Maintainer**: Hung Pham <amser.hoanghung@gmail.com>
([ORCID](https://orcid.org/0000-0002-8271-9355)) \[copyright holder\]

Authors:

- Hung Pham <amser.hoanghung@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-8271-9355)) \[copyright holder\]

Other contributors:

- Francois Husson ([ORCID](https://orcid.org/0000-0002-7271-8877))
  ('missMDA' imputePCA() and MIPCA(), reimplemented in R/pca_imp.R,
  src/pca_kernel.cpp, R/mi_pca_imp.R and src/mipca.cpp) \[contributor,
  copyright holder\]

- Julie Josse ([ORCID](https://orcid.org/0000-0001-9547-891X))
  ('missMDA' imputePCA() and MIPCA(), reimplemented in R/pca_imp.R,
  src/pca_kernel.cpp, R/mi_pca_imp.R and src/mipca.cpp) \[contributor,
  copyright holder\]

- Vincent Audigier ('missMDA' MIPCA() Bayesian sampler, reimplemented in
  R/mi_pca_imp.R and src/mipca.cpp) \[contributor, copyright holder\]

- Posit Software, PBC ('carrier' package code, adapted in R/crate.R and
  tests/testthat/test-crate\*.R) \[copyright holder\]

- Yixuan Qiu (Spectra library, vendored in src/Spectra) \[contributor,
  copyright holder\]

- Gael Guennebaud (Eigen SelfAdjointEigenSolver code, adapted in
  src/Spectra/LinAlg/TridiagEigen.h) \[contributor, copyright holder\]

- Jitse Niesen (Eigen SelfAdjointEigenSolver code, adapted in
  src/Spectra/LinAlg/TridiagEigen.h) \[contributor, copyright holder\]
