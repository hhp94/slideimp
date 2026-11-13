# Vendored files

The multiple-imputation sampler is developed in another repository,
`hhp94/test`, and copied into slideimp as one unit. These files are that
copy:

- `R/mi_pca_imp.R`
- `R/mi_pca_s3.R`
- `R/mi_pca_tune.R`
- `R/mi_pca_flank.R`
- `R/mi_pca_abb.R`
- `src/mipca.cpp`
- `src/mipca_crps.cpp`
- `src/mipca_diag.cpp`
- `tests/testthat/test-mipca-admissible.R`
- `tests/testthat/test-mipca-door-defaults.R`
- `tests/testthat/test-mipca-gram-bound.R`
- `tests/testthat/test-mipca-lifecycle.R`
- `tests/testthat/test-mipca-loop-counts.R`
- `tests/testthat/test-mipca-narrow-join.R`
- `tests/testthat/test-mipca-one-cell.R`
- `tests/testthat/test-mipca-rebase-transform.R`

Each copy lands here as one port commit, and that commit's body names the
source commit as `hhp94/test@<sha>`. Find it with
`git log --grep 'hhp94/test@'`.

A change to these files' code or comments is made in `hhp94/test` first and
reaches slideimp with the next copy. Only their roxygen blocks are written
here. The house rule in `CLAUDE.md` says what that means for a pass over the
package.

A file added to or removed from the copy is added to or removed from this
list in the port commit that brings it.
