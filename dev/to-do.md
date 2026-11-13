# To do

Work that is agreed on but not yet done. An entry stays here until the change
lands; then it is deleted, not marked done - git history is the record of what
happened, this file is the record of what has not.

An entry says what to change and why, enough that a session with no memory of
the conversation can pick it up. If an entry turns out to be wrong, delete it
and say so in chat.

## Draft `NEWS.md` for the release after 1.2.0

`NEWS.md` is written once per release, from the whole diff since the last
release, not one entry per change (see `CLAUDE.md`). Draft the section from
`git diff v1.2.0..HEAD` and `git log v1.2.0..HEAD`, with the maintainer, who
owns its style.

The development-version section already in `NEWS.md` was written one change at
a time, before this rule. Treat it as material for the draft, not as finished
text.

Easy to miss in the diff, noted so the draft does not:

- **Breaking:** the `na_check` argument of `knn_imp()`, `pca_imp()`,
  `slide_imp()` and `group_imp()` is removed, along with the
  `has_remaining_na` attribute it controlled and the note `print()` showed
  from it. `anyNA()` on the result answers the same question.
- Results of `knn_imp()`, `pca_imp()`, `slide_imp()` and `group_imp()` carry
  two new attributes, `mean_filled_cols` and `mean_filled_cells`, that record
  what `post_imp` set to a column mean, and `print()` shows a note for each
  when it is not empty. The help pages gain a "Column mean fill" section and
  describe `colmax` as deciding which columns take part, neighbors included.
- The `post_imp` attribute of a `group_imp()` result now holds the value the
  groups applied: `TRUE` by default, where it was `NULL`, and `NA` when
  `group$parameters` sets it differently between groups or when no group was
  imputed and no value was passed.
- `pca_imp()` on a matrix with no missing values now returns a
  `slideimp_results` with the same attributes as any other call, where it
  returned the bare input. `print()` shows it with the result header, and the
  call now copies the input once to attach the class.
- `sample_na_loc()`, and so `tune_imp()`, no longer draws a column whose
  variance is above zero but below `.Machine$double.eps`, the threshold
  `pca_imp()` uses for eligibility. Such a column was drawn before, and a
  `pca_imp()` repetition that blanked only such columns failed.

## Sampler port: what is still owed

The multiple-imputation sampler (`R/mi_pca_imp.R`, `R/mi_pca_s3.R`,
`R/mi_pca_tune.R`, `R/mi_pca_flank.R`, `R/mi_pca_abb.R`, `src/mipca.cpp`,
`src/mipca_crps.cpp`, `src/mipca_diag.cpp`) is developed in a scratch package
in another repository, `hhp94/test`, and copied here as one unit. The copied
files are listed in `dev/vendored.md`, and the house rule on vendored files
in `CLAUDE.md` applies to them: apart from their roxygen blocks, an edit to
their code or comments is made in `hhp94/test` and arrives with the next
copy. The items below marked "Upstream" are made there. The scratch package
keeps its copies after every port and its gate scripts keep testing them;
the two copies stay in step through `dev/vendored.md` and the export.

Ruled, so not owed (maintainer, 2026-11-13). Do not add these back as items:

- The sampler's refusals stay plain `stop(..., call. = FALSE)`, and its
  messages and warnings stay base `message()` and `warning()`. No cli
  conversion, in the copied files or in the `test-mipca-*` files that pin
  their wording.
- No bytes warning before allocating, and no refusal at `2^31` cells. The
  draws refuse in plain words when R cannot allocate them, naming the size
  in GB and what makes it up. Base R handles the long vector a `2^31`
  refusal was meant to guard.
- `subset` keeps its `NULL` default at `mipca()`, meaning the feature is
  off: draws for every column with a hole. In flank mode it is required,
  and `mipca_flank_plan()` and `mipca_flank_window()` already refuse a
  missing or `NULL` `subset`.
- Every comment in the copied files comes across with the code, written to
  be true inside slideimp.

Owed:

- **The next copy.** Take the export of `hhp94/test` at its latest pushed
  commit, with no change made by hand, and name that commit's sha in the
  commit body. When it runs is the maintainer's call. The first copy
  (`e7f2f0c`, from `hhp94/test@b4cfbf1`) made four changes by hand, and
  `hhp94/test` has since taken all four, so none is needed again: no
  `transform` default at `mipca()`, `mipca_tune()` and `mipca_flank()`;
  `tf_asin()` and `.chk_no_ndiag()` removed; the debug-build line of
  `print.si_tune()` removed. `dqrng (>= 0.4.1)` stays in `LinkingTo`. The
  next copy also brings:
  - the comments;
  - the printing pass: each print shows a header, a line only where
    something waits or was left `NA`, and the next calls;
  - the size refusal;
  - `mipca_tune()` no longer refuses a matrix with a constant column. Its
    holdout draw leaves such a column out, as `sample_na_loc()` does since
    `7a863d0`, so the copy needs `sample_na_loc()` at `7a863d0` or later.
- **The vignette.** `vignettes/articles/multiple-imputation-beta.Rmd` was
  written for the October surface. Its calls were repaired so that it knits
  on the current one, and nothing else was changed. Whether it stays is not
  decided. If it stays, it does not yet cover `abb`, `mipca_complete()` or
  `mipca_flank()`.
- **Roxygen.** Every export of the sampler is a bare `#' @export` and has no
  help page. Write the help pages here, by the rules of `dev/WRITING.md`. The
  pages must say: `subset` is a memory parameter, not a compute one, because
  every missing cell is imputed at every iteration; at `mipca()` `NULL`
  means every column with a hole, and the flank functions require it. A
  `mipca_flank()` result can hold `NA` and the plan says why. A flank window
  is seeded by its target's column index in the plan's matrix, so the same
  target in a plan over fewer columns draws differently.
  `tf_biwhiten(rows = FALSE)`, passed as an object, is the way out when the
  row scale refuses a window. Do not tell users to pin BLAS threads.
- **Columns the core refuses.** The core refuses a column that holds `Inf`,
  is all `NA`, or has one observed value. Real matrices contain these, so a
  wrapper outside the core must carry them round the fit, as `pca_imp()` does
  with its eligible columns. It also owns the near-zero-variance threshold.
  Upstream if the wrapper, or any change to the core's refusal, goes into the
  copied files. Only a wrapper in a file of slideimp's own is written here.
- **`mipca_core()` takes integer `subset` only.** Name resolution moves to
  the callers above it, now that `resolve_subset()` is a sibling internal.
  Upstream.
- **The `Inf` refusal has two shapes.** slideimp reports `Inf` by row and
  column from `check_finite()`. The sampler reports it per column from
  `col_vars()`. Decide whether one message shape wins. A change to the
  sampler's shape is made upstream.
- **Tests not yet written.** The scratch package checks the tuner, the flank
  driver and the ABB path with gate scripts, not testthat files, and the gates
  stay there. No test here runs `mipca_tune()`, `mipca_fit()`,
  `mipca_flank()`, `mipca_complete()` or `abb_cols()`. A test that lands as
  a `test-mipca-*` file is part of the copy, so it is written upstream.
- **Not for the first version.** Eigen's own products in place of the
  sampler's direct BLAS calls, so that parallel chains pay on every BLAS.
  Draws stop being bitwise the BLAS run, so every stored baseline is retaken.
  Upstream.
