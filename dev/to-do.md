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
- **Breaking:** `col_vars()` loses its `cores` argument and runs on one
  thread; a call that passes `cores` is now an error.
- **Breaking:** `mean_imp_col()` is no longer exported. The column means it
  sets, for `post_imp` and for `on_infeasible = "mean"` in every imputer, are
  now summed and divided in long double, so they can differ from 1.2.0's in
  the last bits.
- `pca_imp()` gains `output = c("matrix", "mipca")`. `"mipca"` returns the
  point estimate as a one-draw fit for `mipca_complete()` and `cbind()`, and
  is a new error on an input without unique row and column names. For an
  integer input with no hole, `as.matrix()` of that fit is double where the
  matrix output stays integer.
- `pca_imp()` gains `subset`: the result, in both forms, holds only the
  selected columns, in input order, while every eligible column still takes
  part in the fit. An unknown name or an empty `subset` is an error there,
  where `knn_imp()`'s `subset` still drops the one with a message and returns
  early on the other. The "Column mean fill" section now says
  `mean_filled_cols` counts the result's columns.
- `pca_imp()` no longer writes past a buffer when asked for more components
  than the data has rank (1.2.0 does, in its dsyevr call:
  `dev/eigen_conversion.md`, "The review of the wrapper layer").
- `pca_imp()` under LOBPCG with 17 or more eigenpairs (`ncp` of 16 or more
  when regularized, 17 or more with `method = "EM"`) returns values that
  differ from 1.2.0's by about 1e-8, inside what the fit's own `threshold`
  resolves: dsyevd is given the workspace length of each problem, where
  1.2.0 kept the longest so far (`dev/eigen_conversion.md`, Phase 2,
  step 4). The `n_exact`, `n_lobpcg_ok` and `n_lobpcg_bad` attributes of
  such a result can differ too.
- The `n_iter` attribute of a `pca_imp()` result is an integer; it was a
  double holding a whole number. The same for `info$n_iter` and `prov$iter`
  in the chunks of its `mipca` output (`dev/eigen_conversion.md`, Phase 2,
  step 6).
- `pca_imp()` with `ncp` above 32 returns values that differ from 1.2.0's:
  in the last digits under the exact solver (about 1e-14), by about 1e-8
  under LOBPCG, where the solver counts can differ too. The sum of the
  squared leading singular values is taken by one loop at every length,
  where 1.2.0 gave more than 32 of them to ddot (`dev/eigen_conversion.md`,
  Phase 2, step 6).

## Parts for the other imputers

Decided with the maintainer, 2026-11-17, after `pca_imp()` moved onto parts
(`dev/si_mi_conversion.md`). The parts code is in `R/parts.R`, which holds
`new_part()`, `fill_parts()` for the matrix output, and `parts_to_si_mi()`
and `check_si_mi_names()` for the `"mipca"` output.

- `knn_imp()` goes through the same change: it builds parts, and it gains
  `output = c("matrix", "mipca")`, with the same name rule and the same
  record on every chunk as `pca_imp()`.
- `group_imp()` builds its result from its groups' parts, to save memory:
  no group's full result matrix, and one write into one copy of the input.
- `tune_imp()` and `slide_imp()` keep returning a matrix. Until
  `group_imp()` takes parts, it passes `output = "matrix"` to `pca_imp()`, as
  `slide_imp()` and `tune_imp()` do.

Added 2026-11-17: `subset` narrows the output, package-wide, in input order,
with no flag (`dev/si_mi_conversion.md`, "What was decided"; the doors and
`pca_imp()` do it there, Phases 5 and 6). What is left, delayed by the
maintainer the same day:

- Planned in `dev/si_mi_conversion.md`, "Phase 7", with what was decided
  on 2026-11-18 (the leftover cells' mean, the kernel filling them) and
  what is still to settle before the code.
- `knn_imp()` lands in ONE change (maintainer, 2026-11-17): it builds parts,
  gains `output`, and its `subset` narrows its output through the doors'
  `keep`. Its `subset` resolves as `pca_imp()`'s does, through
  `resolve_narrow_subset()` in `R/utils.R`: an
  unknown name and an empty `subset` are refused, where `resolve_subset()`
  drops the one with a message and returns early on the other. The
  narrowing breaks every caller that passes `subset` today, and the refusals
  may break `slide_imp()`; both are fixed in the same change:
  - `tune_imp()` reads each repetition's estimates at `na_loc` in the full
    matrix's numbering, in both its paths (`R/tune_imp.R`:
    `estimate = imputed_result[na_positions]`), and must rebase them to the
    narrowed columns.
  - `slide_imp()` and `group_imp()` index each window's or group's result
    by its local subset columns, and must read the narrowed result instead.
    Both drop a workaround for `pca_imp()` having no `subset`, by passing it
    instead: `slide_imp()`'s complete-targets check (`R/slide_imp.R`:
    `targets_complete <- imp_method == "pca" &&`) and `group_imp()`'s
    (`R/group_imp.R`: `needs_imp <- vapply(`).
  - `slide_imp(subset = )` and `group_imp(subset = )` narrow their own
    output too, by the same rule.
  - The NEWS entry above gains a **Breaking** line for all of it.
  - `parts_to_si_mi()`'s `prov_extra` is a list beside `parts`, matched to
    them by position and checked only by length, so a part list built in
    another order puts `rank` and `iter` on the wrong chunk without an error.
    With `knn_imp()` a second builder arrives: carry those fields on the
    part instead (`dev/si_mi_conversion.md`, "Review of Phases 1 to 6").

## Kernel cleanup after the Eigen conversion

Agreed 2026-11-18, to run once `dev/eigen_conversion.md` has landed, each as
its own change. Vendored files (`dev/vendored.md`) are out of scope.

- **Investigate dropping `col_min_max()`.** Its one caller is `sim_mat()`,
  rescaling `rnorm()` output to [0, 1] when `beta = TRUE`. Measure whether
  R's own column minima and maxima (`matrixStats` is not a dependency, so
  base R or collapse) cost anything that matters at `sim_mat()`'s sizes,
  and whether the output stays the same under a fixed seed; drop the
  kernel if neither stands in the way.
- **Comment trim pass over `src/`.** Comments that narrate history, repeat
  the code, or explain the Armadillo version a file no longer uses go;
  what stays says why the code is as it is, true for the code as written.
- **Kernel hygiene pass over `src/`.** Hoist loop-invariant work, and move
  checks, allocations and branches that do not depend on the loop index out
  of hot paths. Measure each change against the reference capture and the
  timing driver; no hand unrolling or OpenMP (`dev/CPP.md`). Candidates
  seen in the K-NN kernel while converting it, none measured yet:
  `rows_to_impute` repeats, 0-based, the rows `loc` already holds;
  `obj_masked` and `nmiss_masked` are zero-filled by `std::vector` before
  every cell is written; and `initialize_result()` guards invariants its
  one caller has just established. Seen in the PCA kernel
  (`src/armaSVD.cpp`), none measured apart: `Xhat`, `X_work`, `AA_NxN`, `U`,
  `V`, `Vn` and the eigenpair buffers are zero-filled by `std::vector`
  before every cell that is read is written, where only `fittedX` needs its
  zeros; and `row_w` is copied before it is checked.

## Release condition: pass on the RcppEigen candidate

slideimp 1.2.0 on CRAN does not link to RcppEigen. The sampler port added it
to `LinkingTo`, so the next release is the first that does, whether or not
the Armadillo removal (`dev/eigen_conversion.md`) is in it. RcppEigen is
moving from Eigen 3.4 to Eigen 5, and its maintainer is patching by hand the
reverse dependencies that fail on the candidate
(https://github.com/RcppCore/RcppEigen/issues/151). A release that passes
only on CRAN's RcppEigen would add slideimp to that list.

Before any release that links to RcppEigen, on both CRAN's RcppEigen and the
candidate current on that day (0.3.4.0.2 and 0.4.9.9-2 when this was written,
2026-11-18):

- The whole suite, `check_manual()` and `devtools::check()`. The candidate
  goes in its own private library, selected with `R_LIBS` before R starts;
  `dev/eigen_conversion.md`, "Both RcppEigen versions", says how. Record
  both versions tested in the release commit.
- The sampler's draws under fixed seeds, compared between the two builds:
  worst absolute and worst relative difference, and whether they are
  bitwise identical. If they move beyond what the `test-mipca-*` files
  tolerate, that is a finding for upstream (`hhp94/test`), since the sampler
  is vendored.

State when this was written: the maintainer has compiled the tree against
the candidate. The suite, the check and the draw comparison against it have
not been run.

Nothing is vendored from Eigen and no version is pinned; that was considered
and ruled out (maintainer, 2026-11-18).

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

- **The next copy.** The copied files match `hhp94/test@654aadc`
  (upstream's `vendor_diff.R`, every file MATCH, Spectra included), the sha
  to record upstream as the last port. Take the next copy from its latest
  pushed commit, with no change made by hand, and name that sha in the
  commit body; when it runs is the maintainer's call. Until help pages are
  written into the copied files a whole copy loses nothing. After that,
  use the patch route in upstream's brief (its vendor-to-do, "The next
  port, step by step"), which carries changes since the recorded sha
  without overwriting the roxygen.
- **Makevars: drop OpenMP.** Ruled upstream 2026-11-14: `src/Makevars` and
  `src/Makevars.win` drop `$(SHLIB_OPENMP_CXXFLAGS)` from `PKG_CXXFLAGS` and
  `PKG_LIBS`, since the sampler no longer uses it. This also turns off
  Armadillo's OpenMP in slideimp's own PCA and K-NN code, a cost not yet
  measured here. Measure it on those paths before making the change, as its
  own commit.
- **`USE_FC_LEN_T` before the R headers.** `src/mipca.cpp` and
  `src/mipca_crps.cpp` call BLAS through `R_ext/BLAS.h` and pass `FCONE`,
  but neither defines `USE_FC_LEN_T` before its first R header, and neither
  falls back with `#ifndef FCONE`, as *Writing R Extensions*, "Fortran
  character strings", writes the call (`dev/r_manual/R-exts.md`). The
  manual notes that `USE_FC_LEN_T` is the default from R 4.3.0, the oldest R
  the package supports, so the lengths are passed today; that the calls are
  right on R 4.3.0 itself is not measured here. Ruled 2026-11-18: a file
  that disagrees with the manual is a bug to fix, not a model.
  `src/si_lapack.cpp` has the manual's form. Upstream.
- **The Spectra licence notice.** `src/Spectra/` is vendored, and its
  copyright holders are in `Authors@R`, but no notice file ships its
  licence: `inst/LICENSE.note` covers only the code adapted from carrier.
  Upstream's `inst/COPYRIGHTS` is the text; which file carries it here,
  that note or a new one, is slideimp's choice.
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
  target in a plan over fewer columns draws differently. Do not tell users
  to pin BLAS threads.
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
