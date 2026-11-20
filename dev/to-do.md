# To do

Work that is agreed on but not yet done. An entry stays here until the change
lands; then it is deleted, not marked done - git history is the record of what
happened, this file is the record of what has not.

An entry says what to change and why, enough that a session with no memory of
the conversation can pick it up. If an entry turns out to be wrong, delete it
and say so in chat.

`dev/eigen_conversion.md`, the plan that took Armadillo out of `src/`, was
deleted when its last phase landed. Entries below still cite it by section;
read it with `git show 4f2b669:dev/eigen_conversion.md`.

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
- `pca_imp()` no longer ends the R session when asked for more components
  than the completed data has rank, as on replicate samples. 1.2.0 writes
  past a buffer in its dsyevr call there, and that is a crash, not a silent
  write: built from source on Windows it segfaults on `pca_imp(x, ncp = 5)`,
  every other argument at its default, for 30 samples of which 6 are
  distinct (`dev/review_p2.md`, finding 1; no other platform measured).
  `tests/testthat/test-pca_imp.R` has the regression test.
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
  in the last digits under the exact solver, by about 1e-8 under LOBPCG,
  where the solver counts can differ too. Neither is a bound, since the
  difference grows with the iteration count: under the exact solver, at
  most 8.2e-15 absolute on the plan's cases, and 1.6e-13 absolute, 1.2e-11
  relative, on a fit of 202 iterations (`dev/review_p2.md`, finding 6).
  The sum of the squared leading singular values is taken by one loop at
  every length, where 1.2.0 gave more than 32 of them to ddot
  (`dev/eigen_conversion.md`, Phase 2, step 6).
- `pca_imp()` with `solver = "lobpcg"` on a Gram of order 2 to 4, when one
  eigenpair or every eigenpair of it is asked for, returns values that
  differ from 1.2.0's: by 5.2e-9 absolute at most over the 267 fits measured
  to differ. `n_iter` is the same in each; `n_exact`, `n_lobpcg_ok` and
  `n_lobpcg_bad` differ in most fits that ask for every eigenpair. A product
  of that size goes to BLAS, where Armadillo ran loops of its own
  (`dev/eigen_conversion.md`, Phase 2, step 4, "What changes and is
  reachable"). The default solver takes the exact path at that size and
  returns what it did; `group_imp()` passing `solver = "lobpcg"` to groups
  of two to four features is the likely way in, not measured.
- 1.2.0's error "`eig_sym_sel` failed to converge" is gone, replaced by
  two plain errors (`dev/eigen_conversion.md`, "Phase 4, as carried out").
  Finite input whose squares overflow, unscaled values near 1e160 for one,
  which raised it in 1.2.0, now stops with "the cross-product of the data
  is not finite; values this large overflow when squared, rescale the
  input", before either solver; a test pins its class and wording. A
  failure of LAPACK itself reads "the symmetric eigensolver failed (LAPACK
  dsyevr, info = <n>)"; no input is known to raise it.
- A failed triangular solve inside LOBPCG is no longer an error: that
  outer iteration goes to the exact solver, where 1.2.0 raised Armadillo's
  `std::runtime_error`, "solve(): solution not found". Armadillo's "system
  is singular" line on the console is gone, with the least-squares answer
  it announced (`dev/eigen_conversion.md`, Phase 2, step 4, and "Phase 4, as
  carried out"). No input is known to reach either without an injected
  failure, so this may need no NEWS text.
- `pca_imp()` under `solver = "lobpcg"`, and `"auto"` where it picks
  LOBPCG, returns values that differ in their last digits from the build
  before: the products of the iteration no longer follow Armadillo's choice
  between dgemv and dgemm (`dev/eigen_conversion.md`, "Phase 4, as carried
  out"). Measured worst on the Phase 0 capture: 9.25e-11 absolute, 3.3e-9
  relative on the imputed values. The exact solver is not affected.
- `pca_imp()` under `solver = "lobpcg"`, and `"auto"` where it picks
  LOBPCG, returns values that differ from the build before, and its
  `n_exact`, `n_lobpcg_ok` and `n_lobpcg_bad` attributes can differ:
  `A * Q` is carried through a QR by one dtrsm from the right
  (`si_la::solve_tri_right`), where it was dtrtrs on the transposed
  operands (`dev/pca_imp_structure.md`, "What the wrappers do"). Measured
  on the reference capture, the same under CRAN's RcppEigen and the
  candidate: 5 of 70 cases move, all LOBPCG; worst on the imputed values
  9.26e-10 absolute, 2.88e-8 relative; `criterion_final` by up to 1.07e-5
  relative; no solver count changes. On
  `pca_imp(sim_mat(200, 5000)$input, ncp = 5, solver = "lobpcg")` under
  `set.seed(1)`: 1.2e-8 absolute, 4.4e-8 relative, and one outer iteration
  that went to the exact solver no longer does (`n_exact` 11 to 10,
  `n_lobpcg_bad` 1 to 0). On the 68 edge cases of
  `dev/eigen_p2s4_edges.R`, which ask LOBPCG for up to 41 eigenpairs: 51
  not bitwise, worst on the imputed values 4.12e-8 absolute, 1.27e-6
  relative; `criterion_final` by up to 1.25e-4 relative; the solver counts
  differ in 7 cases.
- `pca_imp()` under either solver returns values that differ in their last
  digits from the build before: the trace of the cross-product, the sum of
  the squared singular values and the sum of the row weights are plain sums,
  no longer in Armadillo's order (`dev/eigen_conversion.md`, "Phase 4, as
  carried out"). Measured worst on the Phase 0 capture: under the exact
  solver 6.0e-15 absolute on the imputed values; under LOBPCG 6.24e-10
  absolute, 3.4e-8 relative. `criterion_final` moves by up to 3.4e-6
  relative.
- `pca_imp()` refuses, with a new error, an input whose eligible columns
  hold more than `.Machine$integer.max` missing values
  (`dev/eigen_conversion.md`, Phase 2, step 6). No test or measurement
  raises it: it needs a matrix of more than that many cells.
- RcppArmadillo has left `LinkingTo` and `-DARMA_64BIT_WORD=1` has left
  `src/Makevars` (`dev/eigen_conversion.md`, Phase 3), so the line already
  in `NEWS.md` saying the package is compiled with 64-bit Armadillo indices
  is no longer true. Building from source no longer needs RcppArmadillo.
- `pca_imp()`'s `criterion_final` is `NA` when the fit converged on its
  objective alone; it was NaN, Inf or a number made by rounding
  (`dev/eigen_conversion.md`, "Phase 4, as carried out").
- The errors raised from C++ are reworded to the *Writing R Extensions* 1.7
  form (`dev/WRITING.md`, section 1): lowercase, no final period, argument
  names in single quotes. Two reach every imputer and `tune_imp()`: an
  infinite value now reads "'obj' holds an infinite value at row <r>,
  column <c>. Infinite values are not supported", and an all-`NA` column,
  "All NA/NaN column detected" in 1.2.0, reads "'obj' has a column with no
  observed value, only NA or NaN. Remove that column before imputation".
  A caller matching the old text breaks.
- **Breaking:** `Depends` is `R (>= 4.4.0)`, raised from 4.3.0. The C++
  needs C++17, which *Writing R Extensions* makes the default from R 4.4.0
  (`dev/eigen_conversion.md`, "The rest of the *Writing R Extensions*
  pass").

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

Agreed 2026-11-18, to run once `dev/eigen_conversion.md` had landed, which it
did on 2026-11-19; each as its own change. Vendored files (`dev/vendored.md`)
are out of scope. The gate of a change is `dev/eigen_p4_gate.sh`: the
reference capture against HEAD and the suite, under CRAN's RcppEigen and the
candidate. Each finding is a proposal with a measurement, put to the
maintainer before it is made; nothing here licenses a hand-rolled kernel
(`dev/CPP.md`) or a call that R's headers do not declare.

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
  zeros; and `row_w` is copied before it is checked. Profile first: after
  the second half of step 4 of Phase 2, `sim_pca_lobpcg` took about 1.5%
  longer than at `ff0e394` over 11 repetitions of the step 4 timing (a
  0.483-0.491 s, b 0.490-0.501 s), and chr21 showed nothing measurable.
  Not profiled. Suspected, an assumption until profiled: the zero-fill of
  every new `Block` in `src/lobpcg_warm.h` and the gathered sub-blocks.
  The copies in
  and out of Armadillo were removed in step 5 and timed 1.003 against step
  4, so they do not account for it, to the driver's resolution of about
  half a percent. Deferred to this pass by the maintainer, 2026-11-18,
  with any `/simplify` calls.
- **How R calls the same LAPACK routines.** Read R's own sources (the
  LAPACK module behind `eigen()`, `qr()`, `svd()` and `La.svd()`), and
  where R and `src/si_lapack.cpp` differ on options, workspace or checks,
  take R's form unless a measurement says otherwise. R's usage is what
  CRAN's machines exercise most. NOT READ yet.
- **Each routine as if chosen fresh.** One symmetric eigensolver or two
  (dsyevr and dsyevd); whether the values-only SVD is needed at all;
  whether a product symmetric by construction should be dsyrk. Fewer
  routines is fewer declarations to depend on in the oldest supported R's
  `R_ext/Lapack.h`, and fewer paths through an external BLAS or LAPACK on
  a check machine. Settled, not to reopen: `syevr_range` asks dsyevr for
  its workspace on every call (maintainer, 2026-11-19).
- **The R callers, now that the return types can change.** Whether the R
  sides of `pca_imp()` and `knn_imp()` are tidied to what the kernels can
  return without Armadillo. Not decided.
- **Symbol visibility.** The DLL exports 433 names, where registration
  with `R_useDynamicSymbols(dll, FALSE)` leaves R needing only
  `R_init_slideimp` (*Writing R Extensions*, "Controlling visibility").
  `PKG_CXXFLAGS = $(CXX_VISIBILITY)` in `src/Makevars` hides the rest on
  Linux and macOS; on Windows `$(CXX_VISIBILITY)` is empty, and the
  manual's route is a `.def` file in `src/` exporting `R_init_slideimp`
  alone. Advice, not a requirement (maintainer, 2026-11-19: fine as it is,
  otherwise for this pass). What it changes on Linux is not measured.

## `criterion_final` below `threshold`: deferred

Deferred by the maintainer, 2026-11-19 (`dev/eigen_conversion.md`, "Phase
4, as carried out", "`criterion_final` as `NA`"): whether `criterion_final`
should be `NA` whenever the objective is below `threshold`, which also
covers a rounding criterion that lands below it. Today it is `NA` only when
the criterion itself is not below `threshold`.

## Test suite: cleanup and tiers for local and CRAN runs

Asked for by the maintainer 2026-11-19. A plan only: every step below is a
proposal, ruled on before a test is moved or dropped. Vendored tests, if
`dev/vendored.md` lists any, are out of scope.

What there is today, counted 2026-11-19: 22 test files, 9115 lines, 279
tests, 2851 expectations. One tier switch, the skip in
`tests/testthat/helper.R` on `MANUAL_TESTS`, used 16 times in 4 files
(`test-tune_imp.R` 7, `test-group_imp.R` 4, `test-pca_imp.R` 4,
`test-crate.R` 1). No test calls `skip_on_cran()`, so CRAN runs everything
that is not manual. Time per file or per test is NOT MEASURED.

1. **Measure first.** Time every test on one machine with a testthat
   reporter that gives times, in a scratch clone, and list the slowest. For
   each test write down in one line what it pins: public API, error
   wording or class, a numeric result, a regression of a named bug, or an
   internal. Measure line coverage of `R/` and `src/` before anything
   moves (`covr` is not a declared dependency, so from a `dev/` script).
2. **Three tiers, proposed.**
   - CRAN: what a caller pins on. The exported API, the error surface, and
     small numeric checks with tolerances, on at most 2 cores, no mirai
     daemons, nothing that depends on a BLAS build's last digits. A time
     budget for the whole tier is set by the maintainer from step 1.
   - Local, run by `devtools::test()`: the CRAN tier and the regression
     tests pinned to bugs, the larger inputs, the internals. Switched with
     testthat's own `skip_on_cran()`.
   - Manual, as today: `MANUAL_TESTS`, through `check_manual()`.
3. **Cleanup, each kind listed for a ruling before any edit.** Tests that
   repeat another; tests that cannot fail (a subset never shown non-empty,
   a comparator never shown to fail); comparisons against a stored
   reference by `identical()`; tests and comments that name Armadillo or
   pin its behaviour; tests of functions no longer exported
   (`test-mean_imp_col.R` is one to look at); files past a thousand lines
   that hold several subjects.
4. **Rules for the pass.** A test is dropped only with a note of which
   test still covers the behaviour, and coverage is measured again after
   and compared with step 1. A regression test pinned to a bug moves tier,
   it is never dropped. The error surface stays in the CRAN tier. One
   commit per file or per kind, the suite run under both RcppEigen
   versions with and without `NOT_CRAN` set, so that both tiers are seen
   to pass.
5. **Write it down.** `.claude/CLAUDE.md`, "Layout", gains the invariant:
   which tier a new test goes in and how each is run.

## Take cli out; audit every message against the R manuals

Planned by the maintainer, 2026-11-19; not started, and its own change, not
part of `dev/eigen_conversion.md`. No entry for it existed before this one.

- cli leaves the package: errors, warnings and messages become base
  `stop()`, `warning()` and `message()`. Until this lands, new code in the
  PCA chain raises plain errors with no class, and existing cli calls are
  left alone.
- In the same pass, every message the package raises, from R and from C++,
  is audited against what the R manuals say about error messages
  (`dev/r_manual/`), and `dev/WRITING.md` is rewritten to that guidance, so
  the two do not disagree afterwards.
- To settle before starting: what becomes of the `slideimp_infeasible`
  class, which `group_imp()` and callers catch by class; the house rule
  "Errors and messages go through cli" in `.claude/CLAUDE.md`, which has to
  change first; and the tests that pin wording.
- The NEWS entry above gains a line for whatever wording or class changes.
- `rlang` goes in the same pass (maintainer, 2026-11-19). `R CMD check`
  notes it in `Imports` but never imported.

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
  goes in its own private library, never over the CRAN one:
  `dev/eigen_install_candidate.R` takes the library as an argument,
  installs from the RcppCore r-universe only when the candidate is missing
  or older, and prints the Eigen version its headers define. A run selects
  it with `R_LIBS` set before R starts; `dev/eigen_p4_gate.sh`, stage
  suite, does that for the suite. "Passes" means the reverse-dependency
  run's standard, the suite and the check, not a compile. Record both
  versions tested in the release commit.
- The reference capture on the BLAS of the rhub Linux images, compared
  against one taken here. Until then the default tolerances of
  `dev/ref_compare.R` are a reading from one BLAS, on one machine.
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
- **`USE_FC_LEN_T` before the R headers.** `src/mipca.cpp` and
  `src/mipca_crps.cpp` call BLAS through `R_ext/BLAS.h` and pass `FCONE`,
  but neither defines `USE_FC_LEN_T` before its first R header, and neither
  falls back with `#ifndef FCONE`, as *Writing R Extensions*, "Fortran
  character strings", writes the call (`dev/r_manual/R-exts.md`). The
  manual notes that `USE_FC_LEN_T` is the default from R 4.3.0, older than
  the oldest R the package supports (4.4.0), so the lengths are passed
  today. Ruled 2026-11-18: a file
  that disagrees with the manual is a bug to fix, not a model.
  `src/si_lapack.cpp` has the manual's form. Fixed upstream in `b44a3b5`,
  after the last copy; the next copy brings it and deletes this item.
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
