# Review: Phase 2, the PCA chain off Armadillo

Checkpoint "Phase 2, before Phase 3" of `dev/eigen_conversion.md`. Range
under review: `59920f4..597d28e` on `si-mi`, the seven commits of steps 1 to
6 (`17e6415`, `905e1ac`, `eade032`, `ff0e394`, `99a8c5e`, `8c15416`,
`597d28e`), working tree clean at the start. `/code-review` at high effort,
one session, no verify pass.

`dev/eigen_conversion.md` was deleted when its last phase landed; read it
with `git show 4f2b669:dev/eigen_conversion.md`. Most scripts this review
ran were retired once that plan closed; read them with
`git show aea6898:dev/review_p2_run.sh`,
`git show aea6898:dev/review_p2_cases.R`,
`git show aea6898:dev/review_p2_compare.R`,
`git show aea6898:dev/review_p2_show.R`,
`git show aea6898:dev/review_p2_tag.sh`,
`git show aea6898:dev/review_p2_suite.sh`,
`git show aea6898:dev/review_p2_fix_check.R`,
`git show aea6898:dev/review_p2_fix_check.sh`,
`git show aea6898:dev/review_p2_fix_suite.R`,
`git show aea6898:dev/review_p2_crash_events.ps1`,
`git show aea6898:dev/probe_p2s4_run.sh`,
`git show aea6898:dev/probe_p2s6_run.sh`,
`git show aea6898:dev/probe_p2s6_niter.R` and
`git show aea6898:dev/audit_eigen_wrappers_chain.sh`. The suite they drove,
`dev/review_p2_suite.R`, stays: `dev/eigen_p4_gate.sh` runs it.

Not a fresh session. The plan asks for one; the maintainer pasted the
command into the session that wrote step 6 and the failure-path half of
step 4, and chose to run it there. What this note says about those two
steps is their author reading them again, and the audit of this note falls
to the same session. The measurements below do not depend on who ran them:
each is a script in `dev/` and a folder of output.

Nothing under `src/` or `R/` was edited, nothing was committed, and nothing
was built in the working tree: every build is a clone in a scratch folder
outside the repository, installed with `dev/eigen_time_install.R` or loaded
there with `load_all1()`. R 4.6.2 ucrt, GCC 14.2.0, R's LAPACK 3.12.0,
RcppEigen 0.3.4.0.2, RcppArmadillo 15.6.2.1 (linked, used by no file at the
head). Run 2026-11-18 and 19. Every figure is a reading on this one machine
with this one BLAS and LAPACK.

Probe scripts, all in `dev/`, each run from the package root:

- `dev/review_p2_run.sh <scratch> <stage>` - every build and run, by stage:
  `build` clones `59920f4` ("old") and `597d28e` ("new") and installs each
  into its own library; `cases` runs `dev/review_p2_cases.R` on both and
  `dev/review_p2_compare.R` on the two outputs; `capture` takes the Phase 0
  capture on a fresh clone of `597d28e` and compares it with the stored
  ones; `edges` and `diag` run the plan's own `dev/probe_p2s4_run.sh` and
  `dev/probe_p2s6_run.sh` with an old commit of the reviewer's choice;
  `time` runs the timing driver, old as a and new as b.
- `dev/review_p2_cases.R <lib> <out dir> <group>` - 462 calls of the public
  `pca_imp()` on one installed build, in seven groups (`tiny`, `tinysweep`,
  `general`, `big`, `validation`, `overflow`, `dup`), each result or error
  in its own file. The `overflow` and `dup` cases each run in a child
  process with a 120 s limit.
- `dev/review_p2_compare.R <old dir> <new dir> [<regex>]` - one line per
  case: identical, identical once integers are read as doubles, or the
  worst absolute and relative difference and the attributes that differ.
- `dev/review_p2_show.R` - the type of each count in a `mipca` fit, what
  changed across the `tinysweep` group, and the names behind the positions
  `dev/ref_compare.R` prints for a diagnostics case.
- `dev/review_p2_tag.sh <scratch> <commit or tag> "<groups>"` - one more
  build and its cases: `v1.2.0`, `905e1ac` and `eade032`.
- `dev/review_p2_suite.sh <scratch>`, with `dev/review_p2_suite.R` - the
  suite gate on a fresh clone of `597d28e`.
- `dev/review_p2_crash_events.ps1` - reads the Windows Application log for
  crash events naming Rterm.exe. It found none (finding 1).

## Summary

| # | Where | Finding | Severity |
|---|---|---|---|
| 1 | `dev/eigen_conversion.md`, Phase 4 and step 3; `tests/testthat/` | The overrun step 3 ended is not silent: on duplicated samples the release, `v1.2.0`, and every commit before step 3 end the R process, 11 of 18 cases, default arguments included. The plan says every such run returned normally and that no test through `pca_imp()` can see it. There is no regression test | high |
| 2 | `src/lobpcg_warm.h`, `matmul_into()` | A Gram of order 2 to 4 under `solver = "lobpcg"` no longer returns what it did: 267 of 280 fits differ, worst 5.2e-9. The map names the cause; no step record lists it as a change, no reminder line carries it, and it was not put to the maintainer | medium |
| 3 | `src/svd_triplet.h`, the stop after `hybrid_topk_eig()` | A Gram holding Inf or NaN is reachable from finite input (values near 1e160 unscaled, 1e200 scaled). The notes leave that unmeasured; no test raises the stop; its text says the solver failed to converge | medium |
| 4 | `src/pca_kernel.cpp`, `sum_two_partial()` | The third two-accumulator sum rests on a ruling given for another sum and is missing from Phase 4's list of what to remove | low |
| 5 | `dev/to-do.md`, the NEWS reminders | Two changes to the error surface have no reminder line: the class and message of a failed triangular solve, and the new refusal of more than `INT_MAX` missing cells | low |
| 6 | `dev/to-do.md`, the reminder for `ncp` above 32 | "About 1e-14" under the exact solver is not a bound: 1.6e-13 absolute and 1.2e-11 relative measured on one fit | low |
| 7 | `src/lobpcg_warm.h`, `lobpcg_detail::Block` | The plan's table names Eigen for four of these files; the phase added none and wrote a matrix type with eleven helpers by hand. No record says the plan changed | low |

Finding 1 is about what is written down, tested and released, not about a
defect in the range: the code of the range is the fix. No finding is a
wrong result at the head.

What holds up, measured or read:

- Both ends of the range install with plain `R CMD INSTALL`: 511
  `-Wignored-attributes` lines each, and one `-Wmaybe-uninitialized` from
  RcppEigen's own `SelfadjointMatrixVector.h`, the same line in both, so not
  from this range. Under `load_all1(timer = FALSE)` the head compiles with
  the 511 and nothing else, as the step records say.
- The capture, retaken on a fresh clone of `597d28e`: 71 cases. Against the
  stored `59920f4` at tolerance 0, every case bitwise and 24 failed, each on
  the type of `n_iter` alone (`dev/probe_p2s6_niter.R`: 24 carry the
  attribute, a double before and an integer after, equal counts, the two
  results `identical()` once it is read as a double). Against the stored
  `p2s6`: 71 of 71 bitwise, none failed. Against `76d4c60` at the defaults:
  those 24 and the 7 deliberate failures of Phase 1. So the whole phase is
  bitwise on the capture, as each step claims of itself.
- The 60 edge cases of `dev/eigen_p2s4_edges.R`, `eade032` against
  `597d28e`, which is steps 4 to 6 in one comparison the plan did not make.
  With the old build left as it was, 6 cases are not bitwise: the five the
  plan gives for the dsyevd workspace (`big_xb_ncp16`, `25`, `17_em`,
  `big_chr21_ncp20`, `30`; worst 6.2e-8 absolute) and `big_xb_ncp40`. With
  the old build given the workspace of each problem's own query
  (`SLIDEIMP_P2S4_FRESH`), one: `big_xb_ncp40`, 8.1e-9 absolute, the rule
  for more than 32 components. Of the 53 results that carry `n_iter`, 52
  are `identical()` once it is read as a double. The event logs: 1207
  `lobpcg_solve()` calls at the head, and the calls end differently in 4 of
  55 cases, `big_xb_ncp40` and the three NaN injections step 4 records.
  Seven cases end in an error in either build, with the classes and
  messages step 4 gives; the old build prints Armadillo's "system is
  singular" line 6 times (the four singular triangles, and two of the NaN
  cases), the head none. Every failure branch of `src/lobpcg_warm.h` is
  taken by injection at the head: 17 failed QRs, 17 failed dsyevd of which
  8 are refusals of NaN, 4 solves by the dgelsd fallback, 2 failed solves,
  2 first Rayleigh-Ritz problems over `X` alone.
- The same probe cannot start from `59920f4`: its instrumented clone of
  that commit dies at `dup_xd_ncp07`. That is finding 1.
- The 49 diagnostics cases of `dev/eigen_p2s6_diag.R`, `59920f4` against
  `597d28e`, the whole phase where the plan compared step 6 alone. Of the
  44 that return a value in both, 25 are bitwise identical. 15 differ in
  one element, `eig_lobpcg_max_rel_res`, by at most 3.3e-24 absolute and
  4.4e-16 relative: every LOBPCG case, and the plain sum of squares in
  `column_norms()` that step 4 says is not bitwise
  (`src/lobpcg_warm.h`: `inline std::vector<double> column_norms(`). The
  other 4 are the cases with 33 and 40 components. The five injected SVD
  failures stop in both builds or in the new one alone, as step 6 says.
- The two boundaries, through `pca_imp()` on installed builds
  (`dev/review_p2_cases.R`, group `big`): `ncp = 32` exact is identical and
  `ncp = 33` is not; `ncp = 15` under LOBPCG is identical and `ncp = 16`,
  17 eigenpairs, is not (9.9e-11 absolute, 3.1e-9 relative, in the range
  step 4 gives).
- 22 calls over the public switches (group `general`): row weights, both
  methods, no scaling, `clamp`, `subset`, `nb.init = 3` with a seed under
  both solvers, `.progress`, an integer input, the auto probe on a 260 x
  300 input, duplicated columns and rows under the exact solver, a matrix
  with no hole, the `mipca` output. 20 are identical once integers are read
  as doubles. The other two are the `mipca` fit, which is `all.equal()` at
  tolerance 0 between the builds and differs in the type of two counts,
  `chunks[[1]]$info$n_iter` and `chunks[[1]]$prov$iter`, double before and
  integer after, as step 6 says.
- 22 calls with arguments the kernel's new `int` parameters could have
  treated differently (group `validation`): `maxiter` and `miniter` of 3e9,
  0 and NA, `ncp` NA and 2.5, the LOBPCG controls at 3e9 and at
  `.Machine$integer.max`, row weights too short, zero, at the floor and
  very unequal, an `Inf` in the input. 14 end in an error with the same
  class and message in both builds; 8 return a value, 7 identical once
  integers are read as doubles and one, `ncp = 39` of 40 columns, with the
  same cells and a `criterion_final` that moved by 2.4e-7 relative, the
  rule for more than 32 components again. `checkmate::assert_int()` refuses
  3e9 before the kernel is reached, as the comment in
  `check_lobpcg_fields()` says.
- Timing, old as a and new as b, the four PCA workloads, 11 repetitions,
  ratios of medians: `chr21_pca_exact` 1.001, `chr21_pca_lobpcg` 1.008,
  `sim_pca_exact` 0.944, `sim_pca_lobpcg` 0.974. Over the whole phase
  nothing is slower: the 1.5% step 4 measured on `sim_pca_lobpcg` is still
  inside the last figure, under the gain of step 6. The products of the
  plan's own step ratios, steps 3 to 6, are close to it: 0.934 and 0.983
  for the two `sim_mat()` workloads.
- The suite on a fresh clone of `597d28e`, manual tests included,
  `load_all1(timer = FALSE)`, the mirai daemons loading that clone
  installed into a private library: 22 files, 276 tests, none failed,
  errored or skipped.
- Armadillo, read in `dev/lib/arma`, version 15.6.2: `trace()` in
  `fn_trace.hpp` and `arrayops::accumulate()` take two partial sums over
  the even and the odd positions and add them; `op_dot::direct_dot()` runs
  `direct_dot_generic()`, the same two partial sums, up to 32 elements and
  BLAS above; `gemm::apply_blas_type()` and `gemv::apply_blas_type()` test
  for square operands up to 4 x 4 before BLAS. `gram_trace()`,
  `sum_two_partial()` and `sum_squares_two_partial()` are those loops.
- Each BLAS call of `src/pca_linalg_utils.h`, `src/gram_ops.h`,
  `src/svd_triplet.h` and `reconstruct_cols()` against the `arma::dgemm_`
  or `arma::dsyrk_` call it replaced: the same routine, transposes,
  dimensions, leading dimensions and offsets. Each expression of
  `src/lobpcg_warm.h` against the Armadillo expression it replaced, for
  the order in which a cell is computed: `P - X * XtP - Ract * RtP` as two
  subtractions cell by cell, `X * VX + pp` as product then sum,
  `0.5 * (M + M.t())` from `M` as it was, the four `-=` between the two
  QRs as one call with `alpha = -1`. The dimensions handed to
  `si_la::qr_econ`, `si_la::solve_tri` and `si_la::syevd` match what
  `src/si_lapack.cpp` reads.
- No narrowing: every cast to `int` in the range follows a check of the
  value (`n_missing`, a row index below `nrow()`, the iteration count below
  `maxiter`), and every offset is a product in `std::size_t`. A `Block`
  that was moved from keeps its dimensions with no cells; none is read
  again before it is assigned. Nothing under `R/`, `tests/` or the tracked
  dev scripts still reads the kernel's old `imputed_values`.
- `air format .` changed nothing tracked. The ASCII check prints
  `.github/workflows/pkgdown.yaml`, which the range does not touch.
  `uv run dev/check_map_anchors.py` reports 82 problems, the known ones in
  `dev/r_manual/R-exts.md` and `dev/review_eigen_wrappers.md`.

Not checked: a second BLAS or LAPACK, Linux, R 4.3.0; the RcppEigen
candidate (the stored `-eigen5` captures were not retaken); the `err_alloc`
stops and the refusal of more than `INT_MAX` missing cells, which nothing
here can raise; the CRAN binary of 1.2.0 (the tag was built from source
here); which status the wrapper returned in finding 3.

## Findings

### 1. The overrun step 3 ended kills the R process, and the notes say it does not

Step 3 moved `eig_sym_sel()` onto the wrapper
(`src/eig_sym_sel.h`: `si_la::syevr_range(PEIG_UPLO,`),
which holds the `n` eigenvalue cells dsyevr documents. Before it, the chain
gave dsyevr `topk` cells, and dsyevr writes past them when the range cuts a
group of equal eigenvalues. The plan treats that as an overrun nobody has
seen do harm:

The three sentences below are quoted from the notes as of `597d28e`; each
was reworded by the fixes ("Audit and fixes" at the end).

- "The review of the wrapper layer" measured it and ends, in
  `dev/eigen_conversion.md`: "10 and 40. Every run returned normally; where
  the cells land was not looked at."
- Phase 4 concludes from that that no test through `pca_imp()` can see it,
  in the same file: "every run measured to write past the buffer returned
  normally".
- The reminder for the NEWS draft reads, in `dev/to-do.md`: "no longer
  writes past a buffer when asked for more components".

READ, in the same section of the plan: those runs were of a clone "whose
`src/eig_sym_sel.h` hands dsyevr an `n`-cell buffer and counts what is
written past `topk`" (`dev/audit_eigen_wrappers_chain.sh`). In that build
nothing is written past any buffer, so "returned normally" says nothing
about the build that has the defect.

MEASURED, on installed builds with no patch
(`bash dev/review_p2_run.sh <scratch> cases dup`,
`bash dev/review_p2_tag.sh <scratch> <ref> dup`). Group `dup` is the six
rank-deficient inputs of the plan's own edge script
(`dev/eigen_p2s4_edges.R`: `xd <- base6[rep(seq_len(6), times = 5), ]`),
built the same way, each under `solver = "lobpcg"` as that script calls
them, `"exact"` and the default `"auto"` with every other argument at its
default: 18 cases, each in its own process.

| build | cases that returned | cases whose process ended with no result |
|---|---|---|
| `v1.2.0`, the release | 7 | 11, child status 5 |
| `59920f4`, the start of the range | 7 | 11, child status 5 |
| `905e1ac`, step 2 | 7 | 11, child status 5 |
| `eade032`, step 3 | 18 | 0 |
| `597d28e`, the head | 18 | 0 |

The 11 are the same in all three: on `xd`, 30 samples of which 6 are
distinct, `ncp` 5, 7 and 12 and `method = "EM"` with `ncp` 7 under the
exact and the default solver (8 cases), and the last three under LOBPCG.
`pca_imp(xd, ncp = 5)`, nothing else passed, is one of them. The 7 that
return are `xd` at `ncp` 5 under LOBPCG held to 15 iterations, and the six
on `xdc`, where the plan counted one cell past. For those 7, `59920f4` and
the head return the same values, identical once `n_iter` is read as a
double.

`xd_ncp05_auto` run three more times on `v1.2.0`, straight from the shell:
"Segmentation fault", exit 139, no output, three of three. The instrumented
clone of `59920f4` that `dev/probe_p2s4_run.sh` builds died too, at
`dup_xd_ncp07`, in both of its runs, which is how this was found.
`dev/review_p2_crash_events.ps1` finds no Application Error event for
Rterm.exe, so the exception code is not confirmed; status 5 is the low
byte of 0xC0000005, an access violation, and that reading is an inference.

What follows from it:

- Step 3 did not tidy a contract. It fixed a crash that the CRAN release
  has, on input a user can have: replicate samples, with `ncp + 1` above
  the rank they leave. The step record, the map's "The eigenvalue buffer of
  dsyevr" and the reminder line all describe less than that.
- A test through `pca_imp()` can see it. `.claude/CLAUDE.md` asks for "the
  regression test pinned to that bug" when a bug is fixed; there is none.
  READ: no PCA test file, and not the helper, builds an input with
  repeated rows (searched for `rbind(` and `rep(` used as an index). The
  input is 30 x 80 and the fit takes 0.01 s. Shown to fail is already
  done: the same call ends the process on `905e1ac`.
- Phase 4 keeps one oracle check alive for this ("guard cells after the
  eigenvalue buffer") on the ground that nothing else can catch it. That
  ground is gone.
- "A patch for CRAN is a separate thread" was said when the overrun had
  been seen to do nothing
  (`dev/eigen_conversion.md`: `Step 3 of Phase 2 ended the overrun on this branch.`).
  The measurement is the reason to look at that thread again.

Only Windows was measured. Where the cells land depends on how the
compiler laid out the frame, so another platform may crash on other input
or on none.

### 2. A Gram of order 2 to 4 under LOBPCG returns other values, and no record lists it

Every product of `src/lobpcg_warm.h` went through one call, to
`si_la::matmul`, removed since in Phase 4, from the function that is now
`product_into()` (`src/lobpcg_warm.h`: `inline void product_into(`),
and the wrapper's header said what it did not reproduce: for square
operands up to 4 x 4 Armadillo runs its own loops.
The map carries the consequence in one sentence, ending "needs a Gram of
at most 4 x 4 under LOBPCG." (in `dev/pca_imp_structure.md` until Phase 4
rewrote "What the wrappers reproduce";
`git show 141a9d0:dev/pca_imp_structure.md` has it).
That is all. Step 4's record lists what changes under
(`dev/eigen_conversion.md`: `What changes, none of it reachable from`)
and this is not there; `dev/to-do.md` has reminder lines for 17 eigenpairs
and for `ncp` above 32 and none for this; no ruling names it. Its gate is
bitwise because neither the capture nor the edge cases hold a Gram that
small.

READ: two shapes take Armadillo's own loops. One eigenpair asked of a Gram
of order 4 or less, where `A * X` is a square matrix times one column
(`method = "EM"`, `ncp = 1`); and every eigenpair asked of one, where `X`
is square (the regularized method with `ncp` one less than the eligible
column count, on a tall input).

MEASURED (`dev/review_p2_cases.R`, groups `tiny` and `tinysweep`, old
against new):

- `tiny`, 12 shapes under `solver = "lobpcg"` with no warm-up and a
  threshold of 1e-12: 7 differ, the 7 inside the rule above except order 2
  with both pairs asked, which is identical. The 4 outside it (order 5;
  order 4 with 2 or 3 pairs) are identical, and so are all 12 under the
  exact solver. So the rule is the boundary.
- `tinysweep`, 40 inputs for each of 9 shapes, `solver = "lobpcg"` and
  every other argument at its default: of the 280 fits inside the rule, 267
  differ (34 to 40 of 40 a shape); of the 80 outside it, none. Over the
  267: worst 5.2e-9 absolute (`w4_em1_s16`, a 4 x 9 input, 690
  iterations), worst 1.0e-9 relative (`w4_em1_s08`), 4 with a cell past
  both 1e-12 absolute and 1e-10 relative, the defaults of
  `dev/ref_compare.R`. `n_iter` and `converged` are the same in all 360.
  Where every pair is asked for, `n_exact` and `n_lobpcg_ok` differ in 71
  of the 80 fits and `n_lobpcg_bad` in 77.

Reach: the caller has to name the solver, since `"auto"` takes the exact
one below order 250 (`R/pca_imp.R`: `(n_gram < 250L || (k_eig / n_gram) > 0.10)`).
`group_imp()` with `solver = "lobpcg"` over groups of two to four features
is the plausible way in. Not measured through it.

The sizes are those of the two changes the maintainer did rule on, and
Phase 4 removes `matmul`'s copy of Armadillo's choice anyway. The finding
is that this one went in without the ruling and the line the plan requires
of every such change. Recommendation: a sentence in step 4's record, a
reminder line, and the maintainer's word; or, if bitwise matters here until
Phase 3, Armadillo's two small-matrix loops reproduced in the wrapper
(NOT READ: the loops), where the oracle's 42 products of that kind would
then hold them to the bit.

### 3. A Gram that is not finite is reachable from finite input

The stop is one line
("top-k symmetric eigensolver failed to converge" in `src/svd_triplet.h`, at
the commit reviewed; Phase 4 has since replaced it with two texts).
Step 3 put a new way into it, the wrapper's refusal of NaN and Inf, and
says of the input that could take it, in `dev/eigen_conversion.md` as of
`597d28e`: "put one there, not measured, and no test does."; the map said
the same, in `dev/pca_imp_structure.md`: "whether any input does is NOT
MEASURED". Both were reworded by the fixes.

MEASURED (`dev/review_p2_cases.R`, group `overflow`: a 30 x 10 input of
ordinary shape times a power of ten, each case in its own process with a
120 s limit, both solvers, old and new):

| input | result, the same in both builds |
|---|---|
| x 1e150, `scale = FALSE` | a value, identical once integers are read as doubles |
| x 1e160, `scale = FALSE` | error: top-k symmetric eigensolver failed to converge |
| x 1e200, `scale = FALSE` | the same error |
| x 1e200, `scale = FALSE`, `ncp = 9` of 10 columns, every eigenpair | the same error |
| x 1e200, `scale = TRUE` | the same error |
| x 1e-200, `scale = FALSE` | error: infeasible, no eligible column (class `slideimp_infeasible`) |

12 of 12 cases returned within the limit in each build, with the same
class and message in both. So the input exists: `check_finite()` passes
it, the squares overflow inside the fit, and scaling does not prevent it at
1e200. The call the wrapper review measured not to return did not show on
these inputs in the old build, the full range included.

Two things follow. The path has a public way in and no test; one line of
`expect_error()` on such an input would pin it, and would have shown
step 3 changing nothing a user sees there. And the message tells a user
with overflowing data that a solver failed to converge. That text is older
than the range and the range did not change it; with the cause now known
it is a candidate for the cli error Phase 4 plans. Not read: whether the
head's stop comes from the refusal or from dsyevr.

### 4. `sum_two_partial()` is kept on another sum's ruling and is not listed for removal

`dev/CPP.md` bans a manual multi-accumulator reduction. The range adds
three: `gram_trace()` (`src/svd_triplet.h`: `inline double gram_trace(`),
ruled in step 5; `sum_squares_two_partial()` in `src/pca_kernel.cpp`,
ruled in step 6; and the sum of the row weights, `sum_two_partial()` in
the same file. (Item 10 of Phase 4 made all three plain sums and removed
the last two functions, so they are named here without a fragment.)
READ: for the third the record does not claim a ruling. It says the
question was put in the same report
(`dev/eigen_conversion.md`: `and stays in Armadillo's order on the same footing.`).
And Phase 4's "Scaffolding to remove" names the first two
(`dev/eigen_conversion.md`: `tail energy, in the order of Armadillo's own`)
and not the third, so the list that phase works from would leave a banned
loop in `src/`.

The measurements behind it are in the record and were not rerun here: a
plain sum leaves the capture bitwise and moves 4 of 5 unequal-weight
diagnostics cases, at most 3.9e-12. Recommendation: the maintainer's word
on this sum by name, and its line in Phase 4. The three loops are one loop
written three times in two files; one helper would make Phase 4's removal
one edit.

### 5. Two changes to the error surface have no reminder line

The plan's rule is that an error changes only where the maintainer
decides, with a line in the NEWS entry of `dev/to-do.md` for each. READ,
that entry has no line for either of these:

- A failed triangular solve
  ("lobpcg_solve: triangular solve failed" in `src/lobpcg_warm.h`, at the
  commit reviewed; Phase 4 has since removed that stop).
  MEASURED by the edge probe's two injections: the condition was
  `std::runtime_error`, "solve(): solution not found", and is
  `Rcpp::exception` with the text above. The change is deliberate and
  recorded in step 4; `.claude/CLAUDE.md` counts wording and class as
  public. No input reaches it without an injection, which may be the
  reason for no line. If so the plan should say that such a change needs
  none.
- More than `INT_MAX` missing cells among the eligible columns
  (`src/pca_kernel.cpp`: `has a row for each and an R matrix holds at most %d rows`),
  a new error. Not raised here. Before the range the same input reached
  Rcpp's conversion of a matrix with that many rows; what that did is not
  read.

Armadillo's line on the console for a singular triangle is gone as well.
Step 4 records it; it was never a condition, so it has no class to pin.

### 6. The reminder for `ncp` above 32 gives a figure that is not a bound

In `dev/to-do.md` as of `597d28e`, reworded by the fixes: "in the last
digits under the exact solver (about 1e-14), by about 1e-8".
The plan's cases gave at most 8.2e-15 absolute and 9.7e-14 relative under
the exact solver, and the diagnostics probe from `59920f4` gives the same
here. MEASURED on another input (`dev/review_p2_cases.R`, `exact_ncp33`:
80 x 70, 500 holes, 202 iterations): 1.6e-13 absolute and 1.2e-11
relative over its 5600 cells, 499 of them changed, with `ncp = 32` on the
same input identical. Still far inside what `threshold` resolves, and under
`dev/ref_compare.R`'s default relative tolerance. For the draft: the figure
grows with the iteration count, so "last digits" is the claim that holds,
not the exponent.

### 7. The plan names Eigen for these files, and the phase wrote a matrix type by hand

"Where Armadillo is used today" gives `src/hybrid_topk_eig.h` and
`src/svd_triplet.h` as becoming Eigen, and `src/pca_kernel.cpp` and
`src/lobpcg_warm.h` as Eigen in part; "What was decided" says where
(`dev/eigen_conversion.md`: `Eigen is used where a small dense step reads more simply than a raw`).
READ: no file of the chain includes an Eigen header at the head. The
kernel holds `std::vector<double>`, and LOBPCG holds its matrices in a
type of its own (`src/lobpcg_warm.h`: `struct Block`) with eleven helpers
around it: `copy_of`, `mul` twice, `tmul`, `sub_mul`, `sub_cells`,
`add_cells`, `mul_plus`, `transpose`, `block_rows`, `put_block` and
`put_block_t`. Each step record says it adds no Eigen; none says the table
changed, or why.

There are reasons one could give: a raw buffer shows the order of every
operation, which the bitwise proof needs, and it keeps the chain clear of
the RcppEigen 0.3.4 to 0.4.9.9 switch. Neither is written down. The cost
is about 180 lines of container code to keep, with Armadillo's member
names (`memptr()`, `n_rows`) for a library that is leaving, and a gate,
"Both RcppEigen versions", that ran at every step on code holding no
Eigen.
Recommendation: one paragraph in the plan saying the chain stays on raw
buffers and what that decides for the table and for Phase 4's hygiene
pass; or, if Eigen is still meant, an item in Phase 4 for `Block`.
`.claude/CLAUDE.md`, "Reach for the solution a dependency publishes", is
the rule this sits beside; RcppEigen is already linked for the sampler.

## Audit and fixes

2026-11-19, the same session, on the maintainer's word to fix all seven as
recommended above. Everything above this heading is the note as written
against `597d28e`, except that six citations of sentences the fixes
reworded are now quotations. Nothing under `src/` or `R/` changed, and
nothing was built in the working tree. Scripts:
`dev/review_p2_fix_check.sh`, which takes the scratch folder and a stage
(`installed`, `clone` or `suite`), with `dev/review_p2_fix_check.R` (one of
the two new tests against one installed build) and
`dev/review_p2_fix_suite.R` (the two PCA test files on a clone).

| # | Outcome |
|---|---|
| 1 | Fixed. The regression test is in `tests/testthat/test-pca_imp.R`; the sentences in the plan, the map and the NEWS reminder say what was measured |
| 2 | Recorded and accepted as it is: step 4's record, the map, a reminder line. No code change |
| 3 | Fixed. A test pins the stop; the plan and the map say it is reachable; Phase 4 has the wording as a choice |
| 4 | Recorded: step 6's record says no word was given by name and that it is accepted since; Phase 4's list names the sum. The three loops were not folded into one helper |
| 5 | Fixed. Two reminder lines, and a third found on the way |
| 6 | Fixed. The reminder gives both readings and says neither is a bound |
| 7 | Recorded and accepted: the chain stays on raw buffers; "What was decided", the table's caption and a Phase 4 item for `Block` say so |

MEASURED, for the two tests
(`tests/testthat/test-pca_imp.R`: `test_that("pca_imp returns on replicate samples with ncp above their rank", {`, `test_that("finite input whose squares overflow stops before the eigensolver", {`),
each run alone in its own process against five installed builds:

| build | replicate samples | overflow |
|---|---|---|
| `v1.2.0` | segmentation fault, exit 139 | 3 expectations pass, then an error: the message differs |
| `59920f4` | segmentation fault, exit 139 | 5 pass |
| `905e1ac` | segmentation fault, exit 139 | 5 pass |
| `eade032` | 14 pass | 5 pass |
| `597d28e` | 14 pass | 5 pass |

So the first test fails on every build with the defect, by ending the
process, and passes from step 3 on. Its four fits raise no warning at the
head. It asserts its own premise: the completed data has lower rank than
the number of eigenpairs asked for, in each of the four calls.

Found on the way: on the overflow input `v1.2.0` stops with "`eig_sym_sel`
failed to converge", where every build of this branch says "top-k
symmetric eigensolver failed to converge". `c7a253b` changed the text,
after the release and before the range, so it is not a finding against the
range; it had no reminder line and has one now. On `v1.2.0` the second
test passes its three premises and fails on the wording, as a test that
pins today's wording should.

The gate, on a clone of `597d28e` holding the working tree's tests,
`load_all1(timer = FALSE)`, manual tests included: the two PCA test files,
39 tests and 352 expectations, none failed, errored or skipped; the whole
suite, the mirai daemons loading the clone installed into a private
library, 22 files, 278 tests, none failed, errored or skipped.

Still open after the fixes: whether 1.2.x is patched for CRAN, which
finding 1 reopens and the maintainer decides; which status the wrapper
returns on the overflow input (by reading, the refusal of a triangle that
is not finite, before dsyevr); and everything under "Not checked" above.
