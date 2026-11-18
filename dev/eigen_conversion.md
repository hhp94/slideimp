# Removing Armadillo from src/

The plan for taking RcppArmadillo out of the package's own C++. Work happens
on the `si-mi` branch, beside the work of `dev/si_mi_conversion.md`
(maintainer, 2026-11-18); there is no separate branch for it. Delete this
file when Phase 3 lands; git history keeps it.

"The start commit" below is the head of `si-mi` when Phase 0 takes its
reference capture. It is not `e5b3cf7`, where the capture scripts were first
run: `si-mi` has changed what `pca_imp()` returns since, on purpose, so the
reference for this plan is captured again, under its own label.

## Why

To cost the CRAN team less. A package that links to RcppArmadillo is
rebuilt and checked whenever RcppArmadillo is submitted, and a translation
unit holding both Armadillo and Eigen is among the heavier compiles on every
check flavor. RcppEigen stays, since the vendored Spectra sampler needs it.
For the same reason, the change must not cost CRAN more than it saves: a
platform-specific failure or a resubmission undoes the point, so careful,
verified steps win over speed.

It ships as its own patch release, the one after 1.2.0, and only after a full
local CI pass on the rhub Linux images.

RcppEigen is itself moving from Eigen 3.4 to Eigen 5
(https://github.com/RcppCore/RcppEigen/issues/151), and its maintainer is
patching the reverse dependencies that fail on the candidate one by one.
slideimp 1.2.0 on CRAN does not link to RcppEigen; the release that first
does must not become one of those failures. So everything here is built and
checked against both RcppEigen versions (see "Both RcppEigen versions").

## What was decided

- Code that uses Armadillo only to hold numbers - containers, indexing,
  element loops - becomes plain Rcpp: raw pointers into R memory, `R_xlen_t`
  or `std::size_t` indices, the same loops.
- Every BLAS and LAPACK call goes through R's own headers, `R_ext/BLAS.h` and
  `R_ext/Lapack.h`, with `F77_CALL` and `FCONE`, the way `src/mipca.cpp`
  already does. That covers the calls the PCA chain makes today through
  Armadillo's wrappers (`dgemm_`, `dsyrk_`, `dsyevr_`, `dsyevd_`) and the
  LAPACK routines behind Armadillo's own algebra (`qr_econ`, `solve`, `inv`,
  `svd`). This is the interface Writing R Extensions documents for BLAS and
  LAPACK, not a non-API entry point.
- Those calls live in one place, the wrapper layer (see "The wrapper
  layer"), never at the call sites. The hard part of this plan is reading
  what an Armadillo call does and reproducing it; the layer puts all of it
  in one small file with its own test, and leaves the rest of the PCA chain
  as substitution.
- Armadillo stays available to dev scripts, and only to them, as the oracle
  the wrappers are tested against (see "The Armadillo oracle"). Nothing from
  it is vendored and the package does not link to it after Phase 3.
- Accepted (maintainer, 2026-11-18): calling R's BLAS and LAPACK directly
  puts slideimp among the packages CRAN notifies when R changes that
  interface, a change RcppArmadillo's maintainer absorbs for its users
  today. `src/mipca.cpp` already carries that for BLAS; this plan adds
  LAPACK. With the wrapper layer such a change is made in one file.
- Eigen is used where a small dense step reads more simply than a raw
  buffer: workspace matrices, `Eigen::Map` over a buffer, element-wise
  expressions. A matrix product that Armadillo sent to BLAS still goes to
  BLAS, never to Eigen's own product, so the arithmetic stays the same.
- The OpenMP flags in `src/Makevars` and `src/Makevars.win` are the
  maintainer's to change and are not part of this plan.

## Out of scope

The vendored files (`dev/vendored.md`): the sampler and `src/Spectra/`. They
use no Armadillo; `src/mipca.cpp` says it must never include it.

## Rules that apply to every step

- Read `dev/CPP.md` before the first edit, and the structure map of the chain
  being changed before touching it: `dev/knn_imp_structure.md` for the K-NN
  chain, `dev/pca_imp_structure.md` for the PCA chain. A step that changes
  cited code fixes the map in the same commit, and
  `uv run dev/check_map_anchors.py` passes before it is committed.
- What R sees does not change. Each internal export returns the same R type,
  the same dims and the same attributes as at the start commit; the harness
  (Phase 0) records them on the unchanged build, so this is checked, not
  remembered. Tidying an R caller - dropping the `[1, ]` in `R/col_vars.R`,
  or the `as.integer()` in `R/mat_miss.R` - is a separate change, later, and
  the maintainer's call.
- 64-bit indexing stays. `-DARMA_64BIT_WORD=1` is in the Makevars so that a
  matrix of more than 2^31 cells indexes correctly; the replacement indexes in
  `R_xlen_t` or `std::size_t`, and never narrows to `int` where the Armadillo
  code used `arma::uword`.
- K-NN workers run on RcppThread. A worker reads raw pointers taken before
  the parallel region and never touches an Rcpp object or the R API.
- One commit per step. Before each: `air format .`, the whole suite with
  `load_all1(timer = FALSE)` and `check_manual()`, and the Phase 0 compare
  script against the reference. The timer must be off because the sampler
  refuses more than one chain in a `-DLOC_TIMER` build.
- No line numbers in this file or the maps; cite code by fragment.
- A step that adds or changes Eigen code is also built, and the suite and
  the compare script run, against the RcppEigen candidate before it is
  committed (next section). A failure there is fixed in that step, in code
  that compiles on both versions, never with a version switch.

## Both RcppEigen versions

Two versions, as read on 2026-11-18; both move, so look them up again
before relying on these numbers:

- CRAN: RcppEigen 0.3.4.0.2, Eigen 3.4. What users and CRAN build with
  today.
- Candidate: RcppEigen 0.4.9.9-5, Eigen 5.0.1, from the RcppCore
  r-universe that issue 151 names (it named 0.4.9.9-2 when the issue was
  opened). What they build with after the switch.

Nothing is vendored and nothing is pinned: the package links to whichever
RcppEigen is installed, and must pass on both. The candidate is installed
into its own private library, never over the CRAN one, by
`dev/eigen_install_candidate.R`, which takes the library as an argument,
installs only when it is missing or older than the r-universe's, and prints
the Eigen version the headers define. A run selects it with `R_LIBS` set
before R starts, the way the timing driver separates its builds. Each
capture records and prints `packageVersion("RcppEigen")` and the library it
came from, so a log says which one it tested.

"Passes" means the reverse-dependency run's standard, not a compile: the
whole suite and, at the end, `devtools::check()`. Compiling is necessary and
proves little.

The breakages issue 151 lists are mostly API: `MappedSparseMatrix`,
`Eigen::all` and `Eigen::last`, functor traits, element-wise `pow()`. The
Eigen this plan adds is `Eigen::Map` over dense buffers and coefficient-wise
expressions, which is the part least likely to move, but that is an
assumption until the candidate build of each step says so.

## The wrapper layer

One header and one source file under `src/`, named when they are written.

- The source file includes R's headers, `R_ext/BLAS.h` and
  `R_ext/Lapack.h`, and nothing from Armadillo or Eigen. The header declares
  plain functions over `double*`, dims and a status: no LAPACK declaration,
  no Armadillo type, no Eigen type. The PCA translation unit can then
  include it beside `RcppArmadillo.h` without the redeclaration warnings
  `src/eig_sym_sel.h` describes, and calls move one at a time. NOT COMPILED:
  this is an assumption until the first wrapper builds without a new
  warning.
- One wrapper per routine the chain calls: `dgemm`, `dsyrk`, `dsyevr`,
  `dsyevd`, and the routines behind `qr_econ`, the triangular `solve` and
  `svd`. Which LAPACK routines those three are is NOT READ; each is written
  down in `dev/pca_imp_structure.md` when it is.
- Each wrapper does its own workspace query, passes Fortran character
  arguments with `FCONE`, and returns LAPACK's `info` as its status. It
  takes and returns plain C++ types and does not call the R API: raising
  the R error stays with the caller, which keeps the failure path the
  Armadillo code had.
- LAPACK takes `int` sizes. A wrapper receives sizes in `R_xlen_t` or
  `std::size_t`, and a size or workspace length that does not fit an `int`
  is reported as a failure, never narrowed.
- It reproduces the Armadillo call it replaces: the same routine, the same
  options, the same behaviour on failure, as READ in armadillo_bits.
- A wrapper is a function call where Armadillo's was inline. Beside a BLAS
  or LAPACK call that should not show; the timing driver measures it.
- Both Makevars already link `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`, read on
  2026-11-18, so nothing new has to link.

## The Armadillo oracle

Dev only, run by hand, like `dev/check_impute_knn.R`: once RcppArmadillo
leaves `DESCRIPTION`, nothing under `tests/` may call it. It stays after
Phase 3, so a wrapper changed later can still be checked against the call it
replaced.

- `dev/oracle_install_armadillo.R` installs RcppArmadillo from CRAN into a
  private library, given as an argument, when it is missing or older than
  CRAN's; `--check` reports both versions and installs nothing. Nothing from
  the library enters the repository.
- A C++ file in dev/, compiled by `Rcpp::sourceCpp()` against that library,
  pulls in the wrapper source and exposes each pair to R: the wrapper, and
  the Armadillo call it replaces. Armadillo and `R_ext/Lapack.h` in one file
  give the redeclaration warnings; in this file they are accepted.
- An R script runs every pair on ordinary inputs and on degenerate ones:
  singular, rank-deficient, one column, one row, and whatever else the
  chain can pass. Per routine it prints the worst absolute and the worst
  relative difference, with the cells the relative figure was taken over,
  asserts on the maximum, and asserts that the pair agree on whether the
  call failed. It asserts each input set is non-empty.
- It must be shown to fail: with a deliberately wrong wrapper (the other
  triangle, say) it reports that routine.
- The latest RcppArmadillo, not a pinned one. The script prints
  `packageVersion("RcppArmadillo")`, and `dev/pca_imp_structure.md` records
  the version the wrappers were READ against. A difference that appears
  after an update is checked against what Armadillo changed first.
- The scripts are tracked, each with its `!` line in `dev/.gitignore`; the
  library is regenerable and is not.

## Who does which work

Split by kind of work, not by file. Opus does what is substitution; Fable
does what depends on reading exactly what an Armadillo or LAPACK call does,
or on chasing a numerical difference: a mistake there compiles, can pass a
loose tolerance, and shows up later as a failure on a platform that was not
tested. When a piece is in doubt, it goes to Fable.

| Work | Model | Why |
|---|---|---|
| Phase 0, what is still owed | Opus | Scripting against rules this file already states |
| The wrapper layer and the oracle | Fable | Reading armadillo_bits, the R interface (`F77_CALL`, `FCONE`, `int` sizes, workspace queries) and the failure paths |
| Phase 1, steps 1 and 2 | Opus | Storage only; the compare script settles bitwise-ness and the RNG question |
| Phase 1, step 3 (the K-NN chain) | Opus | Storage only, but threaded: the worker rule and 64-bit indexing above |
| Phase 2, every file: containers, and each Armadillo call replaced by its wrapper | Opus | Substitution, once the wrapper is tested |
| Phase 2, the call sites in `src/lobpcg_warm.h` and `src/armaSVD.cpp` where a factorization can fail | Fable | What the code does next must match what Armadillo's flag or fallback led to; LOBPCG iterates, so a last-bit change can grow |
| Any step whose results move beyond the noise floor | Fable | Finding which call moved them before going on |
| Phase 3 | Opus | Mechanical |

The wrapper layer touches only new files, so it is written while Phase 1
runs. Phase 1's K-NN chain is what the rest of `dev/si_mi_conversion.md`
waits on ("To settle before the code" there), and is not held behind the
wrappers. Once it is in, that plan's K-NN steps and this plan's Phase 2 run
on the same branch: Phase 2 stays under `src/` in the PCA chain, and a
session of one plan does not edit a file the other is changing.

## Review checkpoints

The gate on every commit is mechanical: the suite and the compare script.
Review is for what those cannot see: a failure path no capture exercises, an
index narrowed to `int`, a worker touching an Rcpp object, a workspace
length that overflows on a large input. So it is not done per commit, but
at four points, each before later work builds on what it covers:

| After | Reviewer | Looks for |
|---|---|---|
| Phase 0 | A fresh Fable session | That the compare can fail (`dev/ref_compare_selftest.R`, rerun), that the capture reaches every export in the table below, and that each solver path is actually taken by some case |
| The wrapper layer and the oracle, before any call site uses them | A fresh Fable session | Each wrapper against armadillo_bits, read independently; the oracle run, and shown to fail |
| Phase 1 | Either model, fresh | `int` narrowing, the worker rule, the RNG draw order |
| Phase 2, before Phase 3 | Fable for the failure-path call sites, either model for the rest | The whole diff of the PCA chain |

A reviewer is a fresh session: a different context matters more than a
different model. It gets the diff and the READ notes in
`dev/pca_imp_structure.md`, not the author's account of why the work is
right. It writes a review note in dev/, with its probe scripts. The
maintainer starts it with `/code-review`, at an effort they choose, with a
message saying what the checkpoint looks for and where the note lands.

The session that did the work up to the checkpoint is the auditor of that
note. `/code-review` starts without the author's context and can load
context unrelated to the change, so a finding may rest on code or rules
that do not apply here, and the reviewer has no verify pass of its own. The
auditor checks each finding against the code and this plan: it reruns the
claims the reviewer measured, measures the ones it left to inspection, and
gives a verdict per finding in chat, with corrections. What the audit
upholds is fixed, on the maintainer's word, before the next phase starts,
so each phase begins from reviewed code.

## Where Armadillo is used today

Written at `e5b3cf7`; on 2026-11-18 the only file under `src/` changed since
was the vendored `src/mipca.cpp`. "Storage" means containers, indexing and element loops only.

| File | What it uses Armadillo for | Becomes |
|---|---|---|
| `src/matrix_checks.h` | storage | plain Rcpp |
| `src/mat_stats.cpp` | storage | plain Rcpp |
| `src/sample_each_rep_cpp.cpp` | storage, and `arma::randu()` | plain Rcpp, `R::runif()` |
| `src/imputed_value.h` | storage, `arma::Mat<mask_t>` | plain Rcpp |
| `src/imputed_value.cpp` | storage | plain Rcpp |
| `src/impute_knn_brute.cpp` | storage | plain Rcpp |
| `src/pca_linalg_utils.h` | `dgemm_` through Armadillo | R's BLAS |
| `src/gram_ops.h` | `dsyrk_`, `dgemm_` through Armadillo | R's BLAS |
| `src/eig_sym_sel.h` | its own `dsyevr_` declaration, typed with `arma::blas_int` | R's LAPACK |
| `src/lobpcg_warm.h` | `qr_econ`, `solve` on a triangle, `norm`, `dsyevd_`, `dsyevr_`, `fliplr`, `reverse`, `find` | R's LAPACK, Eigen where simpler |
| `src/hybrid_topk_eig.h` | storage | Eigen |
| `src/svd_triplet.h` | storage, `trace` | Eigen |
| `src/armaSVD.cpp` | storage, `svd`, `find`, `cumsum`, `dgemm_` | R's BLAS and LAPACK, Eigen |

The PCA files are ONE translation unit: `src/armaSVD.cpp` includes
`src/svd_triplet.h`, which pulls in the rest. `src/RcppExports.cpp` already
includes `RcppArmadillo.h` and `RcppEigen.h` together, so a translation unit
holding both during the switch is not new.

## Phase 0: a reference to compare against

The scripts exist and are tracked: `dev/ref_capture.R`,
`dev/ref_compare.R` and `dev/ref_compare_selftest.R`, shared with the
`si-mi` branch. Their headers are the reference for how to run them; what
they were written to do:

- A capture script. On one build it runs every export whose signature uses
  Armadillo, and the R entry points that reach them (`knn_imp()`,
  `pca_imp()` with each solver, `col_vars()`, `mean_imp_col()`,
  `mat_miss()`, `sim_mat()`, the `tune_imp()` NA placement), on
  `sim_mat()` output and on `dev/GSE286313.EPICv1_chr21.qs2`, under fixed
  seeds. It saves every result, with its type, dims and attributes, to an
  ignored folder under dev/. That folder is regenerable data: rerun the
  script at the start commit.
- Its cases take every path of the PCA kernel that Phase 2 rewrites: exact
  and LOBPCG, the kernel's own auto choice (`pca_imp()` demotes auto to
  exact below 250 rows or eligible columns, so this needs a larger input),
  LOBPCG failing and falling back to exact, EM and regularized, tall and
  wide, with and without the Gram cache, row weights, random starts, and a
  run stopped before it converges. They raise each Inf error with the cell
  it names off the first row and column, and run `col_vars()` and
  `mean_imp_col_internal()` on two threads. The auto choice is timed, so it
  can differ between builds; that case is captured beside the same input
  under each forced solver, so a flip can be read against the one it now
  matches.
- A compare script. It reads a capture and checks a build against it. For a
  double result it reports `max(abs(a - b))` and `max(abs(a - b) / m)`, with
  `m = max(abs(a), abs(b))` per cell (over the cells where `m != 0`, saying
  how many), asserts the worst case against a tolerance, and reports whether
  the two are bitwise identical, byte for byte, without asserting it.
  Dividing by the larger of the two makes the relative figure the same
  whichever capture comes first; dividing by one side lets a cell that side
  holds as 0 pass while the other holds a value under the abs tolerance.
  Integer results, index maps and the NA placement use `identical()`. Type,
  dims and attributes must match exactly; attribute names are compared as a
  set, not in order. Messages are not captured, only warnings and errors.

The reference capture is labelled `76d4c60`, the head of `si-mi` when it
was taken (source commit `603fe33`), and every compare in this plan is
against it. It was taken in a clone of that commit outside the working
tree, since other sessions build here, and copied into `dev/ref/`. A capture
is regenerable: check out `76d4c60` in a clone and rerun.

The noise floor comes first: capture the reference build twice, in two R
processes, and compare them. A threaded BLAS need not be deterministic, and
no tolerance is chosen tighter than what the unchanged build reproduces
against itself. The header of `dev/ref_compare.R` records the result and the
default tolerances set from it, on one machine with one BLAS. Run
2026-11-18 at `76d4c60`, after the Phase 0 review added cases: 75 cases,
all bitwise identical (worst absolute and relative difference 0), and
`dev/ref_compare_selftest.R 76d4c60` caught the 13 planted changes at
tolerance 0 and the 12 expected at the defaults, and reported a 0 turned
into -0 as passing but not bitwise. The self-test fails on the compare as
it was before the review, which divided by the new capture only. A capture
of the same build with every export traced is bitwise identical too, and
showed each path above taken by at least one case: the kernel's auto mode
once (it chose LOBPCG), 12 failed LOBPCG iterations, 3 calls without the
Gram cache (1 tall, 2 wide), 1 unconverged run.

Then the same reference build against the candidate RcppEigen, compared to
the capture. At the start commit the package's own code holds no Eigen, so the
exports captured here should agree to the noise floor; a larger difference
is a finding to explain before Phase 1. This also proves the candidate
library and the `R_LIBS` selection work before any step depends on them.
Run 2026-11-18 against 0.4.9.9-5 (Eigen 5.0.1), as capture
`76d4c60-eigen5`: every translation unit compiled against the candidate's
headers, and all 75 cases are bitwise identical to `76d4c60`, the auto
choice included. The only compiler difference is the count of Eigen's own
`-Wignored-attributes` warnings, which both versions emit (511 and 543).

Every capture here is a `load_all1()` build, which compiles at `-O0` with
`NDEBUG` undefined, so Eigen's own asserts are on in the captures; the
timing builds below are installed at the user's `-O3`. Each comparison is
between builds of one kind.

Timing, in the same phase. Install the reference build into its own library,
and each converted build into another, with `dev/eigen_time_install.R`; it
writes the source commit and the RcppEigen and RcppArmadillo versions the
build compiled against into the library. `dev/eigen_time_builds.R`
alternates Rscript processes between the two libraries
(`dev/eigen_time_worker.R` is one process, and stops if slideimp loaded from
the wrong library) and reports medians over interleaved repetitions of
`pca_imp()` (both solvers) and `knn_imp()` on the chr21 data and on a larger
`sim_mat()`, 200 x 5000 by default. Each process makes one untimed call,
since the first call in a process takes 1.4 to 2.5 times a later one on
the chr21 PCA workloads, a cost common to both builds that would pull every
ratio towards 1. It then runs its call 5 times in one `replicate()` and
reports the total, since one PCA call on chr21 takes under 0.1s. Elapsed
is read from `Sys.time()`, which steps in about 2 us here; `proc.time()`
steps in 10 ms on Windows, a few percent of a chr21 PCA total. Both
builds are the same flags; `load_all1()` builds are not mixed with
installed ones. Removing Armadillo also removes the OpenMP Armadillo turns
on by itself in an OpenMP build, so this is the measurement of what that
costs or saves.

The driver's own noise, read 2026-11-18 on one machine with the BLAS and
OpenMP threads at their defaults: the reference installed twice and timed
against itself, 5 repetitions, 165s in all, before the warm-up call and
the finer clock. Five workloads gave a b/a ratio of medians between 0.97
and 1.05, but every time was a multiple of 10 ms, and the 1.043 of
`chr21_pca_lobpcg` is one tick (0.23 s against 0.24 s). `chr21_knn` gave
1.32, from a stretch of slow runs on one side in repetitions 1 to 3
whichever ran first. The two DLLs differ only in their header timestamp
and checksum, which is measured; that the machine caused the stretch is an
assumption, not measured.

Read again the same day with the current worker, same machine and
settings, 5 repetitions, 194s in all: all six workloads gave a b/a ratio of
medians between 0.988 and 1.008. Single runs still stray (one `chr21_knn`
run took 0.90 s against a median of 0.61 s), which the medians absorb. This
is the band a converted build's ratio is read against; a ratio outside it
is timed again before it counts as a finding.

## Phase 1: storage-only files to plain Rcpp

Expected to be bitwise unchanged, since the loops keep their order; the
compare script asserts on a tolerance anyway and reports bitwise-ness.

1. `src/matrix_checks.h`, and `src/mat_stats.cpp` with `col_min_max()`,
   `col_vars_internal()`, `mean_imp_col_internal()`, `check_inf()`,
   `check_finite()`, `col_miss_internal()`, `row_miss_internal()`.
   `src/imputed_value.h` includes `src/matrix_checks.h`, so its helpers take
   a pointer and dims, callable from Armadillo and Rcpp callers alike while
   the K-NN chain still holds `arma::mat`. `arma_uword_bytes()` stays until
   Phase 3.
2. `src/sample_each_rep_cpp.cpp`. Each `arma::randu()` becomes
   `R::runif(0, 1)`. READ: RcppArmadillo's rng/Alt_R_RNG.h defines
   `randu_val()` as `Rf_runif(0, 1)`. NOT READ: that the scalar
   `arma::randu()` reaches it. The compare script settles it: the sampled
   positions under one `set.seed()` must be `identical()`. The `umat` return
   comes back to R as a double matrix today, and the R side coerces it to
   integer (a comment in `tests/testthat/test-tune_imp.R` says why); the
   replacement keeps the double return, so nothing on the R side moves.
3. The K-NN chain: `src/imputed_value.h`, `src/imputed_value.cpp`,
   `src/impute_knn_brute.cpp`. `arma::datum::inf` and `arma::datum::nan`
   become `std::numeric_limits`; `arma::Mat<mask_t>` becomes a flat buffer.
   The early-abandon chunking in `calc_distance_raw` is algorithmic and
   stays (`dev/CPP.md`).

## Phase 2: the PCA chain to R's BLAS and LAPACK

The wrapper layer and its oracle come first, and pass their review
checkpoint before any file below calls a wrapper. Then leaves first, in this
order, each step compiling and passing on its own:

1. `src/pca_linalg_utils.h`
2. `src/gram_ops.h`
3. `src/eig_sym_sel.h`
4. `src/lobpcg_warm.h`
5. `src/hybrid_topk_eig.h` and `src/svd_triplet.h`
6. `src/armaSVD.cpp`, whose `pca_imp_internal_cpp()` is the export

For each Armadillo algebra call, first READ in RcppArmadillo's
armadillo_bits which LAPACK routine it calls, with which options, and what
it does on failure (READ: the general `solve()` warns and falls back to an
approximate solution on a singular or near-singular system; whether the
triangular path `src/lobpcg_warm.h` takes does the same is NOT READ;
`qr_econ()` returns a flag the caller tests), then call the same routine
with the same options and keep the failure path. Write down what was read
in `dev/pca_imp_structure.md` as the step lands. That reading is done once,
in the wrapper; a step below replaces the call with the wrapper and keeps
what the surrounding code does when it fails.

Each call moves to its wrapper in the step that converts its file.
`src/eig_sym_sel.h` explains why it declares `dsyevr_` itself:
`R_ext/Lapack.h` redeclares symbols Armadillo has already declared, with
different types, and warns in every translation unit that includes both.
The wrapper's source file is its own translation unit and its header
declares no LAPACK symbol, so the PCA translation unit never includes
`R_ext/Lapack.h`; that self-declaration goes in step 3, and step 6 drops
`RcppArmadillo.h` with the last Armadillo type.

If the separate source file turns out not to keep the warnings away, the
fallback is the order this plan had before: steps 1 to 5 replace the
containers while the calls keep Armadillo's declarations, and step 6 moves
every call at once.

Tolerances here come from Phase 0's noise floor, not from a guess. LOBPCG
iterates, so a last-bit difference early can move the converged answer by
more than one ulp; if a step moves results by more than the noise floor,
find out which call did it before going on.

## Phase 3: remove the dependency

- In `DESCRIPTION`, `RcppArmadillo` leaves `LinkingTo`.
- In `src/Makevars` and `src/Makevars.win`, `-DARMA_64BIT_WORD=1` leaves
  `PKG_CPPFLAGS`. The OpenMP flags are left as the maintainer has them.
- `arma_uword_bytes()` in `src/mat_stats.cpp`, its comment on the macro, and
  the test in `tests/testthat/test-mat_miss.R` that pins it are deleted.
- `devtools::document()` regenerates `src/RcppExports.cpp` without
  `RcppArmadillo.h`.
- `.claude/CLAUDE.md`, Layout, says `src/` is Rcpp with RcppArmadillo; it
  says what `src/` is built on instead. `git grep -i armadillo` over the
  tracked notes finds anything else.
- `devtools::check()`, once against CRAN's RcppEigen and once against the
  candidate current on that day, looked up again and not taken from this
  file. Record both versions in the commit body. Local CI on the Linux
  images is the maintainer's call.

## Not settled yet

- Whether the noise floor measured on one machine holds on the BLAS of the
  rhub Linux images. Until a capture is compared there, the default
  tolerances are a reading from one BLAS.
- Whether the R callers are tidied once the return types are free to change.
- What the timings say. They are recorded here when each phase has been
  run, as readings on one machine, not as numbers to pin. Phase 0's reading
  is the driver's noise band, under "Phase 0".
