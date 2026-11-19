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
  `R_ext/Lapack.h`, with `F77_CALL` and `FCONE`, as *Writing R Extensions*
  gives them. No file already in `src/` is the model: where one disagrees
  with the manual, that is a bug to fix (maintainer, 2026-11-18). That
  covers the calls the PCA chain makes today through
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
- The C++ is written for itself, not to keep its R callers unchanged
  (maintainer, 2026-11-18). An internal export takes and returns the type R
  naturally holds for each value: data as `Rcpp::NumericMatrix`, read in
  place; counts, indices and settings as `int` or `Rcpp::IntegerVector`; a
  result per column or per row as a vector, not a 1 x p or n x 1 matrix. No
  step reproduces how RcppArmadillo converted an argument (an index vector
  read as doubles and cast to `arma::uword`), the shape an Armadillo type
  came back to R in, or Armadillo's error text. A conversion or allocation
  in C++ whose only purpose is to spare an R caller a change is not made;
  the R caller changes, in the same commit.
- `std::size_t` or `R_xlen_t` is for offsets computed in C++, such as
  `col * n_rows`, so a matrix past 2^31 cells indexes correctly. A value
  that crosses from R is bounded by `nrow()` or `ncol()`, which R holds as
  `int`, and is `int`.
- An export checks the values it indexes with, as `src/imputed_value.h`
  checks its column groups, rather than trusting its caller: an index out
  of range reads past the matrix instead of failing.
- What the exported functions return, and their errors, change only where
  the maintainer decides, and each such change gets a line in the NEWS
  entry of `dev/to-do.md`. The capture compares them with the start commit.
  Its cases that call an internal export directly fail the compare when the
  export's interface changed; the step names each, and shows its values
  unchanged.
- Parallelism is spent in the K-NN kernel only; the column statistics in
  `src/mat_stats.cpp` run on one thread (maintainer, 2026-11-18). K-NN
  workers run on RcppThread. A worker reads raw pointers taken before the
  parallel region and never touches an Rcpp object or the R API.
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

One header and one source file under `src/`, which are `src/si_lapack.h`
and `src/si_lapack.cpp`. "The wrapper layer and the oracle, as written"
below records what they hold and how they were checked.

R's own guidance governs every wrapper, and wins over Armadillo's habits
and over convenience. It is read from *Writing R Extensions*, kept as
`dev/r_manual/R-exts.md`: the R-devel manual (4.7.0, under development,
2026-11-15), converted from its HTML with pandoc. R-devel is what CRAN
checks against, so it is the stricter reading. The sections that bind
the wrappers:

- "Using `Makevars`": `$(LAPACK_LIBS)` comes before `$(BLAS_LIBS)` and
  `$(FLIBS)`, in that order.
- "Fortran character strings": `USE_FC_LEN_T` is defined before any R
  header is included (and so before `Rcpp.h`, which includes them), and
  every Fortran character argument carries its length through `FCONE`.
- "Numerical analysis subroutines": the routines declared in
  `R_ext/BLAS.h` and `R_ext/Lapack.h` are not part of the official API,
  though unlikely to change, and routines have been added to them (most
  recently in R 4.5.0). A wrapper calls only a routine that R's headers
  declare in every R the package supports (`DESCRIPTION` says `R (>= 4.3.0)`),
  checked against R's NEWS, never one reached by declaring it ourselves.

Where the manual, the installed headers and this plan disagree, the work
stops and says so.

- The source file includes R's headers, `R_ext/BLAS.h` and
  `R_ext/Lapack.h`, and nothing from Armadillo or Eigen. The header declares
  plain functions over `double*`, dims and a status: no LAPACK declaration,
  no Armadillo type, no Eigen type. The PCA translation unit can then
  include it beside `RcppArmadillo.h` without the redeclaration warnings
  `src/eig_sym_sel.h` describes, and calls move one at a time. Measured, in
  "The wrapper layer and the oracle, as written".
- One wrapper per routine the chain calls: `dgemm`, `dsyrk`, `dsyevr`,
  `dsyevd`, and the routines behind `qr_econ`, the triangular `solve` and
  `svd`. READ: those are dgeqrf with dorgqr; dtrtrs with dtrcon, and dgelsd
  as Armadillo's fallback; and dgesdd. `dev/pca_imp_structure.md`, "What the
  wrappers reproduce", has each.
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
  give the redeclaration warnings; in this file they are accepted. As
  written they do not arise: the wrapper source is a unit of its own.
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

## The wrapper layer and the oracle, as written

As committed in `79cce0f`, 2026-11-18, before its review; no file of the
PCA chain calls a wrapper. Every figure is a reading on one machine with
one BLAS. The review changed the two eigen wrappers and the oracle, and
the oracle's and the self-test's counts below are the ones of that commit:
"The review of the wrapper layer" has what changed and the counts since.

- `src/si_lapack.h` declares nine functions in namespace `si_la`: `gemm`,
  `gemv`, `syrk`, `matmul`, `syevr_range`, `syevd`, `qr_econ`, `solve_tri`
  and `svd_values`. Sizes are `std::size_t`; the status is 0, LAPACK's
  positive `info`, or one of five codes of the wrapper's own. An illegal
  argument is refused before the Fortran routine is entered, since under R
  xerbla raises an R error from inside the library. Nothing calls the R
  API or throws; a failed allocation is a status.
- Two more than the list above. `*` on two `arma::mat` is not always dgemm:
  READ, it is dgemv when the product has one row or one column, which the
  products in `src/lobpcg_warm.h` are when one eigenpair is asked for, and
  `+=` and `-=` of a product are one call with `beta = 1`. `gemv` is that
  routine, and `matmul` chooses between the two as Armadillo does, so the
  choice is tested in one place and not rewritten at each call site.
- Not wrapped. `arma::norm(, "fro")` reaches dnrm2, and `arma::dot()`
  reaches ddot above 32 elements (READ). Ruled 2026-11-18: no wrapper for
  either, and a plain loop at each of the chain's three sites, in the step
  that converts the file: the dot of the top singular values in
  `src/armaSVD.cpp`, and the Frobenius norm of the projected block and the
  infinity norm of `A` in `src/lobpcg_warm.h`. None is in a hot path by
  operation count (not timed): each is lower order than the products beside
  it. The Frobenius norm will not be bitwise equal to dnrm2, which scales
  as it sums; the block is a projection of orthonormal vectors, so a plain
  sum of squares cannot overflow, and the value only feeds a threshold
  test. Armadillo's infinity norm calls no BLAS (READ). A difference shows
  in the reference capture at its tolerance.
- Where a wrapper is not Armadillo's call, each recorded in
  `dev/pca_imp_structure.md`: products of square operands up to 4 x 4 go
  to BLAS, where Armadillo runs loops of its own; `A.t() * A` on one object
  is not sent to dsyrk; the dgelsd fallback of the triangular solve is
  sized without ilaenv, whose routine name is a Fortran string longer than
  the one character `FCONE` covers; that fallback is reported by a flag,
  where Armadillo prints a warning; and an empty triangular system sets no
  flag. Against the chain's own dsyevr call, which is not Armadillo's:
  `syevr_range` asks for the workspace and allocates it on every call,
  where `src/eig_sym_sel.h` asks once and keeps the buffers. Not timed.
  Since the review, two more against the chain's calls: `syevr_range`
  holds the `n` eigenvalue cells dsyevr documents, where the chain gives
  it `topk`, and `syevr_range` and `syevd` refuse NaN and Inf in the
  triangle LAPACK reads, where the chain tests nothing.
- Where it is not this plan's wording: `svd_values` asks dgesdd for its
  workspace only from 1024 cells up, and `qr_econ` sizes dorgqr from
  dgeqrf's answer, because Armadillo does and a different length can put
  LAPACK on another code path. `dev/oracle_run.R` is given the library as
  an argument and does not need `R_LIBS`. `Rcpp::sourceCpp()` builds
  `si_lapack.cpp` as a translation unit of its own, so the oracle's unit
  never includes `R_ext/Lapack.h` and there is no redeclaration warning to
  accept.
- R 4.3.0 is not checked. R's NEWS says 4.5.0 extended the routines
  `R_ext/Lapack.h` declares without listing them, so the NEWS cannot
  settle which of the eight LAPACK routines used here R 4.3.0 declares.
  Ruled 2026-11-18: write against R-devel and the installed headers (R
  4.6.1 declares all of them); the check on the Linux images settles 4.3.0.
- Compiled outside the working tree, `-Wall -pedantic`: 511
  `-Wignored-attributes` warnings at `7d0c176` and 511 with the two files
  added, nothing else. The assumption under "The wrapper layer" is
  measured: `src/armaSVD.cpp`, made in the scratch copy to include
  `si_lapack.h` beside `RcppArmadillo.h` and call a wrapper, compiles with
  the same 511.
- Capture `wrap`, of a build that carries the new object: 71 cases, every
  one bitwise identical to `p1r`, and against `76d4c60` the 7 deliberate
  failures `p1r` has. The suite was not run for this step: no test reaches
  a wrapper, and nothing under `R/` or `tests/` changed.
- The oracle, against RcppArmadillo 15.6.2-1: 627 cases, all passing.
  Every pair that enters the same routine is bitwise identical: `gemm`
  154, `gemv` 46, `syrk` 119, `matmul` 83 and one that both sides refuse,
  `syevd` 46, `qr_econ` 17, `solve_tri` 36 (16 solved exactly, 12 by
  Armadillo's fallback, 4 failed on both sides, 4 empty), 4 in the form
  `src/lobpcg_warm.h` writes, `svd_values` 14 and 2 failed on both sides.
  The 42 products Armadillo runs in its own loops differ from BLAS by at
  most 8.9e-16 absolute and 3.8e-15 relative, 18 of them bitwise.
  `syevr_range` against `eig_sym()`, another algorithm: at most 1.7e-13
  absolute and 7.6e-10 relative over the cells above 1e-6 of the largest,
  34 cases. 29 refusals return the code expected.
- `dev/oracle_selftest.R` breaks one wrapper at a time in a copy, twelve
  ways (the other triangle in `syrk`, `syevr_range`, `syevd` and
  `solve_tri` among them); the oracle named the broken routine each time
  and passed on the unbroken source.
- Found on the way. `src/mipca.cpp` and `src/mipca_crps.cpp` do not define
  `USE_FC_LEN_T` before their R headers; they are vendored, so the fix is
  an upstream item in `dev/to-do.md`. `uv run dev/check_map_anchors.py`
  fails on `dev/r_manual/R-exts.md`, which cites paths of R's own source
  tree, and passes on every other note.

## The review of the wrapper layer

The second checkpoint of "Review checkpoints", on `79cce0f`, 2026-11-18. A
fresh session wrote `dev/review_eigen_wrappers.md` with its probe scripts
(`dev/review_eigen_wrappers_mutants.R`, `dev/review_eigen_wrappers_probe.R`
and its `.cpp`, `dev/review_eigen_wrappers_run.R`,
`dev/review_eigen_wrappers_procs.R`). The author's session audited it: every measured claim was rerun and
came back as the note gives it, and all eight findings were upheld. Fixed
the same day on the maintainer's word, on the working tree over `79cce0f`.
Every figure is a reading on one machine: R 4.6.1, R's LAPACK 3.12.0.

What was wrong in the wrappers, and is fixed:

- `syevr_range` passed the caller's `w`, documented as `iu - il + 1`
  cells, to dsyevr, which documents `W` as `n` cells. For a proper subset
  dsyevr stores every eigenvalue of its search interval there before it
  drops the surplus, so a range that cuts a group of equal eigenvalues
  writes past `iu - il + 1`: 34 cells over the reviewer's cases, 17 of
  them in one call. The wrapper now holds the `n` cells and copies the
  count out.
- The chain has the same sizing today, and 1.2.0 has it too (READ at the
  tag). Measured through `pca_imp()` by
  `dev/audit_eigen_wrappers_chain.sh`, on a clone of `79cce0f` whose
  `src/eig_sym_sel.h` hands dsyevr an `n`-cell buffer and counts what is
  written past `topk`: on 30 samples of which 6 are distinct, holes
  included, `ncp = 5` wrote 1 cell past in 1 of 5 calls, `ncp = 7` up to 4,
  `ncp = 12` up to 16, `ncp = 20` up to 9, and `method = "EM"` with
  `ncp = 7` up to 6; on 10 columns of which 3 are identical, `ncp = 8`
  wrote 1. Nothing was written on random data, on the duplicated data
  while `ncp + 1` stayed within its rank, or on the chr21 file at `ncp` 2,
  10 and 40. Every run returned normally; where the cells land was not
  looked at. `dev/to-do.md` carries it until step 3 of Phase 2.
- Neither eigen wrapper tested for NaN or Inf, as the chain's calls do
  not. dsyevr asked for every eigenpair of a matrix holding NaN did not
  return within 60 s, 4 of 4 cases at orders 6 and 60; a proper subset
  came back as `err_count`; `syevd` returned 0 with NaN eigenvalues in 12
  of 20 cases. Both now refuse NaN and Inf in the triangle LAPACK reads,
  before the routine is entered, and accept anything in the other
  triangle. `pca_imp()` asks for every eigenpair when there are fewer
  eligible columns than rows, `ncp` is that column count less one, and the
  method is the regularized one (measured: 41 calls for 6 of 6 on a 40 x 6
  input). Its input is tested for NaN and Inf, so only the
  arithmetic of the fit could put one in a Gram: not measured.

What was wrong in the oracle, and is fixed:

- Class "same" compared at a tolerance and only reported bitwise identity,
  so a wrapper that sent a product to another BLAS routine passed. A case
  of that class that is not bitwise identical now fails;
  `--same-by-tolerance` is the way out on a BLAS that does not repeat
  itself.
- `syevr_range` was compared only with `eig_sym()`, another algorithm, at
  1e-11 and 1e-7. It is now also compared, as `syevr_chain`, with dsyevr
  called as `src/eig_sym_sel.h` calls it, to the bit, vectors included.
- No input lay near a rule copied from Armadillo, and several failure
  paths had no case. Added: ranges that cut a group of equal eigenvalues,
  with guard cells after `w` and `z`; NaN and Inf in either triangle;
  a row times a column of 1000 terms, where the two dgemv routes differ
  in 184 of 200 random products; triangles whose condition estimate lies
  within 1000 of `DBL_EPSILON` on either side, checked on what Armadillo
  reports; a singular value between `DBL_EPSILON` and `n * DBL_EPSILON`;
  the dgelsd fallback at orders 30 and 130; Inf inside a triangle; the
  leading dimension of a transposed operand; an allocation that cannot
  succeed; and the status code a refusal carries.
- `dev/pca_imp_structure.md` left `auxlib::eig_sym_dc()` NOT READ. It is
  read: it refuses NaN and Inf in the upper triangle before LAPACK.

Counts since. The oracle: 748 cases, all passing, every class "same" case
bitwise identical: `gemm` 154, `gemv` 46, `syrk` 119, `matmul` 95 and one
that both sides refuse, `syevr_chain` 48 and 8 that both sides refuse,
`syevd` 48, `qr_econ` 17, `solve_tri` 52 (18 solved exactly, 24 by
Armadillo's fallback, 6 failed on both sides, 4 empty), 4 in the form
`src/lobpcg_warm.h` writes, `svd_values` 14 and 2 failed on both sides.
The 42 products Armadillo runs in its own loops: at most 8.9e-16 absolute
and 7.9e-15 relative, 15 bitwise. `syevr_range` against `eig_sym()`: at
most 5.7e-13 absolute and 4.0e-10 relative over the cells above 1e-6 of
the largest, 48 cases. dsyevr wrote past the count in 10 cases, 55 cells,
and the wrapper left every guard cell alone. 42 refusals return the code
expected. The build outside the tree has the same 511
`-Wignored-attributes` and nothing else; capture `wrapfix` is 71 cases
bitwise identical to `p1r`.

What the oracle does not catch. The reviewer's sixteen wrong wrappers: 14
passed the oracle of `79cce0f`, 3 pass this one.

- `svd_values` asking dgesdd for its optimum below 1024 cells too, and
  `qr_econ` giving dorgqr `max(m, n)` cells: no case shows either rule
  changing a result, so the oracle cannot say they are needed. That is
  evidence for Phase 4, which drops both.
- `syevr_range` without its eigenpair count check. With NaN and Inf now
  refused first, no input known reaches `err_count`; the check stays as a
  guard with no case.
- The workspace is not such a rule for dsyevr: a wrapper that ignores
  LAPACK's answer and passes `26 n` changes the eigenvectors at orders 60
  and 300, which `syevr_chain` catches.

`dev/oracle_selftest.R` now breaks one wrapper at a time 28 ways, the
twelve of `79cce0f` and sixteen for the checks above: the oracle named the
routine aimed at each time and passed on the unbroken source. Thirteen of
the reviewer's sixteen are among what it catches.

Not done, and why. The scan for NaN and Inf is not timed. The Armadillo
version of `src/eig_sym_sel.h` is not patched: a patch for CRAN is a
separate thread (maintainer, 2026-11-18), and the conversion does not wait
for it. Step 3 of Phase 2 ended the overrun on this branch. R 4.3.0 is not checked, as before.

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
| Phase 4 | Fable | Each item changes what a kernel does on purpose, with no Armadillo to compare against |

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
  `mean_imp_col_internal()` on two threads (Phase 1 step 1 removed both).
  The auto choice is timed, so it
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
   `col_vars_internal()`, `check_inf()`, `check_finite()`,
   `col_miss_internal()`, `row_miss_internal()`. `src/imputed_value.h`
   includes `src/matrix_checks.h`, so its helpers take a pointer and dims,
   callable from Armadillo and Rcpp callers alike while the K-NN chain still
   holds `arma::mat`. `mean_imp_col_internal()`, which no R code calls
   (`dev/si_mi_conversion.md`), is deleted rather than ported. RcppThread
   and `cores` leave the file, and `col_vars()` loses its `cores` argument.
   `arma_uword_bytes()` moves to `src/armaSVD.cpp`, so `src/mat_stats.cpp`
   includes no Armadillo.
2. `src/sample_each_rep_cpp.cpp`. Each `arma::randu()` becomes
   `R::runif(0, 1)`. READ: RcppArmadillo's rng/Alt_R_RNG.h defines
   `randu_val()` as `Rf_runif(0, 1)`. NOT READ: that the scalar
   `arma::randu()` reaches it. The compare script settles it: the sampled
   positions under one `set.seed()` must be `identical()`. It takes integer
   vectors and returns an integer matrix, so `sample_na_loc()` drops its
   `storage.mode()` coercion.
3. The K-NN chain: `src/imputed_value.h`, `src/imputed_value.cpp`,
   `src/impute_knn_brute.cpp`. `arma::datum::inf` becomes
   `std::numeric_limits`, and the NaN left in a cell the kernel cannot fill
   becomes R's `NA` (step 3 below); `arma::Mat<mask_t>` becomes a flat
   buffer. The early-abandon chunking in `calc_distance_raw` is algorithmic
   and stays (`dev/CPP.md`).

Step 1 was first committed as `1ac016a` with the Armadillo interfaces
reproduced: the same R types and shapes, `col_idx` cast from doubles, and
`cores`. The maintainer rejected that on 2026-11-18 (the type rules under
"Rules that apply to every step"), and steps 1 and 2 were redone together,
on the working tree over `1ac016a`.

- `src/mat_stats.cpp` includes `Rcpp.h` and nothing of Armadillo or
  RcppThread. Every export takes `obj`, a `const Rcpp::NumericMatrix &`.
  `col_vars_internal()` returns a double vector, `col_miss_internal()` and
  `row_miss_internal()` integer vectors; `R/col_vars.R` drops `[1, ]` and
  `R/mat_miss.R` drops `as.integer()`. `col_min_max()` keeps its 2 x p
  matrix and refuses a matrix with no rows in its own words.
- `col_vars()` loses `cores`, a breaking change noted in `dev/to-do.md`; the
  test that compared `cores = 1` with `cores = 4` goes with it.
  `mean_imp_col_internal()` is deleted. `arma_uword_bytes()` moves to
  `src/armaSVD.cpp`.
- `sample_each_rep_cpp()` takes integer vectors and an `int max_attempts`.
  `sample_na_loc()` already builds the vectors as integers, so they are not
  converted; `max_attempts` is, at the call, since its default is the
  double `100` and `checkmate::assert_count()`, which bounds it to an int,
  leaves its type alone. It returns an integer matrix; `sample_na_loc()` drops its
  `storage.mode()` line. It refuses an index or budget outside its range,
  NA included, and lengths that disagree, before it indexes with them; a
  new test in `tests/testthat/test-tune_imp.R` raises each refusal. Every
  `arma::randu()` is `R::runif(0, 1)`, and the NOT READ above is settled by
  measurement: the positions under one `set.seed()` are identical, so the
  draws are R's.
- `dev/ref_capture.R` drops the four cases whose call no longer exists
  (`col_vars_x2_cores2` and the three `int_mean_imp_col*`). Capture `p1s2`
  against `76d4c60`: 71 cases. Every case through an exported function is
  bitwise identical, `na_loc_x1`, `na_loc_x1_num`, `tune_x1_knn` and
  `tune_x1_pca` among them. Seven fail on purpose: the four that are
  missing, and `int_col_vars`, `int_col_miss` and `int_row_miss`, whose
  container changed. Read as doubles, those three hold identical values
  (worst absolute and relative error 0 over 300, 300 and 60 cells).
- `dev/eigen_p1_edges.R` reduces each result to its values as doubles and
  passes only arguments both builds take. Run on a clone of the start
  commit and on the step, compared at tolerance 0: 50 cases, 48 pass, all
  8 sampler cases among them. Two pass but are not bitwise: `col_min_max()`
  on a column whose minimum or maximum is a zero held as both 0 and -0 can
  return the other sign, since Armadillo kept one running value over even
  rows and one over odd. Its one caller, `sim_mat()`, rescales `rnorm()`
  output. Two fail on purpose: `col_min_max()` on no rows, where the
  Armadillo build raised `min(): object has no elements` for a 0 x 2
  matrix and returned a 2 x 0 one for 0 x 0. `sim_mat()` requires two rows.
- Not built against the RcppEigen candidate: neither step adds Eigen.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading the step installed into a private library: 22
  files, 275 tests, none failed, errored or skipped.

Step 3, on the working tree over `12ac4bd`:

- The three files include `Rcpp.h`, RcppThread and the standard library,
  nothing of Armadillo. `impute_knn_brute()` takes `obj` as a
  `const Rcpp::NumericMatrix &`, `k`, `method` and `cores` as `int` and the
  three groups as integer vectors, and returns `list(loc, value)`: `loc` an
  integer matrix of each missing cell's 1-based row and column, `value` a
  double vector. It returned an n x 3 double matrix holding the indices as
  doubles. `knn_imp()` writes `obj[imputed$loc] <- imputed$value`, and its
  `cbind()` of the two index columns goes.
- Where the kernel leaves a cell it now leaves R's `NA`, not NaN: `value` is
  filled with `NA_REAL` on the main thread, and a worker writes only where
  the weight total is positive. `knn_imp()` drops its pass turning NaN into
  NA. This departs from item 3 as first written, which had
  `arma::datum::nan` become `std::numeric_limits`: under the type rules the
  kernel returns the value R holds for a missing cell.
- Every index that crosses from R, and every row and column index derived
  from one, is `int`; column offsets are taken in `size_t` by `col_ptr()`
  (`src/imputed_value.h`), and so are the result offsets. More than
  2^31 - 1 missing cells in the imputed columns is refused, since `loc` is
  an R matrix. `validate_knn_inputs()` refuses a negative or NA group index
  in its own words; as `arma::uvec` it wrapped to a huge index, refused as
  out of bounds with that number in the message. `tests/testthat/test-knn_imp.R`
  raises both, and `k` of 0 and NA, and pins `loc` and the `NA`.
- The worker lambda captures plain pointers, `std` containers and scalars by
  name, so a use of an Rcpp object inside it does not compile; Rcpp's
  `operator[]` can call R's `warning()`. The masked copy and the mask are
  `std::vector`s, zero-filled before the copy writes every cell (listed for
  the hygiene pass in `dev/to-do.md`, with two other candidates).
- `stop_on_inf()` takes a `const Rcpp::NumericMatrix &`, since no Armadillo
  caller is left.
- Capture `p1s3` against `76d4c60`: 71 cases and the same 7 deliberate
  failures as `p1s2`. Every K-NN case, through `knn_imp()`, `slide_imp()`,
  `group_imp()` and `tune_imp()`, is bitwise identical.
- `dev/eigen_p1_edges.R` gains 20 cases calling the kernel directly, each
  reduced to (row, col, value) as doubles with `NA`: 5, 16, 17 and 40 rows
  around the chunk size, both metrics, `dist_pow` 0, 1 and 2, the largest
  `k`, two threads, a zero distance, no group 3, a target that shares no
  observed row with any candidate and rows no chosen neighbor observes (all
  9 cells left `NA` on both builds), `knn_imp()` on that input with and
  without `post_imp`, integer input with ties, and two refusals both builds
  word alike. 70 cases at tolerance 0: all 20 pass and are bitwise; the 2
  failures are `col_min_max()` on no rows, as in step 2.
- Not built against the RcppEigen candidate: the step adds no Eigen.
- Timed against the start commit, both installed by
  `dev/eigen_time_install.R` with the same flags, on the two K-NN workloads
  only, 5 repetitions a run, BLAS and OpenMP threads at their defaults; read
  2026-11-18 on one machine, two runs of about 145 s. `chr21_knn` gave a b/a
  ratio of medians of 0.945 and 0.932, the step faster in 9 of the 10 paired
  repetitions, on totals near 0.6 s for 5 calls. `sim_knn` gave 0.971 and
  1.050, faster in 4 of 10: outside the band on opposite sides, so no
  difference shows there, and the machine was noisier than when the band
  was read. The reference is the start commit, so this covers steps 1 to 3:
  `knn_imp()` also calls `check_finite()` and `mat_miss()`, which steps 1
  and 2 changed. Where the `chr21_knn` gain comes from is not measured. The
  PCA workloads were not timed for any step of this phase, though
  `pca_imp()` calls `col_vars()`, `check_finite()` and `mat_miss()`.
- The suite with the manual tests, as for steps 1 and 2: 22 files, 275
  tests, none failed, errored or skipped.

The Phase 1 review checkpoint (`dev/ad45ec5..5412876.md`) found no int
narrowing, no Rcpp object or R call inside a worker, and the `R::runif`
draws in the same order and number as Armadillo's. Its other findings were
fixed on the working tree over `5412876`:

- `col_ptr()`, the column-index refusal and a non-negative-count refusal
  live in `src/matrix_checks.h`, which states the int and `size_t` rule
  once. Every column offset in `src/mat_stats.cpp`,
  `src/sample_each_rep_cpp.cpp` and the K-NN chain goes through
  `col_ptr()`, which takes `size_t`, so an `int` argument is widened before
  the multiply. The K-NN and sampler index refusals share one wording; the
  K-NN one changed from "out-of-bounds column index".
- `sample_each_rep_cpp()` refuses a column listed twice in `pool_idx_in`,
  which could fill two slots from one column and return a cell twice.
  `sample_na_loc()` never passes one.
- The K-NN progress bar counts a column with no neighbor, so it reaches its
  total.
- `tests/testthat/test-tune_imp.R` pins the positions one seeded
  `sample_each_rep_cpp()` call draws, read from the reference build, so a
  change to the draws fails the suite and not only the capture.
- `README.Rmd` no longer credits RcppArmadillo and RcppThread to
  `col_vars()`.
- Capture `p1r` against `76d4c60`: the same 7 deliberate failures as `p1s3`,
  everything else bitwise. Edges at tolerance 0: 70 cases, the 2
  `col_min_max()` failures and now `knn_idx_out`, whose message is the new
  shared wording. The suite with the manual tests: 22 files, 276 tests, none
  failed, errored or skipped.

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
approximate solution on a singular or near-singular system, and so does the
triangular path `src/lobpcg_warm.h` takes; `qr_econ()` returns a flag the
caller tests), then call the same routine
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

Step 3 also ends a defect the review of the wrapper layer found in the
chain as it stands: dsyevr is given `topk` eigenvalue cells where it
documents `n`, and writes past them when the range cuts a group of equal
eigenvalues (`dev/to-do.md`, "The eigenvalue buffer of the chain's dsyevr
call"). The wrapper holds the `n` cells, so for finite input the step
changes no result and stops the overrun.

If the separate source file turns out not to keep the warnings away, the
fallback is the order this plan had before: steps 1 to 5 replace the
containers while the calls keep Armadillo's declarations, and step 6 moves
every call at once.

Tolerances here come from Phase 0's noise floor, not from a guess. LOBPCG
iterates, so a last-bit difference early can move the converged answer by
more than one ulp; if a step moves results by more than the noise floor,
find out which call did it before going on.

Step 1, on the working tree over `59920f4`:

- `src/pca_linalg_utils.h` includes `src/si_lapack.h`, `<cmath>` and
  `<cstddef>`, nothing of Armadillo or R. Every helper takes raw buffers and
  `std::size_t` dimensions; the callers in `src/svd_triplet.h`,
  `src/armaSVD.cpp`, `src/hybrid_topk_eig.h`, `src/eig_sym_sel.h` and
  `src/lobpcg_warm.h` pass `memptr()` and the dimensions of their
  `arma::mat`, which they keep until their own step.
- `gemm_nn()` and `gemm_tn_block()` call `si_la::gemm` with the arguments
  they passed to `arma::dgemm_`: the same routine, transposes, leading
  dimensions and offsets. Not `si_la::matmul`, which would send the
  one-column product of `ncp = 1` to dgemv. `gemm_nn()` no longer sizes its
  output: `src/svd_triplet.h` calls `U.set_size()` with the dimensions it
  had. `gemm_tn_block()` takes C as `a_cols x b_cols`, so its leading
  dimension is `a_cols` where it was `C.n_rows`; both callers size `V` as
  `Xhat.n_cols x ncp`, so the two are equal.
- Both return the wrapper's status, `[[nodiscard]]`, and refuse operands
  whose dimensions do not fit (`err_arg`); a non-zero status is an
  `Rcpp::stop()` at each of the three call sites, "internal error: dgemm
  refused its arguments". `arma::blas_int` narrowed a dimension silently.
  Kernel setup refuses dimensions past the `int` range first, so the error
  is not expected to be reachable, and no test raises it.
- `canonicalize_signs()` returns at once on zero rows, where it read
  `col[0]` of an empty column. Not reachable: `ncp + 1 <= min(n_rows,
  n_eligible_cols)` is checked at setup.
- Captures, each with the capture script of its own commit except
  `76d4c60`, retaken here with the script as of `ad45ec5` as Phase 0 took
  it: `76d4c60` (75 cases) and `59920f4` (71), in worktrees of those commits,
  and `p2s1` on this step. `p2s1` against `59920f4` at tolerance 0: 71
  cases, none failed, every one bitwise. Against `76d4c60` at the defaults:
  the 7 deliberate failures of Phase 1, everything else bitwise. Read
  2026-11-18 on R 4.6.2, RcppEigen 0.3.4.0.2, RcppArmadillo 15.6.2-1.
- The step adds no Eigen, so the plan does not ask for the candidate; it
  was run anyway, being cheap. Capture `p2s1-eigen5`, the same source built
  against RcppEigen 0.4.9.9-5 (Eigen 5.0.1) in a private library: 71 cases,
  every one bitwise identical to `p2s1`.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading the step installed into a private library: 22
  files, 276 tests, none failed, errored or skipped.

Step 2, on the working tree over `17e6415`:

- `src/gram_ops.h` includes `src/si_lapack.h`, `src/pca_linalg_utils.h`,
  `<cstddef>`, `<cstring>` and `<vector>`, nothing of Armadillo or R.
  `GramWorkspace` holds its dimensions as `std::size_t`, where
  `arma::blas_int` narrowed them silently. `GramCache` holds the fixed Gram
  and the wide case's scaled changing columns in `std::vector<double>`,
  zero-filled; `set_size` left the second uninitialized, and
  `impute_restandardize()` writes every cell of it before a Gram reads it.
  The helpers take raw buffers; `src/armaSVD.cpp` and `src/svd_triplet.h`
  pass `memptr()` of the `arma::mat` they keep until their own step.
  `copy_upper_cols()` takes the two leading dimensions it read from its two
  matrices, and `gram_changing_block()` takes the row count it read from
  `X_work`. `form_weighted_gram()` takes the cache `const`: it never wrote
  it.
- Every dsyrk and dgemm goes through `si_la::syrk` and `si_la::gemm` with
  the arguments it passed to `arma::dsyrk_` and `arma::dgemm_`.
  `syrk_upper()`, `gram_cache_init()` and `form_weighted_gram()` return the
  wrapper's status, `[[nodiscard]]`, and a non-zero one is an
  `Rcpp::stop()` at the two call sites: "internal error: dsyrk refused its
  arguments" after `gram_cache_init()`, "internal error: dsyrk or dgemm
  refused its arguments" after `form_weighted_gram()`. Not expected to be
  reachable, for the reason step 1 gives, and no test raises it.
- The timer build and the diagnostics build compile, with no warning from
  the chain's own files.
- Phase 0's inputs are built to take all four branches of
  `form_weighted_gram()`: `x1` wide and `xt2` tall with complete columns
  (`dev/ref_capture.R`: `xt2: tall, with complete columns, for LOBPCG.`),
  `xw` wide and `xt` tall with a hole in every column, which the script
  asserts. Capture `p2s2` against `p2s1` at tolerance 0: 71 cases, none
  failed, every one bitwise. Against `76d4c60` at the defaults: the same 7
  deliberate failures as `p2s1`, everything else bitwise. Read 2026-11-18
  on R 4.6.2, RcppEigen 0.3.4.0.2, RcppArmadillo 15.6.2-1.
- Capture `p2s2-eigen5`, the same source against RcppEigen 0.4.9.9-5: 71
  cases, every one bitwise identical to `p2s2`. The step adds no Eigen; run
  anyway, as for step 1.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading the step installed into a private library: 22
  files, 276 tests, none failed, errored or skipped.

Step 3, on the working tree over `905e1ac`:

- `src/eig_sym_sel.h` includes `src/si_lapack.h`, `src/pca_linalg_utils.h`,
  `<algorithm>`, `<cstddef>` and `<utility>`, nothing of Armadillo or R,
  and no longer declares `dsyevr_`. `eig_sym_sel()` calls
  `si_la::syevr_range` with the routine, options and range of the declared
  call; the wrapper's `isuppz` of `2 * topk` and its workspace, the larger
  of the query and `26 n` and `10 n`, are the lengths the chain used (READ
  in both). The wrapper asks for the workspace and allocates it on every
  call, where `EigSymWorkspace::init()` asked once and kept it. The struct
  kept nothing else, so it is now `EigSymRange`: `n`, `topk` and `il..iu`
  as `std::size_t`, where `arma::blas_int` narrowed them. Its `init()`
  checks `1 <= topk <= n` and nothing more, and the stop after it in
  `src/armaSVD.cpp` reads "internal error: eigenpair count out of range",
  where it read "eig workspace query failed". Neither is expected to be
  reachable after the argument checks, and no test raises it.
- `eig_sym_sel()` takes raw buffers and the two dimensions of `A`, which it
  checks against `n` as before, and returns the wrapper's status,
  `[[nodiscard]]`. `src/hybrid_topk_eig.h` sizes `eigvals` and `eigvecs`
  before the call, as `eig_sym_sel()` did, and keeps the failure path: a
  non-zero status is `false`, and `SVD_triplet()` stops with "top-k
  symmetric eigensolver failed to converge". `HybridEigContext` holds a
  `const EigSymRange *`.
- What changes. dsyevr now has the `n` eigenvalue cells it documents, so a
  range that cuts a group of equal eigenvalues no longer writes past
  `eigvals`. That ends the defect `dev/to-do.md` carried: its entry is
  deleted, its NEWS line moved to the list for the draft, and its demand
  for a guard-cell check moved to Phase 4 above. A Gram holding NaN or Inf
  in its upper triangle is refused before dsyevr is entered and stops with
  the message above, where the chain's call came back with a wrong count
  (the same stop) or did not return. Only the arithmetic of the fit could
  put one there, not measured, and no test does.
- `dev/review_eigen_wrappers.md` cited the line that sized the eigenvalue
  buffer; it now names that line at `79cce0f`.
- The timer build and the diagnostics build compile, with no warning from
  the chain's own files. Capture `p2s3` against `p2s2` at tolerance 0: 71
  cases, none failed, every one bitwise. Against `76d4c60` at the
  defaults: the same 7 deliberate failures as `p2s2`, everything else
  bitwise. Capture `p2s3-eigen5`, the same source against RcppEigen
  0.4.9.9-5: 71 cases, every one bitwise identical to `p2s3`. Read
  2026-11-18 on R 4.6.2, RcppEigen 0.3.4.0.2, RcppArmadillo 15.6.2-1.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading the step installed into a private library: 22
  files, 276 tests, none failed, errored or skipped.
- Timed, since the workspace query and allocation on every call, and the
  scan for NaN and Inf, are new work on the exact path: the timing driver,
  `905e1ac` as a and this step as b, both installed by
  `dev/eigen_time_install.R` from worktrees outside the tree, the four PCA
  workloads, 5 repetitions. b/a ratios of medians: `chr21_pca_exact`
  0.998, `chr21_pca_lobpcg` 1.010, `sim_pca_exact` 1.001,
  `sim_pca_lobpcg` 1.007. The 1.010 is outside the driver's band of 0.988
  to 1.008, so the two chr21 workloads were timed again, 11 repetitions:
  1.000 for LOBPCG, 1.007 for exact, both inside. No cost shows at these
  sizes (Gram of order 200 for `sim_mat()`); larger orders not timed.

## Phase 3: remove the dependency

- In `DESCRIPTION`, `RcppArmadillo` leaves `LinkingTo`.
- In `src/Makevars` and `src/Makevars.win`, `-DARMA_64BIT_WORD=1` leaves
  `PKG_CPPFLAGS`. The OpenMP flags are left as the maintainer has them.
- `arma_uword_bytes()` in `src/armaSVD.cpp`, its comment on the macro, and
  the test in `tests/testthat/test-mat_miss.R` that pins it are deleted. If
  a Phase 2 step takes Armadillo out of `src/armaSVD.cpp` first, they go in
  that step.
- `devtools::document()` regenerates `src/RcppExports.cpp` without
  `RcppArmadillo.h`.
- `.claude/CLAUDE.md`, Layout, says `src/` is Rcpp with RcppArmadillo; it
  says what `src/` is built on instead. `git grep -i armadillo` over the
  tracked notes finds anything else.
- `devtools::check()`, once against CRAN's RcppEigen and once against the
  candidate current on that day, looked up again and not taken from this
  file. Record both versions in the commit body. Local CI on the Linux
  images is the maintainer's call.

## Phase 4: off Armadillo's choices, before the release

Written down 2026-11-18 at the maintainer's request; the items are the
author's proposals, none ruled yet and none measured unless it says so.

Phases 2 and 3 reproduce Armadillo so that the conversion can be shown
bitwise. Several things the wrappers do exist only for that proof. They are
scaffolding: once Phase 3 has passed, each is removed or replaced by what
*Writing R Extensions* (`dev/r_manual/R-exts.md`) and the LAPACK
documentation give, with no regard for what Armadillo did. The package that
goes to CRAN should carry R's choices, not a second library's.

The gate changes with this phase. The oracle compares against Armadillo, so
a wrapper that leaves Armadillo's choice fails it by design. From here the
gate is the Phase 0 reference capture at its tolerance, the suite, and a
test of the wrapper's own written for the new behaviour, which must be shown
to fail. Each item is its own commit, with worst absolute and worst relative
difference against the capture and a timing from the driver.

Scaffolding to remove, in `src/si_lapack.cpp` unless said otherwise:

- `matmul`, which picks dgemv or dgemm as Armadillo's `*` does. Each call
  site names the routine it wants.
- `svd_values` asking dgesdd for its workspace only from 1024 cells up, and
  `qr_econ` sizing dorgqr from dgeqrf's answer. Every routine is asked for
  its own workspace, and the length used is the larger of the answer and
  the documented minimum. The review of the wrapper layer found no case in
  which either rule changes a result ("The review of the wrapper layer",
  what the oracle does not catch); the length does matter to dsyevr, whose
  vectors change at orders 60 and 300 when its answer is ignored.
- `qr_econ` taking the full-QR branch for square and wide input, if no call
  site passes that shape. Count the shapes the chain passes first.
- The dgelsd fallback in `solve_tri`, with dtrcon and the `approx` flag.
  Armadillo goes on with a minimum-norm answer when the triangle is singular
  to working precision. Without it a singular triangle is a status, and the
  caller in `src/lobpcg_warm.h` decides: restart, or hand over to the dense
  solver. This also deletes the workspace floors computed by hand because
  ilaenv cannot be called under the `FCONE` rule, the weakest code in the
  file. It changes results wherever the fallback fires today; how often real
  data reaches it is not measured, and is measured before this is decided.
- `syevd`, if no call site is left on it after Phase 2.
- The oracle (`dev/oracle_pairs.cpp`, `dev/oracle_run.R`,
  `dev/oracle_selftest.R`, `dev/oracle_install_armadillo.R`), their `!`
  lines, and `dev/lib/`, once the last item above is gone and nothing is
  left to compare. One check of the oracle's outlives it: guard cells after
  the eigenvalue buffer of `syevr_range`, on a range that cuts a group of
  equal eigenvalues, shown to fail when the wrapper hands dsyevr the
  caller's buffer. A test through `pca_imp()` cannot see that overrun:
  every run measured to write past the buffer returned normally. The
  section "What the wrappers reproduce" in
  `dev/pca_imp_structure.md` is rewritten to say what the wrappers do, not
  what they copy.

Choices to make in the same phase:

- Refuse NaN and Inf before every LAPACK call, not only where Armadillo
  does (the SVD and the fallback). The scan is O(n^2) against a
  factorization that is O(n^3); what a LAPACK build does with non-finite
  input is not something to depend on. Ruled for `syevr_range` and `syevd`
  on 2026-11-18, after the review measured a call that does not return,
  and done there; `qr_econ` and the exact path of `solve_tri` are still
  open, and the scan is not timed.
- Nothing compiled writes to the console. Armadillo printed its warnings to
  the error stream; with it gone, a failure is a status, raised as one cli
  error on the R side, on the main thread, never from a worker.
- `-DUSE_FC_LEN_T` in `PKG_CPPFLAGS`, which the manual allows in place of
  the define in each file. It would cover the two sampler files
  (`dev/to-do.md`, "Sampler port: what is still owed") without an edit to
  vendored code. Not compiled.
- As few routines as the chain needs. Each one dropped is one less
  declaration to depend on in the oldest supported R's `R_ext/Lapack.h`,
  and one less path through an external BLAS or LAPACK on a check machine.

To investigate, beyond the list above. Armadillo's choices were made for a
general matrix library, not for this chain under R; with it gone, each one
the code inherited is open again, as long as the answer stays inside R's
guide:

- Read how R's own sources call the same routines (the LAPACK module behind
  `eigen()`, `qr()`, `svd()` and `La.svd()`), and where R and the wrapper
  differ on options, workspace or checks, take R's form unless a
  measurement says otherwise. R's usage is what CRAN's machines exercise
  most. NOT READ yet.
- Go through every remaining routine and option in the chain and ask what
  it would be if chosen fresh: one symmetric eigensolver or two, whether
  the values-only SVD is needed at all, whether a product that is
  symmetric by construction should be dsyrk, whether the workspace of
  `syevr_range` is kept across calls, as `src/eig_sym_sel.h` kept it before
  step 3 of Phase 2.
- Go through *Writing R Extensions* once more against the finished `src/`,
  not only the sections the wrapper layer needed: registration, visibility,
  interrupts, allocation, the portability list. Record each rule and where
  it is met, as the wrapper step did.
- Each finding is a proposal with a measurement, put to the maintainer
  before it is made. Nothing here licenses a hand-rolled kernel
  (`dev/CPP.md`) or a call R's headers do not declare.

## Not settled yet

- Whether the noise floor measured on one machine holds on the BLAS of the
  rhub Linux images. Until a capture is compared there, the default
  tolerances are a reading from one BLAS.
- Whether the R callers are tidied once the return types are free to change.
- What the timings say. They are recorded here when each phase has been
  run, as readings on one machine, not as numbers to pin. Phase 0's reading
  is the driver's noise band, under "Phase 0".
