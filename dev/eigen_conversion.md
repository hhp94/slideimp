# Removing Armadillo from src/

The plan for taking RcppArmadillo out of the package's own C++. Work happens
on the `si-mi` branch, beside the work of `dev/si_mi_conversion.md`
(maintainer, 2026-11-18); there is no separate branch for it. Delete this
file when Phase 4 lands, not Phase 3 (maintainer, 2026-11-19); git history
keeps it.

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
  As carried out, Phase 2 put no Eigen in the chain: the kernel holds
  `std::vector<double>` and raw pointers, and LOBPCG holds its matrices in
  a type of its own, `lobpcg_detail::Block`. The author's reasons, neither
  measured: a raw buffer shows the order of every operation, which the
  bitwise comparison with Armadillo rests on, and the chain stays clear of
  the switch between RcppEigen versions. No step record said the plan had
  changed, which the review of Phase 2 found. Accepted (maintainer,
  2026-11-19, with the recommendations of `dev/review_p2.md` as a whole):
  the chain stays on raw buffers, and Phase 4 has the item for `Block`.
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
  in the reference capture at its tolerance. For the dot of the top
  singular values a plain loop did show there, so step 6 of Phase 2 keeps
  the order of Armadillo's own loop, by a later ruling recorded with that
  step; no wrapper, as ruled here.
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
  10 and 40. Every run of that clone returned normally, which says nothing
  of a build with the defect: the clone's buffer held the `n` cells, so
  nothing was written past it. The review of Phase 2 ran the builds as they
  are, and there the write ends the R process ("The review of Phase 2"
  below). `dev/to-do.md` carried it until step 3 of Phase 2.
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
for it. Step 3 of Phase 2 ended the overrun on this branch. That ruling was
given when the overrun had been seen to do nothing. Measured since, it is a
segfault in 1.2.0 on replicate samples ("The review of Phase 2" below), so
whether 1.2.x is patched is open again and is the maintainer's to decide.
R 4.3.0 is not checked, as before.

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
The last column is the plan as written. Where it names Eigen, Phase 2 used
raw buffers instead ("What was decided", as carried out).

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
  for a guard-cell check moved to Phase 4 above. The defect was a crash,
  not a silent write, which this step did not know: on replicate samples
  every build before this step ends the R process, and this one returns
  ("The review of Phase 2" below, which added the regression test). A Gram
  holding NaN or Inf in its upper triangle is refused before dsyevr is
  entered and stops with the message above, where the chain's call came
  back with a wrong count (the same stop) or did not return. Only the
  arithmetic of the fit could put one there. It does, on finite input
  whose squares overflow, and both ends of Phase 2 stop the same way there;
  measured and tested by the same review.
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

Step 4, the calls that can fail, on the working tree over `eade032`,
committed as `ff0e394`. The containers, products and norms of
`src/lobpcg_warm.h` stayed Armadillo's in that half; the second half, after
"Found on the way" below, replaced them.

- `src/lobpcg_warm.h` includes `src/si_lapack.h`. Its QR, its triangular
  solve and its dsyevd go through `si_la::qr_econ`, `si_la::solve_tri` and
  `si_la::syevd`. `LOBPCGState` loses the `work` and `iwork` buffers it kept
  for dsyevd, and `small_eig_desc()` and `rr_two_block()` the arguments that
  carried them.
- One function, `la_failed()`, sorts a status for the three sites. A
  positive status, which is LAPACK's `info`, and `err_nonfinite` are a
  numerical failure, and take the path the Armadillo call's failure took.
  `err_alloc` stops with "could not allocate the workspace", where a failed
  allocation was `std::bad_alloc` from Armadillo or from the buffers.
  `err_arg` and `err_size` stop with "internal error: ... refused its
  arguments", as in steps 1 to 3. No test raises either stop.
  `dev/pca_imp_structure.md` has the table of what each site does next
  (Layer 6) and what was READ ("What the wrappers reproduce").
- QR: a failure resets the block and returns rank 0, as `arma::qr_econ()`
  returning false did, and the caller's branch is the one it was.
- Triangular solve: written as the oracle's `pair_solve_chain` writes it,
  the two transposes formed and then `si_la::solve_tri('L', ...)`. A
  triangle singular to working precision is still solved by the dgelsd
  fallback and the iteration still goes on with that answer; the `approx`
  flag is not acted on. Armadillo's line on the console, "warning: solve():
  system is singular; attempting approx solution", is gone and nothing
  replaces it. A non-zero status is an R error, as Armadillo's throw was:
  "lobpcg_solve: triangular solve failed", an `Rcpp::exception`, where it
  was "solve(): solution not found", a `std::runtime_error`. The call site
  also stops if `A*Q` and `Q` differ in columns, which Armadillo's
  conformance check covered and a raw pointer does not.
- dsyevd: a non-zero status is `false`, so `LOBPCG_INITIAL_RR_FAILED` or
  `LOBPCG_RR_FAILED`, as `info != 0` was.

The gate, read 2026-11-18 on R 4.6.2, RcppEigen 0.3.4.0.2, RcppArmadillo
15.6.2-1, R's LAPACK 3.12.0:

- The timer build and the diagnostics build compile, with no warning from
  the chain's own files (511 `-Wignored-attributes` each, nothing else).
- Capture `p2s4` against `p2s3` at tolerance 0: 71 cases, none failed, every
  one bitwise. Against `76d4c60` at the defaults: the same 7 deliberate
  failures as `p2s3`, everything else bitwise. Capture `p2s4-eigen5`, the
  same source against RcppEigen 0.4.9.9-5: 71 cases, every one bitwise
  identical to `p2s4`.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading the step installed into a private library: 22
  files, 276 tests, none failed, errored or skipped.
- Timed, since the three wrappers allocate on every call and `solve_tri`
  copies its right-hand side: the timing driver, `eade032` as a and this
  step as b, both installed by `dev/eigen_time_install.R` from clones
  outside the tree, the four PCA workloads, 5 repetitions. b/a ratios of
  medians: `chr21_pca_exact` 1.008, `chr21_pca_lobpcg` 1.003,
  `sim_pca_exact` 1.000, `sim_pca_lobpcg` 1.003, all inside the driver's
  band. No cost shows at these sizes; the workloads ask LOBPCG for few
  eigenpairs, and more were not timed.
- `uv run dev/check_map_anchors.py`: the problems left are the ones in
  `dev/r_manual/R-exts.md` and the two stale citations of
  `dev/review_eigen_wrappers.md`.

What the capture does not reach. `dev/probe_p2s4_run.sh` clones the commit
before the step and the step outside the tree, has `dev/probe_p2s4_patch.R`
write one line per QR, solve, dsyevd and `lobpcg_solve()` call into each,
and runs the Phase 0 capture and `dev/eigen_p2s4_edges.R` on both;
`dev/probe_p2s4_summary.R` reads the lines. The instrumented clones'
captures are bitwise identical to `p2s3` and `p2s4`.

- In the capture, 7 of the 71 cases reach `lobpcg_solve()`: 197 calls, 185
  converged and 12 stopped at `maxiter`, asking for 3, 4 or 6 eigenpairs,
  with 5734 QRs, 2410 triangular solves and 1402 dsyevd of order 6 to 18.
  None of the three calls failed in either build, no QR dropped a column,
  no solve took the dgelsd fallback (smallest condition estimate 2.0e-9),
  and no first Rayleigh-Ritz was over `X` alone. So the capture takes none
  of the failure branches this step touched, and its bitwise result says
  nothing about them.
- `dev/eigen_p2s4_edges.R` adds 60 cases, compared between the two clones
  at tolerance 0. 20 of them inject nothing: LOBPCG with 12 to 41
  eigenpairs on `xb` and chr21, samples and columns that are copies of each
  other with more components than rank, and settings that stall LOBPCG or
  stop it early. 17 of the 20 reach `lobpcg_solve()` (the exact control
  does not, one fit converges inside the warm-up, and one input never
  yields a seed): 628 calls in each build, which in the new one end
  converged 530 times, at `maxiter` 31, stalled 33 and on the divergence
  guard 34, with 18658 QRs, 7908 solves and 4600 dsyevd of order 2 to 123.
  A QR dropped columns 36 times, all on chr21 at `ncp` 30, which takes the
  branch that recomputes `A * Q`. None of the three calls failed, no solve
  took the fallback, and no first Rayleigh-Ritz was over `X` alone; the
  copies converge in every LOBPCG call they reach. So real input took none
  of the failure branches either.
- The other 40 make a call fail on its k-th use, through an environment
  variable the instrumented clones read. A failed QR, 17 cases: while
  seeding, the state is left unseeded and the next iteration is exact; at
  the entry of `lobpcg_solve()`, the first Rayleigh-Ritz is over `X` alone
  and that call runs to `maxiter`; on the residuals,
  `LOBPCG_RESIDUAL_COLLAPSED`; on `P`, `LOBPCG_P_COLLAPSED`. A failed
  dsyevd, 9 cases: `LOBPCG_INITIAL_RR_FAILED`
  on the first, `LOBPCG_RR_FAILED` after. A triangle whose last pivot is set
  to 0, 4 cases: the old build prints Armadillo's warning 4 times, the new
  one reports `approx` 4 times, and all 4 solves then end `LOBPCG_STALLED`
  and go to the exact solver. In all 30 the two builds end every
  `lobpcg_solve()` call alike and return bitwise identical results. A
  failed solve, 2 cases: both builds raise an R error, with the two
  messages above.
- NaN, which `si_la::syevd` refuses and the old call did not, 8 cases. In a
  small Gram of order 8 to 12 the old dsyevd returned `info` of 7, 10 or 11,
  so both builds fail that Rayleigh-Ritz and agree: 3 cases bitwise, and 2
  with the NaN in the Gram itself, which the exact solver then refuses in
  both with "top-k symmetric eigensolver failed to converge". At order 34
  the old dsyevd returned `info = 0` and the builds part. With the NaN in
  the small Gram only (2 cases), the old build ends in "solve(): solution
  not found"; the new one fails that Rayleigh-Ritz, takes the exact solver
  for the iteration and returns a finite result. With the NaN in the Gram
  (1 case), the old `lobpcg_solve()` reports convergence after one
  iteration, since `rn.max()` skips NaN, and `pca_imp()` stops later with
  "the fit's mse is not finite"; the new one stops with "top-k symmetric
  eigensolver failed to converge". No case without an injection puts NaN or
  Inf in a Gram, in the capture or in the edges.

What step 4 changes that the capture cannot see: the workspace dsyevd is
given. `small_eig_desc()` kept the longest workspace so far on the state,
and 1.2.0 does (READ at the tag); `si_la::syevd` gives each problem the
length of its own query. `dev/pca_imp_structure.md`, "The two eigen
wrappers against the chain's calls", has why the length matters: dsyevd's
eigenvectors, not its eigenvalues, depend on it for orders 34 to 79
(`dev/probe_p2s4_syevd_lwork.R`: 92 of 92 calls of those orders differ by
up to 6.9e-16 when the workspace is that of a problem twice the order, 0 of
186 at other orders, and the control with the same length is bitwise, 278
of 278). The oracle's cases with a longer workspace are of order 30, 90 and
200, on either side of that range, which is why it reports the pair as the
same.

- Through `pca_imp()` that is LOBPCG asked for 17 eigenpairs or more: `ncp`
  16 when regularized, 17 with EM. Of the 9 edge cases with many
  eigenpairs, 5 are not bitwise between the builds: `xb` at `ncp` 16 and 25,
  `xb` with EM at 17, chr21 at 20 and 30. Every imputed cell of those
  results differs and no other does. Worst over a result: 6.2e-8 absolute
  and 1.0e-6 relative (`xb` at `ncp` 25, the relative figure over the
  119634 of its 120000 cells that hold a nonzero value). The compare
  script's own worst relative figure, 1.1e-4, is on the `criterion_final`
  attribute. The fits take the same number of iterations in both builds,
  while individual `lobpcg_solve()` calls end differently, so in 3 of the
  5 the `n_exact`, `n_lobpcg_ok` and `n_lobpcg_bad` attributes differ. For
  scale, on that input the exact solver and LOBPCG differ by 1.9e-7 and
  3.4e-6 in the old build, and `threshold` is 1e-6. The 4 that stay
  bitwise: `ncp` 11 on both inputs (12 eigenpairs: orders 34 to 36 occur,
  but the longest workspace held is too short to change dormql's
  blocking), `xb` at `ncp` 40, and the exact control.
- It is the only difference between the builds. With
  `SLIDEIMP_P2S4_FRESH` set, the old clone gives dsyevd the length of its
  own query on every call and is otherwise the code of `eade032`: of its 60
  edge cases against the new clone's, at tolerance 0, the 53 that return a
  value in both builds are bitwise identical and 2 raise the same error in
  both. The 5 that fail are the 2 failed solves and the 3 NaN cases at
  order 34 above, which differ in the error raised and not in a value.
  Without the variable 10 fail: those 5, and the 5 cases above.
- Ruled (maintainer, 2026-11-18): the wrapper's length stays. The buffers
  were kept to save an allocation, which the timing above does not show,
  and their length was never chosen for its effect on the result; with the
  length of each problem's own query the eigenvectors depend on the matrix
  alone. So results for 17 eigenpairs or more move from step 3's, and from
  1.2.0's, by the amounts above, and `dev/to-do.md` has a line for the NEWS
  draft. Keeping step 3's bits would have taken a wrapper argument for the
  workspace length and the largest order so far on the state, scaffolding
  for Phase 4 to remove.

Found on the way, and not changed:

- The oracle never runs `si_la::syevd` against the old call at an order
  between 34 and 79 with a longer workspace. With such a case its class
  "same" would fail, by the measurement above. After the ruling, the cases
  it does run with a longer workspace compare the wrapper with a call the
  chain no longer makes.
- `lobpcg_solve()` takes its largest residual norm with `rn.max()`, which
  returns `-Inf` for a vector of NaN (READ: `op_max::direct_max()` compares
  with `>` from the most negative value; NOT READ: that `Mat::max()`
  reaches it, which the injected case measures), so a residual that is all
  NaN passes the convergence test. Since this step no route to it is known: a
  Gram or a small Gram holding NaN is refused at the Rayleigh-Ritz before
  the loop.
- For Phase 4, on the dgelsd fallback: no input here reached it without an
  injection, in 2410 solves of the capture and 7908 of the edge cases
  that inject nothing; the smallest condition estimate among them is
  1.3e-14 (`xb` at `ncp` 40), against the `DBL_EPSILON` of 2.2e-16 that
  triggers it. In the 4 cases that force it, every solve that took the
  fallback then stalled.

Step 4, the containers and products, on the working tree over `ff0e394`:

- `src/lobpcg_warm.h` holds its matrices in `lobpcg_detail::Block`: a
  `std::vector<double>`, column-major and unpadded, with `std::size_t`
  dimensions, zero-filled when sized. Its members carry Armadillo's names
  (`n_rows`, `n_cols`, `memptr()`, `(i, j)`), so `la_failed()`, the three
  `si_la` call sites and every failure branch of the first half read as
  they did, unchanged to the character. Every product goes through
  `si_la::matmul`, from one function, `matmul_into()`; a non-zero status
  stops with "internal error: dgemm or dgemv refused its arguments", not
  expected to be reachable, as in steps 1 to 3. No wrapper was added.
- Each Armadillo expression is written the way Armadillo evaluated it
  (`dev/pca_imp_structure.md`, "Products of two `arma::mat`", has what was
  READ). `-=` of a product whose operands do not include its target is one
  call with `alpha = -1` and `beta = 1` (`sub_mul()`): the four updates of
  `Pact` and `APact` between their two QRs. When the target is inside the
  product, as in `V -= X * (X.t() * V)` at the three projections onto the
  complement of `X`, the product is formed and subtracted cell by cell
  (READ: `Glue::is_alias()` asks the nested product too). A product inside
  a sum or difference is formed first and added cell by cell, left to right:
  `build_momentum_P()`, the projections of `P` and `AP`, `pp`, `app`, and
  `X * VX + pp`. `M = 0.5 * (M + M.t())` reads every cell of `M` as it was
  (READ: the transpose's proxy has `has_subview`, so the assignment goes
  through a temporary). `find`, `cols()` of an index vector, `fliplr`,
  `reverse`, `head` and the sub-block copies are loops that move values.
- The norms (`dev/pca_imp_structure.md`, "Calls no wrapper stands in for").
  The infinity norm of `A` is bitwise Armadillo's: row sums in column
  order, then the largest taken with `>` from `-Inf`, as READ in `op_sum`
  and `op_max`. The Frobenius norm is the plain sum of squares the ruling
  above asked for. The ruling named three sites; this is a fourth: the
  residual norms, `sqrt(sum(square(R), 0))`, which Armadillo took through
  `op_dot::direct_dot()`, two partial sums up to 32 rows and ddot above
  (READ). Treated as the Frobenius norm, a plain sum of squares and no
  wrapper, so its last bits are not ddot's. The norms decide convergence,
  the active columns, the stall count and the divergence guard, and are
  `max_rel_res`, which only the diagnostics build returns
  (`eig_lobpcg_max_rel_res`). A last-bit change moves a decision only when
  a norm lies within rounding of the tolerance or of the best so far; in
  the capture and the 60 edge cases none did (below). `rn.max()` reached
  `op_max::direct_max()` (READ: `Mat::max()`, which the first half left NOT
  READ), and `max_skip_nan()` does what it did, so a residual of NaN still
  passes the convergence test as "Found on the way" above says.
- What changes, none of it reachable from `pca_imp()`. A block with more
  columns than rows stops in `qr_orthonormalize()` with "internal error: a
  LOBPCG block has more columns than rows", where the loop over the
  diagonal of `R` read past the cells `R` has (NOT READ: whether
  Armadillo's bounds check was on); every block has at most `k` columns and
  `n` rows, and `k <= n` holds since the state comes from `eig_sym_sel()`.
  Blocks of different sizes added, or a sub-block out of range, stop with an
  internal error where Armadillo's conformance check threw. The QR's `Q` and
  `R` are zero-filled where `arma::fill::none` left them; the wrapper writes
  every cell before anything reads one.
- What changes and is reachable. Found by the review of Phase 2; this
  record did not have it until then. A product whose matrix operand is
  square and 4 x 4 or less goes to BLAS through `si_la::matmul`, where
  Armadillo ran loops of its own (`dev/pca_imp_structure.md`, "Two things
  Armadillo does that `si_la::matmul` does not"). In the chain that is a
  Gram of order 2 to 4 under `solver = "lobpcg"` with one eigenpair asked
  for, or every one. Measured, `59920f4` against `597d28e`
  (`dev/review_p2_cases.R`, groups `tiny` and `tinysweep`): of 280 such fits
  at the default controls 267 differ, by at most 5.2e-9 absolute and 1.0e-9
  relative on the result's cells, with `n_iter` and `converged` the same in
  each; where every eigenpair is asked for, `n_exact` and `n_lobpcg_ok`
  differ in 71 of 80 fits and `n_lobpcg_bad` in 77. The 80 fits outside
  the rule and every fit under the exact solver are identical. The gate
  below is bitwise because neither the capture nor the edge cases hold a
  Gram that small. The default solver takes the exact path below order 250,
  so the caller has to name LOBPCG. Accepted (maintainer, 2026-11-19, with
  the recommendations of `dev/review_p2.md` as a whole): it stays, with a
  line in `dev/to-do.md` for the NEWS draft. Phase 4 takes Armadillo's
  choice of routine out of `matmul` and leaves these products on BLAS.
- Armadillo left in the file, all of it at the boundary with
  `src/hybrid_topk_eig.h`, which step 5 converts: `#include
  <RcppArmadillo.h>`; `LOBPCGState`'s `X` and `P` as `arma::mat`, which
  `hybrid_topk_eig_impl()` copies into `eigvecs`, resets in
  `choose_auto_exact()` and tests with `seeded()`; `LOBPCGResult::eigvals`
  as an `arma::vec`, which it moves into its own; the Gram as a
  `const arma::mat &`, read in place through `memptr()`; `X_curr` and
  `X_prev` of `seed_lobpcg_state()` and `refresh_momentum_P()` as
  `arma::mat`, which is how the caller holds them; and `from_arma()` and
  `to_arma()`, which copy each of them into a `Block` on the way in and out
  of one on the way out. Step 5 changes those types and removes the copies.
  `dev/probe_p2s4_patch.R` writes `ctx.lobpcg_state.X.n_cols` into
  `src/hybrid_topk_eig.h`, which step 5 keeps or updates.
- The probe, to compare this half with `ff0e394`. `dev/probe_p2s4_patch.R`
  chose its patch by mode, and its "old" patch was for the file of
  `eade032`. It now chooses by what the file holds: the calls through
  Armadillo, or the calls through `si_la` with either container, the two
  anchors that differ (the recomputed `A * Q` and the transpose of `R`)
  taken from whichever is present. Applied to the file of `eade032`,
  `ff0e394` and this tree, each patches cleanly. `dev/probe_p2s4_run.sh`
  makes the `edges-fresh` run only when the old file has the hook for it,
  and `dev/probe_p2s4_summary.R` reads each log by its own format, so
  either build may be on either side.

The gate, read 2026-11-18 on R 4.6.2, RcppEigen 0.3.4.0.2, RcppArmadillo
15.6.2-1, R's LAPACK 3.12.0:

- The timer build and the diagnostics build compile, with 511
  `-Wignored-attributes` lines each and no other compiler warning. R
  notes the incomplete final line of `~/.R/Makevars.win`, a file of the
  machine's, not of the package.
- Capture `p2s4b` against `p2s4` at tolerance 0: 71 cases, none failed,
  every one bitwise. Against `76d4c60` at the defaults: the same 7
  deliberate failures as `p2s4` (`col_vars_x2_cores2`, `int_col_miss`,
  `int_col_vars`, `int_mean_imp_col`, `int_mean_imp_col_cores2`,
  `int_mean_imp_col_inf_error`, `int_row_miss`), everything else bitwise.
  Capture `p2s4b-eigen5`, the same source against RcppEigen 0.4.9.9-5: 71
  cases, every one bitwise identical to `p2s4b`.
- `bash dev/probe_p2s4_run.sh <scratch> ff0e394`, `ff0e394` as old and this
  tree as new. Each instrumented clone's capture is bitwise identical to
  `p2s3` and `p2s4`. The 60 edge cases at tolerance 0: none failed, every
  one bitwise, and the 5 that end in an error raise the same class and
  message in both, 3 "top-k symmetric eigensolver failed to converge" and 2
  "lobpcg_solve: triangular solve failed". The two builds' event logs are
  identical byte for byte, 75097 lines for the edge cases and 9743 for the
  capture. So in the edge cases both make the same 1207 `lobpcg_solve()`
  calls, asking for 1 to 41 eigenpairs, and end each alike: 1074 converged,
  37 stalled, 35 at `maxiter`, 34 on the divergence guard, 11
  `LOBPCG_RR_FAILED`, 6 `LOBPCG_INITIAL_RR_FAILED`, 6
  `LOBPCG_RESIDUAL_COLLAPSED` and 4 `LOBPCG_P_COLLAPSED`. They make 43966
  QRs, 17 failed by injection and 36 dropping columns; 19031 triangular
  solves, 4 by the dgelsd fallback, all injected; and 10743 dsyevd calls of
  order 2 to 123, 17 failed: 9 by an injected `info` and 8 refused for an
  injected NaN.
  2 first Rayleigh-Ritz problems were over `X` alone. Every failure branch
  of the file is taken and comes out as in `ff0e394`, and one eigenpair
  takes the dgemv routes.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading this tree installed into a private library: 22
  files, 276 tests, none failed, errored or skipped.
- Timed, since every block is now allocated, zero-filled and, at the
  boundary, copied: the timing driver, `ff0e394` as a and this tree as b,
  both installed by `dev/eigen_time_install.R` from clones outside the
  tree, the four PCA workloads, 5 repetitions. b/a ratios of medians:
  `chr21_pca_exact` 1.006, `chr21_pca_lobpcg` 1.023, `sim_pca_exact` 1.005,
  `sim_pca_lobpcg` 1.013. The two LOBPCG workloads are outside the driver's
  band of 0.988 to 1.008, so they were timed again, 11 repetitions:
  `chr21_pca_lobpcg` 1.004, inside, and `sim_pca_lobpcg` 1.015, with the
  ranges of its runs all but apart (a 0.483 to 0.491 s, b 0.490 to
  0.501 s). So the LOBPCG path costs about 1.5% more on `sim_mat()` input
  (Gram of order 200) and nothing measurable on chr21. Where it goes is not
  profiled; the new work is the zero-fill, the copies of the state and the
  sub-blocks and transposes written out.
- `uv run dev/check_map_anchors.py`: the problems left are the ones in
  `dev/r_manual/R-exts.md` and the two stale citations of
  `dev/review_eigen_wrappers.md`, 82 in all. The ASCII check of
  `.claude/CLAUDE.md` names one file, `.github/workflows/pkgdown.yaml`,
  whose emoji came in with `0b1e87c` and which this step does not touch.

Step 5, on the working tree over `99a8c5e`:

- `src/hybrid_topk_eig.h` and `src/svd_triplet.h` include `<Rcpp.h>`, for
  `Rcpp::stop`, and nothing of Armadillo; `src/lobpcg_warm.h` the same, so
  the boundary step 4 left there is gone. `src/armaSVD.cpp`, which keeps its
  `arma::mat` until step 6, now includes `RcppArmadillo.h` itself, ahead of
  `src/svd_triplet.h`: it had it only through that header, and
  RcppArmadillo stops the build when `Rcpp.h` came first (READ:
  `RcppArmadillo.h`, the `#error` on `Rcpp_hpp`).
- `lobpcg_detail::Block` moves ahead of `LOBPCGState`, which holds `X` and
  `P` as `Block`s, and gains `is_empty()`, so `seeded()` and
  `hybrid_topk_eig_impl()` read as they did. `LOBPCGResult` holds the
  eigenvalues as a `std::vector<double>`. `lobpcg_solve()` takes the Gram
  as a raw buffer and its two dimensions, which it checks as before.
  `seed_lobpcg_state()` and `refresh_momentum_P()` take `Block`s.
  `from_arma()` and `to_arma()` are gone, and with them the copies step 4
  made at the boundary: `lobpcg_solve()` copies the state into its working
  blocks and moves them back on success, and `hybrid_topk_eig_impl()` copies
  each solve's eigenvectors into one `Block`, which seeds or refreshes the
  state and then becomes `X_prev`.
- `hybrid_topk_eig()` takes raw buffers: the Gram and its two dimensions,
  and the `topk` eigenvalue and `n x topk` eigenvector cells of the range.
  On a LOBPCG success it copies into them, after a check that the state
  still has that shape (it keeps the shape of the block it was seeded from)
  and stops with "internal error: the LOBPCG state does not match the
  eigenpair range" otherwise; not expected to be reachable, and no test
  raises it. Before, the result was resized to whatever the state held.
- `SVD_triplet()` takes raw buffers and reads every dimension from
  `GramWorkspace` and the eigenpair range. It no longer sizes its outputs:
  `src/armaSVD.cpp` sizes `U`, `V`, `Vn`, `eigvecs`, `eigvals`, `vs_top` and
  `d_inv` once, before the EM loop, with the dimensions `SVD_triplet()` gave
  them on every call, so `V`'s fixed rows stay unwritten in the wide case as
  before. `X_work`, `GramWorkspace` and `GramCache` are taken `const`: it
  never wrote them.
- No wrapper was added. The one Armadillo call, `arma::trace()`, is
  `gram_trace()`. READ (`fn_trace.hpp`): two partial sums, over the even and
  the odd diagonal cells, added at the end. That is the multi-accumulator
  reduction `dev/CPP.md` bans, so it was measured both ways:
  - In Armadillo's order (capture `p2s5a`, and the tree as it stands): every
    compare below is bitwise.
  - As a plain sum in order (capture `p2s5-plaintrace`, all else the same),
    against `p2s5a` at tolerance 0: 28 of the 71 cases not bitwise, worst on
    the imputed values 6.5e-10 absolute and 3.1e-8 relative
    (`pca_x1_lobpcg`, over the 17733 of its cells that hold a nonzero
    value), worst on `criterion_final` 9.8e-7 relative. Against `p2s4b` at
    the defaults, 22 fail. The probe's 60 edge cases against `99a8c5e`: 52
    not bitwise; worst on the imputed values 6.5e-8 absolute and 1.1e-6
    relative (`big_xb_ncp25`, over its 119634 nonzero cells). The event
    logs part: in 41 of the 55 cases that reach `lobpcg_solve()` some call
    ends after another number of iterations, mostly one more or one fewer,
    and in some with another status (`big_xb_ncp40` and `ctl_x1_tight`:
    converged in one build, at `maxiter` in the other); in 6
    (`big_xb_ncp16`, `25` and `40`, `ctl_x1_tight`, `inj_xb16_syevdnan_01`
    and `02`) that changes which solver an outer iteration used, so
    `n_exact`, `n_lobpcg_ok` and `n_lobpcg_bad` differ; and
    `criterion_final` moves by 0.12 in `dup_xd_ncp12` and is NaN in one build
    and finite in the other in `dup_xd_ncp07_em`, both inputs with copied
    columns, where the criterion is a ratio of two objectives near zero.
    The trace feeds the tail energy, so `sigma2`, the shrinkage and the
    objective, each iteration. Not chased further: a move beyond the noise
    floor is Fable's ("Who does which work").
  - Ruled (maintainer, 2026-11-18): Armadillo's order stays, so the step is
    bitwise. It is scaffolding against `dev/CPP.md`, listed for Phase 4.
- `dev/probe_p2s4_run.sh` copied only `src/lobpcg_warm.h` into its "new"
  clone; it now copies every file under `src/` the working tree changes
  against HEAD. `dev/probe_p2s4_patch.R` finds the `lobpcg_solve()` call in
  `src/hybrid_topk_eig.h` by its form, the Gram an `arma::mat` before this
  step or a raw buffer from it, and writes the injected NaN into `A(0, 0)` or
  `A[0]` to match; applied to the files of `eade032`, `99a8c5e` and this
  tree, each patches cleanly.

The gate, on the tree as it stands (the trace in Armadillo's order), read
2026-11-18 on R 4.6.2, RcppEigen 0.3.4.0.2, RcppArmadillo 15.6.2-1, R's
LAPACK 3.12.0:

- The timer build and the diagnostics build compile, with 511
  `-Wignored-attributes` lines each and no other compiler warning.
- Capture `p2s5` against `p2s4b` at tolerance 0: 71 cases, none failed,
  every one bitwise, and the same against `p2s5a`, which differs from this
  tree only in a comment. Against `76d4c60` at the defaults: the same 7
  deliberate failures as `p2s4b`, everything else bitwise. Capture
  `p2s5-eigen5`, the same source against RcppEigen 0.4.9.9-5: 71 cases,
  every one bitwise identical to `p2s5`.
- `bash dev/probe_p2s4_run.sh <scratch> 99a8c5e`, `99a8c5e` as old and this
  tree as new. Each instrumented clone's capture is bitwise identical to
  `p2s3` and `p2s4`; the 60 edge cases at tolerance 0: none failed, every
  one bitwise, and the 5 that end in an error raise the same class and
  message in both. The two builds' event logs are identical byte for byte,
  75097 lines for the edge cases and 9743 for the capture, with the counts
  step 4 gives: the same 1207 `lobpcg_solve()` calls ending the same way,
  and every failure branch of `src/lobpcg_warm.h` taken by injection. The
  script exits 1 only because its last line tests for the `edges-fresh`
  run, which `99a8c5e` has no hook for.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading this tree installed into a private library: 22
  files, 276 tests, none failed, errored or skipped.
- Timed: the timing driver, `99a8c5e` as a and this tree as b, both
  installed by `dev/eigen_time_install.R` from clones outside the tree, the
  four PCA workloads, 5 repetitions. b/a ratios of medians:
  `chr21_pca_exact` 0.978, `chr21_pca_lobpcg` 0.990, `sim_pca_exact` 0.994,
  `sim_pca_lobpcg` 1.003. The first is outside the driver's band of 0.988 to
  1.008, on the fast side; the two chr21 workloads timed again, 11
  repetitions: 0.996 and 0.990, both inside. No cost and no gain shows; in
  particular the boundary copies this step removes do not account for the
  cost step 4 measured on `sim_pca_lobpcg` (Phase 4, "To investigate").
- `uv run dev/check_map_anchors.py`: 82 problems, the known ones in
  `dev/r_manual/R-exts.md` and `dev/review_eigen_wrappers.md`. The ASCII
  check of `.claude/CLAUDE.md` names only `.github/workflows/pkgdown.yaml`,
  as in step 4.

Step 6, on the working tree over `8c15416`. One session did both shares of
"Who does which work": the file's one call that can fail is in its
diagnostics build, and the rest is the substitution around it.

- `src/armaSVD.cpp` includes `Rcpp.h`, `src/svd_triplet.h`,
  `src/matrix_checks.h` and the standard library, nothing of Armadillo. No
  file under `src/` uses Armadillo now; `src/RcppExports.cpp`, generated,
  still includes `RcppArmadillo.h`, since `DESCRIPTION` links to it until
  Phase 3. `arma_uword_bytes()`, its comment and the test in
  `tests/testthat/test-mat_miss.R` that pinned `-DARMA_64BIT_WORD=1` are
  deleted, as Phase 3 asked of whichever step came first.
- The export, under the type rules. `pca_imp_internal_cpp()` takes `obj` as
  a `const Rcpp::NumericMatrix &`, read in place; `eligible_idx` as an
  integer vector of 0-based columns; `ncp`, `init`, `maxiter`, `miniter`,
  `warmup_iters`, `lobpcg_maxiter` and `trace_iter` as `int`; `row_w` as a
  `const Rcpp::NumericVector &`, copied once into the vector the kernel
  normalizes. `pca_imp()` converts `ncp`, `maxiter` and `miniter` at the
  call, since `checkmate::assert_int()` leaves their type alone and their
  defaults are doubles; the others it already holds as integers. It
  returns the imputed cells as `loc`, an integer matrix with columns `row`
  and `col`, 1-based in `obj`, and `value`, a double vector, where it
  returned `imputed_values`, an n x 3 double matrix holding the indices as
  doubles: the form `impute_knn_brute()` took in Phase 1. `pca_imp()` hands
  both to `new_part()` as they are; its `cbind()` and `as.integer()` of the
  two index columns go, and `clamp` is applied to `value`. The diagnostics
  build returns `obj_hist` and `subspace_cos_hist` as vectors, where
  Armadillo's `vec` came back as an n x 1 matrix.
- Ruled (maintainer, 2026-11-18, during the step): `n_iter` is an integer,
  in the kernel and in what `pca_imp()` returns, and so is any other count
  that should be one. The `n_iter` attribute of the matrix output, and
  `info$n_iter` and `prov$iter` of the `mipca` output, were doubles holding
  a whole number; `pca_imp_info()` defaults it to `0L`, and three tests that
  pinned `0` pin `0L`. No other count crossed as a double: `n_exact`,
  `n_lobpcg_ok`, `n_lobpcg_bad` and `solver_chosen` were `int` already.
  `dev/to-do.md` has the line for the NEWS draft.
- What the kernel refuses. A count below its range, `NA_integer_` included
  (it is `INT_MIN`): `maxiter < 1`, `miniter < 0`, `ncp < 1`, `init < 0`,
  `trace_iter < 0`, and under LOBPCG `lobpcg_maxiter < 1` and
  `warmup_iters < 0`. As `arma::uword` a negative count wrapped to a huge
  one, which was refused by another check or not at all (`init`,
  `trace_iter`). A column index outside `obj`, by `stop_on_bad_col_index()`
  of `src/matrix_checks.h`, in the wording the K-NN kernel and the sampler
  share, where it read "eligible_idx contains out-of-bounds column
  indices". More missing cells than `INT_MAX`, since `loc` has a row for
  each (not reachable in a test). Gone, because an `int` cannot fail them:
  the comparisons of `nrX`, `n_elig`, `ncp`, `k` and `n_gram` with the
  `blas_int` range, and "LOBPCG iteration controls exceed int range". The
  two workspace guards stay. `pca_imp()` validates every one of these
  before the kernel sees it, so none changes an error a user gets; a new
  test in `tests/testthat/test-pca_imp.R` calls the kernel directly, pins
  the types of `loc`, `value` and `n_iter` against `which(is.na(x))`, and
  raises each refusal but the one for `INT_MAX` cells.
- Storage. Every matrix and vector is a `std::vector<double>`, column-major
  and unpadded, zero-filled when sized, where `arma::mat` left `Xhat`,
  `X_work`, `AA_NxN`, `U`, `V`, `Vn` and the eigenpair buffers as allocated
  (`fittedX` was zero-filled before too). A count is held in `std::size_t`
  once checked, so `j * nrX` is computed in `std::size_t`; row indices are
  `int`, so `miss_rows_flat` takes 4 bytes a missing cell where it took 8;
  offsets into it are `std::size_t`. `initialize_missing_start()` no longer
  allocates `fittedX`. `find`, `elem`, `cumsum`, `head`, `tail`, the
  element-wise expressions behind the column means and scales, and the
  tests of `row_w` are loops; `dev/pca_imp_structure.md`, "Calls no wrapper
  stands in for", says which were READ.
- BLAS. `reconstruct_cols()` calls `si_la::gemm` with the arguments it
  passed to `arma::dgemm_`, and a non-zero status stops with "internal
  error: dgemm refused its arguments", as in steps 1 to 3.
- The call that can fail: in the diagnostics build, `arma::svd(V_prev.t() *
  V_new)`, once per iteration from the second on. The product is formed by
  `si_la::matmul`, which sends it where `*` did, and its singular values
  taken by `si_la::svd_values`. A positive status or `err_nonfinite` stops
  with "pca_imp diagnostics: the SVD of the subspace product failed", an
  `Rcpp::exception`, where `arma::svd()` threw `std::runtime_error`, "svd():
  decomposition failed": an R error both ways. `err_alloc` stops with
  "could not allocate the workspace of dgesdd", and a refusal of the
  arguments with an internal error. No capture or test reaches the failure,
  and no route to it is known: `si_la::syevr_range` refuses a Gram holding
  NaN before the eigenvectors exist. `dev/probe_p2s6_run.sh` forces it
  (below). No wrapper was added.
- Two sums, kept in Armadillo's order. Both are the reduction `dev/CPP.md`
  bans, two partial sums over the even and the odd positions added at the
  end, and both were measured as a plain sum in order, in scratch clones
  that differ from this tree in that one loop:
  - `arma::accu(row_w)`, the sum the row weights are divided by, is
    `sum_two_partial()`. As a plain sum the capture is bitwise identical,
    all 71 cases, and so are the 20 edge cases of `dev/eigen_p2s4_edges.R`
    that inject nothing: `pca_imp()` passes equal weights by default, which
    sum exactly in either order, and the capture's two `row.w` cases happen
    not to move. With unequal weights it moves: of the 5 such cases of
    `dev/eigen_p2s6_diag.R`, 4 are not bitwise and 2 fail the default
    tolerance. Worst on the imputed values 3.9e-12 absolute and 5.3e-11
    relative (`xw_ncp2_roww_runif_lobpcg`, over its 720 cells); under the
    exact solver 1.3e-15 and 2.2e-13 (the two `x1` cases, over 1800 cells
    each).
  - `arma::dot()` of the leading singular values with themselves, which the
    tail energy subtracts from the trace, is `sum_squares_two_partial()`.
    "Not wrapped" above ruled a plain loop here, before step 5 measured
    the trace. As a plain sum in order, against `p2s6` at tolerance 0: 15 of
    the 71 cases not bitwise, 10 failing the default tolerance, most on
    `criterion_final` (worst 9.9e-8 relative); worst on the imputed values
    1.7e-11 absolute and 1.1e-10 relative (`pca_x1_lobpcg`, over the 17733
    of its cells that hold a nonzero value). On the 20 edge cases that
    inject nothing: 14 not bitwise, worst on the imputed values 6.5e-8
    absolute and 1.1e-6 relative (`big_xb_ncp25`, over its 119634 nonzero
    cells); `n_exact`, `n_lobpcg_ok` and `n_lobpcg_bad` differ in 4, and
    `criterion_final` is NaN in one build and finite in the other in 4 of
    the cases with copied samples. The same picture as the trace, which
    is the other term of the same difference. So the order of Armadillo's
    own loop is kept, on the footing of `gram_trace()`, and listed with it
    for Phase 4. Ruled (maintainer, 2026-11-18, with the figures above):
    Armadillo's order stays, in place of the plain loop "Not wrapped"
    names. The row-weight sum was put to the maintainer in the same report
    and stays in Armadillo's order on the same footing. That was the
    author's reading, with no word given on that sum by name, as the review
    of Phase 2 found. Accepted since (maintainer, 2026-11-19, with the
    recommendations of `dev/review_p2.md` as a whole): `sum_two_partial()`
    stays through Phases 2 and 3 and is on Phase 4's list with the other
    two.
- What is not bitwise: more than 32 components. Armadillo's `dot()` ran its
  own loop up to 32 values and gave a longer vector to ddot (READ:
  `op_dot::direct_dot()`). "Not wrapped" rules ddot out, so
  `sum_squares_two_partial()` sums every length the same way, and a fit
  with `ncp` above 32 no longer matches step 5 to the bit. Of the edge
  cases, `big_xb_ncp40` is the one with that many: imputed values differ by
  8.1e-9 absolute and 2.8e-7 relative over its 119634 nonzero cells,
  `criterion_final` by 6.6e-7 relative, and the three solver counts differ.
  Of the diagnostics cases, the 4 with `ncp` of 33 and 40: under the exact
  solver the imputed cells differ by at most 8.2e-15 absolute and 9.7e-14
  relative, under LOBPCG at `ncp` 40 by 4.0e-8 and 3.7e-7. It is the only
  cause: a scratch clone that calls ddot above 32 values, and is this tree
  otherwise, is bitwise identical to step 5 on all 9 many-eigenpair edge
  cases and on every diagnostics case that returns a value. For scale, the
  dsyevd workspace of step 4 moved such fits by 6.2e-8. `dev/to-do.md` has
  a line for the NEWS draft.

The gate, read 2026-11-18 on R 4.6.2, RcppEigen 0.3.4.0.2, RcppArmadillo
15.6.2-1 (still linked, used by no file), R's LAPACK 3.12.0:

- The timer build and the diagnostics build compile, with 511
  `-Wignored-attributes` lines each and no other compiler warning.
- Capture `p2s6` against `p2s5` at tolerance 0: 71 cases, every one bitwise,
  and 24 failed on purpose: each `pca_imp()` case, on the type of its
  `n_iter` attribute and nothing else. `dev/probe_p2s6_niter.R` reads the
  two captures: of the 69 cases with a value in both, 24 carry `n_iter`, a
  double before and an integer after in all 24, the counts equal in all 24,
  and the two results `identical()` in all 24 once the count is read as a
  double on both sides. Against `76d4c60` at the defaults: those 24 and the
  7 deliberate failures of Phase 1, everything else bitwise. Capture
  `p2s6-eigen5`, the same source against RcppEigen 0.4.9.9-5: 71 cases,
  none failed, every one bitwise identical to `p2s6`. The step adds no
  Eigen; run anyway.
- `bash dev/probe_p2s4_run.sh <scratch> 8c15416`, step 5 as old and this
  tree as new; the script now copies the changed files under `R/` into its
  new clone as well as those under `src/`. Each instrumented clone's
  capture is bitwise identical to the capture of the build it stands for,
  the new one's but for `n_iter`. The 60 edge cases at tolerance 0: the 5
  that end in an error raise the same class and message in both; of the 55
  that return a value, 54 are `identical()` once `n_iter` is read as a
  double, and the one left is `big_xb_ncp40`, above. The two builds' event
  logs are identical line for line outside that case, 69674 lines each, so
  every failure branch of `src/lobpcg_warm.h` is still taken by injection
  and comes out as in step 5; within it the old build logs 5423 events and
  the new 5752.
- `bash dev/probe_p2s6_run.sh <scratch> 8c15416`: the kernel called
  directly in the diagnostics build, on a clone of each side instrumented
  by `dev/probe_p2s6_patch.R`, 49 cases of `dev/eigen_p2s6_diag.R`, each
  result reduced to a form both kernels give. Of the 44 that return a
  value in both, 40 are bitwise identical at tolerance 0, histories and
  solver logs included: wide and tall inputs, with and without complete
  columns, 1 to 32 components (1 takes the dgemv route of the product),
  both solvers, no scaling, the EM method, a random start, unequal row
  weights. The 4 that are not are the ones with 33 and 40 components,
  above. The SVD made to fail: with a NaN written into the product (3
  cases), both builds stop, the old with "svd(): decomposition failed" and
  the new with its own message; with the status of `si_la::svd_values` set
  to 1 (2 cases), the new build stops the same way, and the old has no
  such way in, so it runs on. A case whose injection names an SVD the call
  never reaches is bitwise identical.
- The suite with the manual tests, under `load_all1(timer = FALSE)`, the
  mirai daemons loading this tree installed into a private library: 22
  files, 276 tests, none failed, errored or skipped. One test left
  (`arma_uword_bytes()`) and one came (the kernel's types and refusals).
- Timed, since every buffer is now zero-filled and the return changed
  shape: the timing driver, `8c15416` as a and this tree as b, both
  installed by `dev/eigen_time_install.R` from clones outside the tree, the
  four PCA workloads. b/a ratios of medians over 5 repetitions:
  `chr21_pca_exact` 0.996, `chr21_pca_lobpcg` 0.981, `sim_pca_exact` 0.951,
  `sim_pca_lobpcg` 0.951; over 11: 0.999, 0.992, 0.934 and 0.956. The two
  chr21 workloads are inside the driver's band of 0.988 to 1.008 on the
  second run. The two `sim_mat()` workloads are faster, with the ranges of
  their runs apart (`sim_pca_exact`: a 0.385 to 0.391 s, b 0.358 to
  0.372 s; `sim_pca_lobpcg`: a 0.506 to 0.515 s, b 0.482 to 0.496 s). Where
  the gain comes from is not profiled. A guess, not measured: `pca_imp()`
  no longer builds an integer map and a value vector out of an n x 3 double
  matrix, three copies of the imputed cells, and those inputs have many
  such cells for their size.
- `uv run dev/check_map_anchors.py`: 82 problems, the known ones in
  `dev/r_manual/R-exts.md` and `dev/review_eigen_wrappers.md`; a third
  citation of that note, and one of `dev/si_mi_conversion.md`, named lines
  this step rewrote and now name what replaced them. The ASCII check of
  `.claude/CLAUDE.md` names only `.github/workflows/pkgdown.yaml`, as in
  step 4.

## The review of Phase 2

The fourth checkpoint of "Review checkpoints", on `59920f4..597d28e`, run
2026-11-18 and 19. `dev/review_p2.md` has the findings; its probe scripts
are `dev/review_p2_run.sh`, `dev/review_p2_cases.R`,
`dev/review_p2_compare.R`, `dev/review_p2_show.R`, `dev/review_p2_tag.sh`,
`dev/review_p2_suite.sh` with `dev/review_p2_suite.R`, and
`dev/review_p2_crash_events.ps1`. Not a fresh session: the maintainer ran
`/code-review` in the session that wrote step 6 and the failure-path half
of step 4, and the same session audited the note. What it says of those two
steps is their author reading them again; its measurements are scripts with
saved output. Seven findings, none a wrong result at `597d28e`. All seven
were fixed or recorded on 2026-11-19 on the maintainer's word, as the note
recommends, with no change under `src/` or `R/`. Every figure is a reading
on one machine: Windows, R 4.6.2, R's LAPACK 3.12.0.

What the notes had wrong, and is fixed:

- The overrun step 3 ended is a crash. On 30 samples of which 6 are
  distinct, 11 of 18 calls end the R process on `v1.2.0` built from source,
  on `59920f4` and on `905e1ac`, `pca_imp(x, ncp = 5)` at its defaults
  among them; `eade032` and `597d28e` return on all 18. This file said
  every such run returned normally, from runs of a clone that held the full
  buffer, and concluded that no test through `pca_imp()` could see it.
  Those sentences are corrected where they stood, and the regression test
  exists
  (`tests/testthat/test-pca_imp.R`: `test_that("pca_imp returns on replicate samples with ncp above their rank", {`).
  Shown to fail: `dev/review_p2_fix_check.sh` ran that test alone against
  five installed builds, and the process ended with a segmentation fault,
  exit 139, on `v1.2.0`, `59920f4` and `905e1ac`, and passed, 14
  expectations, on `eade032` and `597d28e`. Only Windows was measured.
- A Gram that is not finite is reachable from finite input: unscaled
  values near 1e160, or near 1e200 with `scale = TRUE`, stop with "top-k
  symmetric eigensolver failed to converge" under either solver, the same
  at both ends of the phase. Step 3 and the map called it unmeasured. A
  test now pins the class and the wording
  (`tests/testthat/test-pca_imp.R`: `test_that("finite input whose squares overflow stops before the eigensolver", {`).
  Phase 4 renamed the test when the text changed.
  Run the same way it passes on the four builds of this branch and fails
  on `v1.2.0` on the wording alone, which read "`eig_sym_sel` failed to
  converge" until `c7a253b`; `dev/to-do.md` has a line for that. By reading
  alone, the stop under the exact solver is the wrapper's refusal: the
  triangle is tested before dsyevr is entered. Which status came back is
  not measured.
- The reminder for `ncp` above 32 gave "about 1e-14" under the exact
  solver. A fit of 202 iterations moved by 1.6e-13 absolute and 1.2e-11
  relative, so the line now says the figure grows with the iteration count.

What was recorded, with no change to the code:

- A Gram of order 2 to 4 under `solver = "lobpcg"` returns other values
  than the Armadillo build. Step 4's record has it now ("What changes and
  is reachable"), with a line in `dev/to-do.md`. Accepted as it is.
- `sum_two_partial()` is named in step 6's record and on Phase 4's list.
- `dev/to-do.md` has lines for the class and message of a failed triangular
  solve and for the refusal of more than `INT_MAX` missing values, which
  had none.
- The chain stays on raw buffers where "What was decided" and "Where
  Armadillo is used today" name Eigen; both say so now, and Phase 4 has the
  item for `lobpcg_detail::Block`.

The gate after the fixes, a clone of `597d28e` holding the working tree's
tests, `load_all1(timer = FALSE)`, manual tests included: the two PCA test
files, 39 tests and 352 expectations, none failed, errored or skipped.
The whole suite on such a clone, the mirai daemons loading it installed
into a private library: 22 files, 278 tests, none failed, errored or
skipped, two tests more than step 6.

Left open. Whether 1.2.x is patched for CRAN, which the measurement above
reopens and the maintainer decides. Another platform, another BLAS or
LAPACK, R 4.3.0. The three two-partial-sum loops were not folded into one
helper: Phase 4 removes them, and an edit to them now would need the
bitwise gate again.

## Phase 3: remove the dependency

- In `DESCRIPTION`, `RcppArmadillo` leaves `LinkingTo`.
- In `src/Makevars` and `src/Makevars.win`, `-DARMA_64BIT_WORD=1` leaves
  `PKG_CPPFLAGS`. The OpenMP flags are left as the maintainer has them.
- `arma_uword_bytes()` in `src/armaSVD.cpp`, its comment on the macro, and
  the test in `tests/testthat/test-mat_miss.R` that pinned it went in step 6
  of Phase 2, which took Armadillo out of that file.
- `devtools::document()` regenerates `src/RcppExports.cpp` without
  `RcppArmadillo.h`.
- `.claude/CLAUDE.md`, Layout, says `src/` is Rcpp with RcppArmadillo; it
  says what `src/` is built on instead. `git grep -i armadillo` over the
  tracked notes finds anything else.
- `devtools::check()`, once against CRAN's RcppEigen and once against the
  candidate current on that day, looked up again and not taken from this
  file. Record both versions in the commit body. Local CI on the Linux
  images is the maintainer's call.

As carried out, 2026-11-19, on `f999310`. Every build ran on scratch clones
by `bash dev/eigen_p3_gate.sh <scratch> <stage>`, which says what each
stage does; `dev/eigen_p3_check.R` runs `document()` and `check()` there.

- `DESCRIPTION`, `src/Makevars` and `src/Makevars.win` as listed above. With
  the macro gone neither Makevars assigns `PKG_CPPFLAGS`; `load_all1()`
  still appends its macros through the user Makevars, and a
  `load_all1(timer = TRUE, debug = TRUE)` build shows `-DLOC_TIMER` and
  `-DPCA_IMP_DIAGNOSTICS=1` on all 12 compile lines and no `ARMA` on any.
- `devtools::document()` took the `#include <RcppArmadillo.h>` line out of
  `src/RcppExports.cpp` and changed nothing else: not `man/`, not
  `NAMESPACE`, not `R/RcppExports.R`. `RcppEigen.h` stays included, since
  RcppEigen stays in `LinkingTo`. The build gives 511
  `-Wignored-attributes` lines, as step 6 of Phase 2 did, and no other
  compiler warning.
- Versions, looked up that day: CRAN RcppEigen 0.3.4.0.2 (Eigen 3.4); the
  candidate 0.4.9.9-5 (Eigen 5.0.1), still the newest on the r-universe.
- Capture: `f999310` and the Phase 3 tree, each on both RcppEigen versions,
  71 cases each. The Phase 3 tree against `f999310` at tolerance 0: no case
  failed and every one bitwise, on both versions. Against `76d4c60` at the
  defaults: 31 of 75 failed, the 24 `n_iter` cases of step 6 of Phase 2 and
  the 7 deliberate failures of Phase 1, every compared value bitwise.
- Suite, manual tests included, `load_all1(timer = FALSE)`, the mirai
  daemons loading the clone installed into a private library: 22 files,
  278 tests, 2844 expectations, none failed, errored, skipped or warned, on
  both versions.
- `devtools::check()` on both versions, `00install.out` confirming which
  RcppEigen include folder each used: 0 errors, 1 warning, 1 note, the same
  on both. The warning is "Undocumented code objects", 18 exports of the
  vendored sampler (`mipca()` and its helpers); the note is `rlang` in
  `Imports` and not imported from, kept on purpose because cli needs it
  (maintainer, 2026-11-19). Phase 3 changes no R code, `NAMESPACE`
  or `Imports`, so neither comes from it; not measured on `f999310`.
- `.claude/CLAUDE.md`, Layout, says what `src/` is built on now, and
  `dev/pca_imp_structure.md` and the comment in `load_all1()` no longer
  name the macro. `NEWS.md` still says the package is compiled with 64-bit
  Armadillo indices; `dev/to-do.md` has the reminder for the draft.
- This file stays until Phase 4 lands (maintainer, 2026-11-19).

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
- In `src/svd_triplet.h`, `gram_trace()` in Armadillo's order, two partial
  sums, which `dev/CPP.md` bans; kept through Phases 2 and 3 by ruling
  (step 5 of Phase 2). A plain sum in order moves results by the amounts
  step 5 measured, which put 22 of the capture's cases past its default
  tolerance, most on `criterion_final` alone.
- In `src/armaSVD.cpp`, `sum_squares_two_partial()`, the other term of the
  tail energy, in the order of Armadillo's own `dot()` loop; kept on the
  same footing (step 6 of Phase 2, which has what a plain sum moves). The
  two go together: one commit, one measurement.
- In `src/armaSVD.cpp`, `sum_two_partial()`, the sum of the row weights in
  the order of `arma::accu()`, the third loop of the same kind (step 6 of
  Phase 2: a plain sum leaves the capture bitwise and moves 4 of the 5
  diagnostics cases with unequal weights, by at most 3.9e-12). Left off
  this list until the review of Phase 2 found it missing.
- The oracle (`dev/oracle_pairs.cpp`, `dev/oracle_run.R`,
  `dev/oracle_selftest.R`, `dev/oracle_install_armadillo.R`), their `!`
  lines, and `dev/lib/`, once the last item above is gone and nothing is
  left to compare. One check of the oracle's outlives it: guard cells after
  the eigenvalue buffer of `syevr_range`, on a range that cuts a group of
  equal eigenvalues, shown to fail when the wrapper hands dsyevr the
  caller's buffer. A test through `pca_imp()` sees that overrun only where
  it does harm: `tests/testthat/test-pca_imp.R` has one, on replicate
  samples, which ends the test process on a Windows build with the defect
  ("The review of Phase 2" above). Where the cells land depends on how the
  compiler laid out the frame, so on another platform the same build may
  pass it; the guard cells do not depend on that. The
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
  error on the R side, on the main thread, never from a worker. One such
  error has a known cause to name: "top-k symmetric eigensolver failed to
  converge" is what finite input whose squares overflow gets, values near
  1e160 unscaled ("The review of Phase 2" above), where nothing failed to
  converge and the Gram was refused as not finite. The wrapper's status
  tells the two apart; the text and class are the maintainer's to choose,
  and a test pins today's.
- `lobpcg_detail::Block` in `src/lobpcg_warm.h`, a matrix type written by
  hand with eleven helpers, about 180 lines, carrying Armadillo's member
  names (`memptr()`, `n_rows`). The chain stays on raw buffers ("What was
  decided", as carried out), so the type stays; the hygiene pass decides
  its names and whether any helper is one of several of the same shape. It
  is not replaced by Eigen without a measurement that the order of every
  operation is kept.
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
- The cost of the second half of step 4 of Phase 2. Against `ff0e394`,
  `sim_pca_lobpcg` took about 1.5% longer over 11 repetitions of
  the step 4 timing (a 0.483-0.491 s, b 0.490-0.501 s); chr21 showed
  nothing measurable. Not profiled. Suspected, an assumption until
  profiled: the zero-fill of every new `Block` in `src/lobpcg_warm.h`, the
  copies of the state and the Gram's eigenvectors in and out through
  `from_arma()`/`to_arma()`, the gathered sub-blocks, and the transposes
  written out around the triangular solve. Deferred by the maintainer on
  2026-11-18 to this phase, with the hygiene pass and any `/simplify`
  calls; profile before changing anything. Step 5 removed the copies
  through `from_arma()`/`to_arma()`, and `sim_pca_lobpcg` timed 1.003
  against step 4, so they do not account for it, to the driver's
  resolution of about half a percent.
- Go through *Writing R Extensions* once more against the finished `src/`,
  not only the sections the wrapper layer needed: registration, visibility,
  interrupts, allocation, the portability list. Record each rule and where
  it is met, as the wrapper step did.
- Each finding is a proposal with a measurement, put to the maintainer
  before it is made. Nothing here licenses a hand-rolled kernel
  (`dev/CPP.md`) or a call R's headers do not declare.

### Phase 4, as carried out

Each entry is one commit, ruled by the maintainer before it was made. The
gate is `dev/eigen_p4_gate.sh`: a clone of HEAD and a clone with the change,
captured in one run on one machine, so tolerance 0 is asked where an item
is expected to be bitwise.

**Workspace queries and the shape of `qr_econ`** (ruled 2026-11-19).

- `si_la::svd_values` asks dgesdd for its workspace at every size, and is
  compiled only when `PCA_IMP_DIAGNOSTICS` is set: its one caller is in that
  build. The installed DLL imports 10 BLAS and LAPACK routines where it
  imported 11, dgesdd gone (objdump of the DLL the suite stage installed).
- `si_la::qr_econ` asks dgeqrf and dorgqr each for its own workspace, with
  `max(1, n)` as the floor of both, and keeps one form: `m >= n`, factored
  in `q`. `m < n` is `err_arg`. Its one call site stops before the call on
  a block with more columns than rows
  (`src/lobpcg_warm.h`: `internal error: a LOBPCG block has more columns than rows`).
  A square block took Armadillo's full-QR branch before and takes the one
  form now.
- `-DUSE_FC_LEN_T` in `PKG_CPPFLAGS`: not made. The manual says the macro is
  the default from R 4.3.0, which is what `DESCRIPTION` requires, so the
  flag would do nothing on a supported R.
- Measured on R 4.6.1, `0c76579` against the change: the Phase 0 capture is
  bitwise on all 71 cases under RcppEigen 0.3.4.0.2 and under 0.4.9.9.5,
  and the 49 cases of `dev/eigen_p2s6_diag.R`, the diagnostics build, are
  bitwise (worst absolute and relative difference 0 in all three). The
  suite with manual tests passes under both: 278 tests, 0 failed, 0
  skipped. Not measured: whether any case passes a square block to
  `qr_econ`; no timing, since the change adds one workspace query per QR
  and removes none of the work.
- The oracle and `dev/review_eigen_wrappers_probe.cpp` call `svd_values`
  and `qr_econ` on wide input; neither was rerun, and both now need
  `-DPCA_IMP_DIAGNOSTICS=1` to compile. No user-visible change, so no NEWS
  reminder.

**The triangular solve is plain dtrtrs** (ruled 2026-11-19).

- `si_la::solve_tri` is dtrtrs and nothing else: dtrcon, the dgelsd
  fallback, the `rcond` and `approx` outputs and the workspace floors
  computed by hand are gone. R's own `qr.coef()` solves a triangle the same
  way. A positive status makes `qr_orthonormalize_update_matvec()` return
  false, which its two callers already treat as `LOBPCG_P_COLLAPSED`, so
  that outer iteration goes to the exact solver; the R error "lobpcg_solve:
  triangular solve failed" no longer exists.
- How often the fallback fired, the measurement this item waited for: at
  `0c76579`, instrumented, no solve of the Phase 0 capture took it, and at
  `cefd88b` 4 of 19164 in `dev/eigen_p2s4_edges.R` did, all four in cases
  that set a pivot to 0 by injection. The smallest estimate
  among the others was 1.5e-14, in `big_xb_ncp40`, against the threshold of
  2.2e-16. Not measured on real data beyond the chr21 cases of the capture.
- Measured on R 4.6.1, `cefd88b` against the change: the capture is bitwise
  on 71 cases under both RcppEigen versions, and the diagnostics build on
  49. Of the 60 edge cases, 58 give bitwise the same result. The other two,
  `inj_solve_01` and `inj_solve_02`, where the solve is made to return 1,
  ended in the error and now return a result with no NA or NaN, after
  `P_COLLAPSED` at iteration 0. The three `inj_singular` cases under
  LOBPCG alone took the fallback and stalled (`STALLED` at iteration 13, 13
  and 7); now dtrtrs reports the zero and they end in `P_COLLAPSED` at
  iteration 0, 0 and 2, with the same imputed values to the bit, since both
  ways end in the exact solver. The Gram-NaN cases end in the same error as
  before. The suite with manual tests passes under both RcppEigen versions.
- The DLL imports 8 BLAS and LAPACK routines: dtrcon and dgelsd are gone.
- `dev/probe_p2s4_patch.R` instruments either form of the call; for the
  plain one it logs `rcond -1`, there being no estimate.
- NaN or Inf in the triangle or the right-hand side is no longer refused
  here; the Gram scan and the residual test of the next items are what
  stand in front of it.

**One test for NaN and Inf, two plain errors, a NaN residual fails**
(ruled 2026-11-19).

- `SVD_triplet()` tests the upper triangle of the Gram where it is formed
  (`src/gram_ops.h`: `[[nodiscard]] inline bool gram_upper_finite(`) and
  stops before either solver. No scan was added to `qr_econ` or
  `solve_tri`. The wrappers `syevr_range` and `syevd` keep theirs, so the
  exact path scans the Gram twice.
- `hybrid_topk_eig()` returns the status of the exact solver in place of a
  flag, and `SVD_triplet()` picks the text. "top-k symmetric eigensolver
  failed to converge" is gone. Finite input whose squares overflow gets
  "the cross-product of the data is not finite; values this large overflow
  when squared, rescale the input"; LAPACK's own failure gets "the
  symmetric eigensolver failed (LAPACK dsyevr, info = <n>)". Both are plain
  `Rcpp::stop()` errors with no class of their own, raised on the main
  thread: nothing in the PCA chain starts a thread (READ: RcppThread is
  included by the K-NN kernel alone). The maintainer plans to take cli out
  of the package (`dev/to-do.md`), so nothing here goes through it.
- In `lobpcg_solve()` a residual norm that is NaN or Inf ends the solve as
  `LOBPCG_BAD_SEED`, the status of the divergence guard, before the
  convergence test.
- Measured on R 4.6.1, `6b25996` against the change: the capture is bitwise
  on 71 cases under both RcppEigen versions, the diagnostics build on 49,
  and 57 of 60 edge cases. The other three are the cases that write a NaN
  into the Gram inside `lobpcg_solve()`; they raised the old text and raise
  the overflow text, by way of the wrapper's `err_nonfinite`. No
  `lobpcg_solve()` call ends differently in any of the 57 logged cases. The
  suite with manual tests passes under both versions, with the overflow
  test pinned to the new text.
- Timing, changed over base, 5 calls per run. The first two runs were
  taken on a busy machine (one run of 11.8 s against a median of 0.97) and
  settle nothing: over 5 repetitions chr21 exact 1.164, chr21 lobpcg 0.985,
  sim exact 1.066, sim lobpcg 0.945; over 11, chr21 exact 1.016 and sim
  exact 1.203. A third, over 11 repetitions: chr21 exact 1.013
  (a 0.329-0.445 s, b 0.326-0.368 s), sim exact 0.998 (a 0.902-1.711 s,
  b 0.900-2.352 s). The scan is not visible at the driver's resolution.
- The 60 edge cases exercised neither the LAPACK-failure text, which needs
  dsyevr to return `info > 0`, nor the residual test: the Gram-NaN
  injections end at the first Rayleigh-Ritz, where `si_la::syevd` refuses
  the small Gram, and never reach a residual. Two injections were added
  for them (maintainer, 2026-11-19), with 8 cases, 68 in all:
  `eiginfo` makes the exact solver report `info = 1`, and `resnan` sets the
  first column of `A * X` to NaN before the residuals of one iteration.
  Measured, `6b25996` against `f7afaac`, each injection logged as hit in
  all its 4 cases: the 4 `eiginfo` cases, under either solver, raised the
  old text and raise "the symmetric eigensolver failed (LAPACK dsyevr,
  info = 1)". The 4 `resnan` cases end the solve as `BAD_SEED` at the
  injected iteration, where the build before ended it as `RR_FAILED` at
  the same iteration: there the NaN went one step further, into the small
  Gram, which `si_la::syevd` refused. Both send the outer iteration to the
  exact solver and the results are bitwise the same. So the residual test
  is shown to fire, and this injection does not show a case it alone
  catches; the one known, a NaN that passes a dsyevd without the wrapper's
  scan (Phase 2, step 4), cannot arise while the wrapper scans. Neither
  check has a test in `tests/`: both need an instrumented build.

**`matmul` removed; a block of one column goes to dgemv** (item 9;
maintainer, 2026-11-19). `si_la::matmul`, which kept the order of
Armadillo's tests (a result of one row first, then a right operand of one
column, then dgemm), is gone. The LOBPCG products go through one function,
`product_into()` in `src/lobpcg_warm.h`, with one rule of its own: a right
operand of one column is a vector and the product is dgemv on the left
operand; any other is dgemm. `si_la::gemv` stays for it. The product of the
diagnostics build in `src/armaSVD.cpp` is `si_la::gemm`.

- Two forms were built and measured against `6927d58`, each in one run of
  `dev/eigen_p4_gate.sh` (R 4.6.1; the figures are the same under both
  RcppEigen versions). A: dgemm for every product, `si_la::gemv` removed
  too. B: the rule above. Neither is bitwise:

  | | A, dgemm always | B, one column to dgemv |
  | --- | --- | --- |
  | capture, 71 cases, not bitwise | 6 | 6 |
  | worst on the imputed values | 1.19e-10 abs, 5.3e-9 rel | 9.25e-11 abs, 3.3e-9 rel |
  | worst on `criterion_final` | 2.65e-7 rel | 6.94e-8 rel |
  | diagnostics, 49 cases, not bitwise (past the defaults) | 21 (16) | 21 (16) |
  | edges, 68 cases, not bitwise | 52 | 51 |
  | worst on the imputed values | 8.16e-8 abs, 2.4e-6 rel | 6.57e-8 abs, 1.5e-6 rel |
  | edges where `n_exact`, `n_lobpcg_ok`, `n_lobpcg_bad` differ | 10 | 10 |
  | suite with manual tests, both RcppEigen | 278 pass | 278 pass |

  The six capture cases are the same in both, all under LOBPCG:
  `pca_chr21_lobpcg`, `pca_x1_em_lobpcg`, `pca_x1_lobpcg` (the worst),
  `pca_xb_auto`, `pca_xb_lobpcg`, `pca_xt2_lobpcg`; the relative figures on
  the imputed values are over the cells that hold a nonzero value. No case
  under the exact solver moves. In the diagnostics compare the largest
  figures (A 2.05e-2 abs, B 1.54e-4 abs, both near 1 relative) are in one
  per-iteration vector of the kernel's log, not identified further; the
  imputed values there move as in the capture.
- Why B moves at all, when it calls dgemv as Armadillo did for one column:
  the capture's fits ask for more than one component, so the products that
  differ are those B no longer sends where Armadillo's order did, a result
  of one row and several columns (dgemv on the right operand then, dgemm
  now) and a 1 x 1 result (the other operand as the vector). Which sites
  form them was not traced.
- Over the 1316 calls of `lobpcg_solve()` in the edge cases, base, A, B:
  converged 1178, 1174, 1174; stalled 38, 45, 42; at `maxiter` 34, 32, 34;
  `BAD_SEED` 33, 32, 33. The other statuses are the same.
- Timing, with one component, where every product is matrix times vector
  (`pca_lobpcg1`, added to `dev/eigen_time_worker.R`), A as a and B as b,
  11 repetitions: chr21 1.011 (a 0.150-0.215 s, b 0.150-0.231 s), sim 0.986
  (a 0.418-0.606 s, b 0.417-0.737 s). A first run of 5 on a busier machine
  gave 0.914 and 0.706 with a range of 0.504-2.965 s on a, and is not
  relied on. The gate's own four workloads, base against B, 5
  repetitions: 1.014, 0.909, 0.999, 1.000. So B is not measurably faster
  than A here; it was chosen as the convention, dgemv for matrix times
  vector, not for speed or for bits.
- The move is of the size item 10 makes and is absorbed by the reference
  capture that follows it (item 10b). `dev/to-do.md` has the reminder line.
- The oracle no longer compiles: `dev/oracle_pairs.cpp` calls
  `si_la::matmul`. Item 14 removes it.
- Identified during item 10: the per-iteration vector with the largest
  figures in the diagnostics compare is `eig_lobpcg_max_rel_res`, the
  residual LOBPCG stopped at, a number near 1e-9 that is noise below the
  tolerance it passed.

**The trace and two sums as plain sums** (item 10; maintainer, 2026-11-19).
The three reductions kept in Armadillo's order, two partial sums over the
even and the odd positions, are plain sums in order: `gram_trace()` in
`src/svd_triplet.h`, and in `src/armaSVD.cpp` the sum of the squared
leading singular values and the sum of the row weights, both now loops at
their one call site (`sum_squares_two_partial()` and `sum_two_partial()`
are gone). `dev/CPP.md` bans the form they had.

- Measured against `dcc1364`, the item 9 build, in one run of
  `dev/eigen_p4_gate.sh` (R 4.6.1; the same figures under both RcppEigen
  versions). Relative figures on the imputed values are over the cells
  that hold a nonzero value.

  | | not bitwise | worst on the imputed values | worst on `criterion_final` |
  | --- | --- | --- | --- |
  | capture, 71 cases | 27 (21 past the defaults) | 6.24e-10 abs, 3.4e-8 rel | 3.41e-6 rel |
  | diagnostics, 49 cases | 49 (19 past the defaults) | 2.48e-8 abs, 1.1e-7 rel | |
  | edges, 68 cases | 61 (60 past the defaults) | 7.56e-8 abs, 1.2e-6 rel | see below |

  Suite with manual tests: 278 pass under both RcppEigen versions.
- In the capture the worst case on both figures is `pca_x1_lobpcg`. The 21
  past the defaults are all on `criterion_final`, 3 of them on the imputed
  values too. Under the exact solver the imputed values move by at most
  6.0e-15 absolute and `criterion_final` by at most 1.4e-9 relative; the
  larger figures are all under LOBPCG, where the changed tail energy
  reaches the iteration through the next warm start.
- Edges: `n_exact`, `n_lobpcg_ok` or `n_lobpcg_bad` differ in 10 cases;
  `n_iter` and `converged` in none. Over the calls of `lobpcg_solve()`,
  base and item 10: 1316 and 1318 calls; converged 1174, 1172; `BAD_SEED`
  33, 37; the other endings the same. The worst case on the imputed
  values is `big_xb_ncp25`.
- `criterion_final` on the rank-deficient edge cases (duplicated columns,
  `ncp` at or past the rank) is not a number to read: there the tail
  energy is a difference of two equal sums, and the criterion moved from
  NaN to 0.085 (`dup_xd_ncp07`), from 1.7e-15 to 0.5 (`dup_xd_ncp07_em`),
  from 0.013 to 0.496 (`dup_xd_ncp12`) and from NaN to Inf
  (`dup_xdc_ncp08`), while the imputed values of those cases agree to
  5.5e-13 absolute and `converged` is `TRUE` on both sides. It was NaN
  before this item too, so the defect is older than the change; it is
  reported to the maintainer, not fixed here. The maintainer has since
  ruled on it: "`criterion_final` as `NA`" below.
- Timing, 5 repetitions, base as a: 1.298, 0.792, 1.209, 1.060 with ranges
  as wide as 1.1-7.1 s on one side, a busy machine. The change is three
  scalar loops outside every BLAS call, so no effect is expected and this
  run cannot show one either way; it was not repeated.
- `dev/to-do.md` has the reminder line. The reference capture that follows
  (item 10b) absorbs this move and item 9's.

**A new reference capture, `aa1ceee`** (item 10b; 2026-11-19). Taken in a
scratch clone of `aa1ceee`, the item 10 commit, and copied into `dev/ref/`
on this machine with a second capture, `aa1ceee-b`, and one against the
candidate RcppEigen, `aa1ceee-eigen5`. Like every capture it is ignored
and regenerable: clone `aa1ceee` and rerun `dev/ref_capture.R`. From here
an item of Phase 4 that should not move results is compared against it;
`dev/eigen_p4_gate.sh` compares a clone of HEAD with the working tree, so
it needs no label.

- 71 cases. Noise floor: `aa1ceee` against `aa1ceee-b` at tolerance 0, all
  bitwise identical. Against the candidate at tolerance 0: all bitwise
  identical. `dev/ref_compare_selftest.R aa1ceee` caught the 13 planted
  changes at tolerance 0 and the 12 expected at the defaults, and reported
  the 0 turned into -0 as passing but not bitwise.
- The whole move since `76d4c60`, the capture of the Armadillo build (75
  cases): worst on the imputed values 5.97e-10 absolute and 3.1e-8
  relative, worst on `criterion_final` 3.48e-6 relative, all three in
  `pca_x1_lobpcg`; 27 cases not bitwise. The other differences are ones
  ruled on before: `n_iter` an integer where it was a double (24 cases),
  and 4 cases of functions since removed (`col_vars_x2_cores2`,
  `int_mean_imp_col`, `int_mean_imp_col_cores2`,
  `int_mean_imp_col_inf_error`).

**`criterion_final` as `NA`** (maintainer, 2026-11-19: "option 1", then
"Code option 1"). Committed at the maintainer's request with the test as
it now stands NOT RUN.

- The ruling: a fit that stopped because its objective fell below
  `threshold`, and not because its criterion did, reports `criterion_final`
  as `NA`. The fit itself is unchanged.
- What the objective is in the six rank-deficient edge cases (diagnostics
  build of `894826d`, `threshold` 1e-6, 15 iterations each). All six
  stopped on the objective alone:

  | case | last two objectives | criterion from them |
  |---|---|---|
  | `dup_xd_ncp05` | 0, 0 | NaN |
  | `dup_xd_ncp07` | 1.60e-13, 1.46e-13 | 0.085 |
  | `dup_xd_ncp12` | 2.91e-14, 1.47e-14 | 0.496 |
  | `dup_xd_ncp07_em` | 2.84e-14, 1.42e-14 | 0.5 |
  | `dup_xdc_ncp05` | 1.00e-12, 1.29e-13 | 0.871 |
  | `dup_xdc_ncp08` | 0, 1.78e-15 | Inf |

  In `dup_xdc_ncp05` the objective falls steadily from 1.9e-2, so its
  ratio is a real one; in the others it is rounding.
- The kernel edit is one assignment after the loop in `src/armaSVD.cpp`.
  Gate against `894826d` (scratch clone, R 4.6.1, both RcppEigen versions):
  capture 71 of 71 bitwise, diagnostics 49 of 49 bitwise, edges 68 of 68
  bitwise on everything but `criterion_final`, which is `NA` in the six
  cases of the table and unchanged in the other 62.
- The suite: 278 of 279 passed under both versions. The one failure was
  the new test, which expected `NA` at `ncp = 5` on the replicate-sample
  input under the default settings and got 1.1e-16. There both arms of the
  stop rule fired: rounding put the ratio below `threshold`. Measured on
  that input over `ncp` 4, 5, 6, 7, 12, the three solvers and both
  methods: `NA` everywhere from `ncp = 5` up, except 1.1e-16 at `ncp = 5`
  with `method = "regularized"`; at `ncp = 4`, below the rank, a criterion
  as before.
- So the rule as ruled leaves one gap: a criterion that is rounding is
  still reported when it lands below `threshold`. `NA` whenever the
  objective is below `threshold` would close it; that is a different rule
  and is put to the maintainer, not coded.
- The test was rewritten to what holds on every build: a converged fit
  reports `NA` or a criterion below `threshold`, never NaN or a larger
  number. That version was not run before the commit; the suite run was
  stopped at the maintainer's request.
- Run since, on R 4.6.2, in a clone of `e86d532` through the `suite` stage
  of `dev/eigen_p4_gate.sh`: 279 tests, 0 failed, 0 skipped, 2856
  expectations passed, under RcppEigen 0.3.4.0.2 and under 0.4.9.9.5.
- Shown to fail: the test file of `e86d532` on a clone of `894826d`, the
  build before the rule (`dev/criterion_na_check.R`), fails at all four
  `ncp`, with `criterion_final` 1, NaN, Inf and 0.934 at `ncp` 5, 6, 7 and
  12; the other 8 expectations pass. `dev/to-do.md` has the reminder line.

**Rulings on the last items** (maintainer, 2026-11-19: "agree all" to the
author's recommendations). Recorded here; the first three are carried out
below, "`Block` off Armadillo's names", "The wrapper layer's own test" and
"The oracle removed".

- `lobpcg_detail::Block` is renamed off Armadillo's member names. The
  updates that alias their own operand stay as they are until a profile
  shows they cost something.
- The test of the wrapper layer's own, guard cells after the eigenvalue
  buffer of `syevr_range` included, is a script in `dev/` that compiles
  through `Rcpp::sourceCpp()` and is run by hand. No test-only export is
  added to the package.
- The oracle, its four files, their `!` lines and `dev/lib/arma` are
  removed; it has not compiled since `matmul` went. "What the wrappers
  reproduce" in `dev/pca_imp_structure.md` is rewritten to say what the
  wrappers do.
- `syevr_range` goes on asking dsyevr for its workspace on every call.
- The pass through *Writing R Extensions* against the finished `src/` is
  its own item. `src/loc_timer.h` keeps its macro names.
- Not ruled at first: whether the OpenMP flag stays in `src/Makevars`. The
  author gave no recommendation. Ruled since (maintainer, 2026-11-19,
  confirming the ruling of 2026-11-14 recorded upstream): both Makevars drop
  it. Carried out below, "The OpenMP flag dropped".

**`Block` off Armadillo's names** (the first ruling above). In
`src/lobpcg_warm.h` the members of `lobpcg_detail::Block` are `nrow`,
`ncol`, `cells`, `data()`, `empty()` and `clear()`, where they were
Armadillo's `n_rows`, `n_cols`, `mem`, `memptr()`, `is_empty()` and
`reset()`; `src/hybrid_topk_eig.h` follows. A rename only: no line changes
what it does, and no helper was merged or removed.

- `dev/probe_p2s4_patch.R` writes its anchors and the lines it inserts in
  the names the file holds (`nm()`), so one script still patches a clone
  from before the rename and one after.
- Measured on the second machine (R 4.6.2, R's BLAS and LAPACK replaced by
  OpenBLAS 0.3.34), `e86d532` against the change, in one run of
  `dev/eigen_p4_gate.sh`: the capture bitwise on 71 of 71 cases under both
  RcppEigen versions, the diagnostics build on 49 of 49, the edges on 68 of
  68, at tolerance 0. The probe made the same 17 insertions in both edge
  clones, and their event logs are byte for byte the same, 83276 lines. The
  suite with manual tests: 279 tests, 0 failed, 0 skipped, under both
  versions. No timing: the rename adds no work.
- The earlier entries of Phase 4 were measured on the first machine, R
  4.6.1. On the second, HEAD's edges end 1187 of 1316 `lobpcg_solve()`
  calls as `CONVERGED` and 26 as `BAD_SEED`, against 1172 of 1318 and 37 in
  item 10's record: the last digits of another BLAS. A gate compares two
  builds on one machine, so it is not affected; a stored capture compared
  across the two machines would be. The captures of those entries, `aa1ceee`
  among them, are in `dev/ref/` on the first machine only.

**The wrapper layer's own test** (the second ruling above).
`dev/si_lapack_test.R` compiles `dev/si_lapack_test.cpp` beside
`src/si_lapack.h` and `src/si_lapack.cpp` through `Rcpp::sourceCpp()` and
checks every function of the header against what its comment says it does.
Run by hand from the package root, like the oracle; nothing under `R/`,
`src/` or `tests/` changes, and nothing is exported for it.

- The references are R's own: products from R's own loops
  (`options(matprod = "internal")`, which calls no BLAS), `eigen()`,
  `svd()`, `qr()`, `backsolve()` and `forwardsolve()`. Armadillo is not
  loaded. Eigenvectors are checked by their residual and orthonormality,
  since a range that cuts a group of equal eigenvalues has no unique
  answer, and `r` of `qr_econ` up to the sign of each row.
- Checked with `identical()`: 16 guard cells after every buffer a wrapper
  writes, and n more after the eigenvalues and eigenvectors of
  `syevr_range`; the rows past the last of an output with a larger leading
  dimension; the triangle `syrk` does not name; and every output of a call
  that is refused or writes nothing. A size past `int` is asked for with a
  buffer of one cell, which the wrapper must refuse before it reads one.
- Both units are compiled with `-DPCA_IMP_DIAGNOSTICS=1`, so `svd_values`
  is tested too, and linked with `$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)` as
  the package is. Without RcppArmadillo's plugin nothing else adds them: the
  first build failed to link.
- A probe in the same file calls dsyevr as `syevr_range` does, with W given
  n cells, and counts what it wrote past `iu - il + 1`, so the guard cells
  after `w` are shown to have something to catch. The run fails if no case
  makes it write any.
- Measured on the second machine (R 4.6.2, R's BLAS and LAPACK replaced by
  OpenBLAS 0.3.34), on the source at `2bc9d59`: 238 cases, all passing
  (`gemm` 40, `gemv` 20, `syrk` 40, `syevr_range` 61, `syevd` 29,
  `qr_econ` 11, `solve_tri` 23, `svd_values` 14). Worst absolute
  difference over its scale 5.8e-15 (`syevr_range`, order 300), worst
  relative 6.7e-13 over the cells above 1e-6 of the scale (`gemm`, 40 x
  30), worst residual or orthonormality 8.3e-15. The bounds are 1e-12,
  1e-7 and 1e-12: the first and last about a hundredfold above the
  reading, the relative one set by its floor, where one unit in the last
  place of the scale is 2.2e-10. `svd_values` is bitwise on all 9 values,
  R's `svd()` calling dgesdd too. dsyevr wrote past `iu - il + 1` in 10 of
  the 12 ranges that cut a group, 39 cells in all.
- Shown to fail: `dev/si_lapack_selftest.R` breaks one wrapper at a time
  in a copy, 24 ways, and the test named the routine aimed at each time,
  and no other; it passes on the unbroken source. The ruling's case:
  `syevr_range` handing dsyevr the caller's `w` and copying nothing, as the
  chain did before step 3 of Phase 2, fails on exactly those 10 ranges and
  on the guard cells after `w` alone, the first `iu - il + 1` cells being
  right. The others: the other triangle in `syevr_range`, `syevd`, `syrk`
  and `solve_tri`; the NaN scan dropped, or made to read the other
  triangle; one eigenvalue too few copied; `il > iu` let through; alpha
  where beta belongs in `gemm` and `gemv`; C or y left alone where an
  empty product must scale it; `syrk` scaling whole columns when `k` is 0;
  `ldc`, or `lda` of a transposed operand, unchecked; `m < n` let through
  in `qr_econ`; ones below the diagonal of `r`; the transpose solved; a
  zero pivot reported as success; the sign of a singular value; the wrong
  refusal code; `'C'` taken as a transpose. Where an argument check was
  dropped, OpenBLAS returned at once without a word, so dgemm left a
  status of 0 and LAPACK its negative info; under R's own BLAS and LAPACK
  xerbla raises an R error, which the test records as a failing case (not
  run).
- Not covered. `err_count` and `err_alloc`: no input reaches them. A size
  check let through is not among the mutations: the size cases hand the
  wrapper one cell, so a broken copy crashes the child instead of failing.
  A workspace that ignores LAPACK's answer, or a nonzero abstol, leaves the
  eigenpairs right to the tolerance; the oracle caught both by comparing
  bits with the chain's own dsyevr call, which no longer exists apart from
  the wrapper, so nothing here does.
- Not run on the first machine or on R's own BLAS, and under one
  RcppEigen only: nothing in either file includes Eigen. The oracle stays
  until the next ruling removes it.

**The oracle removed** (the third ruling above). `dev/oracle_pairs.cpp`,
`dev/oracle_run.R`, `dev/oracle_selftest.R` and
`dev/oracle_install_armadillo.R` are deleted with their `!` lines in
`dev/.gitignore`, and `dev/lib/arma`, the private RcppArmadillo library
they ran against, is deleted on the second machine. It is ignored, so git
holds none of it, and a copy on the first machine is not touched by this.
Nothing was built or run: the oracle has not compiled since `97525a2`, and
the one line under `src/` that changes is a comment.

- The four files are the same from `59920f4` to `141a9d0`, the last commit
  that has them: `git show 141a9d0:dev/oracle_pairs.cpp`,
  `git show 141a9d0:dev/oracle_run.R`, `git show 141a9d0:dev/oracle_selftest.R`
  and `git show 141a9d0:dev/oracle_install_armadillo.R` print them. The
  version `dev/review_eigen_wrappers.md` reviewed is `79cce0f`'s, and that
  note now says so, which also takes its two stale anchors into the oracle
  off the checker's list.
- "What the wrappers reproduce" in `dev/pca_imp_structure.md` is now "What
  the wrappers do": what each wrapper refuses before its routine, where it
  calls none, its workspace, and the reasons for the eigenvalue buffer of
  `syevr_range`, the tests for NaN and Inf and the workspace lengths, with
  the measurements behind them. What Armadillo did, as read in
  `armadillo_bits`, is gone from the map, and so are "Products of two
  `arma::mat`" and the Armadillo side of "Calls no wrapper stands in for";
  `git show 141a9d0:dev/pca_imp_structure.md` has the section as it was.
  Layer 6 says how the LOBPCG products are formed without pointing at
  Armadillo's rules, and `src/si_lapack.h` points at the new section and at
  `dev/si_lapack_test.R` in place of the oracle.
- The scripts that ran the oracle stay, as records, and each says in its
  header where the oracle's files are now: `dev/probe_p2s4_syevd_lwork.R`,
  `dev/review_eigen_wrappers_mutants.R`, `dev/review_eigen_wrappers_probe.R`,
  `dev/review_eigen_wrappers_probe.cpp` and `dev/review_eigen_wrappers_run.R`.
  The earlier entries of this note name the oracle as it was and are not
  changed.
- What the oracle checked and nothing checks now, from "The wrapper
  layer's own test": a workspace that ignores LAPACK's answer, and a
  nonzero `abstol`.

**The OpenMP flag dropped** (the last ruling above). `src/Makevars` and
`src/Makevars.win` no longer pass `$(SHLIB_OPENMP_CXXFLAGS)`, in
`PKG_CXXFLAGS` or in `PKG_LIBS`. The two entries in `dev/to-do.md` that
covered it are deleted, and `dev/CPP.md` and `dev/knn_imp_structure.md` no
longer say the flags are there.

- What the flag still reached, read with Rtools' `nm` on the gate's two
  builds (second machine, R 4.6.2): before the change, `mipca.o` alone
  calls the OpenMP runtime (`GOMP_parallel`, `omp_get_max_threads`,
  `omp_get_num_threads`, `omp_get_thread_num`), from Eigen's threaded
  matrix product, and the DLL holds 122 `GOMP_` and `omp_` symbols. After
  it, no object and no DLL symbol names either. With Armadillo gone,
  slideimp's own code reached none of it.
- `dev/eigen_p4_gate.sh`, stages capture and suite: at tolerance 0, 71
  cases, 0 failed and 0 not bitwise, under CRAN's RcppEigen and the
  candidate. The suite with manual tests, under each: 279 tests, 0 failed,
  0 skipped, 2856 expectations.
- The sampler's draws are not in the capture, so whether Eigen's product
  on one thread changes their bits is not measured here. The suite's
  sampler tests pass.

## Not settled yet

- Whether the noise floor measured on one machine holds on the BLAS of the
  rhub Linux images. Until a capture is compared there, the default
  tolerances are a reading from one BLAS.
- Whether the R callers are tidied once the return types are free to change.
- What the timings say. They are recorded here when each phase has been
  run, as readings on one machine, not as numbers to pin. Phase 0's reading
  is the driver's noise band, under "Phase 0".