# Removing Armadillo from src/

The plan for taking RcppArmadillo out of the package's own C++. Work happens
on the `eigen` branch, which starts from `e5b3cf7` on `dev` (the sampler
port). Delete this file when the branch merges; git history keeps it.

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

## What was decided

- Code that uses Armadillo only to hold numbers - containers, indexing,
  element loops - becomes plain Rcpp: raw pointers into R memory, `R_xlen_t`
  or `std::size_t` indices, the same loops.
- Every BLAS and LAPACK call goes through R's own headers, `R_ext/BLAS.h` and
  `R_ext/Lapack.h`, with `F77_CALL` and `FCONE`, the way `src/mipca.cpp`
  already does. That covers the calls the PCA chain makes today through
  Armadillo's wrappers (`dgemm_`, `dsyrk_`, `dsyevr_`, `dsyevd_`) and the
  LAPACK routines behind Armadillo's own algebra (`qr_econ`, `solve`, `inv`,
  `svd`).
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
  the same dims and the same attributes as on `e5b3cf7`; the harness
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

## Which model does which step

Opus does most of the plan. Fable takes the steps where the work is reading
exactly what an Armadillo or LAPACK call does, or chasing a numerical
difference: a mistake there compiles, can pass a loose tolerance, and shows
up later as a failure on a platform that was not tested. When a step is in
doubt, it goes to Fable.

| Step | Model | Why |
|---|---|---|
| Phase 0: capture and compare scripts, noise floor, timing driver | Opus | Scripting against rules this file already states |
| Phase 1, steps 1 and 2 | Opus | Storage only; the compare script settles bitwise-ness and the RNG question |
| Phase 1, step 3 (the K-NN chain) | Opus | Storage only, but threaded: check the worker rule and 64-bit indexing above |
| Phase 2, steps 1 to 3 | Opus | One BLAS or LAPACK call per wrapper, arguments in plain sight |
| Phase 2, step 4 (`src/lobpcg_warm.h`) | Fable | What `qr_econ`, the triangular `solve` and `norm` call in armadillo_bits and how each fails; LOBPCG iterates, so a last-bit change can grow |
| Phase 2, step 5 | Opus | Storage and `trace` |
| Phase 2, step 6 (`src/armaSVD.cpp`) | Fable | Which LAPACK routine `svd` calls and its fallback, and the switch of every declaration to `R_ext/BLAS.h` and `R_ext/Lapack.h` with `FCONE`, which is where other platforms warn |
| Any step whose results move beyond the noise floor | Fable | Finding which call moved them before going on |
| Phase 3 | Opus | Mechanical |

After Phase 2, before Phase 3, each model reviews the other's Phase 2 work
cold: Fable reads the diffs of steps 1 to 3 and 5, Opus reads steps 4 and 6.
The reviewer gets the diff and the READ notes added to
`dev/pca_imp_structure.md`, not the author's account of why it is right.

## Where Armadillo is used today

At `e5b3cf7`. "Storage" means containers, indexing and element loops only.

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

Before any source edit, while the branch is still `e5b3cf7`, two scripts in
dev/, tracked (each gets its `!` line in `dev/.gitignore` in the commit that
adds it):

- A capture script. On one build it runs every export whose signature uses
  Armadillo, and the R entry points that reach them (`knn_imp()`,
  `pca_imp()` with each solver, `col_vars()`, `mean_imp_col()`,
  `mat_miss()`, `sim_mat()`, the `tune_imp()` NA placement), on
  `sim_mat()` output and on `dev/GSE286313.EPICv1_chr21.qs2`, under fixed
  seeds. It saves every result, with its type, dims and attributes, to an
  ignored folder under dev/. That folder is regenerable data: rerun the
  script at `e5b3cf7`.
- A compare script. It reads a capture and checks a build against it. For a
  double result it reports `max(abs(a - b))` and `max(abs(a - b) / abs(b))`
  (over the cells where `b != 0`, saying how many), asserts the worst case
  against a tolerance, and reports whether the two are bitwise identical
  without asserting it. Integer results, index maps and the NA placement use
  `identical()`. Type, dims and attributes must match exactly.

The noise floor comes first: capture the reference build twice, in two R
processes, and compare them. A threaded BLAS need not be deterministic, and
no tolerance is chosen tighter than what the unchanged build reproduces
against itself.

Timing, in the same phase. Install the reference build into its own library
in a scratch folder, and each converted build into another; a driver
alternates Rscript processes between the two libraries and reports medians
over interleaved repetitions of `pca_imp()` (both solvers) and `knn_imp()`
on the chr21 data and on a larger `sim_mat()`. Both builds are the same
flags; `load_all1()` builds are not mixed with installed ones. Removing
Armadillo also removes the OpenMP Armadillo turns on by itself in an OpenMP
build, so this is the measurement of what that costs or saves.

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

Leaves first, in this order, each step compiling and passing on its own:

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
in `dev/pca_imp_structure.md` as the step lands.

Declarations switch last. `src/eig_sym_sel.h` explains why it declares
`dsyevr_` itself: `R_ext/Lapack.h` redeclares symbols Armadillo has already
declared, with different types, and warns in every translation unit that
includes both. So steps 1 to 5 replace the containers while the BLAS and
LAPACK calls keep Armadillo's declarations, and the commit that removes the
last Armadillo type from the translation unit (step 6) moves every call to
`R_ext/BLAS.h` and `R_ext/Lapack.h` and drops `RcppArmadillo.h` from it.

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
- `devtools::check()`. Local CI on the Linux images is the maintainer's
  call.

## Not settled yet

- The tolerances, until Phase 0 has measured the noise floor.
- Whether the R callers are tidied once the return types are free to change.
- What the timings say. They are recorded here when Phase 0 and each phase
  have been run, as readings on one machine, not as numbers to pin.
