# Review: the wrapper layer and the oracle

Checkpoint "The wrapper layer and the oracle, before any call site uses
them" of `dev/eigen_conversion.md`. Range under review: `7d0c176..79cce0f`
on `si-mi` (one commit, `79cce0f`), working tree clean at the start.
`/code-review` at high effort, one session, no verify pass: the session
that wrote the wrappers audits this note.

Nothing under `src/` and no oracle file was edited, nothing was committed,
and nothing was built in the working tree: every build was staged in a
scratch folder outside the repository. R 4.6.1 ucrt, GCC 14.2.0, R's LAPACK
3.12.0 (`La_version()`), `R_LIBS_USER` set to this machine's shared library,
Rcpp 1.1.2, RcppArmadillo 15.6.2.1 in `dev/lib/arma`. Run 2026-11-18. Every
figure below is a reading on this one machine with this one BLAS and LAPACK.

The oracle and `dev/lib/arma` were removed after `141a9d0` (Phase 4 of
`dev/eigen_conversion.md`). The version reviewed here is `79cce0f`'s, which
`git show 79cce0f:dev/oracle_pairs.cpp`, `git show 79cce0f:dev/oracle_run.R`,
`git show 79cce0f:dev/oracle_selftest.R` and
`git show 79cce0f:dev/oracle_install_armadillo.R` print; the anchors below
into those files are to that version.

Probe scripts, all in `dev/` (ignored, like the earlier review notes), each
run from the package root with `R_LIBS_USER` set as above:

- `dev/review_eigen_wrappers_run.R <lib> <out dir> check oracle selftest` -
  runs `dev/oracle_install_armadillo.R --check`, `dev/oracle_run.R --verbose`
  and `dev/oracle_selftest.R` unchanged, each in a child process, and keeps
  the output.
- `dev/review_eigen_wrappers_mutants.R <lib> <out dir>` - sixteen wrong
  wrappers the self-test does not try, by the self-test's own mechanism
  (a copy of `src/si_lapack.cpp` changed by textual replacement, the oracle
  run on the copy).
- `dev/review_eigen_wrappers_probe.R <lib> <out dir>` with
  `dev/review_eigen_wrappers_probe.cpp` - direct calls of the wrappers the
  oracle does not make, each group in a child process with a time limit.
- `dev/review_eigen_wrappers_procs.R` - lists R processes, to confirm the
  cases that never return left nothing running. They did not.

## Summary

| # | Where | Finding | Severity |
|---|---|---|---|
| 1 | `src/si_lapack.h`, `syevr_range` | `w` is documented as `iu - il + 1` cells; dsyevr writes up to `n`. Measured: 17 doubles past the end on a 30 x 30 Gram of rank 5 asked for 8 pairs | high |
| 2 | `src/si_lapack.cpp`, `syevr_range` | With NaN in the matrix and every eigenpair asked for, the call does not return (60 s limit hit in 4 of 4 cases, n = 6 and n = 60) | high |
| 3 | `src/si_lapack.cpp`, `syevd`, `syevr_range` | NaN or Inf in the matrix comes back as status 0 with NaN results: `syevd` in 12 of 20 cases, `syevr_range` in 1 of 16 that returned | medium |
| 4 | `dev/oracle_run.R`, class "same" | Bitwise identity is counted and never asserted, so a wrapper that calls another routine passes: `matmul` sent wholly to dgemm passes with 25 of 83 cases no longer bitwise | medium |
| 5 | `dev/oracle_pairs.cpp`, `pair_syevr` | `syevr_range` is checked only against another algorithm at 1e-11 and 1e-7, so its options and workspace are not pinned: `abstol = 1e-12` passes | medium |
| 6 | `dev/oracle_run.R`, inputs | No input lies near a rule the wrappers copy from Armadillo, so moving the rule passes: the rcond threshold by 1000 either way, the `n * eps` cut, the NaN test, the 1024-cell rule, the order of the one-row and one-column tests | medium |
| 7 | `dev/oracle_pairs.cpp`, `wrapper_refusals`, and the eigen pairs | Failure paths with no case: `err_count`, `err_alloc`, the leading-dimension checks of a transposed operand, which code a failure carries, non-finite input to the two eigen wrappers, a range that cuts equal eigenvalues | medium |
| 8 | `dev/pca_imp_structure.md`, "What the wrappers reproduce" | `auxlib::eig_sym_dc()` is left NOT READ; it refuses NaN and Inf in the upper triangle, which is why the oracle's syevr pair cannot agree on failure for such input | low |

14 of the 16 extra wrong wrappers pass the oracle (findings 4 to 7). The
answer to "can the oracle pass with a wrong wrapper the self-test does not
try" is yes, measured.

What holds up, measured or read:

- `dev/oracle_install_armadillo.R dev/lib/arma --check`: library 15.6.2.1,
  CRAN source 15.6.2.1, Windows binary 15.6.2.1, "up to date, nothing
  installed".
- `dev/oracle_run.R dev/lib/arma --verbose`: exit 0, "ORACLE PASSED", 627
  cases, 16 s. Bitwise identical: gemm 148 of 148, gemm_block 6 of 6, gemv
  46 of 46, syrk 116 of 116, syrk_block 3 of 3, matmul 83 of 83, syevd 46
  of 46, qr_econ 17 of 17, solve_tri 32 of 32, solve_tri_chain 4 of 4,
  svd_values 14 of 14; matmul_tiny 18 of 42 (worst abs 8.88e-16, rel
  3.84e-15); syevr_range 10 of 34 (worst abs 1.71e-13, rel 7.57e-10). These
  are the figures the plan records.
- `dev/oracle_selftest.R dev/lib/arma`: 12 of 12 mutations caught, each
  naming the routine aimed at, the unbroken source passes, exit 0,
  "SELFTEST PASSED", 307 s.
- Each wrapper against `armadillo_bits` 15.6.2, read before the READ notes
  were opened: routine, options, workspace and failure agree with
  `src/si_lapack.cpp` and with "What the wrappers reproduce", except
  findings 1 and 8. Detail under "The independent read".
- No argument reaches xerbla. Read: every option character, dimension pair,
  leading dimension and workspace minimum the eleven routines test is
  tested first. Run: 2950 calls over every option and every dimension from
  0 to 4, leading dimensions passed as the row count (0 included): all
  returned 0, no R error, no line from the BLAS. Eight calls the oracle
  does not try return what they should: the four leading-dimension checks
  of a transposed operand in `gemm`, `gemv` and `syrk`, an empty inner
  dimension with a leading dimension of 0, `err_alloc` from `solve_tri`,
  and the boundaries n = 32767 for `syevd` and n = 82595525 for
  `syevr_range`.
- No int narrowing: every size goes through `to_int`, every workspace
  length through `work_len` or `fits_int` in double, and every buffer
  length (`2 * count`, `3 * n`, `8 * p`, `n * n`, `n * nrhs`, `m * n`) is
  computed in `std::size_t`. `syevd` refuses n >= 32767, where
  `1 + 6 n + 2 n^2` leaves int; Armadillo computes that length in int and
  overflows there (READ: `auxlib::eig_sym_dc()`).
- The dgelsd fallback above the oracle's n = 8: 40 singular or
  near-singular systems, n from 8 to 400, both triangles, are bitwise
  identical to `arma::solve()` with the same flag. The hand-computed
  workspace floors do not change the answer at those sizes.
- `dev/r_manual/R-exts.md`: "Fortran character strings" (the define before
  the first R header, `FCONE` once per character argument: 2, 1, 2, 3, 2,
  0, 0, 3, 3, 0, 1 for dgemm, dgemv, dsyrk, dsyevr, dsyevd, dgeqrf, dorgqr,
  dtrtrs, dtrcon, dgelsd, dgesdd, matching `FCLEN` in the installed
  headers), "Calling C from Fortran and vice versa" (`F77_CALL` only),
  "Portable C and C++ code" (a standard header for every name used, `std::`
  throughout, no R header inside `extern "C"`). The argument lists match
  the installed `R_ext/BLAS.h` and `R_ext/Lapack.h` in order and constness,
  dgelsd's non-const sizes included.
- The installed `R_ext/Lapack.h` marks routines "Added in R 4.4.0" and
  "Added in R 4.5.0" in comments. None of the eight used here carries such
  a mark: six are in the main alphabetical part, dgelsd and dgesdd in the
  block headed "Other double precision and double complex Lapack routines
  provided by libRlapack". That is a reading of the 4.6.1 header, not a
  check of R 4.3.0.
- `air format .` changed nothing tracked. The ASCII check prints
  `.github/workflows/pkgdown.yaml`, which this range does not touch; the
  nine files of the range are ASCII. `uv run dev/check_map_anchors.py`
  reports 80 problems, every one in `dev/r_manual/R-exts.md`, as the plan
  says. Run on this note it resolves all nine fragments and reports only
  the four review scripts, which git ignores.

Not checked: R 4.3.0 (left to the Linux images, as ruled); a second BLAS
or LAPACK; timing; whether `pca_imp()` can today reach the inputs of
findings 1 to 3 was read, not run.

## The independent read

Read in `dev/lib/arma/RcppArmadillo/include/armadillo_bits`, version 15.6.2,
then compared with `src/si_lapack.cpp` and with the READ notes.

| wrapper | what was read | agrees |
|---|---|---|
| `gemm`, `syrk` | the chain calls `arma::dgemm_` and `arma::dsyrk_` directly (`src/gram_ops.h`, `src/pca_linalg_utils.h`, `src/armaSVD.cpp`), block forms with a leading dimension larger than the block included; the wrapper passes every argument through | yes |
| `gemv`, `matmul` | `glue_times::apply()` and `apply_inplace_plus()`: empty operand, then one row (dgemv 'T' on B), then one column, then dsyrk for one object, then dgemm; `+=` and `-=` are one call with beta 1. `gemm::apply_blas_type()` and `gemv::apply_blas_type()` run their own loops for square operands up to 4 x 4. `syrk` uses its own loops for a vector or 48 cells or fewer | yes, with the two differences the header states |
| `syevr_range` | no Armadillo counterpart; the chain's own call in `src/eig_sym_sel.h` passes `jobz = 'V'`, `range = 'I'`, `abstol = 0`, `isuppz` of `2 * topk`, workspace the larger of the query and `26 n` and `10 n`, failure on `info != 0` or a wrong count | yes; finding 1 is in both |
| `syevd` | against `small_eig_desc()` in `src/lobpcg_warm.h`: query, then the solve with a buffer at least as long as the query asked; failure on `info != 0`; no finiteness check | yes |
| `qr_econ` | `auxlib::qr_econ()` and `auxlib::qr()`: tall input factors in `Q`, square and wide input in `R` with the first `m` columns copied to `Q`; one workspace, dgeqrf's optimum raised to `max(1, m, n)`; dorgqr with `k = min(m, n)` and not queried | yes |
| `solve_tri` | `glue_solve_tri_default::apply()`, `auxlib::solve_trimat_rcond()`, `auxlib::rcond_trimat()`, `auxlib::solve_approx_svd()`: dtrtrs, dtrcon in the 1-norm with `3 n` and `n` workspace, fallback when the solve failed or the estimate is below `DBL_EPSILON` or NaN, dgelsd on the triangle alone with `rcond = n * DBL_EPSILON`, non-finite refused in the fallback only | yes |
| `svd_values` | `auxlib::svd_dc()`: empty first, non-finite refused, dgesdd `jobz = 'N'`, `iwork` of `8 min(m, n)`, workspace `3 min + max(max, 7 min)` raised to the query, which is made from 1024 cells | yes |

Two things the read found that the notes do not say; neither changes a
result measured here:

- `auxlib::eig_sym_dc()`, the call the oracle compares `syevr_range` with,
  refuses a matrix with NaN or Inf in its upper triangle before LAPACK, and
  asks dsyevd for its workspace only from n = 32 (finding 8).
- Armadillo's dgelsd fallback takes `liwork` as the larger of the query's
  answer and a floor from ilaenv's SMLSIZ; the wrapper's floor uses the
  level count for SMLSIZ = 1, which is never smaller. Both are upper bounds
  on what dgelsd needs.

## Findings

### 1. `syevr_range` gives dsyevr a `w` shorter than LAPACK requires

The header says `w` receives `iu - il + 1` eigenvalues
(`src/si_lapack.h`: `w receives iu - il + 1 eigenvalues, ascending`), and
the wrapper passes the caller's `w` straight to dsyevr
(`src/si_lapack.cpp`: `isuppz.data(), work.data(), &lwork,`). LAPACK's own
description of dsyevr gives `W` the dimension `N`, with the first `M`
elements holding the answer (from the comment block of dsyevr in reference
LAPACK, which is not in this repository). For `range = 'I'` on a proper
subset dsyevr calls dstebz, which stores every eigenvalue it finds in its
search interval in `W` and then discards the ones outside `il..iu`. When
the range cuts a group of eigenvalues bisection cannot separate, it finds
more than `iu - il + 1`.

Measured (`dev/review_eigen_wrappers_probe.R`, group `w_length`): `w` is
given `n` cells, the cells past `iu - il + 1` start as a sentinel, and the
probe counts how many the call changed. Status was 0 and the returned
values agreed with `eigen()` to 4.6e-13 in every case.

| matrix | n | pairs | cells written past `w[count - 1]` |
|---|---|---|---|
| random Gram, top 1, 3 and 5; the oracle's two degenerate cases | 10 to 60 | no group cut | 0 |
| identity | 6 | 5..6 | 4, with either triangle |
| `diag(1,1,2,2,5,5)` | 6 | 6..6 | 1, with either triangle |
| `diag(1,1,2,2,5,5)` | 6 | 4..6 | 1, with either triangle |
| Gram of rank 4 | 10 | 5..10 | 3 (`'U'`), 0 (`'L'`) |
| Gram of rank 4 | 10 | 6..10 | 0 (`'U'`), 1 (`'L'`) |
| 30 rows of which 6 distinct, columns centred | 30 | 26..30, 25..30 | 0 |
| the same | 30 | 23..30 | 17 |
| exactly rank 3, 20 x 50 | 20 | 16..20 | 1 |

34 cells in all, counting each triangle's run. A caller that sizes `w` as the header says has that many
doubles written past its buffer. The last two rows are what `pca_imp()`
would pass when it is asked for more components than the data has rank.

The chain's present code has the same sizing (`eigvals.set_size(ws.topk);`
in `src/eig_sym_sel.h` as of `79cce0f`, gone since step 3 of Phase 2), so
this is not introduced by the wrapper; it is carried over and written into
a contract.
READ, not run through `pca_imp()`.

The oracle cannot see it: its pair allocates the same short buffer
(`dev/oracle_pairs.cpp`: `arma::vec w(count);`) and its one case with
repeated eigenvalues keeps the range off the repeat
(`dev/oracle_run.R`: `"repeated, top pair")`).

Recommendation: the wrapper owns an `n`-cell `W` and copies `iu - il + 1`
values out, so no call site can get it wrong; the oracle gets ranges that
cut a repeat, with guard cells after `w`.

### 2. `syevr_range` does not return on NaN when every eigenpair is asked for

Measured (`dev/review_eigen_wrappers_probe.R --group=nonfinite_syevr`, one
child process per case, 60 s limit; an ordinary call at these sizes returns
at once):

| matrix | pairs `n-1..n` | pairs `1..n` |
|---|---|---|
| NaN on the diagonal, n = 6 | -1004 (`err_count`) | did not return |
| NaN off the diagonal, n = 6 | -1004 | did not return |
| Inf on the diagonal, n = 6 | -1004 | 0, with 6 NaN eigenvalues |
| -Inf off the diagonal, n = 6 | -1004 | -1004 |
| all NaN, n = 6 | -1004 | -1004 |
| NaN on the diagonal, n = 60 | -1004 | did not return |
| NaN off the diagonal, n = 60 | -1004 | did not return |
| Inf on the diagonal, n = 60 | -1004 | -1004 |
| -Inf off the diagonal, n = 60 | -1004 | -1004 |
| all NaN, n = 60 | -1004 | -1004 |

A subset range fails cleanly in all ten cases. The full range, which
dsyevr solves by another routine, hangs in four and returns success with
NaN in one. A hang inside LAPACK cannot be interrupted from R.

The chain can ask for the full range: `k` is `ncp + 1` when regularized
(`src/armaSVD.cpp`: `const std::size_t k = regularized ? n_cp + 1 : n_cp;`)
and `ncp` may be as large as the eligible column count less one
(`R/pca_imp.R`: `max_ncp_cols <- n_elig - 1L`), so `k` can equal the order
of the Gram. READ, not run. Whether a NaN can reach the Gram in
`pca_imp()` is not measured; the chain's own dsyevr call has no check
either, so the wrapper is no worse than what it replaces.

The plan already lists the cure as an unruled Phase 4 proposal ("Refuse NaN
and Inf before every LAPACK call"). The measurement is the reason to rule
on it before `src/eig_sym_sel.h` moves to the wrapper: it is the only
failure found here that is a hang and not a status.

### 3. Non-finite input comes back as status 0 with NaN results

Measured (group `nonfinite_syevd`): `si_la::syevd` on the same five kinds
of matrix, both triangles, n = 6 and n = 60, 20 calls. Status 0 with NaN
eigenvalues in 12: all ten at n = 60 and two at n = 6. A positive `info`
(5 or 2) in the other eight. Every call returned.

| matrix | n = 6, `'U'` | n = 6, `'L'` | n = 60, both triangles |
|---|---|---|---|
| NaN on the diagonal | 5 | 5 | 0, 59 of 60 values NaN |
| NaN off the diagonal | 5 | 5 | 0, 59 of 60 values NaN |
| Inf on the diagonal | 0, 6 of 6 values NaN | 0, 6 of 6 values NaN | 0, 60 of 60 values NaN |
| -Inf off the diagonal | 2 | 5 | 0, 60 of 60 values NaN |
| all NaN | 5 | 5 | 0, 60 of 60 values NaN |

`syevr_range` did the same once among the 16 calls that returned (finding
2's table: Inf on the diagonal, n = 6, full range).

The header's convention is "0 is success; anything else is a failure", and
a call site written to it will take NaN eigenpairs as an answer. Same root
as finding 2 and same cure; the chain's present dsyevd call behaves the
same way (READ: `small_eig_desc()` tests only `info`).

### 4. The oracle passes a wrapper that calls another routine

Class "same" is for pairs that enter the same routine with the same
arguments, and its header says a difference there "is a finding in itself".
But the comparison is a tolerance
(`dev/oracle_run.R`: `same = list(abs = 1e-13, rel = 1e-12, rel_floor = 0),`)
and the bitwise count is printed only
(`dev/oracle_run.R`: `# bitwise identical. Bitwise identity is reported, never asserted.`).

Measured (`dev/review_eigen_wrappers_mutants.R`):

- `matmul_never_gemv`, every product sent to dgemm: "ORACLE PASSED", with
  matmul at 58 of 83 bitwise where the unbroken source has 83 of 83.
- `matmul_order`, the one-column test before the one-row test: "ORACLE
  PASSED". Group `order` of the probe: for a 1 x k row times a k x 1
  column the two dgemv routes differ in 87 of 200 random products at
  k = 5 and 184 of 200 at k = 1000, by at most 1.8e-15 and 1.3e-13. No
  tolerance that lets a BLAS repeat itself loosely can tell them apart.

Both sides are built in one process, the case `.claude/CLAUDE.md` names as
the one where `identical()` stays correct. Recommendation: class "same"
fails when a case is not bitwise identical; if a threaded BLAS is the
worry, that is a reason to say which cases may differ, not to accept all
of them.

### 5. `syevr_range` is compared only with another algorithm

`pair_syevr` compares dsyevr with `arma::eig_sym()`, which is dsyevd, at
1e-11 absolute and 1e-7 relative. Anything that moves the answer by less
passes. Measured: `syevr_abstol` (`abstol = 1e-12` in place of 0) passes
with no bitwise count moved; `work_min_only` (LAPACK's answer to every
workspace query ignored, the documented minimum used) passes too.

The exact reference exists: the chain's own dsyevr call. The pairs unit
includes Armadillo and not `R_ext/Lapack.h`, so it can declare `dsyevr_`
as `src/eig_sym_sel.h` does and call it with that file's arguments. The
pair is then class "same", and `abstol`, `range`, `isuppz` and the
workspace are pinned to the bit. The eig_sym comparison can stay as a
second, independent check.

### 6. No input lies near a rule the wrappers copy

Each of these changes a rule the wrapper takes from Armadillo, and each
passes with no bitwise count moved:

| mutant | rule moved |
|---|---|
| `solve_threshold_up`, `solve_threshold_down` | the fallback threshold, `DBL_EPSILON`, times 1000 and divided by 1000 |
| `solve_no_nan_test` | the NaN test on the estimate |
| `solve_sv_cut` | dgelsd's cut, `n * DBL_EPSILON`, made `DBL_EPSILON` |
| `svd_always_query` | the 1024-cell rule |
| `qr_dorgqr_query` | dorgqr's workspace, dgeqrf's optimum, made `max(m, n)` |
| `work_min_only` | every workspace, the query's answer, made the documented minimum |

The first three are catchable: no `solve_tri` case has an estimate within
three orders of magnitude of the threshold on either side, and none has a
singular value between `eps` and `n * eps` times the largest. The last
three may be equivalent on this LAPACK at any size the oracle could try;
if so, the rule is scaffolding the oracle cannot justify, which is worth
knowing before Phase 4 decides what to remove. The plan gives "a different
length can put LAPACK on another code path" as the reason for two of them;
no case shows a path that differs.

### 7. Failure paths with no case

- `err_count`: no case, and `syevr_no_count` (the check removed) passes.
  The self-test's `syevr_count` proves the oracle sees a failure invented
  on good input, not that it sees the check missing.
- `err_alloc`: no case in the oracle. The probe's `solve_tri` with both
  sizes at `INT_MAX` returns it.
- The leading-dimension checks of a transposed operand: the refusals try
  `lda`, `ldb` and `ldc` for `'N'` only. `gemm_lda_trans` and
  `syrk_lda_trans` (those checks removed for `'T'`) pass. The code under
  review does make the checks (probe, group `refusals`).
- Which code a failure carries: a failed case is recorded as "both
  failed" and nothing else is compared, so `solve_code` (`err_arg`
  returned for `err_nonfinite`) passes. The probe reads -1003 for the
  three non-finite refusals.
- Non-finite input to `syevd` and `syevr_range`, and a range that cuts
  equal eigenvalues: no case; findings 1 to 3 are what such cases show.
- A dgelsd fallback above n = 8: no case in the oracle. The probe's 40
  agree to the bit.

Two extra mutants were caught: `syevr_isuppz_half` (the process died,
exit 884, with no verdict line) and `svd_no_finite_check` (the Inf case;
the NaN case alone would pass, since dgesdd returns a negative `info` for
NaN and both sides then count as failed).

### 8. `auxlib::eig_sym_dc()` is left NOT READ

"What the wrappers reproduce" says of the oracle's reference for
`syevr_range`: "`auxlib::eig_sym_dc()` calls `lapack::syevd`; the rest of
it is NOT READ". The rest refuses a matrix with NaN or Inf in its upper
triangle and returns false. The oracle requires both sides to agree on
failure, and for non-finite input they do not: in the probe the oracle's
own pair gave wrapper status 0 against "Armadillo failed TRUE" (Inf on the
diagonal, n = 6, full range), and the pair cannot be called at all for the
inputs of finding 2. A non-finite case added to the oracle as it stands
would fail for that reason and not for a wrong wrapper.

## The extra mutants, all sixteen

`dev/review_eigen_wrappers_mutants.R dev/lib/arma <out dir>`. The unbroken
source passes with the bitwise counts above.

| mutant | what is wrong | oracle | bitwise counts that moved |
|---|---|---|---|
| `work_min_only` | every workspace is the documented minimum | passed | none |
| `syevr_abstol` | `abstol = 1e-12` | passed | none |
| `syevr_no_count` | the eigenpair count is not checked | passed | none |
| `syevr_isuppz_half` | `isuppz` at half its length | process died, exit 884 | - |
| `matmul_order` | one column tested before one row | passed | matmul_tiny 19 of 42 |
| `matmul_never_gemv` | every product to dgemm | passed | matmul 58 of 83, matmul_tiny 17 of 42 |
| `svd_always_query` | dgesdd asked below 1024 cells too | passed | none |
| `svd_no_finite_check` | NaN and Inf not refused | FAILED: svd_values (the Inf case) | svd_values 14 of 15 |
| `solve_code` | `err_arg` returned for `err_nonfinite` | passed | none |
| `solve_threshold_up` | fallback below `1000 * DBL_EPSILON` | passed | none |
| `solve_threshold_down` | fallback below `DBL_EPSILON / 1000` | passed | none |
| `solve_no_nan_test` | the estimate is not tested for NaN | passed | none |
| `solve_sv_cut` | dgelsd cuts at `DBL_EPSILON` | passed | none |
| `gemm_lda_trans` | `lda`, `ldb` unchecked for a transposed operand | passed | none |
| `syrk_lda_trans` | `lda` unchecked for `'T'` | passed | none |
| `qr_dorgqr_query` | dorgqr given `max(m, n)` | passed | none |
