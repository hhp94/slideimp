# pca_imp() structural map

A descriptive map of the `pca_imp()` call chain, from the R entry point down to BLAS/LAPACK.
Citations are `path:line`. When this file and the code disagree, the code is right.

## Call chain overview

```
R/pca_imp.R:362            pca_imp()                       user-facing entry
  |
  +- R/pca_imp.R:389       check_finite(obj)            -> src/mat_stats.cpp:198
  +- R/pca_imp.R:438       new_lobpcg_control()            R/pca_imp.R:168
  |    +- R/pca_imp.R:94   lobpcg_control()                exported constructor
  +- R/pca_imp.R:457,551   mat_miss()                   -> src/mat_stats.cpp:242,269
  +- R/pca_imp.R:468       col_vars()                   -> src/mat_stats.cpp:30
  +- R/pca_imp.R:581       for (i in seq_len(nb.init))     R-level restart loop
  |    |
  |    +- R/pca_imp.R:592  pca_imp_internal_cpp()
  |         +- R/RcppExports.R:4        .Call(`_slideimp_pca_imp_internal_cpp`)
  |         +- src/RcppExports.cpp:17   SEXP wrapper (BEGIN_RCPP/END_RCPP)
  |         +- src/armaSVD.cpp:306      pca_imp_internal_cpp()   C++ entry symbol
  |              |
  |              +- :463  pass 1: per-column NA counts
  |              +- :512  pass 2: centered moments + Xhat fill
  |              +- :631  initialize_missing_start()   armaSVD.cpp:177
  |              +- :671  gram_cache_init()            gram_ops.h:63  -> dsyrk_
  |              +- :745  EM outer loop
  |                   +- :754 impute_restandardize<>() armaSVD.cpp:32
  |                   +- :771 SVD_triplet()            svd_triplet.h:62
  |                   |    +- form_weighted_gram()     gram_ops.h:156 -> dsyrk_/dgemm_
  |                   |    +- hybrid_topk_eig()        hybrid_topk_eig.h:441
  |                   |    |    +- lobpcg_solve()      lobpcg_warm.h:367 -> dsyevd_
  |                   |    |    +- eig_sym_sel()       eig_sym_sel.h:129  -> dsyevr_
  |                   |    +- gemm_nn/gemm_tn()        pca_linalg_utils.h:85,122 -> dgemm_
  |                   +- :818 reconstruct_cols()       armaSVD.cpp:276 -> dgemm_
  |                   +- :825 missing_weighted_residual_ss()  armaSVD.cpp:236
  |                   +- :864 convergence test
  |              +- :907  postprocessing, imputed_values triplets
  |              +- :1010 Rcpp::List back to R
  |
  +- R/pca_imp.R:656       clamp
  +- R/pca_imp.R:667       write imputed triplets back into obj
  +- R/pca_imp.R:676       post_imp: mean_imp_col() on the columns that can
  |                        still hold NA, skipped when none  -> src/mat_stats.cpp:79
  +- R/pca_imp.R:688       new_slideimp_results()          R/utils.R:5
```

C++ files in the chain: 3 translation units (`src/armaSVD.cpp`, `src/RcppExports.cpp`,
`src/mat_stats.cpp`) and 9 headers (`svd_triplet.h`, `gram_ops.h`, `eig_sym_sel.h`,
`hybrid_topk_eig.h`, `lobpcg_warm.h`, `pca_linalg_utils.h`, `matrix_checks.h`,
`loc_timer.h`, `par_for.h`). Never touched: `find_windows.cpp`, `find_windows_flank.cpp`,
`impute_knn_brute.cpp`, `imputed_value.cpp/.h`, `sample_each_rep_cpp.cpp`.

## Layer 1: R entry point (R/pca_imp.R)

**Control-object constructors.** `resolve_clamp()` (`:10-55`) normalizes `clamp` to `NULL`
or an unnamed length-2 numeric; length 1 is a distinct error (`:15`) from other bad lengths
(`:23`), and `NA` bounds are rejected at `:31`. `lobpcg_control()` (`:94-112`) is the
exported constructor, defaulting to `warmup_iters = 10L`, `tol = 1e-9`, `maxiter = 20`; it
records which arguments were passed explicitly in an `explicit` attribute (`:95-99`,
attached `:110`) and stamps class `"slideimp_lobpcg_control"` (`:109`).
`new_lobpcg_control()` (`:168-257`) is the internal normalizer: `solver = "exact"`
short-circuits to `lobpcg_control(maxiter = 0L)` (`:177-178`), overriding anything the user
passed; `NULL` yields defaults (`:183-184`); a named list is validated against
`formals(lobpcg_control)` (`:193`, `:217-231`) then funneled through `do.call` (`:232`).
An object already carrying the class is re-validated rather than trusted (`:194-215`), by the
same `check_lobpcg_fields()` (`:131-152`) that `lobpcg_control()` calls (`:101`). Hand-built
control objects with no `explicit` attribute have every field treated as explicit
(`:236-239`). A final guard rejects `maxiter == 0` under `"auto"`/`"lobpcg"` (`:250-255`).

**Argument handling and gating.** `pca_imp()` is defined at `:362-706`.
`check_finite(obj)` (`:389`) is the missing-value gate for the whole chain: it rejects
`Inf`/`-Inf` anywhere and any all-`NA`/`NaN` column (`src/mat_stats.cpp:198`). The C++
kernel does no Inf scan of its own; the comment at `src/armaSVD.cpp:326-328` records that
this is deliberate, since the kernel is entered once per `nb.init` on an unchanged `obj`.
`ncp` is bounded by `ncol(obj) - 1L` (`:393`); `row.w` accepts a numeric vector, `"n_miss"`,
or `NULL` (`:399-410`); `seed` is guarded against `seed * (nb.init - 1)` exceeding int range
(`:418-424`); `.progress` becomes `trace_iter` (10 or 0) at `:454`, a C++-side print period.
`cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)` (`:457`), and a matrix with no missingness returns
unchanged at `:460-463` (as a bare matrix, not a `slideimp_results`). `eligible` (`:469-470`)
combines the `colmax` missing-rate ceiling with a non-degenerate-variance test from
`col_vars(obj)` (`:468`). Feasibility caps (`:474-504`) require `ncp` to be at most
`min(nrow(obj) - 2L, n_elig - 1L)`; both abort paths carry class `slideimp_infeasible`
(`:486`, `:502`), which `group_imp()` and `slide_imp()` catch, and a third such abort
(`:507-515`) fires when every column with missingness is ineligible.

**solver = "auto" policy (R half).** `:517-538` decides whether the C++ auto probe is worth
running at all. `k_eig` is `ncp + 1` for `"regularized"` and `ncp` for `"EM"` (`:522`);
`n_gram <- min(nrow(obj), n_elig)` (`:523`) mirrors the Gram dimension the backend will use.
`auto_force_exact` (`:525-526`) demotes `"auto"` to exact when `n_gram < 250` or
`k_eig / n_gram > 0.10`. Otherwise, and only when the user did not set `warmup_iters`
explicitly, `warmup_iters` is raised to `min(50, max(10, ceiling(1.5 * k_eig)))` (`:528-532`).
`solver_code` (`:534-538`) is the integer handed to C++: `0` exact, `1` lobpcg, `2` auto.

**Restart loop and result assembly.** `:581-635`. Each iteration re-seeds R's RNG with
`seed * (i - 1L)` (`:583`, matching `missMDA`) and calls the kernel with `init = 0L` for the
first restart and `init = i` afterwards (`:599`); `init == 0` means mean imputation (zeros
in centered space), anything else a Gaussian random start. After the first call under
`solver = "auto"`, the returned `solver_chosen` code is collapsed to a forced `0`/`1` and
locked in `locked_solver` (`:612-620`) so later restarts skip the probe; the restart with
the lowest `mse` wins (`:622-634`). `clamp` is applied to the value column of the triplet
matrix at `:656-661` (imputed values only), and `:663-667` writes the triplets back through
a two-column index matrix, the indices from C++ already being 1-based original-matrix
indices. `post_imp = TRUE` collects the columns that can still hold NA -- ineligible
columns that had any, plus the columns of NA rows in the triplet -- and runs
`mean_imp_col(obj, subset = na_cols)` only when that set is non-empty (`:676-685`), so
an all-eligible input pays no extra copy here. `new_slideimp_results()`
(`R/utils.R:5`) forces `has_remaining_na` before attaching the class (a classed object
routes `anyNA()` through `any(is.na(x))`, an `n x p` logical), then sets class and
standard attributes; `:696-704` adds `solver_requested`, `solver_chosen`, `solver_probed`,
`n_iter`, `criterion_final`, `converged`, `n_exact`, `n_lobpcg_ok`, `n_lobpcg_bad`.
`slide_imp()` pins every later window to `attr(x, "solver_chosen")`, but only from a window
whose `solver_probed` is `TRUE`. A size-demoted window, or one whose probe never finished,
leaves the next window on `"auto"` (`R/slide_imp.R:511-537`).

## Layer 2: the Rcpp boundary

`R/RcppExports.R:4-6` is the generated shim, calling
`.Call(`_slideimp_pca_imp_internal_cpp`, ...)` with 16 arguments.
`src/RcppExports.cpp:16` declares the C++ signature; `:17-40` is the `RcppExport` SEXP
wrapper, converting each argument through `Rcpp::traits::input_parameter<>` (`:21-36`). An
`Rcpp::RNGScope` (`:20`) makes `Rcpp::rnorm()` inside the kernel use R's RNG stream.
Registration is `src/RcppExports.cpp:192` in `CallEntries`, with `R_init_slideimp` at
`:208-211` and `R_useDynamicSymbols(dll, FALSE)`. The C++ entry symbol is
`pca_imp_internal_cpp` (`src/armaSVD.cpp:306`, marked `// [[Rcpp::export]]` at `:305`).

`obj` arrives as `const arma::mat&`, so Rcpp constructs a read-only Armadillo view over the
SEXP rather than copying. `row_w` is taken **by value** (`src/armaSVD.cpp:316`) because the
kernel normalizes it in place (`:460`).

## Layer 3: kernel setup (src/armaSVD.cpp:306-743)

- Argument sanity (`:336-377`): `maxiter >= 1`, `miniter <= maxiter`, `ncp >= 1`, `solver`
  in 0..2, non-empty and in-bounds `eligible_idx`, `ncp < min(nrX, n_elig)`. `:380-383`
  guards `nrX * n_elig` overflow. `:389-431` guards the `blas_int` range: problem
  dimensions, the `26 * n` dsyevr workspace minimum (`:400`), and the `~3k` Rayleigh-Ritz
  dsyevd workspace LOBPCG would need (`:410-419`).
- Orientation (`:385-387`): `tall = (nrX >= n_elig)`; `k = regularized ? ncp + 1 : ncp`
  eigenpairs are requested; `n_gram = tall ? n_elig : nrX`, so the eigenproblem is always
  over the smaller dimension.
- Pass 1 (`:463-475`) counts NAs per eligible column. `:477-490` splits eligible columns
  into fixed (no NA) and changing, building `perm` so that after permutation the changing
  columns are exactly positions `[n_fixed, n_elig)`. Every later loop derives its local
  index as `n_fixed + ci`, which is why no column-index vector is threaded through.
- `:496-502` builds a CSR-style pair: `miss_rows_offsets` (length `n_mc + 1`) and
  `miss_rows_flat` (length `n_missing`). `for_each_missing_in_col()` (`:12-24`) is the only
  iterator over them.
- Pass 2 (`:512-574`) fills `Xhat` (`:504`) and accumulates per column the observed weight
  `denom`, the reference mean `mu0`, and centered moments `obs_sum_c` / `obs_sumsq_c`.
  Centering against `mu0` before squaring is what keeps the later variance updates away from
  `E[x^2] - mu^2` cancellation (`:551-552`). `:578-585` rejects any eligible column with
  zero observed weight. `:590-601` computes `mean_p` and, when `scale`, the per-column scale
  `et` floored at `PCA_TOL`; `:603-629` centers (and scales) `Xhat` in place.
  `initialize_missing_start()` (`:177-234`, called at `:631`) allocates `fittedX` once and
  writes starting values at missing positions only: zeros for `init == 0`, `Rcpp::rnorm`
  draws otherwise (`:215`).

Workspaces allocated once, outside the loop (`:639-684`): `sqrt_row_w` / `inv_sqrt_row_w`
(the weighted space is `diag(sqrt(w)) * Xhat`, and `U` comes back in it); `GramWorkspace`
holding fixed dsyrk/dgemm dimension arguments (`gram_ops.h:14-34`); `EigSymWorkspace`, which
runs the dsyevr workspace query once (`eig_sym_sel.h:63-120`); `AA_NxN` at `n_gram x n_gram`
and `X_work` at `nrX x n_elig` in the tall case only (`:662-666`); `GramCache`, precomputing
the fixed-column Gram contribution (`:670-674`); and `HybridEigContext` (`:676-684`),
carrying the solver mode, warmup count, LOBPCG options and persistent LOBPCG state across
all outer iterations.

## Layer 4: the EM outer loop (src/armaSVD.cpp:745-902)

One iteration:

1. `Rcpp::checkUserInterrupt()` every 5 iterations (`:747-750`).
2. `impute_restandardize<Scale>()` (`:32-175`, called at `:752-767`) -- the fused heart of
   the loop. It copies the previous fit into `Xhat` at missing positions, updates each
   changing column's mean (and variance when scaling) incrementally from `mu0` plus centered
   contributions, rescales the column in one pass, and simultaneously refreshes the
   sqrt-weighted copy of that column that the Gram step consumes (`:131-138`, `:155-162`).
   Only changing columns are touched.
3. `SVD_triplet()` (`:771`) -- see Layer 5.
4. Shrinkage (`:778-804`): `tail` is the clamped residual energy
   `trace - sum of top ncp squared singular values` (`:783`); under `regularized`, `sigma2`
   is that tail scaled by `scale_factor / denom_sigma`, multiplied by `coeff_ridge`, and
   capped at the `(ncp+1)`-th squared singular value (`:786-787`);
   `lambda_shrinked[i] = v - sigma2 / v` (`:797`).
5. Reconstruction (`:806-820`). `scale_rows_cols_inplace()` (`pca_linalg_utils.h:42`) fuses
   the row de-weighting of `U` with the shrunk column scaling, mutating `U` in place;
   `reconstruct_cols()` (`:276-303`) then does one `dgemm_` writing only the changing block
   of `fittedX` (`:818`). The fixed block is never read inside the loop.
6. Objective (`:822-836`): `objective_full` minus the weighted residual sum of squares over
   missing cells only (`missing_weighted_residual_ss()`, `:236-272`), floored at 0.

**Convergence** (`:864-887`): the criterion is `abs(1 - objective / old)` (`:864`),
reproducing `missMDA::imputePCA()`. It is tested only once `nb_iter >= miniter` (`:868`),
and stopping requires `criterion < threshold || objective < threshold` (`:870`). Hitting
`maxiter` first sets `converged = false` and nothing else (`:876-887`). The kernel stays silent
because a warning raised from C++ can longjmp past the frame's destructors, so `pca_imp()` warns
once, through `cli::cli_warn()`, for the winning restart only (`R/pca_imp.R:648-654`). The loop
is `for (arma::uword nb_iter = 1;; ++nb_iter)` (`:745`), whose only exits are at `:896-901`.
`trace_iter > 0` prints on
iteration 1, every `trace_iter` iterations, and on the last (`:889-894`); print precision is
restored by the `PrecGuard` RAII struct (`:699-711`).

**Postprocessing** (`:907-976`): `reconstruct_cols(U, V, fittedX, 0, n_fixed)` (`:911`)
completes the fixed block the loop skipped, purely so the SSE covers the whole matrix. The
SSE loop (`:914-946`) computes the full residual in original units (times `et^2` per column)
then subtracts the missing-cell contribution, giving an observed-cells-only `mse` (`:948`),
the `nb.init` selection statistic. `imputed_values` (`:950-974`) is an `n_missing x 3`
matrix of `(row_1based, orig_col_1based, value)`, where the value comes from the final
`Xhat` at the missing position, de-standardized -- not from `fittedX` (`:960-962`). The
ordinary-build return list (`:1010-1019`) carries `imputed_values`, `mse`, `n_iter`,
`criterion_final`, `converged`, `solver_chosen`, `n_exact`, `n_lobpcg_ok`, `n_lobpcg_bad`.
`solver_chosen` comes from `HybridEigContext::solver_chosen_code()`
(`hybrid_topk_eig.h:239-262`): 0 forced exact, 1 forced lobpcg, 2 auto never decided, 3 auto
chose exact, 4 auto chose lobpcg. Under `solver = "auto"` only (`R/pca_imp.R:612`), `:617`
collapses 1 and 4 to "lobpcg" and everything else to "exact"; the string is set at `:686`.

## Layer 5: Gram formation and factor recovery

`SVD_triplet()` (`src/svd_triplet.h:62-125`) never forms an SVD directly. It builds the
weighted Gram into `AA_NxN` via `form_weighted_gram()` (`:83`), takes
`trace_val = arma::trace(AA_NxN)` (`:86`) for the tail energy, calls `hybrid_topk_eig()`
(`:89`, failure being an `Rcpp::stop` at `:91`), converts eigenvalues to singular values and
inverses with a `PCA_TOL` floor (`:95-111`), and recovers the other factor with one `dgemm_`
(`:113-124`): tall gives eigenvectors as `V` and `U = X_work * Vn` (`gemm_nn`, `:117`), wide
gives them as `U` and `V = Xhat' * Vn` (`gemm_tn`, `:122`). `split_eigvec_block()` (`:15-54`)
writes the plain and `d_inv`-scaled copies in a single pass.

`form_weighted_gram()` (`src/gram_ops.h:156-265`) has four branches:

| case | how the Gram is built |
| --- | --- |
| tall, cache active | copy cached fixed block (`:179`), `dgemm_` cross block (`:191`), `dsyrk_` changing-changing block (`:207`) |
| tall, no cache | one `dsyrk_` over all of `X_work` (`:218`) |
| wide, cache active | copy cached fixed Gram (`:226`), `dsyrk_` with `beta = 1` adding the changing columns (`:238`) |
| wide, no cache | `dsyrk_` on `Xhat` (`:248`), then scale entries by `sw[i] * sw[j]` (`:250-262`) |

Only the upper triangle is written, matching dsyevr's `uplo = 'U'` (`eig_sym_sel.h:42`). The
LOBPCG path mirrors it down first (`hybrid_topk_eig.h:382`).

## Layer 6: the eigensolver

**Routing** (`src/hybrid_topk_eig.h`). `hybrid_topk_eig()` (`:441-458`) is a timing wrapper:
while the auto probe is undecided (`auto_timing_active()`, `:121-125`) it brackets the real
call with a `steady_clock` measurement and feeds it to `record_auto_probe_time()`
(`:182-202`). `hybrid_topk_eig_impl()` (`:265-435`) picks the path at `:290-373`:

- `Exact` -- always dsyevr.
- `Lobpcg` -- dsyevr while `outer_iter < warmup_iters` (`:290`, `:307-310`) or while there
  is no seed, then LOBPCG.
- `Auto` (`:320-372`) -- exact for at least `max(warmup_iters, auto_min_exact_iter)`
  iterations (`:344`, target at `:116-119`), then `HYB_AUTO_N_PROBE_ITER_DEFAULT = 5` LOBPCG
  probes (`:354-357`), then decides. `maybe_decide_auto()` (`:142-180`) picks LOBPCG only
  when its mean time is below `(1 - 0.10) * exact_mean` (`:163`, margin at `:60`); otherwise
  it picks exact and drops the LOBPCG state so exact mode carries no seed-maintenance
  overhead (`:127-135`).

On LOBPCG success (`:386-403`) the eigenvalues are moved out, `eigvecs` is copied from the
committed state, and `refresh_momentum_P()` rebuilds the search direction from the two most
recent eigenblocks. On failure `n_lobpcg_bad` increments (`:405`) and control falls through
to `eig_sym_sel()` (`:408`); after any exact solve, `seed_lobpcg_state()` re-seeds the warm
start when `keep_lobpcg_seed_state()` allows (`:424-432`).

**Exact path** (`src/eig_sym_sel.h`). `dsyevr_` is declared by hand (`:21-35`) rather than
pulled from `<R_ext/Lapack.h>`; the comment at `:10-20` records that the R header clashes
with Armadillo's own declarations. Hidden Fortran character-length arguments are passed
explicitly (the trailing `1, 1, 1`). `EigSymWorkspace::init()` (`:63-120`) runs the
`lwork = -1` query once and floors the result at the documented `26 * n` / `10 * n` minima
(`:105-114`). `eig_sym_sel()` (`:129-177`) requests eigenpairs `il..iu` with `range = 'I'`
(`:41`, `:72-73`), returning exactly the top `k`. dsyevr returns ascending order, so
`:162-173` reverses eigenvalues and swaps whole eigenvector columns in place, then
`canonicalize_signs()` (`pca_linalg_utils.h:15`) applies the sklearn `svd_flip` convention.
`A` is destroyed by the call (`:126-127`), which is safe because the Gram is rebuilt every
outer iteration.

**LOBPCG path** (`src/lobpcg_warm.h`). This is a **warm-start-only** solver:
`lobpcg_solve()` (`:367`) hard-errors via `Rcpp::stop` on an unseeded state (`:375-376`), so
callers must check `LOBPCGState::seeded()` (`:73-77`) first. The state carries `X` (`n x k`,
orthonormal, sign-canonicalized approximate top-k eigenvectors) and `P` (`n x k`, the search
direction). The header's "P invariant" note (`:21-25`) states that the committed `P` is
deliberately *not* orthonormal and not orthogonal to `X`; entry-time cleanup (`:390-391`)
re-projects and re-orthonormalizes it. Both matvecs are recomputed fresh at entry
(`:397-402`) because `A` changes between outer calls. The initial Rayleigh-Ritz is over the
joint basis `[X, P]` (`rr_two_block()`, `:275-310`), falling back to `X` alone when `P` is
unusable (`:434-448`).

The iteration (`:459-620`) is standard three-block LOBPCG over `[X, R_active, P_active]`.
Besides convergence it exits on `maxiter` (`:459`), a divergence guard at relative residual
> 10 (`:474-484`), a stall counter over `stall_window = 5` (`:486-496`), residual collapse
after QR (`:507-513`), `P` collapse (`:547-552`, `:562-567`), and `LOBPCG_RR_FAILED` when
`rr_two_block` or `small_eig_desc` fails (`:519-525`, `:595-600`). `MATVEC_REFRESH = 20` (`:452`)
periodically recomputes `AX`/`AP` from scratch (`:529-533`, `:615-619`) to shed accumulated
rounding. State is committed only on success (`:625-631`). Small Rayleigh-Ritz problems go
through `small_eig_desc()` (`:203-250`), which calls `dsyevd_` twice (query `:221`, solve
`:239`) with buffers cached on the `LOBPCGState`. Warm-start construction is two entry
points over one helper, `build_momentum_P()` (`:145-174`), which forms
`P = X_prev - X_curr (X_curr' X_prev)`: `seed_lobpcg_state()` (`:318-337`) after an exact
solve, with no norm floor; `refresh_momentum_P()` (`:345-359`) after a LOBPCG success, with
a floor of `1e-10 * sqrt(k)` below which the post-RR `P` is kept instead.

**Where BLAS/LAPACK is entered.**

| symbol | site | purpose |
| --- | --- | --- |
| `dsyrk_` | `gram_ops.h:54`, `:207`, `:238` | Gram formation |
| `dgemm_` | `gram_ops.h:191` | fixed-changing Gram cross block |
| `dgemm_` | `pca_linalg_utils.h:112`, `:151` | `gemm_nn` / `gemm_tn` factor recovery |
| `dgemm_` | `armaSVD.cpp:296` | `reconstruct_cols`, strided column block |
| `dsyevr_` | `eig_sym_sel.h:86`, `:145` | workspace query, exact top-k eigen |
| `dsyevd_` | `lobpcg_warm.h:221`, `:239` | LOBPCG Rayleigh-Ritz |
| Armadillo | `lobpcg_warm.h:91` (`qr_econ`), `:192` (`solve`), `:397` (`A * X`) | LOBPCG orthonormalization and matvecs |

The exact path enters LAPACK once per outer iteration with an `n_gram x n_gram` problem. The
LOBPCG path never touches `n_gram`-sized LAPACK: its dense work is `n_gram x k` matvecs and
QRs plus a `3k x 3k` dsyevd.

## Data structures and ownership

- `obj` crosses as `const arma::mat&` and is never mutated in C++; the write-back happens in
  R (`R/pca_imp.R:667`). `row_w` crosses by value (`src/armaSVD.cpp:316`) and is normalized
  in place (`:460`).
- `Xhat` (`:504`) is the working copy -- permuted column order, centered, scaled, mutated in
  place every outer iteration. `fittedX` is allocated once (`:188`) and overwritten by dgemm
  each iteration. `U` is mutated in place by `scale_rows_cols_inplace()` (`:811`), so after
  the reconstruct step it no longer holds the orthonormal left factor.
- `AA_NxN` is destroyed by every eigen call: dsyevr overwrites it (`eig_sym_sel.h:126-127`),
  and the LOBPCG path mirrors its upper triangle down first (`hybrid_topk_eig.h:382`).
- `HybridEigContext` carries state across outer iterations: `lobpcg_state` (including its
  dsyevd `work`/`iwork` buffers), `X_prev`, `auto_choice`, the auto counters and timing sums,
  and `n_exact`/`n_lobpcg_ok`/`n_lobpcg_bad` (`hybrid_topk_eig.h:65-87`). When auto picks
  exact, only `lobpcg_state.X`, `lobpcg_state.P` and `X_prev` are reset (`:132-134`).
- Return crosses as an `Rcpp::List` (`src/armaSVD.cpp:1010`), copying the triplet matrix
  into R. Only missing cells are transported, never the full imputed matrix. The
  `mat_stats.cpp` entry points instead return by value: `mean_imp_col()` an `n x p` matrix
  (`:79`, `:120`) and `col_vars()` a length-`p` row vector (`:30`), each copied into a SEXP.
- `PCA_TOL = 1e-15` (`src/pca_linalg_utils.h:9`) is the shared floor in `armaSVD.cpp` and
  `svd_triplet.h`: scale floors (`armaSVD.cpp:124`, `:600`), inverse row-weight floor
  (`:644`), shrinkage guard (`:797`), singular-value inverse guard (`svd_triplet.h:109`).
  LOBPCG has its own floors (`lobpcg_warm.h:40`, `:350-351`), and R applies a `1e-10`
  row-weight floor (`R/pca_imp.R:405`, `:553`) and a `.Machine$double.eps` variance floor (`:470`).
- `src/matrix_checks.h` is not in `armaSVD.cpp`'s include chain. It is included by
  `src/mat_stats.cpp:4` and `src/imputed_value.h:9`, and its only function `stop_on_inf()`
  (`matrix_checks.h:8-29`) is reached from `src/mat_stats.cpp:193` via `check_inf()`, which `pca_imp()` hits
  indirectly through
  `mean_imp_col()` (`R/mean_imp_col.R:38`), and only when `post_imp` has columns left
  to fill.

## Compile-time switches

**`-DARMA_64BIT_WORD=1`** -- set in `src/Makevars:1` and `src/Makevars.win:1`. Makes
`arma::uword` 64 bits; without it a matrix past `2^32 - 1` elements wraps silently.
`arma_uword_bytes()` (`src/mat_stats.cpp:181-184`) exists only so the test suite can assert
the flag survived.

**`-DLOC_TIMER`** (`src/loc_timer.h`) -- when defined (`:3`) the header pulls in
`rcpptimer.h` and the `LOC_TIMER_*` macros become real: `LOC_TIMER_OBJ` declares an
`Rcpp::Timer`, `LOC_TIC`/`LOC_TOC` bracket a phase, and `LOC_TIMER_PARAM`/`LOC_TIMER_ARG`
thread the timer through signatures. When undefined, everything collapses to `((void)0)` and
the parameter macros expand to nothing (`:26-34`) -- which is why `SVD_triplet()`'s
signature ends with `LOC_TIMER_PARAM(timer)` (`svd_triplet.h:80`) and its call site ends
with `LOC_TIMER_ARG(pca_imp_gram)` (`armaSVD.cpp:775`). Timed phases:
`pca_imp_internal_cpp_total` (`armaSVD.cpp:325`), `pass1_count_missing` (`:463`),
`pass2_scan` (`:512`), `gram_cache_init` (`:669`), `restandardize` (`:751`), `svd` (`:770`),
`post_svd` (`:778`), `reconstruct` (`:806`), `objective` (`:822`), `postprocessing`
(`:907`), and inside `SVD_triplet` the `form_gram`, `eig` and `recover_factors` phases
(`svd_triplet.h:82`, `:88`, `:113`). `rcpptimer` is deliberately absent from `LinkingTo`, so
its include directory is added by hand in `load_all1()` (`R/dev-utils.R:11-18`).

**`-DPCA_IMP_DIAGNOSTICS`** (`src/armaSVD.cpp`) -- defaulted to 0 at `:1-3`, set by
`load_all1(debug = TRUE)` (`R/dev-utils.R:20`). It gates three things:
`hyb_ctx.init_logs(maxiter)` (`:682-684`), which allocates the five per-iteration log
vectors in `HybridEigContext` (`hybrid_topk_eig.h:97-105`) -- without it they stay empty and
`log_at()` (`hybrid_topk_eig.h:271-286`) is a no-op because its bounds check fails; the
history buffers `eigval_hist`, `obj_hist`, `subspace_cos_hist` (`:738-743`) and their update
block (`:838-854`), which adds an `arma::svd` of `V_prev' * V_new` per iteration to track
subspace rotation; and a wider return list (`:985-1008`). The counters `n_exact`,
`n_lobpcg_ok`, `n_lobpcg_bad` are always maintained and returned; only the per-iteration
logs are gated.

## Parallelism

- **Inside the PCA kernel: none.** `src/armaSVD.cpp` does not include `RcppThread.h`, and
  neither does any header it pulls in. All intra-call parallelism comes from whatever BLAS
  the R installation links.
- **RcppThread** is reached through `par_for()` (`src/par_for.h:33-55`), which both
  `mat_stats.cpp` kernels call (`:39` `col_vars_internal`, `:127` `mean_imp_col_internal`),
  as does the K-NN kernel (`src/impute_knn_brute.cpp:482`, not in this chain). With
  `n_threads <= 1` it runs the body inline and never starts an RcppThread pool
  (`par_for.h:40-47`), and `pca_imp()` reaches both kernels with the default `cores = 1`
  (`R/col_vars.R:37`, `R/mean_imp_col.R:34`). `src/Makevars:3` links `RcppThread::LdFlags()`; `src/Makevars.win:3`
  does not.
- **mirai** is used from R only, one level above `pca_imp()`: `R/group_imp.R:901`
  (`mirai::mirai_map` over groups) and `R/tune_imp.R:948` (over tuning replicates). Both
  share the input matrix through `bigmemory` shared memory rather than serializing it
  (`group_imp.R:831-839`, `tune_imp.R:894-895`) and ship the worker body as a `.crate()`
  closure (`R/crate.R:35`). `slide_imp()` runs its windows sequentially
  (`R/slide_imp.R:450`).
- **BLAS threading** is pinned inside the mirai workers, not in the kernel:
  `RhpcBLASctl::blas_set_num_threads(1)` and `RhpcBLASctl::omp_set_num_threads(1)` at
  `R/group_imp.R:851-852` and `R/tune_imp.R:906-907`, gated by `pin_blas`.
  `check_pin_blas()` (`R/utils.R:120-131`) errors when `pin_blas = TRUE` without
  `RhpcBLASctl` installed, and otherwise emits a tip when more than one BLAS thread is
  detected.
