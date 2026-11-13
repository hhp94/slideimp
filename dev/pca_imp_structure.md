# pca_imp() structural map

A descriptive map of the `pca_imp()` call chain, from the R entry point down to BLAS/LAPACK.
Code is cited by file and by a fragment of one line that occurs exactly once in that file,
never by line number; for example (`R/pca_imp.R`: `pca_imp <- function(`) is the definition
of `pca_imp()`. `uv run dev/check_map_anchors.py` checks every citation, and its header gives
the citation forms. When this file and the code disagree, the code is right.

## Call chain overview

```
R/pca_imp.R                pca_imp()                       user-facing entry
  |
  +- R/pca_imp.R           check_finite(obj)            -> src/mat_stats.cpp
  +- R/pca_imp.R           new_lobpcg_control()            R/pca_imp.R
  |    +- R/pca_imp.R      lobpcg_control()                exported constructor
  +- R/pca_imp.R           mat_miss()                   -> src/mat_stats.cpp
  |    +- R/pca_imp.R      no missing value: return pca_imp_results()
  +- R/pca_imp.R           col_vars()                   -> src/mat_stats.cpp
  +- R/pca_imp.R           for (i in seq_len(nb.init))     R-level restart loop
  |    |
  |    +- R/pca_imp.R      pca_imp_internal_cpp()
  |         +- R/RcppExports.R          .Call(`_slideimp_pca_imp_internal_cpp`)
  |         +- src/RcppExports.cpp      SEXP wrapper (BEGIN_RCPP/END_RCPP)
  |         +- src/armaSVD.cpp          pca_imp_internal_cpp()   C++ entry symbol
  |              |
  |              +- pass 1: per-column NA counts
  |              +- pass 2: centered moments + Xhat fill
  |              +- initialize_missing_start()         src/armaSVD.cpp
  |              +- gram_cache_init()                  src/gram_ops.h -> dsyrk_
  |              +- EM outer loop
  |                   +- impute_restandardize<>()      src/armaSVD.cpp
  |                   +- SVD_triplet()                 src/svd_triplet.h
  |                   |    +- form_weighted_gram()     src/gram_ops.h -> dsyrk_/dgemm_
  |                   |    +- hybrid_topk_eig()        src/hybrid_topk_eig.h
  |                   |    |    +- lobpcg_solve()      src/lobpcg_warm.h -> dsyevd_
  |                   |    |    +- eig_sym_sel()       src/eig_sym_sel.h -> dsyevr_
  |                   |    +- gemm_nn/gemm_tn_block()  src/pca_linalg_utils.h -> dgemm_
  |                   +- reconstruct_cols()            src/armaSVD.cpp -> dgemm_
  |                   +- missing_weighted_residual_ss()  src/armaSVD.cpp
  |                   +- convergence test
  |              +- postprocessing, imputed_values triplets
  |              +- Rcpp::List back to R
  |
  +- R/pca_imp.R           clamp
  +- R/pca_imp.R           write imputed triplets back into obj
  +- R/pca_imp.R           post_imp: mean_imp_col() on the columns that can
  |                        still hold NA, skipped when none  -> src/mat_stats.cpp
  +- R/pca_imp.R           pca_imp_results()
       +- R/utils.R        new_slideimp_results()
```

C++ files in the chain: 3 translation units (`src/armaSVD.cpp`, `src/RcppExports.cpp`,
`src/mat_stats.cpp`) and 9 headers (`svd_triplet.h`, `gram_ops.h`, `eig_sym_sel.h`,
`hybrid_topk_eig.h`, `lobpcg_warm.h`, `pca_linalg_utils.h`, `matrix_checks.h`,
`loc_timer.h`, `par_for.h`). Never touched: `find_windows.cpp`, `find_windows_flank.cpp`,
`impute_knn_brute.cpp`, `imputed_value.cpp/.h`, `sample_each_rep_cpp.cpp`.

## Layer 1: R entry point (`R/pca_imp.R`)

**Control-object constructors.** `resolve_clamp()` normalizes `clamp` to `NULL`
or an unnamed length-2 numeric; length 1 is a distinct error from other bad lengths,
and `NA` bounds are rejected. `lobpcg_control()` is the
exported constructor, defaulting to `warmup_iters = 10L`, `tol = 1e-9`, `maxiter = 20`; it
records which arguments were passed explicitly in an `explicit` attribute
and stamps class `"slideimp_lobpcg_control"`.
`new_lobpcg_control()` is the internal normalizer: `solver = "exact"`
short-circuits to `lobpcg_control(maxiter = 0L)`, overriding anything the user
passed; `NULL` yields defaults; a named list is validated against
`formals(lobpcg_control)` then funneled through `do.call`.
An object already carrying the class is re-validated rather than trusted
(`R/pca_imp.R`: `if (inherits(x, "slideimp_lobpcg_control")) {`), by the
same `check_lobpcg_fields()` that `lobpcg_control()` calls. Hand-built
control objects with no `explicit` attribute have every field treated as explicit
(`R/pca_imp.R`: `if (is.null(explicit)) {`). A final guard rejects `maxiter == 0` under
`"auto"`/`"lobpcg"`.

**Argument handling and gating.** `pca_imp()` is defined at
(`R/pca_imp.R`: `pca_imp <- function(`).
`check_finite(obj)` is the missing-value gate for the whole chain: it rejects
`Inf`/`-Inf` anywhere and any all-`NA`/`NaN` column
(`src/mat_stats.cpp`: `void check_finite(`). The C++ kernel does no Inf scan of its own; the
comment (`src/armaSVD.cpp`: `No Inf scan here.`) records that
this is deliberate, since the kernel is entered once per `nb.init` on an unchanged `obj`.
`ncp` is bounded by `ncol(obj) - 1L`; `row.w` accepts a numeric vector, `"n_miss"`,
or `NULL`; `seed` is guarded against `seed * (nb.init - 1)` exceeding int range;
`.progress` becomes `trace_iter` (10 or 0), a C++-side print period.
`cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)`, and a matrix with no missingness returns
unchanged, through the same `pca_imp_results()` as the main path (`R/pca_imp.R`:
`return(pca_imp_results(obj, post_imp = post_imp, solver = solver))`), so the result
has the same class and attributes either way. `eligible`
combines the `colmax` missing-rate ceiling with a non-degenerate-variance test from
`col_vars(obj)`. Feasibility caps (`R/pca_imp.R`: `max_ncp_rows <- nrow(obj) - 2L`)
require `ncp` to be at most
`min(nrow(obj) - 2L, n_elig - 1L)`; both abort paths carry class `slideimp_infeasible`,
which `group_imp()` and `slide_imp()` catch, and a third such abort
(`R/pca_imp.R`: `if (!any(cmiss[eligible] > 0L)) {`) fires when every column with
missingness is ineligible.

**solver = "auto" policy (R half).** The block (`R/pca_imp.R`: `# solver = "auto" policy.`)
decides whether the C++ auto probe is worth running at all.
`k_eig` is `ncp + 1` for `"regularized"` and `ncp` for `"EM"`;
`n_gram <- min(nrow(obj), n_elig)` mirrors the Gram dimension the backend will use.
`auto_force_exact` demotes `"auto"` to exact when `n_gram < 250` or
`k_eig / n_gram > 0.10`. Otherwise, and only when the user did not set `warmup_iters`
explicitly, `warmup_iters` is raised to `min(50, max(10, ceiling(1.5 * k_eig)))`.
`solver_code` is the integer handed to C++: `0` exact, `1` lobpcg, `2` auto.

**Restart loop and result assembly** (`R/pca_imp.R`: `for (i in seq_len(nb.init)) {`). Each
iteration re-seeds R's RNG with
`seed * (i - 1L)` (matching `missMDA`) and calls the kernel with `init = 0L` for the
first restart and `init = i` afterwards; `init == 0` means mean imputation (zeros
in centered space), anything else a Gaussian random start. After the first call under
`solver = "auto"`, the returned `solver_chosen` code is collapsed to a forced `0`/`1` and
locked in `locked_solver` so later restarts skip the probe; the restart with
the lowest `mse` wins. `clamp` is applied to the value column of the triplet
matrix (imputed values only), and the block (`R/pca_imp.R`: `imp_indices <- cbind(`) writes
the triplets back through
a two-column index matrix, the indices from C++ already being 1-based original-matrix
indices. `post_imp = TRUE` collects the columns that can still hold NA -- ineligible
columns that had any, plus the columns of NA rows in the triplet -- and runs
`mean_imp_col(obj, subset = na_cols)` only when that set is non-empty, so
an all-eligible input pays no extra copy here. `pca_imp_results()`
(`R/pca_imp.R`: `pca_imp_results <- function(`) builds the result on both paths. It calls
`new_slideimp_results()` (`R/utils.R`: `new_slideimp_results <- function(`), which strips
its input to a bare matrix, then sets class and standard attributes, and adds
`solver_requested`, `solver_chosen`, `solver_probed`, `n_iter`, `criterion_final`,
`converged`, `n_exact`, `n_lobpcg_ok`, `n_lobpcg_bad`. On the no-missing-values exit its
defaults stand in for the solve that never ran: `solver_chosen = NA`, `solver_probed = FALSE`,
`n_iter = 0`, `converged = NA`.
`slide_imp()` pins every later window to `attr(x, "solver_chosen")`, but only from a window
whose `solver_probed` is `TRUE`. A size-demoted window, or one whose probe never finished,
leaves the next window on `"auto"` (`R/slide_imp.R`: `pca_solver_current == "auto" &&`).
`pca_imp()` has no `subset`,
so `slide_imp()` does not call it at all for a window whose requested columns hold no missing
value (`R/slide_imp.R`: `targets_complete <- imp_method == "pca" &&`): such a window can
neither be probed nor be infeasible.
`group_imp()` does the same for a group whose feature columns are complete
(`R/group_imp.R`: `cmiss <- unname(mat_miss(`, `needs_imp <- vapply(`).

## Layer 2: the Rcpp boundary

`R/RcppExports.R` is the generated shim, calling
`.Call(`_slideimp_pca_imp_internal_cpp`, ...)` with 16 arguments.
`src/RcppExports.cpp` declares the C++ signature; the `RcppExport` SEXP wrapper
(`src/RcppExports.cpp`: `RcppExport SEXP _slideimp_pca_imp_internal_cpp(`) converts each
argument through `Rcpp::traits::input_parameter<>`. An
`Rcpp::RNGScope` makes `Rcpp::rnorm()` inside the kernel use R's RNG stream.
Registration is in `CallEntries` (`src/RcppExports.cpp`: `{"_slideimp_pca_imp_internal_cpp"`),
with `R_init_slideimp` and `R_useDynamicSymbols(dll, FALSE)`. The C++ entry symbol is
`pca_imp_internal_cpp` (`src/armaSVD.cpp`: `Rcpp::List pca_imp_internal_cpp(`), marked
`// [[Rcpp::export]]`.

`obj` arrives as `const arma::mat&`, so Rcpp constructs a read-only Armadillo view over the
SEXP rather than copying. `row_w` is taken **by value** because the
kernel normalizes it in place (`src/armaSVD.cpp`: `arma::rowvec row_w,`, `row_w /= row_w_sum;`).

## Layer 3: kernel setup (`src/armaSVD.cpp`)

- Argument sanity (`src/armaSVD.cpp`: `if (maxiter == 0)`): `maxiter >= 1`, `miniter <= maxiter`,
  `ncp >= 1`, `solver` in 0..2, non-empty and in-bounds `eligible_idx`, `ncp < min(nrX, n_elig)`.
  A further check (`src/armaSVD.cpp`: `if (n_elig != 0 && nrX >`) guards `nrX * n_elig`
  overflow. The block (`src/armaSVD.cpp`: `BLAS/LAPACK integer range guard`) guards the
  `blas_int` range: problem
  dimensions, the `26 * n` dsyevr workspace minimum, and the `~3k` Rayleigh-Ritz
  dsyevd workspace LOBPCG would need.
- Orientation: `tall = (nrX >= n_elig)`; `k = regularized ? ncp + 1 : ncp`
  eigenpairs are requested; `n_gram = tall ? n_elig : nrX`, so the eigenproblem is always
  over the smaller dimension.
- Pass 1 (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "pass1_count_missing")`) counts NAs per
  eligible column. The next block (`src/armaSVD.cpp`: `fixed_cols = arma::find(col_nmiss == 0)`)
  splits eligible columns
  into fixed (no NA) and changing, building `perm` so that after permutation the changing
  columns are exactly positions `[n_fixed, n_elig)`. Every later loop derives its local
  index as `n_fixed + ci`, which is why no column-index vector is threaded through.
- The block (`src/armaSVD.cpp`: `arma::uvec miss_rows_offsets(`) builds a CSR-style pair:
  `miss_rows_offsets` (length `n_mc + 1`) and
  `miss_rows_flat` (length `n_missing`). `for_each_missing_in_col()` is the only
  iterator over them.
- Pass 2 (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "pass2_scan")`) fills `Xhat` and
  accumulates per column the observed weight
  `denom`, the reference mean `mu0`, and centered moments `obs_sum_c` / `obs_sumsq_c`.
  Centering against `mu0` before squaring is what keeps the later variance updates away from
  `E[x^2] - mu^2` cancellation. A check (`src/armaSVD.cpp`: `if (!(denom[j] > 0.0))`) rejects
  any eligible column with
  zero observed weight. The kernel then computes `mean_p` and, when `scale`, the per-column
  scale `et` floored at `PCA_TOL`; it then centers (and scales) `Xhat` in place.
  `initialize_missing_start()` allocates `fittedX` once and
  writes starting values at missing positions only: zeros for `init == 0`, `Rcpp::rnorm`
  draws otherwise.

Workspaces allocated once, outside the loop: `sqrt_row_w` / `inv_sqrt_row_w`
(the weighted space is `diag(sqrt(w)) * Xhat`, and `U` comes back in it); `GramWorkspace`
holding fixed dsyrk/dgemm dimension arguments (`src/gram_ops.h`: `struct GramWorkspace`);
`EigSymWorkspace`, which runs the dsyevr workspace query once
(`src/eig_sym_sel.h`: `bool init(arma::uword dim_n`); `AA_NxN` at `n_gram x n_gram`
and `X_work` at `nrX x n_elig` in the tall case only; `GramCache`, precomputing
the fixed-column Gram contribution; and `HybridEigContext`,
carrying the solver mode, warmup count, LOBPCG options and persistent LOBPCG state across
all outer iterations.

## Layer 4: the EM outer loop (`src/armaSVD.cpp`)

One iteration:

1. `Rcpp::checkUserInterrupt()` every 5 iterations.
2. `impute_restandardize<Scale>()` -- the fused heart of
   the loop. It copies the previous fit into `Xhat` at missing positions, updates each
   changing column's mean (and variance when scaling) incrementally from `mu0` plus centered
   contributions, rescales the column in one pass, and simultaneously refreshes the
   sqrt-weighted copy of that column that the Gram step consumes.
   Only changing columns are touched.
3. `SVD_triplet()` -- see Layer 5.
4. Shrinkage (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "post_svd")`): `tail` is the clamped
   residual energy `trace - sum of top ncp squared singular values`; under `regularized`,
   `sigma2` is that tail scaled by `scale_factor / denom_sigma`, multiplied by `coeff_ridge`,
   and capped at the `(ncp+1)`-th squared singular value;
   `lambda_shrinked[i] = v - sigma2 / v`.
5. Reconstruction (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "reconstruct")`).
   `scale_rows_cols_inplace()` (`src/pca_linalg_utils.h`: `void scale_rows_cols_inplace(`)
   fuses the row de-weighting of `U` with the shrunk column scaling, mutating `U` in place;
   `reconstruct_cols()` then does one `dgemm_` writing only the changing block
   of `fittedX`. The fixed block is never read inside the loop, and in the wide
   case the fixed rows of `V` it would need are not even computed there (Layer 5).
6. Objective (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "objective")`): `objective_full` minus
   the weighted residual sum of squares over
   missing cells only (`missing_weighted_residual_ss()`), floored at 0.

**Convergence** (`src/armaSVD.cpp`: `const double criterion =`): the criterion is
`abs(1 - objective / old)`, reproducing `missMDA::imputePCA()`. It is tested only once
`nb_iter >= miniter`, and stopping requires `criterion < threshold || objective < threshold`.
Hitting `maxiter` first sets `converged = false` and nothing else. The kernel stays silent
because a warning raised from C++ can longjmp past the frame's destructors, so `pca_imp()` warns
once, through `cli::cli_warn()`, for the winning restart only
(`R/pca_imp.R`: `if (!isTRUE(best_diag$converged)) {`). The loop
is `for (arma::uword nb_iter = 1;; ++nb_iter)`, whose only exits are in the block
(`src/armaSVD.cpp`: `if (stop_now)`).
`trace_iter > 0` prints on
iteration 1, every `trace_iter` iterations, and on the last; print precision is
restored by the `PrecGuard` RAII struct.

**Postprocessing** (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "postprocessing")`): in the wide
case `gemm_tn_block(Xhat, 0, n_fixed, Vn, V)` first fills the fixed rows of `V` that
`SVD_triplet()` leaves unwritten, from the
last iteration's `Vn` and the fixed columns of `Xhat`, which nothing writes after
`pass2_scan`; the tall case returns `V` complete. `reconstruct_cols(U, V, fittedX, 0, n_fixed)`
then completes the fixed block of `fittedX` the loop skipped, purely so the SSE
covers the whole matrix. Both calls no-op when `n_fixed == 0`. The fixed rows of `V` reach R
only through `mse`; `imputed_values` never reads them. The
SSE loop (`src/armaSVD.cpp`: `double sse = 0.0;`) computes the full residual in original
units (times `et^2` per column)
then subtracts the missing-cell contribution, giving an observed-cells-only `mse`,
the `nb.init` selection statistic. `imputed_values` is an `n_missing x 3`
matrix of `(row_1based, orig_col_1based, value)`, where the value comes from the final
`Xhat` at the missing position, de-standardized -- not from `fittedX`. The
ordinary-build return list carries `imputed_values`, `mse`, `n_iter`,
`criterion_final`, `converged`, `solver_chosen`, `n_exact`, `n_lobpcg_ok`, `n_lobpcg_bad`.
`solver_chosen` comes from `HybridEigContext::solver_chosen_code()`
(`src/hybrid_topk_eig.h`: `int solver_chosen_code() const`): 0 forced exact, 1 forced lobpcg,
2 auto never decided, 3 auto chose exact, 4 auto chose lobpcg. Under `solver = "auto"` only
(`R/pca_imp.R`: `if (solver == "auto" && i == 1L) {`), `pca_imp()`
collapses 1 and 4 to "lobpcg" and everything else to "exact"; the string is set at
(`R/pca_imp.R`: `solver_chosen = if (resolved_solver_code == 1L)`).

## Layer 5: Gram formation and factor recovery

`SVD_triplet()` (`src/svd_triplet.h`: `inline void SVD_triplet(`) never forms an SVD directly.
It builds the
weighted Gram into `AA_NxN` via `form_weighted_gram()`, takes
`trace_val = arma::trace(AA_NxN)` for the tail energy, calls `hybrid_topk_eig()`
(failure being an `Rcpp::stop`), converts eigenvalues to singular values and
inverses with a `PCA_TOL` floor, and recovers the other factor with one `dgemm_`
(`src/svd_triplet.h`: `LOC_TIC(timer, "recover_factors")`): tall gives eigenvectors as `V`
and `U = X_work * Vn` (`gemm_nn`); wide
gives them as `U` and forms only the changing rows of `V`,
`V[n_fixed.., ] = Xhat[, n_fixed..]' * Vn`, as one strided-block `dgemm_` in
`gemm_tn_block()` (`src/pca_linalg_utils.h`: `inline void gemm_tn_block(`), reading `n_fixed`
from the
Gram cache the way `form_weighted_gram()` does. `V` is sized once (`set_size`
with unchanged dimensions does not reallocate) and its fixed rows stay unwritten
until postprocessing fills them from the same `Vn` (Layer 4). `split_eigvec_block()`
writes the plain and `d_inv`-scaled copies in a single pass.

`form_weighted_gram()` (`src/gram_ops.h`: `inline void form_weighted_gram(`) has four branches:

| case | how the Gram is built |
| --- | --- |
| tall, cache active | copy cached fixed block, `dgemm_` cross block, `dsyrk_` changing-changing block |
| tall, no cache | one `dsyrk_` over all of `X_work` |
| wide, cache active | copy cached fixed Gram, `dsyrk_` with `beta = 1` adding the changing columns |
| wide, no cache | `dsyrk_` on `Xhat`, then scale entries by `sw[i] * sw[j]` |

Only the upper triangle is written, matching dsyevr's `uplo = 'U'`
(`src/eig_sym_sel.h`: `PEIG_UPLO = 'U'`). The
LOBPCG path mirrors it down first (`src/hybrid_topk_eig.h`: `mirror_upper_to_lower(A)`).

## Layer 6: the eigensolver

**Routing** (`src/hybrid_topk_eig.h`). `hybrid_topk_eig()` is a timing wrapper:
while the auto probe is undecided (`auto_timing_active()`) it brackets the real
call with a `steady_clock` measurement and feeds it to `record_auto_probe_time()`.
`hybrid_topk_eig_impl()` picks the path:

- `Exact` -- always dsyevr.
- `Lobpcg` -- dsyevr while `outer_iter < warmup_iters` or while there
  is no seed, then LOBPCG.
- `Auto` (`src/hybrid_topk_eig.h`: `else // HybridSolverMode::Auto`) -- exact for at least
  `max(warmup_iters, auto_min_exact_iter)`
  iterations, then `HYB_AUTO_N_PROBE_ITER_DEFAULT = 5` LOBPCG
  probes, then decides. `maybe_decide_auto()` picks LOBPCG only
  when its mean time is below `(1 - 0.10) * exact_mean`
  (`src/hybrid_topk_eig.h`: `HYB_AUTO_MARGIN_DEFAULT = 0.10`); otherwise
  it picks exact and drops the LOBPCG state so exact mode carries no seed-maintenance
  overhead (`src/hybrid_topk_eig.h`: `void choose_auto_exact()`).

On LOBPCG success (`src/hybrid_topk_eig.h`: `if (lob_res.converged)`) the eigenvalues are
moved out, `eigvecs` is copied from the
committed state, and `refresh_momentum_P()` rebuilds the search direction from the two most
recent eigenblocks. On failure `n_lobpcg_bad` increments and control falls through
to `eig_sym_sel()`; after any exact solve, `seed_lobpcg_state()` re-seeds the warm
start when `keep_lobpcg_seed_state()` allows.

**Exact path** (`src/eig_sym_sel.h`). `dsyevr_` is declared by hand rather than
pulled from `<R_ext/Lapack.h>`; the comment above it records that the R header clashes
with Armadillo's own declarations. Hidden Fortran character-length arguments are passed
explicitly (the trailing `1, 1, 1`). `EigSymWorkspace::init()` runs the
`lwork = -1` query once and floors the result at the documented `26 * n` / `10 * n` minima.
`eig_sym_sel()` requests eigenpairs `il..iu` with `range = 'I'`,
returning exactly the top `k`. dsyevr returns ascending order, so
it reverses eigenvalues and swaps whole eigenvector columns in place, then
`canonicalize_signs()` (`src/pca_linalg_utils.h`: `void canonicalize_signs(`) applies the
sklearn `svd_flip` convention.
`A` is destroyed by the call, which is safe because the Gram is rebuilt every
outer iteration.

**LOBPCG path** (`src/lobpcg_warm.h`). This is a **warm-start-only** solver:
`lobpcg_solve()` hard-errors via `Rcpp::stop` on an unseeded state, so
callers must check `LOBPCGState::seeded()` first. The state carries `X` (`n x k`,
orthonormal, sign-canonicalized approximate top-k eigenvectors) and `P` (`n x k`, the search
direction). The header's "P invariant" note states that the committed `P` is
deliberately *not* orthonormal and not orthogonal to `X`; entry-time cleanup
(`src/lobpcg_warm.h`: `bool have_P = qr_project_out_qr(`)
re-projects and re-orthonormalizes it. Both matvecs are recomputed fresh at entry
(`src/lobpcg_warm.h`: `arma::mat AX = A * X;`) because `A` changes between outer calls. The
initial Rayleigh-Ritz is over the
joint basis `[X, P]` (`rr_two_block()`), falling back to `X` alone when `P` is
unusable (`src/lobpcg_warm.h`: `fall back to X-only initial RR`).

The iteration (`src/lobpcg_warm.h`: `for (; it < opt.maxiter; ++it)`) is standard three-block
LOBPCG over `[X, R_active, P_active]`.
Besides convergence it exits on `maxiter`, a divergence guard at relative
residual > 10, a stall counter over `stall_window = 5`, residual collapse
after QR, `P` collapse, and `LOBPCG_RR_FAILED` when
`rr_two_block` or `small_eig_desc` fails. `MATVEC_REFRESH = 20`
periodically recomputes `AX`/`AP` from scratch to shed accumulated
rounding. State is committed only on success (`src/lobpcg_warm.h`: `if (out.converged)`).
Small Rayleigh-Ritz problems go
through `small_eig_desc()`, which calls `dsyevd_` twice (query, then solve)
with buffers cached on the `LOBPCGState`. Warm-start construction is two entry
points over one helper, `build_momentum_P()`, which forms
`P = X_prev - X_curr (X_curr' X_prev)`: `seed_lobpcg_state()` after an exact
solve, with no norm floor; `refresh_momentum_P()` after a LOBPCG success, with
a floor of `1e-10 * sqrt(k)` below which the post-RR `P` is kept instead.

**Where BLAS/LAPACK is entered.**

| symbol | site | purpose |
| --- | --- | --- |
| `dsyrk_` | `syrk_upper()`, `form_weighted_gram()` (`src/gram_ops.h`) | Gram formation |
| `dgemm_` | `form_weighted_gram()` (`src/gram_ops.h`) | fixed-changing Gram cross block |
| `dgemm_` | `gemm_nn()`, `gemm_tn_block()` (`src/pca_linalg_utils.h`) | `gemm_nn` / `gemm_tn_block` factor recovery; the block form also fills the fixed rows of `V` after the loop (`src/armaSVD.cpp`: `gemm_tn_block(Xhat, 0, n_fixed, Vn, V)`) |
| `dgemm_` | `reconstruct_cols()` (`src/armaSVD.cpp`) | `reconstruct_cols`, strided column block |
| `dsyevr_` | `EigSymWorkspace::init()`, `eig_sym_sel()` (`src/eig_sym_sel.h`) | workspace query, exact top-k eigen |
| `dsyevd_` | `small_eig_desc()` (`src/lobpcg_warm.h`) | LOBPCG Rayleigh-Ritz |
| Armadillo | `qr_econ`, `solve`, `A * X` (`src/lobpcg_warm.h`) | LOBPCG orthonormalization and matvecs |

The exact path enters LAPACK once per outer iteration with an `n_gram x n_gram` problem. The
LOBPCG path never touches `n_gram`-sized LAPACK: its dense work is `n_gram x k` matvecs and
QRs plus a `3k x 3k` dsyevd.

## Data structures and ownership

- `obj` crosses as `const arma::mat&` and is never mutated in C++; the write-back happens in
  R (`R/pca_imp.R`: `obj[imp_indices] <- best_imputed[, 3]`). `row_w` crosses by value and
  is normalized in place (`src/armaSVD.cpp`: `arma::rowvec row_w,`, `row_w /= row_w_sum;`).
- `Xhat` is the working copy -- permuted column order, centered, scaled, mutated in
  place every outer iteration. `fittedX` is allocated once (`src/armaSVD.cpp`: `fittedX.zeros(`);
  dgemm overwrites only its changing block each iteration, and the fixed block is
  filled once, after the loop (Layer 5). `U` is mutated in place by `scale_rows_cols_inplace()`, so after
  the reconstruct step it no longer holds the orthonormal left factor.
- `AA_NxN` is destroyed by every eigen call: dsyevr overwrites it
  (`src/eig_sym_sel.h`: `destroys A: dsyevr overwrites`),
  and the LOBPCG path mirrors its upper triangle down first
  (`src/hybrid_topk_eig.h`: `mirror_upper_to_lower(A)`).
- `HybridEigContext` carries state across outer iterations: `lobpcg_state` (including its
  dsyevd `work`/`iwork` buffers), `X_prev`, `auto_choice`, the auto counters and timing sums,
  and `n_exact`/`n_lobpcg_ok`/`n_lobpcg_bad` (`src/hybrid_topk_eig.h`: `LOBPCGState lobpcg_state;`).
  When auto picks
  exact, only `lobpcg_state.X`, `lobpcg_state.P` and `X_prev` are reset.
- `V` is sized once by `SVD_triplet()` (`src/svd_triplet.h`: `V.set_size(Xhat.n_cols, ncp);`).
  In the wide case its rows
  `0 .. n_fixed-1` hold whatever the allocation held until postprocessing writes them
  (`src/armaSVD.cpp`: `gemm_tn_block(Xhat, 0, n_fixed, Vn, V)`); nothing reads them inside the
  loop. `Vn` is written only by
  `split_eigvec_block()` (`src/svd_triplet.h`: `Vn.set_size(n, ncp);`) and read only by the
  factor-recovery
  `dgemm_` calls, which is what lets the post-loop fill reuse the last iteration's copy.
- Return crosses as an `Rcpp::List`, copying the triplet matrix
  into R. Only missing cells are transported, never the full imputed matrix. The
  `mat_stats.cpp` entry points instead return by value: `mean_imp_col()` an `n x p` matrix
  and `col_vars()` a length-`p` row vector, each copied into a SEXP.
- `PCA_TOL = 1e-15` (`src/pca_linalg_utils.h`: `inline constexpr double PCA_TOL`) is the
  shared floor in `armaSVD.cpp` and
  `svd_triplet.h`: scale floors, inverse row-weight floor,
  shrinkage guard, singular-value inverse guard.
  LOBPCG has its own floors
  (`src/lobpcg_warm.h`: `double qr_rel_tol = 1e-12;`, `const double min_temporal_norm =`),
  and R keeps row weights at or above `1e-10`, rejecting a smaller supplied
  `row.w` (`R/pca_imp.R`: `lower = 1e-10`) and flooring the weights it computes
  for `row.w = "n_miss"` (`R/pca_imp.R`: `row.w[row.w < 1e-10] <- 1e-10`), and
  applies a `.Machine$double.eps` variance floor.
- `src/matrix_checks.h` is not in `armaSVD.cpp`'s include chain. It is included by
  `src/mat_stats.cpp` and `src/imputed_value.h`, and its only function `stop_on_inf()`
  is reached from `src/mat_stats.cpp` via `check_inf()`, which `pca_imp()` hits
  indirectly through
  `mean_imp_col()` (`R/mean_imp_col.R`: `check_inf(obj)`), and only when `post_imp` has
  columns left to fill.

## Compile-time switches

**`-DARMA_64BIT_WORD=1`** -- set in `src/Makevars` and `src/Makevars.win`. Makes
`arma::uword` 64 bits; without it a matrix past `2^32 - 1` elements wraps silently.
`arma_uword_bytes()` (`src/mat_stats.cpp`: `int arma_uword_bytes()`) exists only so the test
suite can assert the flag survived.

**`-DLOC_TIMER`** (`src/loc_timer.h`) -- when defined the header pulls in
`rcpptimer.h` and the `LOC_TIMER_*` macros become real: `LOC_TIMER_OBJ` declares an
`Rcpp::Timer`, `LOC_TIC`/`LOC_TOC` bracket a phase, and `LOC_TIMER_PARAM`/`LOC_TIMER_ARG`
thread the timer through signatures. When undefined, everything collapses to `((void)0)` and
the parameter macros expand to nothing -- which is why `SVD_triplet()`'s
signature ends with `LOC_TIMER_PARAM(timer)` and its call site ends
with `LOC_TIMER_ARG(pca_imp_gram)`. Timed phases (`src/armaSVD.cpp`):
`pca_imp_internal_cpp_total`, `pass1_count_missing`,
`pass2_scan`, `gram_cache_init`, `restandardize`, `svd`,
`post_svd`, `reconstruct`, `objective`, `postprocessing`,
and inside `SVD_triplet` (`src/svd_triplet.h`) the `form_gram`, `eig` and `recover_factors`
phases. `rcpptimer` is deliberately absent from `LinkingTo`, so
its include directory is added by hand in `load_all1()`
(`dev/dev-utils.R`: `inc <- system.file("include", package = "rcpptimer")`).

**`-DPCA_IMP_DIAGNOSTICS`** (`src/armaSVD.cpp`) -- defaulted to 0, set by
`load_all1(debug = TRUE)` (`dev/dev-utils.R`: `sprintf("-DPCA_IMP_DIAGNOSTICS=%d"`). It gates
three things:
`hyb_ctx.init_logs(maxiter)`, which allocates the five per-iteration log
vectors in `HybridEigContext` (`src/hybrid_topk_eig.h`: `void init_logs(size_t n)`) -- without
it they stay empty and
`log_at()` (`src/hybrid_topk_eig.h`: `auto log_at = [&](`) is a no-op because its bounds check
fails; the
history buffers `eigval_hist`, `obj_hist`, `subspace_cos_hist` and their update
block, which adds an `arma::svd` of `V_prev' * V_new` per iteration to track
subspace rotation; and a wider return list. The counters `n_exact`,
`n_lobpcg_ok`, `n_lobpcg_bad` are always maintained and returned; only the per-iteration
logs are gated.

## Parallelism

- **Inside the PCA kernel: none.** `src/armaSVD.cpp` does not include `RcppThread.h`, and
  neither does any header it pulls in. All intra-call parallelism comes from whatever BLAS
  the R installation links.
- **RcppThread** is reached through `par_for()` (`src/par_for.h`: `inline void par_for(`),
  which both
  `mat_stats.cpp` kernels call (`col_vars_internal`, `mean_imp_col_internal`),
  as does the K-NN kernel (`src/impute_knn_brute.cpp`, not in this chain). With
  `n_threads <= 1` it runs the body inline and never starts an RcppThread pool,
  and `pca_imp()` reaches both kernels with the default `cores = 1`
  (`R/col_vars.R`: `col_vars <- function(`; `R/mean_imp_col.R`: `mean_imp_col <- function(`).
  `src/Makevars` links `RcppThread::LdFlags()`; `src/Makevars.win`
  does not.
- **mirai** is used from R only, one level above `pca_imp()`: `R/group_imp.R`
  (`mirai::mirai_map` over groups) and `R/tune_imp.R` (over tuning replicates). Both
  share the input matrix through `bigmemory` shared memory rather than serializing it
  and ship the worker body as a `.crate()`
  closure (`R/crate.R`: `.crate <- function(`). `slide_imp()` runs its windows sequentially
  (`R/slide_imp.R`: `for (i in seq_along(start)) {`).
- **BLAS threading** is pinned inside the mirai workers, not in the kernel:
  `RhpcBLASctl::blas_set_num_threads(1)` and `RhpcBLASctl::omp_set_num_threads(1)` in
  `R/group_imp.R` and `R/tune_imp.R`, gated by `pin_blas`.
  `check_pin_blas()` (`R/utils.R`: `check_pin_blas <- function(`) errors when `pin_blas = TRUE`
  without
  `RhpcBLASctl` installed, and otherwise emits a tip when more than one BLAS thread is
  detected.
