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
  +- R/parts.R             check_si_mi_names()             output = "mipca" only
  +- R/pca_imp.R           mat_miss()                   -> src/mat_stats.cpp
  |    +- R/pca_imp.R      no missing value: pca_imp_results(), or one
  |                        part through parts_to_si_mi()
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
  |              +- gram_cache_init()                  src/gram_ops.h -> si_la::syrk
  |              +- EM outer loop
  |                   +- impute_restandardize<>()      src/armaSVD.cpp
  |                   +- SVD_triplet()                 src/svd_triplet.h
  |                   |    +- form_weighted_gram()     src/gram_ops.h -> si_la::syrk/gemm
  |                   |    +- hybrid_topk_eig()        src/hybrid_topk_eig.h
  |                   |    |    +- lobpcg_solve()      src/lobpcg_warm.h -> dsyevd_
  |                   |    |    +- eig_sym_sel()       src/eig_sym_sel.h -> si_la::syevr_range
  |                   |    +- gemm_nn/gemm_tn_block()  src/pca_linalg_utils.h -> si_la::gemm
  |                   +- reconstruct_cols()            src/armaSVD.cpp -> si_la::gemm
  |                   +- missing_weighted_residual_ss()  src/armaSVD.cpp
  |                   +- convergence test
  |              +- postprocessing, loc and value of the imputed cells
  |              +- Rcpp::List back to R
  |
  +- R/pca_imp.R           mse guard, per restart
  +- R/pca_imp.R           clamp
  +- R/pca_imp.R           PCA part from loc and value; over the ineligible
  |                        kept columns, a mean part over those with a hole
  |                        (a left NA part without post_imp) and a part with
  |                        nothing to fill over the rest
  +- R/pca_imp.R           pca_imp_info(): the call's record
  +- output = "matrix"
  |    +- R/parts.R        fill_parts(): one copy of obj's kept columns, the
  |    |                   parts cut to them and written
  |    +- R/pca_imp.R      pca_imp_results()
  |         +- R/utils.R   new_slideimp_results()
  +- output = "mipca"
       +- R/parts.R        parts_to_si_mi(): one slice of obj per cut part
            +- R/mi_pca_s3.R   .si_mi_chunk(), .si_mi_new(), .si_mi_join()
```

C++ files in the chain: 4 translation units (`src/armaSVD.cpp`, `src/RcppExports.cpp`,
`src/mat_stats.cpp`, `src/si_lapack.cpp`) and 9 headers (`svd_triplet.h`, `gram_ops.h`,
`eig_sym_sel.h`, `hybrid_topk_eig.h`, `lobpcg_warm.h`, `pca_linalg_utils.h`, `si_lapack.h`,
`matrix_checks.h`, `loc_timer.h`). `si_lapack.cpp/.h` is the wrapper layer ("What the
wrappers reproduce"), and every BLAS and LAPACK call of the chain goes through it. No file
of the chain uses Armadillo, and since Phase 3 of `dev/eigen_conversion.md` the package
does not link to it: RcppArmadillo is not in `LinkingTo`, so the generated
`src/RcppExports.cpp` does not include `RcppArmadillo.h`.
Never touched: `par_for.h`, `find_windows.cpp`, `find_windows_flank.cpp`,
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
`output` is resolved by `match.arg()`; under `"mipca"`, `check_si_mi_names()`
(`R/parts.R`: `check_si_mi_names <- function(obj, call = parent.frame()) {`) refuses an
input without unique row and column names before anything else runs, the
no-missing-values exit included.
`subset` becomes `keep`, the result's columns, ascending, or `NULL` for every
column (`R/pca_imp.R`: `keep <- resolve_narrow_subset(subset, obj)`), through
`resolve_narrow_subset()` (`R/utils.R`:
`resolve_narrow_subset <- function(subset, obj, call = parent.frame()) {`), which
refuses an unknown name and an empty selection and turns a selection of every
column into `NULL`. Every eligible column still enters the kernel: `subset`
narrows only what is returned.
`cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)`, and when no kept column has a hole
the kept columns return unchanged, with the same record as the main path (`R/pca_imp.R`:
`info <- pca_imp_info(post_imp = post_imp, solver = solver)`): through
`pca_imp_results()` for the matrix output, so the result has the same class and
attributes either way, or as one part over the kept columns, with an empty map, through
`parts_to_si_mi()`. The kernel is not called. `eligible`
combines the `colmax` missing-rate ceiling with a non-degenerate-variance test from
`col_vars(obj)`. Feasibility caps (`R/pca_imp.R`: `max_ncp_rows <- nrow(obj) - 2L`)
require `ncp` to be at most
`min(nrow(obj) - 2L, n_elig - 1L)`; both abort paths carry class `slideimp_infeasible`,
which `group_imp()` and `slide_imp()` catch, and a third such abort
(`R/pca_imp.R`: `if (!any(cmiss[out_cols][eligible[out_cols]] > 0L)) {`) fires when
every kept column with missingness is ineligible.

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
the lowest `mse` wins. Each restart's `mse` must be finite or the call stops as a
package defect (`R/pca_imp.R`: `if (!is.finite(cur_obj)) {`): the kernel's SSE sums
the residual over every cell of each eligible column, holes included, before
subtracting the holes' share (`src/armaSVD.cpp`: `sse -= s * ej2;`), so a NaN or
Inf at any imputed cell surfaces there, and the kernel's values themselves are
never NA past this point. `clamp` is applied to the kernel's `value` (imputed
values only). The result is then built from three parts
(`R/parts.R`), which together cover every column: the PCA part
(`R/pca_imp.R`: `pca_part <- new_part(`) is the kernel's `loc` and `value` as they
are, `loc` being an integer matrix with columns `row` and `col` of 1-based
original-matrix indices, over the eligible columns; the
ineligible columns with a hole are one part, filled through `fallback_part()`
under `post_imp = TRUE` and a `left_na` part under `post_imp = FALSE`; and those
without a hole are a mean part with an empty map. So `fallback_part()` copies
only the columns it fills, and a `left_na` part names only columns that keep a
hole. With `subset`, the PCA part still
covers every eligible column and the doors cut it to `keep`, while the mean parts
cover only the kept ineligible columns, so no mean is computed for a column the
result does not hold. An integer input is coerced to double
first. `mean_filled_cols` is the kept ineligible columns that had a hole, under
`post_imp` only (`R/pca_imp.R`: `mean_cols[cmiss[mean_cols] > 0L]`), numbered as
the result's columns, and `mean_filled_cells` is always empty. `pca_imp_info()`
(`R/pca_imp.R`: `pca_imp_info <- function(`) gathers the call's record, one list
for both outputs: `imp_method`, `fallback`, `post_imp`, `mean_filled_cols`,
`mean_filled_cells`, `solver_requested`, `solver_chosen`, `solver_probed`, `n_iter`,
`criterion_final`, `converged`, `n_exact`, `n_lobpcg_ok`, `n_lobpcg_bad`. On the
no-missing-values exit its defaults stand in for the solve that never ran:
`solver_chosen = NA`, `solver_probed = FALSE`, `n_iter = 0L`, `converged = NA`.
`n_iter` is an integer on both paths, as the kernel returns it.

For `output = "matrix"`, `fill_parts()` (`R/parts.R`: `fill_parts <- function(obj, parts, keep = NULL) {`)
copies `obj` once and writes the parts (`R/pca_imp.R`: `obj <- fill_parts(obj, parts, keep = keep)`),
checking that each filling part fills exactly the holes of its columns, and
`pca_imp_results()` (`R/pca_imp.R`: `pca_imp_results <- function(obj, info) {`) builds the
result. It calls `new_slideimp_results()` (`R/utils.R`: `new_slideimp_results <- function(`),
which strips its input to a bare matrix, then sets class and standard attributes, and
adds the rest of the record as attributes, in the order above.

For `output = "mipca"`, every part takes the whole record as its `info`, and
`parts_to_si_mi()` (`R/parts.R`: `parts_to_si_mi <- function(`)
builds a one-draw `si_mi`: one chunk per part with columns, each a slice of `obj`
(no slice when one part covers every column: `obj` itself, cut to its dims and
dimnames when it carries more), its map rebased to the slice, its values as the
draw, `tf_identity()` as the backward map and the empty diagnostic frame, joined
with the parts' columns as the placement. The record rides on every chunk, which
a `cbind()` keeps, though the join does not renumber its column indices, and the
PCA chunk's `prov` also holds `rank` and `iter`.

Both doors take `keep`, the input columns to return, and cut the parts to it
after validating them (`R/parts.R`: `cut_parts <- function(parts, lins, keep, n_col) {`).
With `keep`, the matrix door copies only the kept columns, and the `si_mi` door
slices each part to its kept columns and places them in the fit's own numbering.

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

The kernel takes what R holds (`dev/eigen_conversion.md`, "Rules that apply to every
step"): `obj` as a `const Rcpp::NumericMatrix &`, read in place through its pointer and
never copied; `eligible_idx` as an integer vector of 0-based columns; every count and
setting as `int`, which `pca_imp()` converts at the call where it holds a double
(`R/pca_imp.R`: `maxiter = as.integer(maxiter),`); and `row_w` as a
`const Rcpp::NumericVector &`, copied once into the kernel's own vector because the
kernel normalizes it (`src/armaSVD.cpp`: `std::vector<double> w(row_w.begin(), row_w.end());`,
`w[i] /= row_w_sum;`). Inside, a count is held in `std::size_t` once checked, so an offset
such as `j * nrX` is computed in `std::size_t`; a row or column index of `obj` is `int`
(`src/armaSVD.cpp`: `Types, for the whole file.`).

## Layer 3: kernel setup (`src/armaSVD.cpp`)

- Argument sanity (`src/armaSVD.cpp`: `if (maxiter < 1)`): `maxiter >= 1`,
  `0 <= miniter <= maxiter`, `ncp >= 1`, `init >= 0`, `trace_iter >= 0`, `solver` in 0..2,
  and under LOBPCG `lobpcg_maxiter >= 1` and `warmup_iters >= 0`; a non-empty
  `eligible_idx` whose every index is a column of `obj`, refused by
  `stop_on_bad_col_index()` in the words the K-NN kernel uses (`src/matrix_checks.h`:
  `static inline void stop_on_bad_col_index(`); `ncp < min(nrX, n_elig)`. `NA_integer_` is
  `INT_MIN`, so each lower bound refuses it; as `arma::uword` a negative count wrapped to
  a huge one. A further check (`src/armaSVD.cpp`: `prevent nrX * n_elig overflow later.`)
  guards `nrX * n_elig`. The block (`src/armaSVD.cpp`: `LAPACK workspace guard`) guards
  the `int` range of the `26 * n` dsyevr workspace minimum and of the `~3k` Rayleigh-Ritz
  dsyevd workspace LOBPCG would need; the dimensions themselves fit an `int`, being R's.
  More missing cells than an R matrix has rows is refused once they are counted
  (`src/armaSVD.cpp`: `The eligible columns hold`), since `loc` has a row for each.
- Orientation: `tall = (nrX >= n_elig)`; `k = regularized ? ncp + 1 : ncp`
  eigenpairs are requested; `n_gram = tall ? n_elig : nrX`, so the eigenproblem is always
  over the smaller dimension.
- Pass 1 (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "pass1_count_missing")`) counts NAs per
  eligible column. The next block (`src/armaSVD.cpp`: `perm lists the complete columns first`)
  splits eligible columns
  into fixed (no NA) and changing, building `perm` so that after permutation the changing
  columns are exactly positions `[n_fixed, n_elig)`. Every later loop derives its local
  index as `n_fixed + ci`, which is why no column-index vector is threaded through.
- The block (`src/armaSVD.cpp`: `std::vector<std::size_t> miss_rows_offsets(n_mc + 1, 0);`)
  builds a CSR-style pair: `miss_rows_offsets` (length `n_mc + 1`, `std::size_t`) and
  `miss_rows_flat` (length `n_missing`, `int` row indices). `for_each_missing_in_col()` is the only
  iterator over them.
- Pass 2 (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "pass2_scan")`) fills `Xhat` and
  accumulates per column the observed weight
  `denom`, the reference mean `mu0`, and centered moments `obs_sum_c` / `obs_sumsq_c`.
  Centering against `mu0` before squaring is what keeps the later variance updates away from
  `E[x^2] - mu^2` cancellation. A check (`src/armaSVD.cpp`: `if (!(denom[j] > 0.0))`) rejects
  any eligible column with
  zero observed weight. The kernel then computes `mean_p` and, when `scale`, the per-column
  scale `et` floored at `PCA_TOL`; it then centers (and scales) `Xhat` in place.
  `fittedX` is allocated once, zero-filled, and `initialize_missing_start()`
  writes starting values at missing positions only: zeros for `init == 0`, `Rcpp::rnorm`
  draws otherwise.

Workspaces allocated once, outside the loop: `sqrt_row_w` / `inv_sqrt_row_w`
(the weighted space is `diag(sqrt(w)) * Xhat`, and `U` comes back in it); `GramWorkspace`
holding fixed dsyrk/dgemm dimension arguments (`src/gram_ops.h`: `struct GramWorkspace`);
`EigSymRange`, which holds the range `il..iu` of the top `k` eigenpairs dsyevr is asked
for (`src/eig_sym_sel.h`: `struct EigSymRange`) and refuses a `k` outside `1..n_gram`, an
`Rcpp::stop` not expected to be reachable after the argument checks; `AA_NxN` at `n_gram x n_gram`
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
   `reconstruct_cols()` then does one dgemm, through `si_la::gemm`, writing only the changing block
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
is `for (std::size_t nb_iter = 1;; ++nb_iter)`, whose only exits are in the block
(`src/armaSVD.cpp`: `if (stop_now)`).
`trace_iter > 0` prints on
iteration 1, every `trace_iter` iterations, and on the last; print precision is
restored by the `PrecGuard` RAII struct.

**Postprocessing** (`src/armaSVD.cpp`: `LOC_TIC(pca_imp_gram, "postprocessing")`): in the wide
case `gemm_tn_block()` over columns `0 .. n_fixed-1` of `Xhat` first fills the fixed rows of `V` that
`SVD_triplet()` leaves unwritten, from the
last iteration's `Vn` and the fixed columns of `Xhat`, which nothing writes after
`pass2_scan`; the tall case returns `V` complete. `reconstruct_cols()` over columns
`0 .. n_fixed-1` then completes the fixed block of `fittedX` the loop skipped, purely so the SSE
covers the whole matrix. Both calls no-op when `n_fixed == 0`. The fixed rows of `V` reach R
only through `mse`; `value` never reads them. The
SSE loop (`src/armaSVD.cpp`: `double sse = 0.0;`) computes the full residual in original
units (times `et^2` per column)
then subtracts the missing-cell contribution, giving an observed-cells-only `mse`,
the `nb.init` selection statistic. The imputed cells come back as `loc`, an
`n_missing x 2` integer matrix with columns `row` and `col`, 1-based in `obj`, and
`value`, a double vector (`src/armaSVD.cpp`: `Rcpp::IntegerMatrix loc(`), in column-major
order; each value comes from the final `Xhat` at the missing position, de-standardized --
not from `fittedX`. The
ordinary-build return list carries `loc`, `value`, `mse`, `n_iter` (an integer),
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
the trace of `AA_NxN` for the tail energy (`gram_trace()`), calls `hybrid_topk_eig()`
(failure being an `Rcpp::stop`), converts eigenvalues to singular values and
inverses with a `PCA_TOL` floor, and recovers the other factor with one dgemm
through `si_la::gemm` (`src/svd_triplet.h`: `LOC_TIC(timer, "recover_factors")`): tall gives
eigenvectors as `V` and `U = X_work * Vn` (`gemm_nn`, dgemm whatever the shape, never dgemv);
wide gives them as `U` and forms only the changing rows of `V`,
`V[n_fixed.., ] = Xhat[, n_fixed..]' * Vn`, as one strided-block dgemm in
`gemm_tn_block()` (`src/pca_linalg_utils.h`: `inline int gemm_tn_block(`), reading `n_fixed`
from the
Gram cache the way `form_weighted_gram()` does. Both helpers take raw buffers and their
dimensions, refuse operands whose dimensions do not fit, and return the wrapper's status; a
non-zero one is an `Rcpp::stop` at the call site
(`src/svd_triplet.h`: `const int st = pca_detail::gemm_nn(`). Every dimension is checked
against the `int` range at kernel setup, so that error is not expected to be reachable.
`SVD_triplet()` takes raw buffers and reads every dimension from `GramWorkspace` and the
eigenpair range; `src/armaSVD.cpp` sizes each buffer once, before the EM loop, with the sizes
the header of `SVD_triplet()` lists (`src/armaSVD.cpp`: `every buffer SVD_triplet() reads or writes is sized here`).
So `V`'s fixed rows stay unwritten until postprocessing fills them from the same `Vn`
(Layer 4). `split_eigvec_block()` writes the plain and `d_inv`-scaled copies in a single
pass.

`form_weighted_gram()` (`src/gram_ops.h`: `inline int form_weighted_gram(`) has four branches,
each entering BLAS through `si_la::syrk` and `si_la::gemm` with the arguments it passed to
`arma::dsyrk_` and `arma::dgemm_`:

| case | how the Gram is built |
| --- | --- |
| tall, cache active | copy cached fixed block, dgemm cross block, dsyrk changing-changing block |
| tall, no cache | one dsyrk over all of `X_work` |
| wide, cache active | copy cached fixed Gram, dsyrk with `beta = 1` adding the changing columns |
| wide, no cache | dsyrk on `Xhat`, then scale entries by `sw[i] * sw[j]` |

`src/gram_ops.h` takes raw buffers and the dimensions `GramWorkspace` holds, as `std::size_t`,
and `GramCache` keeps the fixed Gram and the wide case's scaled changing columns in
`std::vector<double>`. `gram_cache_init()` and `form_weighted_gram()` return the wrapper's
status; a non-zero one is an `Rcpp::stop` at the call site (`src/armaSVD.cpp`:
`if (gram_init_st != si_la::ok)`; `src/svd_triplet.h`: `if (gram_st != si_la::ok)`), not
expected to be reachable for the reason given above.

Only the upper triangle is written, matching dsyevr's `uplo = 'U'`
(`src/eig_sym_sel.h`: `PEIG_UPLO = 'U'`). The
LOBPCG path mirrors it down first (`src/hybrid_topk_eig.h`: `mirror_upper_to_lower(A, a_rows)`).

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

On LOBPCG success (`src/hybrid_topk_eig.h`: `if (lob_res.converged)`) the eigenvalues and
the committed `X` are copied into the caller's buffers, after a check that the state still
has the range's `n x topk` shape, which it keeps from its seed (an internal error
otherwise, not expected to be reachable), and `refresh_momentum_P()` rebuilds the search
direction from the two most recent eigenblocks. On failure `n_lobpcg_bad` increments and control falls through
to `eig_sym_sel()`; after any exact solve, `seed_lobpcg_state()` re-seeds the warm
start when `keep_lobpcg_seed_state()` allows.

**Exact path** (`src/eig_sym_sel.h`). `eig_sym_sel()` takes raw buffers and calls
`si_la::syevr_range`, which enters dsyevr with `range = 'I'` for eigenpairs `il..iu`, the
top `k`, asks for the workspace on every call and floors it at the documented
`26 * n` / `10 * n` minima ("What the wrappers reproduce"). `eigvals` and `eigvecs` are
`src/armaSVD.cpp`'s, sized there before the EM loop (`src/armaSVD.cpp`:
`std::vector<double> eigvecs(n_gram * k);`). dsyevr returns
ascending order, so `eig_sym_sel()` reverses eigenvalues and swaps whole eigenvector columns
in place, then `canonicalize_signs()` (`src/pca_linalg_utils.h`: `void canonicalize_signs(`)
applies the sklearn `svd_flip` convention. A non-zero status, from the wrapper or from an
`A` that is not `n_gram` square, makes `hybrid_topk_eig()` return false, and
`SVD_triplet()` stops (`src/svd_triplet.h`: `top-k symmetric eigensolver failed`).
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
(`src/lobpcg_warm.h`: `Block AX = mul(a, n, n, X);`) because `A` changes between outer calls. The
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
through `small_eig_desc()`, which symmetrizes the projected Gram and calls
`si_la::syevd`; the wrapper asks dsyevd for its workspace and allocates it on
every call, so the `LOBPCGState` holds `X` and `P` and nothing else.
Warm-start construction is two entry
points over one helper, `build_momentum_P()`, which forms
`P = X_prev - X_curr (X_curr' X_prev)`: `seed_lobpcg_state()` after an exact
solve, with no norm floor; `refresh_momentum_P()` after a LOBPCG success, with
a floor of `1e-10 * sqrt(k)` below which the post-RR `P` is kept instead.

**Storage on the LOBPCG path.** The iteration holds its matrices in `lobpcg_detail::Block`
(`src/lobpcg_warm.h`: `struct Block`), a column-major buffer whose members carry
Armadillo's names, so the three calls into `si_la` read as they did on `arma::mat`.
Every product goes through `si_la::matmul` (`src/lobpcg_warm.h`: `inline void matmul_into(`),
written the way Armadillo evaluated each one ("Products of two `arma::mat`" below): a
product inside a sum or difference, or one whose target is among its operands, is formed
first and added cell by cell; `-=` of a product whose operands do not include its target is
one call with `beta = 1` (`src/lobpcg_warm.h`: `inline void sub_mul(`). The boundary with
`src/hybrid_topk_eig.h` holds no Armadillo either: `lobpcg_solve()` takes the Gram as a raw
buffer and its two dimensions and reads it in place (`src/lobpcg_warm.h`:
`inline LOBPCGResult lobpcg_solve(const double *a,`), `LOBPCGState` holds `X` and `P` as
`Block`s, `LOBPCGResult` holds the eigenvalues as a `std::vector<double>`, and
`seed_lobpcg_state()` and `refresh_momentum_P()` take the two eigenvector blocks as `Block`s.
`hybrid_topk_eig_impl()` copies the eigenvectors of each solve into one `Block`, which seeds
or refreshes the state and then becomes `X_prev`.

**Where a factorization can fail on the LOBPCG path.** Three calls return a
`si_la` status, and `la_failed()` (`src/lobpcg_warm.h`: `inline bool la_failed(`)
sorts it for all three. A positive status, which is LAPACK's `info`, and
`err_nonfinite` are numerical failures and take the path below. `err_alloc`
stops with "could not allocate the workspace", and `err_arg` and `err_size`,
which no call here can produce once kernel setup has refused dimensions past
`int`, stop with "internal error". Per site:

| site | call | a numerical failure leads to |
| --- | --- | --- |
| `qr_orthonormalize()` | `si_la::qr_econ` | `V` reset and a rank of 0, the same as a block whose first pivot is 0. The caller decides: at entry, `have_P` is false and the first Rayleigh-Ritz is over `X` alone; on the active residuals, `LOBPCG_RESIDUAL_COLLAPSED`; on `P_active`, `LOBPCG_P_COLLAPSED`; in `build_momentum_P()`, no `P`, so the state is left unseeded after an exact solve, or the post-RR `P` is kept after a LOBPCG success |
| `qr_orthonormalize_update_matvec()` | `si_la::solve_tri` | an R error (`src/lobpcg_warm.h`: `Rcpp::stop("lobpcg_solve: triangular solve failed");`), which ends `pca_imp()`. A triangle that is singular to working precision is not a failure: the wrapper returns the least-squares solution with `approx` set, and the iteration goes on with it |
| `small_eig_desc()` | `si_la::syevd` | false. From `rr_two_block()` or the `X`-only form before the loop, `LOBPCG_INITIAL_RR_FAILED`; inside the loop, `LOBPCG_RR_FAILED` |

Every `LOBPCG_*` status but `LOBPCG_CONVERGED` leaves the state untouched and
sends `hybrid_topk_eig_impl()` to the exact solver on the same Gram, as above.
"What the wrappers reproduce" has what each call did before it moved to its
wrapper, and what was measured.

**Where BLAS/LAPACK is entered.**

| symbol | site | purpose |
| --- | --- | --- |
| dsyrk, through `si_la::syrk` | `syrk_upper()`, `form_weighted_gram()` (`src/gram_ops.h`) | Gram formation |
| dgemm, through `si_la::gemm` | `form_weighted_gram()` (`src/gram_ops.h`) | fixed-changing Gram cross block |
| dgemm, through `si_la::gemm` | `gemm_nn()`, `gemm_tn_block()` (`src/pca_linalg_utils.h`) | `gemm_nn` / `gemm_tn_block` factor recovery; the block form also fills the fixed rows of `V` after the loop (`src/armaSVD.cpp`: `Xhat.data(), nrX, n_elig, 0, n_fixed, Vn.data(),`) |
| dgemm, through `si_la::gemm` | `reconstruct_cols()` (`src/armaSVD.cpp`) | `reconstruct_cols`, strided column block |
| dsyevr, through `si_la::syevr_range` | `eig_sym_sel()` (`src/eig_sym_sel.h`) | exact top-k eigen, workspace query included |
| dsyevd, through `si_la::syevd` | `small_eig_desc()` (`src/lobpcg_warm.h`) | LOBPCG Rayleigh-Ritz, workspace query included |
| dgeqrf and dorgqr, through `si_la::qr_econ` | `qr_orthonormalize()` (`src/lobpcg_warm.h`) | LOBPCG orthonormalization |
| dtrtrs and dtrcon, and dgelsd as the fallback, through `si_la::solve_tri` | `qr_orthonormalize_update_matvec()` (`src/lobpcg_warm.h`) | `A * Q` carried through a QR without a new matvec |
| dgemm or dgemv, through `si_la::matmul` | `matmul_into()` (`src/lobpcg_warm.h`) | LOBPCG matvecs and projections |

The exact path enters LAPACK once per outer iteration with an `n_gram x n_gram` problem. The
LOBPCG path never touches `n_gram`-sized LAPACK: its dense work is `n_gram x k` matvecs and
QRs plus a `3k x 3k` dsyevd.

## What the wrappers reproduce

`src/si_lapack.h` declares, and `src/si_lapack.cpp` defines, R's BLAS and LAPACK as plain
functions over `double*`, dimensions and an int status. So far `src/pca_linalg_utils.h` calls
`si_la::gemm`, `src/gram_ops.h` calls `si_la::syrk` and `si_la::gemm`, `src/eig_sym_sel.h`
calls `si_la::syevr_range`, and `src/lobpcg_warm.h` calls `si_la::qr_econ`,
`si_la::solve_tri`, `si_la::syevd` and `si_la::matmul`; the table above says how each site enters BLAS and
LAPACK. Each wrapper stands in for one
call the chain makes, or made until its file was converted, and this section records what
that call does, so that a call
site can be moved to its wrapper without reading Armadillo again.

READ means read in `armadillo_bits` of the installed RcppArmadillo, version 15.6.2-1
(Armadillo 15.6.2). NOT READ is said where it applies. `dev/oracle_run.R` runs every
wrapper beside the call it stands in for.

**How RcppArmadillo configures Armadillo** (READ: `RcppArmadilloConfig.h`, `config.hpp`,
`typedef_elem.hpp`, `compiler_setup.hpp`). `ARMA_USE_BLAS` and `ARMA_USE_LAPACK` are
defined; `blas_int` is `int`; symbols take a trailing underscore; `ARMA_USE_FORTRAN_HIDDEN_ARGS`
is defined, so every Fortran character argument is followed by a hidden `size_t` length,
passed as a literal 1. `ARMA_WARN_LEVEL` is 2: a warning of level 1 or 2 is written to
`ARMA_CERR_STREAM`, which RcppArmadillo points at `Rcpp::Rcerr`, and is not an R warning;
one of level 3 is dropped. `ARMA_CHECK_NONFINITE` is defined.

| wrapper | stands in for | routine and options | on failure |
| --- | --- | --- | --- |
| `si_la::gemm` | `arma::dgemm_`, called directly | dgemm, every argument passed through | BLAS has no status |
| `si_la::gemv` | no direct call; `*` reaches it (below) | dgemv, unit strides | BLAS has no status |
| `si_la::syrk` | `arma::dsyrk_`, called directly | dsyrk, every argument passed through | BLAS has no status |
| `si_la::matmul` | `*`, `+=` and `-=` of a product of two `arma::mat` | dgemv or dgemm (below) | Armadillo throws on dimensions that do not conform |
| `si_la::syevr_range` | the `dsyevr_` that `src/eig_sym_sel.h` declared itself until it moved to this wrapper | dsyevr, `jobz = 'V'`, `range = 'I'`, `abstol = 0` | `info != 0`, or another number of eigenpairs than asked; no test for NaN or Inf (the wrapper has one, below) |
| `si_la::syevd` | `arma::dsyevd_`, which `small_eig_desc()` called directly until it moved to this wrapper | dsyevd, `jobz = 'V'` | `info != 0`; no test for NaN or Inf (the wrapper has one, below) |
| `si_la::qr_econ` | `arma::qr_econ()` | dgeqrf, then dorgqr | returns false, `Q` and `R` reset |
| `si_la::solve_tri` | `arma::solve()` on `trimatl()` or `trimatu()` | dtrtrs, dtrcon, and dgelsd as the fallback | throws `std::runtime_error` |
| `si_la::svd_values` | `arma::svd()` returning the values alone | dgesdd, `jobz = 'N'` | throws `std::runtime_error` |

The three direct calls are Armadillo's own declarations of the Fortran symbols (READ:
`def_blas.hpp`, `def_lapack.hpp`), with nothing between the call site and BLAS. Armadillo
declares no dsyevr and calls it nowhere, so `si_la::syevr_range` has no Armadillo
counterpart. The oracle compares it with dsyevr called as `src/eig_sym_sel.h` called it
before it moved to the wrapper, to the bit, and with `arma::eig_sym()`, which is dsyevd, to a tolerance (READ:
`auxlib::eig_sym_dc()` returns false, before LAPACK, when the upper triangle holds NaN or
Inf; it asks dsyevd for its workspace only from order 32 up, uses the larger of the answer
and `1 + 6 n + 2 n^2` and `3 + 5 n` computed in `blas_int`, and returns `info == 0`).

**The two eigen wrappers against the chain's calls.** Each enters the same routine with
the same options as the call it stands in for, and differs from it in three places. The
first two are on purpose; the third was found when the dsyevd call moved:

- *The eigenvalue buffer of dsyevr.* LAPACK documents `W` as `n` cells whatever the range.
  For a proper subset dsyevr goes through dstebz, which stores every eigenvalue its search
  interval holds and then drops those outside `il..iu`; when the range cuts a group of
  equal eigenvalues that is more than `iu - il + 1` of them. The chain's own call gave `W`
  exactly `topk` cells, so such a call wrote past the end of `eigvals`, and 1.2.0 does.
  `pca_imp()` reaches it when the completed data is rank
  deficient and `ncp`, plus one when regularized, exceeds its rank: samples that are exact
  copies of each other, holes included, or columns that are. `si_la::syevr_range` holds the
  `n` cells itself and copies `iu - il + 1` out, so its caller's buffer needs no more, and
  `eig_sym_sel()` now calls it. The write was not harmless. Measured on Windows
  (`dev/review_p2.md`, finding 1): on 30 samples of which 6 are distinct, 1.2.0 built from
  source and every commit before the move end the R process with a segmentation fault,
  `pca_imp(x, ncp = 5)` at its defaults included, and the builds after it return. The
  regression test is that input
  (`tests/testthat/test-pca_imp.R`: `test_that("pca_imp returns on replicate samples with ncp above their rank", {`).
- *NaN and Inf.* Neither of the chain's own calls tested its matrix, and LAPACK does not
  refuse such input. Measured on one machine with R's own LAPACK: dsyevr asked for every eigenpair
  of a matrix holding NaN does not return, a proper subset comes back with `info = 0` and
  another number of eigenpairs than asked, and dsyevd can return `info = 0` with NaN
  eigenvalues. Both wrappers test the triangle LAPACK reads, and only that one, before
  the routine is entered, and return `err_nonfinite`. `pca_imp()` refuses NaN and Inf in
  its input (`R/pca_imp.R`: `check_finite(obj)`), so only arithmetic inside the fit could
  put one in a Gram. Finite input whose squares overflow does: unscaled values near 1e160,
  or near 1e200 with `scale = TRUE`, end in the stop below under either solver, the same
  before and after the wrappers, and near 1e150 the fit returns (measured,
  `dev/review_p2.md`, finding 3). A test pins that stop
  (`tests/testthat/test-pca_imp.R`: `test_that("finite input whose squares overflow stops in the eigensolver", {`).
  Which status the wrapper returned there is not measured. On the exact path such a Gram
  is now the stop of a failed eigensolve (Layer 6, "Exact path"), where the chain's own
  call either reported a wrong count, the same stop, or did not return. On the LOBPCG
  path a small Gram holding NaN or Inf is a failed Rayleigh-Ritz, so that outer iteration
  goes to the exact solver. The chain's own dsyevd call did one of three things with it,
  measured by writing a NaN into the matrix (`dev/eigen_conversion.md`, Phase 2, step 4):
  at orders 8 to 12 dsyevd returned `info > 0`, the same failed Rayleigh-Ritz; at order 34
  it returned `info = 0` and the NaN went on, to end either in the error of a failed
  triangular solve or, when the Gram itself held the NaN, in `lobpcg_solve()` reporting
  convergence, because the largest residual norm is taken with a `>` that NaN never
  passes, as it still is (`src/lobpcg_warm.h`: `const double max_res = max_skip_nan(rn.data(), rn.size());`).
- *The workspace of dsyevd.* Until it moved to the wrapper, `small_eig_desc()` kept its
  dsyevd workspace on the `LOBPCGState` and only ever lengthened it, so a Rayleigh-Ritz
  problem that followed a larger one was solved with the larger one's workspace.
  `si_la::syevd` gives dsyevd the length of its own query on every call. The length is
  not neutral. READ, in LAPACK 3.12.0 (`dsyevd.f`, `dormql.f`, fetched from the reference
  repository, not from R's sources): dsyevd keeps `1 + 2 n + n^2` cells and hands the
  rest to dormtr, which is dormql for `uplo = 'U'`; dormql applies the reflectors in
  blocks only when its block size is below `n - 1` and that remainder holds its `T`
  array of 4160 cells and at least two columns a block, and otherwise one at a time. At
  the queried length, `1 + 6 n + 2 n^2`, the remainder is `n^2 + 4 n + 1`. With a block
  size of 32 (NOT READ: ilaenv; it is the value the measured boundaries give) that makes
  the eigenvectors depend on the length for orders 34 to 79 and for no other. Measured,
  `dev/probe_p2s4_syevd_lwork.R`, orders 2 to 140: with the workspace of a problem 1.5,
  2 or 3 times the order, every call of order 34 to 79 differs from the wrapper's and no
  other does; with that of the order plus one, orders 62 to 79; the eigenvalues never
  differ, and the eigenvectors by at most 6.9e-16. The oracle runs the pair with a
  longer workspace only at orders 30, 90 and 200, on either side of that range, and so
  reports them the same. In `pca_imp()` the first Rayleigh-Ritz of every `lobpcg_solve()`
  call after a state's first is of order `2 k` with the workspace of order `3 k` held, so
  the two differ from 17 eigenpairs up, which is `ncp` 16 when regularized
  (`14 k^2 + 10 k >= 4159`). LOBPCG iterates on the difference: step 4 has what it does
  to the result.

`dev/eigen_conversion.md`, "The review of the wrapper layer", has the measurements behind
the first two, and Phase 2, step 4, those behind the third.

**`arma::qr_econ(Q, R, X)`**, which `qr_orthonormalize()` called until it moved to the
wrapper (`src/lobpcg_warm.h`: `if (la_failed(si_la::qr_econ(v_rows, v_cols, V.memptr(), Q.memptr(),`). READ:
`fn_qr.hpp`, `auxlib::qr_econ()`, `auxlib::qr()`. For `m > n` it copies `X` into `Q`, runs
dgeqrf there, copies the upper triangle into an `n x n` `R`, zero below, and runs dorgqr on
`Q` with `k = n`. For `m <= n`, square included, it hands over to the full `qr()`: dgeqrf
in an `m x n` `R`, the first `m` columns copied into an `m x m` `Q`, `R` zeroed below its
diagonal, dorgqr with `k = m`. Both are dgeqrf on the `m x n` matrix and dorgqr on its first
`min(m, n)` columns, which is how the wrapper writes it. One workspace serves both
routines: dgeqrf is asked for its optimum, the answer is raised to `max(1, m, n)`, and
dorgqr is not asked. The wrapper has left Armadillo on two of these since Phase 4 of
`dev/eigen_conversion.md`: it asks each routine for its own workspace, with `max(1, n)`
as the floor of both, and it keeps the `m >= n` form alone, factoring in `Q`, and returns
`err_arg` for `m < n`, a shape no call site passes. An empty `X` succeeds without LAPACK. Nothing checks for NaN or Inf.
On `info != 0` from either routine `auxlib` returns false; `qr_econ()` then resets `Q` and
`R`, warns at level 3, which is dropped, and returns false. The caller tested that flag,
and now tests the wrapper's status the same way (Layer 6, "Where a factorization can fail
on the LOBPCG path"). A dimension past `blas_int` was an exception from Armadillo, so an
R error and not the flag (READ: `arma_assert_blas_size()`, called in `auxlib::qr_econ()`
before LAPACK, throws `std::runtime_error`), and so was a failed allocation (READ:
`arma_stop_bad_alloc()` throws `std::bad_alloc`; NOT READ: that every allocation of the
call reaches it). The wrapper returns them as `err_size` and `err_alloc`, and the call
site stops on both.

**`arma::solve(arma::trimatl(A), B)`**, which `qr_orthonormalize_update_matvec()` called as
`AQ = arma::solve(arma::trimatl(R.t()), AQ.t()).t()` until it moved to the wrapper
(`src/lobpcg_warm.h`: `if (la_failed(si_la::solve_tri('L', Rt.n_rows, Y.n_cols, Rt.memptr(),`). READ:
`fn_solve.hpp`, `glue_solve_tri_default::apply()`, `auxlib::solve_trimat_rcond()`,
`auxlib::rcond_trimat()`, `auxlib::solve_approx_svd()`. It is not a plain triangular
solve:

1. dtrtrs with `trans = 'N'`, `diag = 'N'` on a copy of `B`. `info != 0`, a zero on the
   diagonal, counts as failed, with the estimate left at 0.
2. Otherwise dtrcon estimates the reciprocal condition number in the 1-norm; `info != 0`
   there makes it 0.
3. If the solve failed, or the estimate is below `DBL_EPSILON` or NaN, Armadillo writes a
   level-2 warning ("system is singular; attempting approx solution", or "system is close
   to singular; rcond: ...") and discards the triangular answer. It copies the triangle
   alone into a square matrix and solves by dgelsd, from `B` as given, treating singular
   values under `n * DBL_EPSILON` times the largest as zero: the minimum-norm least-squares
   solution. That fallback refuses NaN or Inf in the triangle or in `B`.
4. If the fallback fails too, `solve()` throws `std::runtime_error`, "solve(): solution not
   found", which Rcpp turns into an R error.

So this call site went on with an approximate `AQ` whenever `R` was singular to working
precision, having printed a line and raised nothing. `si_la::solve_tri` does the same
arithmetic and returns the estimate and a flag in place of the print. For an empty system
Armadillo's estimate stays 0, so it prints the warning and returns the empty result; the
wrapper returns the empty result with the flag unset.

The call site now writes the two transposes out, as the oracle's `pair_solve_chain` does
(READ: `glue_solve_tri_default::apply()` evaluates the triangle's expression into a matrix
with `quasi_unwrap` before dtrtrs, so `R.t()` was formed there too), ignores the flag, and
so goes on with the approximate `AQ` as before, without the printed line. A non-zero
status, which is the fallback failing or meeting NaN or Inf, is an R error, as Armadillo's
throw in item 4 was, in the chain's own words and as an `Rcpp::exception` where
Armadillo's was a `std::runtime_error`.

**`arma::svd(X)`**, which the `PCA_IMP_DIAGNOSTICS` build called as
`arma::svd(V_prev.t() * V_new)` until step 6 of Phase 2, once per iteration from the
second on. READ: `fn_svd.hpp`, `auxlib::svd_dc()`. dgesdd with
`jobz = 'N'` on a copy of `X`. An empty `X` gives an empty result. NaN or Inf in `X` is
refused before LAPACK. The workspace is the larger of
`3 min(m, n) + max(max(m, n), 7 min(m, n))` and dgesdd's own optimum, which is asked for
only when `X` has 1024 cells or more. A refusal or `info != 0` makes `svd()` throw
`std::runtime_error`, "svd(): decomposition failed". The wrapper asks dgesdd for its
optimum at every size since Phase 4 of `dev/eigen_conversion.md`, and is compiled only
when `PCA_IMP_DIAGNOSTICS` is set, so the package as shipped does not reference dgesdd.

The call site now forms the product with `si_la::matmul`, which sends it where `*` did
(dgemv for one component, dgemm otherwise), and hands it to `si_la::svd_values`
(`src/armaSVD.cpp`: `st = si_la::svd_values(n_cp, n_cp, V_cross.data(), V_cross_sv.data());`).
A positive status or `err_nonfinite` is an R error, as the throw was, in the chain's own
words and as an `Rcpp::exception` (`src/armaSVD.cpp`:
`the SVD of the subspace product failed`); `err_alloc` and the wrapper's refusals of its
arguments stop too, the second as an internal error. The smallest singular value is
taken with `std::min_element`, where `s.min()` was (NOT READ: `Mat::min()`); singular
values are not NaN on success. No route to the failure is known: a Gram holding NaN is
refused before the eigenvectors exist. `dev/probe_p2s6_run.sh` forces it in both
versions.

**Products of two `arma::mat`** (READ: `glue_times_meat.hpp`, `mul_gemm.hpp`,
`mul_gemv.hpp`, `mul_syrk.hpp`). `A * B` and `A.t() * B` choose their routine at run time,
in this order, and `si_la::matmul` keeps the order:

| product | routine |
| --- | --- |
| either operand has no cells | no BLAS: the result is zero-filled |
| the result has one row | dgemv, `trans = 'T'`, on `B`, with the row as the vector |
| else `B` has one column | dgemv on `A`, `'N'` for `A * B`, `'T'` for `A.t() * B` |
| else `A.t() * A` on one object | dsyrk into the upper triangle, then mirrored; a vector, or `A` of 48 cells or fewer, by Armadillo's own loops |
| else | dgemm, `'N'` or `'T'` for `A`, `'N'` for `B` |

`C += A * B` and `C -= A * B` do not form the product and then add it: when `C` is not one
of the operands they are one call of the same routine with `beta = 1` and `alpha = 1` or
`-1` (`glue_times::apply_inplace_plus()`). When `C` is an operand, nested products
included (READ: `Glue::is_alias()` and `Op::is_alias()` ask their operands in turn), the
product is formed first and subtracted cell by cell: `V -= X * (X.t() * V)` was that form,
and is written so now (`src/lobpcg_warm.h`: `sub_cells(V, mul(X, tmul(X, V)));`). A product
inside a sum or difference is formed first too, as `P = X_prev - X_curr * (X_curr.t() * X_prev)`
was and is written now (`src/lobpcg_warm.h`: `sub_cells(P, mul(X_curr, tmul(X_curr, X_prev)));`):
READ, no overload in `operator_plus.hpp` or `operator_minus.hpp` names `glue_times`, so the
product is an operand of the sum (NOT READ: the proxy that then evaluates it). The oracle
does not cover those forms; the capture runs them. Three operands written `A * (B * C)` are two products, the
inner one first: `depth_lhs` counts only operands nested on the left (READ), so the
parenthesised product is one operand (NOT READ: `partial_unwrap`, which evaluates it).

Two things Armadillo does that `si_la::matmul` does not:

- For a square matrix operand of 4 x 4 or less (dgemv), or two square operands of the same
  size up to 4 x 4 (dgemm), Armadillo runs loops of its own and calls no BLAS (READ: the
  tests in `gemm::apply_blas_type()` and `gemv::apply_blas_type()`; NOT READ: the loops).
  The wrapper calls BLAS, so such a product can differ in its last bits. In the chain that
  needs a Gram of at most 4 x 4 under LOBPCG. Measured (`dev/review_p2.md`, finding 2): with
  `solver = "lobpcg"` on a Gram of order 2 to 4, when one eigenpair or every eigenpair is
  asked for, 267 of 280 fits differ from the Armadillo build, by at most 5.2e-9 absolute;
  fits outside that rule, and every fit under the exact solver, are identical.
- The dsyrk route for `A.t() * A` on one object. No call site of the chain has that form.

An in-place product with an empty inner dimension is passed on to dgemm with a leading
dimension of 0 (READ from the arguments `gemm::apply_blas_type()` builds; not run, since
BLAS answers that through xerbla). The wrapper scales `C` by `beta` and calls nothing.

**Calls no wrapper stands in for.** Four norms, a dot product and a trace, and the plain
loops that replace them in `src/lobpcg_warm.h` and `src/svd_triplet.h`:

- `arma::norm(P, "fro")` in `build_momentum_P()` was dnrm2 over all cells, recomputed by a
  scaled loop when the result is 0 or not finite (READ: `fn_norm.hpp`,
  `op_norm::vec_norm_2_direct_std()`; NOT READ: the scaled loop). Now a plain sum of
  squares (`src/lobpcg_warm.h`: `!(frobenius_norm(P) > min_projected_fro_norm))`), as ruled
  (`dev/eigen_conversion.md`, "Not wrapped"): not bitwise equal to dnrm2, and only a
  threshold reads it.
- `arma::norm(A, "inf")` was `as_scalar(max(sum(abs(A), 1), 0))` with no BLAS (READ:
  `op_norm::mat_norm_inf()`). `sum(, 1)` starts each row's sum from its first column and adds
  the others in column order, whichever of its two paths it takes (READ:
  `op_sum::apply_mat_noalias()`, `op_sum::apply_proxy_noalias()`), and `max()` starts from
  `-Inf` and keeps a value only when it is `>` the best so far, so NaN is skipped (READ:
  `op_max::direct_max()`, which keeps two running maxima and returns the larger; the maximum
  of values that are not NaN does not depend on that order). `norm_inf()` does the same, from
  a zero sum, which adds nothing to a value that is not negative (`src/lobpcg_warm.h`:
  `row_sum[i] += std::abs(col[i]);`), so it is bitwise equal.
- The residual norms, `arma::sqrt(arma::sum(arma::square(R), 0))`, which no note named before
  the second half of step 4: `sum()` of `square()` of a matrix takes each column's sum with
  `op_dot::direct_dot()` (READ: `op_sum::apply()` for `eOp<T1, eop_square>`,
  `op_sum::apply_mat_square_noalias()`), so two partial sums up to 32 rows and ddot above.
  Now a plain sum of squares (`src/lobpcg_warm.h`: `inline std::vector<double> column_norms(`),
  not bitwise equal to either, on the same footing as the Frobenius norm: the norms feed the
  tests of the iteration and `max_rel_res`, which only the diagnostics build reports
  (`src/armaSVD.cpp`: `Rcpp::Named("eig_lobpcg_max_rel_res")`). `rn.max()` reached
  `op_max::direct_max()` (READ: `Mat::max()`), which `max_skip_nan()` reproduces.
- `arma::dot(vs_top.head(ncp), vs_top.head(ncp))`, the sum of the squared leading
  singular values that the tail energy subtracts from the trace, was Armadillo's own
  loop up to 32 values, two partial sums over the even and the odd positions added at
  the end, and ddot above (READ: `op_dot::apply()`, where two `subview_col` take the
  vector branch, `partial_unwrap` being fast for them; `op_dot::direct_dot()`;
  `op_dot::direct_dot_generic()`). Now `sum_squares_two_partial()` (`src/armaSVD.cpp`:
  `static inline double sum_squares_two_partial(`), which takes the two partial sums at
  every length: bitwise equal up to 32 components, not above, where it stands in for
  ddot. Kept in that order on the footing of `gram_trace()` below
  (`dev/eigen_conversion.md`, step 6 of Phase 2, has the measurements).
- `arma::accu(row_w)`, the sum the row weights are divided by, was two partial sums over
  the even and the odd positions, added at the end (READ: `op_accu_mat::apply()`,
  `arrayops::accumulate()`). Now `sum_two_partial()` (`src/armaSVD.cpp`:
  `static inline double sum_two_partial(`), in that order, so it is bitwise equal. Equal
  weights sum exactly in either order; with unequal ones a plain sum moves results
  (`dev/eigen_conversion.md`, step 6), so it is kept on the same footing.
- The rest of what `src/armaSVD.cpp` asked of Armadillo did no arithmetic that depends
  on an order, and is written as loops: `row_w /= row_w_sum` divided each cell (READ:
  `arrayops::inplace_div_base()`); `cumsum()` of a column kept one running sum (READ:
  `op_cumsum::apply_noalias()`); `is_finite()` tested each cell (READ:
  `arrayops::is_finite()`); the element-wise `/`, `%`, `-`, `+` and `sqrt()` behind the
  column means and scales, `find()`, `elem()`, `any()` and `max()` of an index vector
  (NOT READ; the captures are bitwise across the change).

- `arma::trace(AA_NxN)` in `SVD_triplet()`, the tail energy's trace, was two partial sums,
  one over the even and one over the odd diagonal cells, added at the end (READ:
  `fn_trace.hpp`, `trace()` of a `Base`), no BLAS. Now `gram_trace()` (`src/svd_triplet.h`:
  `static inline double gram_trace(`), which sums in that order, so it is bitwise equal. That
  order is the multi-accumulator reduction `dev/CPP.md` bans, kept by the maintainer's ruling
  because a plain sum moves results through the EM loop (`dev/eigen_conversion.md`, step 5);
  Phase 4 revisits it.

`find`, `fliplr` and `reverse`, which `src/lobpcg_warm.h` called and now writes as loops, move
values and do no arithmetic (NOT READ); the captures of step 4 are bitwise across the change.

## Data structures and ownership

- `obj` crosses as `const Rcpp::NumericMatrix &` and is never mutated in C++; the write-back happens in
  R, through the matrix door (`R/pca_imp.R`: `obj <- fill_parts(obj, parts, keep = keep)`), or not at
  all through the `si_mi` door, whose chunks hold the holes and their values. `row_w` crosses by reference and
  the kernel normalizes its own copy (`src/armaSVD.cpp`: `the kernel's own copy of the weights`).
- `Xhat` is the working copy -- permuted column order, centered, scaled, mutated in
  place every outer iteration. `fittedX` is allocated once (`src/armaSVD.cpp`: `std::vector<double> fittedX(nrX * n_elig);`);
  dgemm overwrites only its changing block each iteration, and the fixed block is
  filled once, after the loop (Layer 5). `U` is mutated in place by `scale_rows_cols_inplace()`, so after
  the reconstruct step it no longer holds the orthonormal left factor.
- `AA_NxN` is destroyed by every eigen call: dsyevr overwrites it
  (`src/eig_sym_sel.h`: `destroys A: dsyevr overwrites`),
  and the LOBPCG path mirrors its upper triangle down first
  (`src/hybrid_topk_eig.h`: `mirror_upper_to_lower(A, a_rows)`).
- `HybridEigContext` carries state across outer iterations: `lobpcg_state` (`X` and `P`),
  `X_prev`, `auto_choice`, the auto counters and timing sums,
  and `n_exact`/`n_lobpcg_ok`/`n_lobpcg_bad` (`src/hybrid_topk_eig.h`: `LOBPCGState lobpcg_state;`).
  When auto picks
  exact, only `lobpcg_state.X`, `lobpcg_state.P` and `X_prev` are reset.
- Every matrix and vector of the kernel is a `std::vector<double>`, column-major and
  unpadded, zero-filled when sized, where `arma::mat` left its cells as allocated.
- `V` is sized once, before the EM loop (`src/armaSVD.cpp`: `std::vector<double> V(n_elig * n_cp);`).
  In the wide case its rows
  `0 .. n_fixed-1` stay zero until postprocessing writes them
  (`src/armaSVD.cpp`: `n_gram, n_cp, V.data());`); nothing reads them inside the
  loop. `Vn` is sized there too (`src/armaSVD.cpp`: `std::vector<double> Vn(n_gram * n_cp);`), written
  only by `split_eigvec_block()` and read only by the
  factor-recovery
  dgemm calls, which is what lets the post-loop fill reuse the last iteration's copy.
- Return crosses as an `Rcpp::List` holding `loc` and `value`, allocated as R objects
  and written in place. Only missing cells are transported, never the full imputed matrix.
  `col_vars()` instead returns by value from `mat_stats.cpp`, a length-`p` vector
  written straight into a SEXP. `mean_imp_col()` reaches no kernel of its own: it builds a mean
  part in R and fills it through the matrix door (`R/parts.R`:
  `fill_parts <- function(obj, parts, keep = NULL) {`), one copy of the matrix.
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
- `src/matrix_checks.h` is included by `src/armaSVD.cpp`, for `col_ptr()` and
  `stop_on_bad_col_index()`, as by
  `src/mat_stats.cpp`, `src/sample_each_rep_cpp.cpp` and `src/imputed_value.h`.
  Besides those it holds `stop_on_inf()`, which the PCA kernel does not call and which
  is reached from `src/mat_stats.cpp` via `check_inf()`, which `mean_imp_col()` calls
  (`R/mean_imp_col.R`: `check_inf(obj)`). `pca_imp()` no longer reaches it: its mean
  fill goes through `fallback_part()` directly, and `check_finite()` has already
  refused `Inf`.

## Compile-time switches

`src/Makevars` and `src/Makevars.win` set no macro: `-DARMA_64BIT_WORD=1`, which made
`arma::uword` 64 bits, left them in Phase 3 of `dev/eigen_conversion.md`. The two below are
added by `load_all1()` (`dev/dev-utils.R`) through the user Makevars.

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
block, which adds the singular values of `V_prev' * V_new` per iteration (`si_la::matmul`, `si_la::svd_values`) to track
subspace rotation; and a wider return list. The counters `n_exact`,
`n_lobpcg_ok`, `n_lobpcg_bad` are always maintained and returned; only the per-iteration
logs are gated.

## Parallelism

- **Inside the PCA kernel: none.** `src/armaSVD.cpp` does not include `RcppThread.h`, and
  neither does any header it pulls in. All intra-call parallelism comes from whatever BLAS
  the R installation links.
- **RcppThread** is reached through `par_for()` (`src/par_for.h`: `inline void par_for(`),
  which only the K-NN kernel calls (`src/impute_knn_brute.cpp`, not in this
  chain); `col_vars()` (`R/col_vars.R`: `col_vars <- function(`), which
  `pca_imp()` does call, runs on one thread and takes no `cores`.
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
