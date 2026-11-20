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
wrappers do"), and every BLAS and LAPACK call of the chain goes through it. No file
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
`pca_imp_internal_cpp` (`src/armaSVD.cpp`: `pca_imp_internal_cpp(const Rcpp::NumericMatrix &obj,`), marked
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
  (`src/armaSVD.cpp`: `has a row for each and an R matrix holds at most %d rows`), since `loc` has a row for each.
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
`criterion_final` is the criterion that stopped a converged fit. When the objective alone
stopped it, the criterion being at or above `threshold`, there is none to report and the
value is `NA` (`src/armaSVD.cpp`:
`(converged && !(criterion < threshold)) ? NA_REAL : criterion;`). The `NA` can withhold a
well-defined ratio, and it also covers one that is rounding noise, as when `ncp` is at or
above the rank of the data from the first iteration and the objective sits at the rounding
floor of the subtraction, with no tolerance needed to tell the two apart.
`solver_chosen` comes from `HybridEigContext::solver_chosen_code()`
(`src/hybrid_topk_eig.h`: `int solver_chosen_code() const`): 0 forced exact, 1 forced lobpcg,
2 auto never decided, 3 auto chose exact, 4 auto chose lobpcg. Under `solver = "auto"` only
(`R/pca_imp.R`: `if (solver == "auto" && i == 1L) {`), `pca_imp()`
collapses 1 and 4 to "lobpcg" and everything else to "exact"; the string is set at
(`R/pca_imp.R`: `solver_chosen = if (resolved_solver_code == 1L)`).

## Layer 5: Gram formation and factor recovery

`SVD_triplet()` (`src/svd_triplet.h`: `inline void SVD_triplet(`) never forms an SVD directly.
It builds the
weighted Gram into `AA_NxN` via `form_weighted_gram()`, stops if its upper triangle
holds NaN or Inf (`src/svd_triplet.h`: `if (!gram_upper_finite(AA_NxN, n))`), the one
such test in the chain and the stop that finite input whose squares overflow gets, takes
the trace of `AA_NxN` for the tail energy (`gram_trace()`), calls `hybrid_topk_eig()`
(a status other than 0 being an `Rcpp::stop`), converts eigenvalues to singular values and
inverses with a `PCA_TOL` floor, and recovers the other factor with one dgemm
through `si_la::gemm` (`src/svd_triplet.h`: `LOC_TIC(timer, "recover_factors")`): tall gives
eigenvectors as `V` and `U = X_work * Vn` (`gemm_nn`, dgemm whatever the shape, never dgemv);
wide gives them as `U` and forms only the changing rows of `V`,
`V[n_fixed.., ] = Xhat[, n_fixed..]' * Vn`, as one strided-block dgemm in
`gemm_tn_block()` (`src/pca_linalg_utils.h`: `gemm_tn_block(const double *A, const std::size_t a_rows,`), reading `n_fixed`
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

`form_weighted_gram()` (`src/gram_ops.h`: `form_weighted_gram(const double *Xhat, const double *sw, const bool tall,`) has four branches,
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
`26 * n` / `10 * n` minima ("What the wrappers do"). `eigvals` and `eigvecs` are
`src/armaSVD.cpp`'s, sized there before the EM loop (`src/armaSVD.cpp`:
`std::vector<double> eigvecs(n_gram * k);`). dsyevr returns
ascending order, so `eig_sym_sel()` reverses eigenvalues and swaps whole eigenvector columns
in place, then `canonicalize_signs()` (`src/pca_linalg_utils.h`: `void canonicalize_signs(`)
applies the sklearn `svd_flip` convention. A non-zero status, from the wrapper or from an
`A` that is not `n_gram` square, is what `hybrid_topk_eig()` returns, and
`SVD_triplet()` stops with a text chosen by it: LAPACK's own `info`
(`src/svd_triplet.h`: `Rcpp::stop("the symmetric eigensolver failed (LAPACK dsyevr, info = %d)",`),
a failed allocation, a non-finite Gram, which says the same as the test before the
solve, or an internal error for the rest. All are plain errors with no class of their own.
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
(`src/lobpcg_warm.h`: `struct Block`), a column-major buffer with its dimensions `nrow` and
`ncol`, its cells in a `std::vector<double>`, and `data()` for the first cell, which is
what the calls into `si_la` are handed.
Every product goes through one function (`src/lobpcg_warm.h`: `inline void product_into(`):
a right operand of one column is a vector and the product is dgemv through `si_la::gemv`,
any other is dgemm through `si_la::gemm`. A product inside a sum or
difference, or one whose target is among its operands, is formed into a block of its own
first and subtracted cell by cell (`src/lobpcg_warm.h`: `sub_cells(V, mul(X, tmul(X, V)));`,
`sub_cells(P, mul(X_curr, tmul(X_curr, X_prev)));`), and of three operands the inner product
is formed first; `-=` of a product whose operands do not include its target is one call with
`beta = 1` (`src/lobpcg_warm.h`: `inline void sub_mul(`). The boundary with
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
| `qr_orthonormalize_update_matvec()` | `si_la::solve_tri` | false, which both callers turn into `LOBPCG_P_COLLAPSED`. The one numerical failure of plain dtrtrs is a zero on the diagonal of `R`, which the pivot test of `qr_orthonormalize()` has already ruled out for a block that kept every column, so no input is known to reach it. Nothing estimates the condition of `R`: a triangle singular to working precision is solved as it is |
| `small_eig_desc()` | `si_la::syevd` | false. From `rr_two_block()` or the `X`-only form before the loop, `LOBPCG_INITIAL_RR_FAILED`; inside the loop, `LOBPCG_RR_FAILED` |

Every `LOBPCG_*` status but `LOBPCG_CONVERGED` leaves the state untouched and
sends `hybrid_topk_eig_impl()` to the exact solver on the same Gram, as above.
"What the wrappers do" has what each wrapper refuses before its routine and what
its status means.

**Where BLAS/LAPACK is entered.**

| symbol | site | purpose |
| --- | --- | --- |
| dsyrk, through `si_la::syrk` | `syrk_upper()`, `form_weighted_gram()` (`src/gram_ops.h`) | Gram formation |
| dgemm, through `si_la::gemm` | `form_weighted_gram()` (`src/gram_ops.h`) | fixed-changing Gram cross block |
| dgemm, through `si_la::gemm` | `gemm_nn()`, `gemm_tn_block()` (`src/pca_linalg_utils.h`) | `gemm_nn` / `gemm_tn_block` factor recovery; the block form also fills the fixed rows of `V` after the loop (`src/armaSVD.cpp`: `pca_detail::gemm_tn_block(Xhat.data(), nrX, n_elig, 0, n_fixed,`) |
| dgemm, through `si_la::gemm` | `reconstruct_cols()` (`src/armaSVD.cpp`) | `reconstruct_cols`, strided column block |
| dsyevr, through `si_la::syevr_range` | `eig_sym_sel()` (`src/eig_sym_sel.h`) | exact top-k eigen, workspace query included |
| dsyevd, through `si_la::syevd` | `small_eig_desc()` (`src/lobpcg_warm.h`) | LOBPCG Rayleigh-Ritz, workspace query included |
| dgeqrf and dorgqr, through `si_la::qr_econ` | `qr_orthonormalize()` (`src/lobpcg_warm.h`) | LOBPCG orthonormalization |
| dtrtrs, through `si_la::solve_tri` | `qr_orthonormalize_update_matvec()` (`src/lobpcg_warm.h`) | `A * Q` carried through a QR without a new matvec |
| dgemm through `si_la::gemm`, or dgemv through `si_la::gemv` for a block of one column | `product_into()` (`src/lobpcg_warm.h`) | LOBPCG matvecs and projections |

The exact path enters LAPACK once per outer iteration with an `n_gram x n_gram` problem. The
LOBPCG path never touches `n_gram`-sized LAPACK: its dense work is `n_gram x k` matvecs and
QRs plus a `3k x 3k` dsyevd.

## What the wrappers do

`src/si_lapack.h` declares, and `src/si_lapack.cpp` defines, R's BLAS and LAPACK as plain
functions over `double*`, `std::size_t` dimensions and an int status. The header's opening
comment has the conventions and each function's comment what it computes; this section has
what each does around its routine, and why. `src/pca_linalg_utils.h` calls `si_la::gemm`,
`src/gram_ops.h` calls `si_la::syrk` and `si_la::gemm`, `src/eig_sym_sel.h` calls
`si_la::syevr_range`, `src/lobpcg_warm.h` calls `si_la::qr_econ`, `si_la::solve_tri`,
`si_la::syevd`, `si_la::gemm` and `si_la::gemv`, and `src/armaSVD.cpp` calls `si_la::gemm`,
and `si_la::svd_values` in the diagnostics build; Layer 6 has where each site enters BLAS and
LAPACK and what it does with a status. `dev/si_lapack_test.R` tests every function against
R's own products, `eigen()`, `svd()`, `qr()`, `backsolve()` and `forwardsolve()`, with guard
cells after every buffer a wrapper writes, and `dev/si_lapack_selftest.R` shows it failing on
a broken wrapper (`dev/eigen_conversion.md`, "The wrapper layer's own test"). READ means read
in the source named; NOT READ is said where it applies.

**Before the routine.** Each wrapper checks its arguments before Fortran is entered. An
option character other than `'N'` or `'T'` for a transpose and `'U'` or `'L'` for a triangle,
a leading dimension below the rows of its buffer, or dimensions the routine would refuse
together are `err_arg`, so no illegal argument reaches xerbla, which under R raises an R error
from inside the library and unwinds no C++ frame. A dimension past `int` is `err_size`, never
narrowed, and so is a workspace whose documented minimum is past `int`, tested before the
routine is asked for its optimum, since LAPACK's own arithmetic for that answer has wrapped by
then (`src/si_lapack.cpp`: `bool fits_int(const double len)`). A failed allocation is
`err_alloc`. A positive status is the routine's own `info`, passed on. Nothing in the file
calls the R API, prints or throws: each call site raises its own error.

| wrapper | routine and options | refused before the routine, beyond the options | calls no routine when |
| --- | --- | --- | --- |
| `si_la::gemm` | dgemm, every argument passed through | `ldc < m`; when `k > 0`, `lda` or `ldb` below the rows its operand is stored with | `m` or `n` is 0, and writes nothing; `k` is 0, and `C` becomes `beta * C` |
| `si_la::gemv` | dgemv, unit strides | `lda < m` | `y` is empty; `x` is empty, and `y` becomes `beta * y`, where dgemv itself would leave `y` alone |
| `si_la::syrk` | dsyrk, into the `uplo` triangle of `C` alone | `ldc < n`; when `k > 0`, `lda` below the rows `A` is stored with | `n` is 0; `k` is 0, and that triangle becomes `beta * C` |
| `si_la::syevr_range` | dsyevr, `jobz = 'V'`, `range = 'I'`, `abstol = 0` | not `1 <= il <= iu <= n`; NaN or Inf in the `uplo` triangle (`err_nonfinite`) | never: `n = 0` has no range |
| `si_la::syevd` | dsyevd, `jobz = 'V'` | NaN or Inf in the `uplo` triangle (`err_nonfinite`) | `n` is 0 |
| `si_la::qr_econ` | dgeqrf, then dorgqr with `k = n` | `m < n` | `n` is 0 |
| `si_la::solve_tri` | dtrtrs, `trans = 'N'`, `diag = 'N'` | nothing | `n` or `nrhs` is 0 |
| `si_la::svd_values` | dgesdd, `jobz = 'N'`; compiled only when `PCA_IMP_DIAGNOSTICS` is set | NaN or Inf anywhere in `a` (`err_nonfinite`) | `min(m, n)` is 0 |

Where a wrapper scales instead of calling, it scales as BLAS defines the result: `beta = 0`
assigns 0, so a NaN or Inf already there does not survive, and `beta = 1` leaves the cells
alone (`src/si_lapack.cpp`: `void scale_range(double *x, const std::size_t len, const double beta)`).
After the routine, `syevr_range` returns `err_count` when dsyevr reports success with another
number of eigenpairs than asked; with NaN and Inf refused first, no input is known to reach it.
On any failure what the outputs hold is unspecified; `dev/si_lapack_test.R` checks that a
refused call leaves them as they were.

**The eigenvalue buffer of `syevr_range`.** LAPACK documents `W` as `n` cells whatever the
range. For a proper subset dsyevr goes through dstebz, which stores every eigenvalue its search
interval holds and then drops those outside `il..iu`; when the range cuts a group of equal
eigenvalues that is more than `iu - il + 1` of them. The wrapper holds the `n` cells itself
and copies `iu - il + 1` out (`src/si_lapack.cpp`:
`std::copy(w_all.begin(), w_all.begin() + count, w);`), so its caller's buffer needs no more.
`pca_imp()` reaches such a range when the completed data is rank deficient and `ncp`, plus one
when regularized, exceeds its rank: samples that are exact copies of each other, holes
included, or columns that are. A call that hands dsyevr a buffer of `iu - il + 1` cells, as
the chain did before it moved to this wrapper and 1.2.0 does, writes past it, and the write
is not harmless. Measured on Windows (`dev/review_p2.md`, finding 1): on 30 samples of which 6
are distinct, 1.2.0 built from source ends the R process with a segmentation fault,
`pca_imp(x, ncp = 5)` at its defaults included. The regression test is that input
(`tests/testthat/test-pca_imp.R`: `test_that("pca_imp returns on replicate samples with ncp above their rank", {`).
It sees the overrun only where it does harm; `dev/si_lapack_test.R` sees it in guard cells
after `w`, on ranges that cut a group, where it also shows dsyevr writing past the count.

**NaN and Inf.** LAPACK does not refuse them. Measured on one machine with R's own LAPACK
(`dev/eigen_conversion.md`, "The review of the wrapper layer"): dsyevr asked for every
eigenpair of a matrix holding NaN does not return, a proper subset comes back with `info = 0`
and another number of eigenpairs than asked, and dsyevd can return `info = 0` with NaN
eigenvalues. So `syevr_range` and `syevd` test the triangle LAPACK reads, and only that one,
before the routine is entered (`src/si_lapack.cpp`:
`bool tri_finite(const char uplo, const std::size_t n, const double *a)`), and `svd_values`
tests all of `a`. `qr_econ` and `solve_tri` test nothing, and NaN or Inf goes through to their
output. In the chain, `pca_imp()` refuses NaN and Inf in its input (`R/pca_imp.R`:
`check_finite(obj)`), and the Gram is tested where it is formed (Layer 5), so finite input
whose squares overflow, unscaled values near 1e160 or near 1e200 with `scale = TRUE`
(`dev/review_p2.md`, finding 3), stops there with a message that names the cause
(`tests/testthat/test-pca_imp.R`: `test_that("finite input whose squares overflow stops before the eigensolver", {`).
A non-finite Gram that reached the exact path would still be refused by `syevr_range`. On the
LOBPCG path a small Gram holding NaN or Inf is a failed Rayleigh-Ritz, and a residual norm
that is NaN or Inf ends the solve as `LOBPCG_BAD_SEED` before the convergence test
(`src/lobpcg_warm.h`: `if (!lobpcg_detail::all_finite(rn))`), since the largest norm is taken
by `max_skip_nan()`, which passes over NaN; either way that outer iteration goes to the exact
solver.

**Workspace.** Every wrapper that needs a workspace asks its routine for it on every call,
with a length of -1, and uses the larger of the answer and the routine's documented minimum
(`src/si_lapack.cpp`: `bool work_len(const double queried, const double minimum, int &out)`):
`26 n` and `10 n` for dsyevr, `1 + 6 n + 2 n^2` and `3 + 5 n` for dsyevd, `max(1, n)` for
dgeqrf and for dorgqr, and `3 min(m, n) + max(max(m, n), 7 min(m, n))` for dgesdd. An answer
that is not usable, 0, negative or NaN, leaves the minimum. Nothing is kept between calls, so
the length is a function of the problem alone. That matters, since the length can change the
result:

- dsyevr: given `26 n` in place of its answer, the eigenvectors change at orders 60 and 300
  (`dev/eigen_conversion.md`, "The review of the wrapper layer").
- dsyevd: READ, in LAPACK 3.12.0 (`dsyevd.f`, `dormql.f`, fetched from the reference
  repository, not from R's sources): dsyevd keeps `1 + 2 n + n^2` cells and hands the rest to
  dormtr, which is dormql for `uplo = 'U'`; dormql applies the reflectors in blocks only when
  its block size is below `n - 1` and that remainder holds its `T` array of 4160 cells and at
  least two columns a block, and otherwise one at a time. At the queried length,
  `1 + 6 n + 2 n^2`, the remainder is `n^2 + 4 n + 1`. With a block size of 32 (NOT READ:
  ilaenv; it is the value the measured boundaries give) that makes the eigenvectors depend on
  the length for orders 34 to 79 and for no other. Measured, `dev/probe_p2s4_syevd_lwork.R`,
  orders 2 to 140: with the workspace of a problem 1.5, 2 or 3 times the order, every call of
  order 34 to 79 differs from the wrapper's and no other does; with that of the order plus
  one, orders 62 to 79; the eigenvalues never differ, and the eigenvectors by at most 6.9e-16.
  The chain once kept the longest dsyevd workspace so far on the LOBPCG state;
  `dev/eigen_conversion.md`, Phase 2, step 4, has what giving each call its own did to
  results.

**`qr_econ`.** One form, `m >= n`: `x` is copied into `q`, dgeqrf factors it there, the upper
triangle it leaves is copied into the `n x n` `r` with zeros below, and dorgqr forms the `n`
columns of `q` with `k = n`. The sign of each column is LAPACK's: the diagonal of `r` is not
made positive. `m < n` is `err_arg`, a shape the one call site stops on before the call
(`src/lobpcg_warm.h`: `internal error: a LOBPCG block has more columns than rows`). The call
(`src/lobpcg_warm.h`: `if (la_failed(si_la::qr_econ(v_rows, v_cols, V.data(), Q.data(),`)
and what a failure there leads to are in Layer 6.

**`solve_tri`.** dtrtrs alone, which is how R's own `qr.coef()` solves a triangle (READ:
`qr_coef_real` in R's `Lapack.c`, dtrtrs with `info != 0` an error). A positive status `i` is
dtrtrs's: diagonal element `i` of the triangle is exactly zero, and `b` is as it was given.
Nothing estimates the condition of the triangle, so one singular to working precision is
solved as it is, and NaN or Inf in either operand goes through to the answer. The one call
site solves `R' Y = (A Q)'` with the two transposes written out (`src/lobpcg_warm.h`:
`if (la_failed(si_la::solve_tri('L', Rt.nrow, Y.ncol, Rt.data(),`).

**`svd_values`.** The diagnostics build only: it is compiled when `PCA_IMP_DIAGNOSTICS` is
set, so the package as shipped does not reference dgesdd. Its one caller forms `V_prev' V_new`
with `si_la::gemm`, once per iteration from the second on, and hands it over
(`src/armaSVD.cpp`: `st = si_la::svd_values(n_cp, n_cp, V_cross.data(), V_cross_sv.data());`).
A positive status or `err_nonfinite` is an R error, as an `Rcpp::exception`
(`src/armaSVD.cpp`: `the SVD of the subspace product failed`); `err_alloc` and a refusal of
the arguments stop too, the second as an internal error. The smallest singular value is taken
with `std::min_element`. No route to the failure is known, since a Gram holding NaN is refused
before the eigenvectors exist; `dev/probe_p2s6_run.sh` forces it.

**Sums outside the wrappers.** The chain's other reductions call no BLAS. Each is a plain loop
in order; `dev/CPP.md` bans the multi-accumulator form some of them once had, and
`dev/eigen_conversion.md` has what each moved when it took this form ("Not wrapped" for the
two norms of LOBPCG, "Phase 4, as carried out" for the three sums of the EM loop):

- The Frobenius norm of `P` in `build_momentum_P()`, a sum of squares, read only by a
  threshold (`src/lobpcg_warm.h`: `!(frobenius_norm(P) > min_projected_fro_norm))`).
- The residual norms of LOBPCG, one sum of squares a column (`src/lobpcg_warm.h`:
  `inline std::vector<double> column_norms(`). They feed the iteration's tests and
  `max_rel_res`, which only the diagnostics build reports (`src/armaSVD.cpp`:
  `Rcpp::Named("eig_lobpcg_max_rel_res")`).
- The infinity norm of the Gram, which scales LOBPCG's tolerance, `norm_inf()`: each row's sum
  of absolute values, from zero in column order (`src/lobpcg_warm.h`:
  `row_sum[i] += std::abs(col[i]);`), then the largest.
- The largest of a set of norms, `max_skip_nan()`: from `-Inf`, keeping a value only when it
  is `>` the best so far, so NaN is skipped and a set of NaN gives `-Inf`.
- The sum of the squared leading singular values that the tail energy subtracts from the
  trace (`src/armaSVD.cpp`: `sum_top_sq += vs_top[cp] * vs_top[cp];`), and the trace itself
  (`src/svd_triplet.h`: `static inline double gram_trace(`).
- The sum the row weights are divided by, taken in the loop that tests their sign
  (`src/armaSVD.cpp`: `row_w_sum += w[i];`).

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
block, which adds the singular values of `V_prev' * V_new` per iteration (`si_la::gemm`, `si_la::svd_values`) to track
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
