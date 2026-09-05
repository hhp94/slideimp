# knn_imp() structural map

A file-and-line map of the K-NN imputation path, from the exported R function
down to the C++ kernel. Descriptive only: what the code is and how it connects.

Orientation: `knn_imp()` documents rows as samples and columns as features, and
the neighbor search runs **between columns** -- a "neighbor" is another column of
`obj`, and the distance loop walks rows. Hence `k` is capped at `ncol(obj) - 1`
(`R/knn_imp.R:84`), and `tests/testthat/test-knn_imp.R:64-66` transposes when
comparing against `impute::impute.knn`, which searches rows.

## Call chain overview

```
R/knn_imp.R:61                      knn_imp()  -- user entry, validation, partition
  +- src/mat_stats.cpp:197          check_finite()  (Inf + all-NA column gate)
  +- R/utils.R:56                   resolve_subset()  -- names/ints -> col indices
  +- R/mat_miss.R:33                mat_miss()
  |    +- src/mat_stats.cpp:241     col_miss_internal() via R/RcppExports.R:44
  +- [abort site 1] R/knn_imp.R:115-123   k > n_elig - 1
  +- [abort site 2] R/knn_imp.R:134-142   all subset cols exceed colmax
  +- R/RcppExports.R:16             impute_knn_brute()  -- .Call shim
  |    +- src/RcppExports.cpp:69    _slideimp_impute_knn_brute (SEXP wrapper)
  |         +- src/impute_knn_brute.cpp:365  impute_knn_brute()   <-- C++ entry
  |              +- src/imputed_value.h:35        validate_knn_inputs()
  |              +- src/matrix_checks.h:8         stop_on_inf()
  |              +- src/impute_knn_brute.cpp:387  copy_with_mask lambda
  |              +- src/imputed_value.cpp:18      initialize_result_matrix()
  |              +- RcppThread::parallelFor       (src/impute_knn_brute.cpp:450)
  |                   +- src/impute_knn_brute.cpp:279  distance_vector()
  |                   |    +- :183  distance_vector_impl<Metric>()
  |                   |         +- :41   calc_distance_raw<Metric, Bound>()
  |                   |         +- :97   calc_distance_raw_complete<Metric, Bound>()
  |                   |         +- :153  insert_before_k()
  |                   |         +- :159  insert_if_better_than_worst()
  |                   +- src/impute_knn_brute.cpp:305  knn_weights()
  |                   +- src/imputed_value.cpp:178     impute_column_values()
  +- R/knn_imp.R:160-163            NaN -> NA, scatter triplets back into obj
  +- R/mean_imp_col.R:34            mean_imp_col()  (only when post_imp = TRUE)
  |    +- src/mat_stats.cpp:190     check_inf()
  |    +- src/mat_stats.cpp:78      mean_imp_col_internal()
  +- R/utils.R:5                    new_slideimp_results()
```

## Layer 1 -- R entry point (`R/knn_imp.R`)

Signature at `R/knn_imp.R:61-72`. Validation block, in order
(`R/knn_imp.R:74-90`): numeric matrix with `>= 1` row and `>= 2` columns
(`:74-81`); `check_finite(obj)` (`:82`); `method` matched against
`"euclidean"` / `"manhattan"` (`:83`); `k` integer in `[1, ncol(obj) - 1]`
(`:84`); `cores >= 1` (`:85`); `colmax` in `[0, 1]` (`:86`); `dist_pow` scalar,
`>= 0`, not infinite (`:88`); the `post_imp`, `.progress` and `na_check` flags
(`:87`, `:89-90`).

Missing-value gating is delegated to C++. `check_finite()`
(`src/mat_stats.cpp:197`) makes one pass per column and refuses twice: any
`Inf`/`-Inf` aborts with the row and column of the offender
(`src/mat_stats.cpp:218-223`), and a column ending the pass with no finite cell
aborts with "All NA/NaN column detected" (`src/mat_stats.cpp:233-236`).
`check_inf()` (`src/mat_stats.cpp:190`) is the weaker Inf-only variant, called
from `mean_imp_col()` at `R/mean_imp_col.R:38` but not from `knn_imp()`. Both
are declared in `R/RcppExports.R:36-42`. `NA` and `NaN` are indistinguishable
throughout; the policy is stated at `R/slideimp-package.R:7-16`.

`subset` (`R/knn_imp.R:92`) goes through `resolve_subset()` (`R/utils.R:56-99`),
which accepts `NULL` (all columns, `:60-61`), a character vector matched against
`colnames(obj)` with unmatched names dropped and reported (`:62-77`), or
integerish indices (`:78-89`). It returns `NULL` for "nothing to do" (`:91-96`),
which becomes an early return of the untouched input (`R/knn_imp.R:93-96`).

Missingness and partitioning. `cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)`
at `R/knn_imp.R:99`; `mat_miss()` (`R/mat_miss.R:33`) calls `col_miss_internal()`
(`src/mat_stats.cpp:241`), which counts `std::isnan` per column and does not
count `Inf` as missing (`src/mat_stats.cpp:256-259`). A second early return fires
if no subset column has any missing value (`R/knn_imp.R:102-107`). `colmax` then
becomes an eligibility mask, `eligible <- miss_rate < min(colmax, 1)`
(`R/knn_imp.R:112`); a column at or above the threshold is dropped from the K-NN
stage entirely -- neither imputed nor usable as a neighbor. Three disjoint index
groups follow at `R/knn_imp.R:125-132`: `grp_impute` (eligible, has missing
values, in `subset` -- these get imputed), `grp_miss_no_imp` (eligible, has
missing values, not in `subset` -- neighbor candidate only), and `grp_complete`
(eligible, fully observed -- neighbor candidate only). `method` is remapped to
the integer the kernel expects, `0` = euclidean, `1` = manhattan
(`R/knn_imp.R:144`). So `subset` narrows only the *target* set: every eligible
column stays a neighbor candidate regardless of `subset`.

`rowmax` is not an argument of `knn_imp()`: it belongs to the tuning path
(`R/tune_imp.R:55`, `:184`), capping injected missingness per row.

After the `.Call`, `R/knn_imp.R:160` maps `NaN` back to `NA_real_` (the kernel
emits `NaN` for cells it could not fill), `R/knn_imp.R:162-163` uses the returned
`(row, col)` pairs as a two-column index matrix to scatter values into `obj`,
`R/knn_imp.R:166-168` optionally runs column-mean fill over `subset`, and
`R/knn_imp.R:170-178` attaches the `slideimp_results` class and attributes via
`new_slideimp_results()` (`R/utils.R:5-19`).

## Layer 2 -- the .Call boundary

`R/RcppExports.R:16-18` defines the R-side shim `impute_knn_brute(obj, k,
grp_impute, grp_miss_no_imp, grp_complete, method, dist_pow, cores, pb)`, whose
body is a single `.Call` to `_slideimp_impute_knn_brute`.

The generated C wrapper is `src/RcppExports.cpp:69-85`; it converts each SEXP
through `Rcpp::traits::input_parameter<...>` (`:73-81`) and calls the C++ symbol
at `:82`. Registration is `:195` inside the `CallEntries` table, installed by
`R_init_slideimp` at `:208`. The exact C++ entry symbol is `impute_knn_brute`,
defined at `src/impute_knn_brute.cpp:365` under the `// [[Rcpp::export]]` marker
on `:364`. It is the only C++ entry point the K-NN stage needs; `check_finite`,
`col_miss_internal`, and `mean_imp_col_internal` are separate `.Call`s made from
R around it.

## Layer 3 -- C++ translation units and headers on the path

Four compiled `.cpp` files are reached: `src/RcppExports.cpp` (generated SEXP
wrappers), `src/impute_knn_brute.cpp` (K-NN kernel and entry point),
`src/imputed_value.cpp` (result layout and accumulation loop), and
`src/mat_stats.cpp` (`check_finite`, `check_inf`, `col_miss_internal`,
`mean_imp_col_internal`).

Project headers reached, following the `#include` graph:

- `src/imputed_value.h` -- included by `src/impute_knn_brute.cpp:1` and
  `src/imputed_value.cpp:1`. Defines `mask_t` / `MaskMat` (`:11-12`),
  `GroupLayout` (`:24-33`), the header-inline `validate_knn_inputs()`
  (`:35-110`), and the two out-of-line declarations (`:116`, `:128`).
- `src/matrix_checks.h` -- pulled in by `src/imputed_value.h:8` and
  `src/mat_stats.cpp:4`. Holds only `stop_on_inf()` (`src/matrix_checks.h:8-29`).
- `src/loc_timer.h` -- included at `src/impute_knn_brute.cpp:2`. Every macro
  compiles to `((void)0)` unless `LOC_TIMER` is defined (`:26-34`), so the
  `LOC_*` calls at `src/impute_knn_brute.cpp:447-448` and `485-486` vanish in a
  normal build.

Plus the external `RcppArmadillo.h` (`src/imputed_value.h:4`) and `RcppThread.h`
(`src/impute_knn_brute.cpp:5`). The remaining `src/` headers (`eig_sym_sel.h`,
`gram_ops.h`, `hybrid_topk_eig.h`, `lobpcg_warm.h`, `pca_linalg_utils.h`,
`svd_triplet.h`) and `.cpp` files (`armaSVD.cpp`, `find_windows.cpp`,
`find_windows_flank.cpp`, `sample_each_rep_cpp.cpp`) serve the PCA, windowing,
and simulation paths and are not reached from `knn_imp()`.

`src/Makevars:1` sets `-DARMA_64BIT_WORD=1`, so `arma::uword` is 64-bit and
kernel index arithmetic does not wrap on matrices past `2^32 - 1` elements;
`src/mat_stats.cpp:180` exposes `arma_uword_bytes()` so the test suite can assert
that flag survived.

## Layer 4 -- kernel setup and result skeleton

1. `validate_knn_inputs()` (`src/imputed_value.h:35`) re-checks everything R
   already checked, so the kernel is safe when called directly: shape (`:44-52`),
   `k >= 1` (`:54-57`), `method` in `{0, 1}` (`:59-62`), `dist_pow` finite and
   non-negative (`:64-67`), every group index in range (`:69-88`), and
   `k <= n_working - 1` for `n_working` = sum of the three group sizes
   (`:90-109`). Then `stop_on_inf(obj)` at `src/impute_knn_brute.cpp:378`
   rescans for infinities.
2. `GroupLayout layout{...}` at `src/impute_knn_brute.cpp:379` is the single
   source of truth for group boundaries: groups 1 and 2 occupy local columns
   `[0, n_masked())` of the working matrices, and `complete_start()`
   (`src/imputed_value.h:31`) is a *virtual* index used as a tag -- a neighbor id
   `>= complete_start()` means "group 3, offset `id - complete_start()` into
   `grp_complete`".
3. Working buffers `obj_masked` (`arma::mat`), `nmiss_masked`
   (`arma::Mat<uint8_t>`) and `n_col_valid` are allocated at
   `src/impute_knn_brute.cpp:383-385` and filled by `copy_with_mask`
   (`:387-405`), which calls `std::isnan` once per cell and reuses the result as
   the mask byte, the count increment, and the zeroing test (`:399-402`). Missing
   cells therefore hold a literal `0.0` in `obj_masked` and are excluded by the
   mask, never by a branch in the distance loop.
4. Group 1 is copied at `src/impute_knn_brute.cpp:408-411`, group 2 at `:425-428`.
   Group 3 is **not** copied: the kernel reads it straight out of `obj` through
   `grp_complete` (`:430`, and the lambdas at `:220-234`).
5. `initialize_result_matrix()` (`src/imputed_value.cpp:18-165`) runs between the
   group-1 and group-2 copies (`src/impute_knn_brute.cpp:416-417`) so the call
   can bail early when there is nothing to impute (`:419-422`).

**Result skeleton.** The return value is a triplet matrix: 1-based row, original
1-based column, value. Missing counts per target column come from `n_col_valid`
rather than a second pass over the mask (`src/imputed_value.cpp:90`), and a
prefix sum builds `col_offsets` (`:71-99`) so each target column owns a
contiguous slice of result rows and threads never collide. The matrix is filled
with `NaN` (`:107-108`), columns 0 and 1 are populated up front (`:152-161`), and
`rows_to_impute_vec[i]` records which rows of target column `i` are missing
(`:112-149`). Empty cases short-circuit at `:33-36` and `:102-105`, and
`Rcpp::stop` consistency guards surround the bookkeeping (`:39-67`, `:78-88`,
`:124-131`, `:137-147`).

## Layer 5 -- neighbor search (`src/impute_knn_brute.cpp:41-302`)

**Metrics.** `EuclideanMetric` (`:18-21`) accumulates `diff * diff`,
`ManhattanMetric` (`:23-26`) accumulates `std::abs(diff)`. The metric is a
template parameter and `distance_vector()` (`:279`) dispatches once per target
column on `method` (`:290-301`), so the inner body is monomorphized.

**Distance over incomplete data.** `calc_distance_raw<Metric, Bound>` (`:41-88`)
handles a pair where either side may be missing. Per row it forms
`valid = target_nmiss[r] & other_nmiss[r]` (`:60`) and adds
`valid * Metric::accumulate(diff)` (`:62`) -- branchless, and correct because
missing values were zeroed during the masked copy. It returns `dist / n_valid`
(`:87`), the mean per-coordinate contribution over rows observed in *both*
columns, or `arma::datum::inf` when the two columns share no observed row
(`:82-85`). No square root is taken for the euclidean case, so the stored
quantity is a mean squared difference.
`calc_distance_raw_complete<Metric, Bound>` (`:97-140`) is the group-3 variant:
the other side is fully observed, so only the target mask matters (`:122`) and
`n_valid` is known up front and passed in (`:104`).

**Pruning.** Both kernels chunk rows in blocks of `GRAIN = 16` (`:29`) and, when
`Bound` is true, test a partial-sum bound at each chunk boundary, returning `inf`
on a prune (`:66-72`, `:125-131`; derivations in the comments at `:34-39`,
`:93-95`). `Bound = false` instantiations pass the default `worst_dist = +inf`
(`:48`, `:104`), making the check vacuous.

**Top-k selection.** `distance_vector_impl` (`:183-276`) keeps a
`std::vector<NeighborInfo>` of capacity `k` (`:201-202`) and runs two phases with
three shared cursors `p1`, `p2`, `c` (`:236-238`) so every candidate column is
visited exactly once:

- *Fill* (`:241-252`): append the first `k` candidates with no pruning via
  `insert_before_k` (`:153-156`), in the order group 3, masked columns below the
  target, masked columns above it. `p2` starts at `index + 1` (`:237`) and the
  `p1` loop stops at `index`, so the target column is never its own neighbor.
  Then one full `std::sort` on those `k` entries (`:254-256`).
- *Replacement* (`:259-273`): resume each cursor, compute the distance with
  `Bound = true` against the current worst (`top_k.back().distance`), and hand it
  to `insert_if_better_than_worst` (`:159-172`), which rejects anything not
  strictly better than the worst (`:161-164`), overwrites the last slot, and
  bubbles it into position (`:168-171`).

The structure is therefore a sorted fixed-size array with insertion, not a heap
and not a full sort of all candidates: exactly one `std::sort`, over `k` elements.

**Ties.** Both comparisons are strict `<` (`:168` for the bubble, and the
`dist >= top_k.back().distance` early return at `:161`), so an equal-distance
candidate never displaces an incumbent: ties resolve in favor of the
first-encountered column under the fixed scan order above, independent of thread
count.

## Layer 6 -- weights and imputation

**Weights.** `knn_weights()` (`src/impute_knn_brute.cpp:305-359`) returns a
vector of ones when `dist_pow == 0` (`:311-314`) -- the plain unweighted mean of
neighbors. Otherwise it finds the smallest strictly positive finite distance
`min_pos` and notes whether any distance is zero (`:316-330`); if no positive
finite distance exists, weights stay equal (`:333-336`). Two live branches: when
some distance is exactly zero (`:338-348`), `floor = min_pos * sqrt(DBL_EPSILON)`
and `w_j = (floor / max(d_j, floor)) ^ dist_pow`, so duplicate columns get
weight 1 and everything else a tiny positive weight; otherwise (`:350-356`)
`w_j = (min_pos / d_j) ^ dist_pow`. Either way the largest weight is 1, and a
neighbor at infinite distance gets weight 0.

**Accumulation.** `impute_column_values()` (`src/imputed_value.cpp:178-269`)
splits the neighbor ids into a masked bucket and a complete bucket by comparing
against `complete_start()` (`:209-219`). Masked neighbors contribute
`w * nmiss[row]` to both the weighted sum and the weight total (`:222-237`), so a
neighbor that is itself missing in a given row drops out for that row only.
Complete neighbors add `w * value` per row and a single shared
`total_complete_w` to every row's weight total (`:243-262`). The final value is
`weighted_sum / weight_total`, or `NaN` when the total is zero (`:264-268`),
which R converts to `NA` at `R/knn_imp.R:160`.

## Feasibility and abort paths

`slideimp_infeasible` is a cli condition class raised by the R layer so the
wrappers can catch a per-group or per-window failure and fall back instead of
dying. `knn_imp()` has exactly two sites:

1. `R/knn_imp.R:115-123` -- `k > n_elig - 1`, where `n_elig` is the number of
   columns passing `colmax` (`R/knn_imp.R:113`); message "`k` (...) exceeds
   usable columns (...)". Exercised at `tests/testthat/test-knn_imp.R:126-128`
   and `tests/testthat/test-group_imp.R:702-704`.
2. `R/knn_imp.R:134-142` -- `length(grp_impute) == 0`, i.e. every subset column
   with missing values also exceeds `colmax`; message "All subset columns with
   missing values exceed `colmax` (...)". Exercised at
   `tests/testthat/test-group_imp.R:706-714`.

Everything else that aborts here is a plain error, not `slideimp_infeasible`:
`check_finite()` for `Inf` and all-`NA` columns
(`tests/testthat/test-knn_imp.R:122-123` and `131-137`), plus the C++ guards in
`validate_knn_inputs()` and `initialize_result_matrix()`. Catch sites are in the
wrappers only (next section); each switches on `on_infeasible` to rethrow, skip
the block unchanged, or fall back to `mean_imp_col()`.

## Data structures and ownership

- The R matrix is deep-copied at the boundary:
  `Rcpp::traits::input_parameter<const arma::mat&>` (`src/RcppExports.cpp:73`)
  resolves to a `ConstReferenceInputParameter` holding an `arma::mat` built by
  Rcpp's `MatrixExporter`, which allocates and copies element by element. The
  kernel's `const arma::mat& obj` references that copy, not R's memory, and is
  never written to.
- The kernel returns a freshly allocated `(n_missing x 3)` `arma::mat` of
  triplets, not an imputed matrix. The only in-place mutation of the user's
  matrix happens in R, at `R/knn_imp.R:163`.
- Groups 1 and 2 are materialized a second time into `obj_masked` /
  `nmiss_masked` (`src/impute_knn_brute.cpp:383-384`) with NaNs zeroed and a
  parallel byte mask; `nmiss_masked` is `arma::Mat<uint8_t>`
  (`src/imputed_value.h:11-12`), one byte per cell rather than a `umat`. Group 3
  is read in place from `obj`, so a matrix dominated by complete columns pays no
  second copy.
- Column indices cross the boundary 0-based (`R/knn_imp.R:150-152` subtracts 1)
  and come back 1-based (`src/imputed_value.cpp:159-160` adds 1).
- Per-thread scratch is small and local: `top_k`
  (`src/impute_knn_brute.cpp:201`), `nn_columns` and `weights` (`:468-473`), and
  the two accumulators in `impute_column_values` (`src/imputed_value.cpp:196-197`).
  `mean_imp_col_internal()` (`src/mat_stats.cpp:78`) is likewise out-of-place: a
  fresh matrix with untouched columns `memcpy`d (`src/mat_stats.cpp:136-140`),
  and `R/mean_imp_col.R:54` reattaches dimnames.

## Parallelism

`RcppThread::parallelFor` at `src/impute_knn_brute.cpp:450-483` runs over
`[0, layout.n_imp)` -- one task per target column. `n_threads` and `n_batches`
are both set to `cores`, floored at 1 (`src/impute_knn_brute.cpp:439-441`, passed
at `:483`), so the index range splits into as many contiguous batches as threads.

All threads read `obj`, `obj_masked`, `nmiss_masked`, and `grp_complete` without
writing them, and each writes only rows `[col_offsets(i), col_offsets(i+1))` of
`result`, a disjoint slice per target column, so no locks are taken. Interruption
is checked with `RcppThread::checkUserInterrupt()` on every fifth index
(`src/impute_knn_brute.cpp:455-458`), which RcppThread routes to the main thread
to unwind the pool. The optional progress bar is an `RcppThread::ProgressBar` in
a `unique_ptr` created only when `pb` is true (`:442-446`) and incremented from
the workers at `:478-481`.

`src/Makevars:2-3` adds the OpenMP flags and links `RcppThread::LdFlags()`;
`src/Makevars.win:2-3` omits the RcppThread link line. `mean_imp_col_internal()`
uses the same pattern (`src/mat_stats.cpp:126-170`), reached when
`post_imp = TRUE`.

## How the wrappers reach knn_imp()

`group_imp()` picks the function object at `R/group_imp.R:811`
(`imp_fn <- if (is_knn_mode) knn_imp else pca_imp`), pre-builds a per-group
parameter list injecting `cores` and a group-local `subset` (`:803-806`), and
invokes it via `do.call()` inside a `tryCatch` -- `:859-882` in the mirai branch,
`:912-932` in the sequential branch. Each group is a column slice
`obj[, indices[[i]]$col_idx]` (`:857`, `:910`), so `knn_imp()` sees a submatrix
and group-local indices.

`slide_imp()` calls `knn_imp()` directly at `R/slide_imp.R:440-451` for each
window `obj[, start[i]:end[i]]` (`R/slide_imp.R:434-435`), passing
`na_check = FALSE`, `.progress = FALSE`, and a window-local `subset`
(`R/slide_imp.R:450`), wrapped in the same `slideimp_infeasible` handler
(`R/slide_imp.R:475-485`). Both wrappers set `na_check = FALSE`
(`R/group_imp.R:800`) and do the final NA accounting themselves.

## Test entry points

`tests/testthat/test-knn_imp.R:1-45` calls `impute_knn_brute()` directly and
asserts the returned `(row, col)` triplet locations match `which(is.na(...))`;
`:77-112` covers `subset` by index and by name with `post_imp` both ways;
`:114-129` covers all-NA rows and columns plus abort site 1; `:131-138` covers
`Inf` rejection; `tests/testthat/test-group_imp.R:700-715` pins the
`slideimp_infeasible` class on both abort sites.
