# knn_imp() structural map

A file-and-line map of the K-NN imputation path, from the exported R function
down to the C++ kernel. Descriptive only: what the code is and how it connects.

Orientation: `knn_imp()` documents rows as samples and columns as features, and
the neighbor search runs **between columns** -- a "neighbor" is another column of
`obj`, and the distance loop walks rows. Hence `k` is capped at `ncol(obj) - 1`
(`R/knn_imp.R:121`), and `dev/probe_impute_knn.R` transposes both sides when
comparing against `impute::impute.knn`, which searches rows. That comparison is
real and asserted -- it agrees to 2.2e-16 absolute and 4.5e-16 relative, not
bitwise. It lives in `dev/` rather than under `tests/` because `impute` is a
Bioconductor package outside `Suggests`, and any `impute::` call under `tests/`
makes `R CMD check` report an unstated dependency regardless of skips; the
reasoning is repeated at `tests/testthat/test-knn_imp.R:59-65` so a future
session does not put it back.

## Call chain overview

```
R/knn_imp.R:98                      knn_imp()  -- user entry, validation, partition
  +- src/mat_stats.cpp:198          check_finite()  (Inf + all-NA column gate)
  +- R/utils.R:62                   resolve_subset()  -- names/ints -> col indices
  +- R/mat_miss.R:33                mat_miss()
  |    +- src/mat_stats.cpp:242     col_miss_internal() via R/RcppExports.R:44
  +- [abort site 1] R/knn_imp.R:175-183   k > n_elig - 1
  +- [abort site 2] R/knn_imp.R:194-202   all subset cols exceed colmax
  +- R/RcppExports.R:16             impute_knn_brute()  -- .Call shim
  |    +- src/RcppExports.cpp:69    _slideimp_impute_knn_brute (SEXP wrapper)
  |         +- src/impute_knn_brute.cpp:397  impute_knn_brute()   <-- C++ entry
  |              +- src/imputed_value.h:36        validate_knn_inputs()
  |              +- src/matrix_checks.h:8         stop_on_inf()
  |              +- src/impute_knn_brute.cpp:419  copy_with_mask lambda
  |              +- src/imputed_value.cpp:18      initialize_result_matrix()
  |              +- par_for()                     (src/impute_knn_brute.cpp:482)
  |                   +- src/impute_knn_brute.cpp:310  distance_vector()
  |                   |    +- :192  distance_vector_impl<Metric>()
  |                   |         +- :46   calc_distance_raw<Metric, Bound>()
  |                   |         +- :102  calc_distance_raw_complete<Metric, Bound>()
  |                   |         +- :158  insert_before_k()
  |                   |         +- :164  insert_if_better_than_worst()
  |                   +- src/impute_knn_brute.cpp:336  knn_weights()
  |                   +- src/imputed_value.cpp:178     impute_column_values()
  +- R/knn_imp.R:220-223            NaN -> NA, scatter triplets back into obj
  +- R/mean_imp_col.R:34            mean_imp_col()  (post_imp, and only if
  |                                 anything is still missing)
  |    +- src/mat_stats.cpp:191     check_inf()
  |    +- src/mat_stats.cpp:79      mean_imp_col_internal()
  +- R/utils.R:5                    new_slideimp_results()
```

## Layer 1 -- R entry point (`R/knn_imp.R`)

Signature at `R/knn_imp.R:98-109`. Validation block, in order
(`R/knn_imp.R:110-132`): numeric matrix with `>= 1` row and `>= 2` columns
(`:111-118`); `check_finite(obj)` (`:119`); `method` matched against
`"euclidean"` / `"manhattan"` (`:120`); `k` integer in `[1, ncol(obj) - 1]`
(`:121`); `cores >= 1` (`:122`); `colmax` in `[0, 1]` (`:123`); `dist_pow` a
finite number `>= 0` (`:125-130`); the `post_imp`, `.progress` and `na_check`
flags (`:124`, `:131-132`). All but two go through `checkmate` with an explicit
`.var.name`, so a rejection names the user's argument. The exceptions are `method`,
matched with `match.arg()`, and `check_finite()`, which raises from C++ through
`Rcpp::stop` (`src/mat_stats.cpp:219-224`). `tests/testthat/test-knn_imp.R:538-551`
pins the `checkmate` path for `dist_pow`.

Missing-value gating is delegated to C++. `check_finite()`
(`src/mat_stats.cpp:198`) makes one pass per column and refuses twice: any
`Inf`/`-Inf` aborts with the row and column of the offender
(`src/mat_stats.cpp:219-224`), and a column ending the pass with no finite cell
aborts with "All NA/NaN column detected" (`src/mat_stats.cpp:234-237`).
`check_inf()` (`src/mat_stats.cpp:191`) is the weaker Inf-only variant, called
from `mean_imp_col()` at `R/mean_imp_col.R:38` but not from `knn_imp()`. Both
are declared in `R/RcppExports.R:36-42`. `NA` and `NaN` are indistinguishable
throughout; the policy is stated at `R/slideimp-package.R:7-16`.

`subset` (`R/knn_imp.R:134`) goes through `resolve_subset()` (`R/utils.R:62-105`),
which accepts `NULL` (all columns, `:66-67`), a character vector matched against
`colnames(obj)` with unmatched names dropped and reported (`:68-83`), or
integerish indices (`:84-95`). It returns `NULL` for "nothing to do" (`:97-102`),
after reporting it itself, which becomes an early return of the unchanged input
(`R/knn_imp.R:135-147`).

Both early returns hand back the input through `new_slideimp_results()`, not as
a bare matrix, so the class and attributes of the result do not depend on which
path the data takes; `tests/testthat/test-knn_imp.R:502-536` pins that, and pins
that the empty-subset path emits exactly one message. `pca_imp()` still returns
a bare matrix on its own no-missing-values path (`R/pca_imp.R:460-463`); see
`dev/to-do.md`.

Missingness and partitioning. `cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)`
at `R/knn_imp.R:150`; `mat_miss()` (`R/mat_miss.R:33`) calls `col_miss_internal()`
(`src/mat_stats.cpp:242`), which counts `std::isnan` per column and does not
count `Inf` as missing (`src/mat_stats.cpp:257-260`). A second early return fires
if no subset column has any missing value (`R/knn_imp.R:152-166`). `colmax` then
becomes an eligibility mask, `eligible <- miss_rate <= min(colmax, 1)`
(`R/knn_imp.R:172`); the test is inclusive, matching `pca_imp()`
(`R/pca_imp.R:469`), so a column strictly above the threshold is dropped from
the K-NN stage entirely -- neither imputed nor usable as a neighbor. Three
disjoint index groups follow at `R/knn_imp.R:185-192`: `grp_impute` (eligible, has missing
values, in `subset` -- these get imputed), `grp_miss_no_imp` (eligible, has
missing values, not in `subset` -- neighbor candidate only), and `grp_complete`
(eligible, fully observed -- neighbor candidate only). `method` is remapped to
the integer the kernel expects, `0` = euclidean, `1` = manhattan
(`R/knn_imp.R:204`). So `subset` narrows only the *target* set: every eligible
column stays a neighbor candidate regardless of `subset`.

`rowmax` is not an argument of `knn_imp()`: it belongs to the tuning path
(`R/tune_imp.R:55`, `:184`), capping injected missingness per row.

After the `.Call`, `R/knn_imp.R:220` maps `NaN` back to `NA_real_` (the kernel
emits `NaN` for cells it could not fill), `R/knn_imp.R:222-223` uses the returned
`(row, col)` pairs as a two-column index matrix to scatter values into `obj`,
`R/knn_imp.R:232-240` optionally runs column-mean fill, and
`R/knn_imp.R:242-250` attaches the `slideimp_results` class and attributes via
`new_slideimp_results()` (`R/utils.R:5-25`).

**Which columns post-imputation touches.** Only two kinds of subset column can
still be missing after the kernel: one `colmax` held out of `grp_impute`, and
one the kernel could not finish because no candidate shared an observed row
with it. Both sets are already in hand at `R/knn_imp.R:233-236` -- the first
from `has_miss_idx` minus `grp_impute`, the second from the `NA`s in the third
result column -- so `mean_imp_col()` is called with just those columns, and not
at all when there are none. That is the usual case, and skipping it avoids a
second full-matrix allocation plus a per-column copy of everything K-NN already
filled. Both branches are pinned at `tests/testthat/test-knn_imp.R:473-500`.

## Layer 2 -- the .Call boundary

`R/RcppExports.R:16-18` defines the R-side shim `impute_knn_brute(obj, k,
grp_impute, grp_miss_no_imp, grp_complete, method, dist_pow, cores, pb)`, whose
body is a single `.Call` to `_slideimp_impute_knn_brute`.

The generated C wrapper is `src/RcppExports.cpp:69-85`; it converts each SEXP
through `Rcpp::traits::input_parameter<...>` (`:73-81`) and calls the C++ symbol
at `:82`. Registration is `:195` inside the `CallEntries` table, installed by
`R_init_slideimp` at `:208`. The exact C++ entry symbol is `impute_knn_brute`,
defined at `src/impute_knn_brute.cpp:397` under the `// [[Rcpp::export]]` marker
on `:396`. It is the only C++ entry point the K-NN stage needs; `check_finite`,
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
  `src/imputed_value.cpp:1`. Defines `mask_t` / `MaskMat` (`:12-13`),
  `GroupLayout` (`:25-34`), the header-inline `validate_knn_inputs()`
  (`:36-160`), and the two out-of-line declarations (`:166`, `:178`).
- `src/matrix_checks.h` -- pulled in by `src/imputed_value.h:9` and
  `src/mat_stats.cpp:4`. Holds only `stop_on_inf()` (`src/matrix_checks.h:8-29`).
- `src/par_for.h` -- included by `src/impute_knn_brute.cpp:6` and
  `src/mat_stats.cpp:5`. One function template: `par_for()` runs the loop body
  inline when `n_threads <= 1` and hands off to `RcppThread::parallelFor()`
  otherwise. The reason is in the header's own comment -- RcppThread's free
  `parallelFor()` sets the global pool's thread count and restores it, and
  quickpool implements a change of count by joining every worker and spawning
  new ones, so a one-thread call rebuilds the pool twice for nothing. That is
  the default configuration of the whole chain: `slide_imp()` runs each window
  at one core and `group_imp()` forces one core per worker under mirai. The
  serial and parallel paths are pinned equal at
  `tests/testthat/test-knn_imp.R:410-431`.
- `src/loc_timer.h` -- included at `src/impute_knn_brute.cpp:2`. Unless `LOC_TIMER`
  is defined, the statement macros compile to `((void)0)` and `LOC_TIMER_PARAM` /
  `LOC_TIMER_ARG` to nothing (`:26-34`), so the
  `LOC_*` calls at `src/impute_knn_brute.cpp:479-480` and `519-520` vanish in a
  normal build.

Plus the external `RcppArmadillo.h` (`src/imputed_value.h:4`) and `RcppThread.h`
(`src/impute_knn_brute.cpp:5`). The remaining `src/` headers (`eig_sym_sel.h`,
`gram_ops.h`, `hybrid_topk_eig.h`, `lobpcg_warm.h`, `pca_linalg_utils.h`,
`svd_triplet.h`) and `.cpp` files (`armaSVD.cpp`, `find_windows.cpp`,
`find_windows_flank.cpp`, `sample_each_rep_cpp.cpp`) serve the PCA, windowing,
and simulation paths and are not reached from `knn_imp()`.

`src/Makevars:1` sets `-DARMA_64BIT_WORD=1`, so `arma::uword` is 64-bit and
kernel index arithmetic does not wrap on matrices past `2^32 - 1` elements;
`src/mat_stats.cpp:181` exposes `arma_uword_bytes()` so the test suite can assert
that flag survived.

## Layer 4 -- kernel setup and result skeleton

1. `validate_knn_inputs()` (`src/imputed_value.h:36-160`) re-checks most of what R
   already checked, so the kernel is safe when called directly. All-`NA` columns are
   left to the zero-`n_valid` guard (`src/impute_knn_brute.cpp:215-218`), `cores` is
   floored rather than validated (`:471`), and the `k` bound applies only when a column
   needs imputing. It checks shape (`:45-53`),
   `k >= 1` (`:55-58`), `method` in `{0, 1}` (`:60-63`), `dist_pow` finite and
   non-negative (`:65-68`), every group index in range (`:70-89`), the three
   groups disjoint (`:91-115`), every `grp_complete` column free of `NaN`
   (`:117-138`), and `k <= n_working - 1` for `n_working` = sum of the three
   group sizes (`:140-159`). The last two are not redundant with the R layer:
   group 3 is read unmasked by `impute_column_values()`, so a `NaN` there would
   flow into the weighted sum and produce a missing value with no error, and a
   column claimed by two groups would be its own neighbor at distance zero.
   `tests/testthat/test-knn_imp.R:433-471` exercises both. Then
   `stop_on_inf(obj)` at `src/impute_knn_brute.cpp:410` rescans for infinities.
2. `GroupLayout layout{...}` at `src/impute_knn_brute.cpp:411` is the single
   source of truth for group boundaries: groups 1 and 2 occupy local columns
   `[0, n_masked())` of the working matrices, and `complete_start()`
   (`src/imputed_value.h:32`) is a *virtual* index used as a tag -- a neighbor id
   `>= complete_start()` means "group 3, offset `id - complete_start()` into
   `grp_complete`".
3. Working buffers `obj_masked` (`arma::mat`), `nmiss_masked`
   (`arma::Mat<uint8_t>`) and `n_col_valid` are allocated at
   `src/impute_knn_brute.cpp:415-417` and filled by `copy_with_mask`
   (`:419-437`), which calls `std::isnan` once per cell and reuses the result as
   the mask byte, the count increment, and the zeroing test (`:431-434`). Missing
   cells therefore hold a literal `0.0` in `obj_masked` and are excluded by the
   mask, never by a branch in the distance loop.
4. Group 1 is copied at `src/impute_knn_brute.cpp:440-443`, group 2 at `:457-460`.
   Group 3 is **not** copied: the kernel reads it straight out of `obj` through
   `grp_complete` (`:462`, and the lambdas at `:239-253`).
5. `initialize_result_matrix()` (`src/imputed_value.cpp:18-165`) runs between the
   group-1 and group-2 copies (`src/impute_knn_brute.cpp:448-449`) so the call
   can bail early when there is nothing to impute (`:451-454`).

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

## Layer 5 -- neighbor search (`src/impute_knn_brute.cpp:46-333`)

**Metrics.** `EuclideanMetric` (`:20-23`) accumulates `diff * diff`,
`ManhattanMetric` (`:25-28`) accumulates `std::abs(diff)`. The metric is a
template parameter and `distance_vector()` (`:310`) dispatches once per target
column on `method` (`:321-332`), so the inner body is monomorphized.

**Distance over incomplete data.** `calc_distance_raw<Metric, Bound>` (`:46-92`)
handles a pair where either side may be missing. Per row it forms
`valid = target_nmiss[r] & other_nmiss[r]` (`:64`) and adds
`valid * Metric::accumulate(diff)` (`:66`) -- branchless, and correct because
missing values were zeroed during the masked copy. It returns `dist / n_valid`
(`:91`), the mean per-coordinate contribution over rows observed in *both*
columns, or `arma::datum::inf` when the two columns share no observed row
(`:86-89`). No square root is taken for the euclidean case, so the stored
quantity is a mean squared difference.
`calc_distance_raw_complete<Metric, Bound>` (`:102-145`) is the group-3 variant:
the other side is fully observed, so only the target mask matters (`:126`) and
`n_valid`, the target's own observed count, is known up front and passed in
(`:107`). It has no zero-`n_valid` branch of its own -- see the guard below.

**Pruning.** Both kernels chunk rows in blocks of `GRAIN = 16` (`:31`) and, when
`Bound` is true, test a partial-sum bound at each chunk boundary, returning `inf`
on a prune (`:70-76`, `:129-135`; derivations in the comments at `:36-43`,
`:97-99`). The check sits under `if constexpr (Bound)`, so it is not compiled
into the `Bound = false` instantiations at all; the `worst_dist = +inf` default
(`:52`, `:108`) only covers a `Bound = true` caller with no worst distance yet.
The pruned and unpruned paths are compared against each other implicitly at
`tests/testthat/test-knn_imp.R:237-325`, where the pure-R reference computes
every distance unpruned and the kernel prunes past the first `k`.

**Top-k selection.** `distance_vector_impl` (`:192-307`) keeps a
`std::vector<NeighborInfo>` of capacity `k` (`:220-221`) and runs two phases with
three shared cursors `p1`, `p2`, `c` (`:255-257`) so every candidate column is
visited exactly once. Before any of it, `:215-218` returns an empty vector when
the target column has no observed row: it then overlaps nothing, and
`calc_distance_raw_complete` would divide by a zero `n_valid`. `check_finite()`
keeps that case away from the R entry point, so the guard exists for direct
callers of the kernel; `tests/testthat/test-knn_imp.R:464-470` is one.

- *Fill* (`:260-271`): append the first `k` candidates with no pruning via
  `insert_before_k` (`:158-161`), in the order group 3, masked columns below the
  target, masked columns above it. `p2` starts at `index + 1` (`:256`) and the
  `p1` loop stops at `index`, so the target column is never its own neighbor.
  Then one `std::stable_sort` on those `k` entries (`:273-279`); why stable is
  under Ties below.
- *Replacement* (`:282-296`): resume each cursor, compute the distance with
  `Bound = true` against the current worst (`top_k.back().distance`), and hand it
  to `insert_if_better_than_worst` (`:164-177`), which rejects anything not
  strictly better than the worst (`:166-169`), overwrites the last slot, and
  bubbles it into position (`:173-176`).
- *Tail drop* (`:298-304`): the result is sorted ascending, so candidates at
  distance `+inf` (no observed row in common with the target) sit at the tail
  and are popped before returning. The caller can therefore receive fewer than
  `k` neighbors, or none at all; when none, the column's result rows keep their
  initial `NaN` (`:496-501`).

The structure is therefore a sorted fixed-size array with insertion, not a heap
and not a full sort of all candidates: exactly one `std::stable_sort`, over `k`
elements.

**Ties.** Among candidates at equal distance, the one met first in scan order
wins. Scan order is the fill order above: group 3 in `grp_complete` order, then
the masked columns in local order (`grp_impute`, then `grp_miss_no_imp`),
skipping the target. Two things make it hold. The fill appends in scan order,
and `std::stable_sort` (`:277-279`) keeps equal distances in that order. After
that both comparisons are strict `<` (`:173` for the bubble, and the
`dist >= top_k.back().distance` early return at `:166`), so a candidate that
ties an incumbent never displaces it and is placed behind it. `top_k` is
therefore always ordered by distance and then scan position, the entry evicted
from `top_k.back()` is the last-met of the tied worst, and the result is the
first `k` of all candidates under that order. That final order is also the
order `impute_column_values()` accumulates in, so it is fixed too. None of it
depends on thread count.

It is not "leftmost column wins": a complete column beats a tied column with
missing values to its left. `knn_ref()` in `tests/testthat/test-knn_imp.R`
breaks ties by column order, and agrees with the kernel only because its
`rnorm` data has no ties; on the small-integer matrices of
`dev/probe_knn_ties.R` the two rules give different values in 2-29% of imputed
cells.

The sort used to be `std::sort`, which leaves equal elements in an unspecified
order, so the standard library chose which tied fill entry was evicted.
libstdc++ insertion-sorts ranges of up to 16, which is stable, so at
`k <= 16` the old build already followed the rule (bitwise at `k` = 5, 10 and
16 in the probe) and its results there are unchanged. From `k = 17` it evicted
tied candidates from the front of the run instead of the back, and changed
4-25% of imputed cells on those matrices (GCC 14.2). libc++, which the macOS
builds use, was not measured. `std::stable_sort` may take a temporary buffer
per target column; an A/B in `dev/probe_knn_ties.R` (Part 3) found no cost
above run-to-run noise. Pinned at `tests/testthat/test-knn_imp.R:327-408`,
which reads the chosen neighbor set back out of the imputed value and fails
against the `std::sort` build at `k = 17` and `k = 32`.

## Layer 6 -- weights and imputation

**Weights.** `knn_weights()` (`src/impute_knn_brute.cpp:336-391`) returns a
vector of ones when `dist_pow == 0` (`:342-345`) -- the plain unweighted mean of
neighbors. Otherwise it finds the smallest strictly positive finite distance
`min_pos` and notes whether any distance is zero (`:347-361`). Infinite
distances never reach this function (Layer 5 drops them), so `min_pos` is
non-finite only when every neighbor is an exact duplicate, and weights then
stay equal (`:365-368`). Two live branches: when some distance is exactly zero
(`:370-381`), `floor = min_pos * sqrt(DBL_EPSILON)` and
`w_j = (floor / max(d_j, floor)) ^ dist_pow`, so duplicate columns get
weight 1 and everything else a tiny positive weight; otherwise (`:382-388`)
`w_j = (min_pos / d_j) ^ dist_pow`. Either way the largest weight is 1.

The `d_j` here is the Layer 5 quantity, so under `method = "euclidean"` it is a
mean SQUARED difference and `dist_pow` is an exponent of `2 * dist_pow` on the
square-rooted euclidean distance. The `min_pos` factor is common to every weight
and cancels out of the weighted average; it only keeps the weights near 1. Both
facts are stated on the help page (`R/knn_imp.R:47-74`) and pinned against a
pure-R reference at `tests/testthat/test-knn_imp.R:237-325`, which reproduces
the kernel to 1e-12 on max absolute and max relative error over both metrics,
three `dist_pow` values, and two candidate layouts.

**Accumulation.** `impute_column_values()` (`src/imputed_value.cpp:178-269`)
splits the neighbor ids into a masked bucket and a complete bucket by comparing
against `complete_start()` (`:209-219`). Masked neighbors contribute
`w * nmiss[row]` to both the weighted sum and the weight total (`:222-237`), so a
neighbor that is itself missing in a given row drops out for that row only.
Complete neighbors add `w * value` per row and a single shared
`total_complete_w` to every row's weight total (`:243-262`). The final value is
`weighted_sum / weight_total`, or `NaN` when the total is zero (`:264-268`),
which R converts to `NA` at `R/knn_imp.R:220`.

That per-row drop-out is user-visible and is stated on the help page
(`R/knn_imp.R:62-67`): neighbors are chosen once per column, but the set that
actually contributes can differ from row to row, and a row where none of the
chosen `k` is observed comes back missing. The "mixed candidates" half of
`tests/testthat/test-knn_imp.R:237-325` is the case that exercises it.

## Feasibility and abort paths

`slideimp_infeasible` is a cli condition class raised by the R layer so the
wrappers can catch a per-group or per-window failure and fall back instead of
dying. `knn_imp()` has exactly two sites:

1. `R/knn_imp.R:175-183` -- `k > n_elig - 1`, where `n_elig` is the number of
   columns passing `colmax` (`R/knn_imp.R:173`); message "`k` (...) exceeds
   usable columns (...)". Exercised at `tests/testthat/test-knn_imp.R:115-118`
   and `tests/testthat/test-group_imp.R:739-742`.
2. `R/knn_imp.R:194-202` -- `length(grp_impute) == 0`, i.e. every subset column
   with missing values also exceeds `colmax`; message "All subset columns with
   missing values exceed `colmax` (...)". Exercised at
   `tests/testthat/test-group_imp.R:744-752` and, for the inclusive boundary,
   `tests/testthat/test-knn_imp.R:121-137`.

Everything else that aborts here is a plain error, not `slideimp_infeasible`:
`check_finite()` for `Inf` and all-`NA` columns
(`tests/testthat/test-knn_imp.R:553-560` and `:112-113`), plus the C++ guards in
`validate_knn_inputs()` (`tests/testthat/test-knn_imp.R:433-471`) and
`initialize_result_matrix()`. Catch sites are in the
wrappers only (next section); each switches on `on_infeasible` to abort (`slide_imp()`
rethrows, `group_imp()` raises a new `cli_abort()` with the original as `parent`), skip
the block unchanged, or fall back to `mean_imp_col()`.

## Data structures and ownership

- The R matrix is **not** copied at the boundary. `input_parameter<const
  arma::mat&>` (`src/RcppExports.cpp:73`) resolves to
  `ConstReferenceInputParameter`, which for a double matrix is
  `ArmaMat_InputParameter<..., false_type>`
  (`RcppArmadillo/interface/RcppArmadilloAs.h:573-585`): it constructs
  `arma::mat(ptr, nrow, ncol, false)` over an `Rcpp::NumericMatrix`, `copy_aux_mem =
  false`. The branch is chosen by the Armadillo element type, so a `double` matrix
  always takes it; `:587-599` is for element types that need a cast. For double-storage
  input the `NumericMatrix` wraps R's own memory, so the kernel's `const arma::mat& obj`
  aliases the R matrix and is never written to. Integer-storage input is coerced to a
  fresh double copy when the `NumericMatrix` is built.
- The kernel returns a freshly allocated `(n_missing x 3)` `arma::mat` of
  triplets, not an imputed matrix. The only in-place mutation of the user's
  matrix happens in R, at `R/knn_imp.R:223`.
- So the real full-matrix copies on this path are both on the R side: R's
  copy-on-modify at the scatter (`R/knn_imp.R:223`), and, when
  `mean_imp_col()` actually has work to do, the second `n x m` allocation in
  `mean_imp_col_internal()` (`src/mat_stats.cpp:120`). The second one is
  skipped whenever nothing is still missing.
- Groups 1 and 2 are materialized a second time into `obj_masked` /
  `nmiss_masked` (`src/impute_knn_brute.cpp:415-416`) with NaNs zeroed and a
  parallel byte mask; `nmiss_masked` is `arma::Mat<uint8_t>`
  (`src/imputed_value.h:12-13`), one byte per cell rather than a `umat`. Group 3
  is read in place from `obj`, so a matrix dominated by complete columns pays no
  second copy.
- Column indices cross the boundary 0-based (`R/knn_imp.R:210-212` subtracts 1)
  and come back 1-based (`src/imputed_value.cpp:152` and `:159` add 1).
- Per-thread scratch is small and local: `top_k`
  (`src/impute_knn_brute.cpp:220`), whatever temporary buffer `std::stable_sort`
  takes for it (`:277`), `nn_columns` and `weights` (`:502-507`), and the two
  accumulators in `impute_column_values` (`src/imputed_value.cpp:196-197`).
  `mean_imp_col_internal()` (`src/mat_stats.cpp:79`) is likewise out-of-place: a
  fresh matrix with untouched columns `memcpy`d (`src/mat_stats.cpp:137-141`),
  and `R/mean_imp_col.R:54` reattaches dimnames.

## Parallelism

`par_for()` at `src/impute_knn_brute.cpp:482-517` runs over
`[0, layout.n_imp)` -- one task per target column. `n_threads` and `n_batches`
are both set to `cores`, floored at 1 (`src/impute_knn_brute.cpp:471-473`, passed
at `:517`), so the index range splits into as many contiguous batches as threads.
At one thread `par_for()` (`src/par_for.h`) runs the body inline instead, for
the reason given in Layer 3.

All threads read `obj`, `obj_masked`, `nmiss_masked`, and `grp_complete` without
writing them, and each writes only rows `[col_offsets(i), col_offsets(i+1))` of
`result`, a disjoint slice per target column, so no locks are taken. Interruption
is checked with `RcppThread::checkUserInterrupt()` on every fifth index
(`src/impute_knn_brute.cpp:487-490`), which RcppThread routes to the main thread
to unwind the pool. The optional progress bar is an `RcppThread::ProgressBar` in
a `unique_ptr` created only when `pb` is true (`:474-478`) and incremented from
the workers at `:512-515`.

`src/Makevars:2-3` adds the OpenMP flags and links `RcppThread::LdFlags()`;
`src/Makevars.win:2-3` omits the RcppThread link line. `col_vars_internal()`
(`src/mat_stats.cpp:30-76`) and `mean_imp_col_internal()`
(`src/mat_stats.cpp:79-174`) go through the same `par_for()`.

## How the wrappers reach knn_imp()

`group_imp()` picks the function object at `R/group_imp.R:811`
(`imp_fn <- if (is_knn_mode) knn_imp else pca_imp`), pre-builds a per-group
parameter list injecting `cores` and a group-local `subset` (`:803-806`), and
invokes it via `do.call()` inside a `tryCatch` -- `:859-882` in the mirai branch,
`:929-949` in the sequential branch. Each group is a column slice
`obj[, indices[[i]]$col_idx]` (`:857`, `:927`), so `knn_imp()` sees a submatrix
and group-local indices.

`slide_imp()` calls `knn_imp()` directly at `R/slide_imp.R:460-471` for each
window `obj[, start[i]:end[i]]` (`R/slide_imp.R:454-455`), passing
`na_check = FALSE`, `.progress = FALSE`, and a window-local `subset`
(`R/slide_imp.R:470`), wrapped in the same `slideimp_infeasible` handler
(`R/slide_imp.R:495-505`). Both wrappers set `na_check = FALSE`
(`R/group_imp.R:800`) and do the final NA accounting themselves.

## Test entry points

`tests/testthat/test-knn_imp.R:1-45` calls `impute_knn_brute()` directly and
asserts the returned `(row, col)` triplet locations match `which(is.na(...))`;
`:67-102` covers `subset` by index and by name with `post_imp` both ways;
`:104-119` covers all-NA rows and columns plus abort site 1; `:327-408` pins
the tie rule; `:553-560` covers `Inf` rejection;
`tests/testthat/test-group_imp.R:738-753` pins the `slideimp_infeasible` class
on both abort sites.
