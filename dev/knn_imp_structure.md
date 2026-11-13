# knn_imp() structural map

A structural map of the K-NN imputation path, from the exported R function
down to the C++ kernel. Descriptive only: what the code is and how it connects.
Code is cited by file and by a fragment of one line that occurs exactly once in
that file, never by line number, e.g. (`R/knn_imp.R`: `knn_imp <- function(`).
`uv run dev/check_map_anchors.py` checks every citation; its header gives the
forms.

Orientation: `knn_imp()` documents rows as samples and columns as features, and
the neighbor search runs **between columns** -- a "neighbor" is another column of
`obj`, and the distance loop walks rows. Hence `k` is capped at `ncol(obj) - 1`
(`R/knn_imp.R`: `upper = ncol(obj) - 1`), and `dev/check_impute_knn.R` transposes
both sides when comparing against `impute::impute.knn`, which searches rows. That
comparison is real and asserted -- it agrees to 2.2e-16 absolute and 4.5e-16
relative, not bitwise. It lives in `dev/` rather than under `tests/` because
`impute` is a Bioconductor package outside `Suggests`, and any `impute::` call
under `tests/` makes `R CMD check` report an unstated dependency regardless of
skips; the reasoning is repeated in a comment (`tests/testthat/test-knn_imp.R`:
`The comparison against impute::impute.knn`) so a future session does not put it
back.

## Call chain overview

```
R/knn_imp.R                         knn_imp()  -- user entry, validation, partition
  +- src/mat_stats.cpp              check_finite()  (Inf + all-NA column gate)
  +- R/utils.R                      resolve_subset()  -- names/ints -> col indices
  +- R/mat_miss.R                   mat_miss()
  |    +- src/mat_stats.cpp         col_miss_internal() via R/RcppExports.R
  +- [abort site 1] R/knn_imp.R     k > n_elig - 1
  +- [abort site 2] R/knn_imp.R     all subset cols exceed colmax
  +- R/RcppExports.R                impute_knn_brute()  -- .Call shim
  |    +- src/RcppExports.cpp       _slideimp_impute_knn_brute (SEXP wrapper)
  |         +- src/impute_knn_brute.cpp  impute_knn_brute()   <-- C++ entry
  |              +- src/imputed_value.h       validate_knn_inputs()
  |              +- src/matrix_checks.h       stop_on_inf()
  |              +- src/impute_knn_brute.cpp  copy_with_mask lambda
  |              +- src/imputed_value.cpp     initialize_result_matrix()
  |              +- src/par_for.h             par_for()
  |                   +- src/impute_knn_brute.cpp  distance_vector()
  |                   |    +- src/impute_knn_brute.cpp  distance_vector_impl<Metric>()
  |                   |         +- src/impute_knn_brute.cpp  calc_distance_raw<Metric, Bound>()
  |                   |         +- src/impute_knn_brute.cpp  calc_distance_raw_complete<Metric, Bound>()
  |                   |         +- src/impute_knn_brute.cpp  insert_before_k()
  |                   |         +- src/impute_knn_brute.cpp  insert_if_better_than_worst()
  |                   +- src/impute_knn_brute.cpp  knn_weights()
  |                   +- src/imputed_value.cpp     impute_column_values()
  +- R/knn_imp.R                    NaN -> NA, scatter triplets back into obj
  +- R/mean_imp_col.R               mean_imp_col()  (post_imp, and only if
  |                                 anything is still missing)
  |    +- src/mat_stats.cpp         check_inf()
  |    +- src/mat_stats.cpp         mean_imp_col_internal()
  +- R/utils.R                      new_slideimp_results()
```

## Layer 1 -- R entry point (`R/knn_imp.R`)

Signature (`R/knn_imp.R`: `knn_imp <- function(`). Validation block, in order
(`R/knn_imp.R`: `checkmate::assert_matrix(`): numeric matrix with `>= 1` row and
`>= 2` columns; `check_finite(obj)`; `method` matched against `"euclidean"` /
`"manhattan"`; `k` integer in `[1, ncol(obj) - 1]`; `cores >= 1`; `colmax` in
`[0, 1]`; `dist_pow` a finite number `>= 0`; the `post_imp` and `.progress`
flags. All but two go through `checkmate` with an explicit `.var.name`, so a
rejection names the user's argument. The exceptions are `method`, matched with
`match.arg()`, and `check_finite()`, which raises from C++ through `Rcpp::stop`
(`src/mat_stats.cpp`: `Infinite value found at row`). A test pins the
`checkmate` path for `dist_pow` (`tests/testthat/test-knn_imp.R`:
`dist_pow is validated by name`).

Missing-value gating is delegated to C++. `check_finite()`
(`src/mat_stats.cpp`: `void check_finite(`) makes one pass per column and
refuses twice: any `Inf`/`-Inf` aborts with the row and column of the offender
(`src/mat_stats.cpp`: `Infinite value found at row`), and a column ending the
pass with no finite cell aborts with "All NA/NaN column detected". `check_inf()`
(`src/mat_stats.cpp`: `void check_inf(`) is the weaker Inf-only variant, called
from `mean_imp_col()` (`R/mean_imp_col.R`: `check_inf(obj)`) but not from
`knn_imp()`. Both are declared in `R/RcppExports.R`. `NA` and `NaN` are
indistinguishable throughout; the policy is stated in the package help
(`R/slideimp-package.R`: `@section Missing values and non-finite input`).

`subset` goes through `resolve_subset()`
(`R/utils.R`: `resolve_subset <- function(`), which accepts `NULL` (all
columns), a character vector matched against `colnames(obj)` with unmatched
names dropped and reported, or integerish indices. It returns `NULL` for
"nothing to do", after reporting it itself, which becomes an early return of the
unchanged input (`R/knn_imp.R`: `if (is.null(subset)) {`).

Both early returns hand back the input through `new_slideimp_results()`, not as
a bare matrix, so the class and attributes of the result do not depend on which
path the data takes; a test pins that (`tests/testthat/test-knn_imp.R`:
`early returns carry the slideimp_results class`), and pins that the
empty-subset path emits exactly one message. `pca_imp()` does the same on its
no-missing-values path. The class costs one copy of the input: a matrix the
caller still holds cannot take an attribute in place (`R/utils.R`:
`Attaching anything to a matrix the caller still holds copies it`).

Missingness and partitioning. `cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)`;
`mat_miss()` (`R/mat_miss.R`: `mat_miss <- function(`) calls
`col_miss_internal()` (`src/mat_stats.cpp`: `arma::uvec col_miss_internal(`),
which counts `std::isnan` per column and does not count `Inf` as missing
(`src/mat_stats.cpp`: `count += std::isnan(p[i]);`). A second early return fires
if no subset column has any missing value (`R/knn_imp.R`:
`if (!any(cmiss[subset] > 0L)) {`). `colmax` then becomes an eligibility mask,
`eligible <- miss_rate <= min(colmax, 1)`; the test is inclusive, matching
`pca_imp()` (`R/pca_imp.R`: `eligible <- miss_rate <= min(colmax, 1)`), so a
column strictly above the threshold is dropped from the K-NN stage entirely --
neither imputed nor usable as a neighbor. Three disjoint index groups follow:
`grp_impute` (eligible, has missing values, in `subset` -- these get imputed),
`grp_miss_no_imp` (eligible, has missing values, not in `subset` -- neighbor
candidate only), and `grp_complete` (eligible, fully observed -- neighbor
candidate only). `method` is remapped to the integer the kernel expects, `0` =
euclidean, `1` = manhattan (`R/knn_imp.R`: `method <- switch(method`). So
`subset` narrows only the *target* set: every eligible column stays a neighbor
candidate regardless of `subset`.

`rowmax` is not an argument of `knn_imp()`: it belongs to the tuning path
(`R/tune_imp.R`: `sample_na_loc <- function(`, `max_row_miss <- max_miss_count(`),
capping injected missingness per row.

After the `.Call`, `knn_imp()` maps `NaN` back to `NA_real_`, since the kernel
emits `NaN` for cells it could not fill (`R/knn_imp.R`:
`imputed_values[is.nan(imputed_values)]`); uses the returned `(row, col)` pairs
as a two-column index matrix to scatter values into `obj` (`R/knn_imp.R`:
`obj[imp_indices] <- imputed_values[, 3]`); optionally runs column-mean fill
(`R/knn_imp.R`: `if (post_imp) {`); and attaches the `slideimp_results` class
and attributes via `new_slideimp_results()` (`R/utils.R`:
`new_slideimp_results <- function(`). The constructor first strips its input to
a bare matrix, so a `slideimp_results` passed back in for a second pass does not carry
its class or its old attributes into the new result.

**Which columns post-imputation touches.** Only two kinds of subset column can
still be missing after the kernel: one `colmax` held out of `grp_impute`, and
one the kernel could not finish because no candidate shared an observed row
with it. Both sets are already in hand (`R/knn_imp.R`:
`na_cols <- sort(unique(c(`) -- the first from `has_miss_idx` minus
`grp_impute`, the second from the `NA`s in the third result column -- so
`mean_imp_col()` is called with just those columns, and not at all when there
are none. That is the usual case, and skipping it avoids a second full-matrix
allocation plus a per-column copy of everything K-NN already filled. Both
branches are pinned by a test (`tests/testthat/test-knn_imp.R`:
`post_imp fills columns K-NN could not finish`).

## Layer 2 -- the .Call boundary

`R/RcppExports.R` defines the R-side shim `impute_knn_brute(obj, k,
grp_impute, grp_miss_no_imp, grp_complete, method, dist_pow, cores, pb)`, whose
body is a single `.Call` to `_slideimp_impute_knn_brute`.

The generated C wrapper is in `src/RcppExports.cpp`; it converts each SEXP
through `Rcpp::traits::input_parameter<...>` and calls the C++ symbol.
Registration is its entry in the `CallEntries` table, installed by
`R_init_slideimp`. The exact C++ entry symbol is `impute_knn_brute`, defined in
`src/impute_knn_brute.cpp` under the `// [[Rcpp::export]]` marker. It is the
only C++ entry point the K-NN stage needs; `check_finite`, `col_miss_internal`,
and `mean_imp_col_internal` are separate `.Call`s made from R around it.

## Layer 3 -- C++ translation units and headers on the path

Four compiled `.cpp` files are reached: `src/RcppExports.cpp` (generated SEXP
wrappers), `src/impute_knn_brute.cpp` (K-NN kernel and entry point),
`src/imputed_value.cpp` (result layout and accumulation loop), and
`src/mat_stats.cpp` (`check_finite`, `check_inf`, `col_miss_internal`,
`mean_imp_col_internal`).

Project headers reached, following the `#include` graph:

- `src/imputed_value.h` -- included by `src/impute_knn_brute.cpp` and
  `src/imputed_value.cpp`. Defines `mask_t` / `MaskMat`, `GroupLayout`, the
  header-inline `validate_knn_inputs()`, and the two out-of-line declarations
  (`src/imputed_value.h`: `arma::mat initialize_result_matrix(`,
  `void impute_column_values(`).
- `src/matrix_checks.h` -- pulled in by `src/imputed_value.h` and
  `src/mat_stats.cpp`. Holds only `stop_on_inf()`.
- `src/par_for.h` -- included by `src/impute_knn_brute.cpp` and
  `src/mat_stats.cpp`. One function template: `par_for()` runs the loop body
  inline when `n_threads <= 1` and hands off to `RcppThread::parallelFor()`
  otherwise. The reason is in the header's own comment -- RcppThread's free
  `parallelFor()` sets the global pool's thread count and restores it, and
  quickpool implements a change of count by joining every worker and spawning
  new ones, so a one-thread call rebuilds the pool twice for nothing. That is
  the default configuration of the whole chain: `slide_imp()` runs each window
  at one core and `group_imp()` forces one core per worker under mirai. The
  serial and parallel paths are pinned equal by a test
  (`tests/testthat/test-knn_imp.R`: `the parallel path gives the same answer`).
- `src/loc_timer.h` -- included by `src/impute_knn_brute.cpp`. Unless
  `LOC_TIMER` is defined, the statement macros compile to `((void)0)` and
  `LOC_TIMER_PARAM` / `LOC_TIMER_ARG` to nothing, so the `LOC_*` calls in
  `src/impute_knn_brute.cpp` vanish in a normal build.

Plus the external `RcppArmadillo.h` (included by `src/imputed_value.h`) and
`RcppThread.h` (by `src/impute_knn_brute.cpp`). The remaining `src/` headers
(`eig_sym_sel.h`, `gram_ops.h`, `hybrid_topk_eig.h`, `lobpcg_warm.h`,
`pca_linalg_utils.h`, `svd_triplet.h`) and `.cpp` files (`armaSVD.cpp`,
`find_windows.cpp`, `find_windows_flank.cpp`, `sample_each_rep_cpp.cpp`) serve
the PCA, windowing, and simulation paths and are not reached from `knn_imp()`.

`src/Makevars` sets `-DARMA_64BIT_WORD=1`, so `arma::uword` is 64-bit and
kernel index arithmetic does not wrap on matrices past `2^32 - 1` elements;
`src/mat_stats.cpp` exposes `arma_uword_bytes()` so the test suite can assert
that flag survived.

## Layer 4 -- kernel setup and result skeleton

1. `validate_knn_inputs()` (`src/imputed_value.h`: `void validate_knn_inputs(`)
   re-checks most of what R already checked, so the kernel is safe when called
   directly. All-`NA` columns are left to the zero-`n_valid` guard
   (`src/impute_knn_brute.cpp`: `if (target_n_valid == 0.0)`), `cores` is
   floored rather than validated (`src/impute_knn_brute.cpp`:
   `cores = std::max(1, cores);`), and the `k` bound applies only when a column
   needs imputing. It checks shape, `k >= 1`, `method` in `{0, 1}`, `dist_pow`
   finite and non-negative, every group index in range, the three groups
   disjoint, every `grp_complete` column free of `NaN`, and `k <= n_working - 1`
   for `n_working` = sum of the three group sizes. The partition and `NaN`
   checks are not redundant with the R layer: group 3 is read unmasked by
   `impute_column_values()`, so a `NaN` there would flow into the weighted sum
   and produce a missing value with no error, and a column claimed by two groups
   would be its own neighbor at distance zero. A test exercises both
   (`tests/testthat/test-knn_imp.R`: `the kernel checks the invariants`). Then
   `stop_on_inf(obj)` in `src/impute_knn_brute.cpp` rescans for infinities.
2. `GroupLayout layout{...}` in `src/impute_knn_brute.cpp` is the single
   source of truth for group boundaries: groups 1 and 2 occupy local columns
   `[0, n_masked())` of the working matrices, and `complete_start()`
   (`src/imputed_value.h`: `arma::uword complete_start() const`) is a *virtual*
   index used as a tag -- a neighbor id `>= complete_start()` means "group 3,
   offset `id - complete_start()` into `grp_complete`".
3. Working buffers `obj_masked` (`arma::mat`), `nmiss_masked`
   (`arma::Mat<uint8_t>`) and `n_col_valid` are allocated
   (`src/impute_knn_brute.cpp`: `arma::mat obj_masked(`) and filled by
   `copy_with_mask`, which calls `std::isnan` once per cell and reuses the
   result as the mask byte, the count increment, and the zeroing test
   (`src/impute_knn_brute.cpp`: `mask_t finite = !std::isnan(v);`). Missing
   cells therefore hold a literal `0.0` in `obj_masked` and are excluded by the
   mask, never by a branch in the distance loop.
4. Groups 1 and 2 are copied (`src/impute_knn_brute.cpp`:
   `copy_with_mask(i, grp_impute(i))`, `copy_with_mask(layout.mni_start() + p`).
   Group 3 is **not** copied: the kernel reads it straight out of `obj` through
   `grp_complete` (`src/impute_knn_brute.cpp`: `group 3: read directly from obj`),
   and in the lambdas (`src/impute_knn_brute.cpp`: `auto complete_dist_unbounded`).
5. `initialize_result_matrix()` (`src/imputed_value.cpp`:
   `arma::mat initialize_result_matrix(`) runs between the group-1 and group-2
   copies so the call can bail early when there is nothing to impute
   (`src/impute_knn_brute.cpp`: `if (result.n_rows == 0)`).

**Result skeleton.** The return value is a triplet matrix: 1-based row, original
1-based column, value. Missing counts per target column come from `n_col_valid`
rather than a second pass over the mask (`src/imputed_value.cpp`:
`n_miss_i = mask_n_rows - n_col_valid(i)`), and a prefix sum builds
`col_offsets` (`src/imputed_value.cpp`: `col_offsets.set_size(n_imp + 1)`) so
each target column owns a contiguous slice of result rows and threads never
collide. The matrix is filled with `NaN` (`src/imputed_value.cpp`:
`result.fill(arma::datum::nan)`), columns 0 and 1 are populated up front
(`src/imputed_value.cpp`: `orig_col_1based = grp_impute(i) + 1`), and
`rows_to_impute_vec[i]` records which rows of target column `i` are missing.
Empty cases short-circuit (`src/imputed_value.cpp`: `if (n_imp == 0)`,
`if (sum_missing == 0)`), and `Rcpp::stop` consistency guards surround the
bookkeeping.

## Layer 5 -- neighbor search (`src/impute_knn_brute.cpp`)

**Metrics.** `EuclideanMetric` accumulates `diff * diff`, `ManhattanMetric`
accumulates `std::abs(diff)`. The metric is a template parameter and
`distance_vector()` dispatches once per target column on `method`, so the inner
body is monomorphized.

**Distance over incomplete data.** `calc_distance_raw<Metric, Bound>` handles a
pair where either side may be missing. Per row it forms
`valid = target_nmiss[r] & other_nmiss[r]` and adds
`valid * Metric::accumulate(diff)` -- branchless, and correct because missing
values were zeroed during the masked copy. It returns `dist / n_valid`, the mean
per-coordinate contribution over rows observed in *both* columns, or
`arma::datum::inf` when the two columns share no observed row
(`src/impute_knn_brute.cpp`: `if (n_valid == 0.0)`). No square root is taken for
the euclidean case, so the stored quantity is a mean squared difference.
`calc_distance_raw_complete<Metric, Bound>` is the group-3 variant: the other
side is fully observed, so only the target mask matters
(`src/impute_knn_brute.cpp`: `target_nmiss[rr] * Metric::accumulate(diff)`) and
`n_valid`, the target's own observed count, is known up front and passed in. It
has no zero-`n_valid` branch of its own -- see the guard below.

**Pruning.** Both kernels chunk rows in blocks of `GRAIN = 16` and, when `Bound`
is true, test a partial-sum bound at each chunk boundary, returning `inf` on a
prune (`src/impute_knn_brute.cpp`: `if (dist > worst_dist * (n_valid`,
`if (dist > unnorm_bound)`), with derivations in the comments
(`src/impute_knn_brute.cpp`: `Bound math (conservative but correct)`,
`Bound math (exact)`). The check sits under `if constexpr (Bound)`, so it is not
compiled into the `Bound = false` instantiations at all; the `worst_dist = +inf`
default only covers a `Bound = true` caller with no worst distance yet. The
pruned and unpruned paths are compared against each other implicitly by a test
(`tests/testthat/test-knn_imp.R`: `the documented weighting rule reproduces`),
where the pure-R reference computes every distance unpruned and the kernel
prunes past the first `k`.

**Top-k selection.** `distance_vector_impl` keeps a `std::vector<NeighborInfo>`
of capacity `k` (`src/impute_knn_brute.cpp`: `std::vector<NeighborInfo> top_k;`)
and runs two phases with three shared cursors `p1`, `p2`, `c` so every candidate
column is visited exactly once. Before any of it, a guard
(`src/impute_knn_brute.cpp`: `if (target_n_valid == 0.0)`) returns an empty
vector when the target column has no observed row: it then overlaps nothing, and
`calc_distance_raw_complete` would divide by a zero `n_valid`. `check_finite()`
keeps that case away from the R entry point, so the guard exists for direct
callers of the kernel; one is a test (`tests/testthat/test-knn_imp.R`:
`all_na$obj[, 1] <- NA`).

- *Fill* (`src/impute_knn_brute.cpp`: `c < layout.n_complete && top_k.size() < k`):
  append the first `k` candidates with no pruning via `insert_before_k`, in the
  order group 3, masked columns below the target, masked columns above it. `p2`
  starts at `index + 1` and the `p1` loop stops at `index`, so the target column
  is never its own neighbor. Then one `std::stable_sort` on those `k` entries;
  why stable is under Ties below.
- *Replacement* (`src/impute_knn_brute.cpp`: `for (; c < layout.n_complete; ++c)`):
  resume each cursor, compute the distance with `Bound = true` against the
  current worst (`top_k.back().distance`), and hand it to
  `insert_if_better_than_worst`, which rejects anything not strictly better than
  the worst, overwrites the last slot, and bubbles it into position.
- *Tail drop* (`src/impute_knn_brute.cpp`: `while (!top_k.empty()`): the result
  is sorted ascending, so candidates at distance `+inf` (no observed row in
  common with the target) sit at the tail and are popped before returning. The
  caller can therefore receive fewer than `k` neighbors, or none at all; when
  none, the column's result rows keep their initial `NaN`
  (`src/impute_knn_brute.cpp`: `if (n_neighbors == 0)`).

The structure is therefore a sorted fixed-size array with insertion, not a heap
and not a full sort of all candidates: exactly one `std::stable_sort`, over `k`
elements.

**Ties.** Among candidates at equal distance, the one met first in scan order
wins. Scan order is the fill order above: group 3 in `grp_complete` order, then
the masked columns in local order (`grp_impute`, then `grp_miss_no_imp`),
skipping the target. Two things make it hold. The fill appends in scan order,
and `std::stable_sort` keeps equal distances in that order. After that both
comparisons are strict `<` -- the bubble (`src/impute_knn_brute.cpp`:
`top_k[i].distance < top_k[i - 1].distance`) and the
`dist >= top_k.back().distance` early return -- so a candidate that ties an
incumbent never displaces it and is placed behind it. `top_k` is therefore
always ordered by distance and then scan position, the entry evicted from
`top_k.back()` is the last-met of the tied worst, and the result is the first
`k` of all candidates under that order. That final order is also the order
`impute_column_values()` accumulates in, so it is fixed too. None of it depends
on thread count.

It is not "leftmost column wins": a complete column beats a tied column with
missing values to its left. `knn_ref()` in `tests/testthat/test-knn_imp.R`
breaks ties by column order, and agrees with the kernel only because its
`rnorm` data has no ties; on small-integer matrices, where exact ties are
common, the two rules gave different values in 2-29% of imputed cells. That
measurement was made by `dev/probe_knn_ties.R`, no longer tracked; it is
recoverable with `git show c8f6d3d:dev/probe_knn_ties.R`.

The sort used to be `std::sort`, which leaves equal elements in an unspecified
order, so the standard library chose which tied fill entry was evicted.
libstdc++ insertion-sorts ranges of up to 16, which is stable, so at
`k <= 16` the old build already followed the rule (bitwise at `k` = 5, 10 and
16 in the probe) and its results there are unchanged. From `k = 17` it evicted
tied candidates from the front of the run instead of the back, and changed
4-25% of imputed cells on those matrices (GCC 14.2). libc++, which the macOS
builds use, was not measured. `std::stable_sort` may take a temporary buffer
per target column; an A/B in the same probe (Part 3) found no cost
above run-to-run noise. Pinned by a test (`tests/testthat/test-knn_imp.R`:
`tied neighbors are settled by scan order`), which reads the chosen neighbor set
back out of the imputed value and fails against the `std::sort` build at
`k = 17` and `k = 32`.

## Layer 6 -- weights and imputation

**Weights.** `knn_weights()` (`src/impute_knn_brute.cpp`:
`static arma::vec knn_weights(`) returns a vector of ones when `dist_pow == 0`
-- the plain unweighted mean of neighbors. Otherwise it finds the smallest
strictly positive finite distance `min_pos` and notes whether any distance is
zero. Infinite distances never reach this function (Layer 5 drops them), so
`min_pos` is non-finite only when every neighbor is an exact duplicate, and
weights then stay equal (`src/impute_knn_brute.cpp`:
`if (!std::isfinite(min_pos))`). Two live branches: when some distance is
exactly zero (`src/impute_knn_brute.cpp`: `if (has_zero)`),
`floor = min_pos * sqrt(DBL_EPSILON)` and
`w_j = (floor / max(d_j, floor)) ^ dist_pow`, so duplicate columns get
weight 1 and everything else a tiny positive weight; otherwise
(`src/impute_knn_brute.cpp`: `std::pow(min_pos / top_k[j].distance`)
`w_j = (min_pos / d_j) ^ dist_pow`. Either way the largest weight is 1.

The `d_j` here is the Layer 5 quantity, so under `method = "euclidean"` it is a
mean SQUARED difference and `dist_pow` is an exponent of `2 * dist_pow` on the
square-rooted euclidean distance. The `min_pos` factor is common to every weight
and cancels out of the weighted average; it only keeps the weights near 1. Both
facts are stated on the help page (`R/knn_imp.R`:
`imputed values are distance-weighted averages`) and pinned against a pure-R
reference by a test (`tests/testthat/test-knn_imp.R`:
`the documented weighting rule reproduces`), which reproduces the kernel to
1e-12 on max absolute and max relative error over both metrics, three
`dist_pow` values, and two candidate layouts.

**Accumulation.** `impute_column_values()` (`src/imputed_value.cpp`:
`void impute_column_values(`) splits the neighbor ids into a masked bucket and a
complete bucket by comparing against `complete_start()`
(`src/imputed_value.cpp`: `if (nn_columns(j) < complete_start)`). Masked
neighbors contribute `w * nmiss[row]` to both the weighted sum and the weight
total (`src/imputed_value.cpp`: `for (arma::uword j : masked_j)`), so a neighbor
that is itself missing in a given row drops out for that row only. Complete
neighbors add `w * value` per row and a single shared `total_complete_w` to
every row's weight total. The final value is `weighted_sum / weight_total`, or
`NaN` when the total is zero (`src/imputed_value.cpp`: `ws_ptr[r] / wt_ptr[r]`),
which R converts to `NA` (`R/knn_imp.R`: `imputed_values[is.nan(imputed_values)]`).

That per-row drop-out is user-visible and is stated on the help page
(`R/knn_imp.R`: `Neighbors are selected once per column`): neighbors are chosen
once per column, but the set that actually contributes can differ from row to
row, and a row where none of the chosen `k` is observed comes back missing. The
"mixed candidates" half of the weighting-rule test
(`tests/testthat/test-knn_imp.R`: `"mixed candidates" = mixed_cand`) is the case
that exercises it.

## Feasibility and abort paths

`slideimp_infeasible` is a cli condition class raised by the R layer so the
wrappers can catch a per-group or per-window failure and fall back instead of
dying. `knn_imp()` has exactly two sites:

1. `k > n_elig - 1` (`R/knn_imp.R`: `if (k > n_elig - 1L) {`), where `n_elig` is
   the number of columns passing `colmax`; message "`k` (...) exceeds usable
   columns (...)". Exercised by tests (`tests/testthat/test-knn_imp.R`:
   `"exceeds usable columns"`; `tests/testthat/test-group_imp.R`:
   `Site 1: k > n_elig - 1`).
2. `length(grp_impute) == 0` (`R/knn_imp.R`: `if (length(grp_impute) == 0L) {`),
   i.e. every subset column with missing values also exceeds `colmax`; message
   "All subset columns with missing values exceed `colmax` (...)". Exercised by a
   test (`tests/testthat/test-group_imp.R`: `Site 2: all subset cols exceed colmax`)
   and, for the inclusive boundary, by another (`tests/testthat/test-knn_imp.R`:
   `colmax is inclusive`).

Everything else that aborts here is a plain error, not `slideimp_infeasible`:
`check_finite()` for `Inf` and all-`NA` columns (`tests/testthat/test-knn_imp.R`:
`Throw on Inf`, `"All NA/NaN"`), plus the C++ guards in `validate_knn_inputs()`
(`tests/testthat/test-knn_imp.R`: `the kernel checks the invariants`) and
`initialize_result_matrix()`. Catch sites are in the
wrappers only (next section); each switches on `on_infeasible` to abort (both raise a
new `cli_abort()` with the original as `parent` and the `slideimp_infeasible` class),
skip the block unchanged, or
fall back to `mean_imp_col()`.

## Data structures and ownership

- The R matrix is **not** copied at the boundary. `input_parameter<const
  arma::mat&>` in `src/RcppExports.cpp` resolves to
  `ConstReferenceInputParameter`, which for a double matrix is
  `ArmaMat_InputParameter<..., false_type>` in RcppArmadillo's
  `RcppArmadillo/interface/RcppArmadilloAs.h`: it constructs
  `arma::mat(ptr, nrow, ncol, false)` over an `Rcpp::NumericMatrix`, `copy_aux_mem =
  false`. The branch is chosen by the Armadillo element type, so a `double` matrix
  always takes it; the `ArmaMat_InputParameter<..., true_type>` specialization is
  for element types that need a cast. For double-storage
  input the `NumericMatrix` wraps R's own memory, so the kernel's `const arma::mat& obj`
  aliases the R matrix and is never written to. Integer-storage input is coerced to a
  fresh double copy when the `NumericMatrix` is built.
- The kernel returns a freshly allocated `(n_missing x 3)` `arma::mat` of
  triplets, not an imputed matrix. The only in-place mutation of the user's
  matrix happens in R (`R/knn_imp.R`: `obj[imp_indices] <- imputed_values[, 3]`).
- So the real full-matrix copies on this path are both on the R side: R's
  copy-on-modify at the scatter, and, when `mean_imp_col()` actually has work to
  do, the second `n x m` allocation in `mean_imp_col_internal()`
  (`src/mat_stats.cpp`: `arma::mat out(n_rows, mat.n_cols, arma::fill::none)`).
  The second one is skipped whenever nothing is still missing.
- Groups 1 and 2 are materialized a second time into `obj_masked` /
  `nmiss_masked` (`src/impute_knn_brute.cpp`: `arma::mat obj_masked(`) with NaNs
  zeroed and a parallel byte mask; `nmiss_masked` is `arma::Mat<uint8_t>`
  (`src/imputed_value.h`: `using mask_t = uint8_t;`), one byte per cell rather
  than a `umat`. Group 3 is read in place from `obj`, so a matrix dominated by
  complete columns pays no second copy.
- Column indices cross the boundary 0-based, R subtracting 1 (`R/knn_imp.R`:
  `grp_impute = as.integer(grp_impute - 1L)`), and come back 1-based along with
  the row indices, C++ adding 1 to both (`src/imputed_value.cpp`: `orig_col_1based = grp_impute(i) + 1`,
  `rows_ref(r) + 1`).
- Per-thread scratch is small and local: `top_k` (`src/impute_knn_brute.cpp`:
  `std::vector<NeighborInfo> top_k;`), whatever temporary buffer
  `std::stable_sort` takes for it, `nn_columns` and `weights`
  (`src/impute_knn_brute.cpp`: `arma::uvec nn_columns(n_neighbors);`), and the
  two accumulators in `impute_column_values` (`src/imputed_value.cpp`:
  `arma::vec weighted_sums(`). `mean_imp_col_internal()` (`src/mat_stats.cpp`:
  `arma::mat mean_imp_col_internal(`) is likewise out-of-place: a fresh matrix
  with untouched columns `memcpy`d, and `mean_imp_col()` reattaches dimnames
  (`R/mean_imp_col.R`: `dimnames(res) <- dimnames(obj)`).

## Parallelism

`par_for()` in `src/impute_knn_brute.cpp` runs over `[0, layout.n_imp)` -- one
task per target column. `n_threads` and `n_batches` are both set to `cores`,
floored at 1 (`src/impute_knn_brute.cpp`: `cores = std::max(1, cores);`), and
passed to the call, so the index range splits into as many contiguous batches as
threads. At one thread `par_for()` (`src/par_for.h`) runs the body inline
instead, for the reason given in Layer 3.

All threads read `obj`, `obj_masked`, `nmiss_masked`, and `grp_complete` without
writing them, and each writes only rows `[col_offsets(i), col_offsets(i+1))` of
`result`, a disjoint slice per target column, so no locks are taken. Interruption
is checked with `RcppThread::checkUserInterrupt()` on every fifth index
(`src/impute_knn_brute.cpp`: `if (i % 5 == 0)`), which RcppThread routes to the
main thread to unwind the pool. The optional progress bar is an
`RcppThread::ProgressBar` in a `unique_ptr` created only when `pb` is true and
incremented from the workers (`src/impute_knn_brute.cpp`: `++(*bar);`).

`src/Makevars` adds the OpenMP flags and links `RcppThread::LdFlags()`;
`src/Makevars.win` omits the RcppThread link line. `col_vars_internal()` and
`mean_imp_col_internal()` in `src/mat_stats.cpp` go through the same
`par_for()`.

## How the wrappers reach knn_imp()

`group_imp()` picks the function object (`R/group_imp.R`:
`imp_fn <- if (is_knn_mode) knn_imp else pca_imp`), pre-builds a per-group
parameter list injecting `cores` and a group-local `subset` (`R/group_imp.R`:
`params <- lapply(iter`), and invokes it via `do.call()` inside a `tryCatch` --
in the mirai branch (`R/group_imp.R`: `c(list(obj = sub_mat), task$params)`) and
in the sequential branch (`R/group_imp.R`: `c(list(obj = sub_mat), params[[i]])`).
Each group is a column slice of its `col_idx` (`R/group_imp.R`:
`src[, task$col_idx, drop = FALSE]`, `obj[, indices[[i]]$col_idx, drop = FALSE]`),
so `knn_imp()` sees a submatrix and group-local indices. Two kinds of group never
reach the backend: one whose feature columns hold no missing value
(`R/group_imp.R`: `cmiss <- unname(mat_miss(`, `needs_imp <- vapply(feat_splits`),
and one too small for any `k`, which `group_imp()` settles through
`on_infeasible` itself (`R/group_imp.R`: `small_now <- too_small & needs_imp`,
`groups too small to impute are settled here`) because `knn_imp()` rejects a
one-column matrix as invalid input, not as `slideimp_infeasible`.

`slide_imp()` calls `knn_imp()` directly (`R/slide_imp.R`:
`if (imp_method == "knn") {`) for each window `obj[, start[i]:end[i]]`
(`R/slide_imp.R`: `window_cols <- start[i]:end[i]`), passing `.progress = FALSE`
and a window-local `subset`, wrapped in the same `slideimp_infeasible` handler
(`R/slide_imp.R`: `slideimp_infeasible = function(e) {`).

## Test entry points

A test calls `impute_knn_brute()` directly and asserts the returned `(row, col)`
triplet locations match `which(is.na(...))` (`tests/testthat/test-knn_imp.R`:
`calculate the missing location correctly`); others cover `subset` by index and
by name with `post_imp` both ways (`tests/testthat/test-knn_imp.R`:
`works with post_imp = FALSE/TRUE`), all-NA rows and columns plus abort site 1
(`tests/testthat/test-knn_imp.R`: `extreme missing columns and rows`), the tie
rule (`tests/testthat/test-knn_imp.R`: `tied neighbors are settled by scan order`),
and `Inf` rejection (`tests/testthat/test-knn_imp.R`: `Throw on Inf`); and one
pins the `slideimp_infeasible` class on both abort sites
(`tests/testthat/test-group_imp.R`: `class is set on both abort sites`).
