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
`impute` is not a declared dependency (`.claude/CLAUDE.md`, "Layout"), as a
comment in the tests also says (`tests/testthat/test-knn_imp.R`:
`The comparison against impute::impute.knn`).

## Call chain overview

```
R/knn_imp.R                         knn_imp()  -- user entry: the column scan, knn_imp_parts(), the door
  +- R/col_scan.R                     col_scan()  -- one pass: holes and Inf per column
  |    +- src/mat_stats.cpp         col_scan_cpp() via R/RcppExports.R
  +- R/col_scan.R                     refuse_nonfinite()  (Inf + all-NA column refusal)
R/knn_imp.R                         knn_imp_parts()  -- validation, partition, the parts
  +- R/utils.R                      resolve_narrow_subset()  -- names/ints -> kept cols
  +- [abort site 1] R/knn_imp.R     k > n_elig - 1
  +- [abort site 2] R/knn_imp.R     every kept col with a hole exceeds colmax
  +- R/RcppExports.R                impute_knn_brute()  -- .Call shim
  |    +- src/RcppExports.cpp       _slideimp_impute_knn_brute (SEXP wrapper)
  |         +- src/knn_kernel.cpp        impute_knn_brute()   <-- C++ entry
  |              +- src/knn_kernel.cpp        validate_knn_inputs()
  |              +- src/knn_kernel.cpp        copy_with_mask lambda
  |              +- src/knn_kernel.cpp        initialize_result()
  |              +- src/par_for.h             par_for()
  |                   +- src/knn_kernel.cpp        distance_vector()
  |                   |    +- src/knn_kernel.cpp        distance_vector_impl<Metric>()
  |                   |         +- src/knn_kernel.cpp        calc_distance_raw<Metric, Bound>()
  |                   |         +- src/knn_kernel.cpp        calc_distance_raw_complete<Metric, Bound>()
  |                   |         +- src/knn_kernel.cpp        insert_if_better_than_worst()
  |                   +- src/knn_kernel.cpp        knn_weights()
  |                   +- src/knn_kernel.cpp        impute_column_values()
  |                   +- src/knn_kernel.cpp        masked_column_mean()  (fill_mean only)
  +- R/parts.R                      new_part()  -- the K-NN parts
  +- R/parts.R                      fallback_part()  (post_imp, columns over colmax)
R/parts.R                           fill_imp_parts() -> fill_parts()  -- in knn_imp()
R/knn_imp.R                         knn_imp_results() -> R/utils.R new_slideimp_results()
```

## Layer 1 -- R entry point (`R/knn_imp.R`)

`knn_imp()` (`R/knn_imp.R`: `knn_imp <- function(`) checks `obj`, takes its
scan and refuses on it, then makes two calls: its arguments and the scan
go to `knn_imp_parts()` (`R/knn_imp.R`: `knn_imp_parts <- function(`), which
does the rest of what this layer describes and returns the parts, and the
door fills them (`R/knn_imp.R`:
`knn_imp_results(fill_imp_parts(res), res$info)`). The internal takes the same
formals, defaults included, with `scan` second and without a default, so
`group_imp()` and `slide_imp()` pass it what they would pass `knn_imp()`, and
the scan of their block; a test holds the formals equal
(`tests/testthat/test-parts.R`: `the parts internals take their wrappers' formals`).
Its errors name `knn_imp_parts()`, except the `subset` refusals, which are
handed its caller's frame. It returns `obj` (the matrix the parts index, the
input as it was), `parts`, `keep`, `info`, and `n_miss`, the scan's hole
counts, which the door reads in place of counting again.

Validation, in order. In `knn_imp()` (`R/knn_imp.R`: `checkmate::assert_matrix(`):
a double matrix with `>= 1` row and `>= 2` columns, so an integer or logical
matrix is refused, never converted; then the scan and its refusal
(`R/knn_imp.R`: `refuse_nonfinite(scan)`). In `knn_imp_parts()`: `method`
matched against `"euclidean"` / `"manhattan"`; `k` integer in
`[1, ncol(obj) - 1]`; `cores >= 1`; `colmax` in `[0, 1]`; `dist_pow` a finite
number `>= 0`; the `post_imp` and `.progress` flags. All but two go through
`checkmate` with an explicit `.var.name`, so a rejection names the user's
argument. The exceptions are `method`, matched with `match.arg()`, and the
scan's refusal, a plain `stop()` (`R/col_scan.R`: `refuse_nonfinite <- function(`).
A test pins the `checkmate` path for `dist_pow` (`tests/testthat/test-knn_imp.R`:
`dist_pow is validated by name`).

Missing-value gating is one pass, the column scan (`R/col_scan.R`:
`col_scan <- function(`), through `col_scan_cpp()` (`src/mat_stats.cpp`:
`Rcpp::List col_scan_cpp(`), which every entry point of the package takes of
its input. It reads each column once and counts, never refuses: per column
the holes (`NA` and `NaN`; `Inf` is not a hole) and whether it holds `Inf` or
`-Inf`, and, when asked, the variance and the holes per row. `knn_imp()` asks
for neither. The policy is applied in R. `refuse_nonfinite()` refuses an `Inf`
anywhere, or a column with no observed value, with one message that names no
cell or column (`R/col_scan.R`:
`'obj' holds infinite values or columns with no observed value`).
`refuse_inf()` (`R/col_scan.R`: `refuse_inf <- function(`) refuses an `Inf`
alone, for the entry points that accept a column with no observed value:
`sample_na_loc()`, and the sampler's `abb_cols()` and flank plan. `NA` and
`NaN` are indistinguishable throughout; the policy is stated in the package
help (`R/slideimp-package.R`: `@section Missing values and non-finite input`).

`subset` goes through `resolve_narrow_subset()` (`R/knn_imp.R`:
`keep <- resolve_narrow_subset(subset, obj, call = parent.frame())`), as in `pca_imp()`: it accepts
`NULL` (all columns), a character vector matched against `colnames(obj)`, or
integerish indices, and returns the selected columns ascending as `keep`, or
`NULL` for every column. Each of these is a cli error that names `knn_imp()`:
an unknown name, a name given twice, names on an input without column names,
and an empty `subset`. Bad numbers go through `checkmate`. `keep` is the columns the
result holds, in input order: the door narrows the output to it (`R/parts.R`).

The early return fires if no column the result holds has a missing value
(`R/knn_imp.R`: `if (!any(cmiss[out_cols] > 0L)) {`). It returns no parts,
and the door hands back those columns as they are (`R/parts.R`:
`fill_imp_parts <- function(res) {`) through `knn_imp_results()`
(`R/knn_imp.R`: `knn_imp_results <- function(obj, info) {`), the input itself
without `subset` and `obj[, keep, drop = FALSE]` with it. The record is
`knn_imp_info()`'s (`R/knn_imp.R`: `knn_imp_info <- function(`), the
same fields the main path sets, so the class and attributes of the result do
not depend on which path the data takes; a test pins that
(`tests/testthat/test-knn_imp.R`:
`early returns carry the slideimp_results class`). Its message is "No missing
values in input. Returning input unchanged." without `subset`, as in
`pca_imp()`, and says the `subset` columns come back unchanged with it. The
class costs one copy of the input: a matrix the caller still holds cannot take
an attribute in place (`R/utils.R`:
`Attaching anything to a matrix the caller still holds copies it`).

Missingness and partitioning. `cmiss` is the scan's hole count
(`R/knn_imp.R`: `cmiss <- scan$n_miss`). Without the variance the scan
counts `std::isnan` per column, with no branch on a hole, and does not count
`Inf` as missing (`src/mat_stats.cpp`: `const int hole = std::isnan(v);`);
`cmiss` decides the early return above. `colmax` then becomes an eligibility mask,
`eligible <- miss_rate <= min(colmax, 1)`; the test is inclusive, matching
`pca_imp()` (`R/pca_imp.R`: `eligible <- miss_rate <= min(colmax, 1)`), so a
column strictly above the threshold is dropped from the K-NN stage entirely --
neither imputed nor usable as a neighbor. Three disjoint index groups follow:
`grp_impute` (eligible, has missing values, in `subset` -- these get imputed),
`grp_miss_no_imp` (eligible, has missing values, not in `subset` -- neighbor
candidate only), and `grp_complete` (eligible, fully observed -- neighbor
candidate only). `method` is remapped to the integer the kernel expects, `0` =
euclidean, `1` = manhattan (`R/knn_imp.R`: `method <- switch(method`). So for
the kernel `subset` narrows only the *target* set: every eligible column stays
a neighbor candidate regardless of `subset`. Every target column has the same
candidates with or without `subset`, and only their scan order moves. Scan
order decides only between tied neighbors (Layer 5, Ties), so on untied data the
narrowed result is the full result's columns, bit for bit, as a test asserts
(`tests/testthat/test-knn_imp-parts.R`:
`subset returns the full output's columns, in input order`).

`rowmax` is not an argument of `knn_imp()`: it belongs to the tuning path
(`R/tune_imp.R`: `sample_na_loc <- function(`, `max_row_miss <- max_miss_count(`),
capping injected missingness per row.

Nothing is converted before the kernel is called. `obj` is double from the
door, so the kernel and the door read the matrix the caller passed, and Rcpp
makes no copy of it (Data structures, below).

**The parts.** The kernel passes `fill_mean = post_imp` (`R/knn_imp.R`:
`fill_mean = post_imp,`), so under `post_imp` it fills each cell it leaves with
its column's observed mean itself (Layer 6). `knn_imp()` then builds parts
(`R/parts.R`), which together cover every column the result holds:

- The K-NN part: the kernel's `loc` as its map and `value` as its values, over
  `grp_impute` and the eligible complete columns kept (`R/knn_imp.R`:
  `complete_kept <- intersect(eligible_complete, out_cols)`), which have nothing
  to fill. `grp_impute` lies inside `keep`, so this part needs no cut.
- Under `post_imp = FALSE`, when the kernel left cells, the columns holding them
  are a second K-NN part (`R/knn_imp.R`: `left_cols <- unique(loc[left, "col"])`),
  `left_na`, whose map holds the cells K-NN did fill, since every column of a
  `left_na` part must keep a hole (`R/parts.R`:
  `full <- which(listed >= n_miss[part$cols])`).
- The columns over `colmax` the result holds (`R/knn_imp.R`:
  `mean_cols <- setdiff(out_cols, eligible_idx)`): a mean part from
  `fallback_part()` under `post_imp` (`R/knn_imp.R`:
  `fallback_part(obj, mean_cols, "mean")`), whose mean is `.colMeans()` of the
  column's observed values, and a `left_na` part with an empty map without it.
  Each has a hole, since its missing rate is above `colmax`, and an observed
  value, since the door refuses a column with none (`refuse_nonfinite()`).

`knn_imp()` hands them to `fill_parts()` through `fill_imp_parts()`, which
copies the kept columns once and writes every part's values;
`group_imp()` and `slide_imp()` take them instead (below). `knn_imp_results()`
attaches the
`slideimp_results` class and attributes through `new_slideimp_results()`
(`R/utils.R`: `new_slideimp_results <- function(`), which first strips its
input to a bare matrix, so a `slideimp_results` passed back in for a second pass
does not carry its class or its old attributes into the new result. Tests:
`tests/testthat/test-knn_imp-parts.R`, and `post_imp fills columns K-NN could
not finish` in `tests/testthat/test-knn_imp.R`.

**What post-imputation records.** Nothing in the values tells a column mean
from a K-NN average afterwards, so the two sets go onto the result as
attributes, in the result's column numbering: the columns over `colmax` as
`mean_filled_cols`, the cells the kernel left (its `left` rows) as
`mean_filled_cells`, a two-column `row`/`col` matrix built by
`mean_filled_cells_matrix()` (`R/utils.R`:
`mean_filled_cells_matrix <- function(`). Both are always present and empty
when `post_imp` filled nothing, early returns included. `slide_imp()` and
`group_imp()` map the record from each window's or group's `info` to the
result's columns through one builder (`R/utils.R`:
`mean_filled_record <- function(info, to) {`): `slide_imp()` hands it the
window's targets' positions in its `subset` (`R/slide_imp.R`:
`mean_filled[[i]] <- mean_filled_record(info, tj)`), and `group_imp()` its
features' columns in the result (`R/group_imp.R`:
`lapply(feat_splits, match, table = result_cols)`). `print()`
reports them. Pinned by tests
(`tests/testthat/test-knn_imp.R`:
`knn_imp records what post_imp set to a column mean`;
`tests/testthat/test-knn_imp-parts.R`: `the record counts the result's columns`).

## Layer 2 -- the .Call boundary

`R/RcppExports.R` defines the R-side shim `impute_knn_brute(obj, k,
grp_impute, grp_miss_no_imp, grp_complete, method, dist_pow, fill_mean, cores,
pb)`, whose body is a single `.Call` to `_slideimp_impute_knn_brute`.
`fill_mean` defaults to `FALSE`, so a direct call without it leaves the cells
it cannot fill `NA`.

The generated C wrapper is in `src/RcppExports.cpp`; it converts each SEXP
through `Rcpp::traits::input_parameter<...>` and calls the C++ symbol.
Registration is its entry in the `CallEntries` table, installed by
`R_init_slideimp`. The exact C++ entry symbol is `impute_knn_brute`, defined in
`src/knn_kernel.cpp` under the `// [[Rcpp::export]]` marker. It is the
only C++ entry point the K-NN stage needs; `col_scan_cpp`, the scan, is
a separate `.Call` made from R before it.

## Layer 3 -- C++ translation units and headers on the path

Three compiled `.cpp` files are reached: `src/RcppExports.cpp` (generated SEXP
wrappers), `src/knn_kernel.cpp` (the whole K-NN kernel) and
`src/mat_stats.cpp` (`col_scan_cpp`, the scan).

`src/knn_kernel.cpp` has no header of its own. It defines `MaskByte`,
`GroupLayout`, `validate_knn_inputs()`, the distance and neighbor code, the
weights, the result layout and the accumulation loop
(`src/knn_kernel.cpp`: `KnnResult initialize_result(`,
`void impute_column_values(`) in one anonymous namespace, so the exported
`impute_knn_brute()` at the end of the file is the only name it gives the
linker.

Project headers reached, following the `#include` graph:

- `src/matrix_checks.h` -- pulled in by `src/knn_kernel.cpp` and
  `src/mat_stats.cpp` and `src/sample_each_rep_cpp.cpp`. Holds `col_ptr()`,
  the column-index refusal that `validate_knn_inputs()` and the sampler share,
  and a refusal of a negative count.
- `src/par_for.h` -- included by `src/knn_kernel.cpp` only. One
  function template: `par_for()` runs the loop body
  inline when `n_threads <= 1` and hands off to `RcppThread::parallelFor()`
  otherwise. The reason is in the header's own comment -- RcppThread's free
  `parallelFor()` sets the global pool's thread count and restores it, and
  quickpool implements a change of count by joining every worker and spawning
  new ones, so a one-thread call rebuilds the pool twice for nothing. That is
  the default configuration of the whole chain: `slide_imp()` runs each window
  at one core and `group_imp()` forces one core per worker under mirai. The
  serial and parallel paths are pinned equal by a test
  (`tests/testthat/test-knn_imp.R`: `the parallel path gives the same answer`).
- `src/loc_timer.h` -- included by `src/knn_kernel.cpp`. Unless
  `LOC_TIMER` is defined, the statement macros compile to `((void)0)` and
  `LOC_TIMER_PARAM` / `LOC_TIMER_ARG` to nothing, so the `LOC_*` calls in
  `src/knn_kernel.cpp` vanish in a normal build.

Plus the external `Rcpp.h` (included by `src/knn_kernel.cpp`) and
`RcppThread.h` (by `src/knn_kernel.cpp`). The remaining `src/` headers
(`eig_sym_sel.h`, `gram_ops.h`, `hybrid_topk_eig.h`, `lobpcg_warm.h`,
`pca_linalg_utils.h`, `svd_triplet.h`) and `.cpp` files (`pca_kernel.cpp`,
`find_windows.cpp`, `find_windows_flank.cpp`, `sample_each_rep_cpp.cpp`) serve
the PCA, windowing, and simulation paths and are not reached from `knn_imp()`.

**Index types.** `k`, `method`, `cores` and the group indices cross from R as
`int`, and every row and column index the kernel derives from them is `int`
too, bounded by `nrow()` or `ncol()`, which R holds as `int`. The offset of
a column into a buffer is taken in `size_t` by `col_ptr()`
(`src/matrix_checks.h`: `inline T *col_ptr(`), so a matrix past `2^31` cells
indexes correctly. The result offsets are taken in `size_t` too
(`src/knn_kernel.cpp`: `std::size_t sum_missing = 0;`). The result is an R matrix, whose dims are
`int`, so more than `2^31 - 1` missing cells in `grp_impute` is refused
(`src/knn_kernel.cpp`: `has a row for each and an R matrix holds at most %d rows`).

## Layer 4 -- kernel setup and result skeleton

1. `validate_knn_inputs()` (`src/knn_kernel.cpp`: `validate_knn_inputs(const Rcpp::NumericMatrix &obj, const int k,`)
   re-checks most of what R already checked, so the kernel is safe when called
   directly. All-`NA` columns are left to the zero-`n_valid` guard
   (`src/knn_kernel.cpp`: `if (target_n_valid == 0.0)`), `cores` is
   floored rather than validated (`src/knn_kernel.cpp`:
   `cores = std::max(1, cores);`), and the `k` bound applies only when a column
   needs imputing. It checks shape, `k >= 1`, `method` in `{0, 1}`, `dist_pow`
   finite and non-negative, every group index in `[0, ncol)` (a negative index
   and `NA`, which is `INT_MIN`, included), the three groups
   disjoint, every `grp_complete` column free of `NaN`, and `k <= n_working - 1`
   for `n_working` = sum of the three group sizes. The partition and `NaN`
   checks are not redundant with the R layer: group 3 is read unmasked by
   `impute_column_values()`, so a `NaN` there would flow into the weighted sum
   and produce a missing value with no error, and a column claimed by two groups
   would be its own neighbor at distance zero. A test exercises both, and the
   index and `k` refusals too (`tests/testthat/test-knn_imp.R`:
   `the kernel checks the invariants`). The kernel does not look for
   infinities: the scan at the door refuses them before the call
   (`src/knn_kernel.cpp`: `No scan for infinite values here`).
2. `GroupLayout layout{...}` in `src/knn_kernel.cpp` is the single
   source of truth for group boundaries: groups 1 and 2 occupy local columns
   `[0, n_masked())` of the working buffers, and `complete_start()`
   (`src/knn_kernel.cpp`: `int complete_start() const`) is a *virtual*
   index used as a tag -- a neighbor id `>= complete_start()` means "group 3,
   offset `id - complete_start()` into `grp_complete`".
3. Working buffers `obj_masked` (`std::vector<double>`), `nmiss_masked`
   (`std::vector<uint8_t>`), both flat and column-major, and `n_col_valid` are
   allocated (`src/knn_kernel.cpp`: `std::vector<double> obj_masked(n_cells);`) and filled by
   `copy_with_mask`, which calls `std::isnan` once per cell and reuses the
   result as the mask byte, the count increment, and the zeroing test
   (`src/knn_kernel.cpp`: `MaskByte finite = !std::isnan(v);`). Missing
   cells therefore hold a literal `0.0` in `obj_masked` and are excluded by the
   mask, never by a branch in the distance loop.
4. Groups 1 and 2 are copied (`src/knn_kernel.cpp`:
   `copy_with_mask(i, grp_impute[i])`, `copy_with_mask(layout.mni_start() + p`).
   Group 3 is **not** copied: the kernel reads it straight out of `obj` through
   `grp_complete` (`src/knn_kernel.cpp`: `group 3: read directly from obj`),
   and in the lambdas (`src/knn_kernel.cpp`: `auto complete_dist_unbounded`).
5. `initialize_result()` (`src/knn_kernel.cpp`:
   `KnnResult initialize_result(`) runs between the group-1 and group-2
   copies so the call can bail early when there is nothing to impute
   (`src/knn_kernel.cpp`: `if (res.value.size() == 0)`).

6. Under `fill_mean`, a `grp_impute` column with no observed value is refused
   after the group-1 copy (`src/knn_kernel.cpp`:
   `'grp_impute' column %d has no observed value, so it has `), since it has no
   mean to fall back on; the scan at the door keeps such a column away from the R
   entry point, so only a direct call meets it.

**Result skeleton.** The return value is `list(loc, value, left)`: `loc` an
integer matrix of 1-based row and original 1-based column, one row per missing
cell of `grp_impute`, with no dimnames: R gives it the column names `row` and
`col` a part's map has (`R/knn_imp.R`:
`dimnames(imputed$loc) <- list(NULL, c("row", "col"))`), since an Rcpp
dimnames assignment is an rchk error in Rcpp's own header; `value`
a double vector of the same length; and `left` the 1-based rows of `loc` K-NN
found no value for. `loc` and `value` are R
objects allocated here, on the main thread, before any worker starts, and so is
a byte per entry that the workers set for the cells they leave
(`src/knn_kernel.cpp`: `res.left.assign(sum_missing, 0);`), collected into
`left` after the loop. Missing
counts per target column come from `n_col_valid` rather than a second pass
over the mask (`src/knn_kernel.cpp`: `n_miss_i = n_rows - n_col_valid[i]`),
and a prefix sum builds `col_offsets` (`src/knn_kernel.cpp`:
`col_offsets[i + 1] = sum_missing;`) so each target column owns a contiguous
slice of `value` and threads never collide. `value` starts as `NA`
(`src/knn_kernel.cpp`: `value = Rcpp::NumericVector(n_total, NA_REAL);`),
`loc` is filled up front, in one scan down each mask column
(`src/knn_kernel.cpp`: `orig_col_1based = grp_impute[i] + 1`). Its first
column is also the one record of which rows of target column `i` are missing:
the worker reads entries `col_offsets[i]` to `col_offsets[i + 1]` of it
through a plain pointer and subtracts the 1 (`src/knn_kernel.cpp`:
`const int row_idx = loc_rows[r] - 1;`).
Empty cases short-circuit (`src/knn_kernel.cpp`: `if (n_imp == 0)`,
`if (sum_missing == 0)`). It checks nothing that its one caller has just
established: `copy_with_mask` adds to `n_col_valid` the same byte it writes to
the mask, so column `i` of the mask holds exactly `n_rows - n_col_valid[i]`
zeros (`src/knn_kernel.cpp`: `so it holds exactly as many zeros as its slice`). Its one
refusal is a result too large for an R matrix, above.

## Layer 5 -- neighbor search (`src/knn_kernel.cpp`)

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
`DIST_INF`, +infinity (`src/knn_kernel.cpp`:
`constexpr double DIST_INF`), when the two columns share no observed row
(`src/knn_kernel.cpp`: `if (n_valid == 0.0)`). No square root is taken for
the euclidean case, so the stored quantity is a mean squared difference.
`calc_distance_raw_complete<Metric, Bound>` is the group-3 variant: the other
side is fully observed, so only the target mask matters
(`src/knn_kernel.cpp`: `target_nmiss[rr] * Metric::accumulate(diff)`) and
`n_valid`, the target's own observed count, is known up front and passed in. It
has no zero-`n_valid` branch of its own -- see the guard below.

**Pruning.** Both kernels chunk rows in blocks of `GRAIN = 16` and, when `Bound`
is true, test a partial-sum bound at each chunk boundary, returning `inf` on a
prune (`src/knn_kernel.cpp`: `if (dist > worst_dist * (n_valid`,
`if (dist > unnorm_bound)`), with derivations in the comments
(`src/knn_kernel.cpp`: `Bound math (conservative but correct)`,
`Bound math (exact)`). The check sits under `if constexpr (Bound)`, so it is not
compiled into the `Bound = false` instantiations at all; their callers leave
`worst_dist` out, since it is never read there. The
pruned and unpruned paths are compared against each other implicitly by a test
(`tests/testthat/test-knn_imp.R`: `the documented weighting rule reproduces`),
where the pure-R reference computes every distance unpruned and the kernel
prunes past the first `k`.

**Top-k selection.** `distance_vector_impl` keeps a `std::vector<NeighborInfo>`
of capacity `k` (`src/knn_kernel.cpp`: `std::vector<NeighborInfo> top_k;`)
and runs two phases with three shared cursors `p1`, `p2`, `c` so every candidate
column is visited exactly once. Before any of it, a guard
(`src/knn_kernel.cpp`: `if (target_n_valid == 0.0)`) returns an empty
vector when the target column has no observed row: it then overlaps nothing, and
`calc_distance_raw_complete` would divide by a zero `n_valid`. The scan
keeps that case away from the R entry point, so the guard exists for direct
callers of the kernel; one is a test (`tests/testthat/test-knn_imp.R`:
`all_na$obj[, 1] <- NA`).

- *Fill* (`src/knn_kernel.cpp`: `c < layout.n_complete && top_k.size() < k`):
  append the first `k` candidates with no pruning, in the
  order group 3, masked columns below the target, masked columns above it. `p2`
  starts at `index + 1` and the `p1` loop stops at `index`, so the target column
  is never its own neighbor. Then one `std::stable_sort` on those `k` entries.
  Ties, below, says why the sort is stable.
- *Replacement* (`src/knn_kernel.cpp`: `for (; c < layout.n_complete; ++c)`):
  resume each cursor, compute the distance with `Bound = true` against the
  current worst (`top_k.back().distance`), and hand it to
  `insert_if_better_than_worst`, which rejects anything not strictly better than
  the worst, overwrites the last slot, and bubbles it into position.
- *Tail drop* (`src/knn_kernel.cpp`: `while (!top_k.empty()`): the result
  is sorted ascending, so candidates at distance `+inf` (no observed row in
  common with the target) sit at the tail and are popped before returning. The
  caller can therefore receive fewer than `k` neighbors, or none at all; when
  none (`src/knn_kernel.cpp`: `if (!top_k.empty())`), every entry of the
  column is left: its bytes are set (`src/knn_kernel.cpp`:
  `left_ptr[offset + r] = 1;`), and its entries of `value` keep their initial
  `NA`, or under `fill_mean` take the column's observed mean
  (`src/knn_kernel.cpp`: `fill_mean ? masked_column_mean(data, col) : 0.0;`).
  The column still advances the progress bar.

The structure is therefore a sorted fixed-size array with insertion, not a heap
and not a full sort of all candidates: exactly one `std::stable_sort`, over `k`
elements.

**Ties.** Among candidates at equal distance, the one met first in scan order
wins. Scan order is the fill order above: group 3 in `grp_complete` order, then
the masked columns in local order (`grp_impute`, then `grp_miss_no_imp`),
skipping the target. Two things make it hold. The fill appends in scan order,
and `std::stable_sort` keeps equal distances in that order. After that both
comparisons are strict `<` -- the bubble (`src/knn_kernel.cpp`:
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
common, the two rules give different values in 2-29% of imputed cells
(measured).

The sort must be stable: `std::sort` leaves equal elements in an unspecified
order, so the standard library would choose which tied fill entry is evicted.
`std::stable_sort` may take a temporary buffer per target column; no cost above
run-to-run noise was measured. Pinned by a test (`tests/testthat/test-knn_imp.R`:
`tied neighbors are settled by scan order`), which reads the chosen neighbor set
back out of the imputed value at `k = 17` and `k = 32`.

## Layer 6 -- weights and imputation

**Weights.** `knn_weights()` (`src/knn_kernel.cpp`:
`std::vector<double> knn_weights(`) returns a vector of ones when `dist_pow == 0`
-- the plain unweighted mean of neighbors. Otherwise it finds the smallest
strictly positive distance `min_pos` and notes whether any distance is
zero. Infinite distances never reach this function (Layer 5 drops them), so
`min_pos` keeps its starting `+inf` only when every neighbor is an exact
duplicate, and weights then stay equal (`src/knn_kernel.cpp`:
`if (min_pos == DIST_INF)`). The rest is one loop with two cases
(`src/knn_kernel.cpp`: `has_zero ? min_pos * std::sqrt(`,
`std::pow(ref / std::max(top_k[j].distance, ref), dist_pow)`): when some
distance is exactly zero, `ref = min_pos * sqrt(DBL_EPSILON)` and
`w_j = (ref / max(d_j, ref)) ^ dist_pow`, so duplicate columns get
weight 1 and everything else a tiny positive weight; otherwise `ref = min_pos`,
which no distance is below, and
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

**Accumulation.** `impute_column_values()` (`src/knn_kernel.cpp`:
`void impute_column_values(`) takes the neighbors as the search returned
them and goes over them twice, the masked ones first and the complete ones
after, telling them apart by comparing the id against `complete_start()`
(`src/knn_kernel.cpp`: `if (local >= complete_start)`). Masked
neighbors contribute `w * nmiss[row]` to both the weighted sum and the weight
total (`src/knn_kernel.cpp`: `double weight = w * valid;`), so a neighbor
that is itself missing in a given row drops out for that row only. Complete
neighbors add `w * value` per row and a single shared `total_complete_w` to
every row's weight total. The final value is `weighted_sum / weight_total`
(`src/knn_kernel.cpp`: `ws_ptr[r] / wt_ptr[r]`), written only where the
total is positive. Elsewhere the entry is left (`src/knn_kernel.cpp`:
`left[r] = 1;`): it keeps the `NA` that `initialize_result()` wrote, so R
receives its own missing value and converts nothing, or under `fill_mean`
takes the column's observed mean, computed at the first such row and only then
(`src/knn_kernel.cpp`: `col_mean = masked_column_mean(d, col);`).

**The mean.** `masked_column_mean()` (`src/knn_kernel.cpp`:
`double masked_column_mean(const KnnData &d, const int i) {`) sums the
target's masked column, whose holes hold `0.0`, in double and in row order, and
divides by its observed count, `n_col_valid`. It is taken by the worker that
owns the column, from the copy the kernel already holds, so it costs one pass
over that column and nothing of the input's size. The R fallback for the
columns over `colmax` takes `.colMeans()`, which sums and divides in long
double, so the two can differ in the last bits; they never fill the same
column. Pinned by tests (`tests/testthat/test-knn_imp.R`:
`fill_mean sets each cell K-NN leaves to its column's observed mean`,
`post_imp fills a cell K-NN leaves with the mean of`). The first reports the
difference from `.colMeans()` and asserts it under 1e-14 absolute and 1e-13
relative.

That per-row drop-out is user-visible and is stated on the help page
(`R/knn_imp.R`: `Neighbors are selected once per column`): neighbors are chosen
once per column, but the set that actually contributes can differ from row to
row, and a row where none of the chosen `k` is observed comes back missing, or
holding its column mean under `post_imp`. The
"mixed candidates" half of the weighting-rule test
(`tests/testthat/test-knn_imp.R`: `"mixed candidates" = mixed_cand`) is the case
that exercises it.

## Feasibility and abort paths

`slideimp_infeasible` is a cli condition class raised by the R layer so the
wrappers can catch a per-group or per-window failure and fall back instead of
dying. `knn_imp()` has exactly two sites:

1. `k > n_elig - 1` (`R/knn_imp.R`: `if (k > n_elig - 1L) {`), where `n_elig` is
   the number of columns passing `colmax`; message "`k` (...) must be below the
   number of usable columns (...)". Exercised by tests
   (`tests/testthat/test-knn_imp.R`:
   `"must be below the number of usable columns"`; `tests/testthat/test-knn_imp.R`:
   `Site 1: k > n_elig - 1`).
2. `length(grp_impute) == 0` (`R/knn_imp.R`: `if (length(grp_impute) == 0L) {`),
   i.e. every subset column with missing values also exceeds `colmax`; message
   "No column to impute is within `colmax` (...)". Exercised by a
   test (`tests/testthat/test-knn_imp.R`: `Site 2: all subset cols exceed colmax`)
   and, for the inclusive boundary, by another (`tests/testthat/test-knn_imp.R`:
   `colmax is inclusive`).

Everything else that aborts here is a plain error, not `slideimp_infeasible`:
the scan's refusal of `Inf` and all-`NA` columns, at the door
(`tests/testthat/test-knn_imp.R`: `Throw on Inf`,
`Behavior with extreme missing columns and rows`), the `subset` and
names refusals above (`tests/testthat/test-knn_imp-parts.R`:
`subset refuses what it cannot select, before the kernel`), plus the C++ guards in `validate_knn_inputs()`
(`tests/testthat/test-knn_imp.R`: `the kernel checks the invariants`),
the `fill_mean` refusal of a column with no observed value, and
the refusal in `initialize_result()` of a result past an R matrix's rows. A
part that fails its checks is a package defect, a plain `stop()` with the
prefix "slideimp internal, part:". Catch sites are in the
wrappers only (next section); each switches on `on_infeasible` to abort, skip the
block, or fall back to column means over the block's targets. An abort raises a new
`cli_abort()` of class `slideimp_infeasible`; `slide_imp()` and the sequential path
of `group_imp()` pass the original as `parent`, while the mirai path of
`group_imp()` gets only the message back from the daemon and raises without one. A
skip adds nothing to `slide_imp()`'s sums and gives the group no part in
`group_imp()`. The mean fallback is one `fallback_part()` over the block's targets
with a hole, in `slide_imp()` and in `impute_group()`.

## Data structures and ownership

- The R matrix is **not** copied at the boundary. `input_parameter<const
  Rcpp::NumericMatrix&>` in `src/RcppExports.cpp` resolves to
  `ConstReferenceInputParameter` (Rcpp's `InputParameter.h`), which builds the
  `NumericMatrix` from the SEXP; that goes through `r_cast<REALSXP>`, which
  returns a double matrix as it is (READ in Rcpp's `vector/Vector.h`). So the
  kernel's `obj` is R's own memory, read in place and never written to.
  Integer-storage input would be coerced to a fresh double copy there, so
  every entry point refuses it at its door, and none converts it (Layer 1).
  That the kernel then makes no further copy is READ, not measured.
- The kernel returns `list(loc, value, left)`, a fresh `n_missing x 2`
  integer matrix, a double vector and the left rows, not an imputed matrix. R
  builds parts from them and holds no copy of the input until a door.
- So the one full-matrix copy on the matrix output is the matrix door's: the
  first write into the input without `subset` (`R/parts.R`:
  `out[at] <- parts[[i]]$values`), or `obj[, keep, drop = FALSE]` of the kept
  columns with it (`R/parts.R`: `out <- obj[, keep, drop = FALSE]`). The door
  also forms each part's linear index once, in double, and counts holes before
  and after the write. Measured by `dev/probe/probe_knn_imp_peak.R` (one run per
  stage, peak working set over a 2000 x 4000 input, `k = 10`, 4 cores): with
  `post_imp` 2.07 x the input; with `subset` over a tenth of the columns at
  most 0.95 x.
- Groups 1 and 2 are materialized a second time into `obj_masked` /
  `nmiss_masked` (`src/knn_kernel.cpp`: `std::vector<double> obj_masked(n_cells);`)
  with NaNs zeroed and a parallel byte mask; `nmiss_masked` holds `MaskByte`
  (`src/knn_kernel.cpp`: `using MaskByte = std::uint8_t;`), one byte per cell.
  Both are zero-filled by `std::vector` before `copy_with_mask` writes every
  cell. Group 3 is read in place from `obj`, so a matrix dominated by
  complete columns pays no second copy.
- Column indices cross the boundary 0-based, R subtracting 1 (`R/knn_imp.R`:
  `grp_impute = grp_impute - 1L,`), and come back 1-based along with
  the row indices, C++ adding 1 to both (`src/knn_kernel.cpp`: `orig_col_1based = grp_impute[i] + 1`,
  `loc_row[w] = r + 1;`).
- Per-thread scratch is small and local: `top_k` (`src/knn_kernel.cpp`:
  `std::vector<NeighborInfo> top_k;`), whatever temporary buffer
  `std::stable_sort` takes for it, `weights`
  (`src/knn_kernel.cpp`: `weights = knn_weights(top_k, dist_pow);`), and the
  two accumulators in `impute_column_values` (`src/knn_kernel.cpp`:
  `std::vector<double> weighted_sums(`). The mean part allocates in
  proportion to the columns it fills: `fallback_part()` copies the kept
  columns over `colmax` out, with their `is.na()` (`R/parts.R`:
  `sub <- obj[, cols, drop = FALSE]`), and keeps the holes and their values.

## Parallelism

`par_for()` in `src/knn_kernel.cpp` runs over `[0, layout.n_imp)` -- one
task per target column. `n_threads` and `n_batches` are both set to `cores`,
floored at 1 (`src/knn_kernel.cpp`: `cores = std::max(1, cores);`), and
passed to the call, so the index range splits into as many contiguous batches as
threads. At one thread `par_for()` (`src/par_for.h`) runs the body inline
instead, for the reason given in Layer 3.

A worker never touches an Rcpp object or the R API: Rcpp's `operator[]` can
call R's `warning()` on an index past the end (Rcpp's `vector/traits.h`,
`check_index`), and R is single-threaded. The loop body's lambda therefore
captures by name only plain pointers taken before the loop, `std` containers,
scalars and `data`, a `KnnData` holding the pointers and counts the kernel
reads (`src/knn_kernel.cpp`: `struct KnnData {`,
`[&data, k, method, dist_pow, fill_mean, value_ptr, left_ptr, loc_row_ptr,`), so a
use of `obj`, `grp_complete` or `value` inside it would not compile; the
functions it calls take the same struct.

All threads read `obj`, `obj_masked`, `nmiss_masked`, `grp_complete` and the
first column of `loc` without writing them, and each writes only entries `[col_offsets[i], col_offsets[i+1])`
of `value` and of the left bytes, a disjoint slice per target column, so no
locks are taken. The bytes are a `std::vector<std::uint8_t>`, whose elements
are distinct memory locations, unlike a `std::vector<bool>`. Interruption
is checked with `RcppThread::checkUserInterrupt()` on every fifth index, `INTERRUPT_EVERY`
(`src/knn_kernel.cpp`: `if (i % INTERRUPT_EVERY == 0)`), which RcppThread routes to the
main thread to unwind the pool. The optional progress bar is an
`RcppThread::ProgressBar` in a `unique_ptr` created only when `pb` is true and
incremented from the workers (`src/knn_kernel.cpp`: `++(*bar);`).

`src/Makevars` links `RcppThread::LdFlags()`; `src/Makevars.win` omits
that link line. Neither passes the OpenMP flags. The K-NN kernel is the
only caller of `par_for()`; the column statistics in `src/mat_stats.cpp`
run on one thread.

## How the wrappers reach knn_imp()

`group_imp()` reaches the internal, not `knn_imp()`: it picks the function
object (`R/group_imp.R`:
`imp_fn <- if (is_knn_mode) knn_imp_parts else pca_imp_parts`), pre-builds a
per-group parameter list injecting a group-local `subset` for both backends
and `cores` for K-NN (`R/group_imp.R`:
`p$subset <- seq_len(feat_lengths[[i]])`), and runs each group through
`impute_group()` (`R/group_imp.R`: `impute_group <- function(`), in both
running modes. It calls the internal via `do.call()` inside a `tryCatch`
(`R/group_imp.R`: `c(list(obj = sub_mat, scan = sub_scan), params)`)
on a column slice of the group's `col_idx` without dimnames
(`R/group_imp.R`: `group_slice <- function(obj, col_idx) {`), features
first, so the internal sees a submatrix and group-local indices. The slice's
scan is the door's, indexed by `col_idx` (`R/col_scan.R`:
`scan_cols <- function(scan, j) {`), so no slice is read to check it. Its parts
are cut to the features and moved into `obj`'s numbering (`R/parts.R`:
`rebase_parts <- function(parts, col_idx, keep = NULL) {`), which sorts them
back into `obj`'s column order, since a group names its features in any
order. A group the backend refuses has no part under `"skip"` and one mean
part under `"mean"`. Every group's parts are written in one `fill_parts()`
at the end (`R/group_imp.R`:
`obj <- fill_parts(obj, parts, keep = result_cols, n_miss = cmiss)`), with the
door's hole counts, so every group reads `obj` as it was passed, aux columns
included, and the result holds only the `subset` columns when `subset` is
given. Under mirai each task carries its slice and its scan by value and
hands back its parts (`R/group_imp.R`:
`sub_mat = group_slice(obj, group_cols[[i]]),`). Two kinds of group
never reach the backend: one whose feature columns hold no missing value
(`R/group_imp.R`: `cmiss <- scan$n_miss`, `needs_imp <- vapply(feat_splits`),
and one too small for any `k`, which `group_imp()` settles through
`on_infeasible` itself (`R/group_imp.R`: `small_now <- too_small & needs_imp`,
`Groups too small to impute are settled here`) because the backend rejects
a one-column matrix as invalid input, not as `slideimp_infeasible`.

`slide_imp()` calls the internal too, `knn_imp_parts()` (`R/slide_imp.R`:
`if (imp_method == "knn") {`), for each window `obj[, start[i]:end[i]]`
(`R/slide_imp.R`: `window_cols <- start[i]:end[i]`), passing the window's
scan, which is the door's scan indexed by the window's columns (`R/slide_imp.R`:
`sub_scan <- scan_cols(scan, window_cols)`), `.progress = FALSE`
and a window-local `subset`, ascending, wrapped in a `slideimp_infeasible`
handler (`R/slide_imp.R`: `slideimp_infeasible = function(e) {`);
`pca_imp_parts()` windows take the same `subset`. Windows overlap, so their
parts are not written as they are: `window_hole_values()` fills a window's
parts into a copy of its targets and reads the value at each hole, and the
values are summed into `acc`, one double per hole of the `subset` columns
(`R/slide_imp.R`: `acc[at] <- acc[at] + window_hole_values(res, holes)`).
A window's targets are consecutive in `subset`, so their holes are one run of
`acc`, and the first window that holds a column keeps the rows of its holes
as integer, since which() returns double positions past 2^31 - 1 cells
(`R/slide_imp.R`: `hole_rows[[i]] <- hole_row_index(at_new, nrow(holes))`).
`slide_parts()` divides each sum by the number of windows that ran over its
column and builds the parts from those rows, a `left_na` one for the columns
where some window left a hole (`R/slide_imp.R`:
`slide_parts <- function(subset, hole_n, acc, n_ran, seen, rows, kind) {`),
and one `fill_parts()` writes them into a copy of the `subset` columns
(`R/slide_imp.R`: `result <- fill_parts(obj, parts, keep = result_cols,`). So
no matrix the size of `obj` is built besides the result, `obj` is read once by
the scan and once by the windows, and no observed cell is written.

`tune_imp()` with `.f = "knn_imp"` or `.f = "slide_imp"` injects the columns
of each repetition's masked cells as `subset`, and reads the estimates at each
cell's column's position in that subset (`R/tune_imp.R`:
`read_loc <- if (is_subset_mode) {`).

## Test entry points

A test calls `impute_knn_brute()` directly and asserts the returned `loc` is
`identical()` to `which(is.na(...), arr.ind = TRUE)` (`tests/testthat/test-knn_imp.R`:
`calculate the missing location correctly`); others cover `subset` by index and
by name with `post_imp` both ways (`tests/testthat/test-knn_imp.R`:
`works with post_imp = FALSE/TRUE`), all-NA rows and columns plus abort site 1
(`tests/testthat/test-knn_imp.R`: `extreme missing columns and rows`), the tie
rule (`tests/testthat/test-knn_imp.R`: `tied neighbors are settled by scan order`),
and `Inf` rejection (`tests/testthat/test-knn_imp.R`: `Throw on Inf`); and one
pins the `slideimp_infeasible` class on both abort sites
(`tests/testthat/test-knn_imp.R`: `class is set on both abort sites`). The
narrowing by `subset` and the parts are in
`tests/testthat/test-knn_imp-parts.R`; `tune_imp()` reading the narrowed result
is in `tests/testthat/test-tune_imp.R`
(`tune_imp() reads knn_imp() estimates at the masked cells`); the kernel's `fill_mean` and the mean
over observed values in `tests/testthat/test-knn_imp.R`
(`leftover_input <- function() {`, an input with both kinds of cell K-NN
leaves).
