# One census at the front door: scoping

Status: scoped 2026-11-20, read-only; no code written, nothing timed. It
lands LAST, after the other entries of `dev/to-do.md` that touch the doors
("Parts for the other imputers", "Take cli out; audit every message against
the R manuals", "Kernel cleanup after the Eigen conversion"). The to-do entry
"One census at the front door" points here.

Goal: every imputation entry point scans its input as few times as possible
before the first real work, and does it the same way everywhere.

Everything below was confirmed by reading the code at `b1857c9`, except
where it says "inferred". Claims marked "to measure" have not been measured.

## 1. What each entry point scans today

"Full" means one pass over every cell.

- `knn_imp()`: `check_finite()` full (`R/knn_imp.R`: `check_finite(obj)`),
  then the per-column count (`R/knn_imp.R`:
  `cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)`). After the work, when
  `post_imp` has columns to fill, `mean_imp_col()` runs `check_inf()` and a
  column count over the whole matrix (`R/mean_imp_col.R`: `check_inf(obj)`,
  `n_miss <- col_miss_internal(obj)`), and `fill_parts()` counts the input's
  columns and then the output's (`R/parts.R`:
  `n_miss <- col_miss_internal(obj)`,
  `if (wrote && any(col_miss_internal(out)[out_col[fill_cols]] > 0)) {`).
- `pca_imp()`: `check_finite()`, the per-column count, `col_vars()`, and
  the per-row count under `row.w = "n_miss"` (`R/pca_imp.R`:
  `check_finite(obj)`, `cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)`,
  `obj_vars <- col_vars(obj)`,
  `n_miss_per_row <- mat_miss(obj, col = FALSE, prop = FALSE)`). After the
  work, `fill_parts()` adds the same two column counts as above.
- `group_imp()`: `check_finite()` and the per-column count over the whole
  matrix (`R/group_imp.R`: `check_finite(obj)`,
  `cmiss <- unname(mat_miss(obj, col = TRUE, prop = FALSE))`). Each group is
  then copied out (`R/group_imp.R`:
  `sub_mat <- obj[, indices[[i]]$col_idx, drop = FALSE]`,
  `sub_mat <- src[, task$col_idx, drop = FALSE]`) and goes through the whole
  door of `knn_imp()` or `pca_imp()` again. `on_infeasible = "mean"` runs
  `mean_imp_col()` on a block copy (`R/group_imp.R`:
  `obj[, f] <- suppressMessages(mean_imp_col(obj[, f, drop = FALSE]))`).
  That `attributes(obj) <-` makes a wrapper, not a copy, is inferred from
  the comment above it (`R/group_imp.R`:
  `attributes(obj) <- list(dim = dim(obj))`).
- `slide_imp()`: `check_finite()` when every column is covered, otherwise
  `check_inf()` full and then `check_finite()` on COPIED chunks of the
  covered columns (`R/slide_imp.R`: `check_finite(obj)`, `check_inf(obj)`,
  `check_finite(obj[, covered_cols[s:e], drop = FALSE])`). Each window is a
  copy (`R/slide_imp.R`: `sub_mat <- obj[, window_cols, drop = FALSE]`); under
  PCA an `anyNA()` on it or on a second copy of its targets (`R/slide_imp.R`:
  `targets_complete <- imp_method == "pca" &&`); then the imputer's door on
  the window. Windows overlap, so a column is scanned once per window that
  covers it. No "nothing to impute" exit at this level.
- `tune_imp()`: `check_finite()` (`R/tune_imp.R`: `check_finite(obj)`). With
  `na_loc = NULL`, `sample_na_loc()` then runs four more full passes, its
  `check_inf()` redundant here (`R/tune_imp.R`: `check_inf(obj)`,
  `col_miss <- unname(mat_miss(obj, col = TRUE))`,
  `row_miss <- unname(mat_miss(obj, col = FALSE))`,
  `pool_var <- unname(col_vars(obj))[pool_idx]`). A supplied `na_loc` costs
  only a gather (`R/tune_imp.R`: `on_hole <- is.na(obj[na_loc[[i]]])`). Each
  task (repetition x parameter set) copies and blanks the input
  (`R/tune_imp.R`: `pre[na_loc[[rep_id]]] <- NA`, `pre <- src[,]`) and runs
  the imputer's whole door on it.
- `sample_na_loc()` on its own: the four passes above.
- `mean_imp_col()` (internal, no longer exported): `check_inf()` and a
  column count, both full.
- The sampler (vendored, section 5): `.nothing_to_impute()` (`anyNA()` with
  an early exit; a whole-matrix `sum()` when complete; with `subset`, a
  loop that copies one column at a time) (`R/mi_pca_imp.R`:
  `.nothing_to_impute <- function(obj, subset) {`). With a fitted transform,
  `.check_admissible()` runs before the fit and again inside `mipca_core()`
  on the same matrix, each a column count plus `col_vars()`
  (`R/mi_pca_imp.R`: `.check_admissible(obj)`, `.check_admissible(Z)$cmiss`).
  `tf_sv()` builds two logical masks for its [0, 1] check (`R/mi_pca_imp.R`:
  `if (any(y < 0, na.rm = TRUE) || any(y > 1, na.rm = TRUE)) {`).
  `mipca_tune()`, `mipca_flank()`, `.abb_split()` and `abb_cols()` add
  their own (section 5).

Full passes at the door, valid double input: `knn_imp()` 2; `pca_imp()` 3,
or 4 with `row.w = "n_miss"`; `tune_imp()` with default `na_loc` 5 plus
every task's door; `group_imp()` 2 plus 2 to 4 per block copy; `slide_imp()`
1 to 2 plus a copy plus 2 to 4 per window copy.

Integer input, inferred from Rcpp's semantics, to measure: each scanner
takes `const Rcpp::NumericMatrix &`, which coerces an integer matrix, so
every scanner and each kernel makes its own full double copy of it
(`src/mat_stats.cpp`: `void check_finite(`).

## 2. The policies

Wording today, from C++ through `Rcpp::stop` (`src/mat_stats.cpp`:
`void check_finite(`):

- Inf: "'obj' holds an infinite value at row %d, column %d. Infinite values
  are not supported". Names a cell: the first in column-major order.
- All-NA column: "'obj' has a column with no observed value, only NA or
  NaN. Remove that column before imputation". Names nothing.

By entry point:

- `knn_imp()`, `pca_imp()`, `group_imp()`, `tune_imp()`: refuse both. The
  first offending column by index wins; within a column, Inf before all-NA.
  Columns over `colmax` are allowed and mean-filled.
- `pca_imp()` past the door: a column with variance below
  `.Machine$double.eps`, NA variance, or missingness over `colmax` is
  ineligible and mean-filled, not refused. A column whose variance overflows
  to Inf counts as eligible and stops later, at the cross-product error.
- `slide_imp()`: Inf anywhere, reported in the caller's coordinates; all-NA
  only in covered columns; Inf over the whole matrix is checked before any
  all-NA, so an earlier all-NA column does not win as it does in
  `knn_imp()`. Tests pin both halves and the cell (`tests/testthat/test-slide_imp.R`:
  `"row 3, column 5"`).
- `sample_na_loc()`: Inf only; all-NA columns get a cap of 0.
- `mean_imp_col()`: Inf only; all-NA columns come back untouched.
- The sampler: plain `stop()`, per column, up to 8 labels per bullet, under
  a lead line: "non-finite variance (an Inf cell, or overflow)", "no
  observed cells", "only one observed cell" (`R/mi_pca_imp.R`:
  `.admissible_census <- function(Z) {`). Constant columns are flagged, not
  refused. Its "nothing to impute" test is `!anyNA(v) && is.finite(sum(v))`
  and returns the empty fit without the admissibility check.

## 3. The candidate census, and whether it carries every policy

One internal C++ pass, in place, that counts and never refuses. Per column:

- `n_miss` (int): `n_rows` less Welford's observed count, which equals
  today's `isnan` count exactly.
- `inf_row` (int): 1-based row of the column's first Inf, 0 for none. Not a
  flag and not a count: only this rebuilds today's cell message, and it is
  recorded during the pass, so the position comes by indexing, not by a
  search afterwards.
- `var` (double): today's Welford loop unchanged, Inf entering the
  arithmetic, NA below two observed values (`src/mat_stats.cpp`:
  `Rcpp::NumericVector col_vars_internal(`). The sampler's non-finite test
  and `pca_imp()`'s eligibility read its bits.
- Optional, behind a flag: per-row missing counts, for
  `row.w = "n_miss"`, `sample_na_loc()` and the sampler's tuner.

It dispatches on `REALSXP` and `INTSXP`, so integer input is read in place;
an integer matrix cannot hold Inf, and `NA_INTEGER` is missing.

What each policy then is, in R:

- First offender for the four main doors:
  `which(inf_row > 0L | n_miss == n)[1]`, Inf message if `inf_row` is set.
- `slide_imp()`: any `inf_row > 0`, then `n_miss[covered] == n`.
- `sample_na_loc()`, `mean_imp_col()`: `inf_row` only.
- The sampler: `n_miss` and `var`, as `.admissible_census()` reads them now.

Where it needs care:

- No early exit, by design. See section 4.
- The sampler's `is.finite(sum(v))`: not carried. Its nearest form,
  `n_miss == 0 & inf_row == 0`, differs on a complete matrix whose grand
  total overflows: today that goes on to a fit, under the census it would
  return the empty fit. The current test also depends on the platform's
  long double. Upstream's decision.
- Blocks: a block's column statistics are the door census indexed by the
  block's columns (same values, same order within a column), so
  `group_imp()` groups and `slide_imp()` windows need no rescan. A block's
  ROW counts cannot be derived, so `row.w = "n_miss"` still costs a block
  pass.
- `tune_imp()` tasks: counts are derivable from the census plus `na_loc`;
  the variance of the blanked copy is not. So one census per repetition,
  shared across parameter sets, not one per task.
- Staleness: a census describes the matrix as passed. `knn_imp()`'s
  `mean_imp_col()` runs after the K-NN write, so it takes its columns by
  construction, not from the door census.
- The 0-based column-index option: no caller outside the sampler needs it.
  `slide_imp()` already reads the whole matrix through `check_inf()`, so one
  full census reads fewer bytes than today. Only the sampler's
  narrow-`subset` nothing-to-impute path would use it.
- min and max: read only by vendored code (`abb_cols()`, `tf_sv()`), and a
  transform's `forward` cannot receive a census without a change to the
  transform contract.

## 4. What changes for users

- The class and call of the two refusals, inferred, to measure: from
  Rcpp's C++ exception classes and the call `check_finite(obj)` to whatever
  R raises from the entry point. The wording can stay verbatim; the tests
  match fixed substrings. Since this lands after "Take cli out", the
  refusals are plain `stop()` by then, and the NEWS reminder line goes in
  that entry's list.
- Refusing Inf per column instead of per cell breaks the pinned cells in
  `tests/testthat/test-slide_imp.R` and changes the error surface; keeping
  `inf_row` avoids both. The help says how to locate the cells
  (`R/slideimp-package.R`: `@section Missing values and non-finite input:`).
- No input newly refused or newly accepted if the policies are copied
  exactly. The one risk is the sampler's overflow edge in section 3.
- Early exits lost: a refused input reads one full pass before it errors,
  where today it stops at the first bad column (slower failure, never
  slower success); and `mipca(obj, subset = few)` with complete subset
  columns would read the whole matrix where it reads only those columns.
  The "no missing values" exits of `knn_imp()` and `pca_imp()` already sit
  after full passes and are unaffected.
- `dev/knn_imp_structure.md` and `dev/pca_imp_structure.md` cite
  `check_finite()`; fix them in the same change and run
  `uv run dev/check_map_anchors.py`.

## 5. The sampler (vendored: upstream, not here)

What it would read, if upstream adopts the census. It would then call a
slideimp internal, so `hhp94/test` carries the same function.

- `.admissible_census()`: `n_miss`, `var`.
- `.mipca_run()` and `mipca_core()`: one census passed into the core, where
  two run today.
- `.nothing_to_impute()`: `n_miss`, `inf_row` over `subset`, after the
  overflow decision.
- `.abb_split()`: `n_miss`, `inf_row` over the listed columns, dropping its
  block copy and mask (`R/mi_pca_abb.R`:
  `nonfin = unname(colSums(is.infinite(blk)) > 0)`).
- `.flank_census()`: the census, and `inf_row` for its `abb` targets,
  dropping a copy (`R/mi_pca_flank.R`:
  `bad[tg_l] <- colSums(is.infinite(obj[, tg_l, drop = FALSE])) > 0`).
- `abb_cols()`: `n_miss`, and min and max if carried (`R/mi_pca_abb.R`:
  `r <- range(obj, na.rm = TRUE)`).
- `mipca_tune()`: `n_miss`, `var`, row counts (`R/mi_pca_tune.R`:
  `cen <- .check_admissible(obj)`, `row_na = .mat_miss_int(X, col = FALSE),`).
  Its `which(is.na(obj))` stays: it needs positions (`R/mi_pca_tune.R`:
  `na0 <- which(is.na(obj))`).

This overlaps the sampler-port item "The `Inf` refusal has two shapes" in
`dev/to-do.md`; settle them together.

## 6. Recommendation, and what is the maintainer's

Recommended:

- One internal C++ census as in section 3. `check_finite()` and
  `check_inf()` go. `mat_miss()` and `col_vars()` stay exported; whether
  `col_vars()` reads the census's `var` instead of its own loop is decided
  by a measurement.
- One R helper per policy, building today's messages from the census.
- `knn_imp()` and `pca_imp()` each split into the exported function and an
  internal core that takes the census; `group_imp()` and `slide_imp()` pass
  the core an indexed census per block, `tune_imp()` one per repetition.
  The exported signatures do not change.

Rejected:

- `check_finite()` returning counts while still refusing in C++: the
  policies differ (covered columns only, all-NA allowed), so policy flags
  would move into C++.
- R-level `is.finite()` or `colSums()` masks: they allocate n x p.
- A census argument on the exported functions: an API change.
- The column-index option and min/max now: nothing outside the sampler
  reads them. Add them when upstream asks.

The maintainer's, with a position on each:

1. One pass per entry point at the cost of early exit: yes. The early exit
   saves only on refused input; the sampler's narrow-subset path keeps its
   cheap pre-check upstream.
2. min and max in the census: not yet.
3. Refusals per column or per cell: per cell for Inf through `inf_row`;
   all-NA per column, which could now name the column, as a deliberate
   wording change.
4. The sampler's overflow semantics, and whether upstream adopts the census
   at all: with upstream.

To measure before building (`dev/`, by hand, minimal runs):

- The census against today's passes, on `dev/GSE286313.EPICv1_chr21.qs2`
  and on `sim_mat()`.
- Whether the added Inf branch slows the Welford loop.
- Whether the integer coercion copy happens.
- The class and call of today's refusals.
