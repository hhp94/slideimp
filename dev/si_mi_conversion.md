# Parts inside, `si_mi` at the door

The plan for the `si-mi` branch, which starts from `e5b3cf7` on `dev`. Delete
this file when the branch merges; git history keeps it.

The Armadillo removal (`dev/eigen_conversion.md`) is done on this branch too
(maintainer, 2026-11-18). Where this file says "the `eigen` branch", it
means that work; no such branch is merged.

## Why

- One internal currency for every imputer and every fallback. An imputer
  returns a PART: the holes it filled and their values. A fallback is another
  part, placed back by column, not a write into the matrix. The end goal is
  `group_imp()` building its result from its groups' parts; this branch does
  `mean_imp_col()` and `pca_imp()` only.
- `pca_imp()` gains `output = "mipca"`: its point estimate as a one-draw
  `si_mi`, the object `mipca()` returns, so that both read through
  `as.matrix()` and `summary()`.
- Memory. Readings on one machine (Windows, process peak working set), one
  2000 x 20000 input (307 MB) with 10% of columns over `colmax`, one run per
  stage, each in a fresh R process (`dev/probe_pca_imp_peak.R`). "Extra" is
  the peak above the input:

  | Stage, at `e5b3cf7` | Extra |
  |---|---|
  | `pca_imp()`, `post_imp = TRUE` | 3.06 x |
  | `pca_imp()`, `post_imp = FALSE` (the kernel and the write-back) | 2.42 x |
  | The mean-fill step alone (write-back, then `mean_imp_col()`) | 2.99 x |
  | Parts, filled into one copy of the input | 1.65 x |
  | Parts sliced into an `si_mi` | 1.88 x |
  | Parts sliced, then `as.matrix()` | 2.92 x |

  The three part rows include the probe's own hole-sized vectors (about
  0.5 x at 11.75% holes). So when any column is mean-filled, the peak of
  `pca_imp()` is the mean-fill step, not the kernel; filling one copy of the
  input from parts leaves the kernel as the peak. Slicing and then calling
  `as.matrix()` gains nothing, which is why the matrix output never goes
  through an `si_mi`.

## What was decided (maintainer, 2026-11-15)

- **`si_mi` is built at the door only.** Inside the package, functions pass
  parts. A part holds the columns of the input it covers and the holes it
  filled, both in the INPUT's numbering, one value per hole, a `kind`
  (`"pca"`, `"mean"`), an `info` list, and a `left_na` flag. It holds no copy
  of the input.
- **Two doors.** The matrix door copies the input once and writes every
  part's values at its holes. The `si_mi` door slices the input once per part
  and builds the object with the vendored constructors.
- **`pca_imp()` takes `output = c("matrix", "mipca")`**, default `"matrix"`.
  `slide_imp()`, `group_imp()` and `tune_imp()` pass `output = "matrix"`
  explicitly. The value names the functions the result works with
  (`mipca_complete()`, `cbind()` with `mipca()` fits), not the class:
  `si_mi` is an internal name and users never type it (maintainer,
  2026-11-17). A character argument, not a logical flag, so a later form,
  or the same door on another imputer, takes the same argument.
- **`output = "mipca"` refuses an input without unique row and column
  names**, by the rule `mipca()` applies (`R/mi_pca_imp.R`:
  `.chk_matrix <- function(obj) {`): the join matches columns and rows by
  name, so names `pca_imp()` made up would join two unnamed fits by
  position, unchecked. The matrix output keeps taking unnamed input; parts
  address cells by position and need no names (maintainer, 2026-11-17).
- **The `slideimp_results` class and its print method stay** for the matrix
  output. The print is a later plan.
- **The mean fill follows the ABB pattern.** The columns `pca_imp()` leaves
  out are one mean part, placed back by column, as `mipca()` places its
  ABB-listed columns. No cell-level fallback in `pca_imp()`: a column the
  kernel returns any `NA` for moves to the mean part whole.
- **`mean_imp_col()` is no longer exported.** CRAN listed no reverse
  dependency of slideimp 1.2.0, in any field.
- **A point estimate and a sampler fit are not joined** only through the
  kept-draw check, so only when the sampler fit holds more than one draw:
  `cbind()` refuses fits with different numbers of kept draws. A one-draw
  point fit joins a one-draw sampler or ABB fit over other columns, and that
  is allowed (ruled in `hhp94/test`, see "Upstream").
- `knn_imp()`, `slide_imp()` and `group_imp()` are out of scope, except that
  they keep returning what they return today.
- **`subset` narrows the output, package-wide** (maintainer, 2026-11-17). A
  call given `subset` returns only those columns, as `mipca(subset = )`
  returns a fit over the requested columns only (`R/mi_pca_s3.R`:
  `.rebase_subset <- function(x, fr) {`). No flag: `subset = NULL` keeps
  every column, so only a caller who passes `subset` sees a change, and the
  full width back is one assignment on their side. The columns come back in
  the input's order, as `mipca()` returns them, whatever order `subset` names
  them in. The record's column indices (`mean_filled_cols`, the `col` of
  `mean_filled_cells`) count in the output's numbering. This branch does it in
  the doors (Phase 5) and in `pca_imp()` (Phase 6). For `knn_imp()` it is a
  breaking change to callers that pass `subset` today (`slide_imp()`,
  `group_imp()`, `tune_imp()`), so it lands with their change, which is
  delayed (`dev/to-do.md`, "Parts for the other imputers").

## Order of work (maintainer, 2026-11-15)

1. The "Upstream" items are investigated and ruled in `hhp94/test` first.
   Done 2026-11-15; see "Upstream".
2. This branch then runs Phase 0, Phase 1 and Phase 2 up to the point where
   Phase 2 would touch `src/mat_stats.cpp`.
3. From there the `eigen` branch runs alongside, since it rewrites that file.
4. Phase 4's door builds its parts with `.si_mi_chunk()`, which the port of
   `hhp94/test@40eae1d` brought, with the new `print.si_mi()`. Done
   2026-11-16, after Phase 3: only `R/mi_pca_s3.R`, `R/mi_pca_abb.R` and
   `tests/testthat/test-mipca-narrow-join.R` changed since the port of
   `hhp94/test@14a2b05` in `e5b3cf7`.

## Rules that apply to every step

- Read `dev/to-do.md` first, and `dev/pca_imp_structure.md` before changing
  `pca_imp()`. A step that changes cited code fixes the map in the same
  commit, and `uv run dev/check_map_anchors.py` passes before it is
  committed.
- The vendored files (`dev/vendored.md`) are not edited. What they need is
  listed under "Upstream" and made in `hhp94/test`.
- One slideimp-owned file is the only caller of the vendored internals the
  `si_mi` door uses (`R/mi_pca_s3.R`: `.si_mi_new <- function(`,
  `.si_mi_join <- function(parts, call, placement = NULL) {`;
  `R/mi_pca_s3.R`: `.diag_frame <- function(map, draws, p, nms, prop_diag) {`;
  `R/mi_pca_s3.R`: `.si_mi_chunk <- function(`), so a port that changes them
  breaks in one place.
- `output = "matrix"` returns what it returns at `e5b3cf7`: values, type,
  dims, class and attributes, checked against the Phase 0 capture, not
  remembered. Phase 3 found the one intended difference it planned, the
  kernel-`NA` case, unreachable, so there is none.
- No change under `src/` unless a step says so: the `eigen` branch rewrites
  `src/mat_stats.cpp` and `src/armaSVD.cpp` (`dev/eigen_conversion.md`).
- Tests follow `.claude/CLAUDE.md`: worst-case absolute and relative error
  against a stored reference, on a tolerance, with bitwise-ness reported;
  `identical()` only for two objects built in one process and for maps; every
  subset a test asserts over is asserted non-empty first.
- One commit per step. Before each: `air format .`, `devtools::document()`
  if roxygen moved, the suite with `load_all1(timer = FALSE)` and
  `check_manual()`, and the Phase 0 compare script.

## Phase 0: a reference

Before any source edit, two scripts in `dev/`, tracked:

- A capture script, run at `e5b3cf7`, saving each result with its type,
  dims and attributes to an ignored folder under `dev/` (regenerable: rerun
  at `e5b3cf7`). It covers `pca_imp()` with `post_imp` on and off, `clamp`,
  `nb.init > 1`, an input with no holes, and inputs with columns over
  `colmax` and constant columns; `mean_imp_col()` with and without `subset`;
  and `knn_imp()`, `slide_imp()` (both methods, `on_infeasible = "mean"`) and
  `group_imp()`, which call `mean_imp_col()` and so move under Phase 2. On
  `sim_mat()` output and on `dev/GSE286313.EPICv1_chr21.qs2`, under fixed
  seeds.
- A compare script, as `dev/eigen_conversion.md` describes its own. The two
  branches capture the same functions, so one pair of scripts can serve both;
  whichever branch lands them first, the other reuses them.

`dev/probe_pca_imp_peak.R` becomes tracked here and is rerun after Phases 2,
3 and 4.

Done 2026-11-15, as `dev/ref_capture.R` and `dev/ref_compare.R`; each
header says how to run it. The capture covers the union of this phase's list
and `dev/eigen_conversion.md`'s, so the `eigen` branch takes the same pair.
What it found on the reference build:

- 58 cases, about 1.5 s. Every case reaches the branch its name claims:
  `mean_filled_cols` non-empty where columns are over `colmax`, the
  `on_infeasible = "mean"` fallback taken in `slide_imp()` and `group_imp()`,
  the `colmax` boundary column, the LOBPCG solver chosen when asked.
- The imputers refuse an all-`NA` column, while `mean_imp_col()` leaves it
  `NA`. So the capture uses chr21 with its 1797 all-`NA` probes dropped
  (71 x 8503, 16 probes over `colmax`) for the imputers, and keeps them for
  `mean_imp_col()`; the refusal itself is two cases.
- Noise floor: two captures in two processes are bitwise identical, every
  case. The default tolerances (abs `1e-12`, rel `1e-10`) are set above that
  floor, and `dev/ref_compare_selftest.R` shows the compare still fails on
  each kind of planted change, a few-ulp change included at tolerance 0.

The reference capture is `dev/ref/e5b3cf7`, made with `load_all1(timer =
FALSE)`. It is regenerable, so it is not tracked: rerun the capture script
on any commit whose package source is `e5b3cf7`'s (its `_meta.rds` says
which source it ran on).

## Phase 1: the part and the matrix door

A new file of slideimp's own, holding:

- The part constructor and its validator. A map row lies inside the input;
  its column is one of the part's columns; the linear indices are strictly
  increasing (column-major); there is one value per hole; no value is `NA`.
  A `left_na` part has an empty map.
- The matrix door: copy the input once, and for each part assert that every
  cell it writes is still `NA` (a value over an observed cell is a bug) before
  writing. Parts cover disjoint columns. Every hole of the input is written,
  except in `left_na` parts.

Tested on its own: a part over no columns, one column, every column; parts
out of column order; a hole written twice and a write over an observed cell,
both refused.

Done 2026-11-15, as `R/parts.R` (`new_part()`, `validate_part()`,
`fill_parts()`) and `tests/testthat/test-parts.R`. Two readings of the text
above:

- "Every hole of the input is written" is checked per part: a part that is
  not `left_na` fills exactly the holes of its own columns. Columns no part
  covers pass through with their holes. Phase 2 needs this, since
  `mean_imp_col(subset = )` is one mean part over the subset and leaves
  every other column as it is. A caller that must fill every column, as
  `pca_imp()` must, passes parts that cover every column; the `si_mi` door
  will check full cover itself, since placement needs it.
- A failed check is a package defect, so it is a plain `stop()` with the
  fixed prefix "slideimp internal, part:" (`dev/WRITING.md`, section 1),
  not cli.

Measured with `tracemem()`: two writing parts copy the input once; a
`left_na` part, or no parts, copy nothing. NOT MEASURED: that
`col_miss_internal()`, which counts the holes, reads the input without a
copy; RcppArmadillo maps a `const arma::mat &` argument onto R's memory, and
the memory probe after Phase 2 would show a copy.

## Phase 2: `mean_imp_col()` on parts

- A mean part: the mean of the observed values of each requested column, read
  from the input one column at a time with no slice, as the value of each of
  its holes. An all-`NA` column stays `NA`, as today's kernel leaves it
  (`mean_imp_col_internal()`, deleted in Phase 1 of
  `dev/eigen_conversion.md`; `git show 76d4c60:src/mat_stats.cpp`).
- `mean_imp_col()` becomes the matrix door over one mean part, with the same
  refusals and messages (`check_inf()`, the subset resolution) and the same
  result.
- The export goes: drop `@export`, mark the page `@noRd`, take it out of the
  helpers listed in `README.Rmd`. Add a line to the NEWS entry in
  `dev/to-do.md`. `slideimp_paper` calls it, in
  `03.MCAR_MAR_MNAR.qmd`; that repository needs `slideimp:::` or its own
  column mean.
- `knn_imp()`, `slide_imp()` and `group_imp()` call it unchanged; the compare
  script checks them against Phase 0.
- Rerun the memory probe: the mean-fill step should fall from about 3 x to
  about 1 x the input.

Open: how the means are computed. Today's kernel sums in double, in row
order (`mean_imp_col_internal()`, since deleted; see above).
R's `sum()` and `mean()` accumulate in long double, so means taken in R can
differ in the last bits, which the compare script would report and the
tolerance would pass. Keeping them bitwise needs a C++ entry that returns
only the means, which touches `src/mat_stats.cpp`, the file the `eigen`
branch rewrites. Measure the R version's difference and its speed on the
chr21 data before choosing. The same measurement settles whether `cores` is
still worth having.

Measured 2026-11-15 (`dev/probe_mean_source.R`), on chr21 with its all-`NA`
probes kept (71 x 10300, 461 columns with a hole and an observed value) and
on a 2000 x 20000 normal matrix with column offsets of sd 100 (11000 such
columns). Each mean against an exact one (pieces on fixed power-of-two
grids, summed exactly; within about one ulp of the true mean):

| Mean | chr21: max ulps, columns > 1 ulp | 2000 x 20000: max ulps, columns > 1 ulp |
|---|---|---|
| today's kernel, double sum | 5, 131 | 335, 8363 (median 4) |
| `.colMeans()`, long double sum and divide | 1, 0 | 1, 0 |
| `sum(na.rm = TRUE) / n` | 1, 0 | 2, 1 |
| `mean()` of the observed values | 1, 0 | 1, 0 |

So the R means are not bitwise the kernel's, and they are the more accurate
ones; against the kernel the worst is rel `6.7e-16` on chr21 and `4.0e-14` on
the large matrix, well inside the compare's `1e-10`. `cores` bought at most
1.8 x (on the large matrix over every column: 307, 198 and 171 ms at 1, 4
and 8 cores), and only on the full-matrix copy the kernel wrote, which the
parts path does not make.

Done 2026-11-15, with the means from `.colMeans()` per column, as
`mean_part()` in `R/mean_imp_col.R`, and `cores` dropped from
`mean_imp_col()` and from its two callers that passed it.

Speed, `bench::mark` medians, `mean_imp_col()` at 1.2.0 (`cores = 1`)
against on parts: chr21, 4.5 ms against 6.1; the large matrix over every
column, 307 ms against 479; over its last 10% of columns (95% holes, 3.8
million of them), 151 ms against 317. Slower by up to 2 x, in a step that is
a fraction of a second beside the imputer that calls it.

Memory, `dev/probe_mean_fill_peak.R`, on the input of
`dev/probe_pca_imp_peak.R` with its last 10% of columns filled (3.8 million
holes), extra peak working set over the input, fresh process per stage:

| Stage | Extra peak | R allocations over the call |
|---|---|---|
| the 1.2.0 kernel | 2.00 x | 1.00 x (plus its Armadillo output, outside R) |
| `mean_imp_col()` on parts | 1.65 to 1.81 x | 1.96 x |
| `mean_part()` alone | 0.52 to 0.62 x | 0.63 x |
| `fill_parts()` alone, part prebuilt | 1.33 x | 1.33 x |

So the fill falls from 2 x to about 1.7 x, not to the 1 x hoped for above.
Past the door's one copy: the part itself (an integer map and a double value
per hole, 16 bytes a hole, 0.2 x here); the validator's vectors the length
of the map (its two columns and the linear index); and what reading one
column at a time leaves behind (each column, its `is.na()`, the hole rows),
which R collects late, so it shows in the peak. The peak moves by up to
0.2 x between runs of the same code, with the timing of R's collections.
What would take it near 1 x is C++ that counts the holes, takes the means
and fills the map in one pass over the input, and writes the copy: `src/`,
so for the `eigen` branch or after it.

Two changes to the door came out of this:

- `validate_part()` checks the map with `min()`, `max()`, `tabulate()` and
  `is.unsorted()` on a double linear index, not with logical vectors of the
  map's length, and returns that index; `fill_parts()` writes with it.
- `fill_parts()` finds a write over an observed cell by counting the holes
  left in the filled columns after the write, a pass over the result,
  instead of reading the input at every cell of the map before it. With the
  counts checked before the write and the map free of repeats, a hole left
  open means some cell of the map was observed; only then does it look for
  which one, for the message.

Readings and choices:

- The choice of R over a C++ entry is open to the maintainer until this
  branch merges: R is more accurate and leaves `src/` alone; C++ keeps every
  stored result bitwise. Switching is the one line in `fallback_part()` that
  takes the mean, plus the entry in `src/`.
- A column with no observed value is left out of the part by
  `mean_imp_col()`, so the door passes it through as it is. `fallback_part()`
  itself does not take one: its mean would be NaN, and the validator refuses
  it.
- The result keeps today's shape: a double matrix with dim and dimnames and
  no other attribute, integer input coerced. The attribute strip comes after
  the fill, so it edits the door's copy in place.
- Revised 2026-11-16, at the maintainer's call: the per-column `mean_part()`
  became `fallback_part(obj, cols, method = c("mean", "median"))`,
  vectorized over its columns. It copies their submatrix out once with its
  `is.na()`, takes the map from `which(arr.ind = TRUE)`, and sets `kind` to
  the method; `part_kinds` gained `"median"`. A copy of the fallback's
  columns is small beside the input, since most columns are filled by PCA
  or K-NN. No caller takes the median yet; exposing it is a separate,
  public choice. Measured (`dev/probe_fallback_part.R`), on chr21 and on the
  2000 x 20000 matrix over every column and over its last 10%: the mean
  part's map and values are `identical()` to the per-column part's; the
  median, `collapse::fmedian()`, is bitwise `stats::median()` on every
  column (chr21 has 272 columns with an even observed count and 8231 odd,
  the large matrix even only); speed is the same (2.2 ms against 2.2 on
  chr21, 81 ms against 94 over the last 10%), and the median part takes 6.2
  and 101 ms against 121 and 75 for `apply(stats::median)`. R allocations
  over the last 10% rise from 0.63 x to 0.96 x the input, and over every
  column of the large matrix to 2.93 x, the case the fallback does not meet.
  The peak (`dev/probe_mean_fill_peak.R`, same input as the table above):
  `mean_imp_col()` 1.58 x, `fallback_part()` alone 0.67 x, the door alone
  unchanged at 1.33 x.
- `col_miss_internal()` reads the input without a copy, which Phase 1 left
  unmeasured: `fill_parts()` alone peaks at 1.33 x, the same as R's own
  allocations over the call, so its two calls allocated nothing of matrix
  size outside R.
- `mean_imp_col_internal()` stays in `src/mat_stats.cpp`, with no R caller
  but the Phase 0 capture. The `eigen` branch can drop it rather than port it;
  Phase 1 of `dev/eigen_conversion.md` did.
- `slideimp_paper` still calls `mean_imp_col()`; not changed from here.

Checked: the suite with `load_all1(timer = FALSE)` and the manual tests,
20 files, 247 tests, none failed or skipped; the `group_imp` and `tune_imp`
tests again with daemons on an install of this checkout
(`dev/test_mirai_with_dev_install.R`), 118 tests, none failed: with the
installed 1.2.0 copy their daemons would have run the old mean fill. The
Phase 0 compare passes 57 of 58 cases: `mean_xa_cores2` now errors on the
removed `cores`, as intended; six mean-fill cases are not bitwise, worst rel
`6.7e-16`; types, dims, classes and attributes match on every case.

The memory probe (`dev/probe_pca_imp_peak.R`) after Phase 2, against the
table under "Why": `pca_imp()` with `post_imp = TRUE` 2.65 x (3.06 x
before), the mean-fill step alone 2.64 x (2.99 x before), `post_imp = FALSE`
unchanged at 2.42 x. The mean-fill step still holds the write-back's copy;
Phase 3 removes it.

## Phase 3: `pca_imp()` on parts, matrix output

- The kernel's `imputed_values` already is a map and its values: row, column
  in the input's numbering, value. READ: the kernel visits the changing
  columns in the order `arma::find()` gives, ascending, through
  `eligible_idx_perm` (`src/armaSVD.cpp`: `arma::uvec eligible_idx_perm = eligible_idx.elem(perm);`).
  NOT ASSERTED anywhere: the part's validator checks it, and a failure is an
  internal error. No kernel change.
- `clamp` is applied to the values before the part is built.
- The parts:
  - the PCA part: the eligible columns and their holes;
  - the mean part: every other column, plus any eligible column the kernel
    returned an `NA` for. Under `post_imp = FALSE` it is a `left_na` part.
- The write-back through a two-column index and the `mean_imp_col()` call
  after it are replaced by the matrix door.
- `mean_filled_cols` is read off the mean part's columns that have holes.
  `mean_filled_cells` stays on the result and is always empty for
  `pca_imp()`.
- ~~Intended difference: a kernel-`NA` cell moves its whole column to the
  mean part.~~ Dropped, see below: the kernel-`NA` branch cannot be reached.
- The early exit for an input with no holes is unchanged.
- `dev/pca_imp_structure.md` is updated: the write-back it cites goes.

Done 2026-11-16, in `R/pca_imp.R`, with no change to `R/parts.R` or `src/`.

The kernel-`NA` branch was dead code, so it was dropped rather than
changed. READ, not measured, since finite input does not make the kernel
return NaN: the kernel's final SSE sums the squared residual over every cell
of each eligible column, holes included, and only then subtracts the holes'
share (`src/armaSVD.cpp`: `sse -= s * ej2;`). So a NaN or Inf at any imputed
cell makes `mse` NaN, and at 1.2.0 `if (cur_obj < init_obj)` then failed
with base R's "missing value where TRUE/FALSE needed" before any value was
written. In its place, each restart's `mse` must be finite or the call stops
as a package defect, a plain `stop()` with the prefix "slideimp internal,
pca_imp:" (`dev/WRITING.md`, section 1). The part's validator still refuses
an `NA` value, as a second net. No NEWS line: no input reaches either.

As built:

- The PCA part covers `eligible_idx`, including eligible columns with no
  hole. The mean part covers every ineligible column, including the
  zero-variance ones with no hole, so the two parts cover every column, as
  the Phase 4 door needs. `fallback_part()` copies those columns out, which
  is small beside the input.
- `check_finite()` refuses an all-`NA` column before the kernel, so every
  ineligible column has a mean, as `fallback_part()` requires.
- An integer input is coerced to double before the door, as the write into
  it used to coerce it. That is a second copy, for integer input only.
- `best_imputed` is dropped once the part is built.

Tests: the post-step regression test now mocks `fallback_part()` and asserts
it is asked for exactly the ineligible columns, and for none on an
all-eligible input. A new test stubs the kernel with NaN, NA, Inf and -Inf
`mse`, on the first and on a later restart, and with a NaN value under a
finite `mse`, which the validator refuses. With the guard disabled the
`mse` case fails on base R's "missing value" error, so the test can fail.

Checked: the suite with `load_all1(timer = FALSE)` and the manual tests, 20
files, 248 tests, none failed; its two `group_imp` errors were daemons
running an older installed copy ("could not find function new_part"), and
the `group_imp` and `tune_imp` tests on an install of this checkout
(`dev/test_mirai_with_dev_install.R`) pass, 118 tests. The Phase 0 compare
gives what Phase 2 gave: 57 of 58 cases, `mean_xa_cores2` erroring on the
removed `cores`; every `pca_*` case bitwise, `post_imp` on and off, clamp,
`nb.init`, the no-holes exit and chr21 with its 16 columns over `colmax`;
the same six mean-fill cases not bitwise, worst rel `6.7e-16`.

The memory probe (`dev/probe_pca_imp_peak.R`) after Phase 3: `pca_imp()`
with `post_imp = TRUE` 2.43 x (2.65 x after Phase 2, 3.06 x at `e5b3cf7`),
the same as `post_imp = FALSE`, 2.43 x. So the mean fill no longer sets the
peak; the kernel does, as "Why" expected of filling one copy from parts.

## Phase 4: the `si_mi` door, `output = "mipca"`

- `output = c("matrix", "mipca")`, checked with `match.arg()`. It goes
  after `clamp` and before `.progress`, the last argument (maintainer,
  2026-11-17).
- Under `output = "mipca"`, an input without unique row and column names is
  a user error, through cli, raised before the fit runs, so no fit is
  wasted on it.
- The door, in the file of Phase 1 or beside it, builds each part's chunk
  with `.si_mi_chunk()`, never a bare `list()`: per part, slice the input to
  the part's columns (no slice when one part covers every column: the
  chunk's `obj` is the input itself, a reference); rebase the map to the
  slice's columns, as an INTEGER matrix with columns `row` and `col`
  (`.si_mi_chunk()` refuses a double map); one draw
  (`matrix(values, ncol = 1)`); an identity `backward`; the empty diagnostic
  frame, as an ABB piece takes; `prov` with `chains = 1`, `ndraws = 1`,
  `thin = 1`, `thin_effective = 1`, `seam_state = "never"`, and the rank and
  iteration count for the PCA part; and `info`, passed as `.si_mi_chunk()`'s
  argument, holding the same fields, under the same names, as the matrix
  output's attributes. `info` lives on the chunk, not on the object:
  `.si_mi_new()` rebuilds the object from its own fields on every join of
  two or more fits, so a field on the object would be dropped. Then
  `.si_mi_new()` per part and `.si_mi_join()` with the parts' columns,
  concatenated, as placement.
- A `left_na` part sets `prov$left_na`, with an empty map and 0 x 1 draws.
- An input with no holes: one part, empty map.
- What `.si_mi_chunk()` at `hhp94/test@40eae1d` requires, from the reply
  upstream (checked there, not here):
  - `prov` needs only `seam_state` and `thin_effective`; every other field
    is passed through, `thin_final` included, and no reader of a finalized
    fit uses `chains`, `ndraws`, `thin`, the rank or the iteration count.
    Upstream names the last two `rank` and `iter`.
  - `nkeep` to `.si_mi_new()` equals the draw columns of every chunk, 1
    here; `1` and `1L` both pass and both join a sampler or ABB fit.
  - `placement` is INTEGER and uses every input column exactly once, so the
    parts partition the columns, as Phase 3's two parts do.
  - `call` is required; `NULL` passes.
  - The input must have column names, or the door stops with an internal
    error. Settled: `output = "mipca"` refuses such input first (see "What
    was decided"), so the internal error stays unreachable from the door.
- Tests: `as.matrix()` of the `mipca` output is `identical()` to the matrix
  output in the same process, on a named input; an input missing row names,
  column names, or holding a repeated one, is refused under
  `output = "mipca"` and imputed under `output = "matrix"`; columns in input order; `post_imp = FALSE`
  leaves the same cells `NA`; `mipca_complete(m = 1)` works and `m = 2`
  refuses; `summary()` has no cells; `cbind()` with a `mipca()` fit of 100
  draws refuses.
- Roxygen, following `dev/WRITING.md`. It adds `@param output` and
  `@returns` for both forms, and says what one draw means:
  `mipca_complete()` takes `m = 1` only, `summary()` has nothing to show,
  and `cbind()` joins only fits of one draw, a one-draw sampler or ABB fit
  included.
- Rerun the memory probe.

Done 2026-11-17, in `R/parts.R` (`parts_to_si_mi()`, `check_si_mi_names()`),
`R/pca_imp.R`, and one line each in `R/slide_imp.R`, `R/group_imp.R` and
`R/tune_imp.R`; no change under `src/` or to a vendored file. Settled with
the maintainer before the code:

- Under `post_imp = FALSE`, the ineligible columns are two parts, not one:
  a `left_na` part over those with a hole and a mean part with an empty map
  over those without. `print.si_mi()` counts a `left_na` chunk's columns as
  "left NA", so a zero-variance column with no hole would have been counted.
  The matrix output is unchanged: `fill_parts()` writes nothing for either.
- Every chunk carries the call's whole record as `info`: the 14 attributes
  of the matrix output after dim, dimnames and class, under their names and
  types, `mean_filled_cols` included as integers in the input's numbering.
  So the record survives a `cbind()` that separates the chunks, though its
  column indices then still count this fit's columns, not the joined fit's
  (see "Review of Phases 1 to 6").
- The stored call is `.scrub_call(match.call(), "pca_imp")`, so a matrix
  passed by value through `do.call()` is not stored in the fit.
- `knn_imp()` goes the same way later, and `group_imp()` will build from
  parts (`dev/to-do.md`, "Parts for the other imputers"), so the door and the
  name check take any parts and say nothing about PCA. The imputer's own
  `prov` fields (the PCA part's `rank` and `iter`) reach the door as
  `prov_extra`, one list per part.

As built, and where the text above turned out short or wrong:

- A part over no columns must be dropped, which the list above does not
  say: `.si_mi_new()` refuses a chunk with no columns, and Phase 3 builds an
  empty mean part whenever every column is eligible. The door drops it, so
  an all-eligible input is one chunk holding the input itself.
- The door calls three vendored functions the rules above do not list:
  `.scrub_call()`, `.prop_diag_default()` (the frame's `prop_diag`, so
  `summary()` agrees it with a `mipca()` partner) and the exported
  `tf_identity()`, whose `backward` is the identity map. `R/parts.R` is
  still their one caller outside the vendored files.
- `prov` also holds `nmiss` and `thin_final = 1`, as an ABB piece's does.
- The door makes `fill_parts()`'s checks (shared as `check_parts()`), and
  two of its own: the parts cover every column, and no map cell is observed.
  The last is checked at birth, reading the input at every cell of the maps,
  a vector the length of the maps, so a defect shows where it was made
  rather than at the first `as.matrix()`, which asserts it again.
- The door takes an integer input, which the no-holes exit passes through
  as it is. `as.matrix()` of that fit is double, where the matrix output
  stays integer: `as.matrix.si_mi()` assigns the (empty) double draws into
  the stored matrix, and R coerces on that assignment. Every input with a
  hole is coerced to double before either door, so this is the one case
  where `as.matrix()` of the fit is not `identical()` to the matrix output.
- The name check mirrors `.chk_matrix()`: a missing, NA or empty name, or a
  repeated one, is refused, and the test asserts `.chk_matrix()` refuses
  each case too. The error names the imputer, not the helper.
- Vendored wording a `pca_imp()` user now meets, reported here and not
  changed: `print.si_mi()` heads its output "si_mi fit", the class name the
  decisions above call internal; `mipca_complete(m = 2)` says "the fit's 1
  kept draws"; `update()` says "Fit again with mipca() and
  `keep_state = TRUE`". The record in `info` has no exported reader: it is
  `fit$chunks[[i]]$info`.

Tests: `tests/testthat/test-pca_imp-mipca.R` (not under the `test-mipca-*`
prefix, which belongs to the vendored copy) covers the list above, plus the
no-holes exit, the all-eligible single chunk, the record on every chunk, the
stored call without the matrix, the refusal reached before the kernel (a
stubbed kernel that errors shows the stub is in place), and a `cbind()` with
a one-draw `mipca()` fit that joins. `tests/testthat/test-parts.R` adds the
door on its own: interleaved and out-of-order parts against `fill_parts()`,
one part over every column, a part over no columns, a `left_na` part, and
each refusal.

Checked: the suite with `load_all1(timer = FALSE)` and the manual tests, 21
files, 260 tests, 2593 expectations, none failed; its two `group_imp` errors
were daemons running an older installed copy ("could not find function
new_part"), and the `group_imp` and `tune_imp` tests on an install of this
checkout (`dev/test_mirai_with_dev_install.R`) pass, 118 tests. The Phase 0
compare gives what Phase 3 gave: 57 of 58 cases, `mean_xa_cores2` erroring on
the removed `cores`; every `pca_*` case bitwise; the same six mean-fill cases
not bitwise, worst rel `6.65e-16`.

The memory probe (`dev/probe_pca_imp_peak.R`, two stages added), extra peak
over the input, one run per stage:

| Stage | Extra peak |
|---|---|
| `pca_imp()`, `post_imp = TRUE` | 2.42 x (2.43 x after Phase 3) |
| `pca_imp()`, `post_imp = FALSE` | 2.43 x (2.43 x) |
| `pca_imp(output = "mipca")` | 2.43 x |
| the same, then `as.matrix()` | 2.78 x |
| by hand: slices into an `si_mi` | 1.89 x (1.88 x under "Why") |
| by hand: slices, then `as.matrix()` | 3.20 x (2.92 x) |

So the `si_mi` door does not move the peak: its slices (0.9 x here, plus the
maps) sit under the kernel's. `as.matrix()` of the fit adds the completed
matrix beside the fit's slices. The hand-built rows move by up to 0.3 x
between this run and the one under "Why", with the timing of R's
collections, as Phase 2 saw.

## Phase 5: narrowing in the doors

The rule is under "What was decided". Narrowing lives in the doors, once, so
every imputer on parts gets it from one place.

- `fill_parts()` and `parts_to_si_mi()` take `keep`: the columns of the input
  to return, ascending, in the input's numbering. `NULL` means every column
  and is today's behaviour, unchanged.
- The parts must cover every column in `keep`. A part may also cover columns
  outside it: the door drops those columns and their holes, so an imputer
  whose kernel returns holes outside `keep` (the PCA kernel, Phase 6) hands
  its part over whole. An imputer still builds nothing it can skip: a mean
  part covers only kept columns, since each mean costs a pass over its
  column.
- The door cuts the parts to `keep` before its checks, so the checks
  (`check_parts()`, the observed-cell check) read only kept columns.
- The matrix door copies the kept columns once (`obj[, keep, drop = FALSE]`,
  or the whole input when `keep` is every column), rebases the maps through
  one inverse table, as `.rebase_subset()` does, and writes the values.
- The `si_mi` door slices each part to its kept columns; placement is their
  positions in `keep`. The no-slice case (the chunk's `obj` is the input
  itself) applies only when one part covers every column and `keep` is every
  column.
- A part left with no kept column is dropped, as a part over no columns is
  today.
- Tests, in `tests/testthat/test-parts.R`: each door with `keep` is
  `identical()`, in one process, to the same door without it followed by
  `[, keep]` (for the `si_mi` door, through `as.matrix()`): interleaved parts,
  a part crossing the edge of `keep`, a part wholly outside it, a `left_na`
  part, `keep` of one column; `keep` of every column is `identical()` to
  `NULL`. Refused: `keep` unsorted, repeated, out of range, or not covered
  by the parts.
- No imputer passes `keep` yet, so the suite and the Phase 0 compare do not
  move.

Done 2026-11-17, in `R/parts.R` only: `check_keep()`, `cut_parts()`, and a
`keep` argument on `check_parts()`, `fill_parts()` and `parts_to_si_mi()`.
No change under `src/`, to a vendored file, or to any imputer.

As built, and where the text above turned out short or wrong:

- The matrix door keeps its own contract: a kept column that no part covers
  comes back as it is, holes included, as a column no part covers does
  without `keep`. Only the `si_mi` door needs the parts to cover every kept
  column. The list above put that need on both doors.
- The checks do not all read only kept columns, as the list above says.
  `validate_part()` and the two-parts check run on the parts as given, before
  the cut, since the cut relies on a valid part. The hole count compares
  each part with its kept columns only, but `col_miss_internal()` takes no
  column subset, so it still reads every column of the input (no copy), and
  `src/mat_stats.cpp` is not touched on this branch. The observed-cell check
  reads only kept cells. A defect in a column outside `keep` is not looked
  for, since nothing of that column is returned.
- A `keep` that names every column is turned into `NULL` on entry, so it
  takes today's path exactly: no copy before the first write, and the
  no-slice case of the `si_mi` door.
- A `keep` that is empty, double, `NA`, out of range, unsorted or repeated
  is a package defect (`part_defect()`), not a user error: Phase 6 resolves
  and refuses a bad `subset` before any door sees it.
- With `keep`, the matrix door's result is `obj[, keep, drop = FALSE]`, so
  it carries the input's dimnames for the kept columns and no other
  attribute; without it, the input's attributes are untouched, as before.
- The `si_mi` door's placement is `match()` of the parts' kept columns in
  `keep`, the fit's own numbering.

Tests, in `tests/testthat/test-parts.R`: both doors against the full fill cut
to `keep`, with the parts given in either order, over five `keep`s (each part
crossing the edge, interleaved, one part wholly inside and the other wholly
outside, one column, only complete columns); `keep` of every column; a
`left_na` part cut; a kept column no part covers (matrix door); the hole
count refusing within `keep` and ignoring a cut-away column; an observed-cell
write found in the narrowed copy and reported in the input's numbering; the
interleaved fit's placement, slice, rebased map and draws; a part wholly
outside `keep` dropped; one part narrowed to a slice; the coverage refusal;
and each malformed `keep` in both doors.

Checked: the suite with `load_all1(timer = FALSE)` and the manual tests, 21
files, 264 tests, 2657 expectations, none failed; its two `group_imp` errors
are again daemons running an older installed copy, and the `group_imp` and
`tune_imp` tests on an install of this checkout pass, 118 tests. The Phase 0
compare gives what Phase 4 gave: 57 of 58 cases, `mean_xa_cores2` erroring on
the removed `cores`; every `pca_*` case bitwise; the same six cases not
bitwise, worst rel `6.65e-16`. No memory probe: no imputer passes `keep`, so
no path the probe measures changed.

## Phase 6: `pca_imp(subset = )`

- `subset`, resolved to ascending column indices from names or positions.
  It goes after `output` and before `.progress`, as Phase 4 placed `output`,
  so no argument but `.progress` moves.
- An unknown name and an empty `subset` are refused, as `mipca()` refuses
  them (`R/mi_pca_imp.R`:
  `.resolve_cols <- function(subset, obj, what = "subset") {`), through cli
  (maintainer, 2026-11-17). `resolve_subset()` instead drops an unknown name
  with a message and turns an empty selection into an early return
  (`R/utils.R`: `resolve_subset <- function(subset, obj, sort = FALSE) {`);
  once the output's width depends on `subset`, a dropped name is a column
  the caller asked for and does not get. `knn_imp()` follows in its own
  change (`dev/to-do.md`), so `resolve_subset()` stays as it is here.
- The kernel still runs on every eligible column. Every hole in every
  eligible column is updated at each EM iteration and feeds the next fit, so
  a column outside `subset` still shapes the result, and leaving it out
  would change it: `subset` is a memory parameter here, not a compute one,
  as at `mipca()`. The kernel's triplets outside `subset` reach the PCA part
  and the door drops them (Phase 5); no change under `src/`. Returning only
  the kept triplets from C++ waits on a measurement: in the probe's design
  they are about 0.07 x the input (9000 holed eligible columns of 100 holes,
  24 bytes a hole, against 307 MB), which is arithmetic from the probe's
  set-up, not a reading.
- The early exit: no hole in a `subset` column returns the kept columns
  without running the kernel, as `knn_imp()` does. Under `"mipca"` that is
  one part over them with an empty map.
- The infeasibility check is restated over `subset`: no eligible `subset`
  column with a hole is `slideimp_infeasible`, as in `knn_imp()` ("All
  subset columns with missing values exceed colmax"). With `subset = NULL`
  the check and its message are today's.
- The mean part covers only the ineligible `subset` columns; under
  `post_imp = FALSE` Phase 4's split into a `left_na` part and an empty mean
  part applies to those.
- `mean_filled_cols` and `mean_filled_cells` count in the output's numbering,
  on the matrix output and in every chunk's `info`. With `subset = NULL`
  that is the input's numbering, so nothing changes.
- Under `"mipca"` the name check reads the whole input, as `mipca()`'s does.
- `subset` is not added to the parameter registry (`R/utils.R`:
  `.param_registry <- list(`), so `tune_imp()`'s `parameters` and
  `group$parameters` cannot set it, as with `output`. No caller passes it
  yet: their change is delayed with `knn_imp()`'s.
- Tests: with `subset`, the matrix output's values and dimnames are
  `identical()` to the full output's `[, subset]` in one process, since the
  same fit ran, and `as.matrix()` of the `"mipca"` output is `identical()` to
  the narrowed matrix output; `subset` of every column gives the full
  output; names and positions give the same result; an unknown name and an
  empty `subset` are refused; an unsorted `subset` comes back in input
  order; the early exit never reaches the kernel (a
  stubbed kernel that errors); the infeasible refusal; the record's indices
  in the output's numbering, on an input where they differ from the input's
  (asserted non-empty).
- Roxygen, following `dev/WRITING.md`. It adds `@param subset` in its own
  words, not inherited from `knn_imp()`, whose `subset` keeps the full width
  until its own change. It says what is returned, in what order, and that
  every eligible column still takes part in the fit.
- The memory probe gains a stage with `subset` over a tenth of the columns,
  for both outputs.
- `dev/to-do.md` gains the NEWS line: `pca_imp()` gains `subset`.

Done 2026-11-17, in `R/pca_imp.R`, `R/utils.R` (`resolve_narrow_subset()`) and
one roxygen line in `R/slideimp-package.R`; no change under `src/`, to a
vendored file, or to any caller.

As built, and where the text above turned out short or wrong:

- The refusals are one new helper beside `resolve_subset()`, written for
  `knn_imp()` to take in its own change: cli errors for an empty `subset`, a
  name `obj` does not have, a name given twice, and names on an input
  without column names; checkmate (`.var.name = "subset"`) for numbers out
  of range, repeated, `NA`, not whole, or not numbers at all. It returns
  `NULL` for a selection of every column, so `subset = 12:1` on 12 columns
  is the call without `subset`, record and message included.
- The early exit keeps its message when `subset` is `NULL` and says "No
  missing values in the `subset` columns. Returning them unchanged." when it
  is not. The infeasibility message names `subset` in the same way; its fix
  line is unchanged.
- The matrix output's early exit is `obj[, keep, drop = FALSE]`, not a
  door, as the full early exit returns `obj` itself: an integer input stays
  integer there, as before.
- An integer input with holes is still coerced to double whole before the
  doors, so with `subset` it is copied whole once more than a double input.
  Not changed: the coercion is Phase 3's, and integer input is rare.
- The "Column mean fill" section, shared by every imputer's help page, now
  says `mean_filled_cols` holds the index "in the result". That is true of
  every imputer, since only `pca_imp()` narrows and the others' results
  number as their inputs do; `knn_imp.Rd`, `slide_imp.Rd`, `group_imp.Rd`
  and `slideimp-package.Rd` change with it.
- The memory probe's first subset (every tenth column from column 1) picked
  only odd columns, which hold no hole in its design, so the only holes in
  `subset` were in columns over `colmax` and `pca_imp()` refused the call as
  infeasible, as the rule above says it should. The stage takes two columns
  in every twenty instead, one with holes and one without.

Tests: `tests/testthat/test-pca_imp-subset.R`. With `subset` (sorted,
unsorted, names, doubles, one column) and both `post_imp`, the matrix output
is `identical()` to the full output's selected columns, and its record to the
full record but `mean_filled_cols`; a `subset` of every column, in any order,
is the call without it; `mean_filled_cols` counts the result's columns (6 in
the input, 2 in the result); `as.matrix()` of the `"mipca"` output is the
narrowed matrix output, with the record on every chunk and the column left
NA; a `subset` without holes never reaches the kernel (a stubbed kernel that
errors, which a `subset` with a hole shows is in place); the infeasible
refusal and its class; each refusal names `pca_imp()`; `tune_imp()` refuses
`subset` in `parameters`.

Checked: the suite with `load_all1(timer = FALSE)` and the manual tests, 22
files, 272 tests, 2738 expectations, none failed; its two `group_imp` errors
are again daemons running an older installed copy, and the `group_imp` and
`tune_imp` tests on an install of this checkout pass, 118 tests. The Phase 0
compare gives what Phases 4 and 5 gave: 57 of 58 cases, `mean_xa_cores2`
erroring on the removed `cores`; every `pca_*` case bitwise; the same six
cases not bitwise, worst rel `6.65e-16`.

The memory probe, extra peak over the input, one run per stage:

| Stage | Extra peak |
|---|---|
| `pca_imp()`, `post_imp = TRUE` | 2.44 x (2.42 x in Phase 4) |
| `pca_imp()`, `post_imp = FALSE` | 2.43 x (2.43 x) |
| `pca_imp(output = "mipca")` | 2.44 x (2.43 x) |
| the same, then `as.matrix()` | 2.79 x (2.78 x) |
| `pca_imp(subset = )`, a tenth of the columns | 2.43 x |
| the same, `output = "mipca"` | 2.43 x |

So `subset` does not move the peak of one call, as expected: the kernel sets
it, and the kernel still runs on every eligible column. What `subset` saves
is the result a call returns and keeps, a tenth of the input here instead of
a whole copy, which lowers the peak of a caller that holds results while it
works. No caller passes `subset` yet.

## Review of Phases 1 to 6 (2026-11-17)

A code review in a fresh session, on `e5b3cf7..3c84803` over `R/` and
`tests/` without the vendored files, set up as the uncommitted diff of a
throwaway worktree. It ran no code. Its findings were then audited here, and
those that can be measured were measured with `dev/audit_pca_imp_cr.R`,
numbered as below. What the audit found, and what was done:

1. **Fixed.** On the one-chunk path of the `si_mi` door (one part over every
   column and no `subset`, and the no-holes exit), the chunk held the input
   itself, class and attributes included, and `as.matrix.si_mi()` starts
   from that matrix there. So re-imputing a `knn_imp()` result returned a
   `slideimp_results` with `imp_method = "knn"`, and a plain matrix kept its
   own attributes, while a fit of two chunks, built from slices, kept
   neither. `part_chunk_fit()` now cuts the input to its dims and dimnames
   when it carries more (`R/parts.R`: `if (length(setdiff(names(attributes(sub)), c("dim", "dimnames")))) {`),
   which copies it only then. Tests: at the door in `test-parts.R`, and in
   `test-pca_imp-mipca.R` over a `slideimp_results` with holes, one without
   (the no-holes exit) and a plain matrix with an attribute.
2. **Wrong, not changed.** The claim was that `fallback_part()`'s map turns
   double past 2^31 cells, since `which()` returns a double index there.
   R 4.6.2's `arrayInd()` ends with `storage.mode(ind) <- "integer"`, and a
   row or column number is always below 2^31, so the map stays integer.
3. **Recorded, not changed.** The record on each chunk keeps
   `mean_filled_cols` in the fit's own numbering, and a `cbind()` does not
   renumber it, so after a join it points at the joined fit's other
   columns. Nothing reads it yet (Phase 4, "As built"). The indices stay, as
   decided for the output's numbering; the Phase 4 text and the comment in
   `pca_imp()` now say they are not renumbered.
4. **Fixed.** The mean fallback took complete columns too, and
   `fallback_part()` copied them, with their `is.na()`, to fill nothing.
   `mean_imp_col()` now hands it only the columns with a hole (`R/mean_imp_col.R`:
   `cols <- cols[n_miss[cols] > 0L & n_miss[cols] < nrow(obj)]`), and
   `pca_imp()` builds three parts in place of Phase 3's two: the ineligible
   columns with a hole, filled or `left_na`, and those without, a part with
   nothing to fill, as `post_imp = FALSE` already did. On 2000 x 5000 with a
   tenth of the columns holed, the extra R-heap peak of a fallback over every
   column drops from 2.52 x the input to 1.26 x, the result identical. The
   `slide_imp()` and `group_imp()` fallbacks, which pass every target column
   of a window or a group, gain the most. Tests: a new one in
   `test-mean_imp_col.R`, and the post-step test in `test-pca_imp.R` now
   asserts the fallback is asked for columns 1 and 2, not 3, which its title
   already said.
5. **Real, not changed.** The doors count holes again (`check_parts()`, and
   `fill_parts()` after the write): a pass is 0.010 s on 2000 x 5000, two of
   them 1.9 % of a 10-iteration `pca_imp()`, and 10 us on a 200 x 100
   window. Also real: an integer input is copied once more than at
   `e5b3cf7`, its extra R-heap peak 6.08 x the double size against 5.07 x
   there, while a double input is unchanged (1.08 x and 1.07 x). Most of
   that 5 x predates this work: each Rcpp call converts its own copy of an
   integer matrix. One conversion right after the no-holes exit would remove
   both; not done, as outside the review.
6. **Fixed.** Three comments still said `pca_imp()` has no `subset`: in its
   `row.w = "n_miss"` branch, and in the complete-targets checks of
   `slide_imp()` and `group_imp()`, which now say the callers do not pass it
   yet.
7. **Deferred to the `knn_imp()` change.** `prov_extra` is matched to the
   parts by position and checked only by length. `dev/to-do.md` carries it.

Every new test fails with its fix undone: 10 failures over the three new
tests, none elsewhere.

Checked: the suite with `load_all1(timer = FALSE)` and the manual tests, 22
files, 274 tests, 2774 expectations, none failed; its two `group_imp`
errors are daemons running an older installed copy ("could not find
function resolve_narrow_subset"), and the `group_imp` and `tune_imp` tests
on an install of this checkout pass, 118 tests. The Phase 0 compare gives
what Phase 6 gave: 57 of 58 cases, `mean_xa_cores2` erroring on the removed
`cores`; every `pca_*` case bitwise; the same six cases not bitwise, worst
rel `6.65e-16`.

## Phase 7: `knn_imp()` on parts, with its callers (plan)

The change `dev/to-do.md` delays ("Parts for the other imputers"):
`knn_imp()` builds parts, gains `output`, and its `subset` narrows its
output. `tune_imp()`, `slide_imp()` and `group_imp()` change with it, since
the narrowing breaks them. One change (maintainer, 2026-11-17), built in the
steps at the end of this section.

Decided in discussion (maintainer, 2026-11-18):

- Under `post_imp = TRUE`, a cell K-NN leaves is filled with the mean of its
  column's OBSERVED values, as a fallback part takes it and as `pca_imp()`
  and ABB fill. Today `knn_imp()` writes the K-NN values first and runs
  `mean_imp_col()` on the result (`R/knn_imp.R`:
  `obj <- mean_imp_col(obj, subset = na_cols)`), so that mean includes the
  values K-NN imputed in the same column, as the help says (`R/slideimp-package.R`:
  `found no value. For such a cell, the column mean includes the values that`).
  Those cells change value. The help sentence changes with them, and
  `dev/to-do.md` gains the NEWS line.
- The kernel fills those cells itself. It leaves a cell in two ways: a whole
  column, when no candidate shares an observed row with it
  (`src/impute_knn_brute.cpp`: `if (n_neighbors == 0)`), and a single cell,
  when none of the chosen neighbors is observed at that row
  (`src/imputed_value.cpp`: `(wt_ptr[r] > 0.0) ? (ws_ptr[r] / wt_ptr[r]) : arma::datum::nan;`).
  Both kinds sit in an imputed column, whose masked copy (holes set to 0) and
  observed count the kernel already holds (`src/impute_knn_brute.cpp`:
  `n_col_valid(local_pos) = valid_count;`), so the mean costs one pass over
  that column's copy, and nothing of the input's size.
- The columns over `colmax` stay in R. The kernel never sees them: they are
  in none of its three groups (imputed; eligible with holes outside
  `subset`; eligible and complete), so filling them there would mean a
  fourth group passed only to take means. They are a whole-column fallback,
  as in `pca_imp()`: a mean part through `fallback_part()` under `post_imp`,
  a `left_na` part with an empty map without it.

What the call builds, by column and `post_imp`, over the kept columns:

| Column | `post_imp = TRUE` | `post_imp = FALSE` |
|---|---|---|
| imputed, every hole given a K-NN value | the K-NN part | the K-NN part |
| imputed, some cells left by K-NN | the same K-NN part, those cells holding the kernel's mean | a second K-NN part, `left_na` in the widened sense below |
| over `colmax` | a mean part (`fallback_part()`) | a `left_na` part, empty map |
| complete | a part with nothing to fill | the same |

Only one cell of that table needs `left_na` widened: K-NN's leftover cells
under `post_imp = FALSE`. The second K-NN part exists for `print()`, which
counts every column of a `left_na` chunk; one `left_na` K-NN part would count
the columns whose holes were all filled.

### The kernel

- `impute_knn_brute()` gains `fill_mean`, R's `post_imp`. When it is TRUE, a
  cell it would leave NaN gets its column's observed mean: the sum of the
  masked column divided by `n_col_valid`, in double. It is taken by the
  worker that owns the column, once, and only for a column that needs it.
  `check_finite()` refuses a column with no observed value, so the count is
  never 0.
- Double, not long double (maintainer, 2026-11-18). The R fallback keeps
  `.colMeans()`, which sums and divides in long double (`R/mean_imp_col.R`:
  `mean = .colMeans(sub, nrow(sub), ncol(sub), na.rm = TRUE),`), so the two
  can differ in the last bits. They never fill the same column: the kernel
  fills cells of imputed columns, R whole columns over `colmax`.
- It also returns which result rows hold a mean, so R can record
  `mean_filled_cells`. Proposed: a byte per result row, written by the worker
  that owns the row (no race), collected after the loop into the 1-based
  rows, beside the triplet matrix in a list. A fourth double column would
  cost 8 bytes a hole where the byte costs 1.
- With `fill_mean = FALSE` the leftover cells stay NaN, as today. R's pass
  turning NaN into NA (`R/knn_imp.R`: `imputed_values[is.nan(imputed_values)] <- NA_real_`)
  is replaced by dropping those rows from the map.
- `dev/CPP.md` applies: plain loops, no OpenMP, nothing pinned. The map
  (`dev/knn_imp_structure.md`, layers 4 to 6 and "Feasibility and abort
  paths") follows the code in the same commit.
- The `eigen` branch rewrites these files (`dev/eigen_conversion.md`, Phase
  1, the K-NN chain) and checks that the kernel's return is unchanged from
  its start commit. This step changes that return, so the two plans meet
  here.

### Parts (`R/parts.R`)

- `"knn"` joins `part_kinds`.
- `left_na` widens to: a hole in the part's columns that its map does not
  name stays NA. So a `left_na` part may hold a map (`validate_part()` now
  refuses one: `if (left_na && n_map > 0L) {`), and today's `left_na` parts,
  with empty maps, are the special case. Every column of a `left_na` part
  keeps at least one hole: per column, the map lists fewer cells than the
  column has holes, checked with one `tabulate()` of the map's columns.
- The checks that no value lands on an observed cell extend to a `left_na`
  part with a map. In `fill_parts()`, by count after the write: the holes left
  in its columns are the holes before, less its map rows. In the `si_mi`
  door, `refuse_observed_cells()` runs over every part with a map.
- Finding 7 of the review: a part gains a seventh field, `prov`, the named
  list of fields its chunk records. `parts_to_si_mi()` loses `prov_extra`,
  and `pca_imp()` puts `rank` and `iter` on its PCA part.
- The `si_mi` door needs nothing new for the widened `left_na`: a `left_na`
  part with a map becomes a chunk with draws for its map and `prov$left_na`,
  and `as.matrix.si_mi()` writes only map cells, so a hole outside the map
  stays NA (`R/mi_pca_s3.R`: `as.matrix.si_mi <- function(x, draw = 1L, ...) {`).
  But `print.si_mi()` counts every column of a `left_na` chunk
  (`R/mi_pca_s3.R`: `function(ch) if (isTRUE(ch$prov$left_na)) ncol(ch$obj) else 0L,`)
  and calls them "not imputed", which a K-NN column that keeps a hole is
  not. That wording is vendored; see "To settle".

### `knn_imp()`

- `output = c("matrix", "mipca")`, after `dist_pow` and before `.progress`,
  as `pca_imp()` placed it before its last argument; `subset` stays where it
  is. Under `"mipca"`, `check_si_mi_names()` runs first, and the call is
  stored through `.scrub_call(match.call(), "knn_imp")`. Neither argument
  enters the parameter registry.
- `subset` resolves through `resolve_narrow_subset()`, so an unknown name
  and an empty `subset` are refused where `resolve_subset()` drops the one
  with a message and returns early on the other (`R/knn_imp.R`:
  `subset <- resolve_subset(subset, obj)`). The doors take it as `keep`. The
  K-NN part lies inside `keep`, since the kernel imputes only `subset`
  columns; the mean part covers only kept columns over `colmax`.
- Early exit, no hole in a kept column: the matrix output is `obj`, or
  `obj[, keep, drop = FALSE]` with `subset`, an integer input staying
  integer, as in `pca_imp()`; the `"mipca"` output is one part over the kept
  columns with an empty map. Its message is today's without `subset`, and
  says the `subset` columns come back unchanged with it, since "Returning
  input unchanged" is no longer true there. The two infeasibility messages
  are unchanged.
- An integer input with holes is coerced to double once, after the early
  exit, and both the kernel and the doors read that copy. Today the kernel
  converts its own copy and the write into `obj` converts another. That the
  kernel then takes the double input without a copy is an assumption, to
  measure.
- The record is built once for both outputs, as `pca_imp_info()` does it,
  and rides on every part as `info` under `"mipca"`. `mean_filled_cols` and
  `mean_filled_cells` count in the output's numbering, the cells from the
  kernel's mean rows.
- Roxygen, following `dev/WRITING.md`. It adds `@param output`, and
  `@param subset` saying what is returned and in what order, `@returns` for
  both forms, and an example. The "Column mean fill" section loses the sentence above and says
  the mean is over the observed values, for every imputer.

### The callers

- `tune_imp()`, `.f = "knn_imp"` and `"slide_imp"`: the result holds only
  the `subset` columns tune_imp injects, so the estimates are read at
  `match()` of each masked cell's column in that `subset`, in both paths.
- `slide_imp()`: each window passes `subset = local_idx` to `knn_imp()` and
  to `pca_imp()` alike, and gets the narrowed result back. The PCA
  complete-targets check goes (`R/slide_imp.R`:
  `targets_complete <- imp_method == "pca" &&`), as `pca_imp()`'s early exit
  now makes it, and the mean fallback returns the window's targets only, so
  every branch hands back the same width. The accumulation reads the window's
  result whole (`R/slide_imp.R`: `imputed_window[, local_idx, drop = FALSE]`)
  and the record is read in its numbering (`R/utils.R`:
  `local_mean_filled <- function(res, keep, global) {`).
- `group_imp()`: PCA groups pass `subset` too, and the workaround goes
  (`R/group_imp.R`: `needs_imp <- vapply(`). A group's features are the first
  columns of its slice (`R/group_imp.R`:
  `p$subset <- indices[[i]]$features_idx_local`), so the narrowed result is
  exactly those columns, in order. The reads take the result whole rather
  than lean on that layout.
- `slide_imp(subset = )` and `group_imp(subset = )` narrow their own output
  by the same rule. `slide_imp()` resolves through `resolve_narrow_subset()`
  (`R/slide_imp.R`: `subset <- resolve_subset(subset, obj, sort = TRUE)`),
  so its empty-`subset` exit becomes a refusal, and it allocates its result
  over the `subset` columns only. `group_imp()` already refuses a name not in
  `obj`.
- `dev/to-do.md`: the NEWS entry gains the **Breaking** line for the
  narrowing and the refusals, and a line for the changed mean.

### Tests

- The kernel: under `fill_mean = TRUE`, each cell it leaves under `FALSE`
  holds `.colMeans()` of its column's observed values, worst-case absolute and
  relative error asserted on a tolerance, bitwise reported (the kernel sums
  in double, `.colMeans()` in long double); its mean rows are
  `identical()` to the NaN rows under `FALSE`. The input has both kinds of
  leftover (a column sharing no observed row with any candidate, and a cell
  where every chosen neighbor is missing), each asserted present.
- The changed mean, pinned: on a column where the mean over observed values
  and the mean over observed and imputed values differ (asserted), the cell
  holds the first.
- Parts: a `left_na` part with a map through both doors; refused, a `left_na`
  part with a column that keeps no hole, and a map cell over an observed
  value in a `left_na` part, in each door; `prov` on the part reaching its
  chunk.
- `knn_imp()`: as Phases 4 and 6 tested `pca_imp()`. With `subset`, the
  matrix output is `identical()` to the full output's `[, subset]` in one
  process (the kernel imputes each column from the same candidates either
  way); `as.matrix()` of the `"mipca"` output is `identical()` to the matrix
  output; the refusals; the record in the output's numbering on an input where
  it differs from the input's; under `post_imp = FALSE`, the leftover cells
  NA in both outputs, with the second K-NN part's chunk flagged.
- The callers: `slide_imp()` and `group_imp()` with `subset`, narrowed,
  `identical()` to the full output's columns; `tune_imp()` scores the same
  cells as before.
- The Phase 0 compare: every case that has no cell K-NN leaves is expected
  unchanged, and the narrowed cases are compared against the reference's
  columns, which `dev/ref_compare.R` learns to do. A case with leftover cells
  changes there and only there; the compare reports which.
- A memory probe of `knn_imp()` before and after, since the write moves from
  `obj[imp_indices] <- ` into `fill_parts()`.

### Steps

Each passes the suite and the Phase 0 compare (differences as stated), with
`air format .`, `devtools::document()` and the citation check.

1. Parts: `"knn"`, the widened `left_na`, the `prov` field. `pca_imp()` moves
   to `prov`. No output moves.
2. The kernel: `fill_mean` and the mean rows. `knn_imp()` takes the kernel's
   means, and only the columns over `colmax` still go through
   `mean_imp_col()`. Only the leftover cells move.
3. `knn_imp()` on parts and `output = "mipca"`, still full width. The matrix
   output is `identical()` to step 2's.
4. The narrowing: `knn_imp(subset = )` and the callers that it breaks.
5. `slide_imp(subset = )` and `group_imp(subset = )` narrow their own output.

### To settle before the code

- Settled (maintainer, 2026-11-18): `eigen` goes first. Step 1 may go
  ahead; steps 2 to 5 wait until `eigen`'s K-NN chain
  (`dev/eigen_conversion.md`, Phase 1) is committed on this branch, so
  `fill_mean` is written in plain Rcpp.
- `print.si_mi()` says "N columns not imputed" of a `left_na` chunk. For a
  K-NN column that keeps a hole, ask upstream for wording true of both cases
  (columns that keep holes), made in `hhp94/test` and arriving with the next
  copy; or accept the line until then.
- Median stays out: no imputer offers it, and `fallback_part()`'s
  `"median"` has no caller, so a median would be a new user-facing option,
  decided on its own.
- Whether the five steps land as five commits or are squashed into the one
  change the ruling names.

## Upstream, in `hhp94/test`

Investigated and ruled in `hhp94/test` on 2026-11-15. The rulings are in
that repository's vendor to-do (dev/vendor-to-do.md there) and, with the
reasoning, its decision record (dev/DECISIONS.md there, "One chunk
constructor; the join stays unguarded"). What was done there reached this
repository with the port of `hhp94/test@40eae1d` (see "Order of work").

1. **The join guard: not done, by ruling.** A one-draw point fit may join a
   one-draw sampler or ABB fit over other columns. No point marker goes in
   `prov`; nothing reads one. The kept-draw check still refuses a point fit
   against any fit of more than one draw.
2. **`print.si_mi()` with one draw: done in `hhp94/test@7b59766`.** It prints
   "1 draw" and "next: as.matrix() for the completed matrix", offers
   `summary()` only when some frame has cells, and its "left NA" line no
   longer names the flank plan. It reads nothing from `info`.
3. **`summary()` over zero cells: not done, by ruling.** It keeps its table
   of `NA` quantiles under a header that says 0 cells.
4. **The chunk contract: done in `hhp94/test@7b59766`, as `.si_mi_chunk()`**
   in `R/mi_pca_s3.R`. It checks a chunk's fields once, at birth, passes
   `prov` whole, and stores `info` on the chunk when given.

## Not settled yet

- Whether `mean_imp_col()`'s means stay in R, as Phase 2 left them, or move
  to a C++ entry that keeps them bitwise the 1.2.0 kernel's (Phase 2).
- The parts file is `R/parts.R` unless the maintainer renames it.
- The print method, for both outputs: a later plan.
- The kernel's own peak, about 2.4 x the input above it. A separate question,
  not this branch.
