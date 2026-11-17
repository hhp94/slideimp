# Parts inside, `si_mi` at the door

The plan for the `si-mi` branch, which starts from `e5b3cf7` on `dev`. Delete
this file when the branch merges; git history keeps it.

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
  (`src/mat_stats.cpp`: `if (n == 0 || n == n_rows)`).
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
order (`src/mat_stats.cpp`: `arma::mat mean_imp_col_internal(const arma::mat &mat,`).
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
  but the Phase 0 capture. The `eigen` branch can drop it rather than port it.
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
