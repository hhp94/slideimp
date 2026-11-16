# Parts inside, `si_mi` at the door

The plan for the `si-mi` branch, which starts from `e5b3cf7` on `dev`. Delete
this file when the branch merges; git history keeps it.

## Why

- One internal currency for every imputer and every fallback. An imputer
  returns a PART: the holes it filled and their values. A fallback is another
  part, placed back by column, not a write into the matrix. The end goal is
  `group_imp()` building its result from its groups' parts; this branch does
  `mean_imp_col()` and `pca_imp()` only.
- `pca_imp()` gains `output = "si_mi"`: its point estimate as a one-draw
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
- **`pca_imp()` takes `output = c("matrix", "si_mi")`**, default `"matrix"`.
  `slide_imp()`, `group_imp()` and `tune_imp()` pass `output = "matrix"`
  explicitly.
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
4. Phase 4's door builds its parts with `.si_mi_chunk()`, which the vendored
   copy does not have: it predates `hhp94/test@7b59766`. Phase 4 waits for
   the next port, which brings it and the new `print.si_mi()`. When that
   port runs is the maintainer's call. Phases 0 to 3 do not depend on it.

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
  `R/mi_pca_s3.R`: `.diag_frame <- function(map, draws, p, nms, prop_diag) {`),
  and `.si_mi_chunk()` once the next port brings it to `R/mi_pca_s3.R`, so a
  port that changes them breaks in one place.
- `output = "matrix"` returns what it returns at `e5b3cf7`: values, type,
  dims, class and attributes, checked against the Phase 0 capture, not
  remembered. The one intended difference is the kernel-`NA` case in Phase 3.
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
- The write-back and the `mean_imp_col()` call
  (`R/pca_imp.R`: `obj[imp_indices] <- best_imputed[, 3]`) are replaced by
  the matrix door.
- `mean_filled_cols` is read off the mean part's columns that have holes.
  `mean_filled_cells` stays on the result and is always empty for
  `pca_imp()`.
- Intended difference, and the NEWS entry in `dev/to-do.md` gets a line for
  it: today a kernel-`NA` cell is mean-filled alone and reported in
  `mean_filled_cells`; after this, its whole column is mean-filled and
  reported in `mean_filled_cols`. A test pins the new behaviour; reaching
  the branch needs a kernel that returns `NA`, so the test may have to stub
  the kernel, as the existing post-step regression test stubs
  `mean_imp_col()`.
- The early exit for an input with no holes is unchanged.
- `dev/pca_imp_structure.md` is updated: the write-back it cites goes.

## Phase 4: the `si_mi` door, `output = "si_mi"`

- `output = c("matrix", "si_mi")`, checked with `match.arg()`. Where it
  goes in the signature is open (before `.progress`, or last).
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
- A `left_na` part sets `prov$left_na`.
- An input with no holes: one part, empty map.
- Tests: `as.matrix()` of the `si_mi` output is `identical()` to the matrix
  output in the same process; columns in input order; `post_imp = FALSE`
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
constructor; the join stays unguarded"). What
was done there reaches this repository with the next port (see "Order of
work").

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

- How `mean_imp_col()`'s means are computed (Phase 2).
- Where `output` sits in the signature, and the name of the parts file.
- The print method, for both outputs: a later plan.
- The kernel's own peak, about 2.4 x the input above it. A separate question,
  not this branch.
