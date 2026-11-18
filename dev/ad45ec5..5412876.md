# Review: eigen_conversion.md Phase 1 (ad45ec5..5412876)

Scope: the RcppArmadillo-to-plain-Rcpp conversion of `src/mat_stats.cpp`,
`src/matrix_checks.h`, `src/sample_each_rep_cpp.cpp`, `src/imputed_value.h`,
`src/imputed_value.cpp`, `src/impute_knn_brute.cpp` and their R callers
(`R/col_vars.R`, `R/knn_imp.R`, `R/mat_miss.R`, `R/tune_imp.R`), plus the
tests the range touches. `src/mipca.cpp` and the files in `dev/vendored.md`
were skipped; the one vendored caller of a changed export,
`R/mi_pca_imp.R`
(`R/mi_pca_imp.R`: `v <- col_vars(Z)`), passes no `cores`, so it is not
broken.

Code is cited by file and a fragment of one line, as `dev/` notes do.

## The three questions asked

### Int narrowing or overflow in index and offset arithmetic

None found. Every column offset is taken in `std::size_t`, either through
`col_ptr()`, then in `src/imputed_value.h`, or by hand (finding 6). Quantities kept as `int` are each bounded by
`nrow()` or `ncol()`:

- `GroupLayout` counts: `validate_knn_inputs()` checks the groups are in
  range and disjoint before the `static_cast<int>` of their sizes, so the
  sum is at most `ncol()`.
- `initialize_result()` sums the missing counts in `size_t` and refuses a
  total over `INT_MAX` before `no_init_matrix(n_total, 2)`. The `loc`
  matrix then holds `2 * n_total` cells, which can pass `INT_MAX`;
  `Rf_allocMatrix` takes that as a long vector, and every write uses
  `offset + r` in `size_t`.
- `sample_each_rep_cpp()` sums `na_per_col` in `size_t` and refuses a total
  over `INT_MAX`; `res[total_na + j]` is `size_t`. Written rows and columns
  are below `nrow()` and `ncol()`, so the `+ 1` cannot overflow.
- `col_miss_internal()` and `row_miss_internal()` count into `int`, bounded
  by `nrow()` and `ncol()`.

Measured, not assumed: `checkmate::assert_count()` in `sample_na_loc()`
refuses `max_attempts = 3e9` ("Must be of type 'count'") and passes
`.Machine$integer.max`, so the new `const int max_attempts` cannot receive a
value that `Rf_coerceVector` would turn into `NA`.

### Rcpp objects or the R API inside the RcppThread worker lambdas

None found. The only worker left in the range is the `par_for()` lambda in
`src/impute_knn_brute.cpp`. It captures by name: raw pointers
(`obj_masked_ptr`, `nmiss_masked_ptr`, `obj_ptr`, `grp_complete_ptr`,
`value_ptr`), scalars, `std::vector`s, the plain `GroupLayout`, and the
`unique_ptr<RcppThread::ProgressBar>`. `value_ptr` is taken from
`value.begin()` on the main thread after `initialize_result()` has assigned
`value`, and `value` stays in scope, so protected, for the whole loop; the
workers write disjoint slices of it. `distance_vector()`, `knn_weights()`
and `impute_column_values()` take pointers and std containers only. The
lambda calls `RcppThread::checkUserInterrupt()` and the progress bar's
`operator++`, both RcppThread's published thread-safe API.
`src/mat_stats.cpp` no longer has a worker.

### Order or number of `R::runif` draws against the Armadillo version

Unchanged. Read in the installed RcppArmadillo:
`RcppArmadillo/rng/Alt_R_RNG.h` defines `arma_rng_alt::randu_val()` as
`return double(::Rf_runif(0, 1));`, and `R::runif(0, 1)` is a call to the same
`Rf_runif(0, 1)`. The diff replaces the four `arma::randu()` sites (the
Fisher-Yates step over the pool, the first protected row, the second
protected row, the commit loop) one for one, at the same places, under the
same conditions, with the same `std::min(..., (size_t)(u * n))` mapping. No
draw was added, removed or moved across a `continue`. Both exports keep
`Rcpp::RNGScope` in `src/RcppExports.cpp`.

## Findings

Ranked most severe first. None is a regression introduced by the range in
the K-NN results or the sampled positions.

All eight were fixed in the commit after `5412876`, as recorded under
"Phase 1" in `dev/eigen_conversion.md`. A line quoted below that the fix
removed is quoted as it stood at `5412876`, not cited.

### 1. Progress bar never reaches 100% when a column has no neighbor (correctness, low)

In `src/impute_knn_brute.cpp`, at the comment `// no candidate shares an
observed row with this column: its`, the worker returns early when `top_k` is empty, before `++(*bar);`. The bar
is built for `layout.n_imp` steps, so every column whose neighbors all share
no observed row with it leaves the bar short of its total, and it never
prints its completion. `group_imp()` defaults to `.progress = TRUE`, so a
user sees this on an ordinary call where post-imputation then fills those
columns. Pre-existing in a lambda this range rewrote. Fix: increment the bar
before the early return, or once at the top of the body.

### 2. `sample_each_rep_cpp()` does not refuse a repeated pool column (correctness, low)

The new guards in `src/sample_each_rep_cpp.cpp`, such as
`check_range(pool_idx_in, obj.ncol(), "pool_idx_in");`, check range but not uniqueness. A repeated column in
`pool_idx_in` can fill two slots: `left_temp` is rebuilt from the unchanged
`obj` column on each visit, so the same `(row, col)` can be returned twice
and `num_na` overcounts the cells actually blanked. `sample_na_loc()`
never passes a repeat (`seq_len()`, `match()` on unique names, or
`assert_integerish(unique = TRUE)`), so this is reachable only by a direct
call. The K-NN kernel refuses its equivalent (`claim_group` in
`src/imputed_value.h`); the sampler's new test says it "checks the indices
it relies on" but does not check this one. Fix: a `seen` pass like
`claim_group`, or state in the header comment that uniqueness is the
caller's.

### 3. No test pins the sampled positions under a fixed seed (test-coverage)

The draw-order invariance this phase depends on is checked only by the
`dev/` capture and compare scripts. The suite
(`tests/testthat/test-tune_imp.R`: `test_that("set.seed makes sample_na_loc reproducible", {`)
compares two runs of the SAME build, so a later edit that reorders or adds an `R::runif` call (Phase 2
or 3, or the kernel hygiene pass in `dev/to-do.md`) still passes
`devtools::test()`. A stored integer map of positions for one small seeded
call would catch it; per `.claude/CLAUDE.md`, `identical()` is correct for
integer maps even against a stored reference.

### 4. README claims `col_vars()` uses RcppArmadillo and RcppThread (documentation)

The README still credits RcppArmadillo and RcppThread, in the `README.Rmd`
line beginning "Fast column-wise variance via". After this range `col_vars_internal()` is plain Rcpp on one thread. The
README is public documentation and now states something false. Fix in
`README.Rmd`, then re-knit.

### 5. `check_range()` duplicates `validate_knn_inputs()`'s `check_group` (reuse)

The two are `check_range()` in `src/sample_each_rep_cpp.cpp` and the
`check_group` lambda in `src/imputed_value.h`. Both loop an `IntegerVector`, refuse a value below 0 or at or above a
bound, and build the same NA-aware message. Written twice in one phase,
they will drift. One helper in `src/matrix_checks.h`, which both files
already reach, with the bound and the name as arguments.

### 6. Column pointers hand-rolled beside the existing `col_ptr()` (reuse)

The helper was `inline T *col_ptr(T *base, const int col, const int n_rows)`
in `src/imputed_value.h`; two of the hand-rolled sites were
`const double *col = x + j * n_rows;` in `src/mat_stats.cpp` and
`const double *candidate = x + static_cast<std::size_t>(col_idx) * n_rows;`
in `src/sample_each_rep_cpp.cpp`. `col_ptr()` exists to take a column offset in `size_t` in one place, yet
`src/mat_stats.cpp` (five sites), `src/matrix_checks.h` and
`src/sample_each_rep_cpp.cpp` compute `x + j * n_rows` by hand, and four
header comments each restate the "size_t, so past 2^31 cells" rule. Today
every hand-rolled site is safe only because both operands happen to be
`size_t`; a later edit that makes either one `int` reintroduces the
overflow `col_ptr()` was written to prevent. Move `col_ptr()` into
`src/matrix_checks.h` (already included by `src/imputed_value.h`), give it
a `size_t` overload, use it at every site, and state the rule once.

### 7. `check_range()`'s upper bound is vacuous for three of its four calls (simplification)

The bound is set in `src/sample_each_rep_cpp.cpp` by
`const long long int_max = std::numeric_limits<int>::max() + 1LL;`. For `na_per_col`, `row_room` and `pool_cap` the bound is `INT_MAX + 1`, so
`v[i] >= bound` can never hold for an `int`: those calls are a
non-negativity check. The variable named `int_max` holds `INT_MAX + 1`,
and the error text reads "outside 0 .. 2147483647" for a check whose
only reachable failure is a negative value or NA. Simpler: a
`check_nonnegative()` for the three, `check_range()` for `pool_idx_in`
only.

### 8. The plan says nothing is converted at the sampler boundary; `max_attempts` is (documentation, low)

The plan said "already builds as integers, so nothing is converted" in the
step 1-2 record of `dev/eigen_conversion.md`. The vectors are integers, but `max_attempts` defaults to `100`, a double,
and `assert_count()` does not convert it, so Rcpp coerces it at the call.
The conversion is harmless (measured above: the range check keeps it
inside `int`), but the note states it as absent. Name `max_attempts` as
the one scalar converted.

## Checked and clean

- `air format --check .` reports nothing.
- The ASCII check reports only `.github/workflows/pkgdown.yaml`, which the
  range does not touch.
- `uv run dev/check_map_anchors.py`: 0 problems.
- Every R caller validates `obj` as a numeric matrix before
  `check_finite()` or `check_inf()`, so a vector or data frame never reaches
  `Rcpp::NumericMatrix`'s "not a matrix" error.
- `knn_imp()` writes `obj[imputed$loc] <- imputed$value`; with no missing
  cell this is a 0 x 2 index and a length-0 value, which R accepts.
- The removed `cores` argument of `col_vars()` is recorded in the NEWS
  entry of `dev/to-do.md`, as `.claude/CLAUDE.md` asks.
- Items already deferred in `dev/to-do.md` (the zero-fill of `obj_masked`,
  `rows_to_impute` repeating `loc`, `initialize_result()`'s guards) are not
  repeated here.
