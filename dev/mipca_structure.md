# mipca() structural map

A structural map of the multiple-imputation sampler, from the exported R
functions down to the C++ data-augmentation loop. Descriptive only: what the
code is and how it connects. Code is cited by file and by a fragment of one
line that occurs exactly once in that file, never by line number, e.g.
(`R/mi_pca_imp.R`: `mipca_core <- function(`). `uv run
dev/check_map_anchors.py` checks every citation; its header gives the forms.

The sampler is `R/mi_pca_*.R`, `src/mipca*.cpp` and
`tests/testthat/test-mipca-*.R`, as `dev/DECISIONS.md` ("Sampler") records;
read that section before changing anything here, since its decisions and
constraints bind the code this map describes.

Orientation: rows are samples and columns are features, as in `pca_imp()`.
The model is a rank-`S` PCA of the whole matrix: every column takes part in
every fit, and every hole is redrawn at every iteration. `subset` and `store`
choose only which holes' draws are KEPT. The sampler is the
data-augmentation scheme of Audigier, Husson and Josse (BayesMIPCA), cited on
the help page (`R/mi_pca_imp.R`: `#' @references Audigier V, Husson F, Josse J (2016).`):
an I-step draws the holes given the current low-rank mean and noise level, a
P-step refits the PCA on the completed matrix and draws a new mean for the
holes.

## Call chain overview

```
R/mi_pca_imp.R      mipca()        -- the door; seed_given = TRUE
R/mi_pca_tune.R     mipca_fit()    -- the same door, settings from a si_tune table
R/mi_pca_tune.R     mipca_tune()   -- .mipca_run() per rank, keep_state = TRUE, store = holdout + diag
R/mi_pca_flank.R    mipca_flank()  -- .mipca_run() per window, dqset.seed() per target
R/mi_pca_imp.R      .mipca_run()   -- every door above lands here
  +- R/mi_pca_check.R       .chk_matrix(), .chk_cells()
  +- R/mi_pca_transform.R   .resolve_transform()  -- a name or an object -> a transform
  +- R/mi_pca_imp.R         .resolve_cols()  -- subset and abb, by name or number
  +- R/mi_pca_abb.R         .abb_split()  -- abb columns out of obj, subset rebased
  +- R/mi_pca_imp.R         .nothing_to_impute() -> .mipca_empty()  (no sampling)
  +- R/col_scan.R           col_scan(), refuse_nonfinite()
  +- R/mi_pca_imp.R         .check_admissible() -> .admissible_scan()
  +- R/mi_pca_transform.R   tf_biwhiten()  -- its fit step, run once per call on obj
  +- R/mi_pca_imp.R         mipca_core()  -- the sampler's R half
  |    +- R/mi_pca_state.R        .state_check_args(), .state_check_match()  (resume only)
  |    +- R/mi_pca_imp.R          .core_build() -> .mipca_build_stream()
  |    |    +- R/RcppExports.R          .build_block()  -- .Call shim, per 1024 columns
  |    |    |    +- src/mipca.cpp           build_block()  -- centre, record holes, dsyrk clean cols
  |    |    +- R/RcppExports.R          .sym_from_lower()
  |    |         +- src/mipca.cpp           sym_from_lower()
  |    +- R/mi_pca_imp.R          .start_fresh()  (or R/mi_pca_state.R .start_resume())
  |    +- R/mi_pca_imp.R          .alloc_draws()
  |    +- R/RcppExports.R         .mipca_da_loop_multi()  -- .Call shim
  |         +- src/RcppExports.cpp     _slideimp_mipca_da_loop_multi (SEXP wrapper)
  |              +- src/mipca.cpp       mipca_da_loop_multi()   <-- C++ entry
  |                   +- src/mipca.cpp       da_ctx_build(), spectra_ncv_choose()
  |                   +- src/mipca.cpp       need_double()
  |                   +- src/mipca.cpp       acc.clone()  -- one dqrng generator per chain
  |                   +- src/mipca.cpp       da_chain_run()  -- each chain in turn
  |                        +- src/mipca.cpp       mt_chisq()  (Student-t noise only)
  |                        +- src/mipca.cpp       sweep_cols()
  |                        +- BLAS dsyrk  -- the Gram
  |                        +- src/mipca.cpp       spectra_solve()
  |                        |    +- src/Spectra/SymEigsSolver.h   SymEigsSolver
  |                        +- BLAS dgemm  -- the other side's factor
  |                        +- src/mipca.cpp       moy_dot()
  +- R/mi_pca_s3.R          .si_mi_pca_new()  -- the resumable fit, an environment
  +- R/mi_pca_s3.R          mipca_finalize()  (keep_state = FALSE)
       +- R/mi_pca_s3.R          .diag_frame() -> .diag_sample(), .cell_diag_read()
       |    +- R/RcppExports.R         mipca_cell_diag()
       |         +- src/mipca_diag.cpp     mipca_cell_diag()  -- Rhat and ESS per cell
       +- R/mi_pca_s3.R          .rebase_subset(), .si_mi_chunk(), .si_mi_new()
       +- R/mi_pca_abb.R         .abb_join() -> .abb_piece()  (abb columns only)
R/mi_pca_s3.R       update.si_mi_pca()  -- mipca_core() with .state, in place
R/mi_pca_s3.R       as.matrix.si_mi()   -- one completed matrix, backward map per chunk
R/mi_pca_s3.R       summary.si_mi(), summary.si_mi_pca()
R/mi_pca_s3.R       cbind.si_mi() -> .si_mi_join()
R/mi_pca_tune.R     .score_holdout() -> crps_all(), pit_all()
  +- R/RcppExports.R        crps_block()
       +- src/mipca_crps.cpp     crps_block()  -- CRPS of each held-out cell
```

## Layer 1: the doors

Four exported functions start a fit, and all four reach `.mipca_run()`
(`R/mi_pca_imp.R`: `.mipca_run <- function(`). They differ in how they seed
and what they ask it to keep.

- `mipca()` (`R/mi_pca_imp.R`: `mipca <- function(`) passes its arguments
  through with `seed_given = TRUE`, so `mipca_core()` seeds dqrng itself
  (Random numbers, below).
- `mipca_fit()` (`R/mi_pca_tune.R`: `mipca_fit <- function(`) reads the
  transform, `abb` and sampler settings off a `si_tune` table, checks the
  matrix against the table's fingerprint (`R/mi_pca_tune.R`:
  `.tune_fp_check <- function(fp, X) {`), takes `nu` from the held-out
  cells of the chosen rank, floored at `.NU_FLOOR` (`R/mi_pca_tune.R`:
  `nu <- max(nu_implied, .NU_FLOOR)`), and calls `.mipca_run()` through
  `do.call()` with `seed_given = TRUE` (`R/mi_pca_tune.R`: `.mipca_run,`).
  The help page promises that `mipca_fit()` returns the same fit as `mipca()`
  for the same settings.
- `mipca_tune()` (`R/mi_pca_tune.R`: `mipca_tune <- function(`) blanks the
  holdout in a copy of `obj`, then runs one arm per rank, each a
  `.mipca_run()` with `keep_state = TRUE` and an explicit `store`
  (`R/mi_pca_tune.R`: `store = st$store,`). It seeds once, before the arms
  (`R/mi_pca_tune.R`: `dq_old <- .seed_dq(seed)`), and the arms do not seed.
  A failed arm is caught and becomes a row with its message
  (`R/mi_pca_tune.R`: `error = identity`).
- `mipca_flank()` (`R/mi_pca_flank.R`: `mipca_flank <- function(`) takes a
  `si_flank_plan` and runs one `.mipca_run()` per target with status
  `"impute"`, on the target's window of columns, keeping only the target
  (`R/mi_pca_flank.R`: `subset = plan$local[i],`). It seeds dqrng per target
  from `seed` and the target's column (`R/mi_pca_flank.R`:
  `dqrng::dqset.seed(c(seed, tg[i]))`), so the core does not seed. Targets
  with another status become `.abb_piece()` or `.mipca_empty()` fits, and
  every piece is joined in target order (`R/mi_pca_flank.R`:
  `fit <- .si_mi_join(fits, call = cl)`). A window's error is re-raised with
  the target's name; the `nu` floor warning is given once at the door and
  muffled per window (`R/mi_pca_flank.R`:
  `if (identical(conditionMessage(w), .nu_floor_msg())) {`).

`update()` on a resumable fit does not go through `.mipca_run()`: it calls
`mipca_core()` directly with the saved state (Resume, below).

## Layer 2: `.mipca_run()` (`R/mi_pca_imp.R`)

In order:

1. `.chk_matrix()` (`R/mi_pca_check.R`: `.chk_matrix <- function(obj) {`)
   requires a double matrix with unique row and column names. It refuses
   anything else through checkmate and never converts it. `.chk_cells()` refuses a matrix with more
   cells than an integer can index (`R/mi_pca_check.R`: `.chk_cells <- function(`).
2. `.resolve_transform()` (`R/mi_pca_transform.R`:
   `.resolve_transform <- function(transform, n) {`) maps `"identity"` and
   `"logit_biwhiten"` to transform objects, the second being
   (`R/mi_pca_transform.R`: `logit_biwhiten = tf_biwhiten(tf_sv(n), rows = TRUE),`),
   and accepts a list with `forward` and `backward`, or with `fit`.
3. `subset` and `abb` through `.resolve_cols()` (`R/mi_pca_imp.R`:
   `.resolve_cols <- function(subset, obj, what = "subset") {`), which refuses
   an unknown or repeated name itself rather than letting `resolve_subset()`
   drop it (`dev/DECISIONS.md`).
4. With `abb`, `.abb_split()` (`R/mi_pca_abb.R`:
   `.abb_split <- function(obj, abb, subset) {`) records the requested `abb`
   columns and how to put them back in order, and `obj` is replaced by a copy
   of the kept columns (`R/mi_pca_imp.R`: `obj <- obj[, sp$kept, drop = FALSE]`).
   When every requested column is in `abb`, no model is fitted.
5. The check for nothing to impute (`R/mi_pca_imp.R`:
   `.nothing_to_impute <- function(obj, subset, scan = NULL) {`). With
   `subset` it reads only those columns; otherwise it runs the scan of every
   column, with the variances, and the fit then reuses that scan. When nothing is to
   impute, `.mipca_empty()` (`R/mi_pca_imp.R`: `.mipca_empty <- function(`)
   builds a fit with no draws and finalizes it, so the result has the same
   class as a sampled one. With `keep_state = TRUE` the call is refused
   instead, since there is nothing to resume.
6. The full scan, `refuse_nonfinite()` (`R/col_scan.R`:
   `refuse_nonfinite <- function(scan, empty_j = NULL) {`), and, for a
   transform with a `fit` step, `.check_admissible()` and the fit itself
   (`R/mi_pca_imp.R`: `transform <- transform$fit(obj, rank)`). The fit step
   of `tf_biwhiten()` (`R/mi_pca_transform.R`: `fit = function(X, rank) {`)
   runs once per call. It estimates column scales, and row scales with
   `rows = TRUE`, from a rank-`S` residual of the observed cells, and returns
   a transform whose `forward` and `backward` close over those scales. Its scales depend on
   `rank`, so each arm of `mipca_tune()` fits its own.
7. `mipca_core()`, then `.si_mi_pca_new()`, then `mipca_finalize()` unless
   `keep_state` (`R/mi_pca_imp.R`:
   `if (keep_state) fit else mipca_finalize(fit, prop_diag = prop_diag)`).

`scan` lets a caller that already scanned a larger matrix hand in its scan
indexed to `obj`, so `obj` is not read again. Only `mipca_flank()` passes one
(`R/mi_pca_flank.R`: `scan = scan_cols(scan, at[wcols])`), and `scan` cannot
be combined with `abb`.

**Admissibility.** `.admissible_scan()` (`R/mi_pca_imp.R`:
`.admissible_scan <- function(Z, scan, cols = NULL) {`) reads, per column,
whether it is empty, has one observed cell, has a variance that overflowed,
or is constant (`R/mi_pca_imp.R`: `const <- fin & v < .Machine$double.eps`).
`.check_admissible()` (`R/mi_pca_imp.R`: `.check_admissible <- function(`)
refuses the first three, with the columns named, and points at `abb` for a
column with one observed cell. It does not refuse a constant column: the
core fits one, and its centred column is zero. The flank plan leaves a
constant column out of every window instead (`R/mi_pca_flank.R`:
`!(scan$empty | scan$single | scan$const | too_miss | is_abb[read])`).

## Layer 3: `mipca_core()` (`R/mi_pca_imp.R`)

`mipca_core()` (`R/mi_pca_imp.R`: `mipca_core <- function(`) is the R half of
the sampler, shared by a fresh fit and a resume. It takes `Z`, already on the
input's scale and already admitted, and `forward`, the map to the model's
scale.

Settings are checked first, before R's random number stream is drawn
(`tests/testthat/test-mipca-refuse-before-draw.R`:
`bad settings are refused before R's stream is drawn`). `rank` must leave
two dimensions spare (`R/mi_pca_check.R`: `if (S + 2L > min(n, p)) {`): the
C++ solves `S + 1` eigenpairs on a Gram of dimension `min(n, p)` (Layer 5),
and Spectra solves fewer pairs than the dimension (`src/mipca.cpp`:
`const int k_solve = S + 1;`). `nu` must exceed 2, and below `.NU_FLOOR` a fresh fit warns
(`R/mi_pca_imp.R`: `if (is.null(.state) && nu < .NU_FLOOR) {`).

Two constants of the model are fixed here:

- `K <- min(n - 1L, p)` (`R/mi_pca_imp.R`: `K <- min(n - 1L, p)`), the
  number of components a centred `n x p` matrix can hold;
- `dfP`, the residual degrees of freedom of a rank-`S` fit with column means
  (`R/mi_pca_imp.R`:
  `dfP <- as.double(n) * p - (p + as.double(S) * (n - 1L + p - S))`).

**The hole enumeration.** Every hole of `Z` has a position `1..nmiss` in
column-major order. `mi_row` holds its row, `mi_col` its column, `col_ptr`
the first position of each column, and `miss_idx_l` its linear index into
`Xd`, the matrix of dirty columns (Layer 4) (`R/mi_pca_imp.R`:
`miss_idx_l <- mi_row + (mi_col_l - 1L) * n`). The cells whose draws are kept,
`store_idx`, are positions in this enumeration: every hole of the `subset`
columns by default (`R/mi_pca_imp.R`:
`store_idx <- sequence(cnt[cols_eff], from = col_ptr[cols_eff])`), or the
caller's `store`, strictly increasing (`R/mi_pca_imp.R`:
`store_idx <- .chk_store(store_use, nmiss)`). `subset` and `store` exclude
each other. `map`, the fit's `(row, col)` matrix of kept cells, is built from
`store_idx`, so row `j` of every draws matrix is cell `map[j, ]`.

The rest runs in this order: the build (Layer 4); seeding (Random numbers);
the start of each chain, from `.start_fresh()` or `.start_resume()`; the
draws; and the `.Call` (`R/mi_pca_imp.R`: `res <- .mipca_da_loop_multi(`).
The draws are one `n_store x iter_sampling` matrix of `NA` per chain,
allocated through `.alloc_draws()` (`R/mi_pca_imp.R`:
`.alloc_draws <- function(expr, cells, chains, draws, fix) {`), which turns a
failed allocation into a refusal that names the size in GB. With `keep_state` it
then assembles the state (Resume, below).

**A fresh start.** `.start_fresh()` (`R/mi_pca_imp.R`:
`.start_fresh <- function(`) gives each chain:

- its own copy of `Xd`: chain 1 takes the built matrix, every other chain a
  fresh one (`R/mi_pca_imp.R`: `function(ch) if (ch == 1L) Xd else Xd + 0`),
  since the C++ writes into it;
- a starting noise level from the observed sum of squares
  (`R/mi_pca_imp.R`: `sigma_in = rep(sqrt(sum(col_ss) / df_obs), chains),`);
- a dispersion per dirty column, `init` times the column's observed standard
  deviation (`R/mi_pca_imp.R`:
  `disp = (init_c * sqrt(col_ss / (nobs - 1L)))[dirty]`), from which the C++
  draws each chain's own starting mean for the holes;
- empty `xtild_miss` and `eig_start`, which the C++ requires on a fresh run.

## Layer 4: the build (`R/mi_pca_imp.R`, `src/mipca.cpp`)

`.core_build()` (`R/mi_pca_imp.R`:
`.core_build <- function(Z, cnt, .state, forward) {`) calls
`.mipca_build_stream()` (`R/mi_pca_imp.R`:
`.mipca_build_stream <- function(X, cnt, forward, resume) {`), which walks
`Z` in blocks of 1024 columns (`R/mi_pca_imp.R`: `block <- 1024L`) and hands
each to `build_block()` (`src/mipca.cpp`:
`void build_block(SEXP src, int src_off, int j1, int nb, int transformed,`).

**The transform, a block at a time.** With the identity transform the whole
matrix goes to C++ with a column offset, and nothing is copied
(`R/mi_pca_imp.R`: `in_place <- identical(forward, identity)`). Otherwise
`forward()` is called on each block with the block's column indices, and
the block it returns is refused unless it is double. `build_block()` takes
`src` as a `SEXP` and refuses anything but `REALSXP` (`src/mipca.cpp`:
`if (TYPEOF(src) != REALSXP)`), so no conversion can happen at the boundary.
It counts the holes of each column again and compares them with the scan
(`src/mipca.cpp`: `if (holes != cnt[j]) {`): a transform that maps an
observed cell to `NA` or a missing cell to a value is refused, with the
column named, as is one that gives a non-finite value
(`tests/testthat/test-mipca-messages.R`:
`a transform that breaks the hole pattern is named with its fix`). This is
the only place the transformed values are checked.

**Dirty and clean columns.** `.split_cols()` decides the layout
(`R/mi_pca_imp.R`: `split <- n <= p && any(cnt == 0L)`):

- split (`n <= p` with at least one complete column): only the columns with
  holes, the dirty ones, go into `Xd`. Each complete, clean column is centred
  into a scratch block and folded into `g_fixed`, the `n x n` Gram of the
  clean columns, by one `dsyrk` per block (`src/mipca.cpp`:
  `(&uplo, &trans, &n, &kc, &one, sc, &n, &one, REAL(g_fixed), &n FCONE FCONE);`).
  Those columns are never touched again: the loop adds `Xd`'s contribution to
  `g_fixed` at every iteration. `sym_from_lower()` (`src/mipca.cpp`:
  `void sym_from_lower(Rcpp::NumericMatrix g) {`) mirrors its lower triangle
  once the build ends.
- otherwise every column is in `Xd`, clean or not, and there is no `g_fixed`.

Every column of `Xd` holds its observed cells minus the column's observed
mean, and 0 in its holes. `build_block()` also writes, per column, the mean
into `et`, the centred sum of squares into `col_ss`, and the rows of its
holes into `mi_row`.

On a resume the build skips the clean columns entirely and takes `g_fixed`
from the state; it still recomputes `Xd` and the means from `Z`, and refuses
if the holes have moved (`R/mi_pca_imp.R`:
`if (!is.null(.state) && !identical(out$mi_row, .state$mi_row)) {`).

## Layer 5: the Rcpp boundary and chain setup

`R/RcppExports.R` defines `.mipca_da_loop_multi()` (`R/RcppExports.R`:
`.mipca_da_loop_multi <- function(`), a single `.Call` to
`_slideimp_mipca_da_loop_multi` in `src/RcppExports.cpp`, which calls
`mipca_da_loop_multi()` (`src/mipca.cpp`: `Rcpp::List mipca_da_loop_multi(`),
exported under the `// [[Rcpp::export(.mipca_da_loop_multi)]]` marker. The
other two exports of `src/mipca.cpp` are `.build_block` and
`.sym_from_lower`, both `void`: they write into the R objects they are
given. `src/mipca.cpp` declares its dependencies on RcppEigen, dqrng and BH
(`src/mipca.cpp`: `// [[Rcpp::depends(RcppEigen, dqrng, BH)]]`).

`mipca_da_loop_multi()` takes one list entry per chain for each piece of
chain state (`Xd`, draws, column means, hole means, eigensolver start,
`sigma`) and one `streams` entry per chain. Each list arrives as a `SEXP`,
is refused unless it is a list (`src/mipca.cpp`:
`static void need_list(SEXP x, const char *what) {`) and is read with
`VECTOR_ELT()`, so no `Rcpp::List` is built from an argument;
`mipca_cell_diag()` takes its `draws` the same way. It checks every shape, refuses a
dispersion on a resume and a saved state on a fresh run, and checks that each
chain's buffers are double (`src/mipca.cpp`:
`static void need_double(SEXP x, const char *what, int chain) {`), so the
`Rcpp::NumericMatrix` it wraps around each is the R object itself and not a
copy. `da_ctx_build()` (`src/mipca.cpp`:
`static da_ctx da_ctx_build(int n, int nd, Rcpp::NumericVector g_fixed,`)
checks the shared inputs and fills a `da_ctx` (`src/mipca.cpp`:
`struct da_ctx {`), which holds raw pointers into the R objects; each chain's
copy of it points at that chain's buffers.

Two choices are made here, once per call:

- **The Gram's side.** `tall` is true when there is no `g_fixed` and
  `n > nd` (`src/mipca.cpp`: `const bool tall = !has_fixed && n > nd;`).
  Then the Gram is `Xd' Xd`, `nd x nd`; otherwise it is `Xd Xd'` plus
  `g_fixed`, `n x n`. `ng` is its dimension, and the eigensolver works on the
  smaller side of the matrix it can see.
- **`ncv`**, the number of Lanczos vectors (`src/mipca.cpp`:
  `static int spectra_ncv_choose(int n, int k, int ncv_req) {`): the
  caller's `spectra_opts(ncv = )` if given, else `2k + 1` for the `k = S + 1`
  pairs solved, at least `SPECTRA_NCV_FLOOR` and at most `ng`.

Each chain gets its own generator, a copy of dqrng's (`src/mipca.cpp`:
`gens.push_back(acc.clone(static_cast<rng_result_t>(streams[c])));`), or the
generator read back from its saved text on a resume (`src/mipca.cpp`:
`is >> *g;`). Then the chains run one after another, on this thread
(`src/mipca.cpp`:
`for (int c = 0; c < nc; ++c) sigma_out[c] = da_chain_run(ctxs[c], *gens[c]);`),
and each generator is written back out as text (`src/mipca.cpp`:
`os << *gens[c];`). The result is a list of `sigma`, `etd`, `xtild_miss` and
`eig_start` per chain and `rng_state`; the draws and each chain's `Xd` have
been written in place.

## Layer 6: one iteration (`src/mipca.cpp`)

`da_chain_run()` (`src/mipca.cpp`: `const int nopp = tall ? n : nd;`) runs
`iter_warmup + iter_sampling * thin` iterations. Its state between
iterations is `Xd` with its current holes, `etd` (the column means on the
model's scale), `xtild` (the current mean of each hole), `sigma`, and
`eig_start`. On a fresh run `xtild` starts as each hole's column dispersion
times a normal draw (`src/mipca.cpp`:
`xtild[k] = cx.disp[mcol[k] - 1] * z(rng_eng);`), that is around its column
mean, and the first eigensolve starts cold. The phases of one iteration, in
the order the code runs them, named as its timer labels them:

1. **istep.** Each hole of `Xd` becomes its mean plus noise. The noise is
   Gaussian
   (`src/mipca.cpp`: `xp[midx[k] - 1] = xtild[k] + sigma * z(rng_eng);`), or
   Student-t with `nu` degrees of freedom scaled to the same variance
   (`src/mipca.cpp`:
   `const double t_scale = t_noise ? std::sqrt((cx.nu - 2.0) / cx.nu) : 0.0;`).
   A Student-t draw is a normal over the root of a chi-square that
   `mt_chisq()` draws by Marsaglia and Tsang's method (`src/mipca.cpp`:
   `static inline double mt_chisq(dqrng::random_64bit_generator &g,`).
2. **store.** After warmup, on every `thin`-th iteration (`src/mipca.cpp`:
   `if (l > warmup && (l - warmup) % thin == 0) {`), each kept hole's value
   plus its column mean goes into the next column of the chain's draws
   (`src/mipca.cpp`: `dcol[j] = xp[midx[sidx[j] - 1] - 1] + etd[scol[j] - 1];`).
   A draw is therefore the value the I-step just drew, on the model's scale,
   not centred. `as.matrix()` maps it back (Layer 7).
3. **sweep.** Each column of `Xd` is centred again, with its new holes, and
   the shift is added to `etd` (`src/mipca.cpp`:
   `static void sweep_cols(double *x, int n, int nd, double *etd) {`).
4. **gram.** With `g_fixed`, its lower triangle is copied in
   (`src/mipca.cpp`:
   `std::copy(gfp + off + c, gfp + off + ng, gram.data() + off + c);`), then
   one `dsyrk` adds `Xd`'s part, on the side `tall` chose (`src/mipca.cpp`:
   `(&uplo, gram_trans, &ng, &nopp, &one, xp, &n, gram_beta, gram.data(),`).
   Only the lower triangle is ever written; Spectra's
   `DenseSymMatProd<double>` reads the lower triangle by default
   (`src/Spectra/MatOp/DenseSymMatProd.h`:
   `template <typename Scalar_, int Uplo = Eigen::Lower, int Flags = Eigen::ColMajor>`).
   The Gram's trace is the total sum of squares
   (`src/mipca.cpp`:
   `for (int i = 0; i < ng; ++i) ss += gram[static_cast<R_xlen_t>(i) * ng + i];`).
5. **solve.** Spectra's `SymEigsSolver` finds the `S + 1` largest
   eigenpairs (`src/mipca.cpp`:
   `Spectra::SymEigsSolver<Spectra::DenseSymMatProd<double>> eigs(op, Sk, ncv);`)
   through `spectra_solve()` (`src/mipca.cpp`:
   `static int spectra_solve(EigsT &eigs, int n, int k, double tol, int maxitr,`).
   After the first iteration of a fresh run, and from the first of a resume,
   it starts from the sum of the previous solve's eigenvectors
   (`src/mipca.cpp`: `for (int i = 0; i < n; ++i) start[i] += col[i];`).
   Two refusals stop the chain: a solve that does not converge
   (`src/mipca.cpp`:
   `Rcpp::stop("the PCA step did not converge at iteration %d. Raise "`), and
   an `S`-th eigenvalue that is not positive or too small next to the first
   (`src/mipca.cpp`:
   `if (evals[S - 1] <= 0.0 || evals[0] / evals[S - 1] > gram_max_cond)`),
   with the bound set in R (`R/mi_pca_imp.R`:
   `.GRAM_MAX_COND <- 1e-8 / .Machine$double.eps`).
6. **backout.** The noise variance is the sum of squares the `S` components
   leave, over `dfP` (`src/mipca.cpp`:
   `const double sigma2 = std::max(ss - sd2, 0.0) / dfP;`). Each component is
   shrunk by `phi_s = max(0, 1 - floor / lambda_s)` (`src/mipca.cpp`:
   `double ph = 1.0 - floor / evals[s];`), where the floor is
   `n p / K` times `sigma2`, capped at the `(S + 1)`-th eigenvalue
   (`src/mipca.cpp`: `const double floor = std::min(floor_raw, evals[S]);`).
   That cap is why `S + 1` pairs are solved. One `dgemm` projects `Xd` onto
   the eigenvectors, giving the factor of the other side (`src/mipca.cpp`:
   `(proj_trans, &trans_n, &nopp, &S, &ng, &one, xp, &n, wmat.data(), &ng,`).
7. **moy.** `moy_dot()` (`src/mipca.cpp`:
   `static void moy_dot(const double *u, const double *v, const double *w,`)
   forms the shrunk rank-`S` reconstruction at the holes only: for hole `k`,
   the sum over `s` of the row factor, the column factor and `phi_s`. It
   transposes both factors first, so the inner loop reads `S` contiguous
   values per hole.
8. **pstep_noise.** Each hole's mean takes one more normal draw, of variance
   `sigma2` times the sum of the `phi_s` over `K` (`src/mipca.cpp`:
   `const double sd_t = std::sqrt(sigma2 * sphi / static_cast<double>(K));`).
   This is the P-step's draw of the mean; the next I-step adds the noise.

Every iteration ends with a progress line when `refresh` asks for one, and
with `Rcpp::checkUserInterrupt()` (`src/mipca.cpp`:
`Rcpp::checkUserInterrupt();`). The progress line is sent through
`Rcpp::message()` (`src/mipca.cpp`:
`tfm::format("Chain %d Iteration: %*d / %d [%3d%%] (%s)", chain_id,`), so
`suppressMessages()` silences it (`tests/testthat/test-mipca-messages.R`:
`progress lines are messages, and suppressMessages() silences them`). After the last iteration the chain copies
out `etd`, `xtild` and `eig_start` and returns `sigma`.

The rank bound is pinned at its edge (`tests/testthat/test-mipca-gram-bound.R`:
`rank min(n, p) - 2 fits and min(n, p) - 1 refuses`).

## Layer 7: from draws to a result

**The resumable fit.** `.si_mi_pca_new()` (`R/mi_pca_s3.R`:
`.si_mi_pca_new <- function(`) stores `obj` (the matrix the model saw, after
any `abb` split), `map`, the per-chain draws, the transform, the state and a
provenance list in an environment (`R/mi_pca_s3.R`:
`e <- new.env(parent = emptyenv())`) of class `si_mi_pca`. Being an
environment, it has reference semantics: `update()` changes it in place, and
`b <- fit` does not copy it. Its `lifecycle` is `"ok"`, `"failed"` or
`"passed"`, and `.life_check()` (`R/mi_pca_s3.R`: `.life_check <- function(x) {`)
refuses every method on a fit that is not `"ok"`, `print()` excepted
(`tests/testthat/test-mipca-lifecycle.R`:
`a fit torn by a failed update() refuses every method but print()`).

**Finalizing.** `mipca_finalize()` (`R/mi_pca_s3.R`:
`mipca_finalize <- function(x, thin = 1L, prop_diag = 0.5) {`) keeps every
`thin`-th draw of each chain (`R/mi_pca_s3.R`:
`keep <- seq.int(1L, nd, by = k)`), computes Rhat and ESS on a sample of the
kept cells, narrows the fit to its `subset` columns, binds the chains' draws
side by side into one matrix, and returns a `si_mi`. It then empties the
resumable fit and marks it `"passed"` (`R/mi_pca_s3.R`:
`x$lifecycle <- "passed"`), so the draws are not held twice.

- The diagnostic sample (`R/mi_pca_s3.R`:
  `.diag_sample <- function(map, p, prop_diag) {`) takes
  `ceiling(prop_diag * h)` of each column's `h` kept cells, evenly spaced
  (`R/mi_pca_s3.R`:
  `k <- min(h, max(1L, as.integer(ceiling(prop_diag * h))))`), so it is the
  same on every call.
- `mipca_cell_diag()` (`src/mipca_diag.cpp`:
  `Rcpp::NumericMatrix mipca_cell_diag(SEXP draws, Rcpp::IntegerVector rows,`)
  reads those rows of every chain in tiles of 256 cells (`R/mi_pca_s3.R`:
  `mipca_cell_diag(draws, srows, 256L)`) and returns, per cell, the larger of
  the split Rhat of the draws and of their distances from the median
  (`src/mipca_diag.cpp`:
  `rhat = (ISNAN(a) || ISNAN(b)) ? NA_REAL : std::max(a, b);`), and the ESS
  over all chains.
- `.rebase_subset()` (`R/mi_pca_s3.R`: `.rebase_subset <- function(x, fr) {`)
  cuts `obj` to the `subset` columns and renumbers `map` and the diagnostic
  cells to match; the transform is wrapped by `.tf_narrow()` (`R/mi_pca_s3.R`:
  `.tf_narrow <- function(tf, subset) {`) so its `backward()` still receives
  the column numbers of the matrix it was fitted on.

**The finalized fit.** A `si_mi` (`R/mi_pca_s3.R`:
`.si_mi_new <- function(chunks, rownames, nkeep, call, placement = NULL) {`)
is a list of chunks, each a matrix of columns with its own `map`, draws,
transform and diagnostics, plus the column offsets of the chunks and, when
the chunks' columns are not in input order, a `placement` permutation. One
`mipca()` fit is one chunk, or two with `abb`; `cbind()` and `mipca_flank()`
join chunks with `.si_mi_join()` (`R/mi_pca_s3.R`:
`.si_mi_join <- function(parts, call, placement = NULL) {`), which requires
the same row names and the same number of kept draws. `as.matrix()`
(`R/mi_pca_s3.R`: `as.matrix.si_mi <- function(x, draw = 1L, ...) {`) fills
each chunk's holes with one column of its draws, through the chunk's
`backward()` (`R/mi_pca_s3.R`:
`out[idx] <- .backward_at(ch$transform, ch$draws[, draw], ch$map)`). Draw
`b` of a joined fit is draw `b` of every chunk.

**`abb` columns.** They never enter the model. At finalize, `.abb_join()`
(`R/mi_pca_abb.R`: `.abb_join <- function(`) has `.abb_piece()`
(`R/mi_pca_abb.R`: `.abb_piece <- function(`) draw each listed column's holes
by approximate Bayesian bootstrap from its own observed cells, seeded by the
fit's seed and the column (`R/mi_pca_abb.R`:
`dqrng::dqset.seed(c(seed, key[i]))`), and joins the piece to the model's
chunk with the placement that restores the input's column order.

**Tuning reads the resumable fit.** Each `mipca_tune()` arm keeps its fit
resumable so it can read the draws and the end-of-chain `sigma` without
finalizing. `.tune_arm_read()` (`R/mi_pca_tune.R`:
`.tune_arm_read <- function(fit, rec, rank, took) {`) pools the chains'
holdout rows (`R/mi_pca_tune.R`: `D <- .pool_rows(fit$draws, rec$rows_ho)`),
maps them back, and scores them by CRPS and PIT (`R/mi_pca_tune.R`:
`crps = crps_all(d, truth),`). `crps_block()` (`src/mipca_crps.cpp`:
`Rcpp::NumericVector crps_block(Rcpp::NumericMatrix d,`) computes each cell's
CRPS from its sorted draws, with the sum over draws computed by one `dgemv`. The
store set of an arm is the holdout cells plus the `n_diag` cells sampled for
Rhat and ESS, merged into one increasing set of hole positions
(`R/mi_pca_tune.R`: `st <- .tune_store(hp$ho_pos, diag_pos)`).

## Resume

A fit made with `keep_state = TRUE` carries a state of class
`si_mi_pca_state` (`R/mi_pca_imp.R`: `class = "si_mi_pca_state"`), versioned
by `.STATE_VERSION` (`R/mi_pca_state.R`: `.STATE_VERSION <- 10L`). It holds
the settings the chains ran under (rank, thin, spectra settings, `nu`), the
row count and holes per column, `mi_row`, `g_fixed`, the transform's
`forward`, the store set, and per chain: the values of its holes in `Xd`
(`R/mi_pca_imp.R`: `holes = xd_list[[ch]][miss_idx_l],`), `etd`, `sigma`,
`xtild_miss`, `eig_start`, the generator as text, and a fingerprint of `Xd`.

`update.si_mi_pca()` (`R/mi_pca_s3.R`: `update.si_mi_pca <- function(`) calls
`mipca_core()` with the state and no warmup (`R/mi_pca_s3.R`:
`iter_warmup = 0L,`). The core then:

1. refuses a state of another version (`R/mi_pca_state.R`:
   `.state_check_version <- function(.state) {`), and any `forward`, `store`
   or other `chains` passed with one (`R/mi_pca_state.R`:
   `.state_check_args <- function(.state, chains, chains_given, forward, store) {`);
2. refuses a matrix whose row count or holes per column changed
   (`R/mi_pca_state.R`:
   `.state_check_match <- function(.state, n, cnt, rank, nu, spectra, Z) {`),
   and, for a split matrix, whose clean columns changed, by the relative
   change of their absolute sums (`R/mi_pca_state.R`:
   `.state_fp_clean_check <- function(.state, Z, n, cnt) {`);
3. rebuilds `Xd` from `Z` with the saved `forward` (Layer 4) and, per chain,
   shifts it from the rebuilt means to the chain's own and puts the chain's
   holes back (`R/mi_pca_state.R`:
   `Xc <- Xd + rep(etd - st_ch$etd, each = n)`), then compares the
   fingerprint, the absolute sums of at most `max(100, 5%)` of the dirty
   columns (`R/mi_pca_state.R`:
   `m <- min(nd, max(100L, as.integer(ceiling(0.05 * nd))))`), within
   `.STATE_FP_TOL` (`R/mi_pca_state.R`: `.STATE_FP_TOL <- 1e-8`);
4. passes `resuming = TRUE`, the saved `xtild_miss`, `eig_start`, `sigma` and
   generators, and no dispersion.

`update()` then binds the new draws to the old, chain by chain. It marks the
fit `"failed"` only after the core has returned (`R/mi_pca_s3.R`:
`object$lifecycle <- "failed"`) and `"ok"` once every chain is bound and the
state replaced: an error inside the core leaves the fit as it was, and one
while binding leaves the chains with different numbers of draws, which the
mark records. A resumed fit of `a + b` draws agrees with one run of `a + b`
under the same seed to rounding, not bitwise: the first `a` draws are
identical, and the rebuild of `Xd` from `Z` and the shift by `etd` change
the last bits of what follows. A test pins it on the split and the tall
layout, Gaussian and Student-t (`tests/testthat/test-mipca-resume.R`:
`update() extends a fit as one longer run would (`).

## Random numbers

All sampling draws come from dqrng; base R's stream is touched only by
`.resolve_seed()` when `seed = NULL` (`R/mi_pca_imp.R`:
`seed <- sample.int(.Machine$integer.max, 1L)`), and by `mipca_tune()`'s
holdout and scoring.

- With `seed_given = TRUE` (`mipca()`, `mipca_fit()`), `mipca_core()` seeds
  dqrng (`R/mi_pca_imp.R`: `dq <- .seed_dq(seed)`) and restores the state it
  found on exit (`R/mi_pca_imp.R`:
  `on.exit(dqrng::dqrng_set_state(dq), add = TRUE)`).
- Otherwise the core does not seed, and the chains copy whatever state dqrng
  holds: `mipca_tune()` seeds once for every arm, `mipca_flank()` once per
  window.
- Chain `c` is a copy of dqrng's generator on stream `c - 1`
  (`R/mi_pca_imp.R`: `streams = seq_len(chains) - 1L,`). dqrng's `clone()`
  copies the generator and sets the stream on the copy, so the sampler never
  advances dqrng's own generator: two calls from one dqrng state draw the
  same numbers. That held when it was measured: the chains differed from
  each other, and the call left the dqrng state unchanged. Nothing between `mipca_tune()`'s
  arms draws from dqrng, so every rank runs from the same generator state.
- On a resume each chain's generator is read back from its saved text, so
  the chain continues the stream it was on.

## Data structures and ownership

- `obj` is double from the door and nothing converts it. The identity build
  reads it in place; a transform's `forward()` gets a copy of each block.
  `.abb_split()` and `mipca_flank()`'s windows copy the columns they keep.
- `build_block()` and `sym_from_lower()` return nothing: they write into
  `Xd`, `g_fixed`, `et`, `col_ss` and `mi_row`, allocated in
  `.mipca_build_stream()` with the types Rcpp wraps without copying. An
  argument of another type would be wrapped around a converted copy, and
  the writes would be lost.
- `mipca_da_loop_multi()` writes into each chain's `Xd` and draws matrix in
  place. Each chain has its own `Xd` (Layer 3), so no two chains share a
  buffer; after the call `xd_list` holds each chain's last state, which
  `keep_state` reads.
- Memory, per call: one `n x nd` matrix and one `n_store x iter_sampling`
  draws matrix per chain; with a split matrix one `n x n` `g_fixed`; per
  chain, while it runs, an `ng x ng` Gram, the `ng x (S + 1)` eigenvectors,
  and vectors of length `nmiss`. `mipca_finalize()` binds the chains' kept
  draws into one matrix before it empties the resumable fit.
- The state holds two vectors of `nmiss` values per chain, `holes` and
  `xtild_miss`, over every hole and not only the kept ones.

## Compile-time switches

- **`LOC_TIMER`**, set by `load_all1()` by default (`dev/dev-utils.R`): the
  phases of Layer 6 are timed under their labels, and Spectra's matrix-vector
  products under `s_matvec` through a wrapper (`src/mipca.cpp`:
  `class TimedSymMatProd {`), in one rcpptimer object
  (`src/mipca.cpp`: `LOC_TIMER_OBJ(v2_da_phases);`). Such a build refuses
  more than one chain (`src/mipca.cpp`: `if (nc > 1)`), so the `test-mipca-*`
  files need `load_all1(timer = FALSE)`.
- `PCA_IMP_DIAGNOSTICS` does not reach `src/mipca*.cpp`.
- No `EIGEN_USE_BLAS` and no OpenMP flag is set in `src/Makevars` or
  `src/Makevars.win`, so Spectra's matrix-vector product is Eigen's own
  single-threaded code. From C++, only these calls reach R's BLAS: the
  `dsyrk` of the build, the `dsyrk` and `dgemm` of each iteration, and the
  `dgemv` of the CRPS.
  Every BLAS call goes through
  `R_ext/BLAS.h` with `USE_FC_LEN_T` and one `FCONE` per character argument.

## Parallelism

The sampler starts no thread of its own and pins no BLAS (`dev/DECISIONS.md`).
Chains run one after another in one call; windows of `mipca_flank()` and
ranks of `mipca_tune()` run one after another in R. The only parallelism is
inside R's BLAS and LAPACK. From C++ it is in the calls listed above, the
Gram's `dsyrk` and the factor's `dgemm` above all; from R it is in the fit
step of `tf_biwhiten()`, through `tcrossprod()`, `eigen()` and `%*%`. Under a
reference BLAS there is none.

## Test entry points

The smokes run each door once (`tests/testthat/test-smoke.R`:
`smoke: mipca() with tf_biwhiten(tf_sv()) and spectra_opts() fits its transform and draws`,
`smoke: mipca_finalize() thins a kept fit with tf_identity()`,
`smoke: mipca_tune() scores each rank and mipca_fit() fits the best`). The
local tier is the `tests/testthat/test-mipca-*.R` files, among them the
door defaults (`tests/testthat/test-mipca-door-defaults.R`:
`mipca()'s sampler defaults are mipca_core()'s`), admissibility
(`tests/testthat/test-mipca-admissible.R`:
`a fully missing column is refused, not silently returned as NaN`), the
lifecycle, the joins (`tests/testthat/test-mipca-narrow-join.R`:
`a request with nothing to impute returns a fit, and it joins`) and the
columns the flank plan reads (`tests/testthat/test-mipca-flank.R`:
`mipca_flank_plan() reads only the columns a window can hold`). Run them
under `load_all1(timer = FALSE)`.
