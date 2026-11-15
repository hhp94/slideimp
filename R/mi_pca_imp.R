# mi_pca_imp.R --------------------------------------------------------------
# Bayesian multiple imputation under a low-rank Gaussian model, in layers:
#   1. mipca_core()      the DA sampler; loop in C++ (src/mipca.cpp)
#   2. mipca()           the entry point: transform, core, fit object
#   3. si_mi_pca / si_mi the fit object, in R/mi_pca_s3.R
# plus the transforms (tf_*) and free scorers crps_all() / pit_all().
# ------------------------------------------------------------------------------

# ---- Layer 1: core sampler ---------------------------------------------------

# The Gram route is the only eigensolve route, and the Gram squares the
# condition number: an eigenvalue ratio above this stops the run.
.GRAM_MAX_COND <- 1e-8 / .Machine$double.eps

# Read by mipca()'s warning, mipca_flank() and mipca_fit()'s clamp: below
# it, the I-step's t is heavy enough to destabilise the Gaussian P-step.
.NU_FLOOR <- 5

# What ncv = NULL means lives in C++, because it needs the rank and the Gram
# side.
#' @export
spectra_opts <- function(tol = 1e-10, ncv = NULL, maxitr = 1000) {
  tol <- as.double(tol)
  if (length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive number.", call. = FALSE)
  }
  maxitr <- .chk_count(maxitr, "maxitr")
  ## 0 is the C++ wire sentinel for "apply the rule". NULL is mapped here,
  ## not in the predicate, because no other count accepts it.
  if (is.null(ncv)) {
    ncv <- 0L
  } else {
    ncv <- .chk_count(ncv, "ncv", note = "or NULL")
  }
  ## Classed so a hand-built list can't pass as spectra_opts() output.
  structure(
    list(tol = tol, ncv = ncv, maxitr = maxitr),
    class = "spectra_opts"
  )
}

# The seam format's version, one owner: .seam_new() writes it and
# .seam_check_version() refuses anything else.
.SEAM_VERSION <- 9L

# The one owner of the seam-version refusal; update.si_mi_pca() and
# mipca_core() both call it.
.seam_check_version <- function(.state) {
  if (
    !inherits(.state, "si_mi_pca_seam") ||
      !identical(.state$seam_version, .SEAM_VERSION)
  ) {
    stop(
      "update() cannot extend this fit: its saved sampler state is not in ",
      "the format this version reads. Fit again with mipca() and ",
      "`keep_state = TRUE`.",
      call. = FALSE
    )
  }
  invisible(.state)
}

# The text of every refusal that finds the matrix a fit carries changed since
# its chains ran. `why` is one sentence saying what changed. The other
# resume checks guard what only mipca_core()'s callers pass, so they are
# internal errors.
.resume_refused <- function(why) {
  paste0(
    "update() cannot extend this fit: the matrix it carries no longer ",
    "matches the one its chains ran on. ",
    why,
    " Fit again with mipca() ",
    "and `keep_state = TRUE`."
  )
}

# The one constructor of the resume seam. Sentinels C++ reads:
# g_fixed is numeric(0) when not split (C++ tests size() > 0); fp_clean is
# NULL when not split.
.seam_new <- function(
  chains,
  iter,
  rank,
  thin,
  spectra,
  nu,
  n,
  cnt,
  mi_row,
  g_fixed,
  fp_clean,
  forward,
  store,
  chain
) {
  structure(
    list(
      seam_version = .SEAM_VERSION,
      chains = chains,
      iter = iter,
      rank = rank,
      thin = thin,
      spectra = spectra,
      nu = nu,
      n = n,
      cnt = cnt,
      mi_row = mi_row,
      g_fixed = g_fixed,
      fp_clean = fp_clean,
      forward = forward,
      store = store,
      chain = chain
    ),
    class = "si_mi_pca_seam"
  )
}

# Refused, not repaired: a wrong position is a wrong CELL, paired with `map`
# by index alone.
.chk_store <- function(store, nmiss) {
  s <- as.integer(store)
  if (length(s) < 1L || anyNA(s)) {
    stop(
      "`store` must be a non-empty integer vector with no NA.",
      call. = FALSE
    )
  }
  ## min/max rather than the first and last element: the sortedness check
  ## below has not run yet, so reading the ends would be assuming its result.
  if (min(s) < 1L || max(s) > nmiss) {
    stop(
      "`store` positions must lie in 1..",
      nmiss,
      " (the matrix has ",
      nmiss,
      " holes).",
      call. = FALSE
    )
  }
  if (length(s) > 1L && min(diff(s)) <= 0L) {
    stop(
      "`store` must be strictly increasing (sorted, no duplicates).",
      call. = FALSE
    )
  }
  s
}

# PASS 1 for fresh and resume alike: one build call, one record. A resume
# takes `forward` and g_fixed from the seam: a dsyrk is not bitwise
# reproducible, so the seam's Gram is the one that continues the run.
.core_build <- function(Z, cnt, .state, forward) {
  forward_used <- if (is.null(.state)) forward else .state$forward
  built <- .mipca_build_stream(Z, cnt, forward_used, resume = !is.null(.state))
  out <- list(
    Xd = built$Xd,
    ## Returned so the seam stores the map actually applied.
    forward_used = forward_used,
    etd = built$et[built$dirty],
    col_ss = built$col_ss,
    mi_row = built$row,
    col_ptr = built$col_ptr,
    dirty = built$dirty,
    g_fixed = if (is.null(.state)) built$g_fixed else .state$g_fixed
  )
  if (!is.null(.state) && !identical(out$mi_row, .state$mi_row)) {
    stop(.resume_refused("Its missing cells have moved."), call. = FALSE)
  }
  out
}

# The fresh start state, the same seven fields as .start_resume(); each chain
# draws its own start in C++, and its first eigensolve starts cold.
# `disp` is narrowed to the block because C++ indexes it by `mi_col_l`.
# ALIASING: chain 1's panel IS `Xd`, not a copy, and C++ mutates it in place.
.start_fresh <- function(
  Xd,
  etd,
  col_ss,
  cnt,
  dirty,
  nmiss,
  n,
  p,
  init_c,
  chains
) {
  nobs <- n - cnt
  if (min(nobs) < 2L) {
    stop(
      "internal error: a column with under two observed cells reached ",
      "the init."
    )
  }
  df_obs <- as.double(n) * p - nmiss - p
  if (df_obs <= 0) {
    ## Each column has two observed cells, so df_obs >= p.
    stop("internal error: too few observed cells to form a starting sigma.")
  }
  list(
    xd_list = lapply(
      seq_len(chains),
      function(ch) if (ch == 1L) Xd else Xd + 0
    ),
    etd_list = rep(list(etd), chains),
    xtild_list = rep(list(numeric(0)), chains),
    eig_start_list = rep(list(numeric(0)), chains),
    sigma_in = rep(sqrt(sum(col_ss) / df_obs), chains),
    rng_state_in = character(0),
    disp = (init_c * sqrt(col_ss / (nobs - 1L)))[dirty]
  )
}

# The resumed start state, the same seven fields as .start_fresh(). Empty
# `disp` is the resume half of a pairing C++ checks against its `resuming`
# flag. `eig_start` is the vector each chain's next eigensolve starts from,
# passed on as the bytes the run left.
.start_resume <- function(
  .state,
  Xd,
  etd,
  cnt,
  dirty,
  miss_idx_l,
  n,
  nms = NULL
) {
  cnt_d <- cnt[dirty]
  list(
    xd_list = lapply(seq_along(.state$chain), function(ch) {
      st_ch <- .state$chain[[ch]]
      ## `etd` is the REBUILT BASE's; st_ch$etd is this chain's accumulated
      ## one. The difference recentres the base onto the chain.
      Xc <- Xd + rep(etd - st_ch$etd, each = n)
      Xc[miss_idx_l] <- st_ch$holes
      .seam_fp_check(st_ch$fp, Xc, cnt_d, dirty, nms)
      Xc
    }),
    etd_list = lapply(.state$chain, `[[`, "etd"),
    xtild_list = lapply(.state$chain, `[[`, "xtild_miss"),
    eig_start_list = lapply(.state$chain, `[[`, "eig_start"),
    sigma_in = vapply(.state$chain, `[[`, 0, "sigma"),
    rng_state_in = vapply(.state$chain, `[[`, "", "rng"),
    disp = numeric(0)
  )
}

# Arguments a resume refuses. `chains_given` is `!missing(chains)` from
# the core's frame. Returns the chain count, which comes from the seam.
.seam_check_args <- function(.state, chains, chains_given, forward, store) {
  .seam_check_version(.state)
  if (chains_given && as.integer(chains) != .state$chains) {
    stop(
      "internal error: chains disagrees with the seam; a resume continues ",
      "every chain of the fit."
    )
  }
  if (!identical(forward, identity)) {
    stop(
      "internal error: `forward` cannot be given on a resume: the seam ",
      "carries the map the fit was built under, and the rebuild applies ",
      "that one."
    )
  }
  if (!is.function(.state$forward)) {
    stop("internal error: the seam carries no transform map.")
  }
  if (!is.null(store)) {
    stop(
      "internal error: `store` cannot be given on a resume: the seam ",
      "carries the store set the fit was built under, and an extension ",
      "keeps those same cells."
    )
  }
  .state$chains
}

# Matrix and settings against the seam, before anything at matrix scale is
# built.
.seam_check_match <- function(.state, n, cnt, rank, nu, spectra, Z) {
  if (.state$n != n || !identical(.state$cnt, cnt)) {
    stop(
      .resume_refused(
        "Its row count or its count of missing cells per column has changed."
      ),
      call. = FALSE
    )
  }
  ## update() passes the seam's own rank, nu and spectra.
  if (.state$rank != rank) {
    stop("internal error: the seam was written at a different rank.")
  }
  if (.state$nu != nu) {
    stop("internal error: the seam was written with a different nu.")
  }
  if (!identical(.state$spectra, spectra)) {
    stop(
      "internal error: the seam was written with different spectra settings."
    )
  }
  ## LAST, and after the census: the clean enumeration this reads is only
  ## meaningful once `cnt` is known to match. The dirty half's fingerprint
  ## is checked later still, per chain, once its panel exists.
  .seam_fp_clean_check(.state, Z, n, cnt)
  invisible(.state)
}

# The dirty half's fingerprint columns: a stride, not a sample, so it draws
# no RNG; the indices ride in the seam.
.seam_fp_cols <- function(nd) {
  m <- min(nd, max(100L, as.integer(ceiling(0.05 * nd))))
  unique(as.integer(round(seq.int(1L, nd, length.out = m))))
}

# Sums of ABSOLUTE values: the sweep recentres every column, so raw column
# sums sit near zero and a relative error on them is meaningless.
.seam_fp_take <- function(Xd, loc, cnt_d) {
  list(
    col = loc,
    sumabs = colSums(abs(Xd[, loc, drop = FALSE])),
    holes = cnt_d[loc]
  )
}

# Relative. A rebuilt panel drifts by rounding, so the fingerprint is
# checked on a tolerance that passes that drift, never on identical().
.SEAM_FP_TOL <- 1e-8

# `nms` is the fitted matrix's column names (NULL: none), so a refusal names
# a column the caller can find; `dirty` alone is a position in that matrix,
# which is the kept columns when the fit took a list.
.seam_fp_check <- function(fp, Xd, cnt_d, dirty, nms = NULL) {
  lab <- function(k) if (is.null(nms)) dirty[k] else nms[dirty[k]]
  cur <- .seam_fp_take(Xd, fp$col, cnt_d)
  if (!identical(cur$holes, fp$holes)) {
    bad <- lab(fp$col[cur$holes != fp$holes])
    stop(
      .resume_refused(paste0(
        "Its missing cells have changed in ",
        length(bad),
        " column",
        .s(length(bad)),
        " [",
        fmt_trunc(bad, 8),
        "]."
      )),
      call. = FALSE
    )
  }
  if (!all(is.finite(fp$sumabs)) || any(fp$sumabs <= 0)) {
    stop("internal error: the seam's fingerprint is not positive and finite.")
  }
  rel <- abs(cur$sumabs - fp$sumabs) / fp$sumabs
  if (max(rel) >= .SEAM_FP_TOL) {
    worst <- which.max(rel)
    ## Missing cells in place, so the values or the transform moved.
    stop(
      .resume_refused(paste0(
        "Its values, or what its transform makes of them, have changed. ",
        "Column ",
        lab(fp$col[worst]),
        " differs most."
      )),
      call. = FALSE
    )
  }
  invisible(NULL)
}

# The clean half's fingerprint: one per seam, on RAW values, by absolute
# column index, over every clean column (a raw read, so cheap).
.seam_fp_clean_take <- function(Z, clean) {
  if (length(clean) == 0L) {
    return(NULL)
  }
  loc <- clean
  ## Zero for an all-zero clean column (constant columns are admitted); the
  ## check divides by pmax(sumabs, double.xmin) so it still compares exactly.
  sa <- colSums(abs(Z[, loc, drop = FALSE]))
  if (!all(is.finite(sa))) {
    stop("internal error: a clean column's fingerprint is not finite.")
  }
  list(col = loc, sumabs = sa)
}

# Same statistic and tolerance as the dirty half. `fp_clean` is NULL exactly
# when not split, so a split seam without one is refused, never skipped.
.seam_fp_clean_check <- function(.state, Z, n, cnt) {
  split <- .split_cols(n, length(cnt), cnt)$split
  if (!split) {
    ## Not split: the per-chain fingerprint already covers every column.
    return(invisible(NULL))
  }
  fp <- .state$fp_clean
  if (
    !is.list(fp) ||
      !is.numeric(fp$sumabs) ||
      !is.integer(fp$col) ||
      length(fp$col) != length(fp$sumabs) ||
      length(fp$col) == 0L
  ) {
    stop(
      "internal error: the seam carries no clean-column fingerprint but ",
      "was built from a split matrix."
    )
  }
  if (max(fp$col) > length(cnt) || any(cnt[fp$col] != 0L)) {
    stop(
      "internal error: the clean fingerprint names a column that is not ",
      "clean in this census."
    )
  }
  cur <- colSums(abs(Z[, fp$col, drop = FALSE]))
  rel <- abs(cur - fp$sumabs) / pmax(fp$sumabs, .Machine$double.xmin)
  if (!all(is.finite(rel)) || max(rel) >= .SEAM_FP_TOL) {
    worst <- which.max(rel)
    ## A resume keeps the clean half's Gram from the seam and never rereads
    ## those columns, so continuing would sample against data this matrix no
    ## longer holds.
    stop(
      .resume_refused(paste0(
        "Its values have changed. Column ",
        .col_label(Z, fp$col[worst]),
        " differs most."
      )),
      call. = FALSE
    )
  }
  invisible(NULL)
}

# Split rule: wide, with at least one clean column; else nd = p and no Gram.
# The which() is one-time DISCOVERY, not a position lookup.
# One owner: the core's resume branch needs the same enumeration.
.split_cols <- function(n, p, cnt) {
  split <- n <= p && any(cnt == 0L)
  list(split = split, dirty = if (split) which(cnt > 0L) else seq_len(p))
}

# The streaming build: block by block straight into the loop's inputs, with
# no whole-matrix work copy. IN PLACE: .build_block() writes Xd, g_fixed, et,
# col_ss, mi_row and the scratch panel; it never writes `X`, which on the
# identity path is the caller's matrix. Only identical(forward, identity)
# takes that path.
.mipca_build_stream <- function(
  X,
  cnt,
  forward = identity,
  block = 1024L,
  resume = FALSE
) {
  n <- nrow(X)
  p <- ncol(X)
  block <- max(1L, as.integer(block))
  ## Integer asserted, not coerced: goes to C++ as an IntegerVector.
  if (
    !is.integer(cnt) ||
      length(cnt) != p ||
      anyNA(cnt) ||
      any(cnt < 0L) ||
      any(cnt > n)
  ) {
    stop("internal error: the census does not describe this matrix.")
  }
  in_place <- identical(forward, identity)
  ## Coerce an integer matrix here once, rather than let Rcpp copy per block.
  if (in_place && !is.double(X)) {
    storage.mode(X) <- "double"
  }

  sp <- .split_cols(n, p, cnt)
  split <- sp$split
  dirty <- sp$dirty
  rm(sp)
  nd <- length(dirty)
  dirty_pos <- integer(p)
  dirty_pos[dirty] <- seq_len(nd)
  col_ptr <- c(1L, 1L + cumsum(cnt))
  nmiss <- col_ptr[p + 1L] - 1L

  ## The loop's inputs, allocated once at final size; C++ writes into them
  ## in place. Nothing else at matrix scale is allocated here.
  Xd <- matrix(0, n, nd)
  ## No Gram and no dsyrk panel on a resume: C++ skips every clean column.
  g_fixed <- if (split && !resume) matrix(0, n, n) else numeric(0)
  et <- numeric(p)
  col_ss <- numeric(p)
  mi_row <- integer(nmiss)
  ## Clean columns' centered panel for dsyrk, reused across blocks.
  scratch <- if (split && !resume) numeric(n * min(block, p)) else numeric(0)

  j1 <- 1L
  while (j1 <= p) {
    j2 <- min(j1 + block - 1L, p)
    nb <- j2 - j1 + 1L
    if (in_place) {
      .build_block(
        X,
        j1 - 1L,
        j1,
        nb,
        0L,
        scratch,
        Xd,
        g_fixed,
        et,
        col_ss,
        mi_row,
        cnt,
        col_ptr,
        dirty_pos,
        if (resume) 1L else 0L
      )
    } else {
      blk <- forward(X[, j1:j2, drop = FALSE], j1:j2)
      if (!is.double(blk)) {
        storage.mode(blk) <- "double"
      }
      .build_block(
        blk,
        0L,
        j1,
        nb,
        1L,
        scratch,
        Xd,
        g_fixed,
        et,
        col_ss,
        mi_row,
        cnt,
        col_ptr,
        dirty_pos,
        if (resume) 1L else 0L
      )
    }
    j1 <- j2 + 1L
  }
  rm(scratch)
  ## Otherwise `g_fixed` stays numeric(0), not NULL: C++ tests size() > 0 at
  ## both entries, so empty IS the no-Gram spelling.
  if (split && !resume) {
    .sym_from_lower(g_fixed)
  }

  list(
    Xd = Xd,
    g_fixed = g_fixed,
    et = et,
    col_ss = col_ss,
    row = mi_row,
    col_ptr = col_ptr,
    dirty = dirty
  )
}

# Every refusal that names columns or cells truncates the list with
# fmt_trunc(), at 8.

# `subset` -> integer column indices, from names or positions, through
# resolve_subset(). `sort = TRUE` is load-bearing: the map is column-major
# and the store index assumes an ascending subset. The refusals run first so
# resolve_subset()'s cli branches are unreachable and every refusal is a
# plain stop().
.resolve_cols <- function(subset, obj, what = "subset") {
  if (is.null(subset)) {
    return(NULL)
  }
  ## `what` names the caller's argument: `subset`, or the `abb` list, which
  ## resolves the same way.
  arg <- paste0("`", what, "`")
  if (length(subset) == 0L) {
    ## An empty LIST is no list: abb_cols() returns integer(0) on a matrix
    ## with nothing to mark, and that fit is the one `abb = NULL` gives. An
    ## empty REQUEST stays an error.
    if (identical(what, "abb")) {
      return(NULL)
    }
    stop(
      arg,
      " selects no columns. Name at least one column of `obj`.",
      call. = FALSE
    )
  }
  if (is.character(subset)) {
    cn <- colnames(obj)
    if (is.null(cn)) {
      stop(
        arg,
        " names columns, but `obj` has no column names. Set ",
        "colnames(obj), or give column numbers.",
        call. = FALSE
      )
    }
    dup <- unique(subset[duplicated(subset)])
    if (length(dup) > 0L) {
      stop(
        arg,
        " names a column more than once: ",
        fmt_trunc(dup, 8),
        ".",
        call. = FALSE
      )
    }
    ## Refused, not dropped: a requested column decides what the fit stores
    ## and returns, and a listed one leaves the sampler.
    bad <- subset[is.na(match(subset, cn))]
    if (length(bad) > 0L) {
      stop(
        arg,
        " names ",
        length(bad),
        " column",
        .s(length(bad)),
        " that ",
        if (length(bad) == 1L) "is" else "are",
        " not in `obj`: ",
        fmt_trunc(bad, 8),
        ".",
        call. = FALSE
      )
    }
  } else {
    ## Same assertion resolve_subset() makes, run first so the refusal is
    ## checkmate's plain condition and not the cli one behind it.
    checkmate::assert_integerish(
      subset,
      lower = 1L,
      upper = ncol(obj),
      any.missing = FALSE,
      unique = TRUE,
      .var.name = what
    )
  }
  got <- resolve_subset(subset, obj, sort = TRUE)
  ## Guards resolve_subset()'s contract (its empty-selection path), not a
  ## caller error.
  if (!is.integer(got) || length(got) != length(subset)) {
    stop("internal error: ", arg, " did not resolve to one index per request.")
  }
  got
}

# Every mat_miss() call in this file goes through here. Asserted, not
# coerced, because C++ reads the census as truth; unnamed, because it is
# positional.
.mat_miss_int <- function(Z) {
  m <- mat_miss(Z, col = TRUE, prop = FALSE)
  if (!is.integer(m) || length(m) != ncol(Z) || anyNA(m)) {
    stop(
      "internal error: mat_miss() did not return one integer count per ",
      "column."
    )
  }
  unname(m)
}

# The admissibility predicate; the response is .check_admissible(), apart
# because the flank driver answers the same classes differently. The sampler
# refuses only what it cannot compute, so there is no `colmax`. `const` is
# reported, not refused: the tuner's holdout draw and the flank's windows
# read it.
.admissible_census <- function(Z) {
  n <- nrow(Z)
  cmiss <- .mat_miss_int(Z)
  v <- col_vars(Z)

  ## Not !is.finite(): that would also sweep in col_vars()'s NA_real_ for
  ## "fewer than two observed cells", a different class.
  nonfin <- is.nan(v) | is.infinite(v)

  empty <- cmiss == n
  single <- cmiss == n - 1L
  if (!identical(unname(is.na(v) & !nonfin), unname(empty | single))) {
    stop(
      "internal error: col_vars() NA set disagrees with the observed counts."
    )
  }
  fin <- !empty & !single & !nonfin
  const <- fin & v < .Machine$double.eps
  list(
    cmiss = cmiss,
    nonfin = nonfin,
    empty = empty,
    single = single,
    const = const
  )
}

# What a refusal calls columns `k` of `Z`: their names, which are right in
# any numbering -- `Z` is the kept columns when a list is in play, and a
# position among those is not one the caller can find. Every door's matrix
# has unique names (.chk_matrix()); a bare index is the fallback for an
# unnamed matrix handed straight to the core.
.col_label <- function(Z, k) {
  nms <- colnames(Z)
  if (is.null(nms)) k else nms[k]
}

# `lead` and `fix` let the tuner say that its holdout, not the matrix, is
# what left a column short.
.check_admissible <- function(
  Z,
  lead = "Some columns of `obj` cannot be imputed.",
  fix = paste0(
    "Remove these columns from `obj`. A column with one observed cell can ",
    "be listed in `abb` instead."
  )
) {
  cen <- .admissible_census(Z)
  nonfin <- cen$nonfin
  empty <- cen$empty
  single <- cen$single

  ## An admissible matrix returns its census: the core takes `cmiss` as
  ## `cnt`, and the tuner's door reads `const`, off the same scan.
  if (!any(nonfin | empty | single)) {
    return(invisible(cen))
  }

  show <- function(mask, what) {
    k <- which(mask)
    if (length(k) == 0L) {
      return(NULL)
    }
    paste0(
      "  ",
      what,
      ": ",
      length(k),
      " column",
      .s(length(k)),
      " [",
      fmt_trunc(.col_label(Z, k), 8),
      "]"
    )
  }

  ## Two observed cells, because each chain starts spread by its column's
  ## observed sd.
  stop(
    lead,
    " Each column needs at least two observed cells and a finite ",
    "variance.\n",
    paste(
      c(
        show(nonfin, "non-finite variance (an Inf cell, or overflow)"),
        show(empty, "no observed cells"),
        show(single, "only one observed cell")
      ),
      collapse = "\n"
    ),
    "\n",
    fix,
    call. = FALSE
  )
}

mipca_core <- function(
  Z,
  rank,
  subset = NULL,
  ## What a fit keeps: a cell-level set of positions in the hole enumeration.
  store = NULL,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  ## Chains run in turn. No BLAS pin: the BLAS keeps whatever threads it has.
  chains = 4L,
  refresh = 0L,
  spectra = spectra_opts(),
  ## Every chain starts itself; `init` is a dispersion multiple.
  init = 1,
  forward = identity,
  ## Applied only when `seed_given`, after every refusal; a NULL seed is one
  ## draw off R's stream, so set.seed() governs the sampler. A resume never
  ## seeds.
  seed = NULL,
  seed_given = FALSE,
  keep_state = FALSE,
  .state = NULL # internal: a keep_state seam; update.si_mi_pca() only
) {
  spectra <- .chk_spectra(spectra)
  if (!is.function(forward)) {
    stop("`forward` must be a function.", call. = FALSE)
  }
  init_c <- .chk_init(init)
  nu <- .chk_nu(nu)
  ## A warning, not a refusal, fresh runs only. The flank driver muffles
  ## this copy by its text, so it must stay .nu_floor_msg().
  if (is.null(.state) && nu < .NU_FLOOR) {
    warning(.nu_floor_msg(), call. = FALSE)
  }
  ## Chain c runs on dqrng stream c - 1, so chain 1 is bitwise the serial run.
  if (!is.null(.state)) {
    chains <- .seam_check_args(
      .state,
      chains,
      !missing(chains),
      forward,
      store
    )
  }
  chains <- .chk_count(chains, "chains")
  refresh <- .chk_refresh(refresh)
  ## Same predicate and `what` strings as mipca_tune()'s door; keep them in
  ## sync. Floors: warmup may be 0, which every resume passes; sampling and
  ## thin are at least 1.
  iter_warmup <- .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  thin <- as.double(.chk_count(thin, "thin"))
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("Z must be a numeric matrix.", call. = FALSE)
  }
  ## No pass over the values here: the census names a +-Inf column and the
  ## build refuses one per cell.

  n <- nrow(Z)
  p <- ncol(Z)
  ## The one rank ceiling (.chk_rank()), before the build, so a refused
  ## shape pays nothing.
  S <- .chk_rank(rank, n, p)

  K <- min(n - 1L, p)
  ## In double from the first product: integer n * p is NA past 2^31 cells.
  dfP <- as.double(n) * p - (p + as.double(S) * (n - 1L + p - S))
  ## At the largest rank .chk_rank() passes, dfP is p - n + 2 when n <= p
  ## and 2 * (n - p) + 2 when n > p, so it is never below 2.
  if (dfP <= 0) {
    stop("internal error: non-positive residual df.")
  }

  ## ---- admissibility ----------------------------------------------------------
  ## PASS 0: the census becomes `cnt`. On a resume the verdict is skipped and
  ## the census must equal the seam's instead.
  cnt <- if (is.null(.state)) {
    .check_admissible(Z)$cmiss
  } else {
    .mat_miss_int(Z)
  }
  ## Reached from mipca() only with `keep_state = TRUE`; without it a matrix
  ## with nothing to impute returns an empty fit before the core.
  if (sum(cnt) == 0L) {
    stop(
      "`obj` has no missing cells, so there is nothing to resume. Drop ",
      "`keep_state = TRUE`.",
      call. = FALSE
    )
  }

  ## ---- resume: the wrong-matrix guard, before anything is built ---------------
  if (!is.null(.state)) {
    .seam_check_match(.state, n, cnt, S, nu, spectra, Z)
  }

  ## ---- pass 1: the streaming build ---------------------------------------------
  ## `Xd` is private to this call: the P-step sweeps it in place.
  bld <- .core_build(Z, cnt, .state, forward)
  ## Unpacked and the record dropped, so later rm() calls can release each
  ## piece once consumed.
  Xd <- bld$Xd
  etd <- bld$etd
  col_ss <- bld$col_ss
  mi_row <- bld$mi_row
  col_ptr <- bld$col_ptr
  dirty <- bld$dirty
  g_fixed <- bld$g_fixed
  forward_used <- bld$forward_used
  rm(bld)
  nmiss <- length(mi_row)
  nd <- length(dirty)
  mi_col <- rep(seq_len(p), cnt)
  ## Block-local column of every hole, by scatter table; a surviving 0 means
  ## the construction broke.
  inv <- integer(p)
  inv[dirty] <- seq_along(dirty)
  mi_col_l <- inv[mi_col]
  if (min(mi_col_l) < 1L) {
    stop("internal error: a missing cell maps to a clean column.")
  }
  rm(inv)
  ## Linear indices in the block's own coordinates.
  miss_idx_l <- mi_row + (mi_col_l - 1L) * n

  ## ---- storage subsetting ---------------------------------------------------
  ## The holed set is cnt > 0L, NOT `dirty`, which the non-split branch
  ## widens to every column. `subset` returns as asked, `cols` as narrowed.
  if (!is.null(store) && !is.null(subset)) {
    stop(
      "`subset` and `store` cannot both be given: `subset` stores every hole ",
      "in the columns it names and `store` names the cells outright.",
      call. = FALSE
    )
  }
  cols_req <- .resolve_cols(subset, Z)
  if (is.null(cols_req)) {
    cols_req <- seq_len(p)
  }
  ## The store set, cell-level. A resume takes it from the seam.
  store_use <- if (is.null(.state)) store else .state$store
  if (is.null(store_use)) {
    ## A mask on `cnt`, not intersect(): positional, and `cols_req` stays
    ## sorted.
    cols_eff <- cols_req[cnt[cols_req] > 0L]
    ## The enumeration groups holes by column, so each retained column is one
    ## contiguous run and the full-width case is exactly seq_len(nmiss).
    store_idx <- sequence(cnt[cols_eff], from = col_ptr[cols_eff])
  } else {
    ## Against THIS matrix's hole count; on resume the census guard has
    ## already required the seam's matrix.
    store_idx <- .chk_store(store_use, nmiss)
    ## unique() on an ascending vector dedupes runs; it is not a search.
    cols_eff <- unique(mi_col[store_idx])
  }
  mi_col_store <- mi_col_l[store_idx]
  map <- cbind(row = mi_row[store_idx], col = mi_col[store_idx])

  ## ---- DA loop ---------------------------------------------------------------
  niter <- iter_warmup + iter_sampling * thin
  ## Cumulative across resumes: the run's whole length.
  iter_total <- (if (is.null(.state)) 0L else .state$iter) + niter

  ## SEEDED HERE, after the last refusal and before the first draw. The
  ## fresh start's two guards sit past this point, unreachable behind the
  ## census. C++ mutates each chain's Xd and draws in place.
  if (seed_given) {
    ## Resolved to a number first, so the seed this run used can be returned.
    seed <- .resolve_seed(seed)
    dq <- .seed_dq(seed)
    on.exit(dqrng::dqrng_set_state(dq), add = TRUE)
  }
  ## Fresh-vs-resume is answered HERE, once. The two constructors return the
  ## same seven fields, so nothing below this line branches on it.
  st <- if (is.null(.state)) {
    .start_fresh(Xd, etd, col_ss, cnt, dirty, nmiss, n, p, init_c, chains)
  } else {
    .start_resume(.state, Xd, etd, cnt, dirty, miss_idx_l, n, nms = colnames(Z))
  }
  xd_list <- st$xd_list
  etd_list <- st$etd_list
  xtild_list <- st$xtild_list
  eig_start_list <- st$eig_start_list
  sigma_in <- st$sigma_in
  rng_state_in <- st$rng_state_in
  disp <- st$disp
  ## Nothing else reads the base and its column means. Dropped because on a
  ## fresh run chain 1's panel IS `Xd`, and a second name on a buffer C++ is
  ## about to mutate must not survive.
  rm(st, Xd, etd, col_ss)
  draws_list <- .alloc_draws(
    lapply(
      seq_len(chains),
      function(ch) matrix(NA_real_, length(store_idx), iter_sampling)
    ),
    cells = length(store_idx),
    chains = chains,
    draws = iter_sampling,
    fix = "Keep fewer cells with `subset`, or ask for fewer `chains` or `iter_sampling`."
  )
  res <- .mipca_da_loop_multi(
    xd_list = xd_list,
    g_fixed = g_fixed,
    etd_list = etd_list,
    xtild_list = xtild_list,
    eig_start_list = eig_start_list,
    sigma_in = sigma_in,
    disp = disp,
    ## Stated, not inferred from the inputs: C++ checks the emptiness of
    ## rng_state_in, disp, xtild and eig_start as pairing guards on this flag.
    resuming = !is.null(.state),
    rng_state_in = rng_state_in,
    miss_idx_l = miss_idx_l,
    mi_row = mi_row,
    mi_col_l = mi_col_l,
    store_idx = store_idx,
    mi_col_store = mi_col_store,
    draws_list = draws_list,
    streams = seq_len(chains) - 1L,
    S = S,
    K = K,
    p = p,
    dfP = as.double(dfP),
    nu = nu,
    ## `iter_*` is cmdstanr's R-facing spelling; the C++ entry keeps its own
    ## warmup/ndraws names.
    warmup = as.integer(iter_warmup),
    ndraws = as.integer(iter_sampling),
    thin = as.integer(thin),
    gram_max_cond = .GRAM_MAX_COND,
    spectra_tol = spectra$tol,
    spectra_ncv = spectra$ncv,
    spectra_maxitr = spectra$maxitr,
    refresh = refresh
  )
  ## Every consumer assumes one row of draws per row of `map`; a mismatch
  ## would read as plausible numbers against wrong cells.
  if (!identical(nrow(draws_list[[1L]]), nrow(map))) {
    stop("internal error: draws rows do not match the map's rows.")
  }
  ## `draws` is ALWAYS the per-chain list, chains == 1 included.
  out <- list(
    draws = draws_list,
    map = map,
    cols = cols_eff,
    subset = cols_req,
    chains = chains,
    iter = iter_total,
    rank = S,
    thin = as.integer(thin),
    nu = nu,
    ## Holes in the WHOLE matrix, against nrow(map)'s stored ones: the two
    ## differ when a partial store set means this fit cannot realize.
    nmiss = nmiss,
    ## The seed this run was given or drew; NULL when it took the ambient
    ## dqrng state. The ABB draws of a listed column are set from it.
    seed = if (seed_given) as.integer(seed) else NULL
  )
  if (keep_state) {
    ## The compact seam (see .seam_new()). The fingerprint is taken on the
    ## same panel the hole values leave.
    fp_loc <- .seam_fp_cols(nd)
    cnt_d_out <- cnt[dirty]
    ## `nd < p` IS the split flag; the which() is the same discovery that
    ## produced `dirty`.
    fp_clean <- if (nd < p) {
      .seam_fp_clean_take(Z, which(cnt == 0L))
    } else {
      NULL
    }
    out$state <- .seam_new(
      chains = chains,
      iter = iter_total,
      rank = S,
      thin = as.integer(thin),
      spectra = spectra,
      nu = nu,
      n = n,
      cnt = cnt,
      mi_row = mi_row,
      g_fixed = g_fixed,
      fp_clean = fp_clean,
      forward = forward_used,
      ## The caller's set, NOT the derived `store_idx`: NULL means "every
      ## hole in `subset`", keeping the default seam small.
      store = store_use,
      chain = lapply(seq_len(chains), function(ch) {
        list(
          holes = xd_list[[ch]][miss_idx_l],
          etd = res$etd[[ch]],
          sigma = res$sigma[[ch]],
          xtild_miss = res$xtild_miss[[ch]],
          eig_start = res$eig_start[[ch]],
          rng = res$rng_state[[ch]],
          fp = .seam_fp_take(xd_list[[ch]], fp_loc, cnt_d_out)
        )
      })
    )
  }
  out
}


# ---- Layer 2: transforms -----------------------------------------------------

# A transform is list(name, forward, backward), or list(name, fit) whose
# fit(X, rank) returns the fitted one. forward(y, cols) maps an n x k
# COLUMN block, `cols` its global indices; backward(z, col, row) maps cells
# back. No door defaults `transform`: the scale the sampler works on is the
# caller's choice.

# The one reader of the transform names. A list is returned untouched, so
# resolving twice (the flank door, then per window) is harmless.
.resolve_transform <- function(transform, n) {
  msg <- paste0(
    "`transform` must be \"identity\", \"logit_biwhiten\", or a transform ",
    "object: list(name, forward, backward), or list(name, fit) for one that ",
    "is fitted first."
  )
  if (is.character(transform)) {
    if (length(transform) != 1L || is.na(transform)) {
      stop(msg, call. = FALSE)
    }
    tf <- switch(
      transform,
      identity = tf_identity(),
      logit_biwhiten = tf_biwhiten(tf_sv(n), rows = TRUE),
      NULL
    )
    if (is.null(tf)) {
      stop(msg, call. = FALSE)
    }
    return(tf)
  }
  ok <- is.list(transform) &&
    is.character(transform$name) &&
    length(transform$name) == 1L &&
    (is.function(transform$fit) ||
      (is.function(transform$forward) && is.function(transform$backward)))
  if (!ok) {
    stop(msg, call. = FALSE)
  }
  transform
}

# The bare `identity` in `forward` is load-bearing: .mipca_build_stream()
# keys its in-place read on identical(forward, identity).
#' @export
tf_identity <- function() {
  list(name = "identity", forward = identity, backward = function(z, ...) z)
}

# Smithson & Verkuilen (2006) squeeze, then logit. `n` is stated, not learned,
# so the map stays stateless. Not an exact round trip at 0/1 (~3.8e-19): do
# not hand clamp_frac() a near-zero relative tolerance on it.
#' @export
tf_sv <- function(n) {
  n <- .chk_count(n, "`n`")
  if (n < 2L) {
    stop("`n` must be at least 2.", call. = FALSE)
  }
  nd <- as.double(n)
  list(
    name = paste0("sv-logit(n=", n, ")"),
    forward = function(y, ...) {
      ## Short-circuited: the fused `y < 0 | y > 1` allocates three logicals.
      if (any(y < 0, na.rm = TRUE) || any(y > 1, na.rm = TRUE)) {
        stop(
          "This transform reads proportions, so every observed cell of ",
          "`obj` must lie in [0, 1]. For data on another scale, use ",
          "`transform = \"identity\"`.",
          call. = FALSE
        )
      }
      stats::qlogis((y * (nd - 1) + 0.5) / nd)
    },
    backward = function(z, ...) {
      pmin(pmax((stats::plogis(z) * nd - 0.5) / (nd - 1), 0), 1)
    }
  )
}

# Biwhitening: tf_biwhiten() and its helpers, after Landa, Zhang & Kluger
# (2022, arXiv:2103.13840). The noise is a scale mixture: a column and a row
# factor, frozen by the fit below, and a cell factor in the sampler. A
# column's scale is its regression residual variance shrunk toward the pool.
# The helpers are limma's trigammaInverse(), fitFDist() (no covariates) and
# squeezeVar() (Smyth 2004, SAGMB 3(1), art 3).
.trigamma_inverse <- function(x) {
  ## Newton on trigamma(y) = x from the large-y asymptote.
  y <- 0.5 + 1 / x
  for (iter in seq_len(50L)) {
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2L)
    y <- y + dif
    if (max(-dif / y) < 1e-8) {
      return(y)
    }
  }
  stop("internal error: the trigamma inverse did not converge.", call. = FALSE)
}

.fit_fdist <- function(s2, df) {
  e <- log(s2) - digamma(df / 2) + log(df / 2)
  emean <- mean(e)
  ## One column has no spread to read a prior off: the pool is that column.
  evar <- if (length(e) > 1L) {
    sum((e - emean)^2) / (length(e) - 1) - mean(trigamma(df / 2))
  } else {
    -Inf
  }
  if (evar > 0) {
    d0 <- 2 * .trigamma_inverse(evar)
    s20 <- exp(emean + digamma(d0 / 2) - log(d0 / 2))
  } else {
    d0 <- Inf
    s20 <- exp(emean)
  }
  list(d0 = d0, s20 = s20)
}

.squeeze_var <- function(s2, df, d0, s20) {
  if (is.infinite(d0)) {
    return(rep(s20, length(s2)))
  }
  (d0 * s20 + df * s2) / (d0 + df)
}

# Residual df every column keeps after its mean and its regressors:
# k_j = min(rank, n_obs_j - 1 - .BIWHITEN_D_MIN) regressors.
.BIWHITEN_D_MIN <- 2L

#' @export
tf_biwhiten <- function(base = tf_identity(), rows = TRUE) {
  if (
    !is.list(base) || !is.function(base$forward) || !is.function(base$backward)
  ) {
    stop(
      "`base` must be a transform: list(name, forward, backward).",
      call. = FALSE
    )
  }
  if (is.function(base$fit)) {
    stop(
      "`base` must be a transform with no `fit` step, such as tf_identity() ",
      "or tf_sv().",
      call. = FALSE
    )
  }
  if (!isTRUE(rows) && !isFALSE(rows)) {
    stop("`rows` must be TRUE or FALSE.", call. = FALSE)
  }
  nm <- paste0(if (rows) "biwhiten(" else "whiten(", base$name, ")")
  list(
    name = nm,
    fit = function(X, rank) {
      rank <- .chk_count(rank, "`rank`")
      ## One argument: tf_identity()'s forward IS base::identity, which takes
      ## exactly one.
      Z <- base$forward(X)
      obs <- !is.na(Z)
      n <- nrow(Z)
      p <- ncol(Z)
      n_obs <- colSums(obs)
      ## pmax: a column with no observed cells has no mean; 0 / 1 = 0 keeps
      ## it finite here and the core's census refuses it by name.
      mu <- colSums(Z, na.rm = TRUE) / pmax(n_obs, 1L)
      ## rep(each =), not sweep(): no full-size array of the statistic.
      Zf <- Z - rep(mu, each = n)
      ## Every pass-local n x p is dropped after its last reader.
      rm(Z)
      Zf[!obs] <- 0 # column-mean fill, on the centred scale
      e <- eigen(tcrossprod(Zf), symmetric = TRUE)
      W <- e$vectors[, seq_len(rank), drop = FALSE]
      rm(e)
      ## Final for fully observed columns (W is orthonormal over all rows);
      ## every other column is overwritten by its own fit below.
      R <- Zf - W %*% crossprod(W, Zf)
      R[!obs] <- 0 # residuals over OBSERVED cells only
      k <- pmax(0L, pmin(rank, n_obs - 1L - .BIWHITEN_D_MIN))
      full <- n_obs == n & k == rank
      ## Leverage accumulated per ROW, never stored n x p; the mean's is
      ## 1 / n_obs_j per observed cell, a full column's basis one rowSums(W^2).
      lev_row <- drop(obs %*% (1 / pmax(n_obs, 1L))) + rowSums(W^2) * sum(full)
      for (j in which(!full)) {
        # one-time enumeration of the columns to solve
        O <- which(obs[, j]) # discovery: nothing else knows where they are
        if (k[j] == 0L) {
          R[O, j] <- Zf[O, j]
          next
        }
        Wo <- W[O, seq_len(k[j]), drop = FALSE]
        Gi <- solve(crossprod(Wo))
        b <- Gi %*% crossprod(Wo, Zf[O, j])
        R[O, j] <- Zf[O, j] - Wo %*% b
        lev_row[O] <- lev_row[O] + rowSums((Wo %*% Gi) * Wo)
      }
      rm(Zf)
      ## Nothing below reads the signed residual.
      R <- R^2
      ## Exact for a fully observed column; on a partly observed one the
      ## exact df is up to one higher, not charged back.
      df <- pmax(n_obs - 1L, 0L) - k
      rss <- colSums(R)
      s2 <- numeric(p)
      pos <- df > 0L
      s2[pos] <- rss[pos] / df[pos]
      ## The prior reads only columns with positive df AND residual (log 0);
      ## every column is then squeezed, df 0 to exactly the pool.
      use <- pos & rss > 0
      pr <- .fit_fdist(s2[use], df[use])
      s <- sqrt(.squeeze_var(s2, df, pr$d0, pr$s20))
      if (!all(is.finite(s)) || any(s <= 0)) {
        stop(
          "tf_biwhiten() could not estimate a finite, positive noise scale ",
          "for every column of `obj`. Use `transform = \"identity\"` instead.",
          call. = FALSE
        )
      }
      d0 <- pr$d0
      s0 <- sqrt(pr$s20)
      t <- rep(1, n)
      if (rows) {
        ## The SAME residual, rescaled by the column pass.
        q <- drop(R %*% (1 / s^2))
        df_i <- rowSums(obs) - lev_row
        ## No df, or a residual of exactly zero: the row stays at the gauge
        ## level 1, outside the geomean.
        has <- df_i > 0 & q > 0
        t[has] <- sqrt(q[has] / df_i[has])
        if (!all(is.finite(t)) || any(t <= 0)) {
          stop(
            "internal error: the row scale is not finite and positive in every row.",
            call. = FALSE
          )
        }
        t[has] <- t[has] / exp(mean(log(t[has]))) # gauge: geomean 1, so `s` keeps the level
      }
      ## The closures capture THIS frame, which serializes with every fit, so
      ## temporaries go by keep-list (all.names: ls() hides dot bindings).
      rm(
        list = setdiff(ls(all.names = TRUE), c("s", "t", "n", "p", "d0", "s0"))
      )
      list(
        name = nm,
        scale = s, # inspectable; the maps below read these same two vectors
        rowscale = t,
        ## Inspectable too, read by nothing: the prior the columns were shrunk
        ## toward.
        d0 = d0,
        s0 = s0,
        forward = function(y, cols) {
          if (length(cols) != ncol(y)) {
            stop(
              "whiten forward: `cols` must index every column of the block.",
              call. = FALSE
            )
          }
          ## Blocks are COLUMN blocks (every row present), which is what lets
          ## the row scale apply here.
          if (NROW(y) != n) {
            stop(
              "whiten forward: the block must carry every row of the fitted matrix.",
              call. = FALSE
            )
          }
          ## Recycling rather than sweep() (see `fit` above); the length-n row
          ## scale recycles down the columns with no rep().
          z <- base$forward(y) / rep(s[cols], each = NROW(y))
          if (rows) {
            z <- z / t
          }
          z
        },
        backward = function(z, col, row = NULL) {
          if (length(col) != NROW(z)) {
            stop(
              "whiten backward: one column index per cell is required.",
              call. = FALSE
            )
          }
          if (any(col < 1L | col > p)) {
            stop(
              "whiten backward: column index outside the fitted matrix.",
              call. = FALSE
            )
          }
          ## REQUIRED, not defaulted: a forgotten `row` scales draws by t_i
          ## silently.
          if (rows) {
            if (is.null(row)) {
              stop(
                "whiten backward: this transform has a row scale; `row` is required.",
                call. = FALSE
              )
            }
            if (length(row) != NROW(z)) {
              stop(
                "whiten backward: one row index per cell is required.",
                call. = FALSE
              )
            }
            if (any(row < 1L | row > n)) {
              stop(
                "whiten backward: row index outside the fitted matrix.",
                call. = FALSE
              )
            }
          }
          ## z has one ROW per cell, so column-major recycling of the per-cell
          ## s[col] * t[row] scales each row.
          base$backward(if (rows) z * (s[col] * t[row]) else z * s[col])
        }
      )
    }
  )
}

# A fraction in (0, 1]. Zero is refused rather than read as "diagnose
# nothing": an empty frame is a frame that cannot fail.
.chk_frac <- function(x, what) {
  v <- as.double(x)
  if (length(v) != 1L || is.na(v) || v <= 0 || v > 1) {
    stop(
      .arg_name(what),
      " must be a single fraction in (0, 1].",
      call. = FALSE
    )
  }
  v
}

# An argument's name as a refusal prints it, in backticks; a name that
# already carries them keeps one pair. Called on the error path only.
.arg_name <- function(what) {
  paste0("`", gsub("`", "", what, fixed = TRUE), "`")
}

# Closes a TERMINAL `...` that S3 consistency forces a method to carry, so
# a typo is refused, not ignored. `...` comes first so only an exact
# `.what =` binds to this frame.
.chk_no_dots <- function(..., .what) {
  if (...length() > 0L) {
    nm <- ...names()
    nm <- if (is.null(nm)) character(0) else nm[nzchar(nm)]
    stop(
      .what,
      " takes no other arguments",
      if (length(nm) > 0L) paste0(" (got ", fmt_trunc(nm, 8), ")") else "",
      ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# One owner for both doors. Row and column names must be unique, so no name
# silently resolves to the first of two matches. Named `obj` because
# checkmate's message reports this frame's name.
.chk_matrix <- function(obj) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    row.names = "unique",
    col.names = "unique"
  )
}

# A single whole number >= lo (optionally <= hi), coerced and returned; a
# value whose integer is not the value is refused, not truncated.
.chk_count <- function(x, what, lo = 1L, hi = NULL, note = NULL) {
  ## `v` (truncated) and `d` (not) must agree.
  num <- is.numeric(x) || is.character(x)
  v <- if (num) suppressWarnings(as.integer(x)) else NA_integer_
  d <- if (num) suppressWarnings(as.double(x)) else NA_real_
  ## `||` short-circuits, which is what keeps the length test ahead of the
  ## scalar ones: a length-2 `x` never reaches is.na() to error there.
  bad <- length(v) != 1L || is.na(v) || is.na(d) || v != d || v < lo
  if (!bad && !is.null(hi)) {
    bad <- v > hi
  }
  if (bad) {
    if (is.null(hi)) {
      stop(
        .arg_name(what),
        " must be a single integer >= ",
        lo,
        if (is.null(note)) "" else paste0(" (", note, ")"),
        ".",
        call. = FALSE
      )
    }
    stop(
      .arg_name(what),
      " must be a single index in [",
      lo,
      ", ",
      hi,
      "].",
      call. = FALSE
    )
  }
  v
}

# THE rank ceiling, one owner; must match its C++ mirror in da_ctx_build().
# The Gram-route bound is Spectra's own: rank + 2 <= min(n, p).
.chk_rank <- function(rank, n, p) {
  S <- .chk_count(rank, "rank")
  if (S + 2L > min(n, p)) {
    ## For the docs: the sampler solves rank + 1 eigenpairs, and the
    ## eigensolver needs that many below the smaller dimension.
    stop(
      "`rank` must be at most the smaller dimension minus two (",
      min(n, p) - 2L,
      " here).",
      call. = FALSE
    )
  }
  S
}

# nu > 2 or the variance-matched t scale is not real. Shared by the core and
# the empty fit.
.chk_nu <- function(nu) {
  if (!is.numeric(nu) || length(nu) != 1L || is.na(nu) || !(nu > 2)) {
    stop(
      "`nu` must be a single number > 2 (Inf for Gaussian noise).",
      call. = FALSE
    )
  }
  as.double(nu)
}

# The below-floor warning's text, one owner: the core raises it, and the
# flank driver raises it once at its door and recognises the core's copy
# per window by this string alone (mipca_flank()).
.nu_floor_msg <- function() {
  paste0(
    "`nu` is below ",
    .NU_FLOOR,
    ", where the sampler can become unstable. ",
    "Use a `nu` of ",
    .NU_FLOOR,
    " or more."
  )
}

# Predicates shared by mipca_core() and the tuner's and flank driver's doors,
# so a door is never looser than the rule it fronts. Each returns the coerced
# value.

# A positive dispersion multiplier, coerced to double.
.chk_init <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`init` must be a single positive number.", call. = FALSE)
  }
  as.double(x)
}

# The flank driver's windowing threshold (the sampler takes no `colmax`).
# NOT coerced, unlike the counts: as.numeric() would silently accept a
# character threshold; this matches pca_imp()'s assert_number acceptance set.
.chk_colmax <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    stop("`colmax` must be a single number in [0, 1].", call. = FALSE)
  }
  x
}

# Coerced FIRST, then checked: it is how a character "3" is accepted and an
# "x" becomes the NA this refuses.
.chk_refresh <- function(x) {
  ## `note` exists for this line: a floor of 0 is the "off" sentinel, which
  ## the constraint alone does not say.
  .chk_count(x, "`refresh`", lo = 0L, note = "0 = no progress lines")
}

.chk_spectra <- function(x) {
  if (!inherits(x, "spectra_opts")) {
    stop("`spectra` must come from spectra_opts().", call. = FALSE)
  }
  x
}

# The seed as a number, one owner: the caller's, checked, or for NULL one
# draw off R's stream -- the tie to set.seed() the sampler, the tuner and the
# ABB draws all share.
.resolve_seed <- function(seed) {
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  .chk_count(seed, "`seed`", lo = 0L)
}

# Seeds dqrng for one call and returns the state found, for the caller's own
# on.exit().
.seed_dq <- function(seed) {
  seed <- .resolve_seed(seed)
  old <- dqrng::dqrng_get_state()
  dqrng::dqset.seed(seed)
  old
}

# (row, col) -> column-major linear index: one owner of the convention every
# leakage guard is built on. unname(): a one-row map's name would survive
# and fail the identical() guards. NA past 2^31 cells, a size the door
# refuses (.chk_cells()).
#
# `to` is a column table: the cell's column is to[col], the column of a
# joined result (.chunk_to()), and the index is a double, so it does not
# overflow where the joined matrix passes 2^31 cells. A column outside the
# table is refused (.to_col()) and so is a row outside 1..n: left to the
# assignment the index feeds, an NA index is taken in silence when one value
# is assigned, and a row past n is a cell of the next column.
.map_lin <- function(map, n, to = NULL) {
  if (is.null(to)) {
    return(unname(map[, "row"] + (map[, "col"] - 1L) * n))
  }
  row <- map[, "row"]
  if (length(row) > 0L && (anyNA(row) || min(row) < 1L || max(row) > n)) {
    stop("internal error: a stored cell's row is outside the matrix.")
  }
  unname(row + (.to_col(to, map[, "col"]) - 1) * as.double(n))
}

# The draws are the one allocation that grows with the request (cells kept x
# chains x draws) rather than with the matrix, and there is no size to refuse
# at in advance: R handles long vectors, so the limit is the machine's
# memory. When R cannot allocate them, its own message names neither the
# object nor the way out, so it is re-raised with both. R's text is kept, so
# an error that is not about memory still shows as itself. `expr` is the
# allocation, evaluated here; `draws` is per chain; `fix` is the next step
# at the calling site.
.alloc_draws <- function(expr, cells, chains, draws, fix) {
  tryCatch(expr, error = function(e) {
    num <- function(x) {
      format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    }
    gb <- 8 * as.double(cells) * chains * draws / 1024^3
    stop(
      "The draws need ",
      format(gb, digits = 3),
      " GB (",
      num(cells),
      " cell",
      .s(cells),
      " x ",
      num(chains),
      " chain",
      .s(chains),
      " x ",
      num(draws),
      " draw",
      .s(draws),
      "), more memory than R could ",
      "allocate. ",
      fix,
      "\nR reported: ",
      conditionMessage(e),
      call. = FALSE
    )
  })
}

# The size one fit takes, refused at the door and before any sampling: the
# hole map's linear index is integer arithmetic (.map_lin()), NA past this.
# The tuner's door makes the same test. The limit is per fit: the joined
# result of cbind() may pass it, and as.matrix() indexes that in double.
.chk_cells <- function(n, p) {
  if (as.double(n) * p > .Machine$integer.max) {
    stop(
      "`obj` has ",
      format(as.double(n) * p, big.mark = ",", scientific = FALSE),
      " cells, more than one fit can index (",
      .Machine$integer.max,
      "). Fit it in parts by column and join them with cbind().",
      call. = FALSE
    )
  }
  invisible(NULL)
}

# The inverse, for a bare linear index: which column it falls in. One owner
# for the reason .map_lin() has one.
.lin_col <- function(lin, n) {
  (lin - 1L) %/% n + 1L
}

# The forward convention, one owner beside the backward one: forward(blk,
# cols) on an n x k COLUMN block. The identity map takes one argument and the
# build stream never calls it (`.mipca_build_stream()`, the in-place path), so
# neither does this; the same identical() is the discriminator in both.
.forward_at <- function(tf, blk, cols) {
  if (identical(tf$forward, identity)) {
    blk
  } else {
    tf$forward(blk, cols)
  }
}

# The back-transform convention, one owner: backward(draws, col, row), col
# FIRST; a swapped pair mis-scales every draw silently. `d` must already be
# the draws matching `rows` (NULL: the whole map).
.backward_at <- function(tf, d, map, rows = NULL) {
  if (is.null(rows)) {
    tf$backward(d, map[, "col"], map[, "row"])
  } else {
    tf$backward(d, map[rows, "col"], map[rows, "row"])
  }
}

# ---- Nothing to impute: the fit that does not sample --------------------------
# A zero-row fit, the identity of the join, returned without sampling; with
# `leave_na`, the flank driver's target carried with its holes NA.

# One-time DISCOVERY. An Inf or an overflowing sum goes down the ordinary
# path, whose census refuses it by name.
.nothing_to_impute <- function(obj, subset) {
  clean <- function(v) !anyNA(v) && (!is.double(v) || is.finite(sum(v)))
  if (is.null(subset)) {
    return(clean(obj))
  }
  for (j in subset) {
    if (!clean(obj[, j])) {
      return(FALSE)
    }
  }
  TRUE
}

# The transform of a fit that never ran, left unfitted. Its own function so
# the closure captures this frame (one string), not a caller's.
.tf_unfitted <- function(name) {
  force(name)
  list(name = name, backward = function(z, ...) z)
}

# The rowless resumable object the core WOULD return, finalized like any
# other. Validates only what it RECORDS: nkeep is chains x iter_sampling and
# cbind matches on it. `prov` says what ran; `nmiss` is NA because this path
# never scans the matrix. `prov$left_na` exists only under `leave_na`.
.mipca_empty <- function(
  obj,
  subset,
  rank,
  transform,
  prop_diag,
  iter_sampling,
  iter_warmup,
  thin,
  nu,
  chains,
  call,
  leave_na = FALSE
) {
  rank <- .chk_count(rank, "rank")
  chains <- .chk_count(chains, "chains")
  .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  thin <- .chk_count(thin, "thin")
  nu <- .chk_nu(nu)
  ## A column with a hole built this way comes back from as.matrix() still
  ## holed, so that has to be what the caller asked for.
  if (!leave_na && !.nothing_to_impute(obj, subset)) {
    stop("internal error: an empty fit was asked for columns that have holes.")
  }
  if (is.function(transform$fit)) {
    transform <- .tf_unfitted(transform$name)
  }
  call$rank <- rank
  prov <- list(
    chains = chains,
    rank = rank,
    thin = thin,
    nu = nu,
    warmup = 0,
    ndraws = iter_sampling,
    iter = 0,
    nmiss = NA_integer_,
    thin_final = 1L,
    seam_state = "never"
  )
  if (leave_na) {
    prov$left_na <- TRUE
  }
  fit <- .si_mi_pca_new(
    obj = obj,
    map = matrix(integer(0), 0L, 2L, dimnames = list(NULL, c("row", "col"))),
    draws = rep(list(matrix(NA_real_, 0L, iter_sampling)), chains),
    cols = integer(0),
    subset = if (is.null(subset)) seq_len(ncol(obj)) else subset,
    transform = transform,
    seam = NULL,
    prov = prov,
    call = call
  )
  mipca_finalize(fit, prop_diag = prop_diag)
}

# The one implementation behind every door that fits. No holdout here: `obj`
# is always the matrix fitted, and the tuner owns the holdout. `keep_state`
# decides the class returned.
.mipca_run <- function(
  obj,
  rank,
  transform,
  subset = NULL,
  ## The columns that leave the sampler for ABB, in `obj`'s numbering
  ## (.abb_split(), R/mi_pca_abb.R). NULL: no list; every column is sampled.
  abb = NULL,
  ## Forwarded to the core untouched (the cell-level store set; not `subset`).
  store = NULL,
  keep_state = FALSE,
  prop_diag = .prop_diag_default(),
  ## Whether the CALLER named `prop_diag`.
  prop_diag_given = FALSE,
  ## Named, not forwarded blind, so `.state` is unreachable from a door.
  ## The defaults repeat the core's; keep them in sync.
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  refresh = 0L,
  spectra = spectra_opts(),
  init = 1,
  ## Handed down to the core (see mipca_core()).
  seed = NULL,
  seed_given = FALSE,
  ## `...` exists ONLY so `call` matches by full name, never positionally.
  ## It carries nothing: the refusal below empties it.
  ...,
  call = NULL
) {
  .chk_no_dots(..., .what = ".mipca_run()")
  if (keep_state && prop_diag_given) {
    stop(
      "`prop_diag` has no effect with `keep_state = TRUE`. Pass it to ",
      "summary() or to mipca_finalize() instead.",
      call. = FALSE
    )
  }
  .chk_matrix(obj)
  .chk_cells(nrow(obj), ncol(obj))
  ## Before the empty-fit return below, which records the transform's name.
  transform <- .resolve_transform(transform, nrow(obj))
  ## Resolved here: .nothing_to_impute() and .mipca_empty() read it first.
  subset <- .resolve_cols(subset, obj)

  ## THE SPLIT. From here down `obj` is the kept columns and `subset` is in
  ## their numbering, so nothing below knows a list exists; `abb_rec` holds
  ## the requested listed columns until they are drawn.
  abb <- .resolve_cols(abb, obj, what = "abb")
  abb_rec <- NULL
  if (!is.null(abb)) {
    if (!is.null(store)) {
      stop("internal error: `store` and `abb` are not given together.")
    }
    sp <- .abb_split(obj, abb, subset)
    abb_rec <- sp$rec
    ## No requested column is a kept one, so no sampler runs.
    if (sp$n_req_kept == 0L) {
      if (keep_state) {
        stop(
          "Every requested column is listed in `abb`, so no model is fitted ",
          "and there is nothing to resume. Drop `keep_state = TRUE`.",
          call. = FALSE
        )
      }
      ## What .mipca_empty() checks on the other no-sampler path, so a bad
      ## argument is refused whichever columns were asked for; and before
      ## the seed is resolved, so a refused call draws nothing.
      .chk_count(rank, "rank")
      .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
      .chk_count(iter_sampling, "iter_sampling")
      .chk_count(thin, "thin")
      .chk_nu(nu)
      abb_rec$seed <- .resolve_seed(if (seed_given) seed else NULL)
      return(.abb_join(
        NULL,
        abb_rec,
        chains = chains,
        ndraws = iter_sampling,
        thin = thin,
        prop_diag = prop_diag,
        call = if (is.null(call)) {
          .scrub_call(match.call(), ".mipca_run")
        } else {
          call
        }
      ))
    }
    obj <- obj[, sp$kept, drop = FALSE]
    subset <- sp$subset
    ## The holes are all on the listed side. Refused HERE: the core's own
    ## refusal says the matrix has no missing values, and the caller's has.
    if (keep_state && .nothing_to_impute(obj, subset)) {
      stop(
        "No requested column outside `abb` has a missing cell, so no model ",
        "is fitted and there is nothing to resume. Drop `keep_state = TRUE`.",
        call. = FALSE
      )
    }
  }

  ## Only where the answer is a finalized fit: a keep_state request and a
  ## `store` set still go to the core.
  if (!keep_state && is.null(store) && .nothing_to_impute(obj, subset)) {
    cl0 <- if (is.null(call)) .scrub_call(match.call(), ".mipca_run") else call
    empty <- .mipca_empty(
      obj,
      subset = subset,
      rank = rank,
      transform = transform,
      prop_diag = prop_diag,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      thin = thin,
      nu = nu,
      chains = chains,
      call = cl0
    )
    if (is.null(abb_rec)) {
      return(empty)
    }
    ## The kept side had nothing to impute; the listed side still draws, at
    ## fit time, since there is no finalize to wait for. .mipca_empty() made
    ## every refusal above, so the seed is resolved after them.
    abb_rec$seed <- .resolve_seed(if (seed_given) seed else NULL)
    return(.abb_join(
      empty,
      abb_rec,
      chains = chains,
      ndraws = iter_sampling,
      thin = thin,
      prop_diag = prop_diag,
      call = cl0
    ))
  }

  ## Fitted once, on the matrix the sampler sees. The rank ceiling and the
  ## census run first, or the fit fails with R's own messages.
  if (is.function(transform$fit)) {
    .chk_rank(rank, nrow(obj), ncol(obj))
    .check_admissible(obj)
    transform <- transform$fit(obj, rank)
  }

  core <- mipca_core(
    obj,
    forward = transform$forward,
    seed = seed,
    seed_given = seed_given,
    rank = rank,
    subset = subset,
    store = store,
    keep_state = keep_state,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    thin = thin,
    nu = nu,
    chains = chains,
    refresh = refresh,
    spectra = spectra,
    init = init
    ## No `.state`: resume's door is update.si_mi_pca().
  )

  cl <- if (is.null(call)) .scrub_call(match.call(), ".mipca_run") else call
  ## Evaluated, off the core's coerced `rank`: provenance, not a promise.
  cl$rank <- core$rank

  draws <- core$draws
  nkeep <- ncol(draws[[1L]])
  ## Off the core's result, not the seam, so it works without a seam.
  prov <- list(
    chains = core$chains,
    rank = core$rank,
    thin = core$thin,
    nu = core$nu,
    ## Derived rather than re-declared, so it reports what RAN:
    ## iter = iter_warmup + iter_sampling * thin.
    warmup = core$iter - core$thin * nkeep,
    ndraws = nkeep,
    iter = core$iter,
    nmiss = core$nmiss,
    thin_final = 1L,
    seam_state = if (keep_state) "kept" else "never",
    ## The seed the chains started from, as the ABB piece records its own
    ## (cbind.si_mi() says why); NULL when the sampler took the ambient
    ## dqrng state, as a flank window does.
    seed = core$seed
  )

  ## The ABB seed is fixed HERE, at fit time, so when the fit is finalized
  ## does not move its draws: the sampler's own seed, or one draw off R's
  ## stream when the sampler took the ambient dqrng state.
  if (!is.null(abb_rec)) {
    abb_rec$seed <- if (is.null(core$seed)) .resolve_seed(NULL) else core$seed
  }

  fit <- .si_mi_pca_new(
    obj = obj,
    map = core$map,
    draws = draws,
    cols = core$cols,
    subset = core$subset,
    transform = transform,
    seam = core$state,
    prov = prov,
    call = cl,
    abb = abb_rec
  )

  ## Every si_mi that sampled is built by mipca_finalize() (.abb_piece()
  ## builds the one that did not). `prop_diag` is forwarded
  ## because the caller cannot reach mipca_finalize() on this path.
  if (keep_state) fit else mipca_finalize(fit, prop_diag = prop_diag)
}

# The public entry point, a closed signature (no `...`). The defaults repeat
# the core's; keep them in sync.
# For the docs: `store` fails here as "unused argument" (the tuner's store
# set).
#' @export
mipca <- function(
  obj,
  rank,
  transform,
  subset = NULL,
  keep_state = FALSE,
  prop_diag = .prop_diag_default(),
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  refresh = 0L,
  spectra = spectra_opts(),
  init = 1,
  seed = NULL,
  ## Last, so the other arguments keep their positions. For the
  ## docs: abb_cols(obj) proposes a list; the result keeps the input's
  ## column order, as it does with no list.
  abb = NULL
) {
  .mipca_run(
    obj,
    rank = rank,
    transform = transform,
    subset = subset,
    abb = abb,
    keep_state = keep_state,
    prop_diag = prop_diag,
    prop_diag_given = !missing(prop_diag),
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    thin = thin,
    nu = nu,
    chains = chains,
    refresh = refresh,
    spectra = spectra,
    init = init,
    seed = seed,
    seed_given = TRUE,
    call = .scrub_call(match.call(), "mipca")
  )
}

# ---- Scoring: free functions over (draws, truth) -----------------------------
# Independent of mipca so a non-mipca imputer can be scored too. Blanking,
# fitting and scoring live in the tuner (.score_holdout(), R/mi_pca_tune.R);
# these only see numbers.

# Ensemble CRPS per cell (row). Arithmetic in C++ (src/mipca_crps.cpp), where
# the row-wise sort is the entire cost; this is the free-function entry.
#' @export
crps_all <- function(d, truth) {
  crps_block(d, truth)
}

# Randomized PIT on the RANK of the truth among the M draws (Czado, Gneiting
# & Held 2009, eq. 1, applied to the rank). Consumes base RNG unless `u` is
# supplied.
#' @export
pit_all <- function(d, truth, u = NULL) {
  ## Refused, not recycled: `d < truth` recycles a short `truth` silently
  ## (crps_all() refuses the same shape in C++).
  if (length(truth) != nrow(d)) {
    stop("`truth` must supply one value per row of `d`.", call. = FALSE)
  }
  if (is.null(u)) {
    u <- stats::runif(length(truth))
  }
  if (length(u) != length(truth)) {
    stop("`u` must supply one jitter per cell.", call. = FALSE)
  }
  (rowSums(d < truth) + u * (rowSums(d == truth) + 1)) / (ncol(d) + 1)
}
