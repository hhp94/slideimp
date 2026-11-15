# ---- The tuner -----------------------------------------------------------
#
# mipca_tune() sweeps RANK ONLY at nu = Inf and returns a table, one row per
# rank: it scores, it does not store fits, and the tuner and the imputer
# separate completely. Each arm fits the blanked matrix through .mipca_run(),
# is scored here, and is dropped. CRPS and PIT are scored on the data scale;
# the implied nu is read on the LATENT scale, where the I-step applies it.
#
# The BASE RNG ORDER is load-bearing: sampler seed -> na_loc draw -> diag draw
# -> PIT jitter -> per-arm band. dqrng is cloned, never advanced; a NULL seed
# is drawn from R's stream, so set.seed() governs the sampler.

.chk_ranks <- function(rank) {
  ## The integer ceiling is part of the predicate: as.integer() would turn a
  ## larger rank into an NA that passes the duplicate check.
  if (
    !is.numeric(rank) ||
      length(rank) < 1L ||
      anyNA(rank) ||
      any(!is.finite(rank)) ||
      any(rank < 1) ||
      any(rank != trunc(rank)) ||
      any(rank > .Machine$integer.max)
  ) {
    stop("`rank` must be a vector of whole numbers >= 1, no NA.", call. = FALSE)
  }
  if (anyDuplicated(rank)) {
    ## Two runs at one rank reproduce bitwise, so a repeat buys nothing.
    stop("`rank` has duplicates. Give each rank once.", call. = FALSE)
  }
  as.integer(rank)
}

# The base stream's position. NULL if nothing has been drawn yet.
.base_state <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    return(NULL)
  }
  get(".Random.seed", envir = globalenv(), inherits = FALSE)
}

# Put the base stream back. NULL means remove the seed rather than set one,
# so an untouched stream stays untouched.
.base_restore <- function(s) {
  if (is.null(s)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
    return(invisible(NULL))
  }
  assign(".Random.seed", s, envir = globalenv())
  invisible(NULL)
}

# ---- The holdout: blanking, positions, scores --------------------------------
# The tuner is the holdout's only user; the fit never learns a cell was held
# out.

# The holdout's whole leakage surface: blank by NA, never pass positions on.
.blank_cells <- function(X, idx) {
  X[idx] <- NA_real_
  X
}

# A supplied (row, col) design -> linear indices in map order; no RNG.
# Messages say "the holdout": mipca_tune() prefixes its own argument name.
.resolve_holdout <- function(obj, holdout) {
  if (
    !is.matrix(holdout) ||
      !is.numeric(holdout) ||
      ncol(holdout) != 2L ||
      nrow(holdout) < 1L
  ) {
    stop(
      "the holdout must be a numeric two-column (row, col) matrix ",
      "with at least one row.",
      call. = FALSE
    )
  }
  if (anyNA(holdout)) {
    stop("the holdout contains NA positions.", call. = FALSE)
  }
  if (any(holdout != trunc(holdout))) {
    stop("holdout positions must be whole numbers.", call. = FALSE)
  }
  hr <- as.integer(holdout[, 1L])
  hc <- as.integer(holdout[, 2L])
  n <- nrow(obj)
  if (any(hr < 1L | hr > n) || any(hc < 1L | hc > ncol(obj))) {
    stop("holdout positions fall outside `obj`.", call. = FALSE)
  }
  idx <- hr + (hc - 1L) * n
  if (anyDuplicated(idx)) {
    stop("the holdout contains duplicate positions.", call. = FALSE)
  }
  if (!all(is.finite(obj[idx]))) {
    stop(
      "holdout positions must be observed, finite cells of `obj`.",
      call. = FALSE
    )
  }
  ## Threshold 1, not the core's 2: .check_admissible() owns admissibility and
  ## re-checks after blanking; this exists only to name the holdout. is.na()
  ## must stay .mat_miss_int()'s predicate (NaN missing, Inf not).
  hit <- tabulate(hc, nbins = ncol(obj))
  hitc <- which(hit > 0L)
  nobs <- n - colSums(is.na(obj[, hitc, drop = FALSE]))
  if (any(nobs - hit[hitc] < 1L)) {
    bad <- .col_label(obj, hitc[nobs - hit[hitc] < 1L])
    stop(
      "the holdout covers every observed cell of ",
      length(bad),
      " column",
      .s(length(bad)),
      " [",
      fmt_trunc(bad, 8),
      "]. Hold out fewer cells ",
      "in ",
      if (length(bad) == 1L) "it" else "them",
      ".",
      call. = FALSE
    )
  }
  sort(idx)
}

# Per-cell scores on data-scale draws, through the free functions. `u` is
# always the caller's: one jitter shared by every arm.
.score_holdout <- function(d, truth, u) {
  pit <- pit_all(d, truth, u = u)
  mu <- rowMeans(d)
  list(
    pit = pit,
    crps = crps_all(d, truth),
    mean = mu,
    var = rowSums((d - mu)^2) / (ncol(d) - 1L)
  )
}

# Holdout positions in the blanked matrix's hole enumeration, by rank
# arithmetic over the disjoint merge, then guarded.
.tune_hole_pos <- function(na0, ho_lin) {
  u <- c(na0, ho_lin)
  o <- order(u)
  holes <- u[o]
  rk <- integer(length(u))
  rk[o] <- seq_along(u)
  ho_pos <- rk[length(na0) + seq_along(ho_lin)]
  if (!identical(holes[ho_pos], ho_lin)) {
    stop("internal error: holdout positions do not index the hole enumeration.")
  }
  list(holes = holes, ho_pos = ho_pos)
}

# The diagnostic sample over ALL holes, not just the holdout's, stratified by
# column missingness so it spans the gradient where slow mixing concentrates.
# Base RNG only. Returns the drawn positions in the hole enumeration,
# ascending.
.tune_diag_draw <- function(holes, n, p, n_diag) {
  colof <- .lin_col(holes, n)
  cnt <- tabulate(colof, nbins = p)
  dirty <- which(cnt > 0L)
  ns <- min(10L, length(dirty))
  ord <- order(cnt[dirty], dirty)
  strat_of_dirty <- integer(length(dirty))
  strat_of_dirty[ord] <- ceiling(ns * seq_along(ord) / length(ord))
  strat_col <- integer(p)
  strat_col[dirty] <- strat_of_dirty
  pools <- split(seq_along(holes), strat_col[colof])
  if (length(pools) != ns) {
    stop("internal error: diag strata do not tile the dirty columns.")
  }
  alloc <- rep(n_diag %/% ns, ns)
  extra <- n_diag %% ns
  if (extra > 0L) {
    top <- ns - seq_len(extra) + 1L
    alloc[top] <- alloc[top] + 1L
  }
  picks <- vector("list", ns)
  for (s in seq_len(ns)) {
    pool <- pools[[s]]
    k <- min(alloc[s], length(pool))
    picks[[s]] <- pool[sample.int(length(pool), k)]
  }
  pos <- sort(unlist(picks, use.names = FALSE))
  if (length(pos) < 1L) {
    stop("internal error: the diagnostic draw is empty.")
  }
  pos
}

# Sorted union of holdout and diag positions; an overlap is stored once and
# read twice. Rows by arithmetic, guarded.
.tune_store <- function(ho_pos, diag_pos) {
  s <- c(ho_pos, diag_pos)
  o <- order(s)
  ss <- s[o]
  first <- c(TRUE, diff(ss) > 0L)
  store <- ss[first]
  slot <- integer(length(s))
  slot[o] <- cumsum(first)
  rows_ho <- slot[seq_along(ho_pos)]
  rows_diag <- slot[length(ho_pos) + seq_along(diag_pos)]
  if (
    !identical(store[rows_ho], ho_pos) ||
      !identical(store[rows_diag], diag_pos)
  ) {
    stop("internal error: store rows do not round-trip their positions.")
  }
  list(store = store, rows_ho = rows_ho, rows_diag = rows_diag)
}

# Robust/sd for a pure t at nu, monotone up in nu, tending to 0.9999848 and
# not 1, because 1.349 is a rounded 2 * qnorm(0.75). NA in, NA out.
.implied_ratio <- function(nu) {
  (2 * stats::qt(0.75, nu) / 1.349) / sqrt(nu / (nu - 2))
}

.robust_width <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) < 4L) {
    return(NA_real_)
  }
  stats::IQR(z) / 1.349
}

.implied_nu <- function(r) {
  if (is.na(r)) {
    return(NA_real_)
  }
  ## Must stay .implied_ratio(nu_max), the bracket's own endpoint: a literal 1
  ## or the analytic limit crashes uniroot.
  nu_max <- 1e7
  if (!is.finite(r) || r >= .implied_ratio(nu_max)) {
    return(Inf)
  }
  if (r <= .implied_ratio(2.001)) {
    return(2.001)
  }
  stats::uniroot(
    function(nu) .implied_ratio(nu) - r,
    c(2.001, nu_max),
    tol = 1e-8
  )$root
}

# The one derivation, z -> robust width -> implied nu, for the arm read, the
# bootstrap and mipca_fit()'s re-read; each re-reads the cells, never the grid.
# `z` is latent-scale: (truth_lat - mean_lat) / sd_lat off the result table.
.tune_nu_of <- function(z) {
  r <- .robust_width(z)
  list(robust_z = r, nu = .implied_nu(r))
}

# Cluster bootstrap over whole holdout columns, the unit the holdout is drawn
# in. Consumes base RNG.
.tune_nu_band <- function(z, hcol, n_boot = 500L) {
  zs <- split(z, hcol)
  b <- vapply(
    seq_len(n_boot),
    function(k) {
      pick <- sample.int(length(zs), length(zs), replace = TRUE)
      .tune_nu_of(unlist(zs[pick], use.names = FALSE))$nu
    },
    0
  )
  ## All-NA when the width was unreadable (< 4 finite z); quantile() on
  ## nothing errors, so return NA rather than crash.
  if (all(is.na(b))) {
    return(c(NA_real_, NA_real_))
  }
  stats::quantile(b, c(0.16, 0.84), na.rm = TRUE, names = FALSE)
}

# min() of an all-NA vector under na.rm is Inf with a warning, not NA. The
# median/mean reductions return NA on their own and stay unguarded.
.na_red <- function(x, f, ...) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  f(x, ..., na.rm = TRUE)
}

# THE TABLE ROW: the only copy of the schema, so error rows and scored rows
# cannot drift apart under rbind(). Every summary column is a reduction of
# its own row's nested tables.
.tune_row <- function(
  rank,
  took,
  error = NA_character_,
  crps = NA_real_,
  implied_nu = NA_real_,
  implied_lo = NA_real_,
  implied_hi = NA_real_,
  robust_z = NA_real_,
  ess_min = NA_real_,
  ess_p10 = NA_real_,
  ess_med = NA_real_,
  ess_ok = NA_real_,
  rhat_med = NA_real_,
  rhat_p90 = NA_real_,
  rhat_max = NA_real_,
  rhat_ok = NA_real_,
  ndraws = NA_integer_,
  n_diag = NA_integer_,
  sigma = list(NULL),
  result = list(NULL),
  diag = list(NULL)
) {
  data.frame(
    rank = rank,
    crps = crps,
    implied_nu = implied_nu,
    implied_lo = implied_lo,
    implied_hi = implied_hi,
    robust_z = robust_z,
    ess_min = ess_min,
    ess_p10 = ess_p10,
    ess_med = ess_med,
    ess_ok = ess_ok,
    rhat_med = rhat_med,
    rhat_p90 = rhat_p90,
    rhat_max = rhat_max,
    rhat_ok = rhat_ok,
    ndraws = ndraws,
    n_diag = n_diag,
    sigma = I(sigma),
    time = took,
    error = error,
    result = I(result),
    diag = I(diag),
    stringsAsFactors = FALSE
  )
}

# ONE arm, read while its fit is alive; nothing returned holds a draw. `rec`
# is the run's record, built once before the arms. Base RNG: the band only.
.tune_arm_read <- function(fit, rec, rank, took) {
  ## LEAKAGE-CLASS guard: the store's holdout rows must address exactly the
  ## recorded cells, or every score below is against the wrong truth.
  mh <- fit$map[rec$rows_ho, , drop = FALSE]
  md <- fit$map[rec$rows_diag, , drop = FALSE]
  if (!identical(.map_lin(mh, rec$n), rec$lin)) {
    stop(
      "internal error: the arm's holdout rows do not address the recorded cells."
    )
  }

  ## Against what was ASKED FOR, never a neighbour arm. On value: the core
  ## hands warmup back double and thin back integer, the reverse of the door.
  same <- function(a, b) {
    if (is.numeric(a) && is.numeric(b)) {
      identical(as.double(a), as.double(b))
    } else {
      identical(a, b)
    }
  }
  got <- list(
    chains = fit$prov$chains,
    iter_warmup = fit$prov$warmup,
    iter_sampling = fit$prov$ndraws,
    thin = fit$prov$thin,
    spectra = fit$seam$spectra
  )
  bad <- names(got)[!mapply(same, got, rec$run[names(got)])]
  if (length(bad) > 0L) {
    stop(
      "internal error: an arm did not run the way the tuner asked (",
      fmt_trunc(bad, 8),
      ")."
    )
  }

  D <- .pool_rows(fit$draws, rec$rows_ho)
  d <- .backward_at(fit$transform, D, fit$map, rec$rows_ho)
  sc <- .score_holdout(d, rec$truth, u = rec$u)
  sdv <- sqrt(sc$var)

  ## nu on the LATENT scale, where the I-step applies it: the truth through
  ## THIS arm's forward map (each arm fits its own), on the holdout columns'
  ## full blocks the record carries, read at the cells by the record's
  ## positions. The positions are arithmetic, so they are guarded: the latent
  ## truth must map back onto the recorded truth, or the read is off another
  ## cell (LEAKAGE-CLASS: a wrong cell reads as a plausible z).
  tl <- .forward_at(fit$transform, rec$blk, rec$ucol)[rec$pos]
  rt <- .backward_at(fit$transform, tl, fit$map, rec$rows_ho)
  if (
    anyNA(tl) ||
      max(abs(rt - rec$truth)) > 1e-8 * max(1, max(abs(rec$truth)))
  ) {
    stop(
      "internal error: the arm's forward map does not return the recorded truth at the holdout cells."
    )
  }
  mu_l <- rowMeans(D)
  sd_l <- sqrt(rowSums((D - mu_l)^2) / (ncol(D) - 1L))
  z <- (tl - mu_l) / sd_l
  nuz <- .tune_nu_of(z)
  band <- .tune_nu_band(z, rec$hcol)

  dg <- .cell_diag_read(fit$draws, rec$rows_diag)
  rhat <- dg[, 1L]
  ess <- dg[, 2L]
  .tune_row(
    rank = rank,
    took = took,
    crps = mean(sc$crps),
    implied_nu = nuz$nu,
    implied_lo = band[[1L]],
    implied_hi = band[[2L]],
    robust_z = nuz$robust_z,
    ## All-NA when the run is too short to estimate ESS; NA, not Inf.
    ess_min = .na_red(ess, min),
    ess_p10 = .na_red(ess, stats::quantile, 0.10, names = FALSE),
    ess_med = stats::median(ess, na.rm = TRUE),
    ess_ok = mean(ess >= rec$ess_bar, na.rm = TRUE),
    rhat_med = stats::median(rhat, na.rm = TRUE),
    rhat_p90 = stats::quantile(rhat, 0.90, na.rm = TRUE, names = FALSE),
    rhat_max = .na_red(rhat, max),
    rhat_ok = mean(rhat <= rec$rhat_bar, na.rm = TRUE),
    ndraws = length(fit$draws) * ncol(fit$draws[[1L]]),
    n_diag = length(rec$rows_diag),
    sigma = list(vapply(fit$seam$chain, `[[`, 0, "sigma")),
    result = list(data.frame(
      col = mh[, "col"],
      row = mh[, "row"],
      truth = rec$truth,
      mean = sc$mean,
      sd = sdv,
      crps = sc$crps,
      pit = sc$pit,
      ## The latent triple mipca_fit() re-reads nu from; the fit cannot
      ## rebuild it, its transform being re-fit on the full matrix.
      truth_lat = tl,
      mean_lat = mu_l,
      sd_lat = sd_l,
      stringsAsFactors = FALSE
    )),
    diag = list(data.frame(
      col = md[, "col"],
      row = md[, "row"],
      rhat = rhat,
      ess = ess,
      stringsAsFactors = FALSE
    ))
  )
}

# The matrix's fingerprint on the table, so the handoff can refuse another
# matrix. It cannot catch a margin-preserving edit, or a value edit under the
# tolerance away from the scored cells. Compiled scans only: no n x p mask.
.tune_fp <- function(X) {
  row_na <- mat_miss(X, col = FALSE, prop = FALSE)
  if (!is.integer(row_na) || length(row_na) != nrow(X) || anyNA(row_na)) {
    stop(
      "internal error: mat_miss() did not return one integer count per row."
    )
  }
  list(
    dim = dim(X),
    dimnames = dimnames(X),
    row_na = unname(row_na),
    col_na = .mat_miss_int(X),
    row_sum = unname(rowSums(X, na.rm = TRUE)),
    col_sum = unname(colSums(X, na.rm = TRUE)),
    ## Not the sums themselves: a near-zero margin would refuse on rounding.
    scale = max(abs(range(X, na.rm = TRUE)))
  )
}

# NA counts exactly; sums on .SEAM_FP_TOL, never bitwise, since a saved table
# may sum in another order on another build.
.tune_fp_check <- function(fp, X) {
  what <- NULL
  if (!identical(dim(X), fp$dim) || !identical(dimnames(X), fp$dimnames)) {
    what <- "its dimensions or dimnames"
  } else {
    cur <- .tune_fp(X)
    if (
      !identical(cur$row_na, fp$row_na) || !identical(cur$col_na, fp$col_na)
    ) {
      what <- "its missing cells"
    } else {
      den <- max(fp$scale, .Machine$double.xmin)
      rel <- c(
        abs(cur$row_sum - fp$row_sum) /
          (pmax(fp$dim[2L] - fp$row_na, 1L) * den),
        abs(cur$col_sum - fp$col_sum) / (pmax(fp$dim[1L] - fp$col_na, 1L) * den)
      )
      if (!all(is.finite(rel)) || max(rel) >= .SEAM_FP_TOL) {
        what <- "its values"
      }
    }
  }
  if (!is.null(what)) {
    stop(
      "`obj` is not the matrix this table was tuned on: ",
      what,
      " differ. ",
      "Pass that matrix, or run mipca_tune() on this one.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @export
mipca_tune <- function(
  obj,
  rank,
  na_loc = NULL,
  ## (columns, rows per column), not a cell total: the band bootstraps
  ## COLUMNS. NULL means the default span below.
  n_cols = NULL,
  n_rows = 4L,
  transform,
  ## A count, not a fraction: storage is what is being rationed here.
  n_diag = 10000L,
  ## The tuner's own run length, not the core's, so a change to the imputer's
  ## defaults cannot move the length the arms are ranked at.
  iter_warmup = 100L,
  iter_sampling = 100L,
  chains = 4L,
  ## Named, not `...`, so a typo is refused at the door. `forward` and
  ## `.state` are deliberately absent.
  thin = 1,
  init = 1,
  refresh = 0L,
  spectra = spectra_opts(),
  ## NULL draws the sampler's seed from R's stream, so set.seed() governs it.
  seed = NULL,
  ## The columns that leave the sampler for ABB, as at mipca(). Last, so it
  ## shifts no other argument's position.
  abb = NULL
) {
  ranks <- .chk_ranks(rank)
  .chk_matrix(obj)
  ## THE LIST. From here down `obj` is the kept columns, so the holdout is
  ## drawn from them and every index below is in their numbering. `full` is
  ## the caller's matrix, for the fingerprint and for a supplied `na_loc`,
  ## which stays in the caller's numbering; `kept` maps a kept column back
  ## to it, for every column this function reports (NULL: no list).
  abb <- .resolve_cols(abb, obj, what = "abb")
  full <- obj
  kept <- NULL
  if (!is.null(abb)) {
    kept <- .abb_kept(ncol(obj), abb)$kept
    obj <- obj[, kept, drop = FALSE]
  }
  ## A design the tuner draws itself is in the kept numbering and is
  ## reported back through `kept`; a supplied one never leaves the caller's.
  drawn <- is.null(na_loc)
  n <- nrow(obj)
  p <- ncol(obj)
  ## The hole enumeration is integer arithmetic on linear indices; so is a
  ## supplied design's, which is resolved against the caller's matrix.
  if (as.double(n) * ncol(full) > .Machine$integer.max) {
    stop(
      "`obj` has ",
      format(as.double(n) * ncol(full), big.mark = ",", scientific = FALSE),
      " cells, more than mipca_tune() can index (",
      .Machine$integer.max,
      "). Tune on part of its columns.",
      call. = FALSE
    )
  }
  ## Resolved once: the table, every arm and the handoff read one transform.
  transform <- .resolve_transform(transform, n)
  n_diag <- .chk_count(n_diag, "n_diag")
  chains <- .chk_count(chains, "chains")
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  ## Not .chk_count: zero warmup is legal at the core, so the floor is 0.
  iter_warmup <- .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  ## A double, the core's own type for thin.
  thin <- as.double(.chk_count(thin, "thin"))
  ## At the door: inside the per-arm tryCatch a bad value costs a row per arm.
  init <- .chk_init(init)
  refresh <- .chk_refresh(refresh)
  spectra <- .chk_spectra(spectra)
  ## At the door too: pre-empts sample_na_loc()'s cli refusals with the plain
  ## one before a generator is touched. The arms repeat it on the blanking.
  cen <- .check_admissible(obj)
  ## Applied below, after the holdout refusals.
  if (!is.null(seed)) {
    seed <- .chk_count(seed, "`seed`", lo = 0L)
  }

  ## ---- the holdout draw, or the caller's ------------------------------------
  ## The one discovery scan; every later index is arithmetic on it.
  na0 <- which(is.na(obj))
  ## Before the branch: the supplied branch reads `n_cols` too, and a bad TYPE
  ## must be refused as such, not as a disagreement.
  if (!is.null(n_cols)) {
    n_cols <- .chk_count(n_cols, "n_cols")
  }
  if (is.null(na_loc)) {
    ## Drawn in COLUMNS: n_rows cells from each of n_cols columns.
    if (is.null(n_cols)) {
      n_cols <- max(50L, min(500L, as.integer(floor(0.10 * p))))
    }
    n_rows <- .chk_count(n_rows, "n_rows")
    ## Ahead of sample_na_loc(), so its cli_abort refusals are unreachable.
    ## The helper never draws a constant column (variance below double.eps,
    ## the census's `const`), so the count leaves those out to match it.
    ## Keep two observed cells, the census's floor; no missingness cap: the
    ## sampler refuses only what it cannot compute.
    nobs <- n - cen$cmiss
    eligible <- sum(nobs >= n_rows + 2L & !cen$const)
    if (eligible < n_cols) {
      stop(
        "Cannot place a holdout over ",
        n_cols,
        " columns at ",
        n_rows,
        " cells each. That needs ",
        n_cols,
        " non-constant columns with at ",
        "least ",
        n_rows + 2L,
        " observed cells, and `obj` has ",
        eligible,
        ". Lower `n_cols` or `n_rows`, or supply `na_loc`.",
        call. = FALSE
      )
    }
    ## Known before the draw, so the cell floor below refuses first.
    n_ho_cells <- n_cols * n_rows
  } else {
    ## missing(), since the default is a number.
    if (!missing(n_rows)) {
      stop(
        "`n_rows` has no effect when `na_loc` is given, because `na_loc` ",
        "places every held-out cell. Drop one of them.",
        call. = FALSE
      )
    }
    if (!is.list(na_loc) || length(na_loc) < 1L) {
      stop(
        "`na_loc` must be a list of one two-column (row, col) matrix, as ",
        "sample_na_loc() returns with `n_reps = 1`.",
        call. = FALSE
      )
    }
    if (length(na_loc) > 1L) {
      stop(
        "`na_loc` holds ",
        length(na_loc),
        " designs, and mipca_tune() scores ",
        "one. Call mipca_tune() once per design. The tables do not combine ",
        "into one grid, because each design scores different cells.",
        call. = FALSE
      )
    }
    n_ho_cells <- NROW(na_loc[[1L]])
  }

  ## A floor against degenerate grids, not a sufficiency claim.
  if (n_ho_cells <= 2L) {
    stop(
      "The holdout has ",
      n_ho_cells,
      " cell",
      .s(n_ho_cells),
      ", and ",
      "mipca_tune() needs more than two to score a rank. Raise `n_cols` or ",
      "`n_rows`, or supply an `na_loc` with more rows.",
      call. = FALSE
    )
  }

  ## Prefixed: .resolve_holdout() says "the holdout"; this argument is na_loc.
  ## `against` is the matrix whose numbering the design is in.
  resolve <- function(design, against) {
    tryCatch(
      .resolve_holdout(against, design),
      error = function(e) {
        stop("`na_loc`: ", conditionMessage(e), call. = FALSE)
      }
    )
  }
  Xb <- NULL
  if (!is.null(na_loc)) {
    ## Resolved against the CALLER's matrix, whose numbering it is in, so
    ## .resolve_holdout() owns every shape refusal with or without a list;
    ## then re-based onto the kept columns, which refuses a listed one.
    ho_lin <- resolve(na_loc[[1L]], full)
    if (!is.null(kept)) {
      ho_lin <- .abb_tune_rebase(ho_lin, full, kept)
    }
    ## After resolve(), which established the two-column shape.
    if (
      !is.null(n_cols) &&
        length(unique(na_loc[[1L]][, 2L])) != n_cols
    ) {
      stop(
        "`n_cols` disagrees with the supplied `na_loc`. Drop one of them.",
        call. = FALSE
      )
    }
    ## The core's census on the blanking, once and named against na_loc,
    ## rather than one error row per arm. Xb is the matrix every arm fits.
    ## `obj` passed the census above, so only the holdout can fail it here.
    Xb <- .blank_cells(obj, ho_lin)
    .check_admissible(
      Xb,
      lead = paste0(
        "With the `na_loc` cells held out, some columns of `obj` cannot be ",
        "imputed."
      ),
      fix = "Hold out fewer cells in these columns."
    )
  }

  ## SEEDED HERE, after every refusal and before the first draw, so a refused
  ## call touches neither generator. RNG order: this file's header.
  if (!is.null(seed)) {
    base_old <- .base_state()
    dq_old <- .seed_dq(seed)
    on.exit(
      {
        dqrng::dqrng_set_state(dq_old)
        .base_restore(base_old)
      },
      add = TRUE
    )
    set.seed(seed)
  } else {
    ## The sampler's seed off the caller's base stream; dqrng put back on
    ## exit, base left advanced by that one draw.
    dq_old <- .seed_dq(NULL)
    on.exit(dqrng::dqrng_set_state(dq_old), add = TRUE)
  }
  if (is.null(na_loc)) {
    ## rowmax = colmax = 1 switch off the helper's 90% caps, whose cli
    ## refusal covers rows and columns the sampler takes.
    na_loc <- sample_na_loc(
      obj,
      n_cols = n_cols,
      n_rows = n_rows,
      n_reps = 1L,
      rowmax = 1,
      colmax = 1
    )
    ## Guards the helper's contract, not the caller.
    if (NROW(na_loc[[1L]]) != n_ho_cells) {
      stop(
        "internal error: the holdout draw placed ",
        NROW(na_loc[[1L]]),
        " cells where ",
        n_ho_cells,
        " were asked for."
      )
    }
    ho_lin <- resolve(na_loc[[1L]], obj)
  }
  hp <- .tune_hole_pos(na0, ho_lin)
  diag_pos <- .tune_diag_draw(hp$holes, n, p, n_diag)
  st <- .tune_store(hp$ho_pos, diag_pos)
  ## The holdout columns' FULL blocks, truth present, for the latent-scale nu
  ## read: forward(blk, cols) maps column blocks, and each arm maps them
  ## through its own fitted transform. unique() is one-time discovery of the
  ## column set; `pos` is the scatter table into the forward block, computed
  ## rather than searched. n x n_ho_cols doubles, held for the run.
  hcol <- .lin_col(ho_lin, n)
  ucol <- unique(hcol)
  inv <- integer(p)
  inv[ucol] <- seq_along(ucol)
  ## THE RUN'S RECORD: every arm is read against this one.
  rec <- list(
    n = n,
    lin = ho_lin,
    truth = obj[ho_lin],
    hcol = hcol,
    ucol = ucol,
    blk = obj[, ucol, drop = FALSE],
    pos = cbind(ho_lin - (hcol - 1L) * n, inv[hcol]),
    rows_ho = st$rows_ho,
    rows_diag = st$rows_diag,
    ## ONE jitter, shared across every arm: PIT deciles paired across arms.
    u = stats::runif(length(ho_lin)),
    ess_bar = 100L * chains,
    rhat_bar = 1.05,
    ## The six settings mipca_fit() carries, in its order; must match its
    ## `given` list. The tuner's `refresh` is not one.
    run = list(
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      thin = thin,
      spectra = spectra,
      init = init
    )
  )
  ## Once; a supplied design's was built at the door.
  if (is.null(Xb)) {
    Xb <- .blank_cells(obj, ho_lin)
  }

  rows <- vector("list", length(ranks))
  for (i in seq_along(ranks)) {
    t0 <- proc.time()[["elapsed"]]
    ## keep_state = TRUE: the Rhat/ESS read needs the chains axis.
    fit <- tryCatch(
      .mipca_run(
        Xb,
        rank = ranks[i],
        transform = transform,
        store = st$store,
        keep_state = TRUE,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        chains = chains,
        thin = thin,
        init = init,
        refresh = refresh,
        spectra = spectra
      ),
      error = identity
    )
    took <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      ## Caught, not fatal: one bad row must not cost the whole grid.
      rows[[i]] <- .tune_row(ranks[i], took, error = conditionMessage(fit))
      next
    }
    rows[[i]] <- .tune_arm_read(fit, rec, ranks[i], took)
    ## Dropped before the next arm runs, so the tuner's peak is one arm.
    rm(fit)
  }

  out <- do.call(rbind, rows)
  ## REPORTED IN THE CALLER'S NUMBERING: every column the table names -- the
  ## scored cells, the diagnostic cells, a design drawn here -- goes back
  ## through `kept`, each index guarded to be a kept position first. An arm
  ## that errored has no tables. A supplied design was never re-based.
  if (!is.null(kept)) {
    back <- function(k) {
      if (
        length(k) > 0L &&
          (anyNA(k) || min(k) < 1L || max(k) > length(kept))
      ) {
        stop("internal error: a reported column is not a kept position.")
      }
      kept[k]
    }
    for (i in seq_len(nrow(out))) {
      if (!is.null(out$result[[i]])) {
        out$result[[i]]$col <- back(out$result[[i]]$col)
        out$diag[[i]]$col <- back(out$diag[[i]]$col)
      }
    }
    if (drawn) {
      na_loc[[1L]][, 2L] <- back(na_loc[[1L]][, 2L])
    }
  }
  ## One attribute: row subsets and rbind() keep it (the FIRST table's), a
  ## column subset drops it. `rank` lets .tune_record() spot a combined table.
  attr(out, "tune") <- list(
    transform = transform,
    run = rec$run,
    ## The caller's matrix, listed columns included: the handoff is given
    ## that one, and fits it with the same list.
    fp = .tune_fp(full),
    na_loc = na_loc,
    rank = ranks,
    n_diag = n_diag,
    ess_bar = rec$ess_bar,
    rhat_bar = rec$rhat_bar,
    ## Carried to mipca_fit(), as the transform is. NULL: no list.
    abb = abb
  )
  ## `slideimp_tbl` second, as tune_imp() carries it.
  class(out) <- c("si_tune", "slideimp_tbl", "data.frame")
  out
}

# The record mipca_fit() hands off from; it does not carry the matrix, which
# is passed again. The %in% is caller validation, not a lookup.
.tune_record <- function(object) {
  tn <- attr(object, "tune")
  if (
    !is.list(tn) ||
      is.null(tn$fp) ||
      is.null(tn$run) ||
      is.null(tn$transform) ||
      is.null(tn$rank)
  ) {
    ## Hint for the docs: a row subset keeps the record, while a column
    ## subset or subset() drops it.
    stop(
      "This table has lost the settings mipca_fit() reads. Pass the table ",
      "mipca_tune() returned, or a row subset of it.",
      call. = FALSE
    )
  }
  r <- object$rank
  if (!is.numeric(r) || anyDuplicated(r) || !all(r %in% tn$rank)) {
    ## rbind() keeps the first table's record only, so rows from another run
    ## would be fitted under it.
    stop(
      "This table holds rows from more than one mipca_tune() call. Pass the ",
      "table from one call.",
      call. = FALSE
    )
  }
  tn
}

# A bare handoff's rank: argmin mean CRPS, ties to the smaller rank. A default,
# not a chooser: rank is the caller's to choose and the tuner's to score.
.tune_argmin_rank <- function(object) {
  live <- is.na(object$error) & is.finite(object$crps)
  if (!any(live)) {
    stop(
      "No rank in this table has holdout scores, so mipca_fit() cannot ",
      "choose one. The `error` column gives the cause for each rank. Fix it, ",
      "then run mipca_tune() again.",
      call. = FALSE
    )
  }
  r <- object$rank[live]
  r[order(object$crps[live], r)[1L]]
}

# The tuner -> imputer handoff, mipca_fit(tuned, obj): the matrix is passed
# again. Prefixed: a bare `fit` would mask generics::fit.
#' @export
mipca_fit <- function(object, ...) {
  UseMethod("mipca_fit")
}

#' @export
mipca_fit.si_tune <- function(
  object,
  ## The matrix tuned on, passed again and checked against the table's
  ## fingerprint. Second, so the journey reads mipca_fit(tuned, X).
  obj,
  rank = NULL,
  ## `...` HERE makes every setting below full-name only; it carries nothing.
  ## A method cannot drop it: the generic has `...`.
  ...,
  ## The six carried settings. No defaults, on purpose: an unnamed one is
  ## CARRIED from the arms, read by missing() below.
  chains,
  iter_warmup,
  iter_sampling,
  thin,
  spectra,
  init,
  ## Not carried; CmdStan's own default.
  refresh = 100L,
  ## Not carried: straight to .mipca_run() at its own defaults.
  keep_state = FALSE,
  prop_diag = .prop_diag_default(),
  ## This fit's seed, as at mipca(); not carried from the tuner.
  seed = NULL,
  ## NULL reads nu off the chosen arm and clamps it up to .NU_FLOOR.
  nu = NULL,
  ## Named only so the refusals below can say why.
  transform,
  holdout,
  abb
) {
  .chk_no_dots(..., .what = "mipca_fit()")
  owned <- c(
    transform = !missing(transform),
    holdout = !missing(holdout),
    abb = !missing(abb)
  )
  if (any(owned)) {
    hit <- names(owned)[owned]
    ## For the docs: the transform and `abb` are carried because rank and nu
    ## were tuned under them, and a holdout is the tuner's alone, because a
    ## production fit must not be conditioned on deliberately degraded data.
    stop(
      "mipca_fit() does not take ",
      paste(sprintf("`%s`", hit), collapse = ", "),
      ". It fits every observed cell, with the transform and `abb` stored in ",
      "the table. To change them, run mipca_tune() again, or call mipca() ",
      "directly.",
      call. = FALSE
    )
  }
  if (missing(obj)) {
    stop(
      "mipca_fit() needs the matrix the table was tuned on, passed again: ",
      "mipca_fit(tuned, obj). The table does not carry it.",
      call. = FALSE
    )
  }
  .chk_matrix(obj)
  tn <- .tune_record(object)
  .tune_fp_check(tn$fp, obj)

  if (is.null(rank)) {
    rank_src <- "argmin CRPS"
    rank <- .tune_argmin_rank(object)
  } else {
    rank_src <- "given"
    rank <- .chk_count(rank, "rank")
  }
  ## By key over a handful of unique ranks; caller validation, not a lookup.
  k <- which(object$rank == rank)
  if (length(k) != 1L) {
    stop(
      "Rank ",
      rank,
      " was not tuned in this table, which has rank",
      .s(length(object$rank)),
      " ",
      paste(sort(object$rank), collapse = ", "),
      ". mipca_fit() reads `nu` ",
      "from a tuned rank. Tune this rank with mipca_tune() on the same ",
      "`na_loc`, or call mipca() with `nu` set.",
      call. = FALSE
    )
  }
  res <- object$result[[k]]
  if (is.null(res) || nrow(res) < 1L) {
    err <- object$error[k]
    stop(
      "Rank ",
      rank,
      " has no holdout scores, so mipca_fit() has no `nu` to ",
      "read. Choose another rank, or call mipca() with `nu` set.",
      if (!is.na(err)) paste0("\nIts `error` entry: ", err) else "",
      call. = FALSE
    )
  }
  ## Exact: values, not arithmetic, so there is no build to forgive.
  if (!identical(obj[cbind(res$row, res$col)], res$truth)) {
    stop(
      "`obj` is not the matrix this table was tuned on: its values at the ",
      "holdout cells differ. Pass that matrix, or run mipca_tune() on this ",
      "one.",
      call. = FALSE
    )
  }

  ## nu re-read off the chosen rank's cells, not looked up (see .tune_nu_of()),
  ## on the latent scale the table carries for it.
  if (is.null(res$truth_lat)) {
    stop(
      "This table lacks the cells that mipca_fit() reads `nu` from. Run ",
      "mipca_tune() again to make a new table.",
      call. = FALSE
    )
  }
  nuz <- .tune_nu_of((res$truth_lat - res$mean_lat) / res$sd_lat)
  rz <- nuz$robust_z
  nu_implied <- nuz$nu
  if (is.null(nu)) {
    nu_src <- "tuner"
    if (is.na(nu_implied)) {
      ## The reading needs 4 finite standardized errors at the holdout cells.
      stop(
        "mipca_fit() cannot read `nu` at rank ",
        rank,
        ": fewer than 4 ",
        "holdout cells give a usable score. Pass `nu`.",
        call. = FALSE
      )
    }
    ## A message, not a warning: both numbers are also kept in the fit's
    ## `prov$tuned` record, so the clamp is never silent.
    nu <- max(nu_implied, .NU_FLOOR)
    if (nu > nu_implied) {
      message(sprintf(
        paste0(
          "The implied nu at rank %d is %.3g, under the floor of %g, so the ",
          "fit uses nu = %g. Pass `nu` to use another value."
        ),
        rank,
        nu_implied,
        nu,
        nu
      ))
    }
  } else {
    nu_src <- "given"
  }

  ## Carried, not defaulted: a named argument wins over the arms'.
  cfg <- tn$run
  ## Must stay in step with mipca_tune()'s `run` list (guarded below);
  ## missing() takes a name, not a string.
  given <- c(
    chains = !missing(chains),
    iter_warmup = !missing(iter_warmup),
    iter_sampling = !missing(iter_sampling),
    thin = !missing(thin),
    spectra = !missing(spectra),
    init = !missing(init)
  )
  if (!identical(names(given), names(cfg))) {
    stop(
      "internal error: mipca_fit()'s overridable settings and the tuner's carried ",
      "configuration have drifted apart."
    )
  }
  carried <- cfg[!given]
  dots <- mget(names(given)[given], envir = environment())

  ## `call` is handed in: match.call() under do.call() would store the whole
  ## matrix. `quote = TRUE` is load-bearing: without it the language object
  ## is evaluated and this method recurses until the stack overflows.
  cl_fit <- .scrub_call(match.call(), "mipca_fit")
  res_fit <- do.call(
    .mipca_run,
    c(
      ## `abb` as the tuner resolved it: indices of this same matrix, which
      ## the fingerprint above has just required.
      list(obj, rank = rank, transform = tn$transform, nu = nu, abb = tn$abb),
      carried,
      dots,
      ## As mipca() sends it: missing() cannot be read across a frame.
      list(
        refresh = refresh,
        keep_state = keep_state,
        prop_diag = prop_diag,
        prop_diag_given = !missing(prop_diag),
        seed = seed,
        seed_given = TRUE,
        call = cl_fit
      )
    ),
    quote = TRUE
  )

  ## Which door, recorded on the fit; summary columns only.
  keep <- c(
    "rank",
    "ndraws",
    "ess_p10",
    "ess_ok",
    "rhat_max",
    "robust_z",
    "implied_nu",
    "implied_lo",
    "implied_hi"
  )
  arm <- as.data.frame(object[k, keep, drop = FALSE])
  class(arm) <- "data.frame"
  rec <- list(
    rank = rank,
    rank_src = rank_src,
    nu = nu,
    nu_src = nu_src,
    nu_implied = nu_implied,
    robust_z = rz,
    n_ho_cells = nrow(res),
    n_ho_cols = length(unique(res$col)),
    carried = names(carried),
    overridden = names(given)[given],
    arms = arm
  )
  ## keep_state = FALSE finalizes into an si_mi whose FIRST chunk is the
  ## sampler's (with a list, the kept columns; the ABB chunk follows it).
  if (inherits(res_fit, "si_mi_pca")) {
    res_fit$prov$tuned <- rec
  } else {
    res_fit$chunks[[1L]]$prov$tuned <- rec
  }
  res_fit
}

# The summary columns under their own names, and the next call. For the
# docs: row k's scored holdout cells are x$result[[k]] and its per-cell
# Rhat and ESS are x$diag[[k]]; `ess_ok` is the share of those cells with
# ESS >= attr(x, "tune")$ess_bar, `rhat_ok` the share with Rhat <=
# attr(x, "tune")$rhat_bar; `time` is elapsed seconds.
#' @export
print.si_tune <- function(x, ...) {
  .chk_no_dots(..., .what = "print()")
  tn <- attr(x, "tune")
  cols <- c(
    "rank",
    "crps",
    "ess_p10",
    "ess_ok",
    "rhat_max",
    "implied_nu",
    "implied_lo",
    "implied_hi",
    "time"
  )
  ## Degrades to the plain table print rather than refusing.
  if (!is.list(tn) || !all(c(cols, "error") %in% names(x))) {
    return(NextMethod())
  }
  ## Columns beside cells: the band's honesty turns on the column count.
  ho <- tn$na_loc[[1L]]
  n_ho_cols <- length(unique(ho[, 2L]))
  cat(sprintf(
    "si_tune: %d rank%s on a %d x %d matrix, holdout %d cells over %d column%s\n",
    nrow(x),
    .s(nrow(x)),
    tn$fp$dim[1L],
    tn$fp$dim[2L],
    nrow(ho),
    n_ho_cols,
    .s(n_ho_cols)
  ))
  ## A failed rank's message is its row's only reading, so the column is
  ## shown when one failed and left out when none did.
  if (any(!is.na(x$error))) {
    cols <- c(cols, "error")
  }
  tab <- as.data.frame(x)[cols]
  class(tab) <- "data.frame"
  print(tab, row.names = FALSE, digits = 4L)
  cat("  next: mipca_fit(x, obj) fits the rank with the lowest crps.\n")
  cat("        mipca_fit(x, obj, rank) fits another rank.\n")
  invisible(x)
}
