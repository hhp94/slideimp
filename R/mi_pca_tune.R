.chk_ranks <- function(rank) {
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
    stop(
      "`rank` has duplicates; identical arms reproduce bitwise, so a ",
      "repeated rank buys nothing.",
      call. = FALSE
    )
  }
  as.integer(rank)
}

.base_state <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    return(NULL)
  }
  get(".Random.seed", envir = globalenv(), inherits = FALSE)
}

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

.blank_cells <- function(X, idx) {
  X[idx] <- NA_real_
  X
}

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
      "holdout positions must be finite observed cells of `obj`; ",
      "a non-finite cell is a hole, not truth.",
      call. = FALSE
    )
  }
  hit <- tabulate(hc, nbins = ncol(obj))
  hitc <- which(hit > 0L)
  nobs <- n - colSums(is.na(obj[, hitc, drop = FALSE]))
  if (any(nobs - hit[hitc] < 1L)) {
    stop(
      "the holdout would blank every observed cell of at least one column.",
      call. = FALSE
    )
  }
  sort(idx)
}

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
  list(pos = pos, taken = lengths(picks), alloc = alloc)
}

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

.tune_nu_of <- function(z) {
  r <- .robust_width(z)
  list(robust_z = r, nu = .implied_nu(r))
}

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
  if (all(is.na(b))) {
    return(c(NA_real_, NA_real_))
  }
  stats::quantile(b, c(0.16, 0.84), na.rm = TRUE, names = FALSE)
}

.na_red <- function(x, f, ...) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  f(x, ..., na.rm = TRUE)
}

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

.tune_arm_read <- function(fit, rec, rank, took) {
  mh <- fit$map[rec$rows_ho, , drop = FALSE]
  md <- fit$map[rec$rows_diag, , drop = FALSE]
  if (!identical(.map_lin(mh, rec$n), rec$lin)) {
    stop(
      "internal error: the arm's holdout rows do not address the recorded cells."
    )
  }

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

  dg <- .cell_diag_read(fit$draws, rec$rows_diag, rec$run$parallel_chains)
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
    scale = max(abs(range(X, na.rm = TRUE)))
  )
}

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
      " differ ",
      "from the fingerprint recorded at tune time. The scores describe that ",
      "matrix; pass it, or re-run mipca_tune() on this one.",
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
  n_cols = NULL,
  n_rows = 4L,
  transform,
  n_diag = 10000L,
  iter_warmup = 100L,
  iter_sampling = 100L,
  chains = 4L,
  thin = 1,
  init = 1,
  refresh = 0L,
  spectra = spectra_opts(),
  parallel_chains = 1L,
  seed = NULL,
  abb = NULL
) {
  ranks <- .chk_ranks(rank)
  .chk_matrix(obj)
  abb <- .resolve_cols(abb, obj, what = "abb")
  full <- obj
  kept <- NULL
  if (!is.null(abb)) {
    kept <- .abb_kept(ncol(obj), abb)$kept
    obj <- obj[, kept, drop = FALSE]
  }
  drawn <- is.null(na_loc)
  n <- nrow(obj)
  p <- ncol(obj)
  if (as.double(n) * ncol(full) > .Machine$integer.max) {
    stop(
      "`obj` has more cells than the tuner's integer indexing can address.",
      call. = FALSE
    )
  }
  transform <- .resolve_transform(transform, n)
  n_diag <- .chk_count(n_diag, "n_diag")
  chains <- .chk_count(chains, "chains")
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  iter_warmup <- .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  thin <- as.double(.chk_count(thin, "thin"))
  init <- .chk_init(init)
  refresh <- .chk_refresh(refresh)
  spectra <- .chk_spectra(spectra)
  parallel_chains <- .chk_count(parallel_chains, "parallel_chains")
  cen <- .check_admissible(obj)
  if (!is.null(seed)) {
    seed <- .chk_count(seed, "`seed`", lo = 0L)
  }

  na0 <- which(is.na(obj))
  if (!is.null(n_cols)) {
    n_cols <- .chk_count(n_cols, "n_cols")
  }
  if (is.null(na_loc)) {
    if (is.null(n_cols)) {
      n_cols <- max(50L, min(500L, as.integer(floor(0.10 * p))))
    }
    n_rows <- .chk_count(n_rows, "n_rows")
    if (any(cen$const)) {
      k <- which(cen$const)
      stop(
        "`obj` has ",
        length(k),
        " constant column",
        .s(length(k)),
        " [",
        fmt_trunc(.col_label(obj, k), 8),
        "]. The sampler fits them, but the default holdout ",
        "draw needs two distinct observed values in every column of the ",
        "matrix; drop them before tuning, or supply `na_loc`.",
        call. = FALSE
      )
    }
    nobs <- n - cen$cmiss
    eligible <- sum(nobs >= n_rows + 2L)
    if (eligible < n_cols) {
      stop(
        "cannot place a holdout over ",
        n_cols,
        " columns at ",
        n_rows,
        " cells each: that needs ",
        n_cols,
        " columns able to give up ",
        n_rows,
        " observed cells (keeping two) and the matrix has ",
        eligible,
        "; lower `n_cols` or `n_rows`, or supply `na_loc`.",
        call. = FALSE
      )
    }
    n_ho_cells <- n_cols * n_rows
  } else {
    if (!missing(n_rows)) {
      stop(
        "`n_rows` has no effect with a supplied `na_loc`: the design's own ",
        "rows say how many cells each column gives up. Drop one.",
        call. = FALSE
      )
    }
    if (!is.list(na_loc) || length(na_loc) < 1L) {
      stop(
        "`na_loc` must be a list of one two-column (row, col) matrix -- the ",
        "shape sample_na_loc() returns at n_reps = 1.",
        call. = FALSE
      )
    }
    if (length(na_loc) > 1L) {
      stop(
        "`na_loc` holds ",
        length(na_loc),
        " designs and the tuner scores one. ",
        "The rep axis is retired: loop over the designs yourself, one ",
        "mipca_tune() call each, and compare the tables -- two designs score ",
        "different cells, so their rows are not one grid.",
        call. = FALSE
      )
    }
    n_ho_cells <- NROW(na_loc[[1L]])
  }

  if (n_ho_cells <= 2L) {
    stop(
      "the holdout has ",
      n_ho_cells,
      " cell(s); the tuner needs more than ",
      "two before a CRPS column means anything. Raise `n_cols` or `n_rows`, ",
      "or supply an `na_loc` with more rows.",
      call. = FALSE
    )
  }

  resolve <- function(design, against) {
    tryCatch(
      .resolve_holdout(against, design),
      error = function(e) {
        stop("na_loc: ", conditionMessage(e), call. = FALSE)
      }
    )
  }
  Xb <- NULL
  if (!is.null(na_loc)) {
    ho_lin <- resolve(na_loc[[1L]], full)
    if (!is.null(kept)) {
      ho_lin <- .abb_tune_rebase(ho_lin, full, kept)
    }
    if (
      !is.null(n_cols) &&
        length(unique(na_loc[[1L]][, 2L])) != n_cols
    ) {
      stop(
        "`n_cols` disagrees with the supplied `na_loc`; drop one.",
        call. = FALSE
      )
    }
    Xb <- .blank_cells(obj, ho_lin)
    tryCatch(
      .check_admissible(Xb),
      error = function(e) {
        stop(
          "na_loc: after blanking the holdout, ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

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
    dq_old <- .seed_dq(NULL)
    on.exit(dqrng::dqrng_set_state(dq_old), add = TRUE)
  }
  if (is.null(na_loc)) {
    na_loc <- sample_na_loc(
      obj,
      n_cols = n_cols,
      n_rows = n_rows,
      n_reps = 1L,
      rowmax = 1,
      colmax = 1
    )
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
  dd <- .tune_diag_draw(hp$holes, n, p, n_diag)
  st <- .tune_store(hp$ho_pos, dd$pos)
  hcol <- .lin_col(ho_lin, n)
  ucol <- unique(hcol)
  inv <- integer(p)
  inv[ucol] <- seq_along(ucol)
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
    u = stats::runif(length(ho_lin)),
    ess_bar = 100L * chains,
    rhat_bar = 1.05,
    run = list(
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      thin = thin,
      spectra = spectra,
      parallel_chains = parallel_chains,
      init = init
    )
  )
  if (is.null(Xb)) {
    Xb <- .blank_cells(obj, ho_lin)
  }

  rows <- vector("list", length(ranks))
  for (i in seq_along(ranks)) {
    t0 <- proc.time()[["elapsed"]]
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
        spectra = spectra,
        parallel_chains = parallel_chains
      ),
      error = identity
    )
    took <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      rows[[i]] <- .tune_row(ranks[i], took, error = conditionMessage(fit))
      next
    }
    rows[[i]] <- .tune_arm_read(fit, rec, ranks[i], took)
    rm(fit)
  }

  out <- do.call(rbind, rows)
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
  attr(out, "tune") <- list(
    transform = transform,
    run = rec$run,
    fp = .tune_fp(full),
    na_loc = na_loc,
    rank = ranks,
    n_diag = n_diag,
    ess_bar = rec$ess_bar,
    rhat_bar = rec$rhat_bar,
    abb = abb
  )
  class(out) <- c("si_tune", "slideimp_tbl", "data.frame")
  out
}

.tune_record <- function(object) {
  tn <- attr(object, "tune")
  if (
    !is.list(tn) ||
      is.null(tn$fp) ||
      is.null(tn$run) ||
      is.null(tn$transform) ||
      is.null(tn$rank)
  ) {
    stop(
      "this table has lost the record mipca_fit() hands off from. A row ",
      "subset keeps it; a column subset or subset() drops it. Hand off from ",
      "the table mipca_tune() returned, or a row subset of it.",
      call. = FALSE
    )
  }
  r <- object$rank
  if (!is.numeric(r) || anyDuplicated(r) || !all(r %in% tn$rank)) {
    stop(
      "this table's rows are not all from the run its record describes: ",
      "rbind() keeps the first table's record only, so rows from another run ",
      "would be fitted under it. Hand off from one run's table.",
      call. = FALSE
    )
  }
  tn
}

.tune_argmin_rank <- function(object) {
  live <- is.na(object$error) & is.finite(object$crps)
  if (!any(live)) {
    stop(
      "no arm recorded any holdout scores, so there is no grid to choose ",
      "from; read `x$error` for why, or pass `rank` explicitly.",
      call. = FALSE
    )
  }
  r <- object$rank[live]
  r[order(object$crps[live], r)[1L]]
}

#' @export
mipca_fit <- function(object, ...) {
  UseMethod("mipca_fit")
}

#' @export
mipca_fit.si_tune <- function(
  object,
  obj,
  rank = NULL,
  ...,
  chains,
  iter_warmup,
  iter_sampling,
  thin,
  spectra,
  parallel_chains,
  init,
  refresh = 100L,
  keep_state = FALSE,
  prop_diag = .prop_diag_default(),
  seed = NULL,
  nu = NULL,
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
    stop(
      "mipca_fit() does not take ",
      paste(sprintf("`%s`", hit), collapse = ", "),
      ". The transform is carried across from the tuner and is yours to ",
      "choose in a direct mipca() call; a holdout is not, because scoring ",
      "belongs to the tuner and a production fit must not be conditioned on ",
      "deliberately degraded data. The `abb` list is carried across with the ",
      "transform: rank and nu were tuned with those columns left out.",
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
  k <- which(object$rank == rank)
  if (length(k) != 1L) {
    stop(
      "rank ",
      rank,
      " is not a tuned arm (the table has ",
      paste(sort(object$rank), collapse = ", "),
      "); nu is read off the ",
      "chosen rank's Gaussian arm, so an untuned rank goes through ",
      "mipca() directly, or through a second tuner call sharing `na_loc`.",
      call. = FALSE
    )
  }
  res <- object$result[[k]]
  if (is.null(res) || nrow(res) < 1L) {
    err <- object$error[k]
    stop(
      "the rank-",
      rank,
      " arm recorded no scores",
      if (!is.na(err)) paste0(" (", err, ")") else "",
      "; nothing to read nu from.",
      call. = FALSE
    )
  }
  if (!identical(obj[cbind(res$row, res$col)], res$truth)) {
    stop(
      "`obj` disagrees with the rank-",
      rank,
      " arm's recorded truth at its ",
      "holdout cells; the scores describe another matrix. Pass the matrix ",
      "this table was tuned on, or re-run mipca_tune().",
      call. = FALSE
    )
  }

  if (is.null(res$truth_lat)) {
    stop(
      "this table carries no latent-scale cells (it predates the latent nu ",
      "read of 2026-10-29); re-run mipca_tune().",
      call. = FALSE
    )
  }
  nuz <- .tune_nu_of((res$truth_lat - res$mean_lat) / res$sd_lat)
  rz <- nuz$robust_z
  nu_implied <- nuz$nu
  if (is.null(nu)) {
    nu_src <- "tuner"
    if (is.na(nu_implied)) {
      stop(
        "the implied nu off the rank-",
        rank,
        " arm is unreadable (fewer ",
        "than 4 finite z); pass `nu` to override the reading.",
        call. = FALSE
      )
    }
    nu <- max(nu_implied, .NU_FLOOR)
    if (nu > nu_implied) {
      message(sprintf(
        "mipca_fit(): nu %.3g read off the rank-%d arm, %g used (the floor; pass `nu` to override).",
        nu_implied,
        rank,
        nu
      ))
    }
  } else {
    nu_src <- "given"
  }

  cfg <- tn$run
  given <- c(
    chains = !missing(chains),
    iter_warmup = !missing(iter_warmup),
    iter_sampling = !missing(iter_sampling),
    thin = !missing(thin),
    spectra = !missing(spectra),
    parallel_chains = !missing(parallel_chains),
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

  cl_fit <- .scrub_call(match.call(), "mipca_fit")
  res_fit <- do.call(
    .mipca_run,
    c(
      list(obj, rank = rank, transform = tn$transform, nu = nu, abb = tn$abb),
      carried,
      dots,
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
  rec$line <- .tuned_line(rec)
  if (inherits(res_fit, "si_mi_pca")) {
    res_fit$prov$tuned <- rec
  } else {
    res_fit$chunks[[1L]]$prov$tuned <- rec
  }
  res_fit
}

.tuned_line <- function(tp) {
  lo <- tp$arms$implied_lo
  hi <- tp$arms$implied_hi
  nu_txt <- if (identical(tp$nu_src, "given")) {
    sprintf("nu %.4g given", tp$nu)
  } else if (tp$nu > tp$nu_implied) {
    sprintf("nu %.3g read, %.4g used", tp$nu_implied, tp$nu)
  } else {
    sprintf("nu %.4g", tp$nu)
  }
  paste0(
    sprintf(
      "  tuned:   rank %d (%s), %s; holdout %d cells over %d columns, implied-nu band %.3g-%.3g\n",
      tp$rank,
      tp$rank_src,
      nu_txt,
      tp$n_ho_cells,
      tp$n_ho_cols,
      .na_red(lo, min),
      .na_red(hi, max)
    ),
    sprintf(
      "           carried from the tuner run: %s%s\n",
      if (length(tp$carried) > 0L) {
        paste(tp$carried, collapse = ", ")
      } else {
        "none"
      },
      if (length(tp$overridden) > 0L) {
        paste0("; overridden here: ", paste(tp$overridden, collapse = ", "))
      } else {
        ""
      }
    )
  )
}

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
    "time",
    "error"
  )
  if (!is.list(tn) || !all(cols %in% names(x))) {
    cat(
      "si_tune: the tuning record or a summary column is missing, so this prints as a plain table\n"
    )
    return(NextMethod())
  }
  cat(sprintf(
    "si_tune: %d arm%s on a %d x %d matrix, %s transform\n",
    nrow(x),
    .s(nrow(x)),
    tn$fp$dim[1L],
    tn$fp$dim[2L],
    tn$transform$name
  ))
  ho <- tn$na_loc[[1L]]
  n_ho_cols <- length(unique(ho[, 2L]))
  cat(sprintf(
    "  holdout: %d cells over %d column%s; Rhat and ESS read on up to %d cells per arm\n",
    nrow(ho),
    n_ho_cols,
    .s(n_ho_cols),
    tn$n_diag
  ))
  cat(sprintf(
    "  essOK:   share of those cells with ESS >= %d; x$rhat_ok is the share with Rhat <= %.2f\n",
    tn$ess_bar,
    tn$rhat_bar
  ))
  if (!is.null(tn$abb)) {
    cat(sprintf(
      "  abb:     %d column%s listed and left out of every arm; mipca_fit() draws %s by ABB\n",
      length(tn$abb),
      .s(length(tn$abb)),
      if (length(tn$abb) == 1L) "it" else "them"
    ))
  }
  cat(sprintf(
    "  %-5s %9s %7s %7s %7s %8s %8s %7s\n",
    "rank",
    "crps",
    "ess10",
    "essOK",
    "rhatMx",
    "impl.nu",
    "band",
    "time.s"
  ))
  for (k in seq_len(nrow(x))) {
    if (!is.na(x$error[k])) {
      cat(sprintf("  %-5d ERROR: %s\n", x$rank[k], x$error[k]))
      next
    }
    cat(sprintf(
      "  %-5d %9.5f %7.0f %7.2f %7.3f %8.2f %3.1f-%-4.1f %7.1f\n",
      x$rank[k],
      x$crps[k],
      x$ess_p10[k],
      x$ess_ok[k],
      x$rhat_max[k],
      x$implied_nu[k],
      x$implied_lo[k],
      x$implied_hi[k],
      x$time[k]
    ))
  }
  cat(
    "  per-cell detail for arm k: x$result[[k]] (scored holdout cells), x$diag[[k]] (Rhat and ESS)\n"
  )
  cat("  next:    mipca_fit(x, obj) fits the rank with the lowest crps\n")
  cat("           mipca_fit(x, obj, rank) fits another rank\n")
  invisible(x)
}
