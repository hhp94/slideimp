.chk_parameters <- function(parameters) {
  if (!is.data.frame(parameters) || nrow(parameters) < 1L) {
    stop(
      "`parameters` must be a data.frame with at least one row.",
      call. = FALSE
    )
  }
  if (!identical(names(parameters), "rank")) {
    stop(
      "`parameters` must have exactly one column, `rank`; nu is not an axis ",
      "here (every arm is Gaussian and `implied_nu` comes back derived).",
      call. = FALSE
    )
  }
  r <- parameters$rank
  if (
    !is.numeric(r) || anyNA(r) || any(!is.finite(r)) ||
      any(r < 1) || any(r != trunc(r))
  ) {
    stop("`parameters$rank` must be whole numbers >= 1, no NA.", call. = FALSE)
  }
  if (anyDuplicated(r)) {
    stop(
      "`parameters$rank` has duplicates; identical arms reproduce bitwise, ",
      "so a repeated row buys nothing.",
      call. = FALSE
    )
  }
  as.integer(r)
}

.dq_state <- function() {
  dqrng::dqrng_get_state()
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

.tune_grid_row <- function(param_set, rank, took, error = NA_character_,
                           ndraws = NA_integer_, n_diag = NA_integer_,
                           ess_min = NA_real_, ess_p10 = NA_real_,
                           ess_med = NA_real_, ess_ok = NA_real_,
                           rhat_med = NA_real_, rhat_p90 = NA_real_,
                           rhat_max = NA_real_, rhat_ok = NA_real_,
                           robust_z = NA_real_, implied_nu = NA_real_,
                           implied_lo = NA_real_, implied_hi = NA_real_,
                           sigma = list(NULL)) {
  data.frame(
    param_set = param_set, rank = rank,
    ndraws = ndraws, n_diag = n_diag,
    ess_min = ess_min, ess_p10 = ess_p10, ess_med = ess_med,
    ess_ok = ess_ok,
    rhat_med = rhat_med, rhat_p90 = rhat_p90, rhat_max = rhat_max,
    rhat_ok = rhat_ok,
    robust_z = robust_z, implied_nu = implied_nu,
    implied_lo = implied_lo, implied_hi = implied_hi,
    sigma = I(sigma),
    time = took, error = error,
    stringsAsFactors = FALSE
  )
}

.tune_arm_read <- function(fit, prov, i, rank, took) {
  rd <- prov$hold
  n <- prov$dim[[1L]]
  mh <- fit$map[rd$rows_ho, , drop = FALSE]
  md <- fit$map[rd$rows_diag, , drop = FALSE]
  if (
    !identical(fit$holdout$rows, rd$rows_ho) ||
      !identical(.map_lin(mh, n), rd$lin)
  ) {
    stop("internal error: the fit's holdout rows disagree with the store construction.")
  }
  if (!identical(fit$holdout$truth, rd$truth)) {
    stop("internal error: the fit scored different truth than the tuner recorded.")
  }

  sc <- fit$holdout
  sdv <- sqrt(sc$var)
  z <- (rd$truth - sc$mean) / sdv
  nuz <- .tune_nu_of(z)
  band <- .tune_nu_band(z, rd$hcol)

  dg <- .cell_diag_read(fit$draws, rd$rows_diag, prov$run$parallel_chains)
  rhat <- dg[, 1L]
  ess <- dg[, 2L]
  list(
    cells = data.frame(
      param_set = i, rank = rank,
      col = mh[, "col"],
      row = mh[, "row"],
      truth = rd$truth, mean = sc$mean, sd = sdv,
      crps = sc$crps, pit = sc$pit,
      stringsAsFactors = FALSE
    ),
    diag = data.frame(
      param_set = i, rank = rank,
      col = md[, "col"],
      row = md[, "row"],
      rhat = rhat, ess = ess,
      stringsAsFactors = FALSE
    ),
    grid = .tune_grid_row(
      param_set = i, rank = rank, took = took,
      ndraws = length(fit$draws) * ncol(fit$draws[[1L]]),
      n_diag = length(rd$rows_diag),
      ess_min = .na_red(ess, min),
      ess_p10 = .na_red(ess, stats::quantile, 0.10, names = FALSE),
      ess_med = stats::median(ess, na.rm = TRUE),
      ess_ok = mean(ess >= prov$ess_bar, na.rm = TRUE),
      rhat_med = stats::median(rhat, na.rm = TRUE),
      rhat_p90 = stats::quantile(rhat, 0.90, na.rm = TRUE, names = FALSE),
      rhat_max = .na_red(rhat, max),
      rhat_ok = mean(rhat <= prov$rhat_bar, na.rm = TRUE),
      robust_z = nuz$robust_z,
      implied_nu = nuz$nu,
      implied_lo = band[[1L]],
      implied_hi = band[[2L]],
      sigma = list(vapply(fit$seam$chain, `[[`, 0, "sigma"))
    )
  )
}

#' @export
mipca_tune <- function(
  obj,
  parameters,
  na_loc = NULL,
  n_cols = NULL,
  n_rows = 4L,
  transform,
  level = c(0.5, 0.9, 0.99),
  n_diag = 10000L,
  iter_warmup = 100L,
  iter_sampling = 100L,
  chains = 4L,
  thin = 1,
  init = 1,
  colmax = 0.95,
  refresh = 0L,
  spectra = spectra_opts(),
  parallel_chains = 1L,
  seed = NULL
) {
  ranks <- .chk_parameters(parameters)
  .chk_matrix(obj)
  n <- nrow(obj)
  p <- ncol(obj)
  if (as.double(n) * p > .Machine$integer.max) {
    stop("`obj` has more cells than the tuner's integer indexing can address.",
         call. = FALSE)
  }
  transform <- .resolve_transform(transform, n)
  n_diag <- .chk_count(n_diag, "n_diag")
  chains <- .chk_count(chains, "chains")
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  iter_warmup <- .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  thin <- as.double(.chk_count(thin, "thin"))
  init <- .chk_init(init)
  colmax <- .chk_colmax(colmax)
  refresh <- .chk_refresh(refresh)
  spectra <- .chk_spectra(spectra)
  lvl <- as.double(level)
  if (length(lvl) < 1L || anyNA(lvl) || any(lvl <= 0 | lvl >= 1)) {
    stop("`level` must be fractions strictly inside (0, 1).", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .chk_count(seed, "`seed`", lo = 0L)
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
    nobs <- n - tabulate(.lin_col(na0, n), nbins = p)
    eligible <- sum(nobs >= n_rows + 2L & (n - nobs + n_rows) <= 0.9 * n)
    if (eligible < n_cols) {
      stop(
        "cannot place a holdout over ", n_cols, " columns at ", n_rows,
        " cells each: that needs ", n_cols, " columns able to give up ",
        n_rows, " observed cells (keeping two, staying under 90% missing) ",
        "and the matrix has ", eligible, "; lower `n_cols` or `n_rows`, ",
        "or supply `na_loc`.",
        call. = FALSE
      )
    }
    na_loc <- sample_na_loc(
      obj,
      n_cols = n_cols,
      n_rows = n_rows,
      n_reps = 1L
    )
  } else {
    if (!is.list(na_loc) || length(na_loc) < 1L) {
      stop(
        "`na_loc` must be a list of one two-column (row, col) matrix -- the ",
        "shape sample_na_loc() returns at n_reps = 1.",
        call. = FALSE
      )
    }
    if (length(na_loc) > 1L) {
      stop(
        "`na_loc` holds ", length(na_loc), " designs and the tuner scores one. ",
        "The rep axis is retired: loop over the designs yourself, one ",
        "mipca_tune() call each, and compare the objects -- their grids ",
        "are not poolable (see c.si_tune()).",
        call. = FALSE
      )
    }
    if (
      !is.null(n_cols) &&
        !all(
          vapply(na_loc, function(a) length(unique(a[, 2L])), 0L) == n_cols
        )
    ) {
      stop("`n_cols` disagrees with the supplied `na_loc`; drop one.",
           call. = FALSE)
    }
  }

  n_ho_cells <- NROW(na_loc[[1L]])
  if (n_ho_cells <= 2L) {
    stop(
      "the holdout has ", n_ho_cells, " cell(s); the tuner needs more than ",
      "two before a CRPS column means anything. Raise `n_cols` or `n_rows`, ",
      "or supply an `na_loc` with more rows.",
      call. = FALSE
    )
  }

  dq_state <- .dq_state()

  narm <- length(ranks)
  fits <- vector("list", narm)
  grid_rows <- vector("list", narm)
  cells_rows <- vector("list", narm)
  diag_rows <- vector("list", narm)

  ho_mat <- na_loc[[1L]]
  ho_lin <- tryCatch(
    .resolve_holdout(obj, ho_mat),
    error = function(e) {
      stop("na_loc: ", conditionMessage(e), call. = FALSE)
    }
  )
  truth <- obj[ho_lin]
  hcol <- .lin_col(ho_lin, n)
  hp <- .tune_hole_pos(na0, ho_lin)
  dd <- .tune_diag_draw(hp$holes, n, p, n_diag)
  st <- .tune_store(hp$ho_pos, dd$pos)
  u_ho <- stats::runif(length(ho_lin))
  ho_rec <- list(
    lin = ho_lin,
    truth = truth,
    hcol = hcol,
    diag_pos = dd$pos,
    diag_taken = dd$taken,
    rows_ho = st$rows_ho,
    rows_diag = st$rows_diag,
    store_n = length(st$store),
    u = u_ho
  )

  prov <- list(
    n_diag_req = n_diag,
    ess_bar = 100L * chains,
    rhat_bar = 1.05,
    run = list(
      chains = chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      thin = thin,
      spectra = spectra,
      parallel_chains = parallel_chains,
      init = init,
      colmax = colmax,
      refresh = refresh
    ),
    dim = dim(obj),
    dimnames = dimnames(obj),
    dq_state = dq_state,
    hold = ho_rec
  )

  for (i in seq_along(ranks)) {
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      .mipca_run(
        obj,
        rank = ranks[i],
        transform = transform,
        store = st$store,
        holdout = ho_mat,
        level = lvl,
        u = u_ho,
        keep_state = TRUE,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        chains = chains,
        thin = thin,
        init = init,
        colmax = colmax,
        refresh = refresh,
        spectra = spectra,
        parallel_chains = parallel_chains
      ),
      error = identity
    )
    took <- proc.time()[["elapsed"]] - t0

    if (inherits(fit, "error")) {
      grid_rows[[i]] <- .tune_grid_row(
        i, ranks[i], took, error = conditionMessage(fit)
      )
      next
    }

    ar <- .tune_arm_read(fit, prov, i, ranks[i], took)
    fits[[i]] <- fit
    cells_rows[[i]] <- ar$cells
    diag_rows[[i]] <- ar$diag
    grid_rows[[i]] <- ar$grid
  }

  e <- new.env(parent = emptyenv())
  .tune_tables_set(
    e,
    do.call(rbind, grid_rows),
    do.call(rbind, cells_rows),
    do.call(rbind, diag_rows)
  )
  e$na_loc <- na_loc
  e$data <- obj
  e$fits <- fits
  e$transform <- transform
  e$level <- lvl
  prov$rng_end <- .base_state()
  e$prov <- prov
  e$call <- list(.scrub_call(match.call(), "mipca_tune"))
  class(e) <- "si_tune"
  e
}

#' @export
c.si_tune <- function(...) {
  stop(
    "tuners cannot be combined. Two runs share a scoring basis only if they ",
    "share the blanked cells and the PIT jitter, and two calls share ",
    "neither -- so an argmin over the joined grid would compare columns that ",
    "mean different things. Add the ranks to ONE tuner instead: ",
    "update(tuner, parameters = data.frame(rank = <the new ones>)), or score ",
    "them all in a single mipca_tune() call.",
    call. = FALSE
  )
}

.rn_auto <- function(x) {
  row.names(x) <- NULL
  x
}

.tune_tables_set <- function(x, grid, cells, diag) {
  x$grid <- .rn_auto(grid)
  x$cells <- .rn_auto(cells)
  x$diag <- .rn_auto(diag)
  x
}

.tune_blocks <- function(tab, g, lens, what) {
  out <- vector("list", nrow(g))
  off <- 0L
  for (k in seq_len(nrow(g))) {
    if (!is.na(g$error[k]) || lens[k] < 1L) {
      next
    }
    blk <- tab[off + seq_len(lens[k]), , drop = FALSE]
    if (!all(blk$param_set == g$param_set[k])) {
      stop("internal error: a ", what, " block does not belong to its grid row.")
    }
    out[[k]] <- blk
    off <- off + lens[k]
  }
  if (off != nrow(tab)) {
    stop("internal error: ", what, " has rows outside every arm's block.")
  }
  out
}

.tune_rng_scope <- function(object, expr) {
  if (is.null(object$prov$dq_state) || is.null(object$prov$rng_end)) {
    stop(
      "this tuner did not record the RNG states a mutation has to restore, ",
      "so it could not reproduce what a fresh call would have drawn; ",
      "re-run mipca_tune().",
      call. = FALSE
    )
  }
  caller_dq <- .dq_state()
  caller_base <- .base_state()
  on.exit({
    dqrng::dqrng_set_state(caller_dq)
    .base_restore(caller_base)
  }, add = TRUE)
  dqrng::dqrng_set_state(object$prov$dq_state)
  .base_restore(object$prov$rng_end)
  res <- force(expr)
  object$prov$rng_end <- .base_state()
  res
}

.tune_grow <- function(object, parameters, cl) {
  ranks_new <- .chk_parameters(parameters)
  g <- object$grid
  if (is.unsorted(g$param_set, strictly = TRUE)) {
    stop("internal error: this tuner's param_set is not strictly ascending.")
  }
  dup <- intersect(ranks_new, unique(g$rank))
  if (length(dup) > 0L) {
    stop(
      "rank ", paste(dup, collapse = ", "), " is already in this grid. ",
      "Re-running a rank on the same cells with the same streams reproduces ",
      "it bitwise, so the second row would be a duplicate rather than ",
      "information.",
      call. = FALSE
    )
  }
  X <- .tune_data(object, "grow")
  cfg <- .tune_run_config(object)

  na0 <- which(is.na(X))
  cell_lens <- rep(length(object$prov$hold$lin), nrow(g))
  invisible(.tune_blocks(object$cells, g, cell_lens, "cells"))
  invisible(.tune_blocks(object$diag, g, g$n_diag, "diag"))
  grid_new <- list()
  cells_new <- list()
  diag_new <- list()
  fits_new <- list()
  ps0 <- max(g$param_set)

  .tune_rng_scope(object, {
  rd <- object$prov$hold
  hp <- .tune_hole_pos(na0, rd$lin)
  st <- .tune_store(hp$ho_pos, rd$diag_pos)
  if (
    length(st$store) != rd$store_n ||
      !identical(st$rows_ho, rd$rows_ho) ||
      !identical(st$rows_diag, rd$rows_diag)
  ) {
    stop(
      "internal error: the rebuilt store set disagrees with the one this ",
      "tuner's arms were run on."
    )
  }
  for (i in seq_along(ranks_new)) {
    ps <- ps0 + i
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      do.call(
        .mipca_run,
        c(
          list(
            X,
            rank = ranks_new[i],
            transform = object$transform,
            store = st$store,
            holdout = object$na_loc[[1L]],
            level = object$level,
            u = rd$u,
            keep_state = TRUE
          ),
          cfg
        )
      ),
      error = identity
    )
    took <- proc.time()[["elapsed"]] - t0
    if (inherits(fit, "error")) {
      grid_new <- c(grid_new, list(
        .tune_grid_row(ps, ranks_new[i], took, error = conditionMessage(fit))
      ))
      cells_new <- c(cells_new, list(NULL))
      diag_new <- c(diag_new, list(NULL))
      fits_new <- c(fits_new, list(NULL))
      next
    }
    ar <- .tune_arm_read(fit, object$prov, ps, ranks_new[i], took)
    grid_new <- c(grid_new, list(ar$grid))
    cells_new <- c(cells_new, list(ar$cells))
    diag_new <- c(diag_new, list(ar$diag))
    fits_new <- c(fits_new, list(fit))
  }
  })

  .tune_tables_set(
    object,
    rbind(g, do.call(rbind, grid_new)),
    rbind(object$cells, do.call(rbind, cells_new)),
    rbind(object$diag, do.call(rbind, diag_new))
  )
  object$fits <- c(object$fits, fits_new)
  object$call <- c(object$call, list(cl))
  invisible(object)
}

#' @export
update.si_tune <- function(object, iter_sampling, parameters, ...) {
  cl <- .scrub_call(match.call(), "update")
  if (...length() > 0L) {
    nm <- ...names()
    nm <- nm[nzchar(nm)]
    stop(
      "update() on a tuner takes no argument beyond `iter_sampling` or ",
      "`parameters`",
      if (length(nm) > 0L) paste0(" (got ", fmt_trunc(nm, 8), ")") else "",
      ".\nEvery arm extends the way its chains sampled, and every grown arm ",
      "runs the way the arms beside it did.",
      call. = FALSE
    )
  }
  have_iter <- !missing(iter_sampling)
  have_par <- !missing(parameters)
  if (have_iter && have_par) {
    stop(
      "update() extends OR grows, not both in one call. Extending re-reads ",
      "every arm at a new length while growing adds arms at the current one, ",
      "so together they would leave the grid half-described by each. Call it ",
      "twice.",
      call. = FALSE
    )
  }
  if (!have_iter && !have_par) {
    stop(
      "update() on a tuner needs either `iter_sampling`, to extend every arm, ",
      "or `parameters`, to add ranks scored on the cells this tuner already ",
      "fixed.",
      call. = FALSE
    )
  }
  if (have_par) {
    return(.tune_grow(object, parameters, cl))
  }
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")

  g <- object$grid
  narm <- nrow(g)
  grid_rows <- lapply(seq_len(narm), function(k) g[k, , drop = FALSE])
  cells_rows <- vector("list", narm)
  diag_rows <- vector("list", narm)

  .tune_rng_scope(object, {
  for (k in seq_len(narm)) {
    fit <- object$fits[[k]]
    if (is.null(fit)) {
      next
    }
    t0 <- proc.time()[["elapsed"]]
    err <- tryCatch(
      {
        update(fit, iter_sampling)
        NULL
      },
      error = identity
    )
    took <- g$time[k] + (proc.time()[["elapsed"]] - t0)
    if (!is.null(err)) {
      grid_rows[[k]] <- .tune_grid_row(
        g$param_set[k], g$rank[k], took, error = conditionMessage(err)
      )
      object$fits[k] <- list(NULL)
      next
    }
    ar <- .tune_arm_read(fit, object$prov, g$param_set[k], g$rank[k], took)
    cells_rows[[k]] <- ar$cells
    diag_rows[[k]] <- ar$diag
    grid_rows[[k]] <- ar$grid
  }
  })

  .tune_tables_set(
    object,
    do.call(rbind, grid_rows),
    do.call(rbind, cells_rows),
    do.call(rbind, diag_rows)
  )
  object$call <- c(object$call, list(cl))
  invisible(object)
}

.tune_argmin_rank <- function(object) {
  cl <- object$cells
  if (is.null(cl) || nrow(cl) < 1L) {
    stop(
      "no arm recorded any holdout scores, so there is no grid to choose ",
      "from; read `x$grid$error` for why, or pass `rank` explicitly.",
      call. = FALSE
    )
  }
  m <- tapply(cl$crps, cl$rank, mean)
  as.integer(names(m)[which.min(m)])
}

.tune_run_config <- function(object) {
  g <- object$grid
  live <- !vapply(object$fits, is.null, NA)
  if (length(live) != nrow(g)) {
    stop(
      "internal error: the tuner holds ", length(live), " arms and ",
      nrow(g), " grid rows."
    )
  }
  if (!identical(live, is.na(g$error))) {
    stop(
      "internal error: the tuner's live arms and its error-free grid rows ",
      "are not the same arms."
    )
  }
  fits <- object$fits[live]
  if (length(fits) < 1L) {
    stop(
      "no arm of this tuner is live, so there is no configuration to carry; ",
      "read `x$grid$error` for why.",
      call. = FALSE
    )
  }
  out <- object$prov$run
  same <- function(a, b) {
    if (is.numeric(a) && is.numeric(b)) {
      identical(as.double(a), as.double(b))
    } else {
      identical(a, b)
    }
  }
  for (f in fits) {
    got <- list(
      chains = f$prov$chains,
      iter_warmup = f$prov$warmup,
      thin = f$prov$thin,
      spectra = f$seam$spectra
    )
    bad <- names(got)[!mapply(same, got, out[names(got)])]
    if (length(bad) > 0L) {
      stop(
        "internal error: an arm did not run the way the tuner asked (",
        fmt_trunc(bad, 8), ")."
      )
    }
  }
  nd <- vapply(fits, function(f) f$prov$ndraws, integer(1L))
  if (length(unique(nd)) != 1L) {
    stop(
      "internal error: the tuner's arms are not a single length, so there ",
      "is no single configuration to carry forward."
    )
  }
  out$iter_sampling <- nd[[1L]]
  out
}

.tune_data <- function(object, verb) {
  X <- object$data
  if (is.null(X)) {
    stop(
      "this tuner no longer carries the matrix it was tuned on, so there is ",
      "nothing to ", verb, "; re-run mipca_tune().",
      call. = FALSE
    )
  }
  if (
    !identical(dim(X), object$prov$dim) ||
      !identical(dimnames(X), object$prov$dimnames)
  ) {
    stop(
      "this tuner's matrix is not the one it was tuned on: dimensions or ",
      "dimnames disagree with what was recorded at tune time. Its scores ",
      "describe a matrix that is gone; re-run mipca_tune().",
      call. = FALSE
    )
  }
  rd <- object$prov$hold
  if (!identical(X[rd$lin], rd$truth)) {
    stop(
      "this tuner's matrix has been edited since it was tuned: values at ",
      "the recorded holdout cells differ from the truth stored at tune ",
      "time. Its scores describe a matrix that is gone; re-run ",
      "mipca_tune().",
      call. = FALSE
    )
  }
  X
}

#' @export
mipca_fit <- function(object, ...) {
  UseMethod("mipca_fit")
}

#' @export
mipca_fit.si_tune <- function(
  object,
  rank = NULL,
  ...,
  chains,
  iter_warmup,
  iter_sampling,
  thin,
  spectra,
  parallel_chains,
  init,
  colmax,
  refresh,
  keep_state = FALSE,
  prop_diag = .prop_diag_default(),
  seed = NULL,
  nu,
  transform,
  holdout,
  data
) {
  .chk_no_dots(..., .what = "mipca_fit()")
  owned <- c(
    nu = !missing(nu),
    transform = !missing(transform),
    holdout = !missing(holdout)
  )
  if (any(owned)) {
    hit <- names(owned)[owned]
    stop(
      "mipca_fit() carries ", paste(sprintf("`%s`", hit), collapse = ", "),
      " across from the tuner. `nu` and `transform` are yours to choose in a ",
      "direct mipca() call; `holdout` is not, because scoring belongs to ",
      "the tuner and a production fit must not be conditioned on ",
      "deliberately degraded data.",
      call. = FALSE
    )
  }
  if (!missing(data)) {
    stop(
      "`data` is no longer an argument to mipca_fit(). A tuner carries the matrix ",
      "it was tuned on and nothing in this package detaches it, so there is ",
      "no re-attach to make; if that matrix is missing or altered, the ",
      "tuner's scores describe something that is gone and the answer is to ",
      "re-run mipca_tune().",
      call. = FALSE
    )
  }
  if (is.null(rank)) {
    rank_src <- "argmin CRPS"
    rank <- .tune_argmin_rank(object)
  } else {
    rank_src <- "given"
    rank <- .chk_count(rank, "rank")
  }
  g <- object$grid
  if (!any(g$rank == rank)) {
    stop(
      "rank ", rank, " is not a tuned arm (the grid has ",
      paste(sort(unique(g$rank)), collapse = ", "), "); nu is read off the ",
      "chosen rank's Gaussian arm, so an untuned rank goes through ",
      "mipca() directly, or through a second tuner call sharing `na_loc`.",
      call. = FALSE
    )
  }

  X <- .tune_data(object, "fit")

  cl <- object$cells
  sel <- if (is.null(cl)) logical(0L) else cl$rank == rank
  if (!any(sel)) {
    errs <- g$error[g$rank == rank]
    errs <- errs[!is.na(errs)]
    stop(
      "the rank-", rank, " arm recorded no scores",
      if (length(errs) > 0L) {
        paste0(" (", paste(unique(errs), collapse = "; "), ")")
      } else {
        ""
      },
      "; nothing to read nu from.",
      call. = FALSE
    )
  }
  z <- (cl$truth[sel] - cl$mean[sel]) / cl$sd[sel]
  nuz <- .tune_nu_of(z)
  rz <- nuz$robust_z
  nu <- nuz$nu
  if (is.na(nu)) {
    stop(
      "the implied nu off the rank-", rank, " arm is unreadable (fewer ",
      "than 4 finite z); fit through mipca() with a nu of your own.",
      call. = FALSE
    )
  }
  nu_implied <- nu
  nu <- max(nu, .NU_FLOOR)

  cfg <- .tune_run_config(object)
  given <- c(
    chains = !missing(chains),
    iter_warmup = !missing(iter_warmup),
    iter_sampling = !missing(iter_sampling),
    thin = !missing(thin),
    spectra = !missing(spectra),
    parallel_chains = !missing(parallel_chains),
    init = !missing(init),
    colmax = !missing(colmax),
    refresh = !missing(refresh)
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
  dq_old <- .seed_dq(seed)
  if (!is.null(dq_old)) {
    on.exit(dqrng::dqrng_set_state(dq_old), add = TRUE)
  }
  res <- do.call(
    .mipca_run,
    c(
      list(X, rank = rank, transform = object$transform, nu = nu),
      carried,
      dots,
      list(
        keep_state = keep_state, prop_diag = prop_diag,
        prop_diag_given = !missing(prop_diag), call = cl_fit
      )
    ),
    quote = TRUE
  )

  keep <- c(
    "param_set", "rank", "ndraws", "ess_p10", "ess_ok", "rhat_max",
    "robust_z", "implied_nu", "implied_lo", "implied_hi"
  )
  rec <- list(
    rank = rank,
    rank_src = rank_src,
    nu = nu,
    nu_implied = nu_implied,
    robust_z = rz,
    n_cells = sum(sel),
    n_ho_cells = length(object$prov$hold$lin),
    n_ho_cols = length(unique(object$prov$hold$hcol)),
    carried = names(carried),
    overridden = names(given)[given],
    arms = g[g$rank == rank & is.na(g$error), keep, drop = FALSE],
    tune_call = object$call
  )
  rec$line <- .tuned_line(rec)
  if (inherits(res, "si_mi_pca")) {
    res$prov$tuned <- rec
  } else {
    res$chunks[[1L]]$prov$tuned <- rec
  }
  res
}

.tuned_line <- function(tp) {
  lo <- tp$arms$implied_lo
  hi <- tp$arms$implied_hi
  nu_txt <- if (tp$nu > tp$nu_implied) {
    sprintf("%.4g (raised from the tuner's %.3g)", tp$nu, tp$nu_implied)
  } else {
    sprintf("%.4g", tp$nu)
  }
  paste0(
    sprintf(
      "  tuned:   rank %d (%s) and nu %s from a tuner run (%d cells over %d columns, band %.3g-%.3g)\n",
      tp$rank, tp$rank_src, nu_txt, tp$n_cells, tp$n_ho_cols,
      .na_red(lo, min), .na_red(hi, max)
    ),
    sprintf(
      "           carried from the arms: %s%s\n",
      if (length(tp$carried) > 0L) paste(tp$carried, collapse = ", ") else "nothing",
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
  g <- x$grid
  where <- if (is.null(x$data)) {
    "on NO MATRIX -- detached"
  } else {
    sprintf("on a %d x %d matrix", nrow(x$data), ncol(x$data))
  }
  cat(sprintf(
    "si_tune: %d arm%s %s (%s)\n",
    nrow(g), .s(nrow(g)), where, x$transform$name
  ))
  n_ho_cols <- length(unique(x$prov$hold$hcol))
  cat(sprintf(
    "  holdout %d cells over %d column%s, diag %d cells/arm; bars: ESS >= %d, Rhat <= %.2f\n",
    length(x$prov$hold$lin),
    n_ho_cols,
    .s(n_ho_cols),
    x$prov$n_diag_req,
    x$prov$ess_bar, x$prov$rhat_bar
  ))
  for (i in seq_along(x$call)) {
    s <- paste(deparse(x$call[[i]], nlines = 3L), collapse = " ")
    if (nchar(s) > 66L) s <- paste0(substr(s, 1L, 63L), "...")
    cat(sprintf("  %-8s %s\n", if (i == 1L) "call:" else "then:", s))
  }
  cat(sprintf(
    "  %-5s %9s %7s %7s %7s %8s %8s %7s\n",
    "rank", "crps", "ess10", "essOK", "rhatMx",
    "impl.nu", "band", "time.s"
  ))
  blocks <- .tune_blocks(
    x$cells, g, rep(length(x$prov$hold$lin), nrow(g)), "cells"
  )
  for (k in seq_len(nrow(g))) {
    if (!is.na(g$error[k])) {
      cat(sprintf("  %-5d ERROR: %s\n", g$rank[k], g$error[k]))
      next
    }
    cr <- blocks[[k]]$crps
    cat(sprintf(
      "  %-5d %9.5f %7.0f %7.2f %7.3f %8.2f %3.1f-%-4.1f %7.1f\n",
      g$rank[k], mean(cr), g$ess_p10[k], g$ess_ok[k],
      g$rhat_max[k], g$implied_nu[k], g$implied_lo[k], g$implied_hi[k],
      g$time[k]
    ))
  }
  cat("  reference semantics: `b <- x` aliases\n")
  cat("  next:    update(x, iter_sampling) extends every arm; update(x, parameters =) adds ranks\n")
  cat("           mipca_fit(x, rank) hands off\n")
  invisible(x)
}
