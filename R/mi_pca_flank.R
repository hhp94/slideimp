.flank_windows <- function(location, targets, window_size) {
  win <- find_windows_flank(location, targets, window_size)
  start <- win$start
  end <- win$end
  local <- win$subset_local
  if (!identical(start + local - 1L, targets) || any(end < targets)) {
    stop("internal error: a flank window does not hold its own target.")
  }
  list(start = start, end = end, local = local)
}

.flank_door <- function(obj, location, window_size, subset, min_window_n,
                        colmax) {
  .chk_matrix(obj)
  checkmate::assert_numeric(
    location,
    len = ncol(obj),
    any.missing = FALSE,
    sorted = TRUE,
    finite = TRUE,
    .var.name = "location"
  )
  checkmate::assert_number(
    window_size,
    lower = .Machine$double.eps,
    finite = TRUE,
    .var.name = "window_size"
  )
  min_window_n <- .chk_count(min_window_n, "`min_window_n`", lo = 2L)
  colmax <- .chk_colmax(colmax)
  if (missing(subset) || is.null(subset)) {
    stop(
      "`subset` must name the target columns: a flank window is built around ",
      "each one, so there is no \"every column\" default.",
      call. = FALSE
    )
  }
  list(
    min_window_n = min_window_n,
    colmax = colmax,
    tg = .resolve_cols(subset, obj)
  )
}

.flank_census <- function(obj, location, window_size, tg, min_window_n,
                          colmax) {
  p <- ncol(obj)
  cen <- .admissible_census(obj, colmax)
  if (any(cen$nonfin)) {
    k <- which(cen$nonfin)
    stop(
      "`obj` has ", length(k), " column", .s(length(k)), " with a non-finite ",
      "variance (an Inf cell, or overflow) [", fmt_trunc(k, 8), "]; remove or ",
      "repair them before imputing.",
      call. = FALSE
    )
  }
  adm_idx <- unname(which(!(cen$empty | cen$single | cen$const | cen$too_miss)))
  inv <- integer(p)
  inv[adm_idx] <- seq_along(adm_idx)
  tg_adm <- inv[tg]
  holes <- cen$cmiss[tg]

  nt <- length(tg)
  start <- rep(NA_integer_, nt)
  end <- rep(NA_integer_, nt)
  window_n <- rep(NA_integer_, nt)
  a_start <- rep(NA_integer_, nt)
  a_end <- rep(NA_integer_, nt)
  local <- rep(NA_integer_, nt)

  q <- which(tg_adm > 0L)
  has_win <- logical(nt)
  if (length(q) > 0L) {
    win <- .flank_windows(location[adm_idx], tg_adm[q], window_size)
    wn <- win$end - win$start + 1L
    big <- wn >= min_window_n
    kept <- q[big]
    has_win[kept] <- TRUE
    a_start[kept] <- win$start[big]
    a_end[kept] <- win$end[big]
    local[kept] <- win$local[big]
    start[kept] <- adm_idx[win$start[big]]
    end[kept] <- adm_idx[win$end[big]]
    window_n[kept] <- wn[big]
  }

  status <- rep("impute", nt)
  status[!has_win] <- "window_small"
  status[tg_adm == 0L] <- "inadmissible"
  status[holes == 0L] <- "clean"

  if (any(status == "window_small") && !any(status == "impute")) {
    stop(
      "All windows have fewer than `min_window_n` (", min_window_n,
      ") columns. Consider increasing `window_size` or decreasing ",
      "`min_window_n`.",
      call. = FALSE
    )
  }

  list(
    targets = data.frame(
      target = tg,
      name = colnames(obj)[tg],
      holes = holes,
      start = start,
      end = end,
      window_n = window_n,
      status = status,
      stringsAsFactors = FALSE
    ),
    adm_idx = adm_idx,
    tg_adm = tg_adm,
    a_start = a_start,
    a_end = a_end,
    local = local
  )
}

#' @export
mipca_flank_plan <- function(
  obj,
  location,
  window_size,
  subset,
  min_window_n,
  colmax = 0.95
) {
  door <- .flank_door(obj, location, window_size, subset, min_window_n, colmax)
  cen <- .flank_census(
    obj, location, window_size, door$tg, door$min_window_n, door$colmax
  )
  structure(
    list(
      targets = cen$targets,
      obj = obj,
      window_size = window_size,
      min_window_n = door$min_window_n,
      colmax = door$colmax,
      adm_idx = cen$adm_idx,
      tg_adm = cen$tg_adm,
      a_start = cen$a_start,
      a_end = cen$a_end,
      local = cen$local,
      call = .scrub_call(match.call(), "mipca_flank_plan")
    ),
    class = "si_flank_plan"
  )
}

.chk_flank_plan <- function(plan) {
  if (!inherits(plan, "si_flank_plan")) {
    stop(
      "`plan` must be what mipca_flank_plan() returns: the matrix, its ",
      "targets and their windows, checked once.",
      call. = FALSE
    )
  }
  nt <- nrow(plan$targets)
  per_target <- c("tg_adm", "a_start", "a_end", "local")
  if (!all(lengths(plan[per_target]) == nt)) {
    stop("internal error: the plan's windows do not line up with its targets.")
  }
  invisible(plan)
}

.flank_window_cols <- function(plan, i) {
  wcols <- plan$adm_idx[plan$a_start[i]:plan$a_end[i]]
  if (!identical(wcols[plan$local[i]], plan$targets$target[i])) {
    stop("internal error: a flank window's local index misses its target.")
  }
  wcols
}

#' @export
print.si_flank_plan <- function(x, n = 10L, ...) {
  .chk_no_dots(..., .what = "print()")
  tg <- x$targets
  nt <- nrow(tg)
  st <- tg$status
  cat(sprintf(
    "si_flank_plan: %d target%s on a %d x %d matrix, nothing fitted\n",
    nt, .s(nt), nrow(x$obj), ncol(x$obj)
  ))
  cat(sprintf(
    "  windows: within %s of each target, at least %d admissible columns (colmax %s)\n",
    format(x$window_size), x$min_window_n, format(x$colmax)
  ))
  cat(sprintf(
    "  status:  %d to impute, %d clean, %d inadmissible, %d in windows under %d columns\n",
    sum(st == "impute"), sum(st == "clean"), sum(st == "inadmissible"),
    sum(st == "window_small"), x$min_window_n
  ))
  k <- min(nt, n)
  print(tg[seq_len(k), , drop = FALSE], row.names = FALSE)
  if (nt > k) {
    cat(sprintf("  ... %d more; x$targets is the whole table\n", nt - k))
  }
  cat("  next:    mipca_flank_window(x, subset) is one target's window, for mipca_tune()\n")
  cat("           mipca_flank(x, rank, seed = ) imputes\n")
  invisible(x)
}

#' @export
mipca_flank_window <- function(plan, subset) {
  .chk_flank_plan(plan)
  if (missing(subset) || is.null(subset)) {
    stop(
      "`subset` must name exactly one target column: the window is that ",
      "target's.",
      call. = FALSE
    )
  }
  tg <- .resolve_cols(subset, plan$obj)
  if (length(tg) != 1L) {
    stop(
      "`subset` must name exactly one target column: the window is that ",
      "target's.",
      call. = FALSE
    )
  }
  inv <- integer(ncol(plan$obj))
  inv[plan$targets$target] <- seq_len(nrow(plan$targets))
  i <- inv[tg]
  if (i == 0L) {
    stop(
      "column ", colnames(plan$obj)[tg], " is not a target of this plan; a ",
      "window exists only for the columns the plan was built around.",
      call. = FALSE
    )
  }
  if (is.na(plan$a_start[i])) {
    stop(
      "target ", plan$targets$name[i], " has no window to hand back: ",
      if (plan$tg_adm[i] == 0L) {
        "the column is inadmissible."
      } else {
        paste0(
          "fewer than `min_window_n` (", plan$min_window_n,
          ") admissible columns lie within `window_size` of it."
        )
      },
      call. = FALSE
    )
  }
  plan$obj[, .flank_window_cols(plan, i), drop = FALSE]
}

#' @export
mipca_flank <- function(
  plan,
  rank,
  seed = NULL,
  transform,
  prop_diag = .prop_diag_default(),
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  parallel_chains = 1L,
  spectra = spectra_opts(),
  init = 1,
  .progress = TRUE
) {
  if (!isTRUE(.progress) && !isFALSE(.progress)) {
    stop("`.progress` must be TRUE or FALSE.", call. = FALSE)
  }
  .chk_flank_plan(plan)
  obj <- plan$obj
  targets <- plan$targets
  tg <- targets$target
  min_window_n <- plan$min_window_n
  n <- nrow(obj)

  rank <- .chk_count(rank, "rank")
  need <- rank + 2L
  if (need > n) {
    stop(
      "rank must be at most ", n - 2L, " with ", n, " rows: the sampler ",
      "solves rank + 1 eigenpairs, and its eigensolver needs that many below ",
      "the smaller dimension.",
      call. = FALSE
    )
  }
  if (min_window_n < need) {
    stop(
      "the plan's `min_window_n` (", min_window_n, ") must be at least ",
      "rank + 2 = ", need, " at rank ", rank, ": a narrower window ",
      "cannot be fitted, and a window that is kept must be one the sampler ",
      "accepts. Rebuild the plan with a larger `min_window_n`.",
      call. = FALSE
    )
  }
  seed <- .chk_count(seed, "`seed`", lo = 0L)
  chains <- .chk_count(chains, "chains")
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  transform <- .resolve_transform(transform, n)

  status <- targets$status
  nt <- length(tg)
  nwin <- sum(status == "impute")
  if (.progress) {
    message(sprintf(
      "mipca_flank: %d of %d target%s to impute (%d clean, %d inadmissible, %d in windows under %d columns)",
      nwin, nt, .s(nt), sum(status == "clean"),
      sum(status == "inadmissible"), sum(status == "window_small"),
      min_window_n
    ))
  }

  cl <- .scrub_call(match.call(), "mipca_flank")
  dq <- dqrng::dqrng_get_state()
  on.exit(dqrng::dqrng_set_state(dq), add = TRUE)

  fits <- vector("list", nt)
  step <- max(1L, nwin %/% 20L)
  kwin <- 0L
  for (i in seq_len(nt)) {
    if (status[i] != "impute") {
      fits[[i]] <- .mipca_empty(
        obj,
        subset = tg[i],
        rank = rank,
        transform = transform,
        prop_diag = prop_diag,
        iter_sampling = iter_sampling,
        iter_warmup = iter_warmup,
        thin = thin,
        nu = nu,
        chains = chains,
        call = cl,
        leave_na = status[i] != "clean"
      )
      next
    }
    kwin <- kwin + 1L
    wcols <- .flank_window_cols(plan, i)
    if (.progress && (kwin == 1L || kwin %% step == 0L || kwin == nwin)) {
      message(sprintf(" window %d of %d", kwin, nwin))
    }
    dqrng::dqset.seed(seed, stream = (i - 1) * chains)
    fits[[i]] <- tryCatch(
      .mipca_run(
        obj[, wcols, drop = FALSE],
        rank = rank,
        transform = transform,
        subset = plan$local[i],
        keep_state = FALSE,
        prop_diag = prop_diag,
        iter_sampling = iter_sampling,
        iter_warmup = iter_warmup,
        thin = thin,
        nu = nu,
        chains = chains,
        parallel_chains = parallel_chains,
        spectra = spectra,
        colmax = plan$colmax,
        init = init,
        call = cl
      ),
      error = function(e) {
        stop(
          "target ", targets$name[i], " (window of ", length(wcols),
          " columns): ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }

  fit <- do.call(cbind, fits)
  fit$call <- cl
  if (!identical(fit$colnames, targets$name)) {
    stop("internal error: the joined fit does not hold its targets in order.")
  }
  fit$targets <- targets
  fit
}
