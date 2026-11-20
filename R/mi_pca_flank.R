.flank_windows <- function(location, targets, window_size) {
  win <- find_windows_flank(location, targets, window_size)
  start <- win$start
  end <- win$end
  local <- win$subset_local
  if (!identical(start + local - 1L, targets) || any(end < targets)) {
    stop("internal error: a flank window does not hold its own target")
  }
  list(start = start, end = end, local = local)
}

.flank_census <- function(
  obj,
  location,
  window_size,
  tg,
  min_window_n,
  colmax,
  abb = NULL
) {
  p <- ncol(obj)
  is_abb <- if (is.null(abb)) logical(p) else .abb_kept(p, abb)$is_abb
  cen <- .admissible_census(obj)
  bad <- cen$nonfin & !is_abb
  tg_l <- tg[is_abb[tg]]
  if (length(tg_l) > 0L) {
    bad[tg_l] <- colSums(is.infinite(obj[, tg_l, drop = FALSE])) > 0
  }
  if (any(bad)) {
    k <- which(bad)
    stop(
      "'obj' has ",
      length(k),
      " column",
      .s(length(k)),
      " with a non-finite ",
      "variance (an Inf cell, or overflow) [",
      fmt_trunc(.col_label(obj, k), 8),
      "]. Remove or repair these columns before imputing",
      call. = FALSE
    )
  }
  too_miss <- cen$cmiss / nrow(obj) > colmax
  adm_idx <- unname(which(
    !(cen$empty | cen$single | cen$const | too_miss | is_abb)
  ))
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
  tg_abb <- is_abb[tg]
  status[tg_abb] <- "abb"
  status[tg_abb & cen$empty[tg]] <- "inadmissible"
  status[holes == 0L] <- "clean"

  n_small <- sum(status == "window_small")
  n_inadm <- sum(status == "inadmissible")
  if (
    (n_small + n_inadm) > 0L &&
      !any(status == "impute" | status == "abb")
  ) {
    small <- if (n_small > 0L) {
      paste0(
        "\n",
        n_small,
        " target",
        if (n_small == 1L) " has" else "s have",
        " a window with fewer than 'min_window_n' (",
        min_window_n,
        ") columns. Increase 'window_size', or lower 'min_window_n'"
      )
    }
    inadm <- if (n_inadm > 0L) {
      paste0(
        "\n",
        n_inadm,
        " target",
        if (n_inadm == 1L) " is" else "s are",
        " inadmissible: constant, with fewer than two observed cells, or ",
        "missing more than 'colmax' of ",
        if (n_inadm == 1L) "its" else "their",
        " cells"
      )
    }
    stop(
      "no requested target with a missing cell can be imputed",
      small,
      inadm,
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
    tg_abb = unname(tg_abb),
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
  colmax = 0.95,
  abb = NULL
) {
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
  min_window_n <- .chk_count(min_window_n, "min_window_n", lo = 2L)
  colmax <- .chk_colmax(colmax)
  if (missing(subset) || is.null(subset)) {
    stop(
      "'subset' must name the target columns: a flank window is built around ",
      "each one, so there is no \"every column\" default",
      call. = FALSE
    )
  }
  tg <- .resolve_cols(subset, obj)
  abb <- .resolve_cols(abb, obj, what = "abb")
  cen <- .flank_census(
    obj,
    location,
    window_size,
    tg,
    min_window_n,
    colmax,
    abb
  )
  structure(
    list(
      targets = cen$targets,
      obj = obj,
      window_size = window_size,
      min_window_n = min_window_n,
      colmax = colmax,
      abb = abb,
      adm_idx = cen$adm_idx,
      tg_adm = cen$tg_adm,
      tg_abb = cen$tg_abb,
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
      "'plan' must be the object mipca_flank_plan() returns",
      call. = FALSE
    )
  }
  nt <- nrow(plan$targets)
  per_target <- c("tg_adm", "tg_abb", "a_start", "a_end", "local")
  if (!all(lengths(plan[per_target]) == nt)) {
    stop("internal error: the plan's windows do not line up with its targets")
  }
  invisible(plan)
}

.flank_window_cols <- function(plan, i) {
  wcols <- plan$adm_idx[plan$a_start[i]:plan$a_end[i]]
  if (!identical(wcols[plan$local[i]], plan$targets$target[i])) {
    stop("internal error: a flank window's local index misses its target")
  }
  wcols
}

#' @export
print.si_flank_plan <- function(x, n = 10L, ...) {
  .chk_no_dots(..., .what = "print()")
  tg <- x$targets
  nt <- nrow(tg)
  cat(sprintf(
    "si_flank_plan: %d target%s on a %d x %d matrix, nothing fitted\n",
    nt,
    .s(nt),
    nrow(x$obj),
    ncol(x$obj)
  ))
  lv <- c("impute", "clean", "abb", "inadmissible", "window_small")
  cnt <- vapply(lv, function(s) sum(tg$status == s), 0L)
  cat(
    "  status: ",
    paste(lv[cnt > 0L], cnt[cnt > 0L], collapse = ", "),
    "\n",
    sep = ""
  )
  k <- min(nt, n)
  print(tg[seq_len(k), , drop = FALSE], row.names = FALSE)
  if (nt > k) {
    cat(sprintf("  ... %d more row%s in x$targets\n", nt - k, .s(nt - k)))
  }
  cat("  next: mipca_flank(x, rank, seed, transform) imputes.\n")
  cat(
    "        mipca_flank_window(x, subset) gives one target's window to tune on.\n"
  )
  invisible(x)
}

#' @export
mipca_flank_window <- function(plan, subset) {
  .chk_flank_plan(plan)
  tg <- if (missing(subset)) NULL else .resolve_cols(subset, plan$obj)
  if (length(tg) != 1L) {
    stop(
      "'subset' must name exactly one target column of the plan",
      call. = FALSE
    )
  }
  inv <- integer(ncol(plan$obj))
  inv[plan$targets$target] <- seq_len(nrow(plan$targets))
  i <- inv[tg]
  if (i == 0L) {
    stop(
      "column ",
      colnames(plan$obj)[tg],
      " is not a target of this plan. ",
      "Name a column from 'plan$targets$name', or build a plan for this ",
      "column with mipca_flank_plan()",
      call. = FALSE
    )
  }
  if (plan$tg_abb[i]) {
    stop(
      "target ",
      plan$targets$name[i],
      " is listed in 'abb', so it has no ",
      "window. Its missing cells are drawn from its own observed cells, and ",
      "no model is fitted for it",
      call. = FALSE
    )
  }
  if (is.na(plan$a_start[i])) {
    if (plan$tg_adm[i] == 0L) {
      stop(
        "target ",
        plan$targets$name[i],
        " has no window, because the ",
        "column is inadmissible",
        call. = FALSE
      )
    }
    stop(
      "target ",
      plan$targets$name[i],
      " has no window: fewer than ",
      "'min_window_n' (",
      plan$min_window_n,
      ") admissible columns lie ",
      "within 'window_size' of it. Rebuild the plan with a larger ",
      "'window_size' or a smaller 'min_window_n'",
      call. = FALSE
    )
  }
  plan$obj[, .flank_window_cols(plan, i), drop = FALSE]
}

#' @export
mipca_flank <- function(
  plan,
  rank,
  seed,
  transform,
  prop_diag = 0.5,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  spectra = spectra_opts(),
  init = 1,
  .progress = TRUE
) {
  if (!isTRUE(.progress) && !isFALSE(.progress)) {
    stop("'.progress' must be TRUE or FALSE", call. = FALSE)
  }
  .chk_flank_plan(plan)
  obj <- plan$obj
  targets <- plan$targets
  tg <- targets$target
  min_window_n <- plan$min_window_n
  n <- nrow(obj)

  rank <- .chk_count(rank, "rank")
  if (min_window_n < rank + 2L) {
    stop(
      "the plan's 'min_window_n' (",
      min_window_n,
      ") must be at least ",
      "'rank + 2' = ",
      rank + 2L,
      " at rank ",
      rank,
      ", because a narrower ",
      "window cannot be fitted. Rebuild the plan with mipca_flank_plan() and ",
      "a 'min_window_n' of at least ",
      rank + 2L,
      ", or lower 'rank'",
      call. = FALSE
    )
  }
  rank <- .chk_rank(rank, n, min_window_n)
  if (missing(seed)) {
    stop(
      "'seed' must be given: each target draws from a generator seeded by ",
      "'seed' and the target's column",
      call. = FALSE
    )
  }
  seed <- .chk_count(seed, "seed", lo = 0L)
  chains <- .chk_count(chains, "chains")
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  iter_warmup <- .chk_count(iter_warmup, "iter_warmup", lo = 0L)
  thin <- as.double(.chk_count(thin, "thin"))
  nu <- .chk_nu(nu)
  spectra <- .chk_spectra(spectra)
  init <- .chk_init(init)
  if (nu < .NU_FLOOR) {
    warning(.nu_floor_msg(), call. = FALSE)
  }
  transform <- .resolve_transform(transform, n)

  status <- targets$status
  nt <- length(tg)
  nwin <- sum(status == "impute")
  if (any(status == "abb")) {
    .chk_abb_draws(
      n,
      chains,
      iter_sampling,
      fix = "Lower 'iter_sampling' or 'chains'"
    )
  }

  cl <- .scrub_call(match.call(), "mipca_flank")
  dq <- dqrng::dqrng_get_state()
  on.exit(dqrng::dqrng_set_state(dq), add = TRUE)

  fits <- vector("list", nt)
  step <- max(1L, nwin %/% 20L)
  kwin <- 0L
  for (i in seq_len(nt)) {
    if (status[i] == "abb") {
      fits[[i]] <- .abb_piece(
        obj,
        cols = tg[i],
        seed = seed,
        key = tg[i],
        chains = chains,
        ndraws = iter_sampling,
        thin = thin,
        prop_diag = prop_diag,
        call = cl
      )
      next
    }
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
      message(sprintf("processing window %d of %d", kwin, nwin))
    }
    dqrng::dqset.seed(c(seed, tg[i]))
    fits[[i]] <- withCallingHandlers(
      tryCatch(
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
          spectra = spectra,
          init = init,
          call = cl
        ),
        error = function(e) {
          stop(
            "target ",
            targets$name[i],
            " (window of ",
            length(wcols),
            " columns): ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      ),
      warning = function(w) {
        if (identical(conditionMessage(w), .nu_floor_msg())) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  fit <- .si_mi_join(fits, call = cl)
  if (!identical(fit$colnames, targets$name)) {
    stop("internal error: the joined fit does not hold its targets in order")
  }
  fit
}
