# ---- The flank driver: one fit per target column, joined ----------------------
#
# mipca_flank_plan() is the pre-flight, run once (printing it is the dry run);
# mipca_flank(plan, rank, seed =) fits each target on its admissible window,
# stores only the target, and cbind()s the one-column fits.
#
# A list of columns for ABB (`abb`, R/mi_pca_abb.R) is given to the plan and
# to nothing else here: listed columns are out of every window, and a listed
# target is drawn by .abb_piece() with no window and no sampler.

# The targets' windows in the ADMISSIBLE column system, through
# find_windows_flank(). `targets` integer, ascending, unique.
.flank_windows <- function(location, targets, window_size) {
  win <- find_windows_flank(location, targets, window_size)
  start <- win$start
  end <- win$end
  local <- win$subset_local
  ## LEAKAGE-CLASS guard: a wrong `local` stores a neighbour's draws under the
  ## target's name.
  if (!identical(start + local - 1L, targets) || any(end < targets)) {
    stop("internal error: a flank window does not hold its own target.")
  }
  list(start = start, end = end, local = local)
}

# The plan's door; returns what it validated and resolved.
.flank_door <- function(
  obj,
  location,
  window_size,
  subset,
  min_window_n,
  colmax,
  abb = NULL
) {
  .chk_matrix(obj)
  ## slide_imp()'s own assertions on these two.
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
    tg = .resolve_cols(subset, obj), # ascending, unique, integer
    ## The list, stated ONCE for the whole flank journey: here.
    abb = .resolve_cols(abb, obj, what = "abb")
  )
}

# One census for every window, each target's window and status. `a_start` /
# `a_end` / `local` index `adm_idx` (the admissible system); the table's
# `start` / `end` are the same bounds as original column positions.
#
# `abb` (resolved, or NULL) is the list of columns that leave the sampler.
# They join the exclusion below, so every window is built from kept columns
# and the window tuned on is the window fitted; a listed TARGET gets no
# window at all and the status `abb`.
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
  ## A kept column is refused on its variance, as at the sampler's door. A
  ## listed column is in no window, so its variance is nobody's business: a
  ## listed NON-target is never read, and a listed target is read by the
  ## bootstrap, which asks only for finite cells (.abb_split() has the same
  ## predicate). Refused here, at the plan, and not when the drive reaches it.
  bad <- cen$nonfin & !is_abb
  tg_l <- tg[is_abb[tg]]
  if (length(tg_l) > 0L) {
    bad[tg_l] <- colSums(is.infinite(obj[, tg_l, drop = FALSE])) > 0
  }
  if (any(bad)) {
    k <- which(bad)
    stop(
      "`obj` has ",
      length(k),
      " column",
      .s(length(k)),
      " with a non-finite ",
      "variance (an Inf cell, or overflow) [",
      fmt_trunc(.col_label(obj, k), 8),
      "]. Remove or repair these columns before imputing.",
      call. = FALSE
    )
  }
  ## Windowing policy over the sampler's floors; the sampler takes no
  ## `colmax`. Must stay the exact negation of pca_imp's
  ## `miss_rate <= min(colmax, 1)`, off the same scanner.
  too_miss <- !(cen$empty | cen$single | cen$nonfin | cen$const) &
    cen$cmiss / nrow(obj) > colmax
  ## unname(): col_vars()'s names ride into which(), and a named position
  ## fails the identical() guards on the name alone.
  adm_idx <- unname(which(
    !(cen$empty | cen$single | cen$const | too_miss | is_abb)
  ))
  ## Scatter table, original column -> admissible position, 0 = inadmissible.
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

  ## "clean" (no hole) wins over the rest; "inadmissible" and "window_small"
  ## targets are carried with their holes left NA.
  status <- rep("impute", nt)
  status[!has_win] <- "window_small"
  status[tg_adm == 0L] <- "inadmissible"
  ## A listed target is drawn from its own observed cells, so for it
  ## admissibility is a donor question only: a constant column, a one-cell
  ## column and one over `colmax` all draw; a column with no observed cell
  ## has nothing to draw from.
  tg_abb <- is_abb[tg]
  status[tg_abb] <- "abb"
  status[tg_abb & cen$empty[tg]] <- "inadmissible"
  status[holes == 0L] <- "clean"

  ## slide_imp()'s refusal, widened to inadmissible targets; a request of clean
  ## targets only is not refused.
  n_small <- sum(status == "window_small")
  n_inadm <- sum(status == "inadmissible")
  if (
    (n_small + n_inadm) > 0L &&
      !any(status == "impute" | status == "abb")
  ) {
    stop(
      "No requested target with a missing cell can be imputed.",
      if (n_small > 0L) {
        paste0(
          "\n",
          n_small,
          " target",
          if (n_small == 1L) " has" else "s have",
          " a window with fewer than `min_window_n` (",
          min_window_n,
          ") columns. Increase `window_size`, or lower `min_window_n`."
        )
      },
      if (n_inadm > 0L) {
        paste0(
          "\n",
          n_inadm,
          " target",
          if (n_inadm == 1L) " is" else "s are",
          " inadmissible: constant, with fewer than two observed cells, or ",
          "missing more than `colmax` of ",
          if (n_inadm == 1L) "its" else "their",
          " cells."
        )
      },
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

# THE PLAN: slide_imp()'s pre-flight, run once and kept; fits nothing, draws
# no RNG. A classed list: nothing mutates it. No `rank`: the plan is what the
# rank is tuned on, so min_window_n >= rank + 2 is the driver's check.
#' @export
mipca_flank_plan <- function(
  obj,
  location,
  window_size,
  subset,
  min_window_n,
  colmax = 0.95,
  ## The columns that leave the sampler for ABB, as at mipca(). On the plan
  ## and nowhere else in the flank journey: the window, the tuner run on it
  ## and the driver all read it from here.
  abb = NULL
) {
  door <- .flank_door(
    obj,
    location,
    window_size,
    subset,
    min_window_n,
    colmax,
    abb
  )
  cen <- .flank_census(
    obj,
    location,
    window_size,
    door$tg,
    door$min_window_n,
    door$colmax,
    door$abb
  )
  structure(
    list(
      targets = cen$targets,
      obj = obj,
      window_size = window_size,
      min_window_n = door$min_window_n,
      colmax = door$colmax,
      abb = door$abb,
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

# Class and per-target lengths only; not a tamper check. The guard that
# matters is per window, in .flank_window_cols().
.chk_flank_plan <- function(plan) {
  if (!inherits(plan, "si_flank_plan")) {
    stop(
      "`plan` must be the object mipca_flank_plan() returns.",
      call. = FALSE
    )
  }
  nt <- nrow(plan$targets)
  per_target <- c("tg_adm", "tg_abb", "a_start", "a_end", "local")
  if (!all(lengths(plan[per_target]) == nt)) {
    stop("internal error: the plan's windows do not line up with its targets.")
  }
  invisible(plan)
}

# Target i's window as original columns, guarded. The ONE cut, shared by the
# driver and mipca_flank_window(), so the window tuned on is the one fitted.
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
  cat(sprintf(
    "si_flank_plan: %d target%s on a %d x %d matrix, nothing fitted\n",
    nt,
    .s(nt),
    nrow(x$obj),
    ncol(x$obj)
  ))
  ## The count per `status` value, spelled as the table spells it, for the
  ## values present: the preview below shows only the first rows.
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

# One target's window as a matrix, the tuner's way in. A clean target has a
# window.
#' @export
mipca_flank_window <- function(plan, subset) {
  .chk_flank_plan(plan)
  if (missing(subset) || is.null(subset)) {
    stop(
      "`subset` must name exactly one target column of the plan.",
      call. = FALSE
    )
  }
  tg <- .resolve_cols(subset, plan$obj)
  if (length(tg) != 1L) {
    stop(
      "`subset` must name exactly one target column of the plan.",
      call. = FALSE
    )
  }
  ## Scatter table: original column -> its plan row, 0 = not a target.
  inv <- integer(ncol(plan$obj))
  inv[plan$targets$target] <- seq_len(nrow(plan$targets))
  i <- inv[tg]
  if (i == 0L) {
    stop(
      "Column ",
      colnames(plan$obj)[tg],
      " is not a target of this plan. ",
      "Name a column from `plan$targets$name`, or build a plan for this ",
      "column with mipca_flank_plan().",
      call. = FALSE
    )
  }
  ## Before the window test: a listed target has none by design, clean or
  ## not, and "inadmissible" would be the wrong reason.
  if (plan$tg_abb[i]) {
    stop(
      "Target ",
      plan$targets$name[i],
      " is listed in `abb`, so it has no ",
      "window. Its missing cells are drawn from its own observed cells, and ",
      "no model is fitted for it.",
      call. = FALSE
    )
  }
  if (is.na(plan$a_start[i])) {
    stop(
      "Target ",
      plan$targets$name[i],
      " has no window",
      if (plan$tg_adm[i] == 0L) {
        ", because the column is inadmissible."
      } else {
        paste0(
          ": fewer than `min_window_n` (",
          plan$min_window_n,
          ") admissible columns lie within `window_size` of it. Rebuild the ",
          "plan with a larger `window_size` or a smaller `min_window_n`."
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
  ## REQUIRED; the refusal below says why.
  seed,
  transform,
  prop_diag = .prop_diag_default(),
  ## mipca()'s defaults, repeated: keep in step with mipca(). No `refresh`
  ## (`.progress` reports windows) and no `colmax` (the plan's).
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
    stop("`.progress` must be TRUE or FALSE.", call. = FALSE)
  }
  .chk_flank_plan(plan)
  obj <- plan$obj
  targets <- plan$targets
  tg <- targets$target
  min_window_n <- plan$min_window_n
  n <- nrow(obj)

  rank <- .chk_count(rank, "rank")
  ## The Gram bound at the door (Spectra's own: rank + 2 <= min(n, p)), so
  ## every kept window is one the sampler accepts; min_window_n stands in
  ## for p.
  if (min_window_n < rank + 2L) {
    stop(
      "The plan's `min_window_n` (",
      min_window_n,
      ") must be at least ",
      "`rank + 2` = ",
      rank + 2L,
      " at rank ",
      rank,
      ", because a narrower ",
      "window cannot be fitted. Rebuild the plan with mipca_flank_plan() and ",
      "a `min_window_n` of at least ",
      rank + 2L,
      ", or lower `rank`.",
      call. = FALSE
    )
  }
  rank <- .chk_rank(rank, n, min_window_n)
  if (missing(seed)) {
    stop(
      "`seed` must be given: each target draws from a generator seeded by ",
      "`seed` and the target's column.",
      call. = FALSE
    )
  }
  seed <- .chk_count(seed, "`seed`", lo = 0L)
  chains <- .chk_count(chains, "chains")
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  ## At the DOOR with the core's predicates: refused once, not per window.
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  ## Not .chk_count's floor: zero warmup is legal at the core.
  iter_warmup <- .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  thin <- as.double(.chk_count(thin, "thin"))
  nu <- .chk_nu(nu)
  spectra <- .chk_spectra(spectra)
  init <- .chk_init(init)
  ## Warned ONCE here; the loop muffles the core's per-window copies.
  if (nu < .NU_FLOOR) {
    warning(.nu_floor_msg(), call. = FALSE)
  }
  ## Resolved once: every window shares the row count.
  transform <- .resolve_transform(transform, n)

  ## The plan's print gives the counts; `.progress` reports windows only.
  status <- targets$status
  nt <- length(tg)
  nwin <- sum(status == "impute")

  cl <- .scrub_call(match.call(), "mipca_flank")
  ## The loop reseeds per window; the caller's generator goes back as found.
  dq <- dqrng::dqrng_get_state()
  on.exit(dqrng::dqrng_set_state(dq), add = TRUE)

  ## One piece per requested target, in request order.
  fits <- vector("list", nt)
  step <- max(1L, nwin %/% 20L)
  kwin <- 0L
  for (i in seq_len(nt)) {
    if (status[i] == "abb") {
      ## A listed target: no window, no sampler. It draws off
      ## dqset.seed(c(seed, its column index)), as at the plain door, so it
      ## does not move with a neighbour's status and costs the same at any
      ## request position.
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
      ## Not sampled: a zero-hole piece, or its holes left NA.
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
      message(sprintf("Processing window %d of %d", kwin, nwin))
    }
    ## Off c(seed, the target's column index), the key a listed target uses,
    ## with the chains on streams 0 .. chains - 1 of it as at the plain door.
    ## So a window's draws move with no other target, requested or not, and
    ## seeding costs the same at any request position.
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
        ## Re-thrown, never skipped: a skipped window would hide the error
        ## behind a partly imputed result.
        error = function(e) {
          stop(
            "Target ",
            targets$name[i],
            " (window of ",
            length(wcols),
            " columns): ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      ),
      ## Matched on .nu_floor_msg()'s exact text, so no other warning is
      ## touched.
      warning = function(w) {
        if (identical(conditionMessage(w), .nu_floor_msg())) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }

  fit <- do.call(cbind, fits)
  ## cbind() records its own call, which under do.call is placeholders.
  fit$call <- cl
  ## Joined column k must be request k, so row k of plan$targets describes
  ## it. The table stays on the plan, not the fit: a join would stack tables
  ## whose positions belong to different matrices.
  if (!identical(fit$colnames, targets$name)) {
    stop("internal error: the joined fit does not hold its targets in order.")
  }
  fit
}
