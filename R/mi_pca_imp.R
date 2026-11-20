.GRAM_MAX_COND <- 1e-8 / .Machine$double.eps

.NU_FLOOR <- 5

#' @export
spectra_opts <- function(tol = 1e-10, ncv = NULL, maxitr = 1000) {
  tol <- as.double(tol)
  if (length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("'tol' must be a single positive number", call. = FALSE)
  }
  maxitr <- .chk_count(maxitr, "maxitr")
  if (is.null(ncv)) {
    ncv <- 0L
  } else {
    ncv <- .chk_count(ncv, "ncv", note = "or NULL")
  }
  structure(
    list(tol = tol, ncv = ncv, maxitr = maxitr),
    class = "spectra_opts"
  )
}

.SEAM_VERSION <- 10L

.seam_check_version <- function(.state) {
  if (
    !inherits(.state, "si_mi_pca_seam") ||
      !identical(.state$seam_version, .SEAM_VERSION)
  ) {
    stop(
      "update() cannot extend this fit: its saved sampler state is not in ",
      "the format this version reads. Fit again with mipca() and ",
      "'keep_state = TRUE'",
      call. = FALSE
    )
  }
  invisible(.state)
}

.resume_refused <- function(...) {
  paste0(
    "update() cannot extend this fit: the matrix it carries no longer ",
    "matches the one its chains ran on. ",
    ...,
    " Fit again with mipca() ",
    "and 'keep_state = TRUE'"
  )
}

.chk_store <- function(store, nmiss) {
  s <- as.integer(store)
  if (length(s) < 1L || anyNA(s)) {
    stop("'store' must be a non-empty integer vector with no NA", call. = FALSE)
  }
  if (min(s) < 1L || max(s) > nmiss) {
    stop(
      "'store' positions must lie in 1..",
      nmiss,
      " (the matrix has ",
      nmiss,
      " holes)",
      call. = FALSE
    )
  }
  if (length(s) > 1L && min(diff(s)) <= 0L) {
    stop(
      "'store' must be strictly increasing (sorted, no duplicates)",
      call. = FALSE
    )
  }
  s
}

.core_build <- function(Z, cnt, .state, forward) {
  forward_used <- if (is.null(.state)) forward else .state$forward
  built <- .mipca_build_stream(Z, cnt, forward_used, resume = !is.null(.state))
  out <- list(
    Xd = built$Xd,
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
  df_obs <- as.double(n) * p - nmiss - p
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

.start_resume <- function(.state, Xd, etd, dirty, miss_idx_l, n, Z) {
  list(
    xd_list = lapply(seq_along(.state$chain), function(ch) {
      st_ch <- .state$chain[[ch]]
      Xc <- Xd + rep(etd - st_ch$etd, each = n)
      Xc[miss_idx_l] <- st_ch$holes
      .seam_fp_check(st_ch$fp, Xc, dirty, Z)
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

.seam_check_args <- function(.state, chains, chains_given, forward, store) {
  .seam_check_version(.state)
  if (chains_given && as.integer(chains) != .state$chains) {
    stop(
      "internal error: 'chains' disagrees with the seam; a resume continues ",
      "every chain of the fit"
    )
  }
  if (!identical(forward, identity)) {
    stop(
      "internal error: 'forward' cannot be given on a resume: the seam ",
      "carries the map the fit was built under, and the rebuild applies ",
      "that one"
    )
  }
  if (!is.function(.state$forward)) {
    stop("internal error: the seam carries no transform map")
  }
  if (!is.null(store)) {
    stop(
      "internal error: 'store' cannot be given on a resume: the seam ",
      "carries the store set the fit was built under, and an extension ",
      "keeps those same cells"
    )
  }
  .state$chains
}

.seam_check_match <- function(.state, n, cnt, rank, nu, spectra, Z) {
  if (.state$n != n || !identical(.state$cnt, cnt)) {
    stop(
      .resume_refused(
        "Its row count or its count of missing cells per column has changed."
      ),
      call. = FALSE
    )
  }
  if (.state$rank != rank) {
    stop("internal error: the seam was written at a different rank")
  }
  if (.state$nu != nu) {
    stop("internal error: the seam was written with a different nu")
  }
  if (!identical(.state$spectra, spectra)) {
    stop("internal error: the seam was written with different spectra settings")
  }
  .seam_fp_clean_check(.state, Z, n, cnt)
  invisible(.state)
}

.seam_fp_cols <- function(nd) {
  m <- min(nd, max(100L, as.integer(ceiling(0.05 * nd))))
  unique(as.integer(round(seq.int(1L, nd, length.out = m))))
}

.seam_fp_take <- function(Xd, loc) {
  list(
    col = loc,
    sumabs = colSums(abs(Xd[, loc, drop = FALSE]))
  )
}

.SEAM_FP_TOL <- 1e-8

.seam_fp_check <- function(fp, Xd, dirty, Z) {
  cur <- .seam_fp_take(Xd, fp$col)
  if (!all(is.finite(fp$sumabs)) || any(fp$sumabs <= 0)) {
    stop("internal error: the seam's fingerprint is not positive and finite")
  }
  rel <- abs(cur$sumabs - fp$sumabs) / fp$sumabs
  if (max(rel) >= .SEAM_FP_TOL) {
    worst <- which.max(rel)
    stop(
      .resume_refused(
        "Its values, or what its transform makes of them, have changed. ",
        "Column ",
        .col_label(Z, dirty[fp$col[worst]]),
        " differs most."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.seam_fp_clean_take <- function(Z, clean) {
  sa <- colSums(abs(Z[, clean, drop = FALSE]))
  if (!all(is.finite(sa))) {
    stop("internal error: a clean column's fingerprint is not finite")
  }
  list(col = clean, sumabs = sa)
}

.seam_fp_clean_check <- function(.state, Z, n, cnt) {
  split <- .split_cols(n, length(cnt), cnt)$split
  if (!split) {
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
      "was built from a split matrix"
    )
  }
  if (max(fp$col) > length(cnt) || any(cnt[fp$col] != 0L)) {
    stop(
      "internal error: the clean fingerprint names a column that is not ",
      "clean in this census"
    )
  }
  cur <- colSums(abs(Z[, fp$col, drop = FALSE]))
  rel <- abs(cur - fp$sumabs) / pmax(fp$sumabs, .Machine$double.xmin)
  if (!all(is.finite(rel)) || max(rel) >= .SEAM_FP_TOL) {
    worst <- which.max(rel)
    stop(
      .resume_refused(
        "Its values have changed. Column ",
        .col_label(Z, fp$col[worst]),
        " differs most."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.split_cols <- function(n, p, cnt) {
  split <- n <= p && any(cnt == 0L)
  list(split = split, dirty = if (split) which(cnt > 0L) else seq_len(p))
}

.mipca_build_stream <- function(X, cnt, forward, resume) {
  n <- nrow(X)
  p <- ncol(X)
  block <- 1024L
  in_place <- identical(forward, identity)
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

  Xd <- matrix(0, n, nd)
  g_fixed <- if (split && !resume) matrix(0, n, n) else numeric(0)
  et <- numeric(p)
  col_ss <- numeric(p)
  mi_row <- integer(nmiss)
  scratch <- if (split && !resume) numeric(n * min(block, p)) else numeric(0)

  j1 <- 1L
  while (j1 <= p) {
    j2 <- min(j1 + block - 1L, p)
    blk <- if (in_place) X else forward(X[, j1:j2, drop = FALSE], j1:j2)
    if (!is.double(blk)) {
      storage.mode(blk) <- "double"
    }
    .build_block(
      blk,
      if (in_place) j1 - 1L else 0L,
      j1,
      j2 - j1 + 1L,
      if (in_place) 0L else 1L,
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
    j1 <- j2 + 1L
  }
  rm(scratch)
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


.resolve_cols <- function(subset, obj, what = "subset") {
  if (is.null(subset)) {
    return(NULL)
  }
  arg <- .arg_name(what)
  if (length(subset) == 0L) {
    if (identical(what, "abb")) {
      return(NULL)
    }
    stop(
      arg,
      " selects no columns. Name at least one column of 'obj'",
      call. = FALSE
    )
  }
  if (is.character(subset)) {
    cn <- colnames(obj)
    dup <- unique(subset[duplicated(subset)])
    if (length(dup) > 0L) {
      stop(
        arg,
        " names a column more than once: ",
        fmt_trunc(dup, 8),
        call. = FALSE
      )
    }
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
        " not in 'obj': ",
        fmt_trunc(bad, 8),
        call. = FALSE
      )
    }
  } else {
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
  if (!is.integer(got) || length(got) != length(subset)) {
    stop("internal error: ", arg, " did not resolve to one index per request")
  }
  got
}

.mat_miss_int <- function(Z, col = TRUE) {
  m <- mat_miss(Z, col = col, prop = FALSE)
  if (
    !is.integer(m) || length(m) != (if (col) ncol(Z) else nrow(Z)) || anyNA(m)
  ) {
    stop(
      "internal error: mat_miss() did not return one integer count per ",
      if (col) "column" else "row"
    )
  }
  unname(m)
}

.admissible_census <- function(Z) {
  n <- nrow(Z)
  cmiss <- .mat_miss_int(Z)
  v <- col_vars(Z)

  nonfin <- is.nan(v) | is.infinite(v)

  empty <- cmiss == n
  single <- cmiss == n - 1L
  if (!identical(unname(is.na(v) & !nonfin), unname(empty | single))) {
    stop(
      "internal error: col_vars() NA set disagrees with the observed counts"
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

.col_label <- function(Z, k) {
  nms <- colnames(Z)
  if (is.null(nms)) k else nms[k]
}

.col_bullet <- function(what, labels) {
  if (length(labels) == 0L) {
    return(NULL)
  }
  paste0(
    "  ",
    what,
    ": ",
    length(labels),
    " column",
    .s(length(labels)),
    " [",
    fmt_trunc(labels, 8),
    "]"
  )
}

.check_admissible <- function(
  Z,
  lead = "some columns of 'obj' cannot be imputed.",
  fix = paste0(
    "Remove these columns from 'obj'. A column with one observed cell can ",
    "be listed in 'abb' instead"
  )
) {
  cen <- .admissible_census(Z)
  nonfin <- cen$nonfin
  empty <- cen$empty
  single <- cen$single

  if (!any(nonfin | empty | single)) {
    return(invisible(cen))
  }

  stop(
    lead,
    " Each column needs at least two observed cells and a finite ",
    "variance.\n",
    paste(
      c(
        .col_bullet(
          "non-finite variance (an Inf cell, or overflow)",
          .col_label(Z, which(nonfin))
        ),
        .col_bullet("no observed cells", .col_label(Z, which(empty))),
        .col_bullet("only one observed cell", .col_label(Z, which(single)))
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
  store = NULL,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  refresh = 0L,
  spectra = spectra_opts(),
  init = 1,
  forward = identity,
  seed = NULL,
  seed_given = FALSE,
  keep_state = FALSE,
  .state = NULL
) {
  spectra <- .chk_spectra(spectra)
  if (!is.function(forward)) {
    stop("'forward' must be a function", call. = FALSE)
  }
  init_c <- .chk_init(init)
  nu <- .chk_nu(nu)
  if (is.null(.state) && nu < .NU_FLOOR) {
    warning(.nu_floor_msg(), call. = FALSE)
  }
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
  iter_warmup <- .chk_count(iter_warmup, "iter_warmup", lo = 0L)
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  thin <- as.double(.chk_count(thin, "thin"))
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("'Z' must be a numeric matrix", call. = FALSE)
  }

  n <- nrow(Z)
  p <- ncol(Z)
  S <- .chk_rank(rank, n, p)

  K <- min(n - 1L, p)
  dfP <- as.double(n) * p - (p + as.double(S) * (n - 1L + p - S))

  cnt <- if (is.null(.state)) {
    .check_admissible(Z)$cmiss
  } else {
    .mat_miss_int(Z)
  }
  if (sum(cnt) == 0L) {
    stop(
      "'obj' has no missing cells, so there is nothing to resume. Drop ",
      "'keep_state = TRUE'",
      call. = FALSE
    )
  }

  if (!is.null(.state)) {
    .seam_check_match(.state, n, cnt, S, nu, spectra, Z)
  }

  bld <- .core_build(Z, cnt, .state, forward)
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
  inv <- integer(p)
  inv[dirty] <- seq_along(dirty)
  mi_col_l <- inv[mi_col]
  if (min(mi_col_l) < 1L) {
    stop("internal error: a missing cell maps to a clean column")
  }
  rm(inv)
  miss_idx_l <- mi_row + (mi_col_l - 1L) * n

  if (!is.null(store) && !is.null(subset)) {
    stop(
      "'subset' and 'store' cannot both be given: 'subset' stores every hole ",
      "in the columns it names and 'store' names the cells outright",
      call. = FALSE
    )
  }
  cols_req <- .resolve_cols(subset, Z)
  if (is.null(cols_req)) {
    cols_req <- seq_len(p)
  }
  store_use <- if (is.null(.state)) store else .state$store
  if (is.null(store_use)) {
    cols_eff <- cols_req[cnt[cols_req] > 0L]
    store_idx <- sequence(cnt[cols_eff], from = col_ptr[cols_eff])
  } else {
    store_idx <- .chk_store(store_use, nmiss)
    cols_eff <- unique(mi_col[store_idx])
  }
  mi_col_store <- mi_col_l[store_idx]
  map <- cbind(row = mi_row[store_idx], col = mi_col[store_idx])

  niter <- iter_warmup + iter_sampling * thin
  iter_total <- (if (is.null(.state)) 0L else .state$iter) + niter

  if (seed_given) {
    seed <- .resolve_seed(seed)
    dq <- .seed_dq(seed)
    on.exit(dqrng::dqrng_set_state(dq), add = TRUE)
  }
  st <- if (is.null(.state)) {
    .start_fresh(Xd, etd, col_ss, cnt, dirty, nmiss, n, p, init_c, chains)
  } else {
    .start_resume(.state, Xd, etd, dirty, miss_idx_l, n, Z)
  }
  xd_list <- st$xd_list
  etd_list <- st$etd_list
  xtild_list <- st$xtild_list
  eig_start_list <- st$eig_start_list
  sigma_in <- st$sigma_in
  rng_state_in <- st$rng_state_in
  disp <- st$disp
  rm(st, Xd, etd, col_ss)
  draws_list <- .alloc_draws(
    lapply(
      seq_len(chains),
      function(ch) matrix(NA_real_, length(store_idx), iter_sampling)
    ),
    cells = length(store_idx),
    chains = chains,
    draws = iter_sampling,
    fix = "Keep fewer cells with 'subset', or ask for fewer 'chains' or 'iter_sampling'."
  )
  res <- .mipca_da_loop_multi(
    xd_list = xd_list,
    g_fixed = g_fixed,
    etd_list = etd_list,
    xtild_list = xtild_list,
    eig_start_list = eig_start_list,
    sigma_in = sigma_in,
    disp = disp,
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
    warmup = as.integer(iter_warmup),
    ndraws = as.integer(iter_sampling),
    thin = as.integer(thin),
    gram_max_cond = .GRAM_MAX_COND,
    spectra_tol = spectra$tol,
    spectra_ncv = spectra$ncv,
    spectra_maxitr = spectra$maxitr,
    refresh = refresh
  )
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
    nmiss = nmiss,
    seed = if (seed_given) as.integer(seed) else NULL
  )
  if (keep_state) {
    fp_loc <- .seam_fp_cols(nd)
    out$state <- structure(
      list(
        seam_version = .SEAM_VERSION,
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
        fp_clean = if (nd < p) .seam_fp_clean_take(Z, which(cnt == 0L)),
        forward = forward_used,
        store = store_use,
        chain = lapply(seq_len(chains), function(ch) {
          list(
            holes = xd_list[[ch]][miss_idx_l],
            etd = res$etd[[ch]],
            sigma = res$sigma[[ch]],
            xtild_miss = res$xtild_miss[[ch]],
            eig_start = res$eig_start[[ch]],
            rng = res$rng_state[[ch]],
            fp = .seam_fp_take(xd_list[[ch]], fp_loc)
          )
        })
      ),
      class = "si_mi_pca_seam"
    )
  }
  out
}


.resolve_transform <- function(transform, n) {
  msg <- paste0(
    "'transform' must be \"identity\", \"logit_biwhiten\", or a transform ",
    "object: list(name, forward, backward), or list(name, fit) for one that ",
    "is fitted first"
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

#' @export
tf_identity <- function() {
  list(name = "identity", forward = identity, backward = function(z, ...) z)
}

#' @export
tf_sv <- function(n) {
  n <- .chk_count(n, "n")
  if (n < 2L) {
    stop("'n' must be at least 2", call. = FALSE)
  }
  nd <- as.double(n)
  list(
    name = paste0("sv-logit(n=", n, ")"),
    forward = function(y, ...) {
      if (any(y < 0, na.rm = TRUE) || any(y > 1, na.rm = TRUE)) {
        stop(
          "this transform reads proportions, so every observed cell of ",
          "'obj' must lie in [0, 1]. For data on another scale, use ",
          "'transform = \"identity\"'",
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

.trigamma_inverse <- function(x) {
  y <- 0.5 + 1 / x
  for (iter in seq_len(50L)) {
    tri <- trigamma(y)
    dif <- tri * (1 - tri / x) / psigamma(y, deriv = 2L)
    y <- y + dif
    if (max(-dif / y) < 1e-8) {
      return(y)
    }
  }
  stop("internal error: the trigamma inverse did not converge")
}

.fit_fdist <- function(s2, df) {
  e <- log(s2) - digamma(df / 2) + log(df / 2)
  emean <- mean(e)
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

.BIWHITEN_D_MIN <- 2L

#' @export
tf_biwhiten <- function(base = tf_identity(), rows = TRUE) {
  if (
    !is.list(base) || !is.function(base$forward) || !is.function(base$backward)
  ) {
    stop(
      "'base' must be a transform: list(name, forward, backward)",
      call. = FALSE
    )
  }
  if (is.function(base$fit)) {
    stop(
      "'base' must be a transform with no 'fit' step, such as tf_identity() ",
      "or tf_sv()",
      call. = FALSE
    )
  }
  if (!isTRUE(rows) && !isFALSE(rows)) {
    stop("'rows' must be TRUE or FALSE", call. = FALSE)
  }
  nm <- paste0(if (rows) "biwhiten(" else "whiten(", base$name, ")")
  list(
    name = nm,
    fit = function(X, rank) {
      rank <- .chk_count(rank, "rank")
      Z <- base$forward(X)
      obs <- !is.na(Z)
      n <- nrow(Z)
      p <- ncol(Z)
      n_obs <- colSums(obs)
      mu <- colSums(Z, na.rm = TRUE) / pmax(n_obs, 1L)
      Zf <- Z - rep(mu, each = n)
      rm(Z)
      Zf[!obs] <- 0
      e <- eigen(tcrossprod(Zf), symmetric = TRUE)
      W <- e$vectors[, seq_len(rank), drop = FALSE]
      rm(e)
      R <- Zf - W %*% crossprod(W, Zf)
      R[!obs] <- 0
      k <- pmax(0L, pmin(rank, n_obs - 1L - .BIWHITEN_D_MIN))
      full <- n_obs == n & k == rank
      lev_row <- drop(obs %*% (1 / pmax(n_obs, 1L))) + rowSums(W^2) * sum(full)
      for (j in which(!full)) {
        O <- which(obs[, j])
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
      R <- R^2
      df <- pmax(n_obs - 1L, 0L) - k
      rss <- colSums(R)
      s2 <- numeric(p)
      pos <- df > 0L
      s2[pos] <- rss[pos] / df[pos]
      use <- pos & rss > 0
      pr <- .fit_fdist(s2[use], df[use])
      s <- sqrt(.squeeze_var(s2, df, pr$d0, pr$s20))
      if (!all(is.finite(s)) || any(s <= 0)) {
        stop(
          "tf_biwhiten() could not estimate a finite, positive noise scale ",
          "for every column of 'obj'. Use 'transform = \"identity\"' instead",
          call. = FALSE
        )
      }
      d0 <- pr$d0
      s0 <- sqrt(pr$s20)
      t <- rep(1, n)
      if (rows) {
        q <- drop(R %*% (1 / s^2))
        df_i <- rowSums(obs) - lev_row
        has <- df_i > 0 & q > 0
        t[has] <- sqrt(q[has] / df_i[has])
        if (!all(is.finite(t)) || any(t <= 0)) {
          stop(
            "internal error: the row scale is not finite and positive in every row"
          )
        }
        t[has] <- t[has] / exp(mean(log(t[has])))
      }
      rm(
        list = setdiff(ls(all.names = TRUE), c("s", "t", "n", "p", "d0", "s0"))
      )
      list(
        name = nm,
        scale = s,
        rowscale = t,
        d0 = d0,
        s0 = s0,
        forward = function(y, cols) {
          if (length(cols) != ncol(y)) {
            stop(
              "whiten forward: 'cols' must index every column of the block",
              call. = FALSE
            )
          }
          if (NROW(y) != n) {
            stop(
              "whiten forward: the block must carry every row of the fitted matrix",
              call. = FALSE
            )
          }
          z <- base$forward(y) / rep(s[cols], each = NROW(y))
          if (rows) {
            z <- z / t
          }
          z
        },
        backward = function(z, col, row = NULL) {
          if (length(col) != NROW(z)) {
            stop(
              "whiten backward: one column index per cell is required",
              call. = FALSE
            )
          }
          if (any(col < 1L | col > p)) {
            stop(
              "whiten backward: column index outside the fitted matrix",
              call. = FALSE
            )
          }
          if (rows) {
            if (is.null(row)) {
              stop(
                "whiten backward: this transform has a row scale; 'row' is required",
                call. = FALSE
              )
            }
            if (length(row) != NROW(z)) {
              stop(
                "whiten backward: one row index per cell is required",
                call. = FALSE
              )
            }
            if (any(row < 1L | row > n)) {
              stop(
                "whiten backward: row index outside the fitted matrix",
                call. = FALSE
              )
            }
          }
          base$backward(if (rows) z * (s[col] * t[row]) else z * s[col])
        }
      )
    }
  )
}

.chk_frac <- function(x, what) {
  v <- as.double(x)
  if (length(v) != 1L || is.na(v) || v <= 0 || v > 1) {
    stop(.arg_name(what), " must be a single fraction in (0, 1]", call. = FALSE)
  }
  v
}

.arg_name <- function(what) {
  paste0("'", what, "'")
}

.chk_no_dots <- function(..., .what) {
  if (...length() > 0L) {
    nm <- ...names()
    nm <- if (is.null(nm)) character(0) else nm[nzchar(nm)]
    if (length(nm) > 0L) {
      stop(
        .what,
        " takes no other arguments (got ",
        fmt_trunc(nm, 8),
        ")",
        call. = FALSE
      )
    }
    stop(.what, " takes no other arguments", call. = FALSE)
  }
  invisible(NULL)
}

.chk_matrix <- function(obj) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    row.names = "unique",
    col.names = "unique"
  )
}

.chk_count <- function(x, what, lo = 1L, hi = NULL, note = NULL) {
  num <- is.numeric(x) || is.character(x)
  v <- if (num) suppressWarnings(as.integer(x)) else NA_integer_
  d <- if (num) suppressWarnings(as.double(x)) else NA_real_
  bad <- length(v) != 1L || is.na(v) || is.na(d) || v != d || v < lo
  if (!bad && !is.null(hi)) {
    bad <- v > hi
  }
  if (bad) {
    if (is.null(hi)) {
      if (!is.null(note)) {
        stop(
          .arg_name(what),
          " must be a single integer >= ",
          lo,
          " (",
          note,
          ")",
          call. = FALSE
        )
      }
      stop(.arg_name(what), " must be a single integer >= ", lo, call. = FALSE)
    }
    stop(
      .arg_name(what),
      " must be a single index in [",
      lo,
      ", ",
      hi,
      "]",
      call. = FALSE
    )
  }
  v
}

.chk_rank <- function(rank, n, p) {
  S <- .chk_count(rank, "rank")
  if (S + 2L > min(n, p)) {
    stop(
      "'rank' must be at most the smaller dimension minus two (",
      min(n, p) - 2L,
      " here)",
      call. = FALSE
    )
  }
  S
}

.chk_nu <- function(nu) {
  if (!is.numeric(nu) || length(nu) != 1L || is.na(nu) || !(nu > 2)) {
    stop(
      "'nu' must be a single number > 2 (Inf for Gaussian noise)",
      call. = FALSE
    )
  }
  as.double(nu)
}

.nu_floor_msg <- function() {
  paste0(
    "'nu' is below ",
    .NU_FLOOR,
    ", where the sampler can become unstable. ",
    "Use a 'nu' of ",
    .NU_FLOOR,
    " or more"
  )
}


.chk_init <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("'init' must be a single positive number", call. = FALSE)
  }
  as.double(x)
}

.chk_colmax <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    stop("'colmax' must be a single number in [0, 1]", call. = FALSE)
  }
  x
}

.chk_refresh <- function(x) {
  .chk_count(x, "refresh", lo = 0L, note = "0 = no progress lines")
}

.chk_spectra <- function(x) {
  if (!inherits(x, "spectra_opts")) {
    stop("'spectra' must come from spectra_opts()", call. = FALSE)
  }
  x
}

.resolve_seed <- function(seed) {
  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1L)
  }
  .chk_count(seed, "seed", lo = 0L)
}

.seed_dq <- function(seed) {
  seed <- .resolve_seed(seed)
  old <- dqrng::dqrng_get_state()
  dqrng::dqset.seed(seed)
  old
}

.map_lin <- function(map, n, to = NULL) {
  if (is.null(to)) {
    return(unname(map[, "row"] + (map[, "col"] - 1L) * n))
  }
  row <- map[, "row"]
  if (length(row) > 0L && (anyNA(row) || min(row) < 1L || max(row) > n)) {
    stop("internal error: a stored cell's row is outside the matrix")
  }
  unname(row + (.to_col(to, map[, "col"]) - 1) * as.double(n))
}

.alloc_draws <- function(expr, cells, chains, draws, fix) {
  tryCatch(expr, error = function(e) {
    num <- function(x) {
      format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
    }
    gb <- 8 * as.double(cells) * chains * draws / 1024^3
    stop(
      "the draws need ",
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

.chk_cells <- function(
  n,
  p,
  what = "one fit",
  fix = "Fit it in parts by column and join them with cbind()"
) {
  if (as.double(n) * p > .Machine$integer.max) {
    stop(
      "'obj' has ",
      format(as.double(n) * p, big.mark = ",", scientific = FALSE),
      " cells, more than ",
      what,
      " can index (",
      .Machine$integer.max,
      "). ",
      fix,
      call. = FALSE
    )
  }
  invisible(NULL)
}

.chk_abb_draws <- function(n, chains, ndraws, fix) {
  m <- as.double(chains) * ndraws
  if (as.double(n) * m > .Machine$integer.max) {
    stop(
      "the 'abb' columns would draw ",
      format(m, big.mark = ",", scientific = FALSE),
      " values per hole from ",
      n,
      " rows, more than an integer can index (",
      .Machine$integer.max,
      "). ",
      fix,
      call. = FALSE
    )
  }
  invisible(NULL)
}

.lin_col <- function(lin, n) {
  (lin - 1L) %/% n + 1L
}

.forward_at <- function(tf, blk, cols) {
  if (identical(tf$forward, identity)) {
    blk
  } else {
    tf$forward(blk, cols)
  }
}

.backward_at <- function(tf, d, map, rows = NULL) {
  if (is.null(rows)) {
    tf$backward(d, map[, "col"], map[, "row"])
  } else {
    tf$backward(d, map[rows, "col"], map[rows, "row"])
  }
}


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

.tf_unfitted <- function(name) {
  force(name)
  list(name = name, backward = function(z, ...) z)
}

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
  .chk_count(iter_warmup, "iter_warmup", lo = 0L)
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  thin <- .chk_count(thin, "thin")
  nu <- .chk_nu(nu)
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

.mipca_run <- function(
  obj,
  rank,
  transform,
  subset = NULL,
  abb = NULL,
  store = NULL,
  keep_state = FALSE,
  prop_diag = 0.5,
  prop_diag_given = FALSE,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  refresh = 0L,
  spectra = spectra_opts(),
  init = 1,
  seed = NULL,
  seed_given = FALSE,
  ...,
  call = NULL
) {
  .chk_no_dots(..., .what = ".mipca_run()")
  cl <- if (is.null(call)) .scrub_call(match.call(), ".mipca_run") else call
  if (!isTRUE(keep_state) && !isFALSE(keep_state)) {
    stop("'keep_state' must be TRUE or FALSE", call. = FALSE)
  }
  if (keep_state && prop_diag_given) {
    stop(
      "'prop_diag' has no effect with 'keep_state = TRUE'. Pass it to ",
      "summary() or to mipca_finalize() instead",
      call. = FALSE
    )
  }
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  .chk_matrix(obj)
  .chk_cells(nrow(obj), ncol(obj))
  transform <- .resolve_transform(transform, nrow(obj))
  subset <- .resolve_cols(subset, obj)

  abb <- .resolve_cols(abb, obj, what = "abb")
  abb_rec <- NULL
  if (!is.null(abb)) {
    if (!is.null(store)) {
      stop("internal error: 'store' and 'abb' are not given together")
    }
    sp <- .abb_split(obj, abb, subset)
    abb_rec <- sp$rec
    .chk_abb_draws(
      nrow(obj),
      .chk_count(chains, "chains"),
      .chk_count(iter_sampling, "iter_sampling"),
      fix = "Lower 'iter_sampling' or 'chains'"
    )
    if (sp$n_req_kept == 0L) {
      if (keep_state) {
        stop(
          "every requested column is listed in 'abb', so no model is fitted ",
          "and there is nothing to resume. Drop 'keep_state = TRUE'",
          call. = FALSE
        )
      }
      .chk_count(rank, "rank")
      .chk_count(chains, "chains")
      .chk_count(iter_warmup, "iter_warmup", lo = 0L)
      .chk_count(iter_sampling, "iter_sampling")
      .chk_count(thin, "thin")
      .chk_nu(nu)
      abb_rec$seed <- .resolve_seed(seed)
      return(.abb_join(
        NULL,
        abb_rec,
        chains = chains,
        ndraws = iter_sampling,
        thin = thin,
        prop_diag = prop_diag,
        call = cl
      ))
    }
    obj <- obj[, sp$kept, drop = FALSE]
    subset <- sp$subset
    if (keep_state && .nothing_to_impute(obj, subset)) {
      stop(
        "no requested column outside 'abb' has a missing cell, so no model ",
        "is fitted and there is nothing to resume. Drop 'keep_state = TRUE'",
        call. = FALSE
      )
    }
  }

  if (!keep_state && is.null(store) && .nothing_to_impute(obj, subset)) {
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
      call = cl
    )
    if (is.null(abb_rec)) {
      return(empty)
    }
    abb_rec$seed <- .resolve_seed(seed)
    return(.abb_join(
      empty,
      abb_rec,
      chains = chains,
      ndraws = iter_sampling,
      thin = thin,
      prop_diag = prop_diag,
      call = cl
    ))
  }

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
  )

  cl$rank <- core$rank

  draws <- core$draws
  nkeep <- ncol(draws[[1L]])
  prov <- list(
    chains = core$chains,
    rank = core$rank,
    thin = core$thin,
    nu = core$nu,
    warmup = core$iter - core$thin * nkeep,
    ndraws = nkeep,
    iter = core$iter,
    nmiss = core$nmiss,
    thin_final = 1L,
    seam_state = if (keep_state) "kept" else "never",
    seed = core$seed
  )

  if (!is.null(abb_rec)) {
    abb_rec$seed <- .resolve_seed(core$seed)
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

  if (keep_state) fit else mipca_finalize(fit, prop_diag = prop_diag)
}

#' @export
mipca <- function(
  obj,
  rank,
  transform,
  subset = NULL,
  keep_state = FALSE,
  prop_diag = 0.5,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  refresh = 0L,
  spectra = spectra_opts(),
  init = 1,
  seed = NULL,
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


#' @export
crps_all <- function(d, truth) {
  crps_block(d, truth)
}

#' @export
pit_all <- function(d, truth, u = NULL) {
  if (length(truth) != nrow(d)) {
    stop("'truth' must supply one value per row of 'd'", call. = FALSE)
  }
  if (is.null(u)) {
    u <- stats::runif(length(truth))
  }
  if (length(u) != length(truth)) {
    stop("'u' must supply one jitter per cell", call. = FALSE)
  }
  (rowSums(d < truth) + u * (rowSums(d == truth) + 1)) / (ncol(d) + 1)
}
