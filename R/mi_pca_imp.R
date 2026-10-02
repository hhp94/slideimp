.GRAM_MAX_COND <- 1e-8 / .Machine$double.eps
.NU_FLOOR <- 5

#' @export
spectra_opts <- function(tol = 1e-10, ncv = NULL, maxitr = 1000) {
  tol <- as.double(tol)
  if (length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("spectra tol must be a single positive number.")
  }
  maxitr <- .chk_count(maxitr, "spectra maxitr")
  if (is.null(ncv)) {
    ncv <- 0L
  } else {
    ncv <- .chk_count(ncv, "spectra ncv", note = "or NULL")
  }
  structure(
    list(tol = tol, ncv = ncv, maxitr = maxitr),
    class = "spectra_opts"
  )
}

.seam_check_version <- function(.state) {
  if (
    !inherits(.state, "si_mi_pca_seam") ||
      !identical(.state$seam_version, 8L)
  ) {
    stop(
      "the fit's seam is not a version-8 si_mi_pca_seam (it predates the ",
      "current format, or was hand-built); rerun with keep_state = TRUE ",
      "on the current code to get a resumable fit.",
      call. = FALSE
    )
  }
  invisible(.state)
}

.seam_new <- function(chains, iter, rank, thin, spectra, nu, n, cnt, mi_row,
                         g_fixed, fp_clean, forward, store, chain) {
  structure(
    list(
      seam_version = 8L,
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

.chk_store <- function(store, nmiss) {
  s <- as.integer(store)
  if (length(s) < 1L || anyNA(s)) {
    stop("`store` must be a non-empty integer vector with no NA.",
         call. = FALSE)
  }
  if (min(s) < 1L || max(s) > nmiss) {
    stop("`store` positions must lie in 1..", nmiss,
         " (the matrix has ", nmiss, " holes).", call. = FALSE)
  }
  if (length(s) > 1L && min(diff(s)) <= 0L) {
    stop("`store` must be strictly increasing (sorted, no duplicates).",
         call. = FALSE)
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
    stop(
      "the matrix's hole rows do not match the seam's; resume refused.",
      call. = FALSE
    )
  }
  out
}

.start_fresh <- function(Xd, etd, col_ss, cnt, dirty, nmiss, n, p, init_c,
                         chains) {
  nobs <- n - cnt
  if (min(nobs) < 2L) {
    stop("internal error: a column with under two observed cells reached ",
         "the init.")
  }
  df_obs <- as.double(n) * p - nmiss - p
  if (df_obs <= 0) {
    stop("too few observed cells to form a starting sigma.", call. = FALSE)
  }
  list(
    xd_list = lapply(
      seq_len(chains),
      function(ch) if (ch == 1L) Xd else Xd + 0
    ),
    etd_list = rep(list(etd), chains),
    xtild_list = rep(list(numeric(0)), chains),
    sigma_in = rep(sqrt(sum(col_ss) / df_obs), chains),
    rng_state_in = character(0),
    disp = (init_c * sqrt(col_ss / (nobs - 1L)))[dirty]
  )
}

.start_resume <- function(.state, Xd, etd, cnt, dirty, miss_idx_l, n) {
  cnt_d <- cnt[dirty]
  list(
    xd_list = lapply(seq_along(.state$chain), function(ch) {
      st_ch <- .state$chain[[ch]]
      Xc <- Xd + rep(etd - st_ch$etd, each = n)
      Xc[miss_idx_l] <- st_ch$holes
      .seam_fp_check(st_ch$fp, Xc, cnt_d, dirty, ch)
      Xc
    }),
    etd_list = lapply(.state$chain, `[[`, "etd"),
    xtild_list = lapply(.state$chain, `[[`, "xtild_miss"),
    sigma_in = vapply(.state$chain, `[[`, 0, "sigma"),
    rng_state_in = vapply(.state$chain, `[[`, "", "rng"),
    disp = numeric(0)
  )
}

.seam_check_args <- function(.state, chains, chains_given, forward, store) {
  .seam_check_version(.state)
  if (chains_given && as.integer(chains) != .state$chains) {
    stop(
      "chains disagrees with the seam; a resume continues every chain ",
      "of the fit.",
      call. = FALSE
    )
  }
  if (!identical(forward, identity)) {
    stop(
      "`forward` cannot be given on a resume: the seam carries the map ",
      "the fit was built under, and the rebuild applies that one.",
      call. = FALSE
    )
  }
  if (!is.function(.state$forward)) {
    stop("the seam carries no transform map; resume refused.", call. = FALSE)
  }
  if (!is.null(store)) {
    stop(
      "`store` cannot be given on a resume: the seam carries the store set ",
      "the fit was built under, and an extension keeps those same cells.",
      call. = FALSE
    )
  }
  .state$chains
}

.seam_check_match <- function(.state, n, cnt, rank, nu, spectra, Z) {
  if (.state$n != n || !identical(.state$cnt, cnt) || .state$rank != rank) {
    stop("seam state does not match this matrix; resume refused.",
         call. = FALSE)
  }
  if (.state$nu != nu) {
    stop("seam was written with a different nu; resume refused.",
         call. = FALSE)
  }
  if (!identical(.state$spectra, spectra)) {
    stop("seam was written with different spectra settings; resume refused.",
         call. = FALSE)
  }
  .seam_fp_clean_check(.state, Z, n, cnt)
  invisible(.state)
}

.seam_fp_cols <- function(nd) {
  m <- min(nd, max(100L, as.integer(ceiling(0.05 * nd))))
  unique(as.integer(round(seq.int(1L, nd, length.out = m))))
}

.seam_fp_take <- function(Xd, loc, cnt_d) {
  list(
    col = loc,
    sumabs = colSums(abs(Xd[, loc, drop = FALSE])),
    holes = cnt_d[loc]
  )
}

.SEAM_FP_TOL <- 1e-8

.seam_fp_check <- function(fp, Xd, cnt_d, dirty, ch) {
  cur <- .seam_fp_take(Xd, fp$col, cnt_d)
  if (!identical(cur$holes, fp$holes)) {
    stop(
      "resume refused: chain ", ch, "'s rebuilt panel disagrees with the ",
      "seam's fingerprint -- the hole pattern differs in sampled column(s) ",
      fmt_trunc(dirty[fp$col[cur$holes != fp$holes]], 8),
      ".",
      call. = FALSE
    )
  }
  if (!all(is.finite(fp$sumabs)) || any(fp$sumabs <= 0)) {
    stop("internal error: the seam's fingerprint is not positive and finite.")
  }
  rel <- abs(cur$sumabs - fp$sumabs) / fp$sumabs
  if (max(rel) >= .SEAM_FP_TOL) {
    worst <- which.max(rel)
    stop(
      "resume refused: chain ", ch, "'s rebuilt panel disagrees with the ",
      "seam's fingerprint (max rel ", format(max(rel), digits = 3),
      " at column ", dirty[fp$col[worst]], "; hole pattern matches, so the ",
      "matrix values or the transform state moved).",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.seam_fp_clean_take <- function(Z, clean) {
  if (length(clean) == 0L) {
    return(NULL)
  }
  loc <- clean
  sa <- colSums(abs(Z[, loc, drop = FALSE]))
  if (!all(is.finite(sa)) || any(sa <= 0)) {
    stop("internal error: a clean column's fingerprint is not positive ",
         "and finite.")
  }
  list(col = loc, sumabs = sa)
}

.seam_fp_clean_check <- function(.state, Z, n, cnt) {
  split <- .split_cols(n, length(cnt), cnt)$split
  if (!split) {
    return(invisible(NULL))
  }
  fp <- .state$fp_clean
  if (!is.list(fp) || !is.numeric(fp$sumabs) || !is.integer(fp$col) ||
      length(fp$col) != length(fp$sumabs) || length(fp$col) == 0L) {
    stop("the seam carries no clean-column fingerprint but was built from a ",
         "split matrix; resume refused.", call. = FALSE)
  }
  if (max(fp$col) > length(cnt) || any(cnt[fp$col] != 0L)) {
    stop("internal error: the clean fingerprint names a column that is not ",
         "clean in this census.")
  }
  cur <- colSums(abs(Z[, fp$col, drop = FALSE]))
  rel <- abs(cur - fp$sumabs) / fp$sumabs
  if (!all(is.finite(rel)) || max(rel) >= .SEAM_FP_TOL) {
    worst <- which.max(rel)
    stop(
      "resume refused: a clean column's values disagree with the seam's ",
      "fingerprint (max rel ", format(max(rel), digits = 3), " at column ",
      fp$col[worst], "). A resume keeps the clean half's Gram from the ",
      "seam and never rereads those columns, so continuing would sample ",
      "against data this matrix no longer holds.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.split_cols <- function(n, p, cnt) {
  split <- n <= p && any(cnt == 0L)
  list(split = split, dirty = if (split) which(cnt > 0L) else seq_len(p))
}

.mipca_build_stream <- function(X, cnt, forward = identity, block = 1024L,
                                resume = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  block <- max(1L, as.integer(block))
  if (!is.integer(cnt) || length(cnt) != p || anyNA(cnt) || any(cnt < 0L) ||
      any(cnt > n)) {
    stop("internal error: the census does not describe this matrix.")
  }
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
    nb <- j2 - j1 + 1L
    if (in_place) {
      .build_block(
        X, j1 - 1L, j1, nb, 0L, scratch,
        Xd, g_fixed, et, col_ss, mi_row, cnt, col_ptr, dirty_pos,
        if (resume) 1L else 0L
      )
    } else {
      blk <- forward(X[, j1:j2, drop = FALSE], j1:j2)
      if (!is.double(blk)) {
        storage.mode(blk) <- "double"
      }
      .build_block(
        blk, 0L, j1, nb, 1L, scratch,
        Xd, g_fixed, et, col_ss, mi_row, cnt, col_ptr, dirty_pos,
        if (resume) 1L else 0L
      )
    }
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

.resolve_cols <- function(subset, obj) {
  if (is.null(subset)) {
    return(NULL)
  }
  if (length(subset) == 0L) {
    stop("`subset` selected no columns.", call. = FALSE)
  }
  if (is.character(subset)) {
    cn <- colnames(obj)
    if (is.null(cn)) {
      stop(
        "`subset` is character but the matrix has no column names.",
        call. = FALSE
      )
    }
    dup <- unique(subset[duplicated(subset)])
    if (length(dup) > 0L) {
      stop("`subset` repeats: ", fmt_trunc(dup, 8), ".", call. = FALSE)
    }
    bad <- subset[is.na(match(subset, cn))]
    if (length(bad) > 0L) {
      stop(
        "`subset` names ", length(bad), " column",
        .s(length(bad)),
        " that do not exist: ", fmt_trunc(bad, 8),
        ".\nA requested column decides what the fit stores and what ",
        "realization returns, so a name that matches nothing is refused ",
        "rather than dropped.",
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
      .var.name = "subset"
    )
  }
  got <- resolve_subset(subset, obj, sort = TRUE)
  if (!is.integer(got) || length(got) != length(subset)) {
    stop("internal error: `subset` did not resolve to one index per request.")
  }
  got
}

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

.admissible_census <- function(Z, colmax) {
  n <- nrow(Z)
  cmiss <- .mat_miss_int(Z)
  v <- col_vars(Z)

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
  too_miss <- fin & !const & cmiss / n > colmax
  list(
    cmiss = cmiss, nonfin = nonfin, empty = empty, single = single,
    const = const, too_miss = too_miss
  )
}

.check_admissible <- function(Z, colmax) {
  cen <- .admissible_census(Z, colmax)
  cmiss <- cen$cmiss
  nonfin <- cen$nonfin
  empty <- cen$empty
  single <- cen$single
  const <- cen$const
  too_miss <- cen$too_miss

  if (!any(nonfin | empty | single | const | too_miss)) {
    return(invisible(cmiss))
  }

  show <- function(mask, what) {
    k <- which(mask)
    if (length(k) == 0L) {
      return(NULL)
    }
    paste0(
      "  ", what, ": ", length(k), " column",
      .s(length(k)),
      " [", fmt_trunc(k, 8), "]"
    )
  }

  stop(
    "inadmissible columns. Every column reaching the sampler needs finite ",
    "observed values, at least two of them, more than one distinct value ",
    "among them, and missingness at most colmax (", format(colmax), ").\n",
    paste(
      c(
        show(nonfin, "non-finite variance (an Inf cell, or overflow)"),
        show(empty, "no observed cells"),
        show(single, "only one observed cell"),
        show(const, "constant among observed cells"),
        show(too_miss, "missingness above colmax")
      ),
      collapse = "\n"
    ),
    "\nDrop or precondition these columns before fitting; the sampler has no ",
    "posterior for a column it cannot fit.",
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
  parallel_chains = 1L,
  refresh = 0L,
  spectra = spectra_opts(),
  colmax = 0.95,
  init = 1,
  forward = identity,
  keep_state = FALSE,
  .state = NULL
) {
  spectra <- .chk_spectra(spectra)
  if (!is.function(forward)) {
    stop("`forward` must be a function.", call. = FALSE)
  }
  init_c <- .chk_init(init)
  if (!is.numeric(nu) || length(nu) != 1L || is.na(nu) || !(nu > 2)) {
    stop("`nu` must be a single number > 2 (Inf for Gaussian noise).",
         call. = FALSE)
  }
  nu <- as.double(nu)
  if (nu < .NU_FLOOR) {
    warning("`nu` below ", .NU_FLOOR, " makes the I-step's t heavy enough to ",
            "destabilise the Gaussian P-step; 6 to 10 is the measured range.",
            call. = FALSE)
  }
  if (!is.null(.state)) {
    chains <- .seam_check_args(
      .state, chains, !missing(chains), forward, store
    )
  }
  chains <- .chk_count(chains, "chains")
  parallel_chains <- .chk_count(parallel_chains, "parallel_chains")
  refresh <- .chk_refresh(refresh)
  iter_warmup <- .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  thin <- as.double(.chk_count(thin, "thin"))
  colmax <- .chk_colmax(colmax)
  if (!is.matrix(Z) || !is.numeric(Z)) {
    stop("Z must be a numeric matrix.")
  }

  n <- nrow(Z)
  p <- ncol(Z)
  S <- .chk_count(rank, "rank", hi = min(n - 1L, p - 1L))

  if (S + 2L > min(n, p)) {
    stop(
      "rank must be at most the smaller dimension minus two (",
      min(n, p) - 2L,
      " here): the sampler solves rank + 1 eigenpairs, and its eigensolver ",
      "needs that many below the smaller dimension."
    )
  }

  K <- min(n - 1L, p)
  dfP <- n * p - (p + S * (n - 1L + p - S))
  if (dfP <= 0L) {
    stop("non-positive residual df; reduce rank.")
  }

  cnt <- if (is.null(.state)) {
    .check_admissible(Z, colmax)
  } else {
    .mat_miss_int(Z)
  }
  if (sum(cnt) == 0L) {
    stop("Z has no missing values.")
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
    stop("internal error: a missing cell maps to a clean column.")
  }
  rm(inv)
  miss_idx_l <- mi_row + (mi_col_l - 1L) * n

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

  st <- if (is.null(.state)) {
    .start_fresh(Xd, etd, col_ss, cnt, dirty, nmiss, n, p, init_c, chains)
  } else {
    .start_resume(.state, Xd, etd, cnt, dirty, miss_idx_l, n)
  }
  xd_list <- st$xd_list
  etd_list <- st$etd_list
  xtild_list <- st$xtild_list
  sigma_in <- st$sigma_in
  rng_state_in <- st$rng_state_in
  disp <- st$disp
  rm(st, Xd, etd, col_ss)
  draws_list <- lapply(
    seq_len(chains),
    function(ch) matrix(NA_real_, length(store_idx), iter_sampling)
  )
  res <- .mipca_da_loop_multi(
    xd_list = xd_list,
    g_fixed = g_fixed,
    etd_list = etd_list,
    xtild_list = xtild_list,
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
    nthreads = min(parallel_chains, chains),
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
  if (!identical(nrow(draws_list[[1L]]), nrow(map))) {
    stop("internal error: draws rows do not match the map's rows.")
  }
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
    nmiss = nmiss
  )
  if (keep_state) {
    fp_loc <- .seam_fp_cols(nd)
    cnt_d_out <- cnt[dirty]
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
      store = store_use,
      chain = lapply(seq_len(chains), function(ch) {
        list(
          holes = xd_list[[ch]][miss_idx_l],
          etd = res$etd[[ch]],
          sigma = res$sigma[[ch]],
          xtild_miss = res$xtild_miss[[ch]],
          rng = res$rng_state[[ch]],
          fp = .seam_fp_take(xd_list[[ch]], fp_loc, cnt_d_out)
        )
      })
    )
  }
  out
}

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
    is.character(transform$name) && length(transform$name) == 1L &&
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
  n <- .chk_count(n, "`n`")
  if (n < 2L) stop("`n` must be at least 2.", call. = FALSE)
  nd <- as.double(n)
  list(
    name = paste0("sv-logit(n=", n, ")"),
    forward = function(y, ...) {
      if (any(y < 0, na.rm = TRUE) || any(y > 1, na.rm = TRUE)) {
        stop("data must lie in [0, 1].")
      }
      stats::qlogis((y * (nd - 1) + 0.5) / nd)
    },
    backward = function(z, ...) {
      pmin(pmax((stats::plogis(z) * nd - 0.5) / (nd - 1), 0), 1)
    }
  )
}

#' @export
tf_biwhiten <- function(base = tf_identity(), rows = TRUE) {
  if (!is.list(base) || !is.function(base$forward) || !is.function(base$backward)) {
    stop("`base` must be a transform: list(name, forward, backward).", call. = FALSE)
  }
  if (is.function(base$fit)) {
    stop("`base` must be a stateless transform.", call. = FALSE)
  }
  if (!isTRUE(rows) && !isFALSE(rows)) {
    stop("`rows` must be TRUE or FALSE.", call. = FALSE)
  }
  nm <- paste0(if (rows) "biwhiten(" else "whiten(", base$name, ")")
  list(
    name = nm,
    fit = function(X, rank) {
      rank <- .chk_count(rank, "`rank`")
      Z <- base$forward(X)
      obs <- !is.na(Z)
      n_obs <- colSums(obs)
      df <- n_obs - rank
      if (any(df < 1L)) {
        stop(
          "whitening needs more than `rank` observed cells in every column; ",
          sum(df < 1L), " column(s) have too few.",
          call. = FALSE
        )
      }
      mu <- colSums(Z, na.rm = TRUE) / n_obs
      Zf <- Z - rep(mu, each = nrow(Z))
      Zf[!obs] <- 0
      e <- eigen(tcrossprod(Zf), symmetric = TRUE)
      W <- e$vectors[, seq_len(rank), drop = FALSE]
      R <- Zf - W %*% crossprod(W, Zf)
      R[!obs] <- 0
      R2 <- R^2
      s <- sqrt(colSums(R2) / df)
      if (!all(is.finite(s)) || any(s <= 0)) {
        stop("whitening scale is not finite and positive in every column.", call. = FALSE)
      }
      n <- nrow(Z)
      t <- rep(1, n)
      if (rows) {
        df_i <- rowSums(obs) - rank
        if (any(df_i < 1L)) {
          stop(
            "the row scale needs more than `rank` observed cells in every row; ",
            sum(df_i < 1L), " row(s) have too few. Use rows = FALSE.",
            call. = FALSE
          )
        }
        t <- sqrt(rowSums(R2 / rep(s^2, each = nrow(R2))) / df_i)
        if (!all(is.finite(t)) || any(t <= 0)) {
          stop("the row scale is not finite and positive in every row.", call. = FALSE)
        }
        t <- t / exp(mean(log(t)))
      }
      rm(list = setdiff(ls(all.names = TRUE), c("s", "t", "n")))
      p <- length(s)
      list(
        name = nm,
        scale = s,
        rowscale = t,
        forward = function(y, cols) {
          if (length(cols) != ncol(y)) {
            stop("whiten forward: `cols` must index every column of the block.")
          }
          if (NROW(y) != n) {
            stop("whiten forward: the block must carry every row of the fitted matrix.")
          }
          z <- base$forward(y) / rep(s[cols], each = NROW(y))
          if (rows) z <- z / t
          z
        },
        backward = function(z, col, row = NULL) {
          if (length(col) != NROW(z)) {
            stop("whiten backward: one column index per cell is required.")
          }
          if (any(col < 1L | col > p)) {
            stop("whiten backward: column index outside the fitted matrix.")
          }
          if (rows) {
            if (is.null(row)) {
              stop("whiten backward: this transform has a row scale; `row` is required.", call. = FALSE)
            }
            if (length(row) != NROW(z)) {
              stop("whiten backward: one row index per cell is required.")
            }
            if (any(row < 1L | row > n)) {
              stop("whiten backward: row index outside the fitted matrix.")
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
    stop(what, " must be a single fraction in (0, 1].", call. = FALSE)
  }
  v
}

.chk_no_dots <- function(..., .what) {
  if (...length() > 0L) {
    nm <- ...names()
    nm <- if (is.null(nm)) character(0) else nm[nzchar(nm)]
    stop(
      .what, " takes no argument beyond the ones its signature names",
      if (length(nm) > 0L) paste0(" (got ", fmt_trunc(nm, 8), ")") else "",
      ".\nIts `...` is there only because R requires an S3 method to carry ",
      "every formal its generic has; nothing downstream reads it.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.chk_no_ndiag <- function(...) {
  if ("ndiag" %in% ...names()) {
    stop(
      "`ndiag` (a per-column count) is now `prop_diag`, a per-column fraction ",
      "in (0, 1]: a count under-weights the columns holding the slow cells. ",
      "prop_diag = 1 diagnoses every hole.",
      call. = FALSE
    )
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
      stop(
        what, " must be a single integer >= ", lo,
        if (is.null(note)) "" else paste0(" (", note, ")"), ".",
        call. = FALSE
      )
    }
    stop(what, " must be a single index in [", lo, ", ", hi, "].",
         call. = FALSE)
  }
  v
}

.chk_init <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`init` must be a single positive number.", call. = FALSE)
  }
  as.double(x)
}

.chk_colmax <- function(x) {
  if (
    !is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1
  ) {
    stop("`colmax` must be a single number in [0, 1].", call. = FALSE)
  }
  x
}

.chk_refresh <- function(x) {
  .chk_count(x, "`refresh`", lo = 0L, note = "0 = no progress lines")
}

.chk_spectra <- function(x) {
  if (!inherits(x, "spectra_opts")) {
    stop("`spectra` must come from spectra_opts().", call. = FALSE)
  }
  x
}

.seed_dq <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  seed <- .chk_count(seed, "`seed`", lo = 0L)
  old <- dqrng::dqrng_get_state()
  dqrng::dqset.seed(seed)
  old
}

.map_lin <- function(map, n) {
  unname(map[, "row"] + (map[, "col"] - 1L) * n)
}

.lin_col <- function(lin, n) {
  (lin - 1L) %/% n + 1L
}

.backward_at <- function(tf, d, map, rows = NULL) {
  if (is.null(rows)) {
    tf$backward(d, map[, "col"], map[, "row"])
  } else {
    tf$backward(d, map[rows, "col"], map[rows, "row"])
  }
}

.blank_cells <- function(X, idx) {
  X[idx] <- NA_real_
  X
}

.resolve_holdout <- function(obj, holdout, subset = NULL) {
  if (
    !is.matrix(holdout) ||
      !is.numeric(holdout) ||
      ncol(holdout) != 2L ||
      nrow(holdout) < 1L
  ) {
    stop(
      "`holdout` must be a numeric two-column (row, col) matrix ",
      "with at least one row.",
      call. = FALSE
    )
  }
  if (anyNA(holdout)) {
    stop("`holdout` contains NA positions.", call. = FALSE)
  }
  if (any(holdout != trunc(holdout))) {
    stop("`holdout` positions must be whole numbers.", call. = FALSE)
  }
  hr <- as.integer(holdout[, 1L])
  hc <- as.integer(holdout[, 2L])
  n <- nrow(obj)
  if (any(hr < 1L | hr > n) || any(hc < 1L | hc > ncol(obj))) {
    stop("`holdout` positions fall outside `obj`.", call. = FALSE)
  }
  idx <- hr + (hc - 1L) * n
  if (anyDuplicated(idx)) {
    stop("`holdout` contains duplicate positions.", call. = FALSE)
  }
  if (!all(is.finite(obj[idx]))) {
    stop(
      "`holdout` positions must be finite observed cells of `obj`; ",
      "a non-finite cell is a hole, not truth.",
      call. = FALSE
    )
  }
  if (!is.null(subset)) {
    ck <- sort(unique(as.integer(subset)))
    if (!all(hc %in% ck)) {
      stop(
        "`holdout` positions must lie inside `subset`; cells outside it ",
        "are never stored and cannot be scored.",
        call. = FALSE
      )
    }
  }
  hit <- tabulate(hc, nbins = ncol(obj))
  hitc <- which(hit > 0L)
  nobs <- n - colSums(is.na(obj[, hitc, drop = FALSE]))
  if (any(nobs - hit[hitc] < 1L)) {
    stop(
      "`holdout` would blank every observed cell of at least one column.",
      call. = FALSE
    )
  }
  sort(idx)
}

.si_mi_pca_new <- function(obj, map, draws, cols, subset, transform,
                              seam, prov, call) {
  e <- new.env(parent = emptyenv())
  e$obj <- obj
  e$map <- map
  e$draws <- draws
  e$cols <- cols
  e$subset <- subset
  e$transform <- transform
  e$holdout <- NULL
  e$seam <- seam
  e$prov <- prov
  e$call <- call
  e$lifecycle <- "ok"
  class(e) <- "si_mi_pca"
  e
}

.scrub_call <- function(cl, fn) {
  cl[[1L]] <- as.name(fn)
  for (i in seq_along(cl)[-1L]) {
    a <- cl[[i]]
    if (is.symbol(a) || is.call(a) || is.null(a)) {
      next
    }
    if (length(a) <= 1L && !is.object(a)) {
      next
    }
    cl[[i]] <- paste0("<", class(a)[1L], ">")
  }
  cl
}

.si_mi_new <- function(chunks, rownames, nkeep, call) {
  width <- vapply(chunks, function(ch) ncol(ch$obj), 0L)
  nms <- unlist(lapply(chunks, function(ch) colnames(ch$obj)), use.names = FALSE)
  coff <- c(0L, cumsum(width))
  if (
    length(chunks) == 0L || min(width) < 1L || !is.character(nms) ||
      length(nms) != coff[length(coff)]
  ) {
    stop("internal error: a chunk's columns and their names do not line up.")
  }
  structure(
    list(
      chunks = chunks,
      rownames = rownames,
      colnames = nms,
      coff = coff,
      nkeep = nkeep,
      call = call
    ),
    class = "si_mi"
  )
}

.pool_rows <- function(draws, rows) {
  first <- draws[[1L]][rows, , drop = FALSE]
  if (length(draws) == 1L) {
    return(first)
  }
  if (!all(vapply(draws, function(m) is.null(dimnames(m)), TRUE))) {
    stop("internal error: pooled draws carry dimnames, which the fill drops.")
  }
  out <- matrix(0, nrow(first), sum(vapply(draws, ncol, 0L)))
  off <- ncol(first)
  out[, seq_len(off)] <- first
  rm(first)
  for (i in seq_along(draws)[-1L]) {
    k <- ncol(draws[[i]])
    out[, off + seq_len(k)] <- draws[[i]][rows, , drop = FALSE]
    off <- off + k
  }
  out
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

.mipca_empty <- function(obj, subset, rank, transform, prop_diag,
                         iter_sampling, iter_warmup, thin, nu, chains, call,
                         leave_na = FALSE) {
  rank <- .chk_count(rank, "rank")
  chains <- .chk_count(chains, "chains")
  .chk_count(iter_warmup, "`iter_warmup`", lo = 0L)
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  thin <- .chk_count(thin, "thin")
  if (!is.numeric(nu) || length(nu) != 1L || is.na(nu) || !(nu > 2)) {
    stop("`nu` must be a single number > 2 (Inf for Gaussian noise).",
         call. = FALSE)
  }
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
    nu = as.double(nu),
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
  store = NULL,
  holdout = NULL,
  level = c(0.5, 0.9),
  u = NULL,
  keep_state = FALSE,
  prop_diag = .prop_diag_default(),
  prop_diag_given = FALSE,
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  parallel_chains = 1L,
  refresh = 0L,
  spectra = spectra_opts(),
  colmax = 0.95,
  init = 1,
  ...,
  call = NULL
) {
  .chk_no_dots(..., .what = ".mipca_run()")
  if (keep_state && prop_diag_given) {
    stop(
      "`prop_diag` has no effect on a resumable fit: it describes a frozen ",
      "diagnostic frame, and keep_state = TRUE returns one that has none.\n",
      "Ask for the resolution where the frame is actually built: ",
      "summary(fit, prop_diag = ) on the live fit, or ",
      "mipca_finalize(fit, prop_diag = ) when you finish it.",
      call. = FALSE
    )
  }
  .chk_matrix(obj)
  transform <- .resolve_transform(transform, nrow(obj))
  subset <- .resolve_cols(subset, obj)

  if (
    !keep_state && is.null(holdout) && is.null(store) &&
      .nothing_to_impute(obj, subset)
  ) {
    return(.mipca_empty(
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
      call = if (is.null(call)) .scrub_call(match.call(), ".mipca_run") else call
    ))
  }

  Xfit <- obj
  ho <- NULL
  if (!is.null(holdout)) {
    ho <- .resolve_holdout(obj, holdout, subset)
    Xfit <- .blank_cells(obj, ho)
  }

  if (is.function(transform$fit)) {
    transform <- transform$fit(Xfit, rank)
  }

  core <- mipca_core(
    Xfit,
    forward = transform$forward,
    rank = rank,
    subset = subset,
    store = store,
    keep_state = keep_state,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    thin = thin,
    nu = nu,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    spectra = spectra,
    colmax = colmax,
    init = init
  )

  cl <- if (is.null(call)) .scrub_call(match.call(), ".mipca_run") else call
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
    seam_state = if (keep_state) "kept" else "never"
  )

  fit <- .si_mi_pca_new(
    obj = obj,
    map = core$map,
    draws = draws,
    cols = core$cols,
    subset = core$subset,
    transform = transform,
    seam = core$state,
    prov = prov,
    call = cl
  )

  if (!is.null(ho)) {
    fit$holdout <- .holdout_attach(fit, obj, ho, level, u)
  }

  if (keep_state) fit else mipca_finalize(fit, prop_diag = prop_diag)
}

#' @export
mipca <- function(
  obj,
  rank,
  transform,
  subset = NULL,
  store = NULL,
  keep_state = FALSE,
  prop_diag = .prop_diag_default(),
  iter_sampling = 100,
  iter_warmup = 100,
  thin = 1,
  nu = Inf,
  chains = 4L,
  parallel_chains = 1L,
  refresh = 0L,
  spectra = spectra_opts(),
  colmax = 0.95,
  init = 1,
  seed = NULL
) {
  dq <- .seed_dq(seed)
  if (!is.null(dq)) {
    on.exit(dqrng::dqrng_set_state(dq), add = TRUE)
  }
  .mipca_run(
    obj,
    rank = rank,
    transform = transform,
    subset = subset,
    store = store,
    keep_state = keep_state,
    prop_diag = prop_diag,
    prop_diag_given = !missing(prop_diag),
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    thin = thin,
    nu = nu,
    chains = chains,
    parallel_chains = parallel_chains,
    refresh = refresh,
    spectra = spectra,
    colmax = colmax,
    init = init,
    call = .scrub_call(match.call(), "mipca")
  )
}

.holdout_attach <- function(fit, X, ho, level, u = NULL) {
  lin <- .map_lin(fit$map, nrow(X))
  rows <- which(!is.na(X[lin]))
  if (!identical(lin[rows], ho)) {
    if (length(rows) < length(ho)) {
      stop(
        "the store set omits ", length(ho) - length(rows), " of the ",
        length(ho), " holdout cells; a holdout can only be scored on cells ",
        "whose draws were kept.",
        call. = FALSE
      )
    }
    stop(
      "holdout cells do not line up with the map; blanking boundary violated."
    )
  }
  truth <- X[ho]
  d <- .backward_at(
    fit$transform, .pool_rows(fit$draws, rows), fit$map, rows
  )
  if (is.null(u)) {
    u <- stats::runif(length(truth))
  }
  c(
    list(rows = rows, truth = truth, level = level, u = u),
    .score_holdout(d, truth, level, u = u)
  )
}

#' @export
mipca_clone <- function(x, ...) {
  UseMethod("mipca_clone")
}

#' @export
mipca_clone.si_mi_pca <- function(x, ...) {
  .chk_no_dots(..., .what = "mipca_clone()")
  e <- list2env(
    as.list.environment(x, all.names = TRUE),
    parent = emptyenv()
  )
  class(e) <- "si_mi_pca"
  e
}

.life_check <- function(x) {
  st <- x$lifecycle
  if (identical(st, "ok")) {
    return(invisible(x))
  }
  if (identical(st, "failed")) {
    stop(
      "this fit is marked failed: an update() threw partway through it, so ",
      "its chains are not a single length and its seam may not match its ",
      "draws. No repair is attempted -- re-run from mipca().",
      call. = FALSE
    )
  }
  if (identical(st, "passed")) {
    stop(
      "this fit was passed on to the si_mi that mipca_finalize() returned, which ",
      "took its draws and its seam; use that object.",
      call. = FALSE
    )
  }
  stop("this object's lifecycle field is not one of ok / failed / passed.",
       call. = FALSE)
}

#' @export
update.si_mi_pca <- function(
  object,
  iter_sampling,
  ...
) {
  .life_check(object)
  seam <- object$seam
  .seam_check_version(seam)

  if (...length() > 0L) {
    nm <- ...names()
    nm <- nm[nzchar(nm)]
    stop(
      "update() takes no argument beyond `iter_sampling`",
      if (length(nm) > 0L) paste0(" (got ", fmt_trunc(nm, 8), ")") else "",
      ".\nAn extension samples the way the chains it extends did, so rank, ",
      "thin, spectra and the column subset come from the seam, and `init` ",
      "describes a start these chains already have.",
      call. = FALSE
    )
  }

  Xfit <- object$obj
  hidx <- NULL
  if (!is.null(object$holdout)) {
    hrows <- object$holdout$rows
    hidx <- .map_lin(object$map[hrows, , drop = FALSE], nrow(Xfit))
    Xfit <- .blank_cells(Xfit, hidx)
  }

  core <- mipca_core(
    Xfit,
    rank = seam$rank,
    subset = object$subset,
    iter_sampling = iter_sampling,
    iter_warmup = 0L,
    thin = seam$thin,
    spectra = seam$spectra,
    nu = seam$nu,
    keep_state = TRUE,
    .state = seam
  )
  if (!identical(core$map, object$map)) {
    stop("resumed map disagrees with the fit's map; resume refused.")
  }

  object$lifecycle <- "failed"
  newd <- core$draws
  core$draws <- NULL
  for (ch in seq_along(object$draws)) {
    object$draws[[ch]] <- cbind(object$draws[[ch]], newd[[ch]])
    newd[ch] <- list(NULL)
  }
  rm(newd)
  object$seam <- core$state
  object$prov$ndraws <- ncol(object$draws[[1L]])
  object$prov$iter <- core$iter

  if (!is.null(object$holdout)) {
    object$holdout <- .holdout_attach(
      object, object$obj, hidx, object$holdout$level, object$holdout$u
    )
  }
  object$lifecycle <- "ok"
  invisible(object)
}

#' @export
mipca_finalize <- function(x, ...) {
  UseMethod("mipca_finalize")
}

#' @export
mipca_finalize.si_mi <- function(x, ...) {
  stop("this fit is already finalized; mipca_finalize() is a one-way transition.")
}

#' @export
update.si_mi <- function(object, ...) {
  states <- vapply(object$chunks, function(ch) ch$prov$seam_state, "")
  if (all(states == "never")) {
    stop(
      "this fit was built with keep_state = FALSE, so no seam was ever ",
      "kept; rerun with keep_state = TRUE to get a resumable fit."
    )
  }
  te <- .prov_agree(object$chunks, "thin_effective")
  if (is.na(te)) {
    stop(
      "this fit was finalized, which drops the seam; its chunks were thinned ",
      "differently, so there is no single spacing to re-run at -- see ",
      "summary() for the per-chunk numbers."
    )
  }
  stop(
    "this fit was finalized, which drops the seam; re-run with thin = ",
    te, " to get this spacing directly."
  )
}

.retire_holdout <- function(x, dl) {
  ho <- x$holdout
  c(
    list(
      pos = x$map[ho$rows, , drop = FALSE],
      truth = ho$truth,
      level = ho$level,
      u = ho$u
    ),
    .score_holdout(
      .backward_at(x$transform, .pool_rows(dl, ho$rows), x$map, ho$rows),
      ho$truth,
      ho$level,
      u = ho$u
    )
  )
}

.tf_narrow <- function(tf, subset) {
  force(tf)
  force(subset)
  list(
    name = tf$name,
    backward = function(z, col, row = NULL) {
      if (any(col < 1L | col > length(subset))) {
        stop("backward: column index outside the finalized chunk.")
      }
      tf$backward(z, subset[col], row)
    }
  )
}

.rebase_subset <- function(x, amap, fr, ho_out) {
  if (identical(x$subset, seq_len(ncol(x$obj)))) {
    return(list(
      obj = x$obj, loc = NULL, diag = fr, holdout = ho_out,
      transform = x$transform
    ))
  }
  inv <- integer(ncol(x$obj))
  inv[x$subset] <- seq_along(x$subset)
  loc <- inv[amap[, "col"]]
  fr$cells$col <- inv[fr$cells$col]
  fr$holdout$col <- inv[fr$holdout$col]
  if (!is.null(ho_out)) {
    ho_out$pos[, "col"] <- inv[ho_out$pos[, "col"]]
  }
  if (min(c(loc, fr$cells$col, fr$holdout$col, ho_out$pos[, "col"], 1L)) < 1L) {
    stop("internal error: a stored hole maps outside the requested subset.")
  }
  obj <- x$obj[, x$subset, drop = FALSE]
  dimnames(obj) <- list(NULL, colnames(obj))
  list(
    obj = obj,
    loc = loc,
    diag = fr,
    holdout = ho_out,
    transform = .tf_narrow(x$transform, x$subset)
  )
}

#' @export
mipca_finalize.si_mi_pca <- function(x, thin = 1L, prop_diag = .prop_diag_default(),
                               num_threads = getOption("mc.cores", 1L), ...) {
  .life_check(x)
  .chk_no_ndiag(...)
  .chk_no_dots(..., .what = "mipca_finalize()")
  k <- .chk_count(thin, "thin")
  nd <- ncol(x$draws[[1L]])
  if (k > nd) {
    stop(
      "thin (", k, ") exceeds the ", nd, " stored draws per chain; ",
      "there would be nothing left to keep."
    )
  }

  dl <- if (k == 1L) {
    x$draws
  } else {
    keep <- seq.int(1L, nd, by = k)
    lapply(x$draws, function(m) m[, keep, drop = FALSE])
  }
  nkeep <- sum(vapply(dl, ncol, 0L))

  alive <- .alive_mask(
    nrow(x$map),
    if (is.null(x$holdout)) integer(0) else x$holdout$rows
  )
  ho_out <- if (is.null(x$holdout)) NULL else .retire_holdout(x, dl)

  fr <- .diag_frame(
    map = x$map,
    draws = dl,
    p = ncol(x$obj),
    nms = colnames(x$obj),
    hrows = x$holdout$rows,
    hscore = ho_out,
    prop_diag = prop_diag,
    num_threads = num_threads
  )

  amap <- if (is.null(alive)) x$map else x$map[alive, , drop = FALSE]
  rb <- .rebase_subset(x, amap, fr, ho_out)

  prov <- x$prov
  prov$thin_final <- k
  prov$thin_effective <- prov$thin * k
  prov$seam_state <- if (identical(prov$seam_state, "never")) {
    "never"
  } else {
    "dropped"
  }

  chunk <- list(
    obj = rb$obj,
    map = if (is.null(rb$loc)) amap else cbind(row = amap[, "row"], col = rb$loc),
    draws = if (is.null(alive)) {
      if (length(dl) == 1L) dl[[1L]] else do.call(cbind, dl)
    } else {
      .pool_rows(dl, alive)
    },
    transform = rb$transform,
    holdout = rb$holdout,
    diag = rb$diag,
    prov = prov,
    call = x$call
  )
  out <- .si_mi_new(
    chunks = list(chunk),
    rownames = rownames(x$obj),
    nkeep = nkeep,
    call = x$call
  )

  x$draws <- NULL
  x$seam <- NULL
  x$lifecycle <- "passed"
  out
}

#' @export
cbind.si_mi <- function(..., deparse.level = 0) {
  parts <- Filter(Negate(is.null), list(...))
  if (!all(vapply(parts, inherits, TRUE, "si_mi"))) {
    stop(
      "cbind() over a finalized fit takes si_mi objects only; a resumable ",
      "si_mi_pca must be finalized first."
    )
  }
  if (length(parts) == 1L) {
    return(parts[[1L]])
  }
  a <- parts[[1L]]
  for (b in parts[-1L]) {
    if (!identical(a$rownames, b$rownames)) {
      stop(
        "chunks disagree on rownames; no reordering or intersection is done."
      )
    }
    if (!identical(a$nkeep, b$nkeep)) {
      stop(
        "chunks disagree on the number of kept draws (", a$nkeep, " vs ",
        b$nkeep, "); min-trimming would falsify the frozen diagnostics."
      )
    }
  }
  nms <- unlist(lapply(parts, `[[`, "colnames"), use.names = FALSE)
  if (anyDuplicated(nms) > 0L) {
    dup <- unique(nms[duplicated(nms)])
    stop(
      "chunks repeat ", length(dup), " column name", .s(length(dup)), " [",
      fmt_trunc(dup, 8), "]; a column may be joined once.",
      call. = FALSE
    )
  }
  out <- .si_mi_new(
    chunks = unlist(lapply(parts, `[[`, "chunks"), recursive = FALSE),
    rownames = a$rownames,
    nkeep = a$nkeep,
    call = .scrub_call(match.call(), "cbind")
  )
  tg <- do.call(rbind, lapply(parts, `[[`, "targets"))
  if (!is.null(tg)) {
    row.names(tg) <- NULL
    out$targets <- tg
  }
  out
}

#' @export
as.matrix.si_mi <- function(x, draw = 1L, ...) {
  .chk_no_dots(..., .what = "as.matrix()")
  draw <- .chk_count(draw, "draw", hi = x$nkeep)
  n <- length(x$rownames)
  nc <- length(x$chunks)
  out <- if (nc == 1L) {
    x$chunks[[1L]]$obj
  } else {
    matrix(NA_real_, n, x$coff[nc + 1L])
  }
  for (i in seq_len(nc)) {
    ch <- x$chunks[[i]]
    w <- x$coff[i + 1L] - x$coff[i]
    if (nrow(ch$obj) != n || ncol(ch$obj) != w) {
      stop("internal error: a chunk does not fit its joined columns.")
    }
    if (nc > 1L) {
      out[, x$coff[i] + seq_len(w)] <- ch$obj
    }
    idx <- .map_lin(ch$map, n) + x$coff[i] * as.double(n)
    if (!all(is.na(out[idx]))) {
      stop("internal error: a stored cell is not a hole after retirement.")
    }
    out[idx] <- .backward_at(ch$transform, ch$draws[, draw], ch$map)
  }
  dimnames(out) <- list(x$rownames, x$colnames)
  out
}

#' @export
as.matrix.si_mi_pca <- function(x, ...) {
  .life_check(x)
  stop(
    "as.matrix() is refused on a resumable fit: draw spacing is only ",
    "appropriate for multiple imputation after mipca_finalize(), which pools the ",
    "chains and thins. Call mipca_finalize() first, or inspect_completion() to ",
    "eyeball one mid-run completion (not for analysis)."
  )
}

#' @export
inspect_completion <- function(x, chain = 1L, draw = 1L) {
  if (!inherits(x, "si_mi_pca")) {
    stop("inspect_completion() is for a resumable si_mi_pca fit.")
  }
  .life_check(x)
  chain <- .chk_count(chain, "chain", hi = length(x$draws))
  draw <- .chk_count(draw, "draw", hi = ncol(x$draws[[1L]]))
  out <- x$obj
  idx <- .map_lin(x$map, nrow(out))
  fill <- is.na(out[idx])
  out[idx[fill]] <- .backward_at(
    x$transform, x$draws[[chain]][fill, draw], x$map, fill
  )
  out
}

#' @export
print.si_mi_pca <- function(x, ...) {
  .chk_no_dots(..., .what = "print()")
  if (!identical(x$lifecycle, "ok")) {
    cat(sprintf(
      "si_mi_pca fit, %s -- every method but print() refuses\n",
      if (identical(x$lifecycle, "failed")) {
        "FAILED partway through an update()"
      } else if (identical(x$lifecycle, "passed")) {
        "PASSED ON to the si_mi that mipca_finalize() returned"
      } else {
        "in an unrecognized lifecycle state"
      }
    ))
    return(invisible(x))
  }
  cat(sprintf("si_mi_pca fit, resumable (%s)\n", x$transform$name))
  cat(sprintf("  obj:     %d x %d matrix\n", nrow(x$obj), ncol(x$obj)))
  cat(sprintf(
    "  draws:   %d chain%s x %d stored holes x %d draws, latent scale\n",
    length(x$draws),
    .s(length(x$draws)),
    nrow(x$map),
    ncol(x$draws[[1L]])
  ))
  cat(sprintf(
    "  columns: %d requested, %d with holes (stored)\n",
    length(x$subset),
    length(x$cols)
  ))
  if (!is.null(x$seam$store)) {
    cat(sprintf(
      "  store:   cell-level, %d of %d holes kept\n",
      nrow(x$map), x$prov$nmiss
    ))
  }
  cat(sprintf(
    "  sampled: rank %d, warmup %d, thin %d, %d iterations so far\n",
    x$prov$rank,
    x$prov$warmup,
    x$prov$thin,
    x$prov$iter
  ))
  if (!is.null(x$holdout)) {
    cat(.holdout_line(x$holdout, "blanked", provisional = TRUE))
  }
  if (!is.null(x$prov$tuned)) {
    cat(x$prov$tuned$line)
  }
  cat("  next:    update(fit, iter_sampling) to extend, mipca_finalize(fit) to finish\n")
  cat("           summary(fit) for the per-column diagnostics behind that\n")
  invisible(x)
}

.s <- function(n) {
  if (isTRUE(n == 1L)) "" else "s"
}

.agree <- function(vals) {
  if (length(vals) > 0L && isTRUE(all(vals == vals[1L]))) {
    vals[1L]
  } else {
    vals[NA_integer_]
  }
}

.fmt_agree <- function(v) {
  if (is.na(v)) "mixed" else format(v)
}

.prov_agree <- function(chunks, nm) {
  .agree(vapply(chunks, function(ch) as.integer(ch$prov[[nm]]), 0L))
}

.holdout_line <- function(ho, verb, provisional = FALSE) {
  sprintf(
    "  holdout: %d cells %s, %scoverage %s\n",
    length(ho$truth),
    verb,
    if (provisional) "PROVISIONAL " else "",
    paste(
      sprintf("%s=%.3f", names(ho$coverage), ho$coverage),
      collapse = ", "
    )
  )
}

#' @export
print.si_mi <- function(x, ...) {
  .chk_no_dots(..., .what = "print()")
  p <- length(x$colnames)
  h <- sum(vapply(x$chunks, function(ch) nrow(ch$map), 0L))
  nch <- .prov_agree(x$chunks, "chains")
  cat(sprintf(
    "si_mi fit, finalized (%d chunk%s)\n",
    length(x$chunks),
    .s(length(x$chunks))
  ))
  cat(sprintf(
    "  obj:     %d x %d, %d imputed cells\n",
    length(x$rownames),
    p,
    h
  ))
  nleft <- sum(vapply(
    x$chunks,
    function(ch) if (isTRUE(ch$prov$left_na)) ncol(ch$obj) else 0L,
    0L
  ))
  if (nleft > 0L) {
    cat(sprintf(
      "  left NA: %d column%s not imputed, holes kept as they were%s\n",
      nleft,
      .s(nleft),
      if (is.null(x$targets)) "" else "; fit$targets says why"
    ))
  }
  cat(sprintf(
    "  draws:   %d kept, pooled from %s chain%s at thin %s (effective %s)\n",
    x$nkeep,
    .fmt_agree(nch),
    .s(nch),
    .fmt_agree(.prov_agree(x$chunks, "thin_final")),
    .fmt_agree(.prov_agree(x$chunks, "thin_effective"))
  ))
  for (ch in x$chunks) {
    if (!is.null(ch$holdout)) {
      cat(.holdout_line(ch$holdout, "retired"))
    }
    if (!is.null(ch$prov$tuned)) {
      cat(ch$prov$tuned$line)
    }
  }
  ncell <- sum(vapply(x$chunks, function(ch) nrow(ch$diag$cells), 0L))
  cat(sprintf(
    "  diag:    %d cell%s, frozen at finalize; summary(fit) to read\n",
    ncell,
    .s(ncell)
  ))
  invisible(x)
}

.prop_diag_default <- function() 0.5

.alive_mask <- function(nrow_map, hrows) {
  if (length(hrows) == 0L) {
    return(NULL)
  }
  alive <- rep(TRUE, nrow_map)
  alive[hrows] <- FALSE
  alive
}

.diag_sample <- function(map, p, alive, prop_diag) {
  cnt <- tabulate(map[, "col"], nbins = p)
  off <- cumsum(c(0L, cnt))
  if (off[p + 1L] != nrow(map)) {
    stop("internal error: map column runs do not tile the map.")
  }
  dirty <- which(cnt > 0L)
  pick <- vector("list", length(dirty))
  for (i in seq_along(dirty)) {
    j <- dirty[i]
    r <- off[j] + seq_len(cnt[j])
    if (!is.null(alive)) {
      r <- r[alive[r]]
    }
    h <- length(r)
    if (h == 0L) {
      next
    }
    k <- min(h, max(1L, as.integer(ceiling(prop_diag * h))))
    pick[[i]] <- r[floor(seq(1, h, length.out = k))]
  }
  ns <- lengths(pick)
  rows <- unlist(pick, use.names = FALSE)
  if (is.null(rows)) {
    rows <- integer(0)
  }
  owner <- rep.int(seq_along(dirty), ns)
  if (length(rows) > 1L && any(diff(rows) <= 0L)) {
    stop("internal error: the diagnostic sample is not ascending.")
  }
  if (!identical(unname(map[rows, "col"]), dirty[owner])) {
    stop("internal error: a sampled cell is outside the column it came from.")
  }
  list(rows = rows)
}

.cell_diag_read <- function(draws, srows, num_threads) {
  res <- mipca_cell_diag(draws, srows, 256L, num_threads)

  if (ncol(draws[[1L]]) %/% 2L <= 5L) {
    res[, 2L] <- NA_real_
  }
  res
}

.diag_frame <- function(map, draws, p, nms, hrows, hscore, prop_diag,
                           num_threads) {
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  num_threads <- .chk_count(num_threads, "num_threads")
  if (is.null(nms) || length(nms) != p) {
    stop("internal error: the diagnostic frame needs one name per column.")
  }
  alive <- .alive_mask(nrow(map), hrows)
  srows <- .diag_sample(map, p, alive, prop_diag)$rows

  S <- ncol(draws[[1L]])
  res <- .cell_diag_read(draws, srows, num_threads)

  cs <- map[srows, "col"]
  cells <- data.frame(
    col = cs,
    name = nms[cs],
    row = map[srows, "row"],
    rhat = res[, 1L],
    ess = res[, 2L],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  hc <- map[hrows, "col"]
  holdout <- data.frame(
    col = hc,
    name = nms[hc],
    row = map[hrows, "row"],
    crps = if (length(hrows) > 0L) hscore$crps else numeric(0),
    pit = if (length(hrows) > 0L) hscore$pit else numeric(0),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  if (nrow(holdout) != length(hrows)) {
    stop("internal error: holdout scores do not line up with the holdout rows.")
  }

  structure(
    list(cells = cells, holdout = holdout),
    prop_diag = prop_diag,
    nchains = length(draws),
    ndraws = S,
    class = "si_mi_summary"
  )
}

#' @export
summary.si_mi_pca <- function(object, prop_diag = .prop_diag_default(),
                              num_threads = getOption("mc.cores", 1L), ...) {
  .life_check(object)
  .chk_no_ndiag(...)
  .chk_no_dots(..., .what = "summary()")
  .diag_frame(
    map = object$map,
    draws = object$draws,
    p = ncol(object$obj),
    nms = colnames(object$obj),
    hrows = object$holdout$rows,
    hscore = object$holdout,
    prop_diag = prop_diag,
    num_threads = num_threads
  )
}

#' @export
summary.si_mi <- function(object, ...) {
  .chk_no_dots(..., .what = "summary()")
  bind <- function(which) {
    do.call(rbind, lapply(seq_along(object$chunks), function(i) {
      d <- object$chunks[[i]]$diag[[which]]
      d$col <- d$col + object$coff[i]
      cbind(chunk = rep.int(i, nrow(d)), d)
    }))
  }
  agree <- function(nm) {
    .agree(unlist(lapply(object$chunks, function(ch) attr(ch$diag, nm))))
  }
  structure(
    list(cells = bind("cells"), holdout = bind("holdout")),
    prop_diag = agree("prop_diag"),
    nchains = agree("nchains"),
    ndraws = agree("ndraws"),
    class = "si_mi_summary"
  )
}

#' @export
print.si_mi_summary <- function(x, n = 6L, ...) {
  .chk_no_dots(..., .what = "print()")
  cells <- x$cells
  ho <- x$holdout
  num <- function(nm) .fmt_agree(attr(x, nm))
  ncol_aff <- length(unique(c(cells$name, ho$name)))
  pd <- attr(x, "prop_diag")
  how <- if (is.na(pd)) {
    paste0(num("prop_diag"), " of each column's holes")
  } else if (pd >= 1) {
    "every hole"
  } else {
    sprintf("%g%% of each column's holes", 100 * pd)
  }
  cat(sprintf(
    "si_mi diagnostics: %d cell%s (%s) over %d affected column%s\n",
    nrow(cells),
    .s(nrow(cells)),
    how,
    ncol_aff,
    .s(ncol_aff)
  ))
  cat(sprintf(
    "  from %s chains x %s draws; Rhat/ESS on holes, CRPS/PIT on %d holdout cells\n",
    num("nchains"),
    num("ndraws"),
    nrow(ho)
  ))
  shape <- function(v, nm) {
    v <- v[!is.na(v)]
    if (length(v) == 0L) {
      cat(sprintf("  %-8s none\n", nm))
      return(invisible(NULL))
    }
    cat(sprintf(
      "  %-8s %s\n",
      nm,
      paste(sprintf("%.3g", stats::fivenum(v)), collapse = " ")
    ))
  }
  cat("  (min, lower hinge, median, upper hinge, max across cells)\n")
  shape(cells$rhat, "rhat")
  shape(cells$ess, "ess")
  shape(ho$crps, "crps")
  shape(ho$pit, "pit")
  bad <- order(cells$rhat, decreasing = TRUE, na.last = NA)
  if (length(bad) > 0L) {
    k <- min(n, length(bad))
    cat(sprintf("  worst %d cell%s by rhat:\n", k, .s(k)))
    print(cells[bad[seq_len(k)], , drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}

.draw_parts <- function(fit) {
  if (inherits(fit, "si_mi_pca")) {
    return(lapply(fit$draws, function(m) {
      list(d = m, tf = fit$transform, col = fit$map[, "col"], row = fit$map[, "row"])
    }))
  }
  if (inherits(fit, "si_mi")) {
    return(lapply(fit$chunks, function(ch) {
      list(d = ch$draws, tf = ch$transform, col = ch$map[, "col"], row = ch$map[, "row"])
    }))
  }
  stop("not a si_mi_pca or si_mi fit.")
}

#' @export
clamp_frac <- function(fit) {
  parts <- .draw_parts(fit)
  hit <- vapply(
    parts,
    function(p) {
      b <- p$tf$backward(p$d, p$col, p$row)
      sum(b == 0) + sum(b == 1)
    },
    0
  )
  sum(hit) / sum(vapply(parts, function(p) length(p$d), 0))
}

#' @export
crps_all <- function(d, truth) {
  crps_block(d, truth)
}

#' @export
pit_all <- function(d, truth, u = NULL) {
  if (is.null(u)) {
    u <- stats::runif(length(truth))
  }
  if (length(u) != length(truth)) {
    stop("`u` must supply one jitter per cell.", call. = FALSE)
  }
  rowMeans(d < truth) + u * rowMeans(d == truth)
}

.pit_inside <- function(pit, level) {
  a <- (1 - level) / 2
  (outer(pit, a, ">=") & outer(pit, 1 - a, "<=")) * 1
}

.score_holdout <- function(d, truth, level, u = NULL) {
  pit <- pit_all(d, truth, u = u)
  cov <- colMeans(.pit_inside(pit, level))
  mu <- rowMeans(d)
  list(
    pit = pit,
    crps = crps_all(d, truth),
    coverage = stats::setNames(cov, sprintf("%g%%", 100 * level)),
    mean = mu,
    var = rowSums((d - mu)^2) / (ncol(d) - 1L)
  )
}
