.GRAM_MAX_COND <- 1e-8 / .Machine$double.eps

.NU_FLOOR <- 5

#' Eigensolver Settings for Bayesian PCA Imputation
#'
#' Settings for the eigensolver that [mipca()], [mipca_tune()], [mipca_fit()]
#' and [mipca_flank()] run at every iteration of the sampler. The defaults
#' rarely need changing.
#'
#' @details
#' The eigensolver is Spectra's. These settings change how it converges.
#'
#' With `ncv = NULL`, `ncv` is `2 * rank + 3`, at least `20` and at most the
#' smaller dimension of `obj` once the `abb` columns are left out. A value given
#' by hand must be more than `rank + 1` and at most that dimension.
#'
#' If a fit stops with "the PCA step did not converge", raise `maxitr` or `ncv`.
#'
#' @param tol Numeric. Convergence tolerance of the eigenvalues, relative to
#'   their size. Defaults to `1e-10`.
#' @param ncv Integer or `NULL`. Number of Lanczos vectors the eigensolver
#'   keeps; more vectors converge in fewer iterations and use more memory.
#'   Defaults to `NULL`.
#' @param maxitr Integer. Maximum number of eigensolver iterations at each
#'   iteration of the sampler. Defaults to `1000`.
#'
#' @returns A `spectra_opts` object, to pass as `spectra`.
#'
#' @seealso [mipca()], whose `spectra` argument this builds.
#'
#' @examples
#' set.seed(1)
#' sim <- sim_mat(n = 50, p = 20, beta = TRUE)
#'
#' # More Lanczos vectors and iterations than the defaults
#' opts <- spectra_opts(ncv = 15, maxitr = 5000)
#' fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", spectra = opts, seed = 1)
#' fit
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

.split_cols <- function(n, p, cnt) {
  split <- n <= p && any(cnt == 0L)
  list(split = split, dirty = if (split) which(cnt > 0L) else seq_len(p))
}

.mipca_build_stream <- function(X, cnt, forward, resume) {
  n <- nrow(X)
  p <- ncol(X)
  block <- 1024L
  in_place <- identical(forward, identity)
  if (!is.double(X)) {
    stop(
      "internal error: the sampler was given a matrix that is not double",
      call. = FALSE
    )
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
    # a transform of the caller's own returns what it returns: it is
    # refused, not converted, as every input is
    if (!is.double(blk)) {
      stop(
        "the transform's 'forward' returned values of type \"",
        typeof(blk),
        "\". It must return double values",
        call. = FALSE
      )
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

# What the sampler reads off a scan of `Z` (R/col_scan.R), taken with the
# variances: per column the holes, whether it is empty or has one observed
# cell, whether its variance is not finite with no Inf cell in it (it
# overflowed), and whether it is constant. An Inf cell is in `has_inf`; every
# entry point refuses one before it asks whether a column can be fitted.
# `cols` names the columns the scan read, when it read only some of `Z`.
.admissible_scan <- function(Z, scan, cols = NULL) {
  n <- nrow(Z)
  cmiss <- scan$n_miss
  v <- scan$var
  p_read <- if (is.null(cols)) ncol(Z) else length(cols)
  if (is.null(v) || length(cmiss) != p_read) {
    stop(
      "internal error: the sampler needs a scan of every column it reads, ",
      "with the variances",
      call. = FALSE
    )
  }

  nonfin <- is.nan(v) | is.infinite(v)

  empty <- cmiss == n
  single <- cmiss == n - 1L
  if (!identical(is.na(v) & !nonfin, empty | single)) {
    stop(
      "internal error: the scan's NA variances disagree with its counts",
      call. = FALSE
    )
  }
  fin <- !empty & !single & !nonfin
  const <- fin & v < .Machine$double.eps
  list(
    cmiss = cmiss,
    has_inf = scan$has_inf,
    overflow = nonfin & !scan$has_inf,
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

# The columns of `Z` the model cannot fit, refused with their labels: a
# variance that overflowed, no observed cell, or one. `scan` is the scan of
# `Z` with the variances. The entry points refuse an Inf cell, and a column
# with no observed cell, in 'obj' before they get here (refuse_nonfinite()),
# so an empty column is refused here only in a matrix the sampler built,
# such as 'obj' with the holdout of mipca_tune() blanked.
.check_admissible <- function(
  Z,
  scan,
  lead = "some columns of 'obj' cannot be imputed.",
  fix = paste0(
    "Remove these columns from 'obj'. A column with one observed cell can ",
    "be listed in 'abb' instead"
  )
) {
  adm <- .admissible_scan(Z, scan)
  if (any(adm$has_inf)) {
    stop(
      "internal error: .check_admissible() was reached with an Inf cell",
      call. = FALSE
    )
  }
  overflow <- adm$overflow
  empty <- adm$empty
  single <- adm$single

  if (!any(overflow | empty | single)) {
    return(invisible(adm))
  }

  stop(
    lead,
    " Each column needs at least two observed cells and a finite ",
    "variance.\n",
    paste(
      c(
        .col_bullet(
          "non-finite variance (overflow)",
          .col_label(Z, which(overflow))
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
  .state = NULL,
  scan = NULL
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
    chains <- .state_check_args(
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
  if (!is.matrix(Z) || !is.double(Z)) {
    stop("internal error: mipca_core() takes a double matrix", call. = FALSE)
  }

  n <- nrow(Z)
  p <- ncol(Z)
  S <- .chk_rank(rank, n, p)

  K <- min(n - 1L, p)
  dfP <- as.double(n) * p - (p + as.double(S) * (n - 1L + p - S))

  # the scan .mipca_run() took of `Z`, or one taken here. A resume refuses
  # nothing: its matrix was admitted when the state was made.
  if (is.null(scan)) {
    scan <- col_scan(Z, var = is.null(.state))
  }
  cnt <- if (is.null(.state)) {
    refuse_nonfinite(scan)
    .check_admissible(Z, scan)$cmiss
  } else {
    scan$n_miss
  }
  if (sum(cnt) == 0L) {
    stop(
      "'obj' has no missing cells, so there is nothing to resume. Drop ",
      "'keep_state = TRUE'",
      call. = FALSE
    )
  }

  if (!is.null(.state)) {
    .state_check_match(.state, n, cnt, S, nu, spectra, Z)
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
  cols_req <- if (is.null(subset)) seq_len(p) else subset
  if (
    !is.integer(cols_req) ||
      length(cols_req) == 0L ||
      anyNA(cols_req) ||
      is.unsorted(cols_req, strictly = TRUE) ||
      cols_req[1L] < 1L ||
      cols_req[length(cols_req)] > p
  ) {
    stop("internal error: 'subset' is not ascending column indices of 'Z'")
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
    fp_loc <- .state_fp_cols(nd)
    out$state <- structure(
      list(
        state_version = .STATE_VERSION,
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
        fp_clean = if (nd < p) .state_fp_clean_take(Z, which(cnt == 0L)),
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
            fp = .state_fp_take(xd_list[[ch]], fp_loc)
          )
        })
      ),
      class = "si_mi_pca_state"
    )
  }
  out
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

.lin_col <- function(lin, n) {
  (lin - 1L) %/% n + 1L
}

# The scan .mipca_run() reads to ask whether there is anything to impute:
# with `subset`, of its columns only, so nothing outside them is read when
# nothing is imputed; otherwise of every column, with the variances, which
# the fit then reuses. A caller's scan of every column, `scan`, is indexed in
# place of a read. An Inf, or a column with no observed cell, among the
# columns read is refused. Returns the scan, with `nothing` TRUE when no
# column read has a hole.
.nothing_to_impute <- function(obj, subset, scan = NULL) {
  read <- if (!is.null(scan)) {
    if (is.null(subset)) scan else scan_cols(scan, subset)
  } else if (is.null(subset)) {
    col_scan(obj, var = TRUE)
  } else {
    col_scan(obj, cols = subset)
  }
  refuse_nonfinite(read)
  read$nothing <- !any(read$n_miss > 0L)
  read
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
    state_status = "never"
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
    state = NULL,
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
  call = NULL,
  scan = NULL
) {
  .chk_no_dots(..., .what = ".mipca_run()")
  # `scan`: the scan of every column of `obj`, with the variances, from
  # a caller that took one of a larger matrix and indexed it (mipca_flank()),
  # so `obj` is not read again. It cannot come with `abb`, which drops
  # columns from `obj`.
  if (!is.null(scan) && (!is.null(abb) || is.null(scan$var))) {
    stop(
      "internal error: .mipca_run() takes a scan with the variances and ",
      "without 'abb'",
      call. = FALSE
    )
  }
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
  }

  # `full` is the scan of every column, with the variances, that the fit
  # reads; without `subset` it is the one the check for nothing to impute
  # took. The check runs where its answer changes the path: without
  # `keep_state` and `store`, where nothing to impute is an empty fit, and
  # with `keep_state` and `abb`, where it is a refusal.
  full <- scan
  check_empty <- if (keep_state) !is.null(abb) else is.null(store)
  read <- if (check_empty) .nothing_to_impute(obj, subset, scan)
  if (is.null(subset) && is.null(full)) {
    full <- read
  }
  if (keep_state && isTRUE(read$nothing)) {
    stop(
      "no requested column outside 'abb' has a missing cell, so no model ",
      "is fitted and there is nothing to resume. Drop 'keep_state = TRUE'",
      call. = FALSE
    )
  }

  if (!keep_state && isTRUE(read$nothing)) {
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

  if (is.null(full)) {
    full <- col_scan(obj, var = TRUE)
  }
  refuse_nonfinite(full)
  if (is.function(transform$fit)) {
    .chk_rank(rank, nrow(obj), ncol(obj))
    .check_admissible(obj, full)
    transform <- transform$fit(obj, rank)
  }

  core <- mipca_core(
    obj,
    scan = full,
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
    state_status = if (keep_state) "kept" else "never",
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
    state = core$state,
    prov = prov,
    call = cl,
    abb = abb_rec
  )

  if (keep_state) fit else mipca_finalize(fit, prop_diag = prop_diag)
}

#' Multiple Imputation by Bayesian PCA
#'
#' Draw multiple imputations for the missing values of a numeric matrix under a
#' low-rank model using the Audigier, Husson and Josse (BayesMIPCA)
#' data-augmentation sampler.
#'
#' @param obj A double matrix with samples in rows and features in columns. Row
#'   and column names are required and must be unique. `NA` and `NaN` values are
#'   imputed.
#' @param rank Integer. Number of principal components in the model. Use
#'   [mipca_tune()] to choose the rank on held-out values.
#' @param transform Character or a transform object. The scale the model is
#'   fitted on: `"identity"`, `"logit_biwhiten"`, or an object from
#'   [tf_identity()], [tf_sv()] or [tf_biwhiten()]. See Details in [mipca()].
#' @param subset Optional character or integer vector. The columns whose draws
#'   are kept and returned, by name or by number, or `NULL` for every column.
#'   Defaults to `NULL`.
#' @param keep_state Logical. If `TRUE`, return a fit that
#'   [`update()`][mipca_finalize] can extend at the cost of memory. Defaults to
#'   `FALSE`.
#' @param prop_diag Numeric scalar between `0` and `1`. Proportion of each
#'   column's missing values on which [`summary()`][summary.si_mi] computes Rhat
#'   and ESS. `0` is refused. Defaults to `0.5`.
#' @param iter_sampling Integer. Number of draws kept per chain. Defaults to
#'   `100`.
#' @param iter_warmup Integer. Number of iterations discarded at the start of
#'   each chain, `0` or more. Defaults to `100`.
#' @param thin Integer. After warmup, each chain keeps one iteration in every
#'   `thin`, so it runs `iter_sampling * thin` iterations to keep
#'   `iter_sampling` draws. Defaults to `1`.
#' @param nu Numeric. Degrees of freedom of the noise drawn for a missing value:
#'   `Inf` for Gaussian noise, or a number above `2` for Student-t noise of the
#'   same variance. Defaults to `Inf`. Can be chosen by [mipca_tune()].
#' @param chains Integer. Number of independent chains. Defaults to `4`.
#' @param refresh Integer. Print a progress line every `refresh` iterations of
#'   each chain, or `0` for none. Defaults to `0`.
#' @param spectra A `spectra_opts` object. Settings for the eigensolver, from
#'   [spectra_opts()].
#' @param init Numeric. How widely each chain's starting values are spread: each
#'   missing value starts at a random draw whose standard deviation is `init`
#'   times its column's observed standard deviation. Must be positive. Defaults
#'   to `1`.
#' @param seed Integer, numeric, or `NULL`. Seed for the sampler. With `NULL`,
#'   the seed is drawn from R's random number stream. Defaults to `NULL`.
#' @param abb Optional character or integer vector. Columns to draw by
#'   approximate Bayesian bootstrap instead of by the model. [abb_cols()]
#'   proposes them. Defaults to `NULL`.
#'
#' @details
#' **The transform.** `"identity"` fits the values as they are.
#' `"logit_biwhiten"` is for proportions such as DNAm beta values, and every
#' observed value must lie in `[0, 1]`; see [tf_sv()] and [tf_biwhiten()]. Draws
#' are returned on the scale of `obj`.
#'
#' **Which columns take part.** Every column of `obj` takes part in the model,
#' and every missing value is drawn at every iteration. `subset` selects only
#' which columns' draws are kept. The result holds the selected columns in the
#' order they have in `obj`.
#'
#' **Size.** The draws hold one number per missing value kept, per chain, per
#' draw. Keep fewer columns with `subset`, or ask for fewer `chains` or
#' `iter_sampling`. If the selected columns have no missing values, the fit is
#' returned without sampling.
#'
#' **Reproducibility.** With the same `seed`, or the same [set.seed()] before a
#' call with `seed = NULL`, the draws are the same.
#'
#' @section Performance:
#' Linking R to a faster BLAS can make the sampler significantly faster:
#' [OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows) on Windows,
#' or Apple's Accelerate on macOS.
#'
#' @returns With `keep_state = FALSE`, a finalized fit of class `si_mi`. It
#'   holds `chains * iter_sampling` draws of every missing value in the selected
#'   columns. Use:
#' - [`as.matrix(fit, draw = b)`][as.matrix.si_mi] for the completed matrix at
#'   draw `b`.
#' - [`summary()`][summary.si_mi] for Rhat and ESS.
#' - [`cbind()`][cbind.si_mi] to join fits made over different columns.
#'
#' With `keep_state = TRUE`, a resumable fit of class `si_mi_pca`. Extend it
#' with [`update()`][mipca_finalize], then finish it with [mipca_finalize()].
#'
#' @seealso
#' - [mipca_tune()] and [mipca_fit()] to choose `rank` and `nu` on held-out
#'   values, then fit.
#' - [mipca_flank_plan()] to fit each target column on a window of its neighbors
#'   by location.
#' - [pca_imp()] for a single imputation.
#'
#' @references Audigier V, Husson F, Josse J (2016). Multiple imputation for
#'   continuous variables using a Bayesian principal component analysis.
#'   *Journal of Statistical Computation and Simulation*, 86(11), 2140-2156.
#'   \doi{10.1080/00949655.2015.1104683}
#'
#' @examples
#' set.seed(1)
#' sim <- sim_mat(n = 50, p = 20, beta = TRUE)
#' fit <- mipca(sim$input, rank = 2, transform = "logit_biwhiten", seed = 1)
#' fit
#' completed <- lapply(c(1, 200, 400), function(b) as.matrix(fit, draw = b))
#' summary(fit)
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

crps_all <- function(d, truth) {
  crps_block(d, truth)
}

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
