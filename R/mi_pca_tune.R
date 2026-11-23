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
    stop("'rank' must be a vector of whole numbers >= 1, no NA", call. = FALSE)
  }
  if (anyDuplicated(rank)) {
    stop("'rank' has duplicates. Give each rank once", call. = FALSE)
  }
  as.integer(rank)
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
      "with at least one row",
      call. = FALSE
    )
  }
  if (anyNA(holdout)) {
    stop("the holdout contains NA positions", call. = FALSE)
  }
  if (any(holdout != trunc(holdout))) {
    stop("holdout positions must be whole numbers", call. = FALSE)
  }
  hr <- as.integer(holdout[, 1L])
  hc <- as.integer(holdout[, 2L])
  n <- nrow(obj)
  if (any(hr < 1L | hr > n) || any(hc < 1L | hc > ncol(obj))) {
    stop("holdout positions fall outside 'obj'", call. = FALSE)
  }
  idx <- hr + (hc - 1L) * n
  if (anyDuplicated(idx)) {
    stop("the holdout contains duplicate positions", call. = FALSE)
  }
  if (!all(is.finite(obj[idx]))) {
    stop(
      "holdout positions must be observed, finite cells of 'obj'",
      call. = FALSE
    )
  }
  hit <- tabulate(hc, nbins = ncol(obj))
  hitc <- which(hit > 0L)
  nobs <- n - col_scan(obj, cols = hitc)$n_miss
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
    stop("internal error: holdout positions do not index the hole enumeration")
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
    stop("internal error: diag strata do not tile the dirty columns")
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
    stop("internal error: the diagnostic draw is empty")
  }
  pos
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
    stop("internal error: store rows do not round-trip their positions")
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
      "internal error: the arm's holdout rows do not address the recorded cells"
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
    spectra = fit$state$spectra
  )
  bad <- names(got)[!mapply(same, got, rec$run[names(got)])]
  if (length(bad) > 0L) {
    stop(
      "internal error: an arm did not run the way the tuner asked (",
      fmt_trunc(bad, 8),
      ")"
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
      "internal error: the arm's forward map does not return the recorded truth at the holdout cells"
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
    ess_min = .na_red(ess, min),
    ess_p10 = stats::quantile(ess, 0.10, na.rm = TRUE, names = FALSE),
    ess_med = stats::median(ess, na.rm = TRUE),
    ess_ok = .na_red(ess >= rec$ess_bar, mean),
    rhat_med = stats::median(rhat, na.rm = TRUE),
    rhat_p90 = stats::quantile(rhat, 0.90, na.rm = TRUE, names = FALSE),
    rhat_max = .na_red(rhat, max),
    rhat_ok = mean(rhat <= rec$rhat_bar, na.rm = TRUE),
    ndraws = length(fit$draws) * ncol(fit$draws[[1L]]),
    n_diag = length(rec$rows_diag),
    sigma = list(vapply(fit$state$chain, `[[`, 0, "sigma")),
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
  scan <- col_scan(X, rows = TRUE)
  list(
    dim = dim(X),
    dimnames = dimnames(X),
    row_na = scan$row_miss,
    col_na = scan$n_miss,
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
      if (!all(is.finite(rel)) || max(rel) >= .STATE_FP_TOL) {
        what <- "its values"
      }
    }
  }
  if (!is.null(what)) {
    stop(
      "'obj' is not the matrix this table was tuned on: ",
      what,
      " differ. ",
      "Pass that matrix, or run mipca_tune() on this one",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Choose the Rank of Bayesian PCA Imputation on Held-Out Values
#'
#' Hide a set of observed values, fit [mipca()] at each candidate rank, and
#' score how well each rank's draws predict the hidden values. Pass the result
#' to [mipca_fit()] to fit the best rank.
#'
#' @inheritParams mipca
#' @param rank Integer vector. The candidate ranks, each a whole number of at
#'   least `1`, given once.
#' @param na_loc A list or `NULL`. The values to hold out, as one two-column
#'   matrix of row and column positions in a list, the form [sample_na_loc()]
#'   returns with `n_reps = 1`. Each must be an observed value of `obj`. With
#'   `NULL`, the positions are drawn from `n_cols` and `n_rows`. Defaults to
#'   `NULL`.
#' @param n_cols Integer or `NULL`. Number of columns that receive held-out
#'   values. With `NULL`, 10% of the columns, at least `50` and at most `500`.
#'   Defaults to `NULL`.
#' @param n_rows Integer. Number of held-out values in each of the `n_cols`
#'   columns. Leave it out when `na_loc` is given. Defaults to `4`.
#' @param n_diag Integer. Number of missing values, held-out values included, on
#'   which Rhat and ESS are computed for each rank. They are sampled across
#'   columns with few and with many missing values. Defaults to `10000`.
#' @param seed Integer, numeric, or `NULL`. Seed for the held-out values and the
#'   sampler. With `NULL`, both are drawn from R's random number stream.
#'   Defaults to `NULL`.
#'
#' @details
#' **The held-out values.** With `na_loc = NULL`, [sample_na_loc()] picks
#' `n_cols` columns and hides `n_rows` observed values in each. A column is
#' eligible if it is not constant and has at least `n_rows + 2` observed values.
#' Every rank is fitted on the same matrix, with the same values hidden, and is
#' scored on those values only. The holdout needs more than two values. The
#' `abb` columns take no part in tuning, and no value in them is held out.
#'
#' **The rank.** `crps` is the continuous ranked probability score of each
#' held-out value against its draws, on the scale of `obj`, averaged over the
#' held-out values. A lower score is better. [mipca_fit()] fits the rank with
#' the lowest `crps` unless it is given another.
#'
#' **The noise.** `implied_nu` is the Student-t degrees of freedom that the
#' spread of the held-out values implies, `Inf` when they look Gaussian.
#' `implied_lo` and `implied_hi` are its bootstrap interval. `implied_nu` is the
#' raw reading. [mipca_fit()] fits with `nu` set to `implied_nu`, or to `5` if
#' the reading is lower.
#'
#' **Failures.** A rank whose fit fails does not stop the call. Its row holds
#' the message in `error`, and its scores are `NA`.
#'
#' **Reproducibility.** With the same `seed`, the held-out values and the table
#' are the same. A given `seed` leaves R's random number stream as it was. With
#' `seed = NULL`, the same [set.seed()] before the call gives the same table.
#'
#' @returns A data frame of class `si_tune`, one row per rank, with columns:
#' - `rank`: the rank.
#' - `crps`: the mean score of the held-out values. Lower is better.
#' - `implied_nu`, `implied_lo`, `implied_hi`: the degrees of freedom the
#'   held-out values imply, and its 68% bootstrap interval.
#' - `robust_z`: the interquartile range of the standardized held-out values,
#'   divided by `1.349`. Near `1` for Gaussian noise.
#' - `ess_min`, `ess_p10`, `ess_med`: the minimum, 10% quantile and median ESS
#'   over the `n_diag` values.
#' - `ess_ok`: the proportion of them with an ESS of at least `100 * chains`.
#' - `rhat_med`, `rhat_p90`, `rhat_max`: the median, 90% quantile and maximum
#'   Rhat.
#' - `rhat_ok`: the proportion with an Rhat of at most `1.05`.
#' - `ndraws`: the number of draws, over all chains.
#' - `n_diag`: the number of values Rhat and ESS were computed on.
#' - `sigma`: a list column, the noise standard deviation at the end of each
#'   chain.
#' - `time`: the seconds the rank took.
#' - `error`: the message of a failed fit, or `NA`.
#' - `result`: a list column, a data frame with one row per held-out value: its
#'   `row` and `col`, the `truth`, the `mean` and `sd` of its draws, its `crps`
#'   and `pit`, and `truth_lat`, `mean_lat` and `sd_lat` on the scale the model
#'   is fitted on.
#' - `diag`: a list column, a data frame of the `row`, `col`, `rhat` and `ess`
#'   of each of the `n_diag` values.
#'
#' The table also carries the transform, the sampler settings, `abb` and the
#' held-out positions, which [mipca_fit()] reads.
#'
#' @seealso
#' - [mipca_fit()] to fit a rank from this table.
#' - [mipca()] to fit with a rank and `nu` of your own.
#' - [sample_na_loc()] to draw held-out positions to pass as `na_loc`.
#'
#' @examples
#' set.seed(1)
#' sim <- sim_mat(n = 50, p = 20, beta = TRUE)
#' tuned <- mipca_tune(
#'   sim$input,
#'   rank = 1:3,
#'   transform = "logit_biwhiten",
#'   n_cols = 10,
#'   seed = 1
#' )
#' tuned
#' fit <- mipca_fit(tuned, sim$input, seed = 1)
#' fit
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
  .chk_cells(
    n,
    ncol(full),
    what = "mipca_tune()",
    fix = "Tune on part of its columns"
  )
  transform <- .resolve_transform(transform, n)
  n_diag <- .chk_count(n_diag, "n_diag")
  chains <- .chk_count(chains, "chains")
  iter_sampling <- .chk_count(iter_sampling, "iter_sampling")
  if (chains * iter_sampling < 2) {
    stop(
      "mipca_tune() needs at least 2 draws per rank to score the holdout, ",
      "and 'chains' x 'iter_sampling' is 1. Raise 'chains' or 'iter_sampling'",
      call. = FALSE
    )
  }
  iter_warmup <- .chk_count(iter_warmup, "iter_warmup", lo = 0L)
  thin <- as.double(.chk_count(thin, "thin"))
  init <- .chk_init(init)
  refresh <- .chk_refresh(refresh)
  spectra <- .chk_spectra(spectra)
  # one scan of 'obj' serves the refusals and, when the holdout is drawn,
  # its sampling, which reads the row counts as well
  obj_scan <- col_scan(obj, var = TRUE, rows = is.null(na_loc))
  refuse_nonfinite(obj_scan)
  scan <- .check_admissible(obj, obj_scan)
  if (!is.null(seed)) {
    seed <- .chk_count(seed, "seed", lo = 0L)
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
    nobs <- n - scan$cmiss
    eligible <- sum(nobs >= n_rows + 2L & !scan$const)
    if (eligible < n_cols) {
      stop(
        "cannot place a holdout over ",
        n_cols,
        " columns at ",
        n_rows,
        " cells each. That needs ",
        n_cols,
        " non-constant columns with at ",
        "least ",
        n_rows + 2L,
        " observed cells, and 'obj' has ",
        eligible,
        ". Lower 'n_cols' or 'n_rows', or supply 'na_loc'",
        call. = FALSE
      )
    }
    n_ho_cells <- n_cols * n_rows
  } else {
    if (!missing(n_rows)) {
      stop(
        "'n_rows' has no effect when 'na_loc' is given, because 'na_loc' ",
        "places every held-out cell. Drop one of them",
        call. = FALSE
      )
    }
    if (!is.list(na_loc) || length(na_loc) < 1L) {
      stop(
        "'na_loc' must be a list of one two-column (row, col) matrix, as ",
        "sample_na_loc() returns with 'n_reps = 1'",
        call. = FALSE
      )
    }
    if (length(na_loc) > 1L) {
      stop(
        "'na_loc' holds ",
        length(na_loc),
        " designs, and mipca_tune() scores ",
        "one. Call mipca_tune() once per design. The tables do not combine ",
        "into one grid, because each design scores different cells",
        call. = FALSE
      )
    }
    n_ho_cells <- NROW(na_loc[[1L]])
  }

  if (n_ho_cells <= 2L) {
    stop(
      "the holdout has ",
      n_ho_cells,
      " cell",
      .s(n_ho_cells),
      ", and ",
      "mipca_tune() needs more than two to score a rank. Raise 'n_cols' or ",
      "'n_rows', or supply an 'na_loc' with more rows",
      call. = FALSE
    )
  }

  resolve <- function(design, against) {
    tryCatch(
      .resolve_holdout(against, design),
      error = function(e) {
        stop("'na_loc': ", conditionMessage(e), call. = FALSE)
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
        "'n_cols' disagrees with the supplied 'na_loc'. Drop one of them",
        call. = FALSE
      )
    }
    Xb <- .blank_cells(obj, ho_lin)
    .check_admissible(
      Xb,
      col_scan(Xb, var = TRUE),
      lead = paste0(
        "with the 'na_loc' cells held out, some columns of 'obj' cannot be ",
        "imputed."
      ),
      fix = "Hold out fewer cells in these columns"
    )
  }

  if (!is.null(seed)) {
    base_old <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
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
    na_loc <- sample_na_loc_core(
      obj,
      scan = obj_scan,
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
        " were asked for"
      )
    }
    ho_lin <- resolve(na_loc[[1L]], obj)
  }
  hp <- .tune_hole_pos(na0, ho_lin)
  diag_pos <- .tune_diag_draw(hp$holes, n, p, n_diag)
  st <- .tune_store(hp$ho_pos, diag_pos)
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
        spectra = spectra
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
        stop("internal error: a reported column is not a kept position")
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
      "this table has lost the settings mipca_fit() reads. Pass the table ",
      "mipca_tune() returned, or a row subset of it",
      call. = FALSE
    )
  }
  r <- object$rank
  if (!is.numeric(r) || anyDuplicated(r) || !all(r %in% tn$rank)) {
    stop(
      "this table holds rows from more than one mipca_tune() call. Pass the ",
      "table from one call",
      call. = FALSE
    )
  }
  tn
}

.tune_argmin_rank <- function(object) {
  live <- is.na(object$error) & is.finite(object$crps)
  if (!any(live)) {
    stop(
      "no rank in this table has holdout scores, so mipca_fit() cannot ",
      "choose one. The 'error' column gives the cause for each rank. Fix it, ",
      "then run mipca_tune() again",
      call. = FALSE
    )
  }
  r <- object$rank[live]
  r[order(object$crps[live], r)[1L]]
}

#' Fit Bayesian PCA Imputation at a Tuned Rank
#'
#' Fit [mipca()] at a rank from a [mipca_tune()] table, with the `nu`, transform
#' and sampler settings of that tuning run.
#'
#' @inheritParams mipca
#' @param object A `si_tune` object. The table [mipca_tune()] returns, or a
#'   subset of its rows.
#' @param obj A double matrix. The matrix `object` was tuned on, passed again.
#'   The table does not carry it.
#' @param rank Integer or `NULL`. A rank in `object`. With `NULL`, the rank with
#'   the lowest `crps`. Defaults to `NULL`.
#' @param ... Not used. Name every argument after `rank`.
#' @param chains Integer. Number of independent chains. Taken from `object` when
#'   left out.
#' @param iter_warmup Integer. Number of iterations discarded at the start of
#'   each chain, `0` or more. Taken from `object` when left out.
#' @param iter_sampling Integer. Number of draws kept per chain. Taken from
#'   `object` when left out.
#' @param thin Integer. After warmup, each chain keeps one iteration in every
#'   `thin`, as in [mipca()]. Taken from `object` when left out.
#' @param spectra A `spectra_opts` object. Settings for the eigensolver, from
#'   [spectra_opts()]. Taken from `object` when left out.
#' @param init Numeric. How widely each chain's starting values are spread, as
#'   in [mipca()]. Taken from `object` when left out.
#' @param refresh Integer. Print a progress line every `refresh` iterations of
#'   each chain, or `0` for none. Defaults to `100`.
#' @param nu Numeric or `NULL`. Degrees of freedom of the noise drawn for a
#'   missing value, as in [mipca()]. With `NULL`, the `implied_nu` of `rank` in
#'   `object`, or `5` if it is lower. Defaults to `NULL`.
#' @param transform,abb Not used. The fit uses the transform and `abb` stored in
#'   `object`. To change them, run [mipca_tune()] again, or call [mipca()].
#'
#' @details
#' **The matrix.** `obj` must be the same matrix `object` was tuned on. Another
#' matrix is refused. The held-out values are not held out in the fit: every
#' observed value of `obj` is used.
#'
#' **The rank and `nu`.** With `rank = NULL`, the rank is the one with the
#' lowest `crps` among the ranks whose fit did not fail. A given `rank` must be
#' in `object`. When the `implied_nu` is below `5`, the fit uses `5` and a
#' message says so. A `nu` that is given is used as it is, with the warning
#' [mipca()] gives below `5`.
#'
#' **The same fit as [mipca()].** The fit is the one [mipca()] returns for the
#' same `obj`, `rank`, `nu`, transform, `abb`, sampler settings and `seed`.
#'
#' @returns As [mipca()]: with `keep_state = FALSE`, a finalized fit of class
#'   `si_mi`. With `keep_state = TRUE`, a resumable fit of class `si_mi_pca`.
#'
#' @seealso
#' - [mipca_tune()] to make the table.
#' - [mipca()] to fit with a rank and `nu` of your own.
#'
#' @examples
#' set.seed(1)
#' sim <- sim_mat(n = 50, p = 20, beta = TRUE)
#' tuned <- mipca_tune(
#'   sim$input,
#'   rank = 1:3,
#'   transform = "logit_biwhiten",
#'   n_cols = 10,
#'   seed = 1
#' )
#'
#' # The rank with the lowest crps
#' fit <- mipca_fit(tuned, sim$input, refresh = 0, seed = 1)
#' fit
#'
#' # Another tuned rank, with more draws per chain
#' fit3 <- mipca_fit(tuned, sim$input, rank = 3, iter_sampling = 200, refresh = 0, seed = 1)
#' @export
mipca_fit <- function(
  object,
  obj,
  rank = NULL,
  ...,
  chains,
  iter_warmup,
  iter_sampling,
  thin,
  spectra,
  init,
  refresh = 100L,
  keep_state = FALSE,
  prop_diag = 0.5,
  seed = NULL,
  nu = NULL,
  transform,
  abb
) {
  if (!inherits(object, "si_tune")) {
    stop(
      "'object' must be the table mipca_tune() returns. Pass it first and the ",
      "matrix second: mipca_fit(tuned, obj)",
      call. = FALSE
    )
  }
  .chk_no_dots(..., .what = "mipca_fit()")
  owned <- c(transform = !missing(transform), abb = !missing(abb))
  if (any(owned)) {
    hit <- names(owned)[owned]
    hit_q <- paste(sprintf("'%s'", hit), collapse = ", ")
    stop(
      "mipca_fit() does not take ",
      hit_q,
      ". It uses the transform and 'abb' stored in the table. To change ",
      "them, run mipca_tune() again, or call mipca() directly",
      call. = FALSE
    )
  }
  if (missing(obj)) {
    stop(
      "mipca_fit() needs the matrix the table was tuned on, passed again: ",
      "mipca_fit(tuned, obj). The table does not carry it",
      call. = FALSE
    )
  }
  .chk_matrix(obj)
  tn <- .tune_record(object)
  .tune_fp_check(tn$fp, obj)

  rank <- if (is.null(rank)) {
    .tune_argmin_rank(object)
  } else {
    .chk_count(rank, "rank")
  }
  k <- which(object$rank == rank)
  if (length(k) != 1L) {
    tuned <- paste(sort(object$rank), collapse = ", ")
    stop(
      "rank ",
      rank,
      " was not tuned in this table, which has rank",
      .s(length(object$rank)),
      " ",
      tuned,
      ". mipca_fit() reads 'nu' ",
      "from a tuned rank. Tune this rank with mipca_tune() on the same ",
      "'na_loc', or call mipca() with 'nu' set",
      call. = FALSE
    )
  }
  res <- object$result[[k]]
  if (is.null(res) || nrow(res) < 1L) {
    err <- object$error[k]
    stop(
      "rank ",
      rank,
      " has no holdout scores, so mipca_fit() has no 'nu' to ",
      "read. Choose another rank, or call mipca() with 'nu' set",
      if (!is.na(err)) ".\nIts 'error' entry: ",
      if (!is.na(err)) err,
      call. = FALSE
    )
  }
  if (!identical(obj[cbind(res$row, res$col)], res$truth)) {
    stop(
      "'obj' is not the matrix this table was tuned on: its values at the ",
      "holdout cells differ. Pass that matrix, or run mipca_tune() on this ",
      "one",
      call. = FALSE
    )
  }

  if (is.null(res$truth_lat)) {
    stop(
      "this table lacks the cells that mipca_fit() reads 'nu' from. Run ",
      "mipca_tune() again to make a new table",
      call. = FALSE
    )
  }
  nuz <- .tune_nu_of((res$truth_lat - res$mean_lat) / res$sd_lat)
  nu_implied <- nuz$nu
  if (is.null(nu)) {
    if (is.na(nu_implied)) {
      stop(
        "mipca_fit() cannot read 'nu' at rank ",
        rank,
        ": fewer than 4 ",
        "holdout cells give a usable score. Pass 'nu'",
        call. = FALSE
      )
    }
    nu <- max(nu_implied, .NU_FLOOR)
    if (nu > nu_implied) {
      message(sprintf(
        "the implied nu at rank %d is %.3g, under the floor of %g, so the fit uses nu = %g. Pass 'nu' to use another value",
        rank,
        nu_implied,
        nu,
        nu
      ))
    }
  }

  cfg <- tn$run
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
      "configuration have drifted apart"
    )
  }
  carried <- cfg[!given]
  dots <- mget(names(given)[given], envir = environment())

  cl_fit <- .scrub_call(match.call(), "mipca_fit")
  do.call(
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
}

#' @rdname mipca_tune
#' @param x A `si_tune` object.
#' @param ... Not used.
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
  if (!is.list(tn) || !all(c(cols, "error") %in% names(x))) {
    return(NextMethod())
  }
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
