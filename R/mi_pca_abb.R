#' Propose Columns to Draw by Approximate Bayesian Bootstrap
#'
#' Find the columns of a matrix of proportions, such as DNAm beta values, whose
#' observed values lie almost all near `0` or almost all near `1`. Pass the
#' result as `abb` to [mipca()], [mipca_tune()] or [mipca_flank_plan()].
#'
#' @param obj A double matrix with samples in rows and features in columns. Row
#'   and column names are required and must be unique. Every observed value must
#'   lie in `[0, 1]`.
#' @param bound Numeric. A value below `bound` counts as near `0`, and a value
#'   above `1 - bound` as near `1`. Must be above `0` and below `0.5`. Defaults
#'   to `0.02`.
#' @param q Numeric. Proportion of a column's observed values that must be near
#'   `0`, or near `1`, for the column to be listed. Must be above `0.5` and at
#'   most `1`. Defaults to `0.95`.
#'
#' @details
#' **The listed columns.** A column listed in `abb` is left out of the model,
#' and its missing values are drawn from its own observed values by approximate
#' Bayesian bootstrap (Rubin and Schenker 1986). Unlike the model, this keeps
#' the exact `0`s and `1`s such columns often hold. Each column is drawn on its
#' own, so its relation to the other columns is not kept. [mipca_tune()] tunes
#' on the other columns, and [`summary()`][summary.si_mi] does not report the
#' listed ones.
#'
#' **Choosing the list.** The defaults list only the most extreme columns. Raise
#' `bound` or lower `q` to list more, or pass any columns to `abb`. A column
#' with no observed values is never listed.
#'
#' @returns An integer vector of column numbers in `obj`, in increasing order.
#'   It is empty when no column is listed, and `abb` takes it as no list.
#'
#' @seealso
#' - [mipca()] to fit with the list.
#' - [mipca_flank_plan()], which leaves the listed columns out of every window.
#'
#' @references Rubin DB, Schenker N (1986). Multiple imputation for interval
#'   estimation from simple random samples with ignorable nonresponse. *Journal
#'   of the American Statistical Association*, 81(394), 366-374.
#'   \doi{10.1080/01621459.1986.10478280}
#'
#' @examples
#' set.seed(1)
#' sim <- sim_mat(n = 50, p = 20, beta = TRUE)
#' obj <- sim$input
#'
#' # Push two columns against 0 and 1
#' obj[, 3] <- obj[, 3] / 100
#' obj[, 4] <- 1 - (1 - obj[, 4]) / 100
#'
#' abb <- abb_cols(obj)
#' abb
#'
#' fit <- mipca(obj, rank = 2, transform = "logit_biwhiten", seed = 1, abb = abb)
#' fit
#' @export
abb_cols <- function(obj, bound = 0.02, q = 0.95) {
  .chk_matrix(obj)
  if (
    !is.numeric(bound) ||
      length(bound) != 1L ||
      is.na(bound) ||
      bound <= 0 ||
      bound >= 0.5
  ) {
    stop("'bound' must be a single number in (0, 0.5)", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0.5 || q > 1) {
    stop("'q' must be a single number in (0.5, 1]", call. = FALSE)
  }
  scan <- col_scan(obj, tails = c(bound, 1 - bound))
  refuse_inf(scan)
  n_obs <- nrow(obj) - scan$n_miss
  # checkmate skips NA and NaN, passes a matrix with no observed cell, and
  # stops at the first cell out of bounds; the range is read only to report
  if (!checkmate::test_numeric(obj, lower = 0, upper = 1)) {
    r <- range(obj, na.rm = TRUE)
    stop(
      "abb_cols() reads a matrix of proportions: every observed cell must ",
      "lie in [0, 1] (this one spans ",
      format(r[1L]),
      " to ",
      format(r[2L]),
      ")",
      call. = FALSE
    )
  }

  need <- pmax(ceiling(q * n_obs - 1e-9), n_obs %/% 2L + 1L)
  m_lo <- scan$n_below >= need
  m_hi <- scan$n_above >= need
  if (any(m_lo & m_hi)) {
    stop("internal error: a column is marked in both tails")
  }
  which(m_lo | m_hi)
}

.abb_col <- function(x, m) {
  m <- .chk_count(m, "m")
  row <- unname(which(is.na(x)))
  h <- length(row)
  nd <- length(x) - h
  if (h < 1L) {
    stop("internal error: ABB was asked for a column with no hole")
  }
  if (nd < 1L) {
    stop("internal error: ABB was asked for a column with no observed cell")
  }
  donors <- x[-row]
  if (!all(is.finite(donors))) {
    stop(
      "internal error: ABB was asked for a column with a non-finite ",
      "observed cell"
    )
  }
  if (as.double(max(nd, h)) * m > .Machine$integer.max) {
    stop("internal error: the ABB index for one column overflows an integer")
  }
  pool <- dqrng::dqsample.int(nd, nd * m, replace = TRUE)
  pick <- dqrng::dqsample.int(nd, h * m, replace = TRUE)
  at <- pick + rep((seq_len(m) - 1L) * nd, each = h)
  if (anyNA(at) || min(at) < 1L || max(at) > nd * m) {
    stop("internal error: an ABB pick falls outside its draw's pool")
  }
  list(row = row, draws = matrix(donors[pool[at]], h, m))
}

.abb_kept <- function(p, abb) {
  if (length(abb) == p) {
    stop(
      "'abb' lists every column, so no column is left to fit. Leave at least ",
      "one column out of 'abb'",
      call. = FALSE
    )
  }
  is_abb <- logical(p)
  is_abb[abb] <- TRUE
  list(is_abb = is_abb, kept = which(!is_abb))
}

.abb_tune_rebase <- function(lin, obj, kept) {
  n <- nrow(obj)
  p <- ncol(obj)
  col <- .lin_col(lin, n)
  inv <- integer(p)
  inv[kept] <- seq_along(kept)
  kc <- inv[col]
  hit <- kc == 0L
  if (any(hit)) {
    bad <- unique(col[hit])
    stop(
      "'na_loc' holds ",
      sum(hit),
      " cell",
      .s(sum(hit)),
      " in ",
      length(bad),
      " column",
      .s(length(bad)),
      " listed in 'abb' [",
      fmt_trunc(colnames(obj)[bad], 8),
      "]. A column listed in 'abb' is ",
      "not fitted, so it cannot be scored. Place the holdout in columns ",
      "that are not listed",
      call. = FALSE
    )
  }
  out <- lin - (col - kc) * n
  if (
    anyNA(out) ||
      is.unsorted(out, strictly = TRUE) ||
      out[1L] < 1L ||
      out[length(out)] > n * length(kept)
  ) {
    stop("internal error: the re-based holdout falls outside the kept matrix")
  }
  out
}

.abb_split <- function(obj, abb, subset) {
  p <- ncol(obj)
  cm <- .abb_kept(p, abb)
  is_abb <- cm$is_abb
  kept <- cm$kept
  req <- if (is.null(subset)) seq_len(p) else subset
  req_abb <- req[is_abb[req]]
  req_kept <- req[!is_abb[req]]
  inv <- integer(p)
  inv[kept] <- seq_along(kept)
  sub_k <- inv[req_kept]
  if (length(sub_k) > 0L && min(sub_k) < 1L) {
    stop("internal error: a requested kept column maps outside the kept set")
  }

  rec <- NULL
  if (length(req_abb) > 0L) {
    # a listed column's holes are drawn from its own observed cells, so it
    # needs one, and no Inf. Only the requested listed columns are read.
    refuse_nonfinite(col_scan(obj, cols = req_abb))
    blk <- obj[, req_abb, drop = FALSE]
    pos <- integer(p)
    pos[req] <- seq_along(req)
    joined <- c(req_kept, req_abb)
    placement <- pos[joined]
    if (!identical(req[placement], joined)) {
      stop(
        "internal error: the placement does not return the requested columns to the input's order"
      )
    }
    rec <- list(obj = blk, key = req_abb, placement = placement, seed = NULL)
  }
  list(
    kept = kept,
    subset = if (length(sub_k) == length(kept)) NULL else sub_k,
    n_req_kept = length(sub_k),
    rec = rec
  )
}

.abb_join <- function(
  kept,
  rec,
  chains,
  ndraws,
  thin,
  thin_final = 1L,
  state_status = "never",
  prop_diag,
  call
) {
  if (is.null(rec) || is.null(rec$seed) || is.null(rec$placement)) {
    stop(
      "internal error: the ABB record is missing, or has no seed or no placement yet"
    )
  }
  pc <- .abb_piece(
    rec$obj,
    cols = seq_len(ncol(rec$obj)),
    seed = rec$seed,
    key = rec$key,
    chains = chains,
    ndraws = ndraws,
    thin = thin,
    thin_final = thin_final,
    state_status = state_status,
    prop_diag = prop_diag,
    call = call
  )
  if (is.null(kept)) {
    return(pc)
  }
  .si_mi_join(list(kept, pc), call = call, placement = rec$placement)
}

.abb_piece <- function(
  obj,
  cols,
  seed,
  key,
  chains,
  ndraws,
  thin = 1L,
  thin_final = 1L,
  state_status = "never",
  prop_diag = 0.5,
  call = NULL
) {
  n <- nrow(obj)
  k <- length(cols)
  if (
    k < 1L ||
      !is.integer(cols) ||
      anyNA(cols) ||
      cols[1L] < 1L ||
      cols[k] > ncol(obj) ||
      is.unsorted(cols, strictly = TRUE)
  ) {
    stop("internal error: the ABB columns are not ascending indices of 'obj'")
  }
  if (
    !is.numeric(key) ||
      length(key) != k ||
      anyNA(key) ||
      min(key) < 1 ||
      max(key) > .Machine$integer.max ||
      any(key != round(key)) ||
      anyDuplicated(key) > 0L
  ) {
    stop(
      "internal error: ABB needs one distinct positive integer key per column"
    )
  }
  key <- as.integer(key)
  if (
    !identical(state_status, "never") && !identical(state_status, "dropped")
  ) {
    stop("internal error: a finalized piece's state is never or dropped")
  }
  seed <- .chk_count(seed, "seed", lo = 0L)
  chains <- .chk_count(chains, "chains")
  ndraws <- .chk_count(ndraws, "ndraws")
  thin <- .chk_count(thin, "thin")
  thin_final <- .chk_count(thin_final, "thin_final")
  m <- chains * ndraws

  sub <- obj[, cols, drop = FALSE]
  dimnames(sub) <- list(NULL, colnames(sub))
  nms <- colnames(sub)

  cnt <- col_scan(sub)$n_miss
  off <- cumsum(c(0L, cnt))
  nh <- off[k + 1L]
  map <- matrix(0L, nh, 2L, dimnames = list(NULL, c("row", "col")))
  draws <- .alloc_draws(
    matrix(NA_real_, nh, m),
    cells = nh,
    chains = chains,
    draws = ndraws,
    fix = "List fewer columns in 'abb', or ask for fewer 'chains' or 'iter_sampling'."
  )

  old <- dqrng::dqrng_get_state()
  on.exit(dqrng::dqrng_set_state(old), add = TRUE)
  for (i in seq_len(k)) {
    if (cnt[i] == 0L) {
      next
    }
    dqrng::dqset.seed(c(seed, key[i]))
    r <- .abb_col(sub[, i], m)
    if (length(r$row) != cnt[i]) {
      stop("internal error: an ABB column's holes disagree with the scan")
    }
    at <- off[i] + seq_len(cnt[i])
    map[at, "row"] <- r$row
    map[at, "col"] <- i
    draws[at, ] <- r$draws
  }
  lin <- .map_lin(map, n, to = seq_len(k))
  if (
    anyNA(draws) || is.unsorted(lin, strictly = TRUE) || !all(is.na(sub[lin]))
  ) {
    stop("internal error: the ABB map is not the listed columns' holes")
  }

  diag <- .diag_frame(
    map = map[integer(0), , drop = FALSE],
    draws = rep(list(matrix(NA_real_, 0L, ndraws)), chains),
    p = k,
    nms = nms,
    prop_diag = prop_diag
  )
  prov <- list(
    thin_effective = thin * thin_final,
    state_status = state_status
  )
  chunk <- .si_mi_chunk(
    obj = sub,
    map = map,
    draws = draws,
    transform = .tf_unfitted("abb"),
    diag = diag,
    prov = prov,
    call = call
  )
  .si_mi_new(
    chunks = list(chunk),
    rownames = rownames(obj),
    nkeep = m,
    call = call
  )
}
