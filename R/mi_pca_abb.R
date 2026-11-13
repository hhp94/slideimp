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
    stop("`bound` must be a single number in (0, 0.5).", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0.5 || q > 1) {
    stop("`q` must be a single number in (0.5, 1].", call. = FALSE)
  }
  n_obs <- nrow(obj) - .mat_miss_int(obj)
  if (any(n_obs > 0L)) {
    r <- range(obj, na.rm = TRUE)
    if (r[1L] < 0 || r[2L] > 1) {
      stop(
        "abb_cols() reads a matrix of proportions: every observed cell must ",
        "lie in [0, 1] (this one spans ",
        format(r[1L]),
        " to ",
        format(r[2L]),
        ").",
        call. = FALSE
      )
    }
  }

  lo <- colSums(obj < bound, na.rm = TRUE)
  hi <- colSums(obj > 1 - bound, na.rm = TRUE)
  need <- pmax(ceiling(q * n_obs - 1e-9), n_obs %/% 2L + 1L)
  has <- n_obs >= 1L
  m_lo <- has & lo >= need
  m_hi <- has & hi >= need
  if (any(m_lo & m_hi)) {
    stop("internal error: a column is marked in both tails.")
  }
  unname(which(m_lo | m_hi))
}

.abb_col <- function(x, m) {
  m <- .chk_count(m, "m")
  row <- unname(which(is.na(x)))
  h <- length(row)
  nd <- length(x) - h
  if (h < 1L) {
    stop("internal error: ABB was asked for a column with no hole.")
  }
  if (nd < 1L) {
    stop(
      "a column listed for ABB has no observed cell, so there is no donor ",
      "to draw from.",
      call. = FALSE
    )
  }
  donors <- x[-row]
  if (!all(is.finite(donors))) {
    stop(
      "a column listed for ABB has a non-finite observed cell.",
      call. = FALSE
    )
  }
  if (as.double(max(nd, h)) * m > .Machine$integer.max) {
    stop("internal error: the ABB index for one column overflows an integer.")
  }
  pool <- dqrng::dqsample.int(nd, nd * m, replace = TRUE)
  pick <- dqrng::dqsample.int(nd, h * m, replace = TRUE)
  at <- pick + rep((seq_len(m) - 1L) * nd, each = h)
  if (anyNA(at) || min(at) < 1L || max(at) > nd * m) {
    stop("internal error: an ABB pick falls outside its draw's pool.")
  }
  list(row = row, draws = matrix(donors[pool[at]], h, m))
}

.abb_kept <- function(p, abb) {
  if (length(abb) == p) {
    stop(
      "`abb` lists every column, so nothing is left for the sampler to fit.",
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
      "`na_loc` holds ",
      sum(hit),
      " cell",
      .s(sum(hit)),
      " in ",
      length(bad),
      " column",
      .s(length(bad)),
      " listed in `abb` [",
      fmt_trunc(colnames(obj)[bad], 8),
      "]. A listed column is not ",
      "fitted, so it has no score to give; the holdout lies in kept ",
      "columns only.",
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
    stop("internal error: the re-based holdout falls outside the kept matrix.")
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
    stop("internal error: a requested kept column maps outside the kept set.")
  }

  rec <- NULL
  if (length(req_abb) > 0L) {
    blk <- obj[, req_abb, drop = FALSE]
    cen <- list(
      empty = .mat_miss_int(blk) == nrow(blk),
      nonfin = unname(colSums(is.infinite(blk)) > 0)
    )
    show <- function(mask, why) {
      if (!any(mask)) {
        return(NULL)
      }
      paste0(
        "  ",
        why,
        ": ",
        sum(mask),
        " column",
        .s(sum(mask)),
        " [",
        fmt_trunc(colnames(blk)[mask], 8),
        "]"
      )
    }
    if (any(cen$empty | cen$nonfin)) {
      stop(
        "`abb` lists columns the bootstrap cannot draw. It draws each hole ",
        "from the column's own finite observed cells, so it needs at least ",
        "one.\n",
        paste(
          c(
            show(cen$empty, "no observed cells"),
            show(cen$nonfin, "a non-finite cell")
          ),
          collapse = "\n"
        ),
        call. = FALSE
      )
    }
    pos <- integer(p)
    pos[req] <- seq_along(req)
    joined <- c(req_kept, req_abb)
    placement <- pos[joined]
    if (!identical(req[placement], joined)) {
      stop(
        "internal error: the placement does not return the requested columns to the input's order."
      )
    }
    rec <- list(
      obj = blk,
      key = req_abb,
      placement = placement,
      nlisted = length(abb),
      seed = NULL
    )
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
  seam_state = "never",
  prop_diag,
  call
) {
  if (is.null(rec) || is.null(rec$seed) || is.null(rec$placement)) {
    stop(
      "internal error: the ABB record is missing, or has no seed or no placement yet."
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
    seam_state = seam_state,
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
  seam_state = "never",
  prop_diag = .prop_diag_default(),
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
    stop("internal error: the ABB columns are not ascending indices of `obj`.")
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
      "internal error: ABB needs one distinct positive integer key per column."
    )
  }
  key <- as.integer(key)
  if (!identical(seam_state, "never") && !identical(seam_state, "dropped")) {
    stop("internal error: a finalized piece's seam is never or dropped.")
  }
  seed <- .chk_count(seed, "`seed`", lo = 0L)
  chains <- .chk_count(chains, "chains")
  ndraws <- .chk_count(ndraws, "ndraws")
  thin <- .chk_count(thin, "thin")
  thin_final <- .chk_count(thin_final, "thin_final")
  m <- chains * ndraws

  sub <- obj[, cols, drop = FALSE]
  dimnames(sub) <- list(NULL, colnames(sub))
  nms <- colnames(sub)

  cnt <- .mat_miss_int(sub)
  off <- cumsum(c(0L, cnt))
  nh <- off[k + 1L]
  map <- matrix(0L, nh, 2L, dimnames = list(NULL, c("row", "col")))
  draws <- matrix(NA_real_, nh, m)

  old <- dqrng::dqrng_get_state()
  on.exit(dqrng::dqrng_set_state(old), add = TRUE)
  for (i in seq_len(k)) {
    if (cnt[i] == 0L) {
      next
    }
    dqrng::dqset.seed(c(seed, key[i]))
    r <- tryCatch(
      .abb_col(sub[, i], m),
      error = function(e) {
        stop("column ", nms[i], ": ", conditionMessage(e), call. = FALSE)
      }
    )
    if (length(r$row) != cnt[i]) {
      stop("internal error: an ABB column's holes disagree with the census.")
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
    stop("internal error: the ABB map is not the listed columns' holes.")
  }

  diag <- .diag_frame(
    map = map[integer(0), , drop = FALSE],
    draws = rep(list(matrix(NA_real_, 0L, ndraws)), chains),
    p = k,
    nms = nms,
    prop_diag = prop_diag,
    num_threads = 1L
  )
  prov <- list(
    chains = chains,
    thin = thin,
    ndraws = ndraws,
    nmiss = nh,
    thin_final = thin_final,
    thin_effective = thin * thin_final,
    seam_state = seam_state,
    abb = TRUE,
    seed = seed,
    key = key
  )
  chunk <- list(
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
