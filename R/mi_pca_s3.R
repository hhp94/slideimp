.si_mi_pca_new <- function(
  obj,
  map,
  draws,
  cols,
  subset,
  transform,
  seam,
  prov,
  call,
  abb = NULL
) {
  e <- new.env(parent = emptyenv())
  e$abb <- abb
  e$obj <- obj
  e$map <- map
  e$draws <- draws
  e$cols <- cols
  e$subset <- subset
  e$transform <- transform
  e$seam <- seam
  e$prov <- prov
  e$call <- call
  e$lifecycle <- "ok"
  class(e) <- "si_mi_pca"
  e
}

.scrub_call <- function(cl, fn) {
  args <- as.list(cl)
  args[[1L]] <- as.name(fn)
  for (i in seq_along(args)[-1L]) {
    a <- args[[i]]
    if (is.symbol(a) || is.call(a) || is.null(a)) {
      next
    }
    if (length(a) <= 1L && !is.object(a)) {
      next
    }
    args[[i]] <- paste0("<", class(a)[1L], ">")
  }
  as.call(args)
}

.si_mi_chunk <- function(
  obj,
  map,
  draws,
  transform,
  diag,
  prov,
  call,
  info = NULL
) {
  if (!is.matrix(obj) || is.null(colnames(obj))) {
    stop("internal error: a chunk's 'obj' is not a matrix with column names")
  }
  if (
    !is.matrix(map) ||
      !is.integer(map) ||
      !identical(colnames(map), c("row", "col"))
  ) {
    stop("internal error: a chunk's map is not an integer (row, col) matrix")
  }
  if (!is.matrix(draws) || nrow(draws) != nrow(map)) {
    stop("internal error: a chunk's draws do not have one row per stored hole")
  }
  if (!is.function(transform$backward)) {
    stop("internal error: a chunk's transform has no backward map")
  }
  if (
    !inherits(diag, "si_mi_summary") ||
      is.null(attr(diag, "prop_diag")) ||
      is.null(attr(diag, "nchains")) ||
      is.null(attr(diag, "ndraws"))
  ) {
    stop("internal error: a chunk's frame is not one .diag_frame() built")
  }
  if (
    !is.character(prov$seam_state) ||
      length(prov$seam_state) != 1L ||
      length(prov$thin_effective) != 1L
  ) {
    stop(
      "internal error: a chunk's prov lacks 'seam_state' or 'thin_effective'"
    )
  }
  ch <- list(
    obj = obj,
    map = map,
    draws = draws,
    transform = transform,
    diag = diag,
    prov = prov,
    call = call
  )
  if (!is.null(info)) {
    ch$info <- info
  }
  ch
}

.si_mi_new <- function(chunks, rownames, nkeep, call, placement = NULL) {
  width <- vapply(chunks, function(ch) ncol(ch$obj), 0L)
  nms <- unlist(
    lapply(chunks, function(ch) colnames(ch$obj)),
    use.names = FALSE
  )
  coff <- c(0L, cumsum(width))
  if (
    length(chunks) == 0L ||
      min(width) < 1L ||
      !is.character(nms) ||
      length(nms) != coff[length(coff)]
  ) {
    stop("internal error: a chunk's columns and their names do not line up")
  }
  if (
    !is.numeric(nkeep) ||
      length(nkeep) != 1L ||
      is.na(nkeep) ||
      nkeep < 1 ||
      nkeep > .Machine$integer.max ||
      nkeep != round(nkeep)
  ) {
    stop("internal error: 'nkeep' is not one positive whole number")
  }
  nkeep <- as.integer(nkeep)
  if (!all(vapply(chunks, function(ch) ncol(ch$draws), 0L) == nkeep)) {
    stop("internal error: a chunk's draws do not hold 'nkeep' columns")
  }
  if (!is.null(placement)) {
    p <- coff[length(coff)]
    if (
      !is.integer(placement) ||
        length(placement) != p ||
        anyNA(placement) ||
        min(placement) < 1L ||
        max(placement) > p ||
        !all(tabulate(placement, nbins = p) == 1L)
    ) {
      stop(
        "internal error: the placement vector is not a permutation of the joined columns"
      )
    }
    if (is.unsorted(placement)) {
      out_nms <- character(p)
      out_nms[placement] <- nms
      nms <- out_nms
    } else {
      placement <- NULL
    }
  }
  structure(
    list(
      chunks = chunks,
      rownames = rownames,
      colnames = nms,
      coff = coff,
      placement = placement,
      nkeep = nkeep,
      call = call
    ),
    class = "si_mi"
  )
}

.chunk_to <- function(x, i) {
  k <- x$coff[i] + seq_len(x$coff[i + 1L] - x$coff[i])
  if (is.null(x$placement)) k else x$placement[k]
}

.to_col <- function(to, col) {
  if (
    length(col) > 0L &&
      (anyNA(col) || min(col) < 1L || max(col) > length(to))
  ) {
    stop("internal error: a stored cell's column is outside its chunk")
  }
  to[col]
}

.pool_rows <- function(draws, rows) {
  first <- draws[[1L]][rows, , drop = FALSE]
  if (length(draws) == 1L) {
    return(first)
  }
  if (!all(vapply(draws, function(m) is.null(dimnames(m)), TRUE))) {
    stop("internal error: pooled draws carry dimnames, which the fill drops")
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

#' @export
mipca_clone <- function(x) {
  if (!inherits(x, "si_mi_pca")) {
    stop(
      "'x' must be a resumable fit, the object mipca() returns with ",
      "'keep_state = TRUE'. A finalized fit (\"si_mi\") is a value: copy it ",
      "with b <- fit",
      call. = FALSE
    )
  }
  .life_check(x)
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
      "this fit is marked failed: an update() call stopped partway through ",
      "it, so its chains may hold different numbers of draws. Fit again ",
      "with mipca()",
      call. = FALSE
    )
  }
  if (identical(st, "passed")) {
    stop(
      "this fit was finalized: mipca_finalize() moved its draws to the ",
      "\"si_mi\" object it returned. Use that object",
      call. = FALSE
    )
  }
  stop(
    "this object is not a fit from mipca(): its 'lifecycle' field is not ",
    "\"ok\", \"failed\" or \"passed\"",
    call. = FALSE
  )
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

  .chk_no_dots(..., .what = "update()")

  core <- mipca_core(
    object$obj,
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
    stop("internal error: the resumed map disagrees with the fit's map")
  }

  object$lifecycle <- "failed"
  newd <- core$draws
  core$draws <- NULL
  for (ch in seq_along(object$draws)) {
    object$draws[[ch]] <- .alloc_draws(
      cbind(object$draws[[ch]], newd[[ch]]),
      cells = nrow(newd[[ch]]),
      chains = length(object$draws),
      draws = ncol(object$draws[[ch]]) + ncol(newd[[ch]]),
      fix = "Ask update() for fewer 'iter_sampling'."
    )
    newd[ch] <- list(NULL)
  }
  rm(newd)
  object$seam <- core$state
  object$prov$ndraws <- ncol(object$draws[[1L]])
  object$prov$iter <- core$iter
  object$lifecycle <- "ok"
  invisible(object)
}

#' @export
mipca_finalize <- function(x, ...) {
  UseMethod("mipca_finalize")
}

#' @export
mipca_finalize.si_mi <- function(x, ...) {
  stop(
    "this fit is already finalized. Pass it to mipca_complete() or ",
    "summary()",
    call. = FALSE
  )
}

#' @export
update.si_mi <- function(object, ...) {
  states <- vapply(object$chunks, function(ch) ch$prov$seam_state, "")
  if (all(states == "never")) {
    stop(
      "update() cannot extend this fit: it holds no sampler state to resume ",
      "from. Fit again with mipca() and 'keep_state = TRUE'",
      call. = FALSE
    )
  }
  te <- .prov_agree(object$chunks, "thin_effective")
  stop(
    "update() cannot extend a finalized fit. Fit again with mipca() and ",
    "'keep_state = TRUE'",
    if (!is.na(te) && te != 1) ", with 'thin = ",
    if (!is.na(te) && te != 1) te,
    if (!is.na(te) && te != 1) "' to keep this fit's thinning",
    call. = FALSE
  )
}

.tf_narrow <- function(tf, subset) {
  force(tf)
  force(subset)
  list(
    name = tf$name,
    backward = function(z, col, row = NULL) {
      if (any(col < 1L | col > length(subset))) {
        stop(
          "backward: column index outside the finalized chunk",
          call. = FALSE
        )
      }
      tf$backward(z, subset[col], row)
    }
  )
}

.rebase_subset <- function(x, fr) {
  if (identical(x$subset, seq_len(ncol(x$obj)))) {
    return(list(obj = x$obj, loc = NULL, diag = fr, transform = x$transform))
  }
  loc <- integer(0)
  if (nrow(x$map) > 0L) {
    inv <- integer(ncol(x$obj))
    inv[x$subset] <- seq_along(x$subset)
    loc <- inv[x$map[, "col"]]
    fr$cells$col <- inv[fr$cells$col]
    if (min(c(loc, fr$cells$col, 1L)) < 1L) {
      stop("internal error: a stored hole maps outside the requested subset")
    }
  }
  obj <- x$obj[, x$subset, drop = FALSE]
  dimnames(obj) <- list(NULL, colnames(obj))
  list(
    obj = obj,
    loc = loc,
    diag = fr,
    transform = .tf_narrow(x$transform, x$subset)
  )
}

#' @export
mipca_finalize.si_mi_pca <- function(x, thin = 1L, prop_diag = 0.5, ...) {
  .life_check(x)
  .chk_no_dots(..., .what = "mipca_finalize()")
  k <- .chk_count(thin, "thin")
  nd <- ncol(x$draws[[1L]])
  if (k > nd) {
    stop(
      "'thin' (",
      k,
      ") is larger than the ",
      nd,
      " draws each chain holds, ",
      "so no draw would be kept. Use a 'thin' of at most ",
      nd,
      call. = FALSE
    )
  }
  if (!is.null(x$abb)) {
    .chk_abb_draws(
      nrow(x$abb$obj),
      length(x$draws),
      ceiling(nd / k),
      fix = "Use a larger 'thin' here"
    )
  }

  dl <- if (k == 1L) {
    x$draws
  } else {
    keep <- seq.int(1L, nd, by = k)
    lapply(x$draws, function(m) m[, keep, drop = FALSE])
  }
  nkeep <- sum(vapply(dl, ncol, 0L))

  fr <- .diag_frame(
    map = x$map,
    draws = dl,
    p = ncol(x$obj),
    nms = colnames(x$obj),
    prop_diag = prop_diag
  )

  rb <- .rebase_subset(x, fr)

  prov <- x$prov
  prov$thin_final <- k
  prov$thin_effective <- prov$thin * k
  prov$seam_state <- if (identical(prov$seam_state, "never")) {
    "never"
  } else {
    "dropped"
  }

  chunk <- .si_mi_chunk(
    obj = rb$obj,
    map = if (is.null(rb$loc)) {
      x$map
    } else {
      cbind(row = x$map[, "row"], col = rb$loc)
    },
    draws = if (length(dl) == 1L) {
      dl[[1L]]
    } else {
      .alloc_draws(
        do.call(cbind, dl),
        cells = nrow(dl[[1L]]),
        chains = length(dl),
        draws = ncol(dl[[1L]]),
        fix = "Pass a larger 'thin' to mipca_finalize(), or fit with fewer 'chains' or 'iter_sampling'."
      )
    },
    transform = rb$transform,
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

  if (!is.null(x$abb)) {
    out <- .abb_join(
      out,
      x$abb,
      chains = length(dl),
      ndraws = ncol(dl[[1L]]),
      thin = prov$thin,
      thin_final = k,
      seam_state = prov$seam_state,
      prop_diag = prop_diag,
      call = x$call
    )
  }

  x$draws <- NULL
  x$seam <- NULL
  x$lifecycle <- "passed"
  out
}

#' @export
cbind.si_mi <- function(..., deparse.level = 1) {
  parts <- Filter(Negate(is.null), list(...))
  if (!all(vapply(parts, inherits, TRUE, "si_mi"))) {
    stop(
      "cbind() joins finalized fits (\"si_mi\") only. Call mipca_finalize() ",
      "on an \"si_mi_pca\" fit first",
      call. = FALSE
    )
  }
  if (length(parts) == 1L) {
    return(parts[[1L]])
  }
  call <- .scrub_call(match.call(), "cbind")
  .si_mi_join(parts, call = call)
}

.si_mi_join <- function(parts, call, placement = NULL) {
  a <- parts[[1L]]
  for (b in parts[-1L]) {
    if (!identical(a$rownames, b$rownames)) {
      stop(
        "the fits have different row names. cbind() joins fits of the same ",
        "rows in the same order",
        call. = FALSE
      )
    }
    if (!isTRUE(a$nkeep == b$nkeep)) {
      stop(
        "the fits hold different numbers of kept draws (",
        a$nkeep,
        " and ",
        b$nkeep,
        "). Fit and finalize every part with the same 'chains', ",
        "'iter_sampling' and 'thin'",
        call. = FALSE
      )
    }
  }
  nms <- unlist(lapply(parts, `[[`, "colnames"), use.names = FALSE)
  if (anyDuplicated(nms) > 0L) {
    dup <- unique(nms[duplicated(nms)])
    stop(
      "the fits repeat ",
      length(dup),
      " column name",
      .s(length(dup)),
      " [",
      fmt_trunc(dup, 8),
      "]. Join each column once",
      call. = FALSE
    )
  }
  part_pl <- lapply(parts, `[[`, "placement")
  if (!all(vapply(part_pl, is.null, TRUE))) {
    if (!is.null(placement)) {
      stop("internal error: a part of the join carries a placement of its own")
    }
    w <- vapply(parts, function(p) length(p$colnames), 0L)
    off <- c(0L, cumsum(w))
    placement <- unlist(
      lapply(seq_along(parts), function(i) {
        off[i] + if (is.null(part_pl[[i]])) seq_len(w[i]) else part_pl[[i]]
      }),
      use.names = FALSE
    )
  }
  .si_mi_new(
    chunks = unlist(lapply(parts, `[[`, "chunks"), recursive = FALSE),
    rownames = a$rownames,
    nkeep = a$nkeep,
    call = call,
    placement = placement
  )
}

#' @export
`[.si_mi` <- function(x, ...) {
  stop(
    "[ does not work on an \"si_mi\" fit. Choose columns at fit time with ",
    "'subset', join fits with cbind(), or subset a completed matrix: ",
    "as.matrix(fit, draw = b)[, j]",
    call. = FALSE
  )
}


#' @export
as.matrix.si_mi <- function(x, draw = 1L, ...) {
  .chk_no_dots(..., .what = "as.matrix()")
  draw <- .chk_count(draw, "draw", hi = x$nkeep)
  n <- length(x$rownames)
  nc <- length(x$chunks)
  whole <- nc == 1L && is.null(x$placement)
  out <- if (whole) {
    x$chunks[[1L]]$obj
  } else {
    matrix(NA_real_, n, x$coff[nc + 1L])
  }
  for (i in seq_len(nc)) {
    ch <- x$chunks[[i]]
    to <- .chunk_to(x, i)
    if (nrow(ch$obj) != n || ncol(ch$obj) != length(to)) {
      stop("internal error: a chunk does not fit its joined columns")
    }
    if (!whole) {
      out[, to] <- ch$obj
    }
    idx <- .map_lin(ch$map, n, to)
    if (!all(is.na(out[idx]))) {
      stop("internal error: a stored cell is not a hole")
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
    "as.matrix() needs a finalized fit. Call mipca_finalize() first. To look ",
    "at one completion of an unfinished fit, use inspect_completion(), but ",
    "not for analysis",
    call. = FALSE
  )
}

#' @export
mipca_complete <- function(fit, m = 20L) {
  if (inherits(fit, "si_mi_pca")) {
    .life_check(fit)
    stop(
      "mipca_complete() takes a finalized fit. Call mipca_finalize() first",
      call. = FALSE
    )
  }
  if (!inherits(fit, "si_mi")) {
    stop(
      "mipca_complete() takes a fit from mipca() or mipca_finalize()",
      call. = FALSE
    )
  }
  m <- .chk_count(m, "m")
  if (m > fit$nkeep) {
    stop(
      "'m' (",
      m,
      ") is more than the fit's ",
      fit$nkeep,
      " kept draws. ",
      "Use an 'm' of at most ",
      fit$nkeep,
      call. = FALSE
    )
  }
  draw <- as.integer(round(seq(1, fit$nkeep, length.out = m)))
  if (
    draw[1L] != 1L || draw[m] > fit$nkeep || is.unsorted(draw, strictly = TRUE)
  ) {
    stop(
      "internal error: the spaced draws are not strictly increasing in ",
      "[1, nkeep]"
    )
  }
  out <- lapply(draw, function(b) as.matrix(fit, draw = b))
  attr(out, "draw") <- draw
  out
}

#' @export
inspect_completion <- function(x, chain = 1L, draw = 1L) {
  if (!inherits(x, "si_mi_pca")) {
    stop(
      "inspect_completion() takes an unfinished fit from mipca() with ",
      "'keep_state = TRUE'. For a finalized fit, use as.matrix()",
      call. = FALSE
    )
  }
  .life_check(x)
  chain <- .chk_count(chain, "chain", hi = length(x$draws))
  draw <- .chk_count(draw, "draw", hi = ncol(x$draws[[1L]]))
  out <- x$obj
  idx <- .map_lin(x$map, nrow(out))
  if (!all(is.na(out[idx]))) {
    stop("internal error: a stored cell is not a hole")
  }
  out[idx] <- .backward_at(x$transform, x$draws[[chain]][, draw], x$map)
  out
}

#' @export
print.si_mi_pca <- function(x, ...) {
  .chk_no_dots(..., .what = "print()")
  if (!identical(x$lifecycle, "ok")) {
    cat(
      "si_mi_pca fit, ",
      if (identical(x$lifecycle, "failed")) {
        "FAILED partway through an update(). Fit it again with mipca() or mipca_fit()."
      } else if (identical(x$lifecycle, "passed")) {
        "PASSED ON to mipca_finalize(). Use the si_mi that it returned."
      } else {
        "in an unrecognized lifecycle state. Fit it again with mipca() or mipca_fit()."
      },
      "\n",
      sep = ""
    )
    return(invisible(x))
  }
  cat(sprintf(
    "si_mi_pca fit, resumable: %d x %d, rank %d, %d chain%s x %d draws\n",
    nrow(x$obj),
    ncol(x$obj),
    x$prov$rank,
    length(x$draws),
    .s(length(x$draws)),
    ncol(x$draws[[1L]])
  ))
  if (!is.null(x$abb)) {
    cat(sprintf(
      "  abb:  %d column%s left out of the fit, drawn by ABB at mipca_finalize()\n",
      ncol(x$abb$obj),
      .s(ncol(x$abb$obj))
    ))
  }
  cat("  next: update() to extend, or mipca_finalize() to finish.\n")
  cat("        summary() shows Rhat and ESS per cell.\n")
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

#' @export
print.si_mi <- function(x, ...) {
  .chk_no_dots(..., .what = "print()")
  h <- sum(vapply(x$chunks, function(ch) nrow(ch$map), 0L))
  cat(sprintf(
    "si_mi fit: %d x %d, %d imputed cell%s, %d draw%s\n",
    length(x$rownames),
    length(x$colnames),
    h,
    .s(h),
    x$nkeep,
    .s(x$nkeep)
  ))
  nleft <- sum(vapply(
    x$chunks,
    function(ch) if (isTRUE(ch$prov$left_na)) ncol(ch$obj) else 0L,
    0L
  ))
  if (nleft > 0L) {
    cat(sprintf(
      "  left NA: %d column%s not imputed; as.matrix() returns %s holes as NA.\n",
      nleft,
      .s(nleft),
      if (nleft == 1L) "its" else "their"
    ))
  }
  ncell <- sum(vapply(x$chunks, function(ch) nrow(ch$diag$cells), 0L))
  nxt <- if (x$nkeep == 1L) {
    "as.matrix() for the completed matrix"
  } else if (x$nkeep < 20L) {
    sprintf(
      "mipca_complete() with `m` at most %d for completed matrices",
      x$nkeep
    )
  } else {
    "mipca_complete() for completed matrices"
  }
  cat(
    "  next: ",
    nxt,
    if (ncell > 0L) ", summary() for Rhat and ESS",
    "\n",
    sep = ""
  )
  invisible(x)
}


.diag_sample <- function(map, p, prop_diag) {
  cnt <- tabulate(map[, "col"], nbins = p)
  off <- cumsum(c(0L, cnt))
  if (off[p + 1L] != nrow(map)) {
    stop("internal error: map column runs do not tile the map")
  }
  dirty <- which(cnt > 0L)
  pick <- vector("list", length(dirty))
  for (i in seq_along(dirty)) {
    j <- dirty[i]
    r <- off[j] + seq_len(cnt[j])
    h <- length(r)
    k <- min(h, max(1L, as.integer(ceiling(prop_diag * h))))
    pick[[i]] <- r[floor(seq(1, h, length.out = k))]
  }
  ns <- lengths(pick)
  rows <- unlist(pick, use.names = FALSE)
  owner <- rep.int(seq_along(dirty), ns)
  if (length(rows) > 1L && any(diff(rows) <= 0L)) {
    stop("internal error: the diagnostic sample is not ascending")
  }
  if (!identical(unname(map[rows, "col"]), dirty[owner])) {
    stop("internal error: a sampled cell is outside the column it came from")
  }
  list(rows = rows)
}

.cell_diag_read <- function(draws, srows) {
  res <- mipca_cell_diag(draws, srows, 256L)

  if (ncol(draws[[1L]]) %/% 2L <= 5L) {
    res[, 2L] <- NA_real_
  }
  res
}

.diag_frame <- function(map, draws, p, nms, prop_diag) {
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  if (is.null(nms) || length(nms) != p) {
    stop("internal error: the diagnostic frame needs one name per column")
  }
  S <- ncol(draws[[1L]])
  if (nrow(map) == 0L) {
    srows <- integer(0)
    res <- matrix(numeric(0), 0L, 2L)
  } else {
    srows <- .diag_sample(map, p, prop_diag)$rows
    res <- .cell_diag_read(draws, srows)
  }
  cs <- map[srows, "col"]
  cells <- structure(
    list(
      col = unname(cs),
      name = nms[cs],
      row = unname(map[srows, "row"]),
      rhat = unname(res[, 1L]),
      ess = unname(res[, 2L])
    ),
    class = "data.frame",
    row.names = .set_row_names(length(srows))
  )

  structure(
    list(cells = cells),
    prop_diag = prop_diag,
    nchains = length(draws),
    ndraws = S,
    class = "si_mi_summary"
  )
}

#' @export
summary.si_mi_pca <- function(object, prop_diag = 0.5, ...) {
  .life_check(object)
  .chk_no_dots(..., .what = "summary()")
  .diag_frame(
    map = object$map,
    draws = object$draws,
    p = ncol(object$obj),
    nms = colnames(object$obj),
    prop_diag = prop_diag
  )
}

#' @export
summary.si_mi <- function(object, ...) {
  .chk_no_dots(..., .what = "summary()")
  cells <- do.call(
    rbind,
    lapply(seq_along(object$chunks), function(i) {
      d <- object$chunks[[i]]$diag$cells
      d$col <- .to_col(.chunk_to(object, i), d$col)
      cbind(chunk = rep.int(i, nrow(d)), d)
    })
  )
  agree <- function(nm) {
    .agree(unlist(lapply(object$chunks, function(ch) attr(ch$diag, nm))))
  }
  structure(
    list(cells = cells),
    prop_diag = agree("prop_diag"),
    nchains = agree("nchains"),
    ndraws = agree("ndraws"),
    class = "si_mi_summary"
  )
}

#' @export
print.si_mi_summary <- function(x, n = 6L, ...) {
  .chk_no_dots(..., .what = "print()")
  n <- .chk_count(n, "n", lo = 0L)
  cells <- x$cells
  num <- function(nm) .fmt_agree(attr(x, nm))
  ncol_aff <- length(unique(cells$name))
  cat(sprintf(
    "si_mi_summary: Rhat and ESS for %d cell%s in %d column%s\n",
    nrow(cells),
    .s(nrow(cells)),
    ncol_aff,
    .s(ncol_aff)
  ))
  cat(sprintf(
    "  prop_diag %s, %s chain%s x %s draw%s\n",
    num("prop_diag"),
    num("nchains"),
    .s(attr(x, "nchains")),
    num("ndraws"),
    .s(attr(x, "ndraws"))
  ))
  q <- function(v) {
    stats::quantile(v, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  }
  print(signif(rbind(rhat = q(cells$rhat), ess = q(cells$ess)), 3L))
  bad <- order(cells$rhat, decreasing = TRUE, na.last = NA)
  if (length(bad) > 0L && n > 0L) {
    k <- min(n, length(bad))
    cat(sprintf("worst %d cell%s by rhat:\n", k, .s(k)))
    print(cells[bad[seq_len(k)], , drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}

.draw_parts <- function(fit) {
  if (inherits(fit, "si_mi_pca")) {
    .life_check(fit)
    return(lapply(fit$draws, function(m) {
      list(
        d = m,
        tf = fit$transform,
        col = fit$map[, "col"],
        row = fit$map[, "row"]
      )
    }))
  }
  if (inherits(fit, "si_mi")) {
    return(lapply(fit$chunks, function(ch) {
      list(
        d = ch$draws,
        tf = ch$transform,
        col = ch$map[, "col"],
        row = ch$map[, "row"]
      )
    }))
  }
  stop("'fit' must be a fit from mipca() or mipca_finalize()", call. = FALSE)
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
