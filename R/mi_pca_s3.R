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
    stop("internal error: a chunk's columns and their names do not line up.")
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
        "internal error: the placement vector is not a permutation of the joined columns."
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
    stop("internal error: a stored cell's column is outside its chunk.")
  }
  to[col]
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

#' @export
mipca_clone <- function(x, ...) {
  UseMethod("mipca_clone")
}

#' @export
mipca_clone.si_mi_pca <- function(x, ...) {
  .life_check(x)
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
  stop(
    "this object's lifecycle field is not one of ok / failed / passed.",
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
    stop(
      "resumed map disagrees with the fit's map; resume refused.",
      call. = FALSE
    )
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
    "this fit is already finalized; mipca_finalize() is a one-way transition.",
    call. = FALSE
  )
}

#' @export
update.si_mi <- function(object, ...) {
  states <- vapply(object$chunks, function(ch) ch$prov$seam_state, "")
  if (all(states == "never")) {
    stop(
      "this fit was built with keep_state = FALSE, so no seam was ever ",
      "kept; rerun with keep_state = TRUE to get a resumable fit.",
      call. = FALSE
    )
  }
  te <- .prov_agree(object$chunks, "thin_effective")
  if (is.na(te)) {
    stop(
      "this fit was finalized, which drops the seam; its chunks were thinned ",
      "differently, so there is no single spacing to re-run at -- see ",
      "summary() for the per-chunk numbers.",
      call. = FALSE
    )
  }
  stop(
    "this fit was finalized, which drops the seam; re-run with thin = ",
    te,
    " to get this spacing directly.",
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
          "backward: column index outside the finalized chunk.",
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
  inv <- integer(ncol(x$obj))
  inv[x$subset] <- seq_along(x$subset)
  loc <- inv[x$map[, "col"]]
  fr$cells$col <- inv[fr$cells$col]
  if (min(c(loc, fr$cells$col, 1L)) < 1L) {
    stop("internal error: a stored hole maps outside the requested subset.")
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
mipca_finalize.si_mi_pca <- function(
  x,
  thin = 1L,
  prop_diag = .prop_diag_default(),
  num_threads = getOption("mc.cores", 1L),
  ...
) {
  .life_check(x)
  .chk_no_dots(..., .what = "mipca_finalize()")
  k <- .chk_count(thin, "thin")
  nd <- ncol(x$draws[[1L]])
  if (k > nd) {
    stop(
      "thin (",
      k,
      ") exceeds the ",
      nd,
      " stored draws per chain; ",
      "there would be nothing left to keep.",
      call. = FALSE
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
    prop_diag = prop_diag,
    num_threads = num_threads
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

  chunk <- list(
    obj = rb$obj,
    map = if (is.null(rb$loc)) {
      x$map
    } else {
      cbind(row = x$map[, "row"], col = rb$loc)
    },
    draws = if (length(dl) == 1L) dl[[1L]] else do.call(cbind, dl),
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
cbind.si_mi <- function(..., deparse.level = 0) {
  parts <- Filter(Negate(is.null), list(...))
  if (!all(vapply(parts, inherits, TRUE, "si_mi"))) {
    stop(
      "cbind() over a finalized fit takes si_mi objects only; a resumable ",
      "si_mi_pca must be finalized first.",
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
        "chunks disagree on rownames; no reordering or intersection is done.",
        call. = FALSE
      )
    }
    if (!identical(a$nkeep, b$nkeep)) {
      stop(
        "chunks disagree on the number of kept draws (",
        a$nkeep,
        " vs ",
        b$nkeep,
        "); min-trimming would falsify the frozen diagnostics.",
        call. = FALSE
      )
    }
  }
  nms <- unlist(lapply(parts, `[[`, "colnames"), use.names = FALSE)
  if (anyDuplicated(nms) > 0L) {
    dup <- unique(nms[duplicated(nms)])
    stop(
      "chunks repeat ",
      length(dup),
      " column name",
      .s(length(dup)),
      " [",
      fmt_trunc(dup, 8),
      "]; a column may be joined once.",
      call. = FALSE
    )
  }
  part_pl <- lapply(parts, `[[`, "placement")
  if (!all(vapply(part_pl, is.null, TRUE))) {
    if (!is.null(placement)) {
      stop("internal error: a part of the join carries a placement of its own.")
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
    "`[` is not defined for a finalized fit (si_mi). Choose columns at fit ",
    "time with `subset =`, join fits with cbind(), or complete first and ",
    "subset the matrix: as.matrix(fit, draw = b)[, j].",
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
      stop("internal error: a chunk does not fit its joined columns.")
    }
    if (!whole) {
      out[, to] <- ch$obj
    }
    idx <- .map_lin(ch$map, n, to)
    if (!all(is.na(out[idx]))) {
      stop("internal error: a stored cell is not a hole.")
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
    "as.matrix() is refused on a resumable fit: it is not finalized (chains ",
    "not pooled, diagnostics not frozen). Call mipca_finalize() first, or ",
    "inspect_completion() to eyeball one mid-run completion (not for analysis).",
    call. = FALSE
  )
}

#' @export
mipca_complete <- function(fit, m = 20L) {
  if (inherits(fit, "si_mi_pca")) {
    .life_check(fit)
    stop(
      "mipca_complete() takes a finalized fit; call mipca_finalize() first.",
      call. = FALSE
    )
  }
  if (!inherits(fit, "si_mi")) {
    stop(
      "mipca_complete() takes a fit from mipca() or mipca_finalize().",
      call. = FALSE
    )
  }
  m <- .chk_count(m, "m")
  if (m > fit$nkeep) {
    stop(
      "m = ",
      m,
      " is more completions than the fit's ",
      fit$nkeep,
      " kept draws.",
      call. = FALSE
    )
  }
  draw <- as.integer(round(seq(1, fit$nkeep, length.out = m)))
  if (
    draw[1L] != 1L || draw[m] > fit$nkeep || is.unsorted(draw, strictly = TRUE)
  ) {
    stop(
      "internal error: the spaced draws are not strictly increasing in ",
      "[1, nkeep]."
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
      "inspect_completion() is for a resumable si_mi_pca fit.",
      call. = FALSE
    )
  }
  .life_check(x)
  chain <- .chk_count(chain, "chain", hi = length(x$draws))
  draw <- .chk_count(draw, "draw", hi = ncol(x$draws[[1L]]))
  out <- x$obj
  idx <- .map_lin(x$map, nrow(out))
  if (!all(is.na(out[idx]))) {
    stop("internal error: a stored cell is not a hole.")
  }
  out[idx] <- .backward_at(x$transform, x$draws[[chain]][, draw], x$map)
  out
}

#' @export
print.si_mi_pca <- function(x, ...) {
  .chk_no_dots(..., .what = "print()")
  if (!identical(x$lifecycle, "ok")) {
    cat(sprintf(
      "si_mi_pca fit, %s\n  every method but print() stops with an error\n",
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
  cat(sprintf("si_mi_pca fit, resumable, %s transform\n", x$transform$name))
  cat(sprintf("  obj:     %d x %d matrix\n", nrow(x$obj), ncol(x$obj)))
  cat(sprintf(
    "  draws:   %d chain%s x %d stored holes x %d draws, on the latent scale\n",
    length(x$draws),
    .s(length(x$draws)),
    nrow(x$map),
    ncol(x$draws[[1L]])
  ))
  cat(sprintf(
    "  columns: %d requested, %d with stored holes\n",
    length(x$subset),
    length(x$cols)
  ))
  if (!is.null(x$abb)) {
    cat(sprintf(
      "  abb:     %d column%s listed, %d requested; not in obj, drawn at mipca_finalize()\n",
      x$abb$nlisted,
      .s(x$abb$nlisted),
      ncol(x$abb$obj)
    ))
  }
  if (!is.null(x$seam$store)) {
    cat(sprintf(
      "  store:   cell-level, draws kept for %d of %d holes\n",
      nrow(x$map),
      x$prov$nmiss
    ))
  }
  cat(sprintf(
    "  sampled: rank %d, warmup %d, thin %d, %d iterations so far\n",
    x$prov$rank,
    x$prov$warmup,
    x$prov$thin,
    x$prov$iter
  ))
  if (!is.null(x$prov$tuned)) {
    cat(x$prov$tuned$line)
  }
  cat(
    "  next:    update(fit, iter_sampling) to extend, mipca_finalize(fit) to finish\n"
  )
  cat(
    "           summary(fit) for per-cell Rhat and ESS, to decide between them\n"
  )
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
      "  left NA: %d column%s not imputed, holes still NA; the flank plan's $targets$status says why\n",
      nleft,
      .s(nleft)
    ))
  }
  nabb <- sum(vapply(
    x$chunks,
    function(ch) if (isTRUE(ch$prov$abb)) ncol(ch$obj) else 0L,
    0L
  ))
  if (nabb > 0L) {
    cat(sprintf(
      "  abb:     %d column%s drawn by ABB from %s own observed cells\n",
      nabb,
      .s(nabb),
      if (nabb == 1L) "its" else "their"
    ))
  }
  cat(sprintf(
    "  draws:   %d kept, pooled from %s chain%s; thin %s at finalize, %s overall\n",
    x$nkeep,
    .fmt_agree(nch),
    .s(nch),
    .fmt_agree(.prov_agree(x$chunks, "thin_final")),
    .fmt_agree(.prov_agree(x$chunks, "thin_effective"))
  ))
  for (ch in x$chunks) {
    if (!is.null(ch$prov$tuned)) {
      cat(ch$prov$tuned$line)
    }
  }
  ncell <- sum(vapply(x$chunks, function(ch) nrow(ch$diag$cells), 0L))
  cat(sprintf(
    "  diag:    Rhat and ESS for %d cell%s, frozen at finalize; summary(fit) shows them\n",
    ncell,
    .s(ncell)
  ))
  invisible(x)
}

.prop_diag_default <- function() 0.5

.diag_sample <- function(map, p, prop_diag) {
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
    h <- length(r)
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

.diag_frame <- function(map, draws, p, nms, prop_diag, num_threads) {
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  num_threads <- .chk_count(num_threads, "num_threads")
  if (is.null(nms) || length(nms) != p) {
    stop("internal error: the diagnostic frame needs one name per column.")
  }
  srows <- .diag_sample(map, p, prop_diag)$rows

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

  structure(
    list(cells = cells),
    prop_diag = prop_diag,
    nchains = length(draws),
    ndraws = S,
    class = "si_mi_summary"
  )
}

#' @export
summary.si_mi_pca <- function(
  object,
  prop_diag = .prop_diag_default(),
  num_threads = getOption("mc.cores", 1L),
  ...
) {
  .life_check(object)
  .chk_no_dots(..., .what = "summary()")
  .diag_frame(
    map = object$map,
    draws = object$draws,
    p = ncol(object$obj),
    nms = colnames(object$obj),
    prop_diag = prop_diag,
    num_threads = num_threads
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
  pd <- attr(x, "prop_diag")
  how <- if (is.na(pd)) {
    "prop_diag differs by chunk"
  } else if (pd >= 1) {
    "every hole"
  } else {
    sprintf("%g%% of each column's holes", 100 * pd)
  }
  cat(sprintf(
    "si_mi_summary: Rhat and ESS for %d cell%s (%s) over %d column%s\n",
    nrow(cells),
    .s(nrow(cells)),
    how,
    ncol_aff,
    .s(ncol_aff)
  ))
  cat(sprintf(
    "  from %s chains x %s draws\n",
    num("nchains"),
    num("ndraws")
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
  bad <- order(cells$rhat, decreasing = TRUE, na.last = NA)
  if (length(bad) > 0L && n > 0L) {
    k <- min(n, length(bad))
    cat(sprintf("  worst %d cell%s by rhat:\n", k, .s(k)))
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
  stop("not a si_mi_pca or si_mi fit.", call. = FALSE)
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
