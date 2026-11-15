# mi_pca_s3.R ---------------------------------------------------------------
# The fit object (Layer 3), its exits (as.matrix(), mipca_complete(),
# inspect_completion(), clamp_frac()) and the per-cell diagnostic frame
# behind summary() (Layer 3b). The sampler, the transforms and the entry
# point mipca() are in R/mi_pca_imp.R.
# ------------------------------------------------------------------------------

# ---- Layer 3: the fit object, in two stages ----------------------------------
# si_mi_pca (resumable, an environment) -> mipca_finalize() -> si_mi (a value
# that cbind() joins). Disjoint classes, no inheritance: a subclass would let
# cbind.si_mi dispatch on a resumable fit. Mutation exists exactly where
# resume exists.

# The resumable stage's one constructor. `map` is never rewritten: update()'s
# wrong-matrix refusal compares it identical(). `parent = emptyenv()` is
# required, or the matrix in this frame rides into every saved fit.
#
# With a list in play (`abb`, .abb_split()), `obj`, `map`, `cols` and
# `subset` are the KEPT columns in their own numbering -- the matrix the
# sampler saw, so update() resumes against it unchanged -- and `abb` holds
# the requested listed columns, their keys, their seed and the placement
# that returns the input's column order, until mipca_finalize() draws them.
# NULL when no column was listed.
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

# The one owner of a stored `$call`. The head is rewritten to `fn` (a method
# or do.call() would record the method or the whole closure), and an argument
# that arrived as a VALUE (do.call() inlines `obj`) becomes a placeholder
# naming its class, so the record never carries the matrix.
.scrub_call <- function(cl, fn) {
  ## Over a list, not the call: `cl[[i]]` on a call walks the pairlist, so
  ## do.call(cbind, fits) over N parts would cost O(N^2) here.
  args <- as.list(cl)
  args[[1L]] <- as.name(fn)
  for (i in seq_along(args)[-1L]) {
    a <- args[[i]]
    ## Language and NULL pass through untouched: those ARE expressions.
    if (is.symbol(a) || is.call(a) || is.null(a)) {
      next
    }
    ## A bare scalar is cheap and informative; anything longer, or anything
    ## carrying a class, is the inlined-value case.
    if (length(a) <= 1L && !is.object(a)) {
      next
    }
    args[[i]] <- paste0("<", class(a)[1L], ">")
  }
  as.call(args)
}

# The finalized stage's one constructor, cbind() included. `rownames` and
# `nkeep` are hoisted out of the chunks; `colnames` and `coff` are derived
# here, chunk i holding joined columns coff[i] + 1 .. coff[i + 1].
#
# THE PLACEMENT VECTOR. `placement[k]` is the OUTPUT column of joined column
# k, the column as.matrix() writes it to and summary() numbers it by. NULL is
# the identity: chunk order is then the column order, as it is on every path
# but the plain door's join with an `abb` list (.abb_join()), whose two
# chunks interleave in the input's order. `colnames` is in OUTPUT order, what
# colnames(as.matrix(fit)) returns. An identity placement is stored as NULL.
# `coff` numbers the JOINED columns and `colnames` the OUTPUT ones; the one
# reader that crosses from the first to the second is .chunk_to().
.si_mi_new <- function(chunks, rownames, nkeep, call, placement = NULL) {
  width <- vapply(chunks, function(ch) ncol(ch$obj), 0L)
  nms <- unlist(
    lapply(chunks, function(ch) colnames(ch$obj)),
    use.names = FALSE
  )
  coff <- c(0L, cumsum(width))
  ## The closed form's own guard: every chunk has columns, every column has a
  ## name, and the last offset is the joined width.
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
    ## A permutation of the joined columns, in O(p): one slot each, each slot
    ## once. A wrong formula stops here instead of returning plausible
    ## columns under the wrong names.
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
    ## An ascending permutation is the identity.
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

# THE OUTPUT COLUMNS of chunk i: its column k is output column to[k]. One
# owner of the step from a chunk's own numbering to the result's:
# as.matrix() sends the columns and the holes through this table and
# summary() its `col`.
.chunk_to <- function(x, i) {
  k <- x$coff[i] + seq_len(x$coff[i + 1L] - x$coff[i])
  if (is.null(x$placement)) k else x$placement[k]
}

# A chunk's own columns `col` read through its table `to` (.chunk_to()),
# guarded first: R's indexing drops a 0 and returns NA past the end, both
# without a word, and a dropped index recycles the rest onto other cells.
.to_col <- function(to, col) {
  if (
    length(col) > 0L &&
      (anyNA(col) || min(col) < 1L || max(col) > length(to))
  ) {
    stop("internal error: a stored cell's column is outside its chunk.")
  }
  to[col]
}

# Rows of the per-chain draws, pooled across chains; one caller, the tuner.
.pool_rows <- function(draws, rows) {
  first <- draws[[1L]][rows, , drop = FALSE]
  if (length(draws) == 1L) {
    return(first)
  }
  ## Preallocate and fill: one chain's subset at a time, where pooling then
  ## subsetting, or subsetting then cbind(), holds the whole block twice. The
  ## fill drops dimnames, so none may exist.
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

# COPY: `b <- fit` is a second name for one environment. Shallow suffices:
# every field is an ordinary copy-on-modify value.
#' @export
mipca_clone <- function(x, ...) {
  UseMethod("mipca_clone")
}

#' @export
mipca_clone.si_mi_pca <- function(x, ...) {
  ## A copy of a failed or passed fit would be refused everywhere anyway;
  ## refused here so the door says so first.
  .life_check(x)
  .chk_no_dots(..., .what = "mipca_clone()")
  ## as.list.environment() explicitly, not as.list(): `x` carries a class
  ## attribute, so plain as.list() would dispatch.
  e <- list2env(
    as.list.environment(x, all.names = TRUE),
    parent = emptyenv()
  )
  class(e) <- "si_mi_pca"
  e
}

# The lifecycle gate. One owner of the refusal, so a failed or passed fit
# gets the same message from every door. print() does not call it (see
# print.si_mi_pca()).
.life_check <- function(x) {
  st <- x$lifecycle
  if (identical(st, "ok")) {
    return(invisible(x))
  }
  if (identical(st, "failed")) {
    stop(
      "This fit is marked failed: an update() call stopped partway through ",
      "it, so its chains may hold different numbers of draws. Fit again ",
      "with mipca().",
      call. = FALSE
    )
  }
  if (identical(st, "passed")) {
    stop(
      "This fit was finalized: mipca_finalize() moved its draws to the ",
      "`si_mi` object it returned. Use that object.",
      call. = FALSE
    )
  }
  ## A hand-built environment, or a state with no branch. Loud, because
  ## silence here would read as permission.
  stop(
    "This object is not a fit from mipca(): its `lifecycle` field is not ",
    "\"ok\", \"failed\" or \"passed\".",
    call. = FALSE
  )
}

# EXTEND every chain as if it never stopped, from the seam; the resume is
# checked on a tolerance plus a fingerprint. MUTATES `object` in place. No
# `keep_state`: mipca_finalize() is how a fit drops its seam.
#' @export
update.si_mi_pca <- function(
  object,
  iter_sampling,
  ...
) {
  .life_check(object)
  seam <- object$seam
  .seam_check_version(seam)

  ## Refused, not forwarded. For the docs: rank, thin, spectra and subset come
  ## from the seam, and `init` describes a start these chains already have.
  .chk_no_dots(..., .what = "update()")

  ## `subset` is the requested one, not the narrowed one.
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
  ## The map is a pure function of (NA pattern, subset), and the core has
  ## already matched the NA pattern against the seam.
  if (!identical(core$map, object$map)) {
    stop("internal error: the resumed map disagrees with the fit's map.")
  }

  ## EVERYTHING BELOW MUTATES IN PLACE: marked `failed` before the first
  ## write and cleared after the last, so a throw in between stays marked.
  object$lifecycle <- "failed"
  ## Chain at a time, NOT Map(cbind, ...), which holds every chain's old, new
  ## and joined panels at once.
  newd <- core$draws
  core$draws <- NULL
  for (ch in seq_along(object$draws)) {
    object$draws[[ch]] <- .alloc_draws(
      cbind(object$draws[[ch]], newd[[ch]]),
      cells = nrow(newd[[ch]]),
      chains = length(object$draws),
      draws = ncol(object$draws[[ch]]) + ncol(newd[[ch]]),
      fix = "Ask update() for fewer `iter_sampling`."
    )
    ## `newd[ch] <- list(NULL)`, never `newd[[ch]] <- NULL`: the latter
    ## DELETES the element and renumbers every chain after it.
    newd[ch] <- list(NULL)
  }
  rm(newd)
  object$seam <- core$state
  object$prov$ndraws <- ncol(object$draws[[1L]])
  object$prov$iter <- core$iter
  object$lifecycle <- "ok"
  ## Invisibly, after rjags::update().
  invisible(object)
}

# FINISH: trim, pool, drop the seam and narrow `obj`, one transition.
# Returns a one-chunk si_mi.
#' @export
mipca_finalize <- function(x, ...) {
  UseMethod("mipca_finalize")
}

#' @export
mipca_finalize.si_mi <- function(x, ...) {
  stop(
    "This fit is already finalized. Pass it to mipca_complete() or ",
    "summary().",
    call. = FALSE
  )
}

# Un-updatable, and the refusal must say WHICH un-updatable. This method has
# to EXIST rather than dispatch to stats::update.default, which would find
# `$call` on this object and silently re-run the whole fit.
#' @export
update.si_mi <- function(object, ...) {
  ## Read across ALL chunks, not off chunk 1: after cbind() the two recoveries
  ## can coexist.
  states <- vapply(object$chunks, function(ch) ch$prov$seam_state, "")
  if (all(states == "never")) {
    stop(
      "update() cannot extend this fit: it was made with ",
      "`keep_state = FALSE`. Fit again with mipca() and `keep_state = TRUE`.",
      call. = FALSE
    )
  }
  ## The parts of a joined fit can be thinned differently, and then no one
  ## `thin` gives this fit's spacing of draws.
  te <- .prov_agree(object$chunks, "thin_effective")
  stop(
    "update() cannot extend a finalized fit. Fit again with mipca() and ",
    "`keep_state = TRUE`",
    if (is.na(te) || te == 1) {
      "."
    } else {
      paste0(", with `thin = ", te, "` to keep this fit's thinning.")
    },
    call. = FALSE
  )
}

# The transform in the CHUNK's column system: untranslated, a rebased index
# passes the range check and mis-scales every draw silently. Backward-only,
# so a stray forward fails loudly. force() both: an unforced promise keeps
# .rebase_subset()'s frame, and with it the whole matrix, in every chunk.
.tf_narrow <- function(tf, subset) {
  force(tf)
  force(subset)
  list(
    name = tf$name,
    backward = function(z, col, row = NULL) {
      ## Past the end of `subset` the read returns NA, silently.
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

# Narrow `obj` to the REQUESTED subset (not the map's columns) and re-base the
# map and the diagnostics through ONE scatter table, so the two cannot
# disagree. The identity selection is skipped: it would be a full copy.
.rebase_subset <- function(x, fr) {
  ## NULL `loc` means "unchanged", so the caller hands the map through.
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
  ## Row names dropped: they are hoisted onto si_mi. The identity branch
  ## keeps them, since stripping would cost the copy it avoids.
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
  ...
) {
  .life_check(x)
  .chk_no_dots(..., .what = "mipca_finalize()")
  k <- .chk_count(thin, "thin")
  nd <- ncol(x$draws[[1L]])
  if (k > nd) {
    stop(
      "`thin` (",
      k,
      ") is larger than the ",
      nd,
      " draws each chain holds, ",
      "so no draw would be kept. Use a `thin` of at most ",
      nd,
      ".",
      call. = FALSE
    )
  }

  ## TRIM, then POOL: spacing is a within-chain property, so a stride taken
  ## after pooling would step across chain boundaries. k == 1 is skipped.
  dl <- if (k == 1L) {
    x$draws
  } else {
    keep <- seq.int(1L, nd, by = k)
    lapply(x$draws, function(m) m[, keep, drop = FALSE])
  }
  ## NOT pooled yet: the frame below needs the chains axis, and the pooled
  ## matrix is built once, at the chunk.
  nkeep <- sum(vapply(dl, ncol, 0L))

  ## FREEZE the diagnostic frame over the trimmed draws that get stored --
  ## the last moment it can be computed, since Rhat is between-chain.
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
  ## 200 draws that came from 2000 is a different object from 200 that came
  ## from 200 -- the effective thin records which.
  prov$thin_effective <- prov$thin * k
  prov$seam_state <- if (identical(prov$seam_state, "never")) {
    "never"
  } else {
    "dropped"
  }

  chunk <- list(
    obj = rb$obj,
    ## No column moved: `x$map` as it stands.
    map = if (is.null(rb$loc)) {
      x$map
    } else {
      cbind(row = x$map[, "row"], col = rb$loc)
    },
    ## Same identity: one chain needs no cbind at all.
    draws = if (length(dl) == 1L) {
      dl[[1L]]
    } else {
      .alloc_draws(
        do.call(cbind, dl),
        cells = nrow(dl[[1L]]),
        chains = length(dl),
        draws = ncol(dl[[1L]]),
        fix = "Pass a larger `thin` to mipca_finalize(), or fit with fewer `chains` or `iter_sampling`."
      )
    },
    transform = rb$transform,
    diag = rb$diag,
    prov = prov,
    call = x$call
  )
  out <- .si_mi_new(
    chunks = list(chunk),
    ## Hoisted onto si_mi: what cbind matches on.
    rownames = rownames(x$obj),
    nkeep = nkeep,
    call = x$call
  )

  ## THE LISTED SIDE, drawn here and not beside the sampler: the kept-draw
  ## count is fixed only now, after update() extended and `thin` trimmed.
  ## The seed was fixed at fit time.
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

  ## PASSED ON, and the draws go with it, so what is left holds them once.
  ## The PEAK is 2x, not 1x: the pooled matrix above is allocated while every
  ## chain still exists, which any pooling into one matrix must do (thin > 1
  ## adds the trimmed copy). Under 2x needs chunks that keep the chains
  ## apart, not a different fill. Marked LAST, after `out` is built, so a
  ## throw above leaves a usable fit.
  x$draws <- NULL
  x$seam <- NULL
  x$lifecycle <- "passed"
  out
}

# JOIN CHUNKS along the variable axis: pure concatenation, no alignment, no
# min-trimming, each column name once.
# For the docs: fits made apart and joined here take different seeds. On one
# seed, two sampled parts that agree on rows, width and hole count share
# their noise, and two listed columns at one column index with one donor
# count draw the same donors. Not refused: on any other shapes one seed
# shares nothing, and each chunk's prov records the seed it started from.
#' @export
cbind.si_mi <- function(..., deparse.level = 0) {
  parts <- Filter(Negate(is.null), list(...))
  if (!all(vapply(parts, inherits, TRUE, "si_mi"))) {
    stop(
      "cbind() joins finalized fits (`si_mi`) only. Call mipca_finalize() ",
      "on an `si_mi_pca` fit first.",
      call. = FALSE
    )
  }
  if (length(parts) == 1L) {
    return(parts[[1L]])
  }
  ## Here, not in the argument: match.call() reads this frame.
  call <- .scrub_call(match.call(), "cbind")
  .si_mi_join(parts, call = call)
}

# THE JOIN, one owner for cbind.si_mi() and the ABB join (.abb_join()): the
# checks, the placement, one construction. `placement`, when given, places
# the parts' columns as concatenated in argument order. It is the caller's
# derivation over parts that hold their own columns in order, so a part that
# already moved a column is refused, not composed with.
.si_mi_join <- function(parts, call, placement = NULL) {
  a <- parts[[1L]]
  for (b in parts[-1L]) {
    ## Not reordered or intersected: the caller's rows are the caller's.
    if (!identical(a$rownames, b$rownames)) {
      stop(
        "The fits have different row names. cbind() joins fits of the same ",
        "rows in the same order.",
        call. = FALSE
      )
    }
    ## Not trimmed to the smaller count: the frozen Rhat and ESS describe
    ## every kept draw.
    if (!identical(a$nkeep, b$nkeep)) {
      stop(
        "The fits hold different numbers of kept draws (",
        a$nkeep,
        " and ",
        b$nkeep,
        "). Fit and finalize every part with the same `chains`, ",
        "`iter_sampling` and `thin`.",
        call. = FALSE
      )
    }
  }
  ## Validation of caller input, once: one hashed pass, not a position
  ## lookup. Hint for the docs: say which two parts share the name.
  nms <- unlist(lapply(parts, `[[`, "colnames"), use.names = FALSE)
  if (anyDuplicated(nms) > 0L) {
    dup <- unique(nms[duplicated(nms)])
    stop(
      "The fits repeat ",
      length(dup),
      " column name",
      .s(length(dup)),
      " [",
      fmt_trunc(dup, 8),
      "]. Join each column once.",
      call. = FALSE
    )
  }
  ## The parts' placement vectors, composed: each part keeps its own column
  ## order and the parts come back in argument order, as base cbind() does.
  ## Part i's columns are shifted by the width before it; NULL reads as the
  ## identity. All NULL leaves `placement` as given: NULL for the flank
  ## driver's join, the record's for the ABB join.
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

# REFUSED: `fit[, j]` reads as a column subset, and base `[` answers it with
# "incorrect number of dimensions" (two indices) or a classless list (one),
# neither of which says what to do. Columns are chosen at fit time
# (`subset =`) or joined (cbind()); nothing here reads an si_mi with `[`.
# `$` and `[[` stay open: every method reads fields through them.
# Hint for the docs: show `as.matrix(fit, draw = b)[, j]`.
#' @export
`[.si_mi` <- function(x, ...) {
  stop(
    "`[` does not work on an `si_mi` fit. Choose columns at fit time with ",
    "`subset`, join fits with cbind(), or subset a completed matrix: ",
    "as.matrix(fit, draw = b)[, j].",
    call. = FALSE
  )
}

# ---- Exits: the backward map lives here, applied exactly once ----------------

# One completed dataset, with the backward map applied HERE, once.
#' @export
as.matrix.si_mi <- function(x, draw = 1L, ...) {
  .chk_no_dots(..., .what = "as.matrix()")
  draw <- .chk_count(draw, "draw", hi = x$nkeep)
  n <- length(x$rownames)
  nc <- length(x$chunks)
  ## One chunk with no column moved IS the result as stored: started from,
  ## not allocated beside. The index below is the same either way.
  whole <- nc == 1L && is.null(x$placement)
  ## Inline, not a helper taking `out`: assigning into it inside a function
  ## copies the joined matrix once per chunk.
  out <- if (whole) {
    x$chunks[[1L]]$obj
  } else {
    matrix(NA_real_, n, x$coff[nc + 1L])
  }
  for (i in seq_len(nc)) {
    ch <- x$chunks[[i]]
    ## Chunk i's column k is output column to[k]. The columns and the holes
    ## go through this ONE table on every fit, so a draw cannot land in
    ## another column than its cell did.
    to <- .chunk_to(x, i)
    if (nrow(ch$obj) != n || ncol(ch$obj) != length(to)) {
      stop("internal error: a chunk does not fit its joined columns.")
    }
    if (!whole) {
      out[, to] <- ch$obj
    }
    idx <- .map_lin(ch$map, n, to)
    ## Every stored cell is a hole of `obj` (nothing is held out). Asserted,
    ## not masked: a draw over an observed value is leakage.
    if (!all(is.na(out[idx]))) {
      stop("internal error: a stored cell is not a hole.")
    }
    out[idx] <- .backward_at(ch$transform, ch$draws[, draw], ch$map)
  }
  ## The primitive, on the private matrix: row names come from the one
  ## hoisted copy, and a narrowed chunk carries none of its own.
  dimnames(out) <- list(x$rownames, x$colnames)
  out
}

# Refused because the fit is not finalized, not on autocorrelation.
#' @export
as.matrix.si_mi_pca <- function(x, ...) {
  ## First, so a passed-on fit is pointed at the object that can answer.
  .life_check(x)
  stop(
    "as.matrix() needs a finalized fit. Call mipca_finalize() first. To look ",
    "at one completion of an unfinished fit, use inspect_completion(), but ",
    "not for analysis.",
    call. = FALSE
  )
}

# m completed datasets for Rubin's rules, at draws spaced evenly over the
# pooled chains.
#' @export
mipca_complete <- function(fit, m = 20L) {
  if (inherits(fit, "si_mi_pca")) {
    .life_check(fit)
    stop(
      "mipca_complete() takes a finalized fit. Call mipca_finalize() first.",
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
  ## Hint for the docs: more draws means update() on the live fit
  ## before mipca_finalize(), or a refit with a larger iter_sampling.
  if (m > fit$nkeep) {
    stop(
      "`m` (",
      m,
      ") is more than the fit's ",
      fit$nkeep,
      " kept draws. ",
      "Use an `m` of at most ",
      fit$nkeep,
      ".",
      call. = FALSE
    )
  }
  ## Closed form, then guarded. m <= nkeep makes the step >= 1, so the
  ## rounded picks are strictly increasing inside [1, nkeep]; a wrong formula
  ## stops here instead of handing out one completion twice.
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

# The debug path: ONE completion off ONE chain of a live fit, on the full
# `obj`, for looking at. The draws are unpooled and unthinned, so these
# completions must not feed pooling. With a list in play `obj` is the kept
# columns: the listed ones are not drawn until mipca_finalize().
#' @export
inspect_completion <- function(x, chain = 1L, draw = 1L) {
  if (!inherits(x, "si_mi_pca")) {
    stop(
      "inspect_completion() takes an unfinished fit from mipca() with ",
      "`keep_state = TRUE`. For a finalized fit, use as.matrix().",
      call. = FALSE
    )
  }
  .life_check(x)
  chain <- .chk_count(chain, "chain", hi = length(x$draws))
  draw <- .chk_count(draw, "draw", hi = ncol(x$draws[[1L]]))
  out <- x$obj
  idx <- .map_lin(x$map, nrow(out))
  ## The same assertion as.matrix.si_mi() makes, for the same reason. Holes
  ## outside the narrowed subset have no stored draws and stay NA.
  if (!all(is.na(out[idx]))) {
    stop("internal error: a stored cell is not a hole.")
  }
  out[idx] <- .backward_at(x$transform, x$draws[[chain]][, draw], x$map)
  out
}

# A bare list holding a multi-GB matrix floods the terminal; these exist so
# it never does, and so the stage is visible without inspecting fields.
#' @export
print.si_mi_pca <- function(x, ...) {
  ## Safe to refuse: no printing route passes an argument, so this refuses
  ## only a caller's own extra one.
  .chk_no_dots(..., .what = "print()")
  ## The one method that does not call .life_check(): a broken object must
  ## not look healthy, so this reports the state instead of refusing.
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
  ## What the caller needs to choose the next call, and no more: the
  ## transform, the run settings and the tuner's record stay in the fields.
  cat(sprintf(
    "si_mi_pca fit, resumable: %d x %d, rank %d, %d chain%s x %d draws\n",
    nrow(x$obj),
    ncol(x$obj),
    x$prov$rank,
    length(x$draws),
    .s(length(x$draws)),
    ncol(x$draws[[1L]])
  ))
  ## Said only when true: the matrix above is the kept columns, and these
  ## are not in it until mipca_finalize() draws them.
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

# The plural of a count. isTRUE(n == 1L), not `==` and not identical():
# .prov_agree() returns NA_integer_ when two chunks disagree, where `if (NA)`
# would error, and a count may arrive as a double 1.
.s <- function(n) {
  if (isTRUE(n == 1L)) "" else "s"
}

# ONE value only when every chunk agrees, NA of the input's type otherwise:
# cbind() checks only rownames and nkeep. Every object-level read of a
# per-chunk field goes through here.
.agree <- function(vals) {
  if (length(vals) > 0L && isTRUE(all(vals == vals[1L]))) {
    vals[1L]
  } else {
    vals[NA_integer_]
  }
}

# How a disagreement prints. "mixed" rather than "NA", because NA in a header
# reads as "unknown" when the truth is "more than one answer".
.fmt_agree <- function(v) {
  if (is.na(v)) "mixed" else format(v)
}

# Per-chunk provenance read as an object-level field, through .agree().
.prov_agree <- function(chunks, nm) {
  .agree(vapply(chunks, function(ch) as.integer(ch$prov[[nm]]), 0L))
}

#' @export
print.si_mi <- function(x, ...) {
  .chk_no_dots(..., .what = "print()")
  h <- sum(vapply(x$chunks, function(ch) nrow(ch$map), 0L))
  ## What the caller needs to use the draws, and no more: chains, thinning,
  ## the ABB columns and the tuner's record stay in the chunks' fields.
  cat(sprintf(
    "si_mi fit: %d x %d, %d imputed cell%s, %d draws\n",
    length(x$rownames),
    length(x$colnames),
    h,
    .s(h),
    x$nkeep
  ))
  ## Said only when true: a column carried along and NOT imputed comes back
  ## from as.matrix() with its holes still NA. Counted off the chunks' own
  ## flag (.mipca_empty(leave_na = TRUE)), which only the flank driver sets;
  ## the flank plan's table says why.
  nleft <- sum(vapply(
    x$chunks,
    function(ch) if (isTRUE(ch$prov$left_na)) ncol(ch$obj) else 0L,
    0L
  ))
  if (nleft > 0L) {
    cat(sprintf(
      "  left NA: %d column%s not imputed. The flank plan's targets$status says why.\n",
      nleft,
      .s(nleft)
    ))
  }
  cat(
    "  next: mipca_complete() for completed matrices, summary() for Rhat and ESS\n"
  )
  invisible(x)
}

# ---- Layer 3b: the per-cell diagnostic frame ---------------------------------
# Stored holes sampled per column, deterministically, as a fraction of each
# column's holes, read on the latent scale by the C++ kernel
# (src/mipca_diag.cpp): classic split-Rhat with the Gelman fold, and one
# Geyer ESS.

# One owner of the default per-column fraction.
.prop_diag_default <- function() 0.5

# WHICH cells get diagnosed: ceiling(prop_diag * h) of each dirty column's h
# holes, evenly spread over its contiguous map run -- arithmetic on a
# tabulate(), guarded. A pure function of (map, prop_diag); not nested
# across prop_diag.
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
    ## ceiling(), so no column is ever DROPPED.
    k <- min(h, max(1L, as.integer(ceiling(prop_diag * h))))
    ## floor(), not round(): at spacing >= 1 floor stays strictly increasing,
    ## where round() can tie two neighbours.
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
  ## unname(): a ONE-row read comes back named ("col"), and identical()
  ## would then refuse any fit whose sample is a single cell.
  if (!identical(unname(map[rows, "col"]), dirty[owner])) {
    stop("internal error: a sampled cell is outside the column it came from.")
  }
  list(rows = rows)
}

# Per-cell Rhat / ESS off the kernel, for an ascending row set the caller
# has fixed. Shared by .diag_frame() and the tuner, which reads its diag
# cells straight off the stored map. Returns column 1 Rhat, column 2 ESS.
.cell_diag_read <- function(draws, srows) {
  res <- mipca_cell_diag(draws, srows, 256L)

  ## Masked here, not in the kernel, which follows posterior: at 11 or fewer
  ## draws per chain its ESS is chains * (draws %/% 2) whatever the data.
  if (ncol(draws[[1L]]) %/% 2L <= 5L) {
    res[, 2L] <- NA_real_
  }
  res
}

# THE FRAME, one row per cell, no reduction: per-cell measurements, not
# verdicts. The one builder for summary() and the freeze; `draws` must still
# carry the chains axis.
.diag_frame <- function(map, draws, p, nms, prop_diag) {
  prop_diag <- .chk_frac(prop_diag, "prop_diag")
  ## Asserted rather than branched on: mipca() requires unique colnames, so
  ## an internal error beats a frame of NA names.
  if (is.null(nms) || length(nms) != p) {
    stop("internal error: the diagnostic frame needs one name per column.")
  }
  ## Ascending by construction, which the kernel's tiled gather wants.
  srows <- .diag_sample(map, p, prop_diag)$rows

  S <- ncol(draws[[1L]])
  res <- .cell_diag_read(draws, srows)

  ## row.names = NULL: with ONE sampled cell, data.frame() would take a name
  ## ("rhat") as the row name.
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
    ## `ndraws`: the draws actually here after thinning/extension.
    ndraws = S,
    class = "si_mi_summary"
  )
}

# ON DEMAND, and nothing is stored on the resumable class: the frame feeds
# the extend-or-finish decision, so it describes the draws as they are NOW.
# Consumes no RNG.
#' @export
summary.si_mi_pca <- function(object, prop_diag = .prop_diag_default(), ...) {
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

# FROZEN: Rhat is not recomputable after pooling. Chunks keep their own
# frames (they may differ in prop_diag); `col` is the OUTPUT column, the one
# as.matrix() gives the cell, read off the same table (.chunk_to()).
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
  ## One value only when every chunk agrees. unlist() rather than a typed
  ## vapply so one helper serves prop_diag (double) and the two integers.
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

# The tables run to tens of thousands of rows: quantiles across cells for
# shape, then the worst cells by Rhat. Display only; `$cells` holds the
# measurements.
#' @export
print.si_mi_summary <- function(x, n = 6L, ...) {
  .chk_no_dots(..., .what = "print()")
  ## 0 is legal: the header and the shape table without the worst cells.
  n <- .chk_count(n, "n", lo = 0L)
  cells <- x$cells
  ## NA means the chunks disagree, printed as "mixed".
  num <- function(nm) .fmt_agree(attr(x, nm))
  ## `name`, not `col`: unique by .chk_matrix() and cbind(), so it is the
  ## identity that survives a join.
  ncol_aff <- length(unique(cells$name))
  ## prop_diag below 1 tells the reader the worst cells are the worst
  ## SEEN, not the true worst.
  cat(sprintf(
    "si_mi_summary: Rhat and ESS for %d cell%s in %d column%s\n",
    nrow(cells),
    .s(nrow(cells)),
    ncol_aff,
    .s(ncol_aff)
  ))
  cat(sprintf(
    "  prop_diag %s, %s chains x %s draws\n",
    num("prop_diag"),
    num("nchains"),
    num("ndraws")
  ))
  ## A cell with too few draws for ESS is NA and left out of its row.
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

# Every stored draw PAIRED with the transform that maps it back, since the
# two stages keep the transform in different places.
.draw_parts <- function(fit) {
  ## `col`/`row` ride with each part: a cell-aware backward (tf_biwhiten)
  ## needs to know which column and sample each draw is.
  if (inherits(fit, "si_mi_pca")) {
    ## The lifecycle gate, here for every reader: a passed fit's draws are
    ## NULL (clamp_frac() would read 0 / 0 = NaN) and a failed one's are torn.
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
  stop("`fit` must be a fit from mipca() or mipca_finalize().", call. = FALSE)
}

# Fraction of draws sitting exactly on the [0, 1] boundary after the backward
# map: draws a bounded map clamped, and genuine 0/1 reproduction. An exit:
# draws are stored latent, so the backward map is applied here, at the door.
#' @export
clamp_frac <- function(fit) {
  parts <- .draw_parts(fit)
  hit <- vapply(
    parts,
    function(p) {
      b <- p$tf$backward(p$d, p$col, p$row)
      ## Two sums, not `b == 0 | b == 1`: the sets are disjoint, so the
      ## count is the same integer, and this drops one draws-scale logical.
      sum(b == 0) + sum(b == 1)
    },
    0
  )
  sum(hit) / sum(vapply(parts, function(p) length(p$d), 0))
}
