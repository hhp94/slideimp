# mi_pca_abb.R --------------------------------------------------------------
# Regioned imputation: the helper that proposes which columns leave the
# sampler (abb_cols()), the per-column approximate Bayesian bootstrap that
# draws their holes (.abb_col()), the finalized piece that carries those
# draws and joins like any other (.abb_piece(), .abb_join()), and what each
# door needs to take the list: the plain door's split (.abb_split()), the
# tuner's holdout re-base (.abb_tune_rebase()), and the complement both read
# (.abb_kept()), which the flank plan reads too. The caller lists the
# columns and ABB draws them; no door decides the list.
# ------------------------------------------------------------------------------

# WHICH columns to list, read off the raw matrix alone. The quantile rule in
# its counting form (type 1, no interpolation): the low tail marks a column
# when at least ceiling(q * n_obs) of its observed cells are below `bound`,
# the high tail when that many are above 1 - `bound`. Both are counts, so the
# rule is the same on any monotone scale and needs no sort. Returns ascending
# column indices in `obj`'s own numbering.
# For the docs: `bound` is where the model stops being trusted, not a
# detection limit; a wider bound takes the unmethylated peak, not a fringe.
#' @export
abb_cols <- function(obj, bound = 0.02, q = 0.95) {
  .chk_matrix(obj)
  ## Under one half, the two bands cannot meet.
  if (
    !is.numeric(bound) ||
      length(bound) != 1L ||
      is.na(bound) ||
      bound <= 0 ||
      bound >= 0.5
  ) {
    stop("`bound` must be a single number in (0, 0.5).", call. = FALSE)
  }
  ## Above one half, or one column could clear the count in both tails.
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q <= 0.5 || q > 1) {
    stop("`q` must be a single number in (0.5, 1].", call. = FALSE)
  }
  n_obs <- nrow(obj) - .mat_miss_int(obj)
  ## Caller input, checked once. A matrix off the [0, 1] scale would have
  ## every negative cell counted as "in the low band". An all-NA matrix has
  ## no range and marks nothing below; the census says so without a mask.
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

  ## Each comparison allocates one logical matrix the size of `obj`.
  lo <- colSums(obj < bound, na.rm = TRUE)
  hi <- colSums(obj > 1 - bound, na.rm = TRUE)
  ## The 1e-9 keeps q * n_obs from tipping over an integer on fp noise. The
  ## floor is what q > 1/2 means in counts -- more than half -- and holds
  ## the nudge to it: for a q within 1e-9 / n_obs of one half the nudged
  ## count would be exactly half, which both tails can clear at once.
  need <- pmax(ceiling(q * n_obs - 1e-9), n_obs %/% 2L + 1L)
  ## A column with no observed cell clears "0 of 0" in both tails and has no
  ## donor to draw from: never marked.
  has <- n_obs >= 1L
  m_lo <- has & lo >= need
  m_hi <- has & hi >= need
  ## The arithmetic the two range checks above buy.
  if (any(m_lo & m_hi)) {
    stop("internal error: a column is marked in both tails.")
  }
  ## One-time DISCOVERY of the list, not a position lookup.
  unname(which(m_lo | m_hi))
}

# ONE listed column's holes, drawn `m` times by the approximate Bayesian
# bootstrap (Rubin and Schenker 1986). Per kept draw: resample the donors
# (the observed cells) with replacement to the donor count, then draw each
# hole from that pool with replacement. The two stages are the method: holes
# of one draw share a pool, which is what carries the uncertainty about the
# column's distribution into the draws.
#
# `x` is the raw column; a hole is is.na() (NaN included), the predicate
# .mat_miss_int() counts by. Returns `row`, the holes in ascending order (the
# map's order within a column), and `draws`, holes x m, on the DATA scale:
# every value is one of the column's donors, so an exact zero stays one.
#
# Consumes dqrng's current generator and neither seeds nor restores it: the
# caller seeds one generator per column and puts the state back. Base R's
# stream is not touched.
.abb_col <- function(x, m) {
  m <- .chk_count(m, "m")
  ## One-time DISCOVERY of the holes. unname(): a column cut from a matrix
  ## carries its row names, and a named `row` fails the map's identical()
  ## guards.
  row <- unname(which(is.na(x)))
  h <- length(row)
  nd <- length(x) - h
  if (h < 1L) {
    stop("internal error: ABB was asked for a column with no hole.")
  }
  if (nd < 1L) {
    stop(
      "A column listed in `abb` has no observed cell to draw its missing ",
      "cells from.",
      call. = FALSE
    )
  }
  donors <- x[-row]
  if (!all(is.finite(donors))) {
    stop(
      "A column listed in `abb` has a non-finite observed cell. Its missing ",
      "cells are drawn from its observed cells, so those must be finite.",
      call. = FALSE
    )
  }
  ## Both index vectors are integer; past this, nd * m would be NA.
  if (as.double(max(nd, h)) * m > .Machine$integer.max) {
    stop("internal error: the ABB index for one column overflows an integer.")
  }
  ## Draw b's pool is pool[(b - 1) * nd + 1 .. b * nd], positions in `donors`.
  pool <- dqrng::dqsample.int(nd, nd * m, replace = TRUE)
  ## Each hole's slot within its own draw's pool, then that slot's position
  ## in `pool`: arithmetic on the block offset, guarded.
  pick <- dqrng::dqsample.int(nd, h * m, replace = TRUE)
  at <- pick + rep((seq_len(m) - 1L) * nd, each = h)
  if (anyNA(at) || min(at) < 1L || max(at) > nd * m) {
    stop("internal error: an ABB pick falls outside its draw's pool.")
  }
  list(row = row, draws = matrix(donors[pool[at]], h, m))
}

# The list's complement, one owner for the plain door and the tuner: which
# of `p` columns are listed (a mask) and which are kept (ascending indices).
# `abb` is resolved.
.abb_kept <- function(p, abb) {
  if (length(abb) == p) {
    stop(
      "`abb` lists every column, so no column is left to fit. Leave at least ",
      "one column out of `abb`.",
      call. = FALSE
    )
  }
  is_abb <- logical(p)
  is_abb[abb] <- TRUE
  ## One-time DISCOVERY of the complement.
  list(is_abb = is_abb, kept = which(!is_abb))
}

# THE RE-BASE at the tuner's door. The tuner runs on the kept columns alone
# -- its holdout lies in them, so rank and nu are tuned for the sampler's
# side -- and a supplied `na_loc` arrives in the caller's numbering. `lin`
# is that design as .resolve_holdout() returns it against the caller's
# matrix `obj`: ascending linear indices, every shape refusal already made,
# by the one function that owns them. Returned in the kept matrix's linear
# indexing. A design that touches a listed column is refused: a listed
# column is not fitted, so it has no score to give.
.abb_tune_rebase <- function(lin, obj, kept) {
  n <- nrow(obj)
  p <- ncol(obj)
  col <- .lin_col(lin, n)
  ## Scatter table, caller's column -> kept position, 0 = listed.
  inv <- integer(p)
  inv[kept] <- seq_along(kept)
  kc <- inv[col]
  hit <- kc == 0L
  if (any(hit)) {
    ## Caller input, once: unique() names the columns, it finds no slot.
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
      "]. A column listed in `abb` is ",
      "not fitted, so it cannot be scored. Place the holdout in columns ",
      "that are not listed.",
      call. = FALSE
    )
  }
  ## Same row, `col - kc` columns to the left. `kept` ascends, so the order
  ## of the cells is the order they came in.
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

# THE SPLIT at the plain door: which columns the sampler keeps, where the
# requested ones fall in the kept numbering, and the record of the requested
# listed columns that rides on the fit until they are drawn. `abb` and
# `subset` are resolved (ascending indices of `obj`; `subset` NULL for all).
#
# KEYS. A listed column draws off dqset.seed(c(seed, key)) with its key its
# index in `obj`, so its draws move neither with `subset` nor with the rest
# of the list, and it draws as it would as a flank target.
.abb_split <- function(obj, abb, subset) {
  p <- ncol(obj)
  cm <- .abb_kept(p, abb)
  is_abb <- cm$is_abb
  kept <- cm$kept
  req <- if (is.null(subset)) seq_len(p) else subset
  req_abb <- req[is_abb[req]]
  req_kept <- req[!is_abb[req]]
  ## Kept numbering by scatter table; a surviving 0 means a listed column
  ## was taken for a kept one.
  inv <- integer(p)
  inv[kept] <- seq_along(kept)
  sub_k <- inv[req_kept]
  if (length(sub_k) > 0L && min(sub_k) < 1L) {
    stop("internal error: a requested kept column maps outside the kept set.")
  }

  rec <- NULL
  if (length(req_abb) > 0L) {
    blk <- obj[, req_abb, drop = FALSE]
    ## Refused HERE, before a sampler runs, not at finalize after it: no
    ## donor, or a non-finite cell. One donor is enough. Read off the CELLS,
    ## the predicate .abb_col() draws by, and not off the variance: a
    ## one-cell column has none, and finite cells can overflow one. NaN is a
    ## hole (is.na()), so the only non-finite donor is an Inf.
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
        "`abb` lists columns that cannot be drawn. A listed column's missing ",
        "cells are drawn from its own observed cells, so it needs at least ",
        "one, and all of them finite.\n",
        paste(
          c(
            show(cen$empty, "no observed cells"),
            show(cen$nonfin, "a non-finite cell")
          ),
          collapse = "\n"
        ),
        "\nRemove these columns from `obj` before fitting.",
        call. = FALSE
      )
    }
    ## THE PLACEMENT (.si_mi_new()). The join holds the requested kept
    ## columns, then the requested listed ones, each ascending; `req`
    ## ascends too, so a column's output position is its rank in `req`. By
    ## scatter table: placement[k] is where joined column k goes. Guarded
    ## here by reading it back -- output column placement[k] is requested
    ## column req[placement[k]], which has to be joined column k itself --
    ## and by the constructor as a permutation. Two ascending runs alone
    ## would pass any interleaving of the two sides.
    pos <- integer(p)
    pos[req] <- seq_along(req)
    joined <- c(req_kept, req_abb)
    placement <- pos[joined]
    if (!identical(req[placement], joined)) {
      stop(
        "internal error: the placement does not return the requested columns to the input's order."
      )
    }
    rec <- list(obj = blk, key = req_abb, placement = placement, seed = NULL)
  }
  list(
    kept = kept,
    ## NULL when every kept column is requested, as it is with no list.
    subset = if (length(sub_k) == length(kept)) NULL else sub_k,
    n_req_kept = length(sub_k),
    rec = rec
  )
}

# THE LISTED SIDE, drawn and joined: the one place a record becomes a piece.
# `rec` is .abb_split()'s record with its seed set (.resolve_seed(), by the
# caller, after its last refusal); `kept` is the kept side's finalized piece,
# or NULL when no kept column was requested. mipca_finalize() calls this
# once the kept-draw count is fixed; .mipca_run() calls it at fit time on
# the two paths where no sampler ran. The listed chunk is stored AFTER the
# kept one, and the record's placement vector (`rec$placement`,
# .abb_split()) puts every column back where the input had it, so the exits
# return the input's column order.
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
  ## cbind()'s join (its checks: row names, kept draws, a name joined once),
  ## under the door's call and the record's placement. `rec$placement` was
  ## derived for two pieces that each hold their columns ascending; the join
  ## refuses a piece that already moved one.
  .si_mi_join(list(kept, pc), call = call, placement = rec$placement)
}

# THE PIECE: a finalized one-chunk si_mi over the listed columns `cols` of
# `obj`, every hole drawn by .abb_col(). It joins a fitted piece through
# cbind.si_mi() and every exit reads it unchanged: the draws are stored on
# the data scale, so its backward map is the identity. One column is the
# flank driver's listed target; several are the plain door's listed side.
#
# Built here and not by mipca_finalize(): there are no chains to trim or
# pool and nothing to diagnose, since ABB draws are independent by
# construction. Its frame is the empty one .diag_frame() gives a map with no
# rows, so its attributes agree with a fitted partner's.
#
# `key` is one positive integer per column of `cols`, by every caller the
# column's index in the caller's matrix. A column draws off
# dqset.seed(c(seed, key)), so its draws depend on (seed, its key, its
# cells, the kept-draw count) and on no other column. dqrng packs the pair
# into one 64-bit seed, `seed` high and `key` low, so for seed >= 1 no key
# gives the scalar seed a sampler's chains start from, and the cost does not
# grow with the key as a stream number's does. The caller's dqrng state is
# put back; base R's stream is not touched. `chains` and `ndraws` (kept
# draws per chain) are the partner's: cbind() matches on their product, and
# summary() reads them.
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
  ## The doors resolve and sort `cols`; this guards their arithmetic.
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
  ## Whole numbers dqset.seed() takes as one 32-bit half of the seed.
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

  ## Row names dropped: they are hoisted onto si_mi, as in .rebase_subset().
  sub <- obj[, cols, drop = FALSE]
  dimnames(sub) <- list(NULL, colnames(sub))
  nms <- colnames(sub)

  ## The map is column-major, so column i's holes are rows off[i] + 1 ..
  ## off[i + 1] of it: arithmetic on the census, guarded below.
  cnt <- .mat_miss_int(sub)
  off <- cumsum(c(0L, cnt))
  nh <- off[k + 1L]
  map <- matrix(0L, nh, 2L, dimnames = list(NULL, c("row", "col")))
  draws <- .alloc_draws(
    matrix(NA_real_, nh, m),
    cells = nh,
    chains = chains,
    draws = ndraws,
    fix = "List fewer columns in `abb`, or ask for fewer `chains` or `iter_sampling`."
  )

  old <- dqrng::dqrng_get_state()
  on.exit(dqrng::dqrng_set_state(old), add = TRUE)
  for (i in seq_len(k)) {
    ## A listed column with no hole is carried and draws nothing.
    if (cnt[i] == 0L) {
      next
    }
    dqrng::dqset.seed(c(seed, key[i]))
    ## Re-thrown with the column's name; the doors refuse a donor-less
    ## column before any draw, so this is the piece's own last word.
    r <- tryCatch(
      .abb_col(sub[, i], m),
      error = function(e) {
        stop("Column ", nms[i], ": ", conditionMessage(e), call. = FALSE)
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
  ## Every stored cell is a hole, in which(is.na()) order, and every one of
  ## them was drawn. The index in double (`to = `): the integer form is NA
  ## where the listed block passes 2^31 cells.
  lin <- .map_lin(map, n, to = seq_len(k))
  if (
    anyNA(draws) || is.unsorted(lin, strictly = TRUE) || !all(is.na(sub[lin]))
  ) {
    stop("internal error: the ABB map is not the listed columns' holes.")
  }

  ## The empty frame, off the one builder: no ABB cell gets an Rhat or ESS.
  diag <- .diag_frame(
    map = map[integer(0), , drop = FALSE],
    draws = rep(list(matrix(NA_real_, 0L, ndraws)), chains),
    p = k,
    nms = nms,
    prop_diag = prop_diag
  )
  prov <- list(
    chains = chains,
    thin = thin,
    ndraws = ndraws,
    nmiss = nh,
    thin_final = thin_final,
    thin_effective = thin * thin_final,
    seam_state = seam_state,
    ## Marks an ABB piece; with `seed` and `key`, what would redraw it.
    abb = TRUE,
    seed = seed,
    key = key
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
