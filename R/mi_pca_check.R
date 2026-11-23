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
    mode = "double",
    row.names = "unique",
    col.names = "unique",
    .var.name = "obj"
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
