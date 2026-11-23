.STATE_VERSION <- 10L

.state_check_version <- function(.state) {
  if (
    !inherits(.state, "si_mi_pca_state") ||
      !identical(.state$state_version, .STATE_VERSION)
  ) {
    stop(
      "update() cannot extend this fit: its saved sampler state is not in ",
      "the format this version reads. Fit again with mipca() and ",
      "'keep_state = TRUE'",
      call. = FALSE
    )
  }
  invisible(.state)
}

.resume_refused <- function(...) {
  paste0(
    "update() cannot extend this fit: the matrix it carries no longer ",
    "matches the one its chains ran on. ",
    ...,
    " Fit again with mipca() ",
    "and 'keep_state = TRUE'"
  )
}

.start_resume <- function(.state, Xd, etd, dirty, miss_idx_l, n, Z) {
  list(
    xd_list = lapply(seq_along(.state$chain), function(ch) {
      st_ch <- .state$chain[[ch]]
      Xc <- Xd + rep(etd - st_ch$etd, each = n)
      Xc[miss_idx_l] <- st_ch$holes
      .state_fp_check(st_ch$fp, Xc, dirty, Z)
      Xc
    }),
    etd_list = lapply(.state$chain, `[[`, "etd"),
    xtild_list = lapply(.state$chain, `[[`, "xtild_miss"),
    eig_start_list = lapply(.state$chain, `[[`, "eig_start"),
    sigma_in = vapply(.state$chain, `[[`, 0, "sigma"),
    rng_state_in = vapply(.state$chain, `[[`, "", "rng"),
    disp = numeric(0)
  )
}

.state_check_args <- function(.state, chains, chains_given, forward, store) {
  .state_check_version(.state)
  if (chains_given && as.integer(chains) != .state$chains) {
    stop(
      "internal error: 'chains' disagrees with the state; a resume continues ",
      "every chain of the fit"
    )
  }
  if (!identical(forward, identity)) {
    stop(
      "internal error: 'forward' cannot be given on a resume: the state ",
      "carries the map the fit was built under, and the rebuild applies ",
      "that one"
    )
  }
  if (!is.function(.state$forward)) {
    stop("internal error: the state carries no transform map")
  }
  if (!is.null(store)) {
    stop(
      "internal error: 'store' cannot be given on a resume: the state ",
      "carries the store set the fit was built under, and an extension ",
      "keeps those same cells"
    )
  }
  .state$chains
}

.state_check_match <- function(.state, n, cnt, rank, nu, spectra, Z) {
  if (.state$n != n || !identical(.state$cnt, cnt)) {
    stop(
      .resume_refused(
        "Its row count or its count of missing cells per column has changed."
      ),
      call. = FALSE
    )
  }
  if (.state$rank != rank) {
    stop("internal error: the state was written at a different rank")
  }
  if (.state$nu != nu) {
    stop("internal error: the state was written with a different nu")
  }
  if (!identical(.state$spectra, spectra)) {
    stop(
      "internal error: the state was written with different spectra settings"
    )
  }
  .state_fp_clean_check(.state, Z, n, cnt)
  invisible(.state)
}

.state_fp_cols <- function(nd) {
  m <- min(nd, max(100L, as.integer(ceiling(0.05 * nd))))
  unique(as.integer(round(seq.int(1L, nd, length.out = m))))
}

.state_fp_take <- function(Xd, loc) {
  list(
    col = loc,
    sumabs = colSums(abs(Xd[, loc, drop = FALSE]))
  )
}

.STATE_FP_TOL <- 1e-8

.state_fp_check <- function(fp, Xd, dirty, Z) {
  cur <- .state_fp_take(Xd, fp$col)
  if (!all(is.finite(fp$sumabs)) || any(fp$sumabs <= 0)) {
    stop("internal error: the state's fingerprint is not positive and finite")
  }
  rel <- abs(cur$sumabs - fp$sumabs) / fp$sumabs
  if (max(rel) >= .STATE_FP_TOL) {
    worst <- which.max(rel)
    stop(
      .resume_refused(
        "Its values, or what its transform makes of them, have changed. ",
        "Column ",
        .col_label(Z, dirty[fp$col[worst]]),
        " differs most."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.state_fp_clean_take <- function(Z, clean) {
  sa <- colSums(abs(Z[, clean, drop = FALSE]))
  if (!all(is.finite(sa))) {
    stop("internal error: a clean column's fingerprint is not finite")
  }
  list(col = clean, sumabs = sa)
}

.state_fp_clean_check <- function(.state, Z, n, cnt) {
  split <- .split_cols(n, length(cnt), cnt)$split
  if (!split) {
    return(invisible(NULL))
  }
  fp <- .state$fp_clean
  if (
    !is.list(fp) ||
      !is.numeric(fp$sumabs) ||
      !is.integer(fp$col) ||
      length(fp$col) != length(fp$sumabs) ||
      length(fp$col) == 0L
  ) {
    stop(
      "internal error: the state carries no clean-column fingerprint but ",
      "was built from a split matrix"
    )
  }
  if (max(fp$col) > length(cnt) || any(cnt[fp$col] != 0L)) {
    stop(
      "internal error: the clean fingerprint names a column that is not ",
      "clean in this scan"
    )
  }
  cur <- colSums(abs(Z[, fp$col, drop = FALSE]))
  rel <- abs(cur - fp$sumabs) / pmax(fp$sumabs, .Machine$double.xmin)
  if (!all(is.finite(rel)) || max(rel) >= .STATE_FP_TOL) {
    worst <- which.max(rel)
    stop(
      .resume_refused(
        "Its values have changed. Column ",
        .col_label(Z, fp$col[worst]),
        " differs most."
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}
