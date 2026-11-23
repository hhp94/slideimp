# CRAN tier: the smokes, and nothing else. This file has no skip, so it is
# the whole of a CRAN run. A smoke calls an exported function on its main
# path, on a small input, on at most 2 cores and with no mirai daemons, and
# asserts the shape, the attributes, that no NA is left where none should be,
# and one value with a tolerance against a value computed without the code
# under test. Every export, and every compiled path, is named by some smoke.
# The error surface is in the local tier.
#
# CRAN runs these on BLAS, LAPACK, compilers and long-double settings that
# cannot be reproduced here, so a tolerance sits far above rounding: 1e-10
# where compiled code is checked against an R reference, 1e-8 where the value
# goes through BLAS or LAPACK. A real fault, a wrong neighbor or a wrong
# component, moves a value by orders of magnitude more. Counts and indices
# are exact, at 0.

# Fixtures -------------------------------------------------------------------

# K-NN by brute force in base R, the rule knn_imp() documents: the distance
# from column j to column c is the mean squared difference over the rows both
# observe, the k nearest are the neighbors, and a missing cell is the plain
# mean of the neighbors observed at its row. dist_pow = 0 only.
smoke_knn_ref <- function(x, k) {
  out <- x
  for (j in which(colSums(is.na(x)) > 0L)) {
    d <- vapply(
      seq_len(ncol(x)),
      function(cc) {
        both <- !is.na(x[, j]) & !is.na(x[, cc])
        if (cc == j || !any(both)) Inf else mean((x[both, j] - x[both, cc])^2)
      },
      numeric(1)
    )
    nb <- order(d)[seq_len(k)]
    for (r in which(is.na(x[, j]))) {
      out[r, j] <- mean(x[r, nb], na.rm = TRUE)
    }
  }
  out
}

# A 10 x 8 matrix with one hole, at [1, c1], and a trap for any method that
# looks past columns 1-4: c5 is an exact copy of c1, observed at row 1, so it
# is c1's nearest neighbor over the whole matrix, while c2, c1 plus a little
# noise, is its nearest inside columns 1-4. Imputed within c1-c4, the hole
# takes c2's value; imputed over all eight columns, it takes c5's.
smoke_trap_input <- function() {
  set.seed(42)
  x <- matrix(
    stats::rnorm(80),
    10,
    8,
    dimnames = list(paste0("r", 1:10), paste0("c", 1:8))
  )
  x[, 2] <- x[, 1] + 0.05 * stats::rnorm(10)
  x[, 5] <- x[, 1]
  x[1, 1] <- NA
  x
}

# An exactly rank-2 12 x 8 matrix and the same with 4 holes. Centered, it is
# still rank 2, so EM PCA with ncp = 2 converges to the true values.
smoke_rank2 <- function() {
  set.seed(1)
  truth <- matrix(stats::rnorm(24), 12, 2) %*%
    matrix(stats::rnorm(16), 2, 8) +
    5
  x <- truth
  x[cbind(c(1, 4, 7, 10), c(2, 3, 5, 8))] <- NA
  list(x = x, truth = truth)
}

# A matrix of proportions, rank 2 on the logit scale plus noise. By default
# the noise level differs by column, from 0.1 to 0.8, so that tf_biwhiten()'s
# prior on the column scales is finite; `noise = c(0.3, 0.3)` gives one level,
# whose held-out values look Gaussian to mipca_tune(). Up to `n_na` holes at
# random cells, and two more in c04.
smoke_beta <- function(n = 20L, p = 15L, n_na = 12L, noise = c(0.1, 0.8)) {
  set.seed(7)
  noise <- rep(seq(noise[1], noise[2], length.out = p), each = n)
  y <- stats::plogis(
    matrix(stats::rnorm(n * 2), n) %*%
      matrix(stats::rnorm(2 * p), 2) +
      noise * matrix(stats::rnorm(n * p), n)
  )
  dimnames(y) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  y[cbind(sample.int(n, n_na, TRUE), sample.int(p, n_na, TRUE))] <- NA
  y[c(3, 9), 4] <- NA
  y
}

# tf_biwhiten(tf_sv(n), rows = TRUE)$fit(y, rank), rebuilt in base R from its
# help page by another route: the sv-logit squeeze by hand, the top-rank left
# singular vectors by svd(), each column's residual and leverage by lm.fit(),
# and the Smyth prior's d0 by uniroot() on trigamma(). Every column must
# have at least rank + 3 observed values.
smoke_biwhiten_ref <- function(y, rank) {
  n <- nrow(y)
  z <- stats::qlogis((y * (n - 1) + 0.5) / n)
  obs <- !is.na(z)
  zf <- sweep(z, 2L, colMeans(z, na.rm = TRUE))
  zf[!obs] <- 0
  w <- svd(zf, nu = rank, nv = 0L)$u
  res <- matrix(0, n, ncol(z))
  lev <- numeric(n)
  for (j in seq_len(ncol(z))) {
    o <- obs[, j]
    f <- stats::lm.fit(w[o, , drop = FALSE], zf[o, j])
    res[o, j] <- f$residuals
    lev[o] <- lev[o] + 1 / sum(o) + rowSums(qr.Q(f$qr)^2)
  }
  df <- colSums(obs) - 1 - rank
  s2 <- colSums(res^2) / df
  e <- log(s2) - digamma(df / 2) + log(df / 2)
  evar <- stats::var(e) - mean(trigamma(df / 2))
  if (evar > 0) {
    d0 <- 2 *
      stats::uniroot(
        function(h) trigamma(h) - evar,
        c(1e-8, 1e8),
        tol = 1e-14
      )$root
    s20 <- exp(mean(e) + digamma(d0 / 2) - log(d0 / 2))
    s <- sqrt((d0 * s20 + df * s2) / (d0 + df))
  } else {
    s <- rep(sqrt(exp(mean(e))), ncol(z))
  }
  q <- colSums(t(res^2) / s^2)
  rs <- sqrt(q / (rowSums(obs) - lev))
  list(scale = s, rowscale = rs / exp(mean(log(rs))))
}

# Each completed matrix of a mipca() fit: the dims and dimnames of `y`, the
# observed cells of `y` as they are, and finite values in [0, 1] in the holes.
# Only these, of the draws: seeded draws move with the BLAS.
smoke_expect_draws <- function(fit, y, ndraws, lo = 0, hi = 1) {
  expect_identical(fit$nkeep, as.integer(ndraws))
  ms <- lapply(seq_len(ndraws), function(b) as.matrix(fit, draw = b))
  same_names <- vapply(ms, function(m) identical(dimnames(m), dimnames(y)), NA)
  expect_true(all(same_names))
  obs <- vapply(ms, function(m) m[!is.na(y)], numeric(sum(!is.na(y))))
  expect_identical(obs, matrix(y[!is.na(y)], nrow(obs), ndraws))
  holes <- vapply(ms, function(m) m[is.na(y)], numeric(sum(is.na(y))))
  expect_gt(length(holes), 0L)
  expect_true(all(is.finite(holes) & holes >= lo & holes <= hi))
  invisible(ms)
}

# Matrix statistics --------------------------------------------------------

test_that("smoke: col_vars() matches var() column by column", {
  x <- matrix(c(1, 4, 2, 8, NA, 3, NA, NA, 5, 0, 1, 1), 4, 3)
  colnames(x) <- c("a", "b", "c")
  got <- col_vars(x)
  expect_identical(names(got), colnames(x))
  # c has one observed value, so its variance is NA in both
  expect_worst_err(
    got,
    apply(x, 2, stats::var, na.rm = TRUE),
    abs_tol = 1e-10,
    rel_tol = 1e-10
  )
})

test_that("smoke: mat_miss() counts NA and NaN by column and by row", {
  x <- matrix(c(1, NA, 3, NaN, 5, 6, NA, NA, 9, Inf, 11, 12), 4, 3)
  dimnames(x) <- list(paste0("r", 1:4), paste0("c", 1:3))
  by_col <- mat_miss(x)
  by_row <- mat_miss(x, col = FALSE)
  expect_type(by_col, "integer")
  expect_identical(names(by_col), colnames(x))
  expect_identical(names(by_row), rownames(x))
  expect_worst_err(by_col, colSums(is.na(x)), abs_tol = 0, rel_tol = 0)
  expect_worst_err(by_row, rowSums(is.na(x)), abs_tol = 0, rel_tol = 0)
  expect_worst_err(
    mat_miss(x, prop = TRUE),
    colMeans(is.na(x)),
    abs_tol = 1e-10,
    rel_tol = 1e-10
  )
})

test_that("smoke: sim_mat() returns the matrix, its metadata and its NA", {
  set.seed(3)
  sim <- sim_mat(n = 10, p = 6, perc_total_na = 0.1, perc_col_na = 0.5)
  expect_s3_class(sim, "slideimp_sim")
  x <- sim$input
  expect_identical(dim(x), c(10L, 6L))
  expect_identical(colnames(x), sim$col_group$feature)
  expect_identical(rownames(x), sim$row_group$sample)
  expect_true(all(x >= 0 & x <= 1, na.rm = TRUE))
  # ceiling(0.1 * 10 * 6) holes, spread over floor(0.5 * 6) columns
  expect_worst_err(sum(is.na(x)), 6, abs_tol = 0, rel_tol = 0)
  expect_identical(sum(colSums(is.na(x)) > 0), 3L)
})

test_that("smoke: sample_na_loc() draws n_rows observed cells in n_cols columns", {
  set.seed(5)
  x <- matrix(stats::rnorm(60), 10, 6)
  x[c(1, 2), 1] <- NA
  x[5, 4] <- NA
  locs <- sample_na_loc(x, n_cols = 3L, n_rows = 2L, n_reps = 2L)
  expect_length(locs, 2L)
  expect_true(all(vapply(locs, is.integer, NA)))
  expect_identical(
    lapply(locs, dimnames),
    rep(list(list(NULL, c("row", "col"))), 2)
  )
  expect_identical(vapply(locs, anyDuplicated, 1L), c(0L, 0L))
  expect_false(anyNA(x[do.call(rbind, locs)]))
  # two cells in each of three distinct columns, in each repetition
  per_col <- unlist(lapply(locs, function(loc) as.vector(table(loc[, "col"]))))
  expect_worst_err(per_col, rep(2, 6), abs_tol = 0, rel_tol = 0)
})

# Single imputation ----------------------------------------------------------

test_that("smoke: knn_imp() on 2 cores matches K-NN by brute force", {
  set.seed(11)
  x <- matrix(
    stats::rnorm(120),
    12,
    10,
    dimnames = list(paste0("r", 1:12), paste0("c", 1:10))
  )
  x[cbind(c(1, 5, 2, 9, 3, 12), c(1, 1, 4, 4, 7, 10))] <- NA
  got <- knn_imp(x, k = 3, cores = 2)
  expect_s3_class(got, "slideimp_results")
  expect_identical(dimnames(got), dimnames(x))
  expect_identical(attr(got, "imp_method"), "knn")
  expect_false(anyNA(got))
  expect_worst_err(
    unclass(got)[is.na(x)],
    smoke_knn_ref(x, 3)[is.na(x)],
    abs_tol = 1e-10,
    rel_tol = 1e-10
  )
})

test_that("smoke: pca_imp() with the exact solver recovers a rank-2 matrix", {
  d <- smoke_rank2()
  got <- pca_imp(
    d$x,
    ncp = 2,
    scale = FALSE,
    method = "EM",
    solver = "exact",
    # miniter = maxiter fixes the iteration count, so the threshold decides
    # only `converged`. On an exactly rank-2 matrix the objective at the end
    # is rounding residue, and a threshold below that residue makes the flag,
    # and the non-convergence warning, turn on its sign, which differs by BLAS.
    threshold = 1e-8,
    maxiter = 300,
    miniter = 300
  )
  expect_s3_class(got, "slideimp_results")
  expect_identical(dim(got), dim(d$x))
  expect_identical(attr(got, "solver_chosen"), "exact")
  expect_gt(attr(got, "n_exact"), 0L)
  expect_identical(attr(got, "n_lobpcg_ok"), 0L)
  expect_false(anyNA(got))
  expect_identical(unclass(got)[!is.na(d$x)], d$x[!is.na(d$x)])
  expect_worst_err(
    unclass(got)[is.na(d$x)],
    d$truth[is.na(d$x)],
    abs_tol = 1e-8,
    rel_tol = 1e-8
  )
})

test_that("smoke: pca_imp() with LOBPCG and lobpcg_control() recovers a rank-2 matrix", {
  d <- smoke_rank2()
  ctrl <- lobpcg_control(warmup_iters = 2L, tol = 1e-12, maxiter = 50)
  expect_s3_class(ctrl, "slideimp_lobpcg_control")
  got <- pca_imp(
    d$x,
    ncp = 2,
    scale = FALSE,
    method = "EM",
    solver = "lobpcg",
    lobpcg_control = ctrl,
    # above the rounding residue, as in the exact-solver smoke
    threshold = 1e-8,
    maxiter = 300,
    miniter = 300
  )
  expect_s3_class(got, "slideimp_results")
  expect_identical(dim(got), dim(d$x))
  expect_identical(attr(got, "solver_chosen"), "lobpcg")
  # LOBPCG ran after the warmup, and its solves converged
  expect_gt(attr(got, "n_lobpcg_ok"), 0L)
  expect_false(anyNA(got))
  expect_worst_err(
    unclass(got)[is.na(d$x)],
    d$truth[is.na(d$x)],
    abs_tol = 1e-8,
    rel_tol = 1e-8
  )
})

# Windows and groups ---------------------------------------------------------

test_that("smoke: slide_imp() imputes each window on its own", {
  x <- smoke_trap_input()
  wins <- slide_imp(
    x,
    location = 1:8,
    window_size = 4,
    min_window_n = 2,
    dry_run = TRUE
  )
  expect_s3_class(wins, "slideimp_tbl")
  expect_identical(as.numeric(wins$start), c(1, 5))
  expect_identical(as.numeric(wins$end), c(4, 8))

  got <- slide_imp(
    x,
    location = 1:8,
    window_size = 4,
    min_window_n = 2,
    k = 1,
    .progress = FALSE
  )
  expect_s3_class(got, "slideimp_results")
  expect_identical(dimnames(got), dimnames(x))
  expect_identical(attr(got, "metacaller"), "slide_imp")
  expect_false(anyNA(got))
  want <- cbind(smoke_knn_ref(x[, 1:4], 1), smoke_knn_ref(x[, 5:8], 1))
  expect_worst_err(unclass(got)[1, 1], want[1, 1], 1e-10, 1e-10)
})

test_that("smoke: slideimp_resolve_group() and prep_groups() build the groups", {
  group <- data.frame(
    feature = c(paste0("c", 1:8), "c1"),
    group = c(rep("g1", 4), rep("g2", 4), "g1")
  )
  resolved <- slideimp_resolve_group(group)
  expect_s3_class(resolved, "data.frame")
  expect_identical(nrow(resolved), 8L)

  prepped <- prep_groups(paste0("c", 1:8), group)
  expect_s3_class(prepped, "slideimp_tbl")
  expect_identical(nrow(prepped), 2L)
  expect_true(all(c("feature", "aux", "parameters") %in% names(prepped)))
  expect_setequal(prepped$feature[[1]], paste0("c", 1:4))
  expect_setequal(prepped$feature[[2]], paste0("c", 5:8))
  expect_worst_err(lengths(prepped$feature), c(4, 4), abs_tol = 0, rel_tol = 0)
})

test_that("smoke: group_imp() imputes each group on its own", {
  x <- smoke_trap_input()
  group <- data.frame(
    feature = paste0("c", 1:8),
    group = rep(c("g1", "g2"), each = 4)
  )
  got <- suppressMessages(group_imp(x, group = group, k = 1, .progress = FALSE))
  expect_s3_class(got, "slideimp_results")
  expect_identical(dimnames(got), dimnames(x))
  expect_false(anyNA(got))
  want <- smoke_knn_ref(x[, 1:4], 1)
  expect_worst_err(unclass(got)[1, 1], want[1, 1], 1e-10, 1e-10)
})

# Tuning -------------------------------------------------------------------

test_that("smoke: tune_imp() scores knn_imp() on the cells it masks", {
  set.seed(13)
  x <- matrix(stats::rnorm(80), 10, 8)
  na_loc <- sample_na_loc(x, n_cols = 4L, n_rows = 2L, n_reps = 1L)
  loc <- na_loc[[1]]
  got <- suppressMessages(tune_imp(
    x,
    parameters = data.frame(k = c(1, 2)),
    .f = "knn_imp",
    na_loc = na_loc,
    .progress = FALSE
  ))
  expect_s3_class(got, "slideimp_tune")
  expect_identical(got$k, c(1, 2))
  expect_true(all(c("param_set", "rep_id", "result", "error") %in% names(got)))
  expect_true(all(is.na(got$error)))
  # the truth is x at the masked cells, and the estimate K-NN by brute force
  # on x with those cells masked
  expect_identical(lapply(got$result, `[[`, "truth"), list(x[loc], x[loc]))
  est <- vapply(got$result, `[[`, numeric(8), "estimate")
  expect_false(anyNA(est))
  masked <- x
  masked[loc] <- NA
  want <- cbind(smoke_knn_ref(masked, 1)[loc], smoke_knn_ref(masked, 2)[loc])
  expect_worst_err(est, want, abs_tol = 1e-10, rel_tol = 1e-10)
})

test_that("smoke: compute_metrics() scores a result by hand", {
  results <- data.frame(param_set = 1L)
  results$result <- list(data.frame(
    truth = c(1, 2, 3, 4),
    estimate = c(1.5, 2, NA, 3)
  ))
  got <- compute_metrics(results)
  expect_s3_class(got, "data.frame")
  expect_identical(got$.metric, c("mae", "rmse"))
  expect_identical(got$n, c(4L, 4L))
  expect_identical(got$n_miss, c(1L, 1L))
  expect_false(anyNA(got$.estimate))
  # over the three cells with an estimate: errors 0.5, 0 and 1
  expect_worst_err(
    got$.estimate,
    c(1.5 / 3, sqrt(1.25 / 3)),
    abs_tol = 1e-10,
    rel_tol = 1e-10
  )
})

# Multiple imputation --------------------------------------------------------

test_that("smoke: mipca() with tf_biwhiten(tf_sv()) and spectra_opts() fits its transform and draws", {
  y <- smoke_beta()
  tf <- tf_biwhiten(tf_sv(nrow(y)), rows = TRUE)
  opts <- spectra_opts(tol = 1e-10, maxitr = 1000)
  expect_s3_class(opts, "spectra_opts")
  fit <- mipca(
    y,
    rank = 2L,
    transform = tf,
    spectra = opts,
    chains = 2L,
    iter_warmup = 5L,
    iter_sampling = 5L,
    seed = 1
  )
  expect_s3_class(fit, "si_mi")
  smoke_expect_draws(fit, y, 10L)

  # the deterministic part: the column and row scales the fit carries
  fitted <- fit$chunks[[1]]$transform
  # the prior on the column scales is finite, so the scales are shrunk
  # towards it, not set equal
  expect_true(is.finite(fitted$d0))
  want <- smoke_biwhiten_ref(y, 2L)
  expect_worst_err(unname(fitted$scale), want$scale, 1e-10, 1e-10)
  expect_worst_err(fitted$rowscale, want$rowscale, 1e-10, 1e-10)
})

test_that("smoke: mipca_finalize() thins a kept fit with tf_identity()", {
  set.seed(17)
  y <- matrix(
    stats::runif(200, 0.2, 0.8),
    20,
    10,
    dimnames = list(sprintf("s%02d", 1:20), sprintf("c%02d", 1:10))
  )
  y[cbind(c(2, 7, 11, 15, 19), c(1, 3, 3, 6, 9))] <- NA
  kept <- mipca(
    y,
    rank = 2L,
    transform = tf_identity(),
    keep_state = TRUE,
    chains = 2L,
    iter_warmup = 4L,
    iter_sampling = 6L,
    seed = 1
  )
  expect_s3_class(kept, "si_mi_pca")
  fit <- mipca_finalize(kept, thin = 2L)
  expect_s3_class(fit, "si_mi")
  # 2 chains x 6 draws, every second kept
  expect_worst_err(fit$nkeep, 6, abs_tol = 0, rel_tol = 0)
  smoke_expect_draws(fit, y, 6L, lo = -Inf, hi = Inf)
})

test_that("smoke: abb_cols() lists the columns near 0 or 1, and mipca() draws them by ABB", {
  y <- smoke_beta()
  y[, 3] <- y[, 3] / 100
  y[, 5] <- 1 - (1 - y[, 5]) / 100
  y[c(1, 2), 3] <- NA
  y[4, 5] <- NA
  abb <- abb_cols(y)
  expect_worst_err(abb, c(3, 5), abs_tol = 0, rel_tol = 0)
  fit <- mipca(
    y,
    rank = 2L,
    transform = "logit_biwhiten",
    abb = abb,
    chains = 2L,
    iter_warmup = 5L,
    iter_sampling = 5L,
    seed = 1
  )
  expect_s3_class(fit, "si_mi")
  ms <- smoke_expect_draws(fit, y, 10L)
  # a column drawn by ABB takes only its own observed values
  from_own <- vapply(
    ms,
    function(m) {
      all(vapply(
        abb,
        function(j) all(m[is.na(y[, j]), j] %in% y[!is.na(y[, j]), j]),
        NA
      ))
    },
    NA
  )
  expect_true(all(from_own))
})

test_that("smoke: mipca_tune() scores each rank and mipca_fit() fits the best", {
  y <- smoke_beta(n = 20L, p = 20L, n_na = 16L, noise = c(0.3, 0.3))
  na_loc <- sample_na_loc(y, n_cols = 4L, n_rows = 2L, n_reps = 1L)
  tuned <- mipca_tune(
    y,
    rank = 1:2,
    na_loc = na_loc,
    transform = "logit_biwhiten",
    n_diag = 10L,
    chains = 2L,
    iter_warmup = 5L,
    iter_sampling = 5L,
    seed = 1
  )
  expect_s3_class(tuned, "si_tune")
  expect_identical(tuned$rank, 1:2)
  expect_true(all(is.na(tuned$error)))
  expect_true(all(is.finite(tuned$crps) & tuned$crps >= 0))
  # each rank is scored on the na_loc cells of y, as they are, and its crps
  # is the mean of their scores
  held <- na_loc[[1]][order(na_loc[[1]][, "col"], na_loc[[1]][, "row"]), ]
  cells <- lapply(tuned$result, function(res) {
    unname(cbind(res$row, res$col)[order(res$col, res$row), ])
  })
  expect_identical(cells, rep(list(unname(held)), 2))
  truth_ok <- vapply(
    tuned$result,
    function(res) identical(res$truth, y[cbind(res$row, res$col)]),
    NA
  )
  expect_true(all(truth_ok))
  expect_worst_err(
    tuned$crps,
    vapply(tuned$result, function(res) mean(res$crps), 1),
    1e-10,
    1e-10
  )

  fit <- mipca_fit(tuned, y, refresh = 0, seed = 1)
  expect_s3_class(fit, "si_mi")
  expect_identical(
    fit$chunks[[1]]$prov$rank,
    tuned$rank[which.min(tuned$crps)]
  )
  smoke_expect_draws(fit, y, fit$nkeep)
})

test_that("smoke: mipca_flank_plan(), mipca_flank_window() and mipca_flank() fit each target on its window", {
  y <- smoke_beta(n = 20L, p = 30L, n_na = 20L)
  y[c(2, 6), 8] <- NA
  y[c(5, 11), 20] <- NA
  location <- seq_len(ncol(y))
  plan <- mipca_flank_plan(
    y,
    location,
    window_size = 5,
    subset = c(8, 20),
    min_window_n = 5
  )
  expect_s3_class(plan, "si_flank_plan")
  tg <- plan$targets
  expect_identical(tg$name, c("c08", "c20"))
  expect_identical(tg$status, c("impute", "impute"))
  # every column within 5 of the target, on both sides
  want <- c(
    range(which(abs(location - 8) <= 5)),
    range(which(abs(location - 20) <= 5))
  )
  expect_worst_err(
    c(tg$start[1], tg$end[1], tg$start[2], tg$end[2]),
    want,
    0,
    0
  )

  win <- mipca_flank_window(plan, 20)
  expect_identical(dim(win), c(20L, 11L))
  expect_identical(colnames(win), sprintf("c%02d", 15:25))

  fit <- mipca_flank(
    plan,
    rank = 2,
    seed = 1,
    transform = "logit_biwhiten",
    chains = 2L,
    iter_warmup = 5L,
    iter_sampling = 5L,
    .progress = FALSE
  )
  expect_s3_class(fit, "si_mi")
  smoke_expect_draws(fit, y[, c(8, 20)], 10L)
})
