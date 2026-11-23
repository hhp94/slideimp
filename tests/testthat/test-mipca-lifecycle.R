skip_on_cran()

# The fit both blocks use. `ballast`, when given, is held in this function's
# frame, the one mipca() is called from.
lifecycle_fit <- function(ballast = NULL) {
  force(ballast)
  set.seed(20260819L)
  n <- 30L
  p <- 60L
  Z <- matrix(rnorm(n * 3L), n) %*%
    matrix(rnorm(3L * p), 3L) +
    0.4 * matrix(rnorm(n * p), n)
  X <- plogis(Z / 2)
  dimnames(X) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  X[cbind(sample.int(n, 60L, TRUE), sample.int(40L, 60L, TRUE))] <- NA_real_
  mipca(
    X,
    rank = 3L,
    transform = "identity",
    keep_state = TRUE,
    chains = 2L,
    iter_sampling = 8L,
    iter_warmup = 10L,
    seed = 1L
  )
}

test_that("a fit torn by a failed update() refuses every method but print()", {
  fit <- lifecycle_fit()
  expect_identical(fit$lifecycle, "ok")

  fit$draws[[2L]] <- fit$draws[[2L]][-1L, , drop = FALSE]
  expect_error(update(fit, 4L))
  expect_identical(
    c(ncol(fit$draws[[1L]]), ncol(fit$draws[[2L]])),
    c(12L, 8L)
  )

  expect_identical(fit$lifecycle, "failed")
  expect_error(update(fit, 4L), "marked failed")
  expect_error(mipca_finalize(fit), "marked failed")
  expect_error(as.matrix(fit), "marked failed")
  expect_error(summary(fit), "marked failed")

  expect_output(print(fit), "FAILED")
})

test_that("a finalized fit and its summary print their counts and Rhat", {
  kept <- lifecycle_fit()
  holes <- colSums(is.na(kept$obj))
  fit <- mipca_finalize(kept, thin = 2L)

  out <- capture_output(ret <- expect_invisible(print(fit)))
  expect_identical(ret, fit)
  # 2 chains x 8 draws, every second kept
  expect_identical(
    strsplit(out, "\n")[[1L]],
    c(
      sprintf(
        "si_mi fit: %d x %d, %d imputed cells, 8 draws",
        nrow(kept$obj),
        ncol(kept$obj),
        sum(holes)
      ),
      "  next: as.matrix() with 'draw' for a completed matrix, summary() for Rhat and ESS"
    )
  )

  s <- summary(fit)
  rhat <- s$cells$rhat
  out <- capture_output(ret <- expect_invisible(print(s, n = 1L)))
  expect_identical(ret, s)
  lines <- strsplit(out, "\n")[[1L]]
  # prop_diag = 0.5 samples ceiling(h / 2) cells of a column with h holes
  expect_identical(
    lines[c(1L, 2L, 6L)],
    c(
      sprintf(
        "si_mi_summary: Rhat and ESS for %d cells in %d columns",
        sum(ceiling(holes / 2)),
        sum(holes > 0)
      ),
      "  prop_diag 0.5, 2 chains x 4 draws",
      "worst 1 cell by rhat:"
    )
  )
  q <- strsplit(trimws(lines[4L]), " +")[[1L]]
  expect_identical(q[1L], "rhat")
  expect_worst_err(
    as.numeric(q[-1L]),
    unname(signif(quantile(rhat, c(0, 0.25, 0.5, 0.75, 1)), 3L)),
    1e-12,
    1e-12
  )
  worst <- s$cells[which.max(rhat), ]
  expect_identical(
    strsplit(trimws(lines[8L]), " +")[[1L]][1:4],
    as.character(c(worst$chunk, worst$col, worst$name, worst$row))
  )
})

test_that("a finalized fit refuses [ and update(), saying what to do instead", {
  msg <- function(expr) tryCatch(expr, error = conditionMessage)
  kept <- lifecycle_fit()
  fit1 <- mipca_finalize(lifecycle_fit())
  fit2 <- mipca_finalize(kept, thin = 2L)
  nostate <- mipca(
    kept$obj,
    rank = 3L,
    transform = "identity",
    chains = 1L,
    iter_sampling = 1L,
    iter_warmup = 2L,
    seed = 1L
  )
  expect_s3_class(nostate, "si_mi")

  expect_identical(
    msg(fit1[, 1L]),
    paste0(
      "[ does not work on an \"si_mi\" fit. Choose columns at fit time with ",
      "'subset', join fits with cbind(), or subset a completed matrix: ",
      "as.matrix(fit, draw = b)[, j]"
    )
  )
  finalized <- paste0(
    "update() cannot extend a finalized fit. Fit again with mipca() and ",
    "'keep_state = TRUE'"
  )
  expect_identical(msg(update(fit1)), finalized)
  expect_identical(
    msg(update(fit2)),
    paste0(finalized, ", with 'thin = 2' to keep this fit's thinning")
  )
  expect_identical(
    msg(update(nostate)),
    paste0(
      "update() cannot extend this fit: it holds no sampler state to resume ",
      "from. Fit again with mipca() and 'keep_state = TRUE'"
    )
  )
})

test_that("a fit does not serialize the frame it was built in", {
  ballast_mb <- 8
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)

  # the ballast sits in the frame mipca() is called from, so a fit that kept
  # that frame would carry it into the file
  ballast <- matrix(runif(1e6L), 1000L)
  expect_gt(as.numeric(object.size(ballast)) / 1e6, ballast_mb - 0.1)
  saveRDS(lifecycle_fit(ballast), path)

  expect_lt(file.size(path) / 1e6, 2)
})
