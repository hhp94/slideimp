test_that("mipca_finalize() narrows the transform to match the re-based columns", {
  n <- 20L
  p <- 30L
  set.seed(9)
  Y <- matrix(stats::runif(n * p), n, p)
  Y <- sweep(Y, 2L, seq(0.01, 0.99, length.out = p), "*")
  Y[sample(seq_along(Y), 90L)] <- NA
  dimnames(Y) <- list(paste0("r", seq_len(n)), paste0("c", seq_len(p)))

  cnt <- colSums(is.na(Y))
  sel <- which(cnt > 0)
  sel <- sel[sel >= 20L][1:5]

  hole <- which(is.na(Y[, sel, drop = FALSE]))

  s <- tf_biwhiten(tf_identity())$fit(Y, 2L)$scale
  expect_gt(length(hole), 0L)
  expect_false(identical(as.integer(sel), seq_along(sel)))
  expect_false(isTRUE(all.equal(s[sel], s[seq_along(sel)])))

  run <- function(tf, cols) {
    set.seed(11)
    mipca(
      Y,
      rank = 2L,
      transform = tf,
      subset = cols,
      chains = 1L,
      iter_warmup = 5L,
      iter_sampling = 5L,
      seed = 11
    )
  }

  transforms <- list(
    identity = tf_identity(),
    `whiten(cols only)` = tf_biwhiten(tf_identity(), rows = FALSE),
    biwhiten = tf_biwhiten(tf_identity(), rows = TRUE)
  )

  for (nm in names(transforms)) {
    narrowed <- as.matrix(run(transforms[[nm]], sel), draw = 1L)[hole]
    full <- as.matrix(run(transforms[[nm]], NULL), draw = 1L)[, sel][hole]

    abs_err <- max(abs(narrowed - full))
    rel_err <- max(abs(narrowed - full) / abs(full))
    expect_equal(abs_err, 0, info = paste(nm, "worst absolute error"))
    expect_equal(rel_err, 0, info = paste(nm, "worst relative error"))
  }
})
