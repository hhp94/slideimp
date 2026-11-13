# Contract: mipca() repeats mipca_core()'s sampler defaults (it takes no
# `...`), and the two spellings must stay identical.

test_that("mipca()'s sampler defaults are mipca_core()'s", {
  nm <- c(
    "iter_sampling",
    "iter_warmup",
    "thin",
    "nu",
    "chains",
    "parallel_chains",
    "refresh",
    "spectra",
    "init"
  )
  # Non-vacuity: a typo'd name would compare NULL with NULL.
  expect_length(nm, 9L)
  expect_true(all(nm %in% names(formals(mipca))))
  expect_true(all(nm %in% names(formals(mipca_core))))

  # Unevaluated: `spectra_opts()` compares as written.
  expect_identical(formals(mipca)[nm], formals(mipca_core)[nm])
})

# Contract: `transform` has no default at a door that fits, so a beta matrix is
# never fitted on the identity scale by leaving the argument out.

test_that("transform has no default at mipca(), mipca_tune() and mipca_flank()", {
  doors <- list(
    mipca = mipca,
    mipca_tune = mipca_tune,
    mipca_flank = mipca_flank
  )
  for (nm in names(doors)) {
    f <- formals(doors[[nm]])
    # Non-vacuity: a door without the argument would deparse NULL, not "".
    expect_true("transform" %in% names(f), info = nm)
    expect_identical(deparse(f[["transform"]]), "", info = nm)
  }

  set.seed(1)
  X <- matrix(stats::rnorm(120), 12L, 10L)
  dimnames(X) <- list(sprintf("r%02d", 1:12), sprintf("c%02d", 1:10))
  X[2L, 3L] <- NA_real_
  # Non-vacuity: the same call fits once the transform is named.
  fit <- mipca(
    X,
    rank = 2L,
    transform = "identity",
    chains = 1L,
    iter_sampling = 4L,
    iter_warmup = 4L,
    seed = 1L
  )
  expect_s3_class(fit, "si_mi")
  expect_error(
    mipca(
      X,
      rank = 2L,
      chains = 1L,
      iter_sampling = 4L,
      iter_warmup = 4L,
      seed = 1L
    ),
    "\"transform\" is missing"
  )
})
