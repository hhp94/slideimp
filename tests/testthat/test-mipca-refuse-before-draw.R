skip_on_cran()

test_that("bad settings are refused before R's stream is drawn", {
  set.seed(1)
  Y <- matrix(stats::rnorm(320), 16, 20)
  Y[sample(seq_along(Y), 30L)] <- NA
  dimnames(Y) <- list(paste0("r", 1:16), paste0("c", 1:20))
  L <- colnames(Y)[colSums(is.na(Y)) > 0L]

  seed_moves <- function(expr) {
    set.seed(7)
    s0 <- get(".Random.seed", envir = globalenv())
    force(expr)
    !identical(s0, get(".Random.seed", envir = globalenv()))
  }

  expect_true(seed_moves(mipca(
    Y,
    2L,
    "identity",
    chains = 1L,
    iter_warmup = 4L,
    iter_sampling = 4L
  )))
  expect_true(seed_moves(mipca(
    Y,
    2L,
    "identity",
    subset = L,
    abb = L,
    chains = 1L,
    iter_sampling = 4L
  )))

  expect_false(seed_moves(expect_error(
    mipca(Y, 2L, "identity", prop_diag = 2),
    "'prop_diag' must be a single fraction in (0, 1",
    fixed = TRUE
  )))
  expect_false(seed_moves(expect_error(
    mipca(Y, 2L, "identity", keep_state = NA),
    "'keep_state' must be TRUE or",
    fixed = TRUE
  )))
  expect_false(seed_moves(expect_error(
    mipca(Y, 2L, "identity", subset = L, abb = L, chains = 0L),
    "'chains' must be a single integer >= 1",
    fixed = TRUE
  )))
})
