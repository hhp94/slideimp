skip_on_cran()

# A fit of `a` draws extended by update() to `a + b` against one run of
# `a + b` draws under the same seed. The first `a` draws come from the same
# run, so they are identical. A resume rebuilds `Xd` from the matrix and
# shifts it by each chain's column means, which moves the last bits of what
# follows, so the rest agree to rounding and are not asserted bitwise.
# Both layouts of the build are covered: split (n <= p with complete columns,
# their Gram fixed once and carried in the state) and tall (n > p).

resume_mat <- function(n, p, holes_in) {
  set.seed(20261122L)
  Z <- matrix(rnorm(n * 3L), n) %*%
    matrix(rnorm(3L * p), 3L) +
    0.4 * matrix(rnorm(n * p), n)
  dimnames(Z) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  Z[cbind(
    sample.int(n, 60L, TRUE),
    sample.int(holes_in, 60L, TRUE)
  )] <- NA_real_
  Z
}

resume_pair <- function(X, nu = Inf, seed = 1L) {
  args <- list(
    X,
    rank = 3L,
    transform = "identity",
    keep_state = TRUE,
    chains = 2L,
    iter_warmup = 10L,
    nu = nu
  )
  one <- do.call(mipca, c(args, list(iter_sampling = 12L, seed = seed)))
  two <- do.call(mipca, c(args, list(iter_sampling = 6L, seed = seed)))
  update(two, iter_sampling = 6L)
  list(one = one, two = two)
}

# The worst absolute and relative difference, over every stored cell and
# every draw after the first `a`, across chains; and whether those are
# bitwise identical, reported and not asserted.
resume_err <- function(one, two, a) {
  later <- seq.int(a + 1L, ncol(one$draws[[1L]]))
  x <- unlist(lapply(two$draws, function(m) m[, later]))
  y <- unlist(lapply(one$draws, function(m) m[, later]))
  list(
    n = length(y),
    abs = max(abs(x - y)),
    rel = max(abs(x - y) / abs(y)),
    bitwise = identical(x, y)
  )
}

for (case in list(
  list(name = "split", n = 30L, p = 60L, holes_in = 40L, nu = Inf),
  list(name = "tall", n = 60L, p = 20L, holes_in = 20L, nu = Inf),
  list(name = "split, Student-t", n = 30L, p = 60L, holes_in = 40L, nu = 7)
)) {
  test_that(
    paste0("update() extends a fit as one longer run would (", case$name, ")"),
    {
      X <- resume_mat(case$n, case$p, case$holes_in)
      # the layout this case claims
      cnt <- colSums(is.na(X))
      expect_identical(
        .split_cols(nrow(X), ncol(X), cnt)$split,
        case$name != "tall"
      )

      r <- resume_pair(X, nu = case$nu)
      one <- r$one
      two <- r$two
      expect_identical(two$map, one$map)
      expect_identical(lapply(two$draws, dim), lapply(one$draws, dim))
      expect_identical(two$prov$iter, one$prov$iter)

      # the first 6 draws of each chain: the same run, so bitwise
      for (ch in seq_along(one$draws)) {
        expect_identical(two$draws[[ch]][, 1:6], one$draws[[ch]][, 1:6])
      }

      e <- resume_err(one, two, a = 6L)
      expect_gt(e$n, 0L)
      message(sprintf(
        "resume (%s): %d values, max abs %.3g, max rel %.3g, bitwise %s",
        case$name,
        e$n,
        e$abs,
        e$rel,
        e$bitwise
      ))
      expect_lt(e$abs, 1e-12)
      expect_lt(e$rel, 1e-10)

      # the comparator can fail: another seed's extension is far outside it
      other <- resume_pair(X, nu = case$nu, seed = 2L)$two
      f <- resume_err(one, other, a = 6L)
      expect_gt(f$abs, 1e-3)
    }
  )
}
