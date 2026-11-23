skip_on_cran()

# A small fit for the sampler's messages. One chain, so it runs under a
# timer build too.
msg_obj <- function() {
  set.seed(20261121L)
  n <- 20L
  p <- 30L
  X <- matrix(rnorm(n * 2L), n) %*%
    matrix(rnorm(2L * p), 2L) +
    0.3 * matrix(rnorm(n * p), n)
  dimnames(X) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  X[cbind(sample.int(n, 30L, TRUE), sample.int(p, 30L, TRUE))] <- NA_real_
  X
}

msg_fit <- function(X, ...) {
  mipca(
    X,
    rank = 2L,
    chains = 1L,
    iter_warmup = 4L,
    iter_sampling = 4L,
    seed = 1L,
    ...
  )
}

test_that("progress lines are messages, and suppressMessages() silences them", {
  X <- msg_obj()
  # Rcout once wrote these lines to stdout, where suppressMessages() could
  # not reach them
  msgs <- character(0)
  out <- capture.output(
    fit <- withCallingHandlers(
      msg_fit(X, transform = "identity", refresh = 1L),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
  )
  expect_identical(out, character(0))
  # one line per iteration (8) and the closing line
  expect_length(msgs, 9L)
  expect_match(msgs[1L], "Chain 1 Iteration:", fixed = TRUE)
  expect_match(msgs[9L], "Chain 1 finished in", fixed = TRUE)

  expect_silent(
    suppressMessages(msg_fit(X, transform = "identity", refresh = 1L))
  )
})

test_that("a 'spectra_opts' object edited by hand is refused in its own terms", {
  X <- msg_obj()
  o <- spectra_opts()
  o$tol <- -1
  expect_error(
    msg_fit(X, transform = "identity", spectra = o),
    "'tol' in spectra_opts() must be a positive number. It is -1",
    fixed = TRUE
  )
  o <- spectra_opts()
  o$maxitr <- 0L
  expect_error(
    msg_fit(X, transform = "identity", spectra = o),
    "'maxitr' in spectra_opts() must be at least 1. It is 0",
    fixed = TRUE
  )
})

test_that("a transform that breaks the hole pattern is named with its fix", {
  X <- msg_obj()
  tf <- function(forward) {
    list(name = "bad", forward = forward, backward = function(z, ...) z)
  }

  # an observed cell sent to NA: reported as NA or NaN, not NaN alone
  to_na <- tf(function(y, ...) {
    y[which(!is.na(y))[1L]] <- NA_real_
    y
  })
  expect_error(
    msg_fit(X, transform = to_na),
    "the transform maps observed values in column 1 to NA or NaN (1 of them)",
    fixed = TRUE
  )

  # every missing cell given a value
  fill <- tf(function(y, ...) {
    y[is.na(y)] <- 0
    y
  })
  expect_error(
    msg_fit(X, transform = fill),
    "A transform must map NA to NA",
    fixed = TRUE
  )
  expect_error(
    msg_fit(X, transform = fill),
    "the transform maps missing cells in column",
    fixed = TRUE
  )

  # an observed cell sent to Inf
  to_inf <- tf(function(y, ...) {
    y[which(!is.na(y))[1L]] <- Inf
    y
  })
  expect_error(
    msg_fit(X, transform = to_inf),
    "to a non-finite value. A transform must give a finite value",
    fixed = TRUE
  )
})
