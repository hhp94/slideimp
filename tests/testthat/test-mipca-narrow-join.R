.nj <- local({
  set.seed(20261001L)
  n <- 30L
  p <- 40L
  Xc <- matrix(stats::rnorm(n * p), n, p)
  dimnames(Xc) <- list(
    sprintf("s%02d", seq_len(n)),
    sprintf("c%02d", seq_len(p))
  )
  X <- Xc
  other <- setdiff(seq_len(p), c(5L, 6L, 12L))
  X[cbind(sample.int(n, 80L, TRUE), sample(other, 80L, TRUE))] <- NA_real_
  X[c(3L, 11L, 19L, 27L), 5L] <- NA_real_
  X[c(2L, 9L, 20L, 28L), 6L] <- NA_real_
  fit <- function(M, cols) {
    mipca(
      M,
      rank = 2L,
      transform = "identity",
      subset = cols,
      chains = 2L,
      iter_sampling = 6L,
      iter_warmup = 4L,
      seed = 1L
    )
  }
  list(X = X, Xc = Xc, fit = fit)
})

test_that("a narrowed chunk does not keep the fit it was finalized from alive", {
  f <- mipca(
    .nj$X,
    rank = 2L,
    transform = "identity",
    subset = "c05",
    keep_state = TRUE,
    chains = 2L,
    iter_sampling = 6L,
    iter_warmup = 4L,
    seed = 1L
  )
  gone <- FALSE
  reg.finalizer(f, function(e) gone <<- TRUE)
  g <- mipca_finalize(f)
  rm(f)
  gc()
  expect_true(gone)
  expect_false(anyNA(as.matrix(g, 1L)))
})

test_that("cbind() refuses a column that is already in the join", {
  a <- .nj$fit(.nj$X, "c05")
  b <- .nj$fit(.nj$X, "c06")
  expect_error(cbind(a, a), "repeat 1 column name [c05]", fixed = TRUE)
  expect_error(cbind(a, b, a), "repeat 1 column name [c05]", fixed = TRUE)
})

test_that("summary() numbers a joined fit's columns as as.matrix() does", {
  j <- cbind(
    .nj$fit(.nj$X, "c05"),
    .nj$fit(.nj$X, "c12"),
    .nj$fit(.nj$X, "c06")
  )
  cells <- summary(j)$cells
  expect_identical(sort(unique(cells$col)), c(1L, 3L))
  expect_identical(colnames(as.matrix(j, 1L))[cells$col], cells$name)
})

test_that("a request with nothing to impute returns a fit, and it joins", {
  X <- .nj$X

  whole <- .nj$fit(.nj$Xc[, 1:10], NULL)
  expect_identical(as.matrix(whole, 1L), .nj$Xc[, 1:10])

  clean <- .nj$fit(X, "c12")
  expect_identical(nrow(clean$chunks[[1L]]$map), 0L)
  expect_identical(as.matrix(clean, 1L), X[, "c12", drop = FALSE])

  a <- .nj$fit(X, "c05")
  b <- .nj$fit(X, "c06")
  r <- as.matrix(cbind(a, clean, b), 3L)
  Xt <- X[, c("c05", "c12", "c06")]
  expect_identical(sum(is.na(Xt)), 8L)
  expect_identical(dimnames(r), dimnames(Xt))
  expect_identical(r[!is.na(Xt)], Xt[!is.na(Xt)])
  expect_false(anyNA(r))
  expect_identical(r[, 1L], as.matrix(a, 3L)[, 1L])
  expect_identical(r[, 3L], as.matrix(b, 3L)[, 1L])
})

test_that("a one-draw part built with a double nkeep joins a sampler fit", {
  X <- .nj$X
  obj <- X[, "c05", drop = FALSE]
  map <- which(is.na(obj), arr.ind = TRUE)
  dimnames(map) <- list(NULL, c("row", "col"))
  expect_gt(nrow(map), 0L)
  ch <- .si_mi_chunk(
    obj = obj,
    map = map,
    draws = matrix(seq_len(nrow(map)) / 10, ncol = 1),
    transform = list(name = "identity", backward = function(z, ...) z),
    diag = .diag_frame(
      map = map[integer(0), , drop = FALSE],
      draws = list(matrix(NA_real_, 0L, 1L)),
      p = 1L,
      nms = "c05",
      prop_diag = 0.5
    ),
    prov = list(seam_state = "never", thin_effective = 1),
    call = NULL
  )
  point <- .si_mi_new(list(ch), rownames(X), nkeep = 1, call = NULL)
  expect_identical(point$nkeep, 1L)

  samp <- mipca(
    X,
    rank = 2L,
    transform = "identity",
    subset = "c06",
    chains = 1L,
    iter_sampling = 1L,
    iter_warmup = 4L,
    seed = 1L
  )
  j <- cbind(point, samp)
  expect_identical(
    as.matrix(j),
    cbind(as.matrix(point), as.matrix(samp))
  )
})
