skip_if_not_manual()

# The PCA kernel, pca_imp_internal_cpp(), on the input families the gate
# scripts of the Eigen conversion fed it by hand (dev/eigen_p2s4_edges.R and
# dev/eigen_p2s6_diag.R at d5a34e5): many eigenpairs under LOBPCG, ncp = 1,
# EM and regularized fits, no scaling, row weights, a random start, a hole in
# every column, one iteration, samples and columns that are copies with more
# components than rank, LOBPCG settings that start it at once, and a seeded
# random sweep over all of them. The settings that stall LOBPCG are in
# test-pca_kernel.R, which runs by default. Those gates also made a QR, a
# triangular solve or an eigensolver call fail on purpose; that needs an
# instrumented build, so it is not carried here.
#
# The oracle is the other solver: each input is fitted once under the exact
# solver and once under LOBPCG, over the same fixed number of EM iterations,
# so both take the same steps and any difference is the eigensolvers'. Over
# every fit in this file the two differ by at most 5.22e-9 absolute and
# 1.97e-6 relative (measured by printing kernel_err() for each), and a fit
# with one component fewer differs from the full one by 5.65e-2 or more on
# the first block's inputs, so the tolerances below sit between what the
# solvers leave and what a wrong answer looks like.

kernel_abs_tol <- 2e-7
kernel_rel_tol <- 1e-5

# Worst absolute and worst relative difference between two fits' imputed
# values, the relative one over the cells whose `want` value is above 1e-12
# in magnitude.
kernel_err <- function(got, want) {
  d <- abs(got$value - want$value)
  nz <- abs(want$value) > 1e-12
  c(
    abs = max(d),
    rel = if (any(nz)) max(d[nz] / abs(want$value[nz])) else 0
  )
}

# Fits `x` under both solvers and checks what must hold whatever the input:
# the same cells imputed, every value finite, the input left as it was, the
# solver counts adding up, and the two solvers agreeing on the worst cell.
# Returns both fits for the caller's own checks.
expect_solvers_agree <- function(case, x, ncp, ...) {
  # Force every argument before either fit. kernel_fit() sets the seed, and
  # an argument still unevaluated there would draw from the RNG after it in
  # the first fit only, so the two fits would start from different random
  # starts under init > 0.
  force(ncp)
  args <- list(...)
  x0 <- x + 0
  ex <- do.call(kernel_fit, c(list(x, ncp, "exact"), args))
  lo <- do.call(kernel_fit, c(list(x, ncp, "lobpcg"), args))
  iters <- args$iters %||% kernel_iters

  expect_identical(x, x0, label = paste(case, "input"))
  expect_gt(length(ex$value), 0L, label = paste(case, "imputed cells"))
  expect_identical(lo$loc, ex$loc, label = paste(case, "imputed cells"))
  expect_true(all(is.finite(ex$value)), label = paste(case, "exact finite"))
  expect_true(all(is.finite(lo$value)), label = paste(case, "LOBPCG finite"))
  expect_identical(ex$n_exact, iters, label = paste(case, "exact count"))
  # every iteration is solved once, by LOBPCG or by the exact solver; a
  # failed LOBPCG solve is counted bad and its exact fallback counted exact
  expect_identical(
    lo$n_exact + lo$n_lobpcg_ok,
    iters,
    label = paste(case, "LOBPCG counts")
  )
  expect_lte(lo$n_lobpcg_bad, lo$n_exact, label = paste(case, "fallbacks"))

  e <- kernel_err(lo, ex)
  label <- sprintf(
    "%s: max abs error %.3g, max rel error %.3g",
    case,
    e[["abs"]],
    e[["rel"]]
  )
  expect_lte(e[["abs"]], kernel_abs_tol, label = label)
  expect_lte(e[["rel"]], kernel_rel_tol, label = label)
  list(exact = ex, lobpcg = lo)
}

test_that("the solvers agree on every input family the gates fed the kernel", {
  set.seed(1)
  wide <- sim_mat(n = 60, p = 300, perc_total_na = 0.1, perc_col_na = 0.5)$input
  set.seed(2)
  tall <- sim_mat(n = 300, p = 60, perc_total_na = 0.1, perc_col_na = 0.5)$input
  set.seed(3)
  holes <- sim_mat(n = 80, p = 120, perc_total_na = 0.1, perc_col_na = 1)$input
  set.seed(4)
  w <- runif(nrow(holes), 0.1, 1)
  set.seed(11)
  big <- sim_mat(n = 120, p = 160, perc_total_na = 0.1, perc_col_na = 0.5)$input
  # premise: a hole in every column
  expect_true(all(colSums(is.na(holes)) > 0L))

  cases <- list(
    list("wide ncp 1", wide, 1L),
    list("wide ncp 2", wide, 2L),
    list("wide ncp 5", wide, 5L),
    list("wide ncp 11", wide, 11L),
    list("wide ncp 17", wide, 17L),
    list("wide ncp 17 EM", wide, 17L, regularized = FALSE),
    list("tall ncp 1", tall, 1L),
    list("tall ncp 11", tall, 11L),
    list("holes ncp 3", holes, 3L),
    list("holes ncp 20", holes, 20L),
    list("holes ncp 3 row weights", holes, 3L, row_w = w),
    list("holes ncp 20 row weights", holes, 20L, row_w = w),
    list("holes ncp 3 no scaling", holes, 3L, scale = FALSE),
    list("holes ncp 20 no scaling", holes, 20L, scale = FALSE),
    list("holes ncp 3 random start", holes, 3L, init = 1L),
    list(
      "holes ncp 20 random start EM",
      holes,
      20L,
      init = 1L,
      regularized = FALSE
    ),
    list("holes one iteration", holes, 3L, iters = 1L, warmup = 0L),
    list("wide no warm-up", wide, 3L, warmup = 0L),
    list("wide warm-up 1", wide, 3L, warmup = 1L),
    list("big ncp 16", big, 16L),
    list("big ncp 25", big, 25L),
    list("big ncp 40", big, 40L),
    list("big ncp 33 EM", big, 33L, regularized = FALSE)
  )

  lobpcg_ran <- character()
  control <- numeric()
  for (cs in cases) {
    args <- cs[-(1:3)]
    fits <- do.call(expect_solvers_agree, c(cs[1:3], args))
    if (fits$lobpcg$n_lobpcg_ok > 0L) {
      lobpcg_ran <- c(lobpcg_ran, cs[[1]])
    }
    # what a wrong answer looks like: one component fewer, exact solver
    if (cs[[3]] >= 2L && (args$iters %||% kernel_iters) > 1L) {
      fewer <- do.call(
        kernel_fit,
        c(list(cs[[2]], cs[[3]] - 1L, "exact"), args)
      )
      control[[cs[[1]]]] <- kernel_err(fewer, fits$exact)[["abs"]]
    }
  }

  # LOBPCG must actually have run, at many eigenpairs too, or the agreement
  # above says nothing about it.
  expect_gte(length(lobpcg_ran), 15L)
  expect_true(all(c("big ncp 25", "big ncp 40") %in% lobpcg_ran))
  # and the tolerance must be able to fail: every one-fewer fit is far
  # outside it.
  expect_gte(length(control), 15L)
  expect_gt(min(control), 1000 * kernel_abs_tol)
})

test_that("copies of samples or columns, with ncp above rank, fit and agree", {
  # 30 samples of which 6 are distinct, so the completed data has rank 6 at
  # most; and 12 columns that are 4 copies of 3. Up to 1.2.0 asking for more
  # components than that ended the R process.
  set.seed(21)
  base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  dup_rows <- base6[rep(seq_len(6), times = 5), ]
  set.seed(22)
  base3 <- matrix(rnorm(40 * 3), 40, 3)
  dup_cols <- base3[, rep(1:3, times = 4)]
  dup_cols[cbind(c(2, 9, 17, 30), c(1, 5, 8, 12))] <- NA
  # premises
  expect_identical(nrow(unique(dup_rows)), 6L)
  expect_true(anyNA(dup_rows))
  expect_identical(qr(base3)$rank, 3L)

  fits <- list(
    rows_5 = expect_solvers_agree("rows ncp 5", dup_rows, 5L, warmup = 1L),
    rows_7 = expect_solvers_agree("rows ncp 7", dup_rows, 7L, warmup = 1L),
    rows_12 = expect_solvers_agree("rows ncp 12", dup_rows, 12L, warmup = 1L),
    cols_5 = expect_solvers_agree("cols ncp 5", dup_cols, 5L, warmup = 1L),
    cols_8 = expect_solvers_agree("cols ncp 8", dup_cols, 8L, warmup = 1L)
  )
  fits$rows_7_em <- expect_solvers_agree(
    "rows ncp 7 EM",
    dup_rows,
    7L,
    warmup = 1L,
    regularized = FALSE
  )
  # LOBPCG must have run, or "agree" is the exact fit against itself: a failed
  # LOBPCG solve falls back to the exact one bit for bit. "cols ncp 8" is
  # left out: every iteration of it is exact, LOBPCG never attempted.
  for (nm in c("rows_5", "rows_7", "rows_12", "cols_5", "rows_7_em")) {
    expect_gt(fits[[nm]]$lobpcg$n_lobpcg_ok, 0L, label = nm)
  }

  # Components past the rank carry nothing: every copy of a sample has its
  # holes in the same cells and is imputed alike, so the completed data keeps
  # 6 distinct samples, rank 5 once centered, and two fits that both ask for
  # at least that impute the same values, under either solver. (Copies of a
  # column have their holes in different cells, so their completed copies
  # differ and no such rule holds.)
  for (solver in c("exact", "lobpcg")) {
    e <- kernel_err(fits$rows_12[[solver]], fits$rows_5[[solver]])
    expect_lte(e[["abs"]], 1e-10, label = paste("rows 12 vs 5", solver))
  }
  # and below the rank they do not, so the check above can fail
  rows_4 <- kernel_fit(dup_rows, 4L, "exact", warmup = 1L)
  e <- kernel_err(rows_4, fits$rows_5$exact)
  expect_gt(e[["abs"]], 1e-3)
})

test_that("the solvers agree over a seeded random sweep of inputs and options", {
  # Each seed draws a shape, a missingness, an ncp up to the kernel's cap and
  # every option. A failure names its seed, so it can be rerun alone.
  n_run <- 0L
  n_lobpcg <- 0L
  for (s in 1:24) {
    set.seed(1000 + s)
    n <- sample(20:120, 1)
    p <- sample(20:160, 1)
    total <- runif(1, 0.02, 0.3)
    col <- runif(1, 0.2, 1)
    # sim_mat() refuses a column share too small for the total; such a draw
    # is skipped, the same draw on every run.
    x <- tryCatch(
      sim_mat(n = n, p = p, perc_total_na = total, perc_col_na = col)$input,
      error = function(e) NULL
    )
    if (is.null(x)) {
      next
    }
    cap <- min(nrow(x) - 2L, length(kernel_eligible(x)) - 1L)
    if (cap < 1L) {
      next
    }
    # drawn here, in this order, under this seed
    ncp <- sample(seq_len(min(cap, 30L)), 1)
    scale <- runif(1) < 0.7
    regularized <- runif(1) < 0.7
    init <- sample(0:1, 1)
    row_w <- if (runif(1) < 0.5) rep(1, n) else runif(n, 0.1, 1)
    warmup <- sample(0:3, 1)
    fits <- expect_solvers_agree(
      sprintf(
        "seed %d (%d x %d, ncp %d, scale %s, regularized %s, init %d, warmup %d)",
        1000 + s,
        n,
        p,
        ncp,
        scale,
        regularized,
        init,
        warmup
      ),
      x,
      ncp,
      scale = scale,
      regularized = regularized,
      init = init,
      row_w = row_w,
      warmup = warmup
    )
    n_run <- n_run + 1L
    n_lobpcg <- n_lobpcg + (fits$lobpcg$n_lobpcg_ok > 0L)
  }
  # the sweep must not have skipped itself away
  expect_gte(n_run, 18L)
  # and LOBPCG must have run in most of it, or the agreement is mostly the
  # exact fit against itself (measured: in 18 of the 22 draws that run)
  expect_gte(n_lobpcg, 15L)
})
