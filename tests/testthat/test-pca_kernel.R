skip_on_cran()

test_that("wide recover_factors block agrees with the stored reference", {
  # SVD_triplet()'s wide branch forms V = Xhat' * Vn only for the changing
  # (NA-containing) columns inside the EM loop, and the fixed columns' rows
  # once after it (src/svd_triplet.h, src/pca_kernel.cpp postprocessing). The
  # fixed rows reach R only through `mse`, the nb.init selection statistic;
  # the imputed values never read them. So both are compared, for a wide
  # input with fixed columns (Gram cache active), a wide input with an NA in
  # every column (no fixed block, cache inactive), and a tall input, whose
  # path is unchanged.
  #
  # The reference was recorded from the code before the change, at 96dec4f,
  # by a generator script committed alongside it in 7d8508c
  # (`git show 7d8508c:dev/gen_pca_recover_block_ref.R`). It is bytes from
  # another build: the comparison is on a tolerance, and bitwise agreement is
  # reported, not asserted.
  ref <- readRDS(test_path("fixtures", "pca_recover_block_ref.rds"))
  claims <- list(
    wide_fixed = list(wide = TRUE, fixed = TRUE),
    wide_all = list(wide = TRUE, fixed = FALSE),
    tall_fixed = list(wide = FALSE, fixed = TRUE)
  )
  expect_setequal(names(ref$cases), names(claims))

  for (nm in names(claims)) {
    case <- ref$cases[[nm]]
    x <- case$x
    n_fixed <- sum(colSums(is.na(x)) == 0L)

    # the input exercises the branch the case is named for
    expect_true(anyNA(x), label = nm)
    expect_identical(ncol(x) > nrow(x), claims[[nm]]$wide, label = nm)
    expect_identical(n_fixed > 0L, claims[[nm]]$fixed, label = nm)

    got <- pca_imp_internal_cpp(
      obj = x,
      eligible_idx = seq_len(ncol(x)) - 1L,
      ncp = ref$ncp,
      scale = TRUE,
      regularized = TRUE,
      threshold = 0,
      init = 0L,
      maxiter = ref$iters,
      miniter = ref$iters,
      row_w = case$row_w,
      coeff_ridge = 1,
      solver = 0L,
      warmup_iters = 0L,
      lobpcg_tol = 1e-8,
      lobpcg_maxiter = 1L
    )
    want <- case$imputed_values

    # structural: the same cells, in the same order. The reference holds
    # each cell as a row of doubles (row, col, value), as the kernel returned
    # them when it was recorded; the kernel now returns the two indices in
    # the integer matrix `loc` and the values in `value`.
    expect_identical(dim(got$loc), c(nrow(want), 2L), label = nm)
    expect_null(dimnames(got$loc), label = nm)
    expect_identical(length(got$value), nrow(want), label = nm)
    expect_identical(unname(got$loc[, 1]), as.integer(want[, 1]), label = nm)
    expect_identical(unname(got$loc[, 2]), as.integer(want[, 2]), label = nm)

    a <- got$value
    b <- want[, 3]
    expect_gt(length(b), 0L)
    err <- abs(a - b)
    # the relative figure is taken over the imputed cells whose reference
    # value is not near zero
    keep <- abs(b) >= 1e-3
    expect_gt(sum(keep), 0L)
    max_abs <- max(err)
    max_rel <- max(err[keep] / abs(b[keep]))
    expect_gt(case$mse, 0)
    mse_abs <- abs(got$mse - case$mse)
    mse_rel <- mse_abs / abs(case$mse)

    cat(sprintf(
      paste0(
        "\n%s: imputed max abs %.3g, max rel %.3g over %d of %d cells; ",
        "mse abs %.3g, rel %.3g; bitwise imputed %s, mse %s\n"
      ),
      nm,
      max_abs,
      max_rel,
      sum(keep),
      length(b),
      mse_abs,
      mse_rel,
      identical(a, b),
      identical(got$mse, case$mse)
    ))

    expect_lt(max_abs, 1e-8, label = nm)
    expect_lt(max_rel, 1e-8, label = nm)
    expect_lt(mse_abs, 1e-8, label = nm)
    expect_lt(mse_rel, 1e-8, label = nm)
  }
})

test_that("the PCA kernel returns integers and refuses what it cannot count with", {
  set.seed(1234)
  x <- sim_mat(20, 12, perc_total_na = 0.1, perc_col_na = 1)$input
  holes <- which(is.na(x), arr.ind = TRUE)
  expect_gt(nrow(holes), 0L)
  args <- list(
    obj = x,
    eligible_idx = seq_len(ncol(x)) - 1L,
    ncp = 2L,
    scale = TRUE,
    regularized = TRUE,
    threshold = 1e-6,
    init = 0L,
    maxiter = 20L,
    miniter = 5L,
    row_w = rep(1, nrow(x)),
    coeff_ridge = 1,
    solver = 0L,
    warmup_iters = 0L,
    lobpcg_tol = 1e-8,
    lobpcg_maxiter = 1L
  )
  got <- do.call(pca_imp_internal_cpp, args)

  # `loc` has one row per hole, 1-based, in column-major order, which is the
  # order which() lists them in. It has no dimnames: pca_imp() sets them in
  # R, since an Rcpp dimnames assignment is an rchk error in Rcpp's header.
  expect_true(is.integer(got$loc))
  expect_null(dimnames(got$loc))
  expect_identical(got$loc[, 1L], unname(holes[, 1]))
  expect_identical(got$loc[, 2L], unname(holes[, 2]))
  expect_true(is.double(got$value))
  expect_null(dim(got$value))
  expect_length(got$value, nrow(holes))
  expect_false(anyNA(got$value))
  expect_true(is.integer(got$n_iter))
  expect_gte(got$n_iter, 5L)

  # a count below its range, NA included; as an unsigned index each of
  # these wrapped to a huge value
  bad <- function(...) {
    do.call(pca_imp_internal_cpp, utils::modifyList(args, list(...)))
  }
  expect_error(bad(maxiter = 0L), "'maxiter' must be >= ", fixed = TRUE)
  expect_error(
    bad(maxiter = NA_integer_),
    "'maxiter' must be >= ",
    fixed = TRUE
  )
  expect_error(bad(miniter = -1L), "'miniter' must be >= ", fixed = TRUE)
  expect_error(bad(miniter = 21L), "'miniter' must be <= ", fixed = TRUE)
  expect_error(bad(ncp = 0L), "'ncp' must be >= ", fixed = TRUE)
  expect_error(bad(ncp = NA_integer_), "'ncp' must be >= ", fixed = TRUE)
  expect_error(bad(init = -1L), "'init' must be >= ", fixed = TRUE)
  expect_error(bad(trace_iter = -1L), "'trace_iter' must be >= ", fixed = TRUE)
  expect_error(
    bad(solver = 1L, lobpcg_maxiter = 0L),
    "'lobpcg_maxiter' must be >= 1 when",
    fixed = TRUE
  )
  expect_error(
    bad(solver = 1L, warmup_iters = -1L),
    "'warmup_iters' must be >= 0 when",
    fixed = TRUE
  )

  # a column index outside the matrix, NA included
  expect_error(
    bad(eligible_idx = c(0L, ncol(x))),
    "eligible_idx' holds 12 at position 2",
    fixed = TRUE
  )
  expect_error(
    bad(eligible_idx = c(0L, NA_integer_)),
    "eligible_idx' holds NA at position 2",
    fixed = TRUE
  )
  expect_error(
    bad(eligible_idx = c(-1L, 0L)),
    "eligible_idx' holds -1 at position 1",
    fixed = TRUE
  )
  expect_error(
    bad(row_w = rep(1, nrow(x) - 1L)),
    "'row_w' length must equal",
    fixed = TRUE
  )
})

test_that("a LOBPCG run that never converges leaves the exact fit untouched", {
  # A tolerance below what double precision reaches, or 2 iterations, fails
  # every LOBPCG solve after the warm-up; each falls back to the exact
  # solver from the same state, so the fit is the exact one. Both pass
  # through BLAS, which need not return the same bits for the same call
  # twice, so the values are asserted to rounding and bitwise reported.
  set.seed(1)
  x <- sim_mat(n = 60, p = 300, perc_total_na = 0.1, perc_col_na = 0.5)$input
  ex <- kernel_fit(x, 3L, "exact")
  # the exact fit it is held to imputes something, and nothing it imputes is
  # missing or infinite
  expect_gt(length(ex$value), 0L)
  expect_true(all(is.finite(ex$value)))
  for (ctrl in list(
    list(tol = 1e-14, lobpcg_maxiter = 40L),
    list(lobpcg_maxiter = 2L)
  )) {
    lo <- do.call(kernel_fit, c(list(x, 3L, "lobpcg"), ctrl))
    expect_identical(lo$n_lobpcg_ok, 0L)
    expect_identical(lo$n_lobpcg_bad, kernel_iters - 2L)
    expect_identical(lo$n_exact, kernel_iters)
    expect_identical(lo$loc, ex$loc)
    expect_worst_err(lo$value, ex$value)
    message(sprintf(
      "failed LOBPCG against exact bitwise: %s",
      identical(lo$value, ex$value)
    ))
  }
})
