#' Ampute NA Given the Output of `sample_each_rep()`
#'
#' Used in testthat only.
#'
#' @param obj Input.
#' @param loc Output of `sample_each_rep()`.
#'
#' @returns `obj` with an `NA` at every cell `loc` names.
#'
#' @noRd
#' @keywords internal
apply_na <- function(obj, loc) {
  obj[cbind(loc[, "row"], loc[, "col"])] <- NA_real_
  obj
}

#' Force pca_imp_internal_cpp() to run a fixed number of iterations.
#' Returns the completed matrix plus LOBPCG/exact counters.
#'
#' @noRd
run_pca_fixed_iters <- function(
  x,
  ctrl = lobpcg_control(),
  ncp = 2L,
  pca_iters = 10L,
  solver = c("exact", "lobpcg", "auto"),
  colmax = 1
) {
  ctrl <- new_lobpcg_control(ctrl, solver = solver)
  solver_code <- switch(solver, exact = 0L, lobpcg = 1L, auto = 2L)
  # mirror pca_imp() eligibility closely enough for fixed-iteration tests:
  # drop columns above colmax and zero-/undefined-variance columns.
  miss_rate <- mat_miss(x, col = TRUE, prop = TRUE)

  cv <- col_vars(x)

  eligible <- miss_rate < min(colmax, 1) &
    !(is.na(cv) | cv < .Machine$double.eps)

  eligible_idx <- as.integer(which(eligible) - 1L)

  if (length(eligible_idx) == 0L) {
    cli::cli_abort("No eligible columns for PCA imputation.")
  }

  cap <- min(nrow(x) - 2L, length(eligible_idx) - 1L)
  if (ncp > cap) {
    cli::cli_abort(
      "{.arg ncp} is too large for the eligible fixed-iteration test problem."
    )
  }

  res <- pca_imp_internal_cpp(
    obj = x,
    eligible_idx = eligible_idx,
    ncp = ncp,
    scale = TRUE,
    regularized = TRUE,
    threshold = 0,
    init = 0L,
    maxiter = pca_iters,
    miniter = pca_iters,
    row_w = rep(1, nrow(x)),
    coeff_ridge = 1,
    solver = solver_code,
    warmup_iters = ctrl$warmup_iters,
    lobpcg_tol = ctrl$tol,
    lobpcg_maxiter = ctrl$maxiter
  )

  x[res$loc] <- res$value

  list(
    mat = x,
    solver = solver,
    solver_chosen_code = res$solver_chosen,
    n_exact = res$n_exact,
    n_lobpcg_ok = res$n_lobpcg_ok,
    n_lobpcg_bad = res$n_lobpcg_bad
  )
}

skip_if_not_manual <- function() {
  skip_if_not(
    identical(Sys.getenv("MANUAL_TESTS"), "true"),
    "Set MANUAL_TESTS=true to run manually-triggered tests"
  )
}

# Asserts that `got` matches `want` cell by cell, on the worst case: the same
# dims and the same NA cells, then max(abs(got - want)) over the cells `want`
# holds a value at, and max(abs(got - want) / abs(want)) over those that are
# also nonzero, each at most its tolerance. Both figures are in the labels, so
# a failure reports them. At least one cell must be compared.
expect_worst_err <- function(got, want, abs_tol = 1e-12, rel_tol = 1e-12) {
  expect_identical(dim(got), dim(want))
  g <- as.vector(got)
  w <- as.vector(want)
  expect_identical(is.na(g), is.na(w))
  ok <- !is.na(w) & !is.na(g)
  expect_gt(sum(ok), 0L)
  d <- abs(g[ok] - w[ok])
  nz <- w[ok] != 0
  max_abs <- max(d)
  max_rel <- if (any(nz)) max(d[nz] / abs(w[ok][nz])) else 0
  label <- sprintf(
    "max abs error %.3g, max rel error %.3g (over %d nonzero cells)",
    max_abs,
    max_rel,
    sum(nz)
  )
  expect_lte(max_abs, abs_tol, label = label)
  expect_lte(max_rel, rel_tol, label = label)
}

kernel_iters <- 10L

kernel_eligible <- function(x) {
  miss <- mat_miss(x, col = TRUE, prop = TRUE)
  cv <- col_vars(x)
  as.integer(which(miss < 1 & !(is.na(cv) | cv < .Machine$double.eps)) - 1L)
}

# One kernel fit at a fixed iteration count. init > 0 draws its start from
# R's RNG, so the seed is set before every call and both solvers start alike.
kernel_fit <- function(
  x,
  ncp,
  solver,
  scale = TRUE,
  regularized = TRUE,
  init = 0L,
  iters = kernel_iters,
  row_w = rep(1, nrow(x)),
  warmup = 2L,
  tol = 1e-10,
  lobpcg_maxiter = 100L
) {
  set.seed(1)
  pca_imp_internal_cpp(
    obj = x,
    eligible_idx = kernel_eligible(x),
    ncp = ncp,
    scale = scale,
    regularized = regularized,
    threshold = 0,
    init = init,
    maxiter = iters,
    miniter = iters,
    row_w = row_w,
    coeff_ridge = 1,
    solver = switch(solver, exact = 0L, lobpcg = 1L),
    warmup_iters = warmup,
    lobpcg_tol = tol,
    lobpcg_maxiter = lobpcg_maxiter
  )
}

# A named 12 x 8 input where K-NN leaves cells, at k = 1. c3 is c2 plus a
# little noise, so each is the other's nearest neighbor, and both are missing
# in row 1, so K-NN leaves row 1 in each. It fills row 5 of c2, where c3 is
# observed, and both holes of c5. c7 is over colmax; c1, c4, c6 and c8 are
# complete.
knn_left_input <- function() {
  set.seed(3)
  n <- 12
  obj <- matrix(
    rnorm(n * 8),
    n,
    8,
    dimnames = list(paste0("r", 1:n), paste0("c", 1:8))
  )
  obj[, 3] <- obj[, 2] + 1e-3 * rnorm(n)
  obj[c(1, 5), 2] <- NA
  obj[1, 3] <- NA
  obj[c(2, 7), 5] <- NA
  obj[2:n, 7] <- NA
  obj
}

# The values, dims and dimnames only, without the class or the record.
bare_knn <- function(x) {
  attributes(x) <- list(dim = dim(x), dimnames = dimnames(x))
  x
}
