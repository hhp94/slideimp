#' Ampute NA Given the Output of `sample_each_rep()`
#'
#' Used in testthat only.
#'
#' @param obj Input.
#' @param loc Output of `sample_each_rep()`.
#'
#' @returns `NULL` invisibly. Called for its side effects (messages or errors).
#'
#' @noRd
#' @keywords internal
apply_na <- function(obj, loc) {
  obj[cbind(loc[, "row"], loc[, "col"])] <- NA_real_
  obj
}

#' Force pca_imp() to Run a Fixed Number of Iterations
#'
#' Used in testthat only.
#'
#' @noRd
#' @keywords internal
run_pca_fixed_iters <- function(x, ctrl, ncp = 3L, pca_iters = 14L) {
  pca_imp(
    x,
    ncp = ncp,
    scale = TRUE,
    method = "regularized",
    coeff.ridge = 1,
    threshold = 0,
    maxiter = pca_iters,
    miniter = pca_iters,
    nb.init = 1,
    seed = 1,
    lobpcg_control = ctrl,
    colmax = 1,
    post_imp = FALSE,
    na_check = TRUE
  )
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
  ctrl <- new_lobpcg_control(
    ctrl,
    ncp = ncp,
    n = nrow(x),
    p = ncol(x),
    solver = solver
  )
  solver_code <- switch(solver,
    exact = 0L,
    lobpcg = 1L,
    auto = 2L
  )
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

  iv <- res$imputed_values

  if (!is.null(iv) && nrow(iv) > 0L) {
    x[cbind(as.integer(iv[, 1]), as.integer(iv[, 2]))] <- iv[, 3]
  }

  list(
    mat = x,
    solver = solver,
    solver_chosen_code = res$solver_chosen,
    n_exact = res$n_exact,
    n_lobpcg_ok = res$n_lobpcg_ok,
    n_lobpcg_bad = res$n_lobpcg_bad
  )
}

max_abs_diff <- function(a, b) {
  stopifnot(!anyNA(a), !anyNA(b))
  max(abs(a - b))
}
