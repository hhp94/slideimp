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
#' Returns the completed matrix plus LOBPCG/dsyevr counters.
#'
#' @noRd
#' @noRd
run_pca_fixed_iters <- function(x, ctrl, ncp = 3L, pca_iters = 14L) {
  # mirror pca_imp() eligibility: drop all-NA columns and zero-variance columns

  miss_rate <- mat_miss(x, prop = TRUE)
  cv <- col_vars(x)
  eligible <- miss_rate < 1 & !(is.na(cv) | cv < .Machine$double.eps)
  eligible_idx <- which(eligible) - 1L

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
    warmup_iters = ctrl$warmup_iters,
    lobpcg_tol = ctrl$tol,
    lobpcg_maxiter = ctrl$maxiter
  )
  iv <- res$imputed_values
  x[cbind(iv[, 1], iv[, 2])] <- iv[, 3]
  list(
    mat = x,
    n_dsyevr = res$n_dsyevr,
    n_lobpcg_ok = res$n_lobpcg_ok,
    n_lobpcg_bad = res$n_lobpcg_bad
  )
}

max_abs_diff <- function(a, b) {
  stopifnot(!anyNA(a), !anyNA(b))
  max(abs(a - b))
}
