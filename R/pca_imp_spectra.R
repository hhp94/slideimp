#' Spectra Eigensolver Backend for PCA Imputation (Dev Bench)
#'
#' @description
#' Experimental developer benchmark of a Spectra (implicitly restarted
#' Lanczos, via `RSpectra`'s bundled headers) eigensolver backend.
#'
#' Runs a plain (non-incremental) version of the [pca_imp()] EM loop with the
#' eigensolver swapped to Spectra: impute, restandardize, weighted Gram,
#' Spectra top-k eig (warm-started with the previous dominant eigenvector),
#' reconstruct, objective. The loss formula and convergence criterion are the
#' same as [pca_imp()]'s backend, so `n_iter` and the final `objective`
#' should roughly match; imputed values are not returned.
#'
#' Build with `load_all1(timer = TRUE)` and compare the rcpptimer block
#' summaries: `pca_imp_spectra` (this function) against `pca_imp_gram` from
#' `pca_imp(..., solver = "lobpcg")` — in particular the `eig` and
#' `form_gram` blocks.
#'
#' @inheritParams pca_imp
#' @param spectra_tol Numeric. Spectra convergence tolerance.
#' @param spectra_maxiter Integer. Maximum Spectra (restart) iterations per
#'   solve.
#' @param spectra_ncv `NULL` or integer. Krylov subspace dimension; `NULL`
#'   uses RSpectra's default `min(n, max(2 * nev + 1, 20))`.
#'
#' @returns A list with `n_iter`, `objective`, `criterion_final`,
#' `converged`, the final top eigenvalues, and per-EM-iteration Spectra
#' diagnostics (`eig_niter`, `eig_nops`, `eig_max_rel_res`).
#'
#' @keywords internal
#' @export
pca_imp_spectra <- function(
  obj,
  ncp = 2,
  scale = TRUE,
  method = c("regularized", "EM"),
  coeff.ridge = 1,
  row.w = NULL,
  threshold = 1e-6,
  maxiter = 1000,
  miniter = 5,
  spectra_tol = 1e-9,
  spectra_maxiter = 1000,
  spectra_ncv = NULL,
  colmax = 0.9
) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    null.ok = FALSE,
    .var.name = "obj"
  )
  check_finite(obj)
  checkmate::assert_int(
    ncp,
    lower = 1L,
    upper = ncol(obj) - 1L,
    .var.name = "ncp"
  )
  method <- match.arg(method)
  checkmate::assert_flag(scale, .var.name = "scale")
  checkmate::assert_number(coeff.ridge, lower = 0, .var.name = "coeff.ridge")
  checkmate::assert(
    checkmate::check_numeric(
      row.w,
      finite = TRUE,
      any.missing = FALSE,
      len = nrow(obj),
      lower = 1e-10
    ),
    checkmate::check_choice(row.w, choices = "n_miss"),
    checkmate::check_null(row.w),
    .var.name = "row.w"
  )
  checkmate::assert_number(threshold, lower = 0, .var.name = "threshold")
  checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
  checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")
  checkmate::assert_number(
    spectra_tol,
    lower = 0,
    finite = TRUE,
    .var.name = "spectra_tol"
  )
  checkmate::assert_int(spectra_maxiter, lower = 1, .var.name = "spectra_maxiter")
  checkmate::assert_int(
    spectra_ncv,
    lower = 1,
    null.ok = TRUE,
    .var.name = "spectra_ncv"
  )
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")

  # eligibility: same policy as pca_imp().
  cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)
  miss_rate <- cmiss / nrow(obj)
  obj_vars <- col_vars(obj)
  eligible <- miss_rate <= min(colmax, 1) &
    !(obj_vars < .Machine$double.eps | is.na(obj_vars))

  n_elig <- sum(eligible)
  cap <- min(nrow(obj) - 2L, n_elig - 1L)
  if (cap < 1L || ncp > cap) {
    cli::cli_abort(
      "PCA is infeasible: ncp ({ncp}) vs cap ({cap}) with {n_elig} eligible columns.",
      class = "slideimp_infeasible"
    )
  }

  eligible_idx <- which(eligible)

  if (is.null(row.w)) {
    row.w <- rep(1, nrow(obj))
  } else if (is.character(row.w) && row.w == "n_miss") {
    n_miss_per_row <- mat_miss(obj, col = FALSE, prop = FALSE)
    row.w <- 1 - (n_miss_per_row / ncol(obj))
    row.w[row.w < 1e-10] <- 1e-10
  }

  pca_imp_spectra_cpp(
    obj = obj,
    eligible_idx = as.integer(eligible_idx - 1L),
    ncp = ncp,
    scale = scale,
    regularized = (method == "regularized"),
    threshold = threshold,
    maxiter = maxiter,
    miniter = miniter,
    row_w = row.w,
    coeff_ridge = coeff.ridge,
    spectra_tol = spectra_tol,
    spectra_maxiter = spectra_maxiter,
    spectra_ncv = if (is.null(spectra_ncv)) 0L else as.integer(spectra_ncv)
  )
}
