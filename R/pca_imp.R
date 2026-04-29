#' Validate Clamp Bounds
#'
#' @param clamp `NULL` or a numeric vector of length 2.
#' @param arg Argument name used in error messages.
#'
#' @return `NULL` or an unnamed numeric vector of length 2.
#'
#' @keywords internal
#' @noRd
resolve_clamp <- function(clamp, arg = "clamp") {
  if (is.null(clamp)) {
    return(NULL)
  }

  if (length(clamp) == 1L) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a numeric vector of length 2: {.code c(lower, upper)}.",
      "i" = "For upper-bound-only clamping, use {.code c(-Inf, upper)}.",
      "i" = "For lower-bound-only clamping, use {.code c(lower, Inf)}."
    ))
  }

  if (length(clamp) != 2L) {
    cli::cli_abort(c(
      "{.arg {arg}} must be a numeric vector of length 2: {.code c(lower, upper)}.",
      "i" = "Use {.code NULL} for no clamping.",
      "i" = "Use {.code c(-Inf, upper)} or {.code c(lower, Inf)} for one-sided clamping."
    ))
  }

  if (anyNA(clamp)) {
    cli::cli_abort(c(
      "{.arg {arg}} cannot contain missing values.",
      "i" = "Use {.code -Inf} or {.code Inf} for an open bound.",
      "i" = "For example, use {.code c(-Inf, 1)} or {.code c(0, Inf)}, not {.code c(NA, 1)} or {.code c(0, NA)}."
    ))
  }

  checkmate::assert_numeric(
    clamp,
    len = 2L,
    any.missing = FALSE,
    finite = FALSE,
    null.ok = FALSE,
    .var.name = arg
  )

  if (clamp[1L] > clamp[2L]) {
    cli::cli_abort(
      "{.arg {arg}} must be ordered as {.code c(lower, upper)} with lower <= upper."
    )
  }

  return(unname(as.numeric(clamp)))
}

#' LOBPCG Eigensolver Control Options
#'
#' Construct a validated list of control options for the LOBPCG eigensolver
#' used by [pca_imp()]. Most users do not need to call this directly.
#'
#' @param warmup_iters Integer. Number of warm-up iterations before the main
#'   LOBPCG solve. Must be non-negative.
#' @param tol Numeric. Convergence tolerance for the LOBPCG eigensolver. Must
#'   be non-negative and finite.
#' @param maxiter Integer. Maximum number of LOBPCG iterations. Must be
#'   non-negative. In [pca_imp()], `maxiter` must be positive when
#'   `solver = "auto"` or `solver = "lobpcg"`. Use `solver = "exact"` to force
#'   the exact solver.
#'
#' @returns A named list of class `"slideimp_lobpcg_control"` containing
#' `warmup_iters`, `tol`, and `maxiter`.
#'
#' @examples
#' set.seed(123)
#' obj <- sim_mat(10, 10)$input
#'
#' # Use all defaults
#' lobpcg_control()
#'
#' # Override a single option
#' lobpcg_control(maxiter = 50)
#'
#' # Force the exact solver from pca_imp()
#' pca_imp(obj, ncp = 2, solver = "exact")
#'
#' # Pass directly to pca_imp()
#' pca_imp(obj, ncp = 2, lobpcg_control = lobpcg_control(tol = 1e-9))
#'
#' # Or use a named list
#' pca_imp(obj, ncp = 2, lobpcg_control = list(maxiter = 50))
#'
#' @export
lobpcg_control <- function(warmup_iters = 10L, tol = 1e-9, maxiter = 20) {
  explicit <- c(
    warmup_iters = !missing(warmup_iters),
    tol = !missing(tol),
    maxiter = !missing(maxiter)
  )

  checkmate::assert_int(
    warmup_iters,
    lower = 0L,
    .var.name = "warmup_iters"
  )
  checkmate::assert_number(
    tol,
    lower = 0,
    finite = TRUE,
    .var.name = "tol"
  )
  checkmate::assert_int(
    maxiter,
    lower = 0L,
    .var.name = "maxiter"
  )
  structure(
    list(
      warmup_iters = as.integer(warmup_iters),
      tol = as.numeric(tol),
      maxiter = as.integer(maxiter)
    ),
    class = "slideimp_lobpcg_control",
    explicit = explicit
  )
}

#' Validate a LOBPCG Control Object
#'
#' @description
#' Internal helper that accepts `NULL`, a `"slideimp_lobpcg_control"` object,
#' or a (partial) named list, and returns a fully validated control object.
#' Used by [pca_imp()] to normalize the `lobpcg_control` argument before
#' dispatching to the C++ backend.
#'
#' @param x `NULL`, a `"slideimp_lobpcg_control"` object, or a named list.
#'
#' @returns A `"slideimp_lobpcg_control"` object.
#'
#' @keywords internal
#' @noRd
new_lobpcg_control <- function(
  x,
  ncp = NULL,
  n = NULL,
  p = NULL,
  solver = c("auto", "exact", "lobpcg")
) {
  solver <- match.arg(solver)
  # force exact solver. This overrides any lobpcg_control supplied.
  if (solver == "exact") {
    return(lobpcg_control(maxiter = 0L))
  }
  # if no explicit LOBPCG control was supplied, use defaults.
  # for solver = "auto", maxiter must remain positive so the C++ backend can
  # actually probe the LOBPCG branch.
  if (is.null(x)) {
    return(lobpcg_control())
  }
  # validate an already created control object.
  if (inherits(x, "slideimp_lobpcg_control")) {
    out <- x
  } else {
    if (!is.list(x)) {
      cli::cli_abort(
        "{.arg lobpcg_control} must be {.code NULL}, a list, or created by {.fn lobpcg_control}."
      )
    }
    if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
      cli::cli_abort("{.arg lobpcg_control} must be a named list.")
    }
    allowed <- names(formals(lobpcg_control))
    unknown <- setdiff(names(x), allowed)
    if (length(unknown) > 0L) {
      cli::cli_abort(c(
        "{cli::qty(length(unknown))}Unknown LOBPCG control option{?s}: {fmt_trunc(unknown, 10)}.",
        "i" = "Allowed options are: {.arg {allowed}}."
      ))
    }
    out <- do.call(lobpcg_control, x)
  }
  control_names <- names(formals(lobpcg_control))
  explicit <- attr(out, "explicit", exact = TRUE)

  if (is.null(explicit)) {
    # old or manually constructed control objects: treat all fields as explicit
    # to avoid silently overriding user intent.
    explicit <- stats::setNames(rep(TRUE, length(control_names)), control_names)
  } else {
    missing_explicit <- setdiff(control_names, names(explicit))
    if (length(missing_explicit) > 0L) {
      explicit[missing_explicit] <- TRUE
    }
    explicit <- explicit[control_names]
  }

  attr(out, "explicit") <- explicit

  if (solver %in% c("auto", "lobpcg") && out$maxiter == 0L) {
    cli::cli_abort(c(
      "{.arg lobpcg_control$maxiter} is 0, but {.arg solver} is {.val {solver}}.",
      "i" = "Use {.code solver = 'exact'} for the exact solver, or set {.arg maxiter} > 0."
    ))
  }
  out
}

#' PCA Imputation for Numeric Matrices
#'
#' Impute missing values in a numeric matrix using regularized or
#' expectation-maximization PCA imputation.
#'
#' @details
#' This algorithm is based on `missMDA::imputePCA()` and is optimized for tall
#' or wide numeric matrices.
#'
#' @param obj A numeric matrix with samples in rows and features in columns.
#' @param ncp Integer. Number of principal components used to predict missing
#'   entries.
#' @param scale Logical. If `TRUE`, columns are scaled to unit variance.
#' @param method Character. PCA imputation method: either `"regularized"` or
#'   `"EM"`.
#' @param coeff.ridge Numeric. Ridge regularization coefficient. Only used when
#'   `method = "regularized"`. Values less than `1` regularize less, moving
#'   closer to EM PCA. Values greater than `1` regularize more, moving closer
#'   to mean imputation.
#' @param row.w Row weights, internally normalized to sum to `1`. Can be:
#'   * `NULL`: all rows are weighted equally.
#'   * A numeric vector of positive weights with length `nrow(obj)`.
#'   * `"n_miss"`: rows with more missing values receive lower weight.
#' @param threshold Numeric. Convergence threshold.
#' @param seed Integer, numeric, or `NULL`. Random seed for reproducibility.
#' @param nb.init Integer. Number of random initializations. The first
#'   initialization is always mean imputation.
#' @param maxiter Integer. Maximum number of iterations.
#' @param miniter Integer. Minimum number of iterations.
#' @param solver Character. Eigensolver selection. One of `"auto"`, `"exact"`,
#'   or `"lobpcg"`. `"exact"` uses the exact solver. `"lobpcg"` uses the
#'   iterative LOBPCG solver with exact fallback. `"auto"` performs a short
#'   timed probe and chooses LOBPCG only if it is clearly faster than the exact
#'   solver. When `nb.init > 1`, the auto choice from the first PCA initialization
#'   is reused for subsequent PCA initializations.
#' @param lobpcg_control A list of LOBPCG eigensolver control options, usually
#'   created by [lobpcg_control()]. A plain named list is also accepted. Ignored
#'   when `solver = "exact"`.
#' @param colmax Numeric scalar between `0` and `1`. Columns with a missing-data
#'   proportion greater than `colmax` are not imputed.
#' @param post_imp Logical. If `TRUE`, replace any remaining missing values
#'   with column means after PCA imputation.
#' @param clamp Optional numeric vector of length 2 giving lower and upper bounds
#'   for PCA-imputed values. Use `NULL` for no clamping. Use `c(0, 1)` for DNA
#'   methylation beta values. Use `c(lb, Inf)` for only lower bound clamping, or
#'   `c(-Inf, ub)` for only upper bound clamping. Clamping is applied only to
#'   values imputed by the PCA step, not to observed values.
#' @param na_check Logical. If `TRUE`, check whether the result still contains
#'   missing values.
#'
#' @returns A numeric matrix of the same dimensions as `obj`, with missing
#' values imputed. The returned object has class `slideimp_results`.
#'
#' @section Performance tips:
#' `pca_imp()` relies heavily on linear algebra. On Windows, the default BLAS
#' shipped with R may be slow for large matrices. Advanced users can replace
#' it with [OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).
#'
#' PCA imputation speed depends on the eigensolver selected by `solver` and the
#' convergence threshold `threshold`. The exact solver is selected with
#' `solver = "exact"`. The iterative LOBPCG solver is selected with
#' `solver = "lobpcg"`. The default, `solver = "auto"`, performs a short timed
#' probe chooses LOBPCG only when it is clearly faster.
#'
#' For large or approximately low-rank genomic matrices, it can be useful to
#' benchmark `solver = "exact"` against `solver = "lobpcg"` on a representative
#' subset, such as chromosome 22, before tuning accuracy-related parameters such
#' as `ncp`, `coeff.ridge`, `window_size`, or `overlap_size`.
#'
#' The default `threshold = 1e-6` is conservative. In many genomic datasets,
#' `threshold = 1e-5` can be faster while giving very similar imputed values.
#' Check this on a representative subset before using the relaxed threshold in a
#' full analysis.
#'
#' See the pkgdown article
#' [Speeding up PCA imputation](https://hhp94.github.io/slideimp/articles/speeding-up-pca-imputation.html)
#' for a full workflow.
#'
#' @references
#' Josse J, Husson F (2013). Handling missing values in exploratory
#' multivariate data analysis methods. *Journal de la SFdS*, 153(2), 79-99.
#'
#' Josse J, Husson F (2016). missMDA: A Package for Handling Missing Values in
#' Multivariate Data Analysis. *Journal of Statistical Software*, 70(1), 1-31.
#' \doi{10.18637/jss.v070.i01}
#'
#' The PCA imputation algorithm is based on the original `missMDA`
#' implementation by Francois Husson and Julie Josse.
#'
#' @examples
#' set.seed(123)
#' obj <- sim_mat(10, 10)$input
#' sum(is.na(obj))
#' obj[1:4, 1:4]
#'
#' # Randomly initialize missing values 5 times. The first initialization is mean imputation.
#' pca_imp(obj, ncp = 2, nb.init = 5, seed = 123)
#'
#' @export
pca_imp <- function(
  obj,
  ncp = 2,
  scale = TRUE,
  method = c("regularized", "EM"),
  coeff.ridge = 1,
  row.w = NULL,
  threshold = 1e-6,
  seed = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  solver = c("auto", "exact", "lobpcg"),
  lobpcg_control = NULL,
  colmax = 0.9,
  post_imp = TRUE,
  na_check = TRUE,
  clamp = NULL
) {
  # pre-conditioning
  checkmate::assert_matrix(obj, mode = "numeric", null.ok = FALSE, .var.name = "obj")
  check_finite(obj)
  checkmate::assert_int(ncp, lower = 1L, upper = ncol(obj) - 1L, .var.name = "ncp")
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
  checkmate::assert_int(nb.init, lower = 1, .var.name = "nb.init")
  checkmate::assert_int(seed, null.ok = TRUE, lower = 0, .var.name = "seed")
  if (!is.null(seed) && nb.init > 1L && seed > 2147483647L / (nb.init - 1L)) {
    stop("`seed` too large")
  }
  checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
  checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")
  # solver resolves
  solver <- match.arg(solver)
  lobpcg_control <- new_lobpcg_control(
    lobpcg_control,
    ncp = ncp,
    n = nrow(obj),
    p = ncol(obj),
    solver = solver
  )

  lobpcg_explicit <- attr(lobpcg_control, "explicit", exact = TRUE)
  warmup_explicit <- isTRUE(lobpcg_explicit[["warmup_iters"]])

  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  checkmate::assert_flag(na_check, .var.name = "na_check")
  clamp <- resolve_clamp(clamp, arg = "clamp")

  # per-column missingnesss
  cmiss <- mat_miss(obj, col = TRUE, prop = FALSE)

  # early exit if no missingness at all
  if (!any(cmiss > 0L)) {
    cli::cli_inform("No missing values in input. Returning input unchanged.")
    return(obj)
  }

  miss_rate <- cmiss / nrow(obj)

  # eligibility: below colmax and non-degenerate variance
  obj_vars <- col_vars(obj)
  eligible <- miss_rate <= min(colmax, 1) &
    !(obj_vars < .Machine$double.eps | is.na(obj_vars))

  n_elig <- sum(eligible)

  max_ncp_rows <- nrow(obj) - 2L
  max_ncp_cols <- n_elig - 1L
  cap <- min(max_ncp_rows, max_ncp_cols)

  if (cap < 1L) {
    cli::cli_abort(
      c(
        "PCA imputation is infeasible with the current data and settings.",
        "i" = "Number of rows: {nrow(obj)}.",
        "i" = "Number of eligible columns: {n_elig}.",
        "i" = "Relax {.arg colmax}, remove zero-variance columns, or use more rows/columns."
      ),
      class = "slideimp_infeasible"
    )
  }

  if (ncp > cap) {
    row_bound <- max_ncp_rows <= max_ncp_cols

    cli::cli_abort(
      c(
        "{.arg ncp} ({ncp}) exceeds the maximum usable components ({cap}).",
        "i" = if (row_bound) {
          "Limited by rows ({nrow(obj)}). Reduce {.arg ncp} or add more rows."
        } else {
          "Limited by eligible columns ({n_elig}). Reduce {.arg ncp} or relax {.arg colmax}."
        }
      ),
      class = "slideimp_infeasible"
    )
  }

  # at least one eligible column must have missingness
  if (!any(cmiss[eligible] > 0L)) {
    cli::cli_abort(
      c(
        "All columns with missing values are ineligible (exceed {.arg colmax} ({colmax}) or have zero variance).",
        "i" = "Relax {.arg colmax} or remove zero-variance columns."
      ),
      class = "slideimp_infeasible"
    )
  }

  # solver = "auto" policy.
  #
  # the eigensystem dimension is the Gram dimension used by the backend.
  # LOBPCG is usually not worth probing for small eigensystems or when the
  # requested rank is a large fraction of the eigensystem dimension.
  k_eig <- as.integer(ncp + if (method == "regularized") 1L else 0L)
  n_gram <- min(nrow(obj), n_elig)

  auto_force_exact <- solver == "auto" &&
    (n_gram < 250L || (k_eig / n_gram) > 0.10)

  if (solver == "auto" && !auto_force_exact && !warmup_explicit) {
    lobpcg_control$warmup_iters <- as.integer(
      min(50L, max(10L, ceiling(1.5 * k_eig)))
    )
  }

  solver_code <- if (auto_force_exact) {
    0L
  } else {
    switch(solver,
      exact = 0L,
      lobpcg = 1L,
      auto = 2L
    )
  }

  # index of eligible columns (1-based)
  eligible_idx <- which(eligible)

  # row weights
  if (is.null(row.w)) {
    row.w <- rep(1, nrow(obj))
  } else if (is.character(row.w) && row.w == "n_miss") {
    # per-row missingness without materializing any boolean subset. We
    # purposefully scan the entire object here since for pca, we don't have the
    # subset argument and if missingness is spreadout through all columns,
    # obj[, eligible_idx] would materialize a full copy in RAM.
    n_miss_per_row <- mat_miss(obj, col = FALSE, prop = FALSE)
    row.w <- 1 - (n_miss_per_row / ncol(obj))
    row.w[row.w < 1e-10] <- 1e-10
  }

  init_obj <- Inf
  best_imputed <- NULL

  # solver_chosen codes:
  #   0 = forced exact
  #   1 = forced lobpcg/hybrid
  #   2 = auto had no reason/opportunity to choose. Treat as exact
  #   3 = auto chose exact
  #   4 = auto chose lobpcg
  locked_solver <- NULL
  resolved_solver_code <- if (isTRUE(auto_force_exact)) {
    0L
  } else {
    switch(solver,
      exact = 0L,
      lobpcg = 1L,
      auto = NA_integer_
    )
  }

  for (i in seq_len(nb.init)) {
    if (!is.null(seed)) {
      set.seed(seed * (i - 1L)) # exactly as missMDA does
    }

    solver_code_i <- if (is.null(locked_solver)) solver_code else locked_solver

    res.impute <- pca_imp_internal_cpp(
      obj = obj,
      eligible_idx = as.integer(eligible_idx - 1L),
      ncp = ncp,
      scale = scale,
      regularized = (method == "regularized"),
      threshold = threshold,
      init = if (i == 1L) 0L else i,
      maxiter = maxiter,
      miniter = miniter,
      row_w = row.w,
      coeff_ridge = coeff.ridge,
      solver = solver_code_i,
      warmup_iters = lobpcg_control$warmup_iters,
      lobpcg_tol = lobpcg_control$tol,
      lobpcg_maxiter = lobpcg_control$maxiter
    )

    # after the auto probe, select a solver for the rest of the inits.
    if (solver == "auto" && i == 1L) {
      chosen <- as.integer(res.impute$solver_chosen)
      if (length(chosen) != 1L || is.na(chosen)) {
        chosen <- 2L
      }
      locked_solver <- if (chosen %in% c(1L, 4L)) 1L else 0L
      resolved_solver_code <- locked_solver
    }

    cur_obj <- res.impute$mse
    if (cur_obj < init_obj) {
      # res.impute$imputed_values: N_miss x 3 matrix of
      # (row_idx, col_idx, value) in original-matrix, 1-based coordinates
      best_imputed <- res.impute$imputed_values
      init_obj <- cur_obj
    }
  }

  # column/row indices returned by backend are already original-matrix indices.
  if (is.null(best_imputed) || nrow(best_imputed) == 0L) {
    cli::cli_abort("Internal error: PCA imputation produced no imputed values.")
  }

  if (!is.null(clamp)) {
    best_imputed[, 3] <- pmin(
      pmax(best_imputed[, 3], clamp[1L]),
      clamp[2L]
    )
  }

  imp_indices <- cbind(
    as.integer(best_imputed[, 1]),
    as.integer(best_imputed[, 2])
  )
  obj[imp_indices] <- best_imputed[, 3]

  if (post_imp) {
    obj <- mean_imp_col(obj)
  }

  solver_chosen <- if (resolved_solver_code == 1L) "lobpcg" else "exact"

  out <- new_slideimp_results(
    obj,
    "pca",
    fallback = FALSE,
    post_imp = post_imp,
    na_check = na_check
  )
  attr(out, "solver_requested") <- solver
  attr(out, "solver_chosen") <- solver_chosen
  return(out)
}
