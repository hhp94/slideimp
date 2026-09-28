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

  check_lobpcg_fields(warmup_iters, tol, maxiter)

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

#' Validate LOBPCG Control Field Values
#'
#' @description
#' Internal helper holding the field rules for a LOBPCG control object, so that
#' [lobpcg_control()] and `new_lobpcg_control()` cannot drift apart. The latter
#' re-checks the fields of an object that already carries the class, because a
#' control object is an ordinary list and nothing stops a caller from assigning
#' into it after construction.
#'
#' @param warmup_iters,tol,maxiter The three control fields.
#' @param prefix Character prepended to each `.var.name`, so the message names
#'   the argument the user actually passed.
#'
#' @returns `NULL`, invisibly. Called for the side effect of aborting.
#'
#' @keywords internal
#' @noRd
check_lobpcg_fields <- function(warmup_iters, tol, maxiter, prefix = "") {
  # `assert_int` rejects anything outside integer range as non-integerish, so
  # the `arma::uword` overflow the C++ backend guards against cannot get past
  # here; no explicit upper bound is needed for it.
  checkmate::assert_int(
    warmup_iters,
    lower = 0L,
    .var.name = paste0(prefix, "warmup_iters")
  )
  checkmate::assert_number(
    tol,
    lower = 0,
    finite = TRUE,
    .var.name = paste0(prefix, "tol")
  )
  checkmate::assert_int(
    maxiter,
    lower = 0L,
    .var.name = paste0(prefix, "maxiter")
  )
  invisible(NULL)
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
  # validate an already created control object. Carrying the class is not
  # evidence the fields are still sound: a control object is an ordinary list,
  # so `ctrl$tol <- -1` or `ctrl$maxiter <- NULL` survives construction and is
  # only noticed further down - a negative `tol` not at all, a NULL field as
  # `missing value where TRUE/FALSE needed`, a negative `maxiter` as an
  # `arma::uword` overflow reported by the C++ backend as an int-range error.
  # Re-check the fields here, against the same rules the constructor uses.
  control_names <- names(formals(lobpcg_control))
  if (inherits(x, "slideimp_lobpcg_control")) {
    checkmate::assert_list(
      x,
      types = c("numeric", "integer"),
      any.missing = FALSE,
      names = "unique",
      .var.name = "lobpcg_control"
    )
    absent <- setdiff(control_names, names(x))
    if (length(absent) > 0L) {
      cli::cli_abort(c(
        "{.arg lobpcg_control} is missing {cli::qty(length(absent))}field{?s} {.field {absent}}.",
        "i" = "Rebuild it with {.fn lobpcg_control} rather than assigning into it."
      ))
    }
    check_lobpcg_fields(
      x$warmup_iters,
      x$tol,
      x$maxiter,
      prefix = "lobpcg_control$"
    )
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
    unknown <- setdiff(names(x), control_names)
    if (length(unknown) > 0L) {
      cli::cli_abort(c(
        "{cli::qty(length(unknown))}Unknown LOBPCG control option{?s}: {fmt_trunc(unknown, 10)}.",
        "i" = "Allowed options are: {.arg {control_names}}."
      ))
    }
    out <- do.call(lobpcg_control, x)
  }
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
#' expectation-maximization (EM) PCA imputation. Supports warm-start LOBPCG with
#' both the previous eigenblock and search direction.
#'
#' @inheritParams knn_imp
#'
#' @param ncp Integer. Number of principal components used to predict missing
#'   entries.
#' @param scale Logical. If `TRUE`, columns are scaled to unit variance.
#' @param method Character. PCA imputation method: either `"regularized"` or
#'   `"EM"`.
#' @param coeff.ridge Numeric. Ridge regularization, used only when
#'   `method = "regularized"`. Values `< 1` move toward EM PCA; values `> 1`
#'   move toward mean imputation.
#' @param row.w Row weights, normalized to sum to `1`. `NULL` (equal weights),
#'   a positive numeric vector of length `nrow(obj)`, or `"n_miss"`
#'   (down-weight rows with more missing values).
#' @param threshold Numeric. Convergence threshold.
#' @param seed Integer, numeric, or `NULL`. Random seed for reproducibility.
#'   Initialization `i` is seeded with `seed * (i - 1)`, exactly as
#'   `missMDA::imputePCA()` does. `seed = 0` therefore seeds every
#'   initialization with `0`, so all random restarts draw the same start and
#'   `nb.init` collapses to two distinct solves: the mean initialization and
#'   one random one. Because the product has to stay within integer range, the
#'   largest admissible value is `.Machine$integer.max %/% (nb.init - 1)`.
#' @param nb.init Integer. Number of random initializations. The first
#'   initialization is always mean imputation.
#' @param maxiter Integer. Maximum number of iterations.
#' @param miniter Integer. Minimum number of iterations. Must be less than or
#'   equal to `maxiter`.
#' @param solver Character. Eigensolver: `"auto"` (default), `"exact"`, or
#'   `"lobpcg"`. `"auto"` runs a short timed probe and picks `"lobpcg"` only
#'   when clearly faster. Consecutive EM calls warm-start LOBPCG with both the
#'   previous eigenblock and search direction. When `nb.init > 1`, the auto
#'   choice from the first init is reused. See Performance tips.
#' @param lobpcg_control A list of LOBPCG eigensolver control options, usually
#'   created by [lobpcg_control()]. A plain named list is also accepted. Its
#'   fields are validated on every call, so an object modified after
#'   construction is checked again rather than trusted. Ignored when
#'   `solver = "exact"`.
#' @param clamp Optional numeric vector `c(lower, upper)` bounding PCA-imputed
#'   values (use `-Inf`/`Inf` for one-sided, `NULL` for none). E.g., `c(0, 1)`
#'   for DNAm beta values. Observed values are not clamped.
#'
#' @returns A numeric matrix of the same dimensions as `obj`, with missing
#' values imputed. The returned object has class `slideimp_results`.
#'
#' @details
#' This function reimplements the PCA imputation method from the `missMDA`
#' package by Francois Husson and Julie Josse, based on Josse and Husson (2016).
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @section PCA Performance tips:
#' Speed comes from three levers: `solver` (through LOBPCG with warm-start),
#' `threshold`, and `scale`. Tune these first, then accuracy parameters
#' (`ncp`, `coeff.ridge`) on a representative subset.
#'
#' **Exact vs. LOBPCG with warm-start.** Whether `"lobpcg"` beats `"exact"`
#' depends on size and low-rankness: `"lobpcg"` is preferred for large, approximately
#' low-rank matrices with small `ncp`, and `"exact"` for small matrices
#' (including `slide_imp()` windows), where it is faster and more robust.
#' Separately, the warm-start makes each successive solve cheap: `pca_imp()`
#' warm-starts LOBPCG with the previous eigenblock and search direction, so once
#' imputed values stabilize, later solves converge in a few iterations. The
#' payoff therefore grows with the number of EM iterations, independent of
#' low-rankness. `solver = "auto"` (default) probes both and is a safe start.
#'
#' **Threshold.** The default `1e-6` is conservative; `1e-5` is often faster
#' with very similar values.
#'
#' **Scale.** For columns on a common scale (e.g., DNAm beta values in
#' `[0, 1]`), `scale = FALSE` can be faster and more accurate.
#'
#' **Parallel and BLAS.** In parallel via `tune_imp()` or `group_imp()` with a
#' multithreaded BLAS, set `pin_blas = TRUE` to avoid thread oversubscription.
#' On Windows, the stock BLAS can be slow. Advanced users can swap in
#' [OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows).
#'
#' See [Speeding up PCA imputation](https://hhp94.github.io/slideimp/articles/speeding-up-pca-imputation.html)
#' for the full workflow.
#'
#' @references
#' Josse J, Husson F (2013). Handling missing values in exploratory
#' multivariate data analysis methods. *Journal de la SFdS*, 153(2), 79-99.
#'
#' Josse J, Husson F (2016). missMDA: A Package for Handling Missing Values in
#' Multivariate Data Analysis. *Journal of Statistical Software*, 70(1), 1-31.
#' \doi{10.18637/jss.v070.i01}
#'
#' @examples
#' set.seed(123)
#' obj <- sim_mat(10, 10)$input
#' sum(is.na(obj))
#' obj[1:4, 1:4]
#'
#' # Randomly initialize missing values 5 times. The first initialization is
#' # mean imputation. Select `ncp` with `tune_imp()`.
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
  clamp = NULL,
  .progress = FALSE
) {
  # pre-conditioning
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
  checkmate::assert_int(nb.init, lower = 1, .var.name = "nb.init")
  # initialization `i` is seeded with `seed * (i - 1)` (missMDA parity, see the
  # `@param seed` note), so the largest value ever handed to set.seed() is
  # `seed * (nb.init - 1)` and THAT is what has to stay inside integer range.
  # Folding the bound into the assertion keeps the message in checkmate's voice
  # and names the argument the caller passed.
  checkmate::assert_int(
    seed,
    null.ok = TRUE,
    lower = 0,
    upper = .Machine$integer.max %/% max(nb.init - 1L, 1L),
    .var.name = "seed"
  )
  checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
  checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")
  # the C++ backend enforces this too, but only as a bare `Rcpp::stop` that
  # names neither argument nor value; validation belongs on the R side.
  if (miniter > maxiter) {
    cli::cli_abort(c(
      "{.arg miniter} must be less than or equal to {.arg maxiter}.",
      "x" = "{.arg miniter} is {miniter} and {.arg maxiter} is {maxiter}.",
      "i" = "The EM loop runs at least {.arg miniter} and at most {.arg maxiter} iterations."
    ))
  }
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
  checkmate::assert_flag(.progress, .var.name = ".progress")
  trace_iter <- if (isTRUE(.progress)) 10L else 0L

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
    switch(solver, exact = 0L, lobpcg = 1L, auto = 2L)
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
  best_diag <- NULL

  # solver_chosen codes:
  #   0 = forced exact
  #   1 = forced lobpcg/hybrid
  #   2 = auto had no reason/opportunity to choose. Treat as exact
  #   3 = auto chose exact
  #   4 = auto chose lobpcg
  #
  # only 3 and 4 are verdicts from a completed probe. `solver_chosen` below
  # says which solver RAN, and reports "exact" for four different reasons:
  # forced, demoted by `auto_force_exact` before the kernel was ever called,
  # code 2, and code 3. A caller that wants to reuse the decision - as
  # `slide_imp()` does across windows - has to be able to tell a verdict from
  # the other three, so `solver_probed` is reported alongside it.
  solver_probed <- FALSE
  locked_solver <- NULL
  resolved_solver_code <- if (isTRUE(auto_force_exact)) {
    0L
  } else {
    switch(solver, exact = 0L, lobpcg = 1L, auto = NA_integer_)
  }

  for (i in seq_len(nb.init)) {
    if (!is.null(seed)) {
      set.seed(seed * (i - 1L)) # exactly as missMDA does
    }

    if (trace_iter > 0L && nb.init > 1L) {
      cli::cli_inform("Initialization {i}/{nb.init}")
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
      lobpcg_maxiter = lobpcg_control$maxiter,
      trace_iter = trace_iter
    )

    # after the auto probe, select a solver for the rest of the inits.
    if (solver == "auto" && i == 1L) {
      chosen <- as.integer(res.impute$solver_chosen)
      if (length(chosen) != 1L || is.na(chosen)) {
        chosen <- 2L
      }
      locked_solver <- if (chosen %in% c(1L, 4L)) 1L else 0L
      resolved_solver_code <- locked_solver
      solver_probed <- chosen %in% c(3L, 4L)
    }

    cur_obj <- res.impute$mse
    if (cur_obj < init_obj) {
      best_imputed <- res.impute$imputed_values
      init_obj <- cur_obj
      best_diag <- list(
        n_iter = res.impute$n_iter,
        criterion_final = res.impute$criterion_final,
        converged = res.impute$converged,
        n_exact = res.impute$n_exact,
        n_lobpcg_ok = res.impute$n_lobpcg_ok,
        n_lobpcg_bad = res.impute$n_lobpcg_bad
      )
    }
  }

  # column/row indices returned by backend are already original-matrix indices.
  if (is.null(best_imputed) || nrow(best_imputed) == 0L) {
    cli::cli_abort("Internal error: PCA imputation produced no imputed values.")
  }

  # the non-convergence warning is raised here, not in the kernel. Warning from
  # C++ goes through Rf_warning, which longjmps out of the EM frame whenever the
  # caller has set options(warn = 2) or established a handler that invokes a
  # restart, and a longjmp skips every C++ destructor in that frame. Raising it
  # from R also means one warning per call rather than one per `nb.init`
  # restart, and it reports the restart whose values were actually kept.
  if (!isTRUE(best_diag$converged)) {
    cli::cli_warn(c(
      "Stopped after {maxiter} iteration{?s} without converging.",
      "i" = "Final criterion {signif(best_diag$criterion_final, 3)} is above {.arg threshold} ({threshold}).",
      "i" = "Increase {.arg maxiter} or relax {.arg threshold}."
    ))
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

  # after the write-back, NA can remain only in the ineligible columns that had
  # any, plus wherever the kernel handed back NA for a cell it was asked to
  # fill. Mean-impute exactly those columns, and skip the call when there are
  # none. An unconditional `mean_imp_col(obj)` copied the matrix twice more -
  # once into the C++ result, once when Rcpp wrapped it - and scanned every
  # column for Inf and for its mean, all to change nothing on an all-eligible
  # input. Same shape as the post step in knn_imp().
  if (post_imp) {
    na_cols <- sort(unique(c(
      which(!eligible & cmiss > 0L),
      as.integer(best_imputed[is.na(best_imputed[, 3]), 2])
    )))
    if (length(na_cols) > 0L) {
      obj <- mean_imp_col(obj, subset = na_cols)
    }
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
  attr(out, "solver_probed") <- solver_probed
  attr(out, "n_iter") <- best_diag$n_iter
  attr(out, "criterion_final") <- best_diag$criterion_final
  attr(out, "converged") <- best_diag$converged
  attr(out, "n_exact") <- best_diag$n_exact
  attr(out, "n_lobpcg_ok") <- best_diag$n_lobpcg_ok
  attr(out, "n_lobpcg_bad") <- best_diag$n_lobpcg_bad
  return(out)
}
