#' Validated Clamp Bounds
#'
#' Checks `clamp` and returns it as `c(lower, upper)`, or `NULL` for no
#' clamping.
#'
#' @param clamp Numeric or `NULL`. The `clamp` argument as the user passed it.
#' @param arg Character. The argument name the error messages give. Defaults
#'   to `"clamp"`.
#'
#' @returns `NULL` or an unnamed numeric vector of length 2.
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
      "i" = "Use {.code c(-Inf, upper)} or {.code c(lower, Inf)} for one-sided clamping."
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
      "i" = "Use {.code -Inf} or {.code Inf} for an open bound, as in {.code c(-Inf, 1)}."
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
      "{.arg {arg}} must be {.code c(lower, upper)} with lower <= upper."
    )
  }

  return(unname(as.numeric(clamp)))
}

#' LOBPCG Eigensolver Control Options
#'
#' Construct a validated list of control options for the LOBPCG eigensolver
#' used by [pca_imp()]. Most users do not need to call it.
#'
#' @param warmup_iters Integer. Number of initial [pca_imp()] iterations that
#'   use the exact eigensolver before LOBPCG starts. Must be non-negative.
#' @param tol Numeric. Convergence tolerance for the LOBPCG eigensolver. Must
#'   be non-negative and finite.
#' @param maxiter Integer. Maximum number of LOBPCG iterations per eigenvalue
#'   solve. Must be non-negative, and positive when [pca_imp()] runs with
#'   `solver = "auto"` or `"lobpcg"`. Use `solver = "exact"`, not
#'   `maxiter = 0`, to select the exact solver.
#'
#' @details
#' Under `solver = "auto"`, [pca_imp()] raises a `warmup_iters` that was not
#' set for a large `ncp`, up to `50`. Set `warmup_iters` to keep a fixed value.
#'
#' @returns A `slideimp_lobpcg_control` object: a named list with
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
#' pca_imp(obj, ncp = 2, lobpcg_control = lobpcg_control(tol = 1e-8))
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

#' LOBPCG Control Field Rules
#'
#' Holds the field rules for a LOBPCG control object, so that
#' [lobpcg_control()] and `new_lobpcg_control()` cannot drift apart. The latter
#' re-checks an object that already carries the class, because a control
#' object is an ordinary list that a caller can assign into after construction.
#'
#' @param warmup_iters,tol,maxiter Numeric. The three control fields, as
#'   [lobpcg_control()] takes them.
#' @param prefix Character. The prefix for each `.var.name`, so the message
#'   names the argument the user passed. Defaults to `""`.
#'
#' @returns `NULL`, invisibly. Called for the side effect of aborting.
#'
#' @keywords internal
#' @noRd
check_lobpcg_fields <- function(warmup_iters, tol, maxiter, prefix = "") {
  # `assert_int` rejects anything outside integer range as non-integerish, so
  # a value too large for the `int` the C++ backend takes cannot get past
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

#' Validated LOBPCG Control Object
#'
#' Turns the `lobpcg_control` argument of [pca_imp()] into a validated
#' `slideimp_lobpcg_control` object. It takes `NULL`, a
#' `slideimp_lobpcg_control` object, or a partial named list. With
#' `solver = "exact"` it returns `lobpcg_control(maxiter = 0L)` and ignores
#' `x`.
#'
#' @param x A `slideimp_lobpcg_control` object, a named list, or `NULL`. The
#'   `lobpcg_control` argument as the user passed it.
#' @param solver Character. One of `"auto"`, `"exact"`, or `"lobpcg"`. The
#'   solver that [pca_imp()] resolved.
#'
#' @returns A `slideimp_lobpcg_control` object.
#'
#' @keywords internal
#' @noRd
new_lobpcg_control <- function(x, solver = c("auto", "exact", "lobpcg")) {
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
  # validate an already created control object. It is an ordinary list, so
  # its fields can be changed after construction (`ctrl$tol <- -1`); re-check
  # them against the constructor's rules.
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
        "i" = "Rebuild it with {.fn lobpcg_control}."
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
        "{.arg lobpcg_control} must be {.code NULL}, a named list, or the output of {.fn lobpcg_control}."
      )
    }
    if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
      cli::cli_abort("{.arg lobpcg_control} must be a named list.")
    }
    unknown <- setdiff(names(x), control_names)
    if (length(unknown) > 0L) {
      cli::cli_abort(c(
        "{length(unknown)} unknown option{?s} in {.arg lobpcg_control}:",
        "x" = "{fmt_trunc(unknown, 10)}",
        "i" = "Allowed: {.arg {control_names}}."
      ))
    }
    out <- do.call(lobpcg_control, x)
  }
  explicit <- attr(out, "explicit", exact = TRUE)

  if (is.null(explicit)) {
    # a control object without the attribute: treat all fields as explicit,
    # so none is overridden.
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

# What pca_imp() reports about a call: the attributes of its result, in the
# order it sets them. Shared by the no-missing-values exit and the main path so that a caller cannot tell from
# the result which one ran: the same fields with the same types. The defaults
# describe a call that ran no solve. No solver was chosen or probed, no
# iteration ran, and `converged` is NA because there was nothing to converge.
# `n_iter` is an integer on both paths: the kernel returns its count as one.
pca_imp_info <- function(
  post_imp,
  solver,
  solver_chosen = NA_character_,
  solver_probed = FALSE,
  diag = list(
    n_iter = 0L,
    criterion_final = NA_real_,
    converged = NA,
    n_exact = 0L,
    n_lobpcg_ok = 0L,
    n_lobpcg_bad = 0L
  ),
  mean_filled_cols = integer(0),
  mean_filled_cells = mean_filled_cells_matrix()
) {
  list(
    imp_method = "pca",
    fallback = FALSE,
    post_imp = post_imp,
    mean_filled_cols = mean_filled_cols,
    mean_filled_cells = mean_filled_cells,
    solver_requested = solver,
    solver_chosen = solver_chosen,
    solver_probed = solver_probed,
    n_iter = diag$n_iter,
    criterion_final = diag$criterion_final,
    converged = diag$converged,
    n_exact = diag$n_exact,
    n_lobpcg_ok = diag$n_lobpcg_ok,
    n_lobpcg_bad = diag$n_lobpcg_bad
  )
}

# The result of pca_imp(): the shared constructor plus the attributes
# only pca_imp() carries, from pca_imp_info().
pca_imp_results <- function(obj, info) {
  out <- new_slideimp_results(
    obj,
    info$imp_method,
    fallback = info$fallback,
    post_imp = info$post_imp,
    mean_filled_cols = info$mean_filled_cols,
    mean_filled_cells = info$mean_filled_cells
  )
  shared <- c(
    "imp_method",
    "fallback",
    "post_imp",
    "mean_filled_cols",
    "mean_filled_cells"
  )
  for (nm in setdiff(names(info), shared)) {
    attr(out, nm) <- info[[nm]]
  }
  out
}

#' PCA Imputation for Numeric Matrices
#'
#' Impute missing values in a numeric matrix with regularized or
#' expectation-maximization (EM) PCA imputation. The eigensolver can be LOBPCG,
#' warm-started from the previous eigenblock and search direction.
#'
#' @inheritParams knn_imp
#'
#' @param ncp Integer. Number of principal components used to impute missing
#'   values.
#' @param scale Logical. If `TRUE`, columns are scaled to unit variance.
#' @param method Character. One of `"regularized"` or `"EM"`. The PCA
#'   imputation method.
#' @param coeff.ridge Numeric. The ridge regularization for
#'   `method = "regularized"`, ignored otherwise. Values below `1` move toward
#'   EM PCA, and values above `1` move toward mean imputation.
#' @param row.w Numeric vector, `"n_miss"`, or `NULL`. The row weights,
#'   normalized to sum to `1`: positive values of length `nrow(obj)`,
#'   `"n_miss"` to down-weight rows with more missing values, or `NULL` for
#'   equal weights.
#' @param threshold Numeric. Convergence threshold.
#' @param seed Integer, numeric, or `NULL`. Random seed for reproducibility.
#'   As in `missMDA::imputePCA()`, initialization `i` is seeded with
#'   `seed * (i - 1)`. So `seed = 0` gives every random initialization the same
#'   start, and `nb.init` then yields two distinct solves: the mean
#'   initialization and one random one. With `nb.init > 1`, `seed` can be at
#'   most `.Machine$integer.max %/% (nb.init - 1)`.
#' @param nb.init Integer. Number of initializations. The first is always mean
#'   imputation, and the others are random.
#' @param maxiter Integer. Maximum number of iterations.
#' @param miniter Integer. Minimum number of iterations, at most `maxiter`.
#' @param solver Character. One of `"auto"`, `"exact"`, or `"lobpcg"`. The
#'   eigensolver. `"auto"` uses `"exact"` for small problems and otherwise runs
#'   a short timed probe, picking `"lobpcg"` only when it is clearly faster.
#'   With `nb.init > 1`, the choice from the first initialization is reused.
#'   See the PCA Performance tips section.
#' @param lobpcg_control A list or `NULL`. LOBPCG control options, from
#'   [lobpcg_control()] or as a plain named list, or `NULL` for the defaults.
#'   Ignored when `solver = "exact"`.
#' @param clamp Numeric vector of length 2 or `NULL`. The bounds
#'   `c(lower, upper)` for PCA-imputed values, such as `c(0, 1)` for DNA
#'   methylation beta values, or `NULL` for none. Use `-Inf` or `Inf` for a
#'   one-sided bound. Observed values are not clamped.
#'
#' @inherit knn_imp return
#'
#' @details
#' `pca_imp()` reimplements the PCA imputation method from the `missMDA`
#' package by Francois Husson and Julie Josse, based on Josse and Husson (2016).
#'
#' A column whose observed values do not vary takes no part in PCA imputation.
#' Column mean fill treats it like a column over `colmax`, in `post_imp` and in
#' the `mean_filled_cols` attribute.
#'
#' With `subset`, the result holds only the selected columns, in the order
#' they have in `obj`. The other columns take part in the imputation as they
#' do without `subset`. So `subset` makes the result smaller, not the imputation
#' faster. A name that is not in `colnames(obj)`, or a `subset` that
#' selects no column, is an error.
#'
#' @inheritSection slideimp-package Column mean fill
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @section PCA Performance tips:
#' Three arguments set the speed: `solver`, `threshold`, and `scale`. Tune
#' these first, then the accuracy arguments `ncp` and `coeff.ridge`, on a
#' representative subset.
#'
#' **Exact vs. LOBPCG.** `"lobpcg"` suits large, approximately low-rank
#' matrices with small `ncp`. `"exact"` suits small matrices, such as
#' `slide_imp()` windows, where it is faster. It is also more reliable: it
#' always converges, to eigenvectors exact to machine precision, while LOBPCG
#' stops at the tolerance set in [lobpcg_control()] and can take many
#' iterations, or reach `maxiter`, when the eigenvalues around the `ncp`-th one
#' are close together. The warm start makes each later LOBPCG solve cheap:
#' `pca_imp()` starts LOBPCG from the previous eigenblock and search direction,
#' so once the imputed values stabilize, each solve takes a few iterations. The
#' gain grows with the number of EM iterations, whatever the rank.
#' `solver = "auto"`, the default, chooses between the two and is a safe start.
#'
#' **Threshold.** The default `1e-6` is conservative. `1e-5` is often faster
#' and gives very similar values.
#'
#' **Scale.** For columns on a common scale, such as DNA methylation beta
#' values in `[0, 1]`, `scale = FALSE` can be faster and more accurate.
#'
#' **BLAS.** On Windows, the BLAS that ships with R can be slow, and
#' [OpenBLAS](https://github.com/david-cortes/R-openblas-in-windows) can
#' replace it.
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
#' @seealso
#' * [group_imp()] to impute within feature groups.
#' * [mipca()] for multiple imputation with PCA.
#'
#' @examples
#' set.seed(123)
#' obj <- sim_mat(10, 10)$input
#' sum(is.na(obj))
#' obj[1:4, 1:4]
#'
#' # Five initializations: the first is mean imputation, the other four are
#' # random. Select `ncp` with `tune_imp()`.
#' pca_imp(obj, ncp = 2, nb.init = 5, seed = 123)
#'
#' # Impute every column, but return only the first three
#' pca_imp(obj, ncp = 2, subset = 1:3)
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
  clamp = NULL,
  subset = NULL,
  .progress = FALSE
) {
  checkmate::assert_matrix(
    obj,
    mode = "double",
    null.ok = FALSE,
    .var.name = "obj"
  )
  # the row counts only for `row.w = "n_miss"`, which pca_imp_parts()
  # validates with the rest of `row.w`
  scan <- col_scan(obj, var = TRUE, rows = identical(row.w, "n_miss"))
  refuse_nonfinite(scan)
  res <- pca_imp_parts(
    obj = obj,
    scan = scan,
    ncp = ncp,
    scale = scale,
    method = method,
    coeff.ridge = coeff.ridge,
    row.w = row.w,
    threshold = threshold,
    seed = seed,
    nb.init = nb.init,
    maxiter = maxiter,
    miniter = miniter,
    solver = solver,
    lobpcg_control = lobpcg_control,
    colmax = colmax,
    post_imp = post_imp,
    clamp = clamp,
    subset = subset,
    .progress = .progress
  )
  pca_imp_results(fill_imp_parts(res), res$info)
}

# pca_imp() without its door, as knn_imp_parts() is knn_imp() without its
# door: `obj` and its `scan`, with the variances, then the formals of
# pca_imp(), and the same list returned. The scan carries the row counts
# when the door took them; a slice's scan does not, and they are counted
# here. The PCA part covers every eligible column, kept or not.
pca_imp_parts <- function(
  obj,
  scan,
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
  clamp = NULL,
  subset = NULL,
  .progress = FALSE
) {
  # pre-conditioning
  if (is.null(scan$var)) {
    stop(
      "internal error: pca_imp_parts() was given a scan without variances",
      call. = FALSE
    )
  }
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
  # `seed * (nb.init - 1)`, which has to stay inside integer range.
  checkmate::assert_int(
    seed,
    null.ok = TRUE,
    lower = 0,
    upper = .Machine$integer.max %/% max(nb.init - 1L, 1L),
    .var.name = "seed"
  )
  checkmate::assert_int(maxiter, lower = 1, .var.name = "maxiter")
  checkmate::assert_int(miniter, lower = 1, .var.name = "miniter")
  # the C++ backend enforces this too, but its message names neither argument
  # nor value.
  if (miniter > maxiter) {
    cli::cli_abort(c(
      "{.arg miniter} must be less than or equal to {.arg maxiter}.",
      "x" = "{.arg miniter} is {miniter} and {.arg maxiter} is {maxiter}."
    ))
  }
  # solver resolves
  solver <- match.arg(solver)
  lobpcg_control <- new_lobpcg_control(lobpcg_control, solver = solver)

  lobpcg_explicit <- attr(lobpcg_control, "explicit", exact = TRUE)
  warmup_explicit <- isTRUE(lobpcg_explicit[["warmup_iters"]])

  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  clamp <- resolve_clamp(clamp, arg = "clamp")
  # the columns the result holds, ascending: NULL for every column. Every
  # eligible column still takes part in the fit, since every hole in every
  # eligible column is updated at each EM iteration and feeds the next one,
  # so `subset` narrows what is returned, not what is computed. Its refusals
  # name the caller, pca_imp() for a user.
  keep <- resolve_narrow_subset(subset, obj, call = parent.frame())
  out_cols <- if (is.null(keep)) seq_len(ncol(obj)) else keep
  checkmate::assert_flag(.progress, .var.name = ".progress")
  trace_iter <- if (isTRUE(.progress)) 10L else 0L

  # per-column missingness, from the scan
  cmiss <- scan$n_miss

  # early exit if no column the result holds has a hole. Those columns come
  # back unchanged, but in the same shape as the main path's result, as
  # @returns promises: fill_imp_parts() takes them from `obj`.
  if (!any(cmiss[out_cols] > 0L)) {
    if (is.null(keep)) {
      cli::cli_inform("No missing values in input. Returning input unchanged.")
    } else {
      cli::cli_inform(
        "No missing values in the {.arg subset} columns. Returning them unchanged."
      )
    }
    return(list(
      obj = obj,
      parts = NULL,
      keep = keep,
      info = pca_imp_info(post_imp = post_imp, solver = solver),
      n_miss = cmiss
    ))
  }

  miss_rate <- cmiss / nrow(obj)

  # eligibility: below colmax and non-degenerate variance
  obj_vars <- scan$var
  eligible <- miss_rate <= min(colmax, 1) &
    !(obj_vars < .Machine$double.eps | is.na(obj_vars))

  n_elig <- sum(eligible)

  max_ncp_rows <- nrow(obj) - 2L
  max_ncp_cols <- n_elig - 1L
  cap <- min(max_ncp_rows, max_ncp_cols)

  if (cap < 1L) {
    cli::cli_abort(
      c(
        "PCA imputation is infeasible: {.arg obj} needs at least 3 rows and 2 eligible columns.",
        "x" = "{.arg obj} has {nrow(obj)} row{?s} and {n_elig} eligible column{?s}.",
        "i" = "Relax {.arg colmax}, or add rows or columns."
      ),
      class = "slideimp_infeasible"
    )
  }

  if (ncp > cap) {
    row_bound <- max_ncp_rows <= max_ncp_cols

    cli::cli_abort(
      c(
        "{.arg ncp} ({ncp}) must be at most the number of usable components ({cap}).",
        "i" = if (row_bound) {
          "The limit is the number of rows ({nrow(obj)}) minus 2. Reduce {.arg ncp} or add rows."
        } else {
          "The limit is the number of eligible columns ({n_elig}) minus 1. Reduce {.arg ncp} or relax {.arg colmax}."
        }
      ),
      class = "slideimp_infeasible"
    )
  }

  # at least one eligible column the result holds must have missingness
  if (!any(cmiss[out_cols][eligible[out_cols]] > 0L)) {
    cli::cli_abort(
      c(
        if (is.null(keep)) {
          "All columns with missing values are ineligible for PCA imputation."
        } else {
          "All {.arg subset} columns with missing values are ineligible for PCA imputation."
        },
        "x" = "Each has a proportion of missing values above {.arg colmax} ({colmax}) or zero variance.",
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
    # per-row missingness over every column of `obj`, eligible or not, from
    # the scan, or counted here for a slice, whose scan cannot carry
    # them. `subset` does not narrow it, since it narrows only the result.
    n_miss_per_row <- scan$row_miss
    if (is.null(n_miss_per_row)) {
      n_miss_per_row <- col_scan(obj, rows = TRUE)$row_miss
    }
    row.w <- 1 - (n_miss_per_row / ncol(obj))
    row.w[row.w < 1e-10] <- 1e-10
  }

  init_obj <- Inf
  best_loc <- NULL
  best_value <- NULL
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
      eligible_idx = eligible_idx - 1L,
      ncp = as.integer(ncp),
      scale = scale,
      regularized = (method == "regularized"),
      threshold = threshold,
      init = if (i == 1L) 0L else i,
      maxiter = as.integer(maxiter),
      miniter = as.integer(miniter),
      row_w = row.w,
      coeff_ridge = coeff.ridge,
      solver = solver_code_i,
      warmup_iters = lobpcg_control$warmup_iters,
      lobpcg_tol = lobpcg_control$tol,
      lobpcg_maxiter = lobpcg_control$maxiter,
      trace_iter = trace_iter
    )
    # the kernel sets no dimnames, so `loc` gets a part map's column names here.
    dimnames(res.impute$loc) <- list(NULL, c("row", "col"))

    # after the auto probe, select a solver for the rest of the inits.
    if (solver == "auto" && i == 1L) {
      chosen <- res.impute$solver_chosen
      locked_solver <- if (chosen %in% c(1L, 4L)) 1L else 0L
      resolved_solver_code <- locked_solver
      solver_probed <- chosen %in% c(3L, 4L)
    }

    # the kernel's mse sums the squared residual over every cell of each
    # eligible column, holes included, before taking the holes' share back
    # out. So a NaN or Inf at any imputed cell makes it non-finite: this is
    # where a broken fit shows, not in the values.
    cur_obj <- res.impute$mse
    if (!is.finite(cur_obj)) {
      stop(
        "slideimp internal, pca_imp: the fit's mse is not finite (",
        cur_obj,
        ") at initialization ",
        i,
        call. = FALSE
      )
    }
    if (cur_obj < init_obj) {
      best_loc <- res.impute$loc
      best_value <- res.impute$value
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
  if (is.null(best_value) || length(best_value) == 0L) {
    cli::cli_abort(
      "slideimp internal, pca_imp: the fit produced no imputed values"
    )
  }

  # the non-convergence warning is raised here, not in the kernel: a warning
  # from C++ can longjmp out of the EM frame (options(warn = 2), a restarting
  # handler) and skip its destructors. It is also one warning per call, for
  # the restart whose values were kept.
  if (!isTRUE(best_diag$converged)) {
    cli::cli_warn(c(
      "Stopped after {maxiter} iteration{?s} without converging.",
      "i" = "Final criterion {signif(best_diag$criterion_final, 3)} is above {.arg threshold} ({threshold}).",
      "i" = "Increase {.arg maxiter} or relax {.arg threshold}."
    ))
  }

  if (!is.null(clamp)) {
    best_value <- pmin(pmax(best_value, clamp[1L]), clamp[2L])
  }

  # three parts, which together cover every column the result holds. The PCA
  # part is the kernel's `loc` and `value`: the eligible columns and every
  # hole in them, `loc` an integer matrix with columns `row` and `col` in the
  # input's 1-based numbering, in the column-major order the part's validator
  # checks. With `subset`, the door cuts it to the kept
  # columns. The other two are the ineligible columns the result holds: those
  # with a hole, under `post_imp` each hole set to its column mean and
  # otherwise left NA, and those without, a part with nothing to fill. No
  # ineligible column is all NA, since the door refuses one
  # (refuse_nonfinite()), so each has a mean. No part holds NA, as the mse check above guarantees for the
  # kernel's values.
  # `eligible` comes from the scan, which carries no names, so these are
  # plain indices, as parts take.
  pca_part <- new_part(
    cols = eligible_idx,
    map = best_loc,
    values = best_value,
    kind = "pca"
  )
  rm(best_loc, best_value)

  # the ineligible columns the result holds. The PCA part keeps every
  # eligible column, kept or not, and the door cuts it to `keep`.
  mean_cols <- which(!eligible)
  if (!is.null(keep)) {
    mean_cols <- intersect(mean_cols, keep)
  }
  holed_cols <- mean_cols[cmiss[mean_cols] > 0L]
  # The ineligible columns with a hole are one part, filled under `post_imp`
  # and `left_na` without it, so a `left_na` part names only columns that
  # keep a hole (check_parts()). Those with none are a part with nothing to
  # fill, so fallback_part() copies only the columns it fills.
  mean_parts <- list(
    if (post_imp) {
      fallback_part(obj, holed_cols, "mean")
    } else {
      new_part(cols = holed_cols, kind = "mean", left_na = TRUE)
    },
    new_part(cols = setdiff(mean_cols, holed_cols), kind = "mean")
  )
  parts <- c(list(pca_part), mean_parts)
  rm(pca_part, mean_parts)

  info <- pca_imp_info(
    post_imp = post_imp,
    solver = solver,
    solver_chosen = if (resolved_solver_code == 1L) "lobpcg" else "exact",
    solver_probed = solver_probed,
    diag = best_diag,
    # what `post_imp` set to a column mean: the ineligible columns that had a
    # hole, numbered as the result's columns. `mean_filled_cells` keeps its
    # empty default, since the PCA part fills every hole of its columns.
    mean_filled_cols = if (!post_imp) {
      integer(0)
    } else if (is.null(keep)) {
      holed_cols
    } else {
      match(holed_cols, keep)
    }
  )

  list(obj = obj, parts = parts, keep = keep, info = info, n_miss = cmiss)
}
