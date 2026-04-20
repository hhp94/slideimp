#' Impute Numeric Matrix with PCA Imputation
#'
#' Impute missing values in a numeric matrix using (regularized) iterative PCA.
#'
#' @details
#' This algorithm is based on the original `missMDA::imputePCA` function and is
#' optimized for tall or wide numeric matrices.
#'
#' @inheritParams knn_imp
#' @param ncp Integer. Number of components used to predict the missing entries.
#' @param scale Logical. If `TRUE` (default), variables are scaled to have
#'   unit variance.
#' @param method Character. Either `"regularized"` (default) or `"EM"`.
#' @param coeff.ridge Numeric. Ridge regularization coefficient (default is 1).
#'   Only used if `method = "regularized"`. Values < 1 regularize less (closer
#'   to EM); values > 1 regularize more (closer to mean imputation).
#' @param row.w Row weights (internally normalized to sum to 1). Can be one of:
#'   * `NULL` (default): All rows weighted equally.
#'   * A numeric vector: Custom positive weights of length `nrow(obj)`.
#'   * `"n_miss"`: Rows with more missing values receive lower weight.
#' @param threshold Numeric. The threshold for assessing convergence.
#' @param seed Integer. Random number generator seed.
#' @param nb.init Integer. Number of random initializations. The first
#'   initialization is always mean imputation.
#' @param maxiter Integer. Maximum number of iterations for the algorithm.
#' @param miniter Integer. Minimum number of iterations for the algorithm.
#'
#' @inherit knn_imp return
#'
#' @references
#' Josse, J. & Husson, F. (2013). Handling missing values in exploratory
#' multivariate data analysis methods. Journal de la SFdS. 153 (2), pp. 79-99.
#'
#' Josse, J. and Husson, F. (2016). missMDA: A Package for Handling Missing
#' Values in Multivariate Data Analysis. Journal of Statistical Software,
#' 70 (1), pp 1-31. \doi{10.18637/jss.v070.i01}
#'
#' @author Francois Husson and Julie Josse (original `missMDA` implementation).
#'
#' @examples
#' obj <- sim_mat(10, 10)$input
#' sum(is.na(obj))
#' obj[1:4, 1:4]
#' # Randomly initialize missing values 5 times (1st time is mean).
#' pca_imp(obj, ncp = 2, nb.init = 5)
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
  colmax = 0.9,
  post_imp = TRUE,
  na_check = TRUE
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
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  checkmate::assert_flag(na_check, .var.name = "na_check")

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
  eligible <- miss_rate < min(colmax, 1) &
    !(obj_vars < .Machine$double.eps | is.na(obj_vars))

  n_elig <- sum(eligible)
  cap <- min(nrow(obj) - 2L, n_elig - 1L)
  row_bound <- nrow(obj) - 2L <= n_elig - 1L
  if (ncp > cap) {
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

  # iterate over initializations; C++ only returns the imputed triples + mse,
  # not the reconstructed matrix.
  init_obj <- Inf
  best_imputed <- NULL
  for (i in seq_len(nb.init)) {
    if (!is.null(seed)) {
      set.seed(seed * (i - 1L)) # exactly as missMDA does
    }
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
      coeff_ridge = coeff.ridge
    )
    cur_obj <- res.impute$mse
    if (cur_obj < init_obj) {
      # res.impute$imputed_values: N_miss x 3 matrix of
      # (row_idx, col_idx, value) in original-matrix, 1-based coordinates
      best_imputed <- res.impute$imputed_values
      init_obj <- cur_obj
    }
  }

  # column/row indices returned by C++ are already original-matrix indices,
  imp_indices <- cbind(best_imputed[, 1], best_imputed[, 2])
  obj[imp_indices] <- best_imputed[, 3]

  if (post_imp) {
    obj <- mean_imp_col(obj)
  }

  return(
    as_slideimp_results(
      obj,
      "pca",
      fallback = FALSE,
      post_imp = post_imp,
      na_check = na_check
    )
  )
}
