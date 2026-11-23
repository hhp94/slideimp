#' K-Nearest Neighbor Imputation for Numeric Matrices
#'
#' Impute missing values in a numeric matrix with K-nearest neighbors (K-NN),
#' searching the whole matrix for neighbors.
#'
#' @param obj A double matrix. Samples are in rows and features in columns.
#' @param k Integer. Number of nearest neighbors for K-NN imputation.
#' @param colmax Numeric scalar between `0` and `1`. Columns with a proportion
#'  of missing values greater than `colmax` take no part in the imputation, as
#'  described under Column mean fill.
#' @param method Character. One of `"euclidean"` or `"manhattan"`. The distance
#'  used to find neighbors, as defined under Details.
#' @param cores Integer. Number of cores for K-NN imputation. Defaults to `1`.
#' @param post_imp Logical. If `TRUE`, column means replace the missing values
#'  that remain after imputation, as described under Column mean fill.
#' @param subset Optional character or integer vector. The columns to impute
#'  and return, by name or by number, or `NULL` for every column. Defaults to
#'  `NULL`.
#' @param dist_pow Numeric. Power used to penalize more distant neighbors in
#'  the weighted average, where `0` gives an unweighted average. The weighting
#'  rule is given under Details in [knn_imp()].
#' @param .progress Logical. If `TRUE`, show imputation progress.
#'
#' @returns A `slideimp_results` object: a numeric matrix with the rows of
#' `obj` and the columns that `subset` selects, with missing values imputed.
#'
#' @details
#' `knn_imp()` imputes each column from the `k` columns nearest to it, found by
#' an exact brute-force search. The distance between the target column and a
#' candidate column is a mean over the rows where both are observed. Candidates
#' that share different numbers of rows with the target therefore stay
#' comparable. `method = "euclidean"` averages squared differences and takes no
#' square root. `method = "manhattan"` averages absolute differences.
#'
#' A candidate that shares no observed row with the target is never a
#' neighbor. A column can therefore have fewer than `k` neighbors, or none, and
#' K-NN does not impute a column with none.
#'
#' A column with few observed values can be selected as a neighbor because it
#' agrees with the target on one or two shared rows. Such a neighbor has no
#' value at most of the rows to impute. A lower `colmax` keeps these columns
#' out of the neighbors.
#'
#' When `dist_pow > 0`, imputed values are distance-weighted averages of the
#' selected neighbors: a neighbor at distance `d` is weighted
#' `(d_min / d)^dist_pow`, where `d_min` is the smallest positive distance
#' among the neighbors chosen for that column. `d_min` cancels out of the
#' average and only keeps the largest weight at `1`.
#'
#' Under `"euclidean"`, `d` is a mean squared difference, so weighting by
#' `1 / d^dist_pow` equals weighting by `1 / d^(2 * dist_pow)` on the usual,
#' square-rooted euclidean distance. A `dist_pow` tuned for one `method`
#' therefore does not carry over to the other.
#'
#' Neighbors are selected once per column. A neighbor that is missing at a row
#' drops out of the average for that row, which is taken over the neighbors
#' observed there. K-NN does not impute a row where none of the neighbors is
#' observed.
#'
#' With `dist_pow > 0`, a neighbor at distance zero agrees with the target on
#' every shared row and would take an infinite weight. It takes weight `1`
#' instead, and the other neighbors take a small positive weight. Exact
#' matches therefore dominate the average, but the other neighbors stay in it.
#' When every neighbor is at distance zero, all weights are `1`.
#'
#' The result holds the `subset` columns in the order they have in `obj`.
#' Every column within `colmax` can still be a neighbor, in `subset` or not. A
#' name that is not in `colnames(obj)`, or a `subset` that selects no column,
#' is an error.
#'
#' @inheritSection slideimp-package Column mean fill
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @section K-NN performance optimization:
#' - Use `subset` when only specific columns need imputation.
#' - Use [group_imp()] or [slide_imp()] for very large matrices.
#'
#' @references
#' Troyanskaya O, Cantor M, Sherlock G, Brown P, Hastie T, Tibshirani R,
#' Botstein D, Altman RB (2001). Missing value estimation methods for DNA
#' microarrays. *Bioinformatics*, 17(6), 520-525.
#' \doi{10.1093/bioinformatics/17.6.520}
#'
#' @seealso [group_imp()] to impute within feature groups.
#'
#' @examples
#' set.seed(123)
#' obj <- sim_mat(20, 20, perc_col_na = 1)$input
#' sum(is.na(obj))
#'
#' # Select `k` with `tune_imp()`.
#' result <- knn_imp(obj, k = 10)
#' result
#'
#' # Impute and return only the first three columns
#' knn_imp(obj, k = 10, subset = 1:3)
#'
#' @export
knn_imp <- function(
  obj,
  k,
  colmax = 0.9,
  method = c("euclidean", "manhattan"),
  cores = 1,
  post_imp = TRUE,
  subset = NULL,
  dist_pow = 0,
  .progress = FALSE
) {
  checkmate::assert_matrix(
    obj,
    mode = "double",
    min.rows = 1,
    min.cols = 2,
    null.ok = FALSE,
    .var.name = "obj"
  )
  scan <- col_scan(obj)
  refuse_nonfinite(scan)
  res <- knn_imp_parts(
    obj = obj,
    scan = scan,
    k = k,
    colmax = colmax,
    method = method,
    cores = cores,
    post_imp = post_imp,
    subset = subset,
    dist_pow = dist_pow,
    .progress = .progress
  )
  knn_imp_results(fill_imp_parts(res), res$info)
}

# knn_imp() without its door: the parts that impute the columns the result
# holds, for a caller that places them itself (group_imp(), slide_imp()).
# `obj` is a double matrix the door has checked, and `scan` its col_scan()
# (or scan_cols() of the matrix it was sliced from), so nothing here reads
# `obj` to check it. Otherwise the same formals as knn_imp(), defaults
# included, so a caller passes it what it would pass knn_imp(); a test holds
# them equal. Returns a list:
# - `obj`: the matrix the parts index, the input as it was;
# - `parts`: the parts, in `obj`'s numbering, which together cover the
#   columns in `keep`, and may cover others; NULL when no kept column has a
#   hole;
# - `keep`: the columns the result holds, ascending, or NULL for every column;
# - `info`: the record knn_imp_results() attaches;
# - `n_miss`: the holes of each column of `obj`, for fill_parts().
knn_imp_parts <- function(
  obj,
  scan,
  k,
  colmax = 0.9,
  method = c("euclidean", "manhattan"),
  cores = 1,
  post_imp = TRUE,
  subset = NULL,
  dist_pow = 0,
  .progress = FALSE
) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_int(k, lower = 1, upper = ncol(obj) - 1, .var.name = "k")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  checkmate::assert_number(
    dist_pow,
    lower = 0,
    finite = TRUE,
    .var.name = "dist_pow"
  )
  checkmate::assert_flag(.progress, .var.name = ".progress")
  # the columns the result holds, ascending: NULL for every column. Only they
  # are imputed, but every eligible column is a candidate neighbor, kept or
  # not. Its refusals name the caller, knn_imp() for a user.
  keep <- resolve_narrow_subset(subset, obj, call = parent.frame())
  out_cols <- if (is.null(keep)) seq_len(ncol(obj)) else keep

  # per-column missingness, from the scan
  cmiss <- scan$n_miss

  # early exit if no column the result holds has a hole. Those columns come
  # back unchanged, but in the same shape as the main path's result, as
  # @returns promises: fill_imp_parts() takes them from `obj`.
  if (!any(cmiss[out_cols] > 0L)) {
    if (is.null(keep)) {
      cli::cli_inform(
        "No missing values in input. Returning input unchanged."
      )
    } else {
      cli::cli_inform(
        "No missing values in the {.arg subset} columns. Returning them unchanged."
      )
    }
    return(list(
      obj = obj,
      parts = NULL,
      keep = keep,
      info = knn_imp_info(post_imp = post_imp),
      n_miss = cmiss
    ))
  }

  miss_rate <- cmiss / nrow(obj)

  # partitioning: determine eligible columns under colmax threshold.
  # Inclusive, matching pca_imp() and the roxygen ("greater than colmax").
  eligible <- miss_rate <= min(colmax, 1)
  n_elig <- sum(eligible)

  if (k > n_elig - 1L) {
    cli::cli_abort(
      c(
        "{.arg k} ({k}) must be below the number of usable columns ({n_elig}).",
        "i" = "Reduce {.arg k} or relax {.arg colmax} to admit more columns."
      ),
      class = "slideimp_infeasible"
    )
  }

  # indices, which the groups and parts take: the scan carries no names
  eligible_idx <- which(eligible)
  has_miss_idx <- which(cmiss > 0L)
  eligible_has_miss <- intersect(eligible_idx, has_miss_idx)
  eligible_complete <- setdiff(eligible_idx, has_miss_idx)

  grp_impute <- sort(intersect(eligible_has_miss, out_cols))
  grp_miss_no_imp <- sort(setdiff(eligible_has_miss, out_cols))
  grp_complete <- sort(eligible_complete)

  if (length(grp_impute) == 0L) {
    cli::cli_abort(
      c(
        "No column to impute is within {.arg colmax} ({colmax}).",
        "i" = "Relax {.arg colmax} to admit columns with more missingness."
      ),
      class = "slideimp_infeasible"
    )
  }

  method <- switch(method, "euclidean" = 0L, "manhattan" = 1L)

  # impute using brute-force K-NN: `loc` holds the (row, col) of each missing
  # cell of grp_impute, in column-major order with columns `row` and `col`,
  # `value` its value, and `left` the rows of `loc` K-NN found no value for:
  # a column no candidate shares an observed row with, and a cell where none
  # of the chosen neighbors is observed. Under `post_imp` the kernel sets each
  # of those to its column's observed mean; otherwise they hold NA.
  imputed <- impute_knn_brute(
    obj = obj,
    k = k,
    grp_impute = grp_impute - 1L,
    grp_miss_no_imp = grp_miss_no_imp - 1L,
    grp_complete = grp_complete - 1L,
    method = method,
    dist_pow = dist_pow,
    fill_mean = post_imp,
    cores = cores,
    pb = .progress
  )
  # the kernel sets no dimnames, so `loc` gets a part map's column names here.
  dimnames(imputed$loc) <- list(NULL, c("row", "col"))
  loc <- imputed$loc
  left <- imputed$left

  # the parts, which together cover every column the result holds. The K-NN
  # parts cover grp_impute, which lies inside the kept columns, and the
  # eligible complete columns kept, which have nothing to fill. Under
  # `post_imp = FALSE` the columns with a cell K-NN left are a part of their
  # own, `left_na`, whose map holds the cells K-NN did fill, so that a
  # `left_na` part names only columns that keep a hole (check_parts()). The
  # columns over `colmax` the result holds are one more part, filled with
  # their means under `post_imp` and `left_na` without it. Each has a hole,
  # since its missing rate is above `colmax`, and an observed value, since
  # the door refuses a column with none (refuse_nonfinite()).
  complete_kept <- intersect(eligible_complete, out_cols)
  mean_cols <- setdiff(out_cols, eligible_idx)
  if (post_imp || length(left) == 0L) {
    knn_parts <- list(new_part(
      cols = sort(c(grp_impute, complete_kept)),
      map = loc,
      values = imputed$value,
      kind = "knn"
    ))
  } else {
    left_cols <- unique(loc[left, "col"])
    in_left_col <- loc[, "col"] %in% left_cols
    filled <- !in_left_col
    filled_left <- in_left_col
    filled_left[left] <- FALSE
    knn_parts <- list(
      new_part(
        cols = sort(c(setdiff(grp_impute, left_cols), complete_kept)),
        map = loc[filled, , drop = FALSE],
        values = imputed$value[filled],
        kind = "knn"
      ),
      new_part(
        cols = left_cols,
        map = loc[filled_left, , drop = FALSE],
        values = imputed$value[filled_left],
        kind = "knn",
        left_na = TRUE
      )
    )
    rm(in_left_col, filled, filled_left)
  }
  mean_part <- if (post_imp) {
    fallback_part(obj, mean_cols, "mean")
  } else {
    new_part(cols = mean_cols, kind = "mean", left_na = TRUE)
  }
  parts <- c(knn_parts, list(mean_part))
  rm(imputed, knn_parts, mean_part)

  # What `post_imp` set to a column mean, recorded on the result since nothing
  # in the values tells a column mean from a K-NN average afterwards: the
  # columns over `colmax` in `mean_filled_cols`, the cells K-NN left in
  # `mean_filled_cells`, both numbered as the result's columns.
  out_col <- function(cols) if (is.null(keep)) cols else match(cols, keep)
  info <- if (post_imp) {
    knn_imp_info(
      post_imp = post_imp,
      mean_filled_cols = out_col(mean_cols),
      mean_filled_cells = mean_filled_cells_matrix(
        row = loc[left, "row"],
        col = out_col(loc[left, "col"]),
        n_row = nrow(obj)
      )
    )
  } else {
    knn_imp_info(post_imp = post_imp)
  }
  rm(loc, left)

  list(obj = obj, parts = parts, keep = keep, info = info, n_miss = cmiss)
}

# What knn_imp() reports about a call: the attributes of its result after
# dim, dimnames and class, in the order it sets them. The defaults describe a
# call that filled nothing with a mean.
knn_imp_info <- function(
  post_imp,
  mean_filled_cols = integer(0),
  mean_filled_cells = mean_filled_cells_matrix()
) {
  list(
    imp_method = "knn",
    fallback = FALSE,
    post_imp = post_imp,
    mean_filled_cols = mean_filled_cols,
    mean_filled_cells = mean_filled_cells
  )
}

# The result of knn_imp(), from knn_imp_info().
knn_imp_results <- function(obj, info) {
  new_slideimp_results(
    obj,
    info$imp_method,
    fallback = info$fallback,
    post_imp = info$post_imp,
    mean_filled_cols = info$mean_filled_cols,
    mean_filled_cells = info$mean_filled_cells
  )
}
