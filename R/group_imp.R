#' Resolution of a Group Specification
#'
#' Convert a group specification to the data frame that [prep_groups()]
#' reads. This S3 generic is exported so that extension packages, such as
#' `slideimp.extra`, can register methods for other classes.
#'
#' @param x An object. The group specification: `slideimp` has a method for
#'   `data.frame`, and the optional `slideimp.extra` package adds one for a
#'   character platform name.
#'
#' @returns A data frame. The `data.frame` method returns `x` with duplicate
#'   rows removed.
#'
#' @examples
#' df <- data.frame(feature = c("cg1", "cg2"), group = c(1, 1))
#' slideimp_resolve_group(df)
#'
#' @export
slideimp_resolve_group <- function(x) {
  UseMethod("slideimp_resolve_group")
}

#' @rdname slideimp_resolve_group
#' @method slideimp_resolve_group data.frame
#' @export
slideimp_resolve_group.data.frame <- function(x) {
  # collapse::funique() does not take list columns, which the list-column form
  # of `group` is made of; the long format has none and takes the fast path.
  if (any(vapply(x, is.list, logical(1)))) unique(x) else collapse::funique(x)
}

#' @rdname slideimp_resolve_group
#' @method slideimp_resolve_group default
#' @export
slideimp_resolve_group.default <- function(x) {
  cli::cli_abort(c(
    "No {.fun slideimp_resolve_group} method for objects of class {.cls {class(x)[1]}}.",
    "i" = "Supply a data frame, or install and load {.pkg slideimp.extra} to pass a platform name such as {.val EPICv2}.",
    ">" = "See the {.pkg slideimp} README for installation instructions."
  ))
}

#' Feature Removal from Auxiliary Columns
#'
#' A grouped set difference with `collapse`: for each group, remove from
#' `aux` the names that are in the same group of `feat`.
#'
#' @param aux A list of character vectors. The auxiliary column names of each
#'   group.
#' @param feat A list of character vectors, the same length as `aux`. The
#'   feature column names of each group.
#' @param iter Integer. The group positions, used as the levels for
#'   [collapse::gsplit()].
#' @returns A list of character vectors of the same length as `aux`, with `feat`
#'   elements removed.
#'
#' @noRd
#' @keywords internal
remove_feat_from_aux <- function(aux, feat, iter) {
  aux_lens <- lengths(aux)
  if (sum(aux_lens) == 0L) {
    return(aux)
  }

  aux_flat <- unlist(aux, use.names = FALSE)
  aux_gid <- rep.int(iter, aux_lens)

  feat_lens <- lengths(feat)
  feat_flat <- unlist(feat, use.names = FALSE)
  feat_gid <- rep.int(iter, feat_lens)

  keep <- is.na(collapse::fmatch(
    list(aux_gid, aux_flat),
    list(feat_gid, feat_flat)
  ))
  gid_f <- factor(aux_gid[keep], levels = iter)
  out <- collapse::gsplit(aux_flat[keep], gid_f)
  names(out) <- NULL
  out
}

#' List Pruning to a Reference Set
#'
#' Keep, in each element of `lst`, only the values that are in `A`.
#'
#' @param lst A list of character vectors. The names to prune, such as the
#'   features or the auxiliary columns of each group.
#' @param A Character. The names to keep.
#' @param iter Integer. The group positions, used as the levels for
#'   [collapse::gsplit()].
#'
#' @returns A list of character vectors the same length as `lst`. Each element
#'   holds only the values of its `lst` element that are in `A`.
#'
#' @noRd
#' @keywords internal
prune_to_A <- function(lst, A, iter) {
  lens <- lengths(lst)
  if (sum(lens) == 0L) {
    return(lst)
  }
  flat <- unlist(lst, use.names = FALSE)
  keep <- !is.na(collapse::fmatch(flat, A))
  gid <- factor(rep.int(iter, lens), levels = iter)
  out <- collapse::gsplit(flat[keep], gid[keep])
  names(out) <- NULL
  out
}

# Drops the groups that `drop` marks from prep_groups()'s table, with a
# message that names them by label, or by row position where `group` has no
# labels. `why` is the reason, as cli markup, pasted after "dropped: " into
# the message template, so that cli renders its markup.
drop_groups <- function(group, drop, why) {
  has_label <- "group" %in% names(group)
  dropped <- if (has_label) as.character(group$group[drop]) else which(drop)
  cli::cli_inform(c(
    "!" = paste0("{length(dropped)} group{?s} dropped: ", why),
    "i" = if (has_label) {
      "{cli::qty(length(dropped))}Dropped group{?s}: {fmt_trunc(dropped)}"
    } else {
      "Dropped group indices: {fmt_trunc(dropped)}"
    }
  ))
  group <- group[!drop, , drop = FALSE]
  rownames(group) <- NULL
  group
}

# A group's slice of `obj`, the columns `col_idx`, with no dimnames: the
# backends address columns by number, and a slice sent to a mirai daemon
# would carry the row names with every group.
group_slice <- function(obj, col_idx) {
  sub_mat <- obj[, col_idx, drop = FALSE]
  dimnames(sub_mat) <- NULL
  sub_mat
}

# One group of group_imp(), on its slice: `sub_mat` holds the group's
# features as its first columns and its aux columns after them, and
# `col_idx` gives each slice column's column in group_imp()'s `obj`.
# `sub_scan` is the slice's scan, group_imp()'s indexed by `col_idx`, so the
# slice is not read to check it. `feat_out` gives each feature's column in
# group_imp()'s result. `imp_fn` is knn_imp_parts() or pca_imp_parts(), and
# `params` its arguments, with the features as `subset`. Returns a list:
# - `parts`: the group's parts over its features, in `obj`'s numbering;
# - `fallback`: TRUE when the backend found the group infeasible;
# - `mean_filled`: what `post_imp` set to a column mean, numbered as the
#   result's columns, from mean_filled_record(); NULL for a fallback, which
#   records nothing.
# A group the backend finds infeasible goes through `on_infeasible`: under
# "skip" it has no part, so its holes stay NA, and under "mean" one mean
# part over its features with a hole. Under "error" the condition is
# returned in place of the list, for the caller to raise.
#
# It runs on the mirai daemons too, where it reaches the backend through the
# package namespace the daemon has loaded.
impute_group <- function(
  sub_mat,
  sub_scan,
  col_idx,
  feat_out,
  params,
  imp_fn,
  on_infeasible
) {
  res <- tryCatch(
    suppressMessages(do.call(
      imp_fn,
      c(list(obj = sub_mat, scan = sub_scan), params)
    )),
    slideimp_infeasible = function(e) e
  )
  if (inherits(res, "slideimp_infeasible")) {
    if (on_infeasible == "error") {
      return(res)
    }
    parts <- list()
    if (on_infeasible == "mean") {
      feat <- seq_along(feat_out)
      holed <- feat[sub_scan$n_miss[feat] > 0L]
      parts <- list(fallback_part(sub_mat, holed, "mean"))
    }
    return(list(
      parts = rebase_parts(parts, col_idx),
      fallback = TRUE,
      mean_filled = NULL
    ))
  }

  # the backend's result columns are the features, in order: `subset` is
  # the first columns of the slice, and `res$keep` is NULL only when the
  # slice holds nothing else. So its `t`-th result column is `feat_out[t]`.
  list(
    parts = rebase_parts(res$parts, col_idx, res$keep),
    fallback = FALSE,
    mean_filled = mean_filled_record(res$info, feat_out)
  )
}

#' Group Preparation for Imputation
#'
#' Normalize and validate a grouping specification for use with [group_imp()].
#'
#' `prep_groups()` converts a long-format or list-column `group` to a
#' `slideimp_tbl`, and pads small groups when `min_group_size` is above `0`.
#'
#' @inheritParams group_imp
#' @param obj_cn Character. The column names of the data matrix,
#'   usually `colnames(obj)`.
#' @param subset Optional character vector. The features to impute, by name,
#'   or `NULL` for every grouped feature. A group's other features become its
#'   auxiliary columns, and a group left with no feature is dropped with a
#'   message.
#' @param min_group_size Integer. Minimum number of columns per group,
#'   features and auxiliary columns together. Groups smaller than this are
#'   padded with columns sampled at random from the other grouped columns in
#'   `obj_cn`. If `0`, no padding is performed.
#' @param seed Integer, numeric, or `NULL`. Random seed used when padding small
#'   groups.
#'
#' @details
#' Let \eqn{A} be `obj_cn` and \eqn{B} be the union of all feature and
#' auxiliary names in `group`. When `allow_unmapped = FALSE`, `prep_groups()`
#' requires \eqn{A \subseteq B}: every name in `obj_cn` must appear in
#' `group`.
#'
#' Names in \eqn{B} but not in \eqn{A}, such as probes dropped in quality
#' control, are removed from each group. Groups left with no features are
#' dropped with a message.
#'
#' @returns A data frame of class `slideimp_tbl` containing:
#' * `group`: original group labels, if provided, or sequential group labels.
#' * `feature`: a list-column of character vectors containing feature names.
#' * `aux`: a list-column of character vectors containing auxiliary names.
#' * `parameters`: a list-column of per-group configuration lists.
#'
#' @seealso [group_imp()] to impute with the prepared groups.
#'
#' @examples
#' sim <- sim_mat(n = 10, p = 20)
#' prepped <- prep_groups(colnames(sim$input), sim$col_group)
#' prepped
#'
#' @export
prep_groups <- function(
  obj_cn,
  group,
  subset = NULL,
  min_group_size = 0,
  allow_unmapped = FALSE,
  seed = NULL
) {
  feature <- NULL
  # pre-conditioning ---
  checkmate::assert_character(
    obj_cn,
    min.len = 1,
    unique = TRUE,
    any.missing = FALSE,
    .var.name = "obj_cn"
  )
  group <- slideimp_resolve_group(group)
  checkmate::assert_data_frame(group, min.rows = 1, .var.name = "group")
  checkmate::assert_names(
    colnames(group),
    must.include = "feature",
    .var.name = "group"
  )
  checkmate::assert_int(min_group_size, lower = 0, .var.name = "min_group_size")
  checkmate::assert_flag(allow_unmapped, .var.name = "allow_unmapped")
  checkmate::assert_number(seed, null.ok = TRUE, .var.name = "seed")
  checkmate::assert_character(
    subset,
    min.len = 1,
    any.missing = FALSE,
    unique = TRUE,
    null.ok = TRUE,
    .var.name = "subset"
  )

  # Step 1: Normalize input format ----
  if ("group" %in% names(group) && inherits(group$feature, "character")) {
    if (anyNA(group$group)) {
      cli::cli_abort(c(
        "{.code NA} is not allowed in {.code group$group}.",
        "i" = "Give every feature a label in {.code group$group}."
      ))
    }

    checkmate::assert_character(
      group$feature,
      any.missing = FALSE,
      unique = TRUE,
      .var.name = "group$feature"
    )

    if (is.numeric(group$group) && !checkmate::test_integerish(group$group)) {
      cli::cli_warn(c(
        "{.code group$group} has non-integer numeric values.",
        "i" = "Each distinct value forms one group.",
        "i" = "Use integer or character labels."
      ))
    }

    group <- collapse::fsummarize(
      collapse::fgroup_by(group, group),
      feature = list(feature)
    )
  } else if (inherits(group$feature, "character")) {
    cli::cli_abort(c(
      "{.arg group} has a character {.field feature} column but no {.field group} column.",
      "i" = "Add a {.field group} column, or make {.field feature} a list-column with one element per group."
    ))
  }

  checkmate::assert_list(
    group$feature,
    types = "character",
    min.len = 1,
    unique = TRUE,
    .var.name = "group$feature"
  )

  # Normalize aux - always present as a list of character vectors
  if ("aux" %in% names(group)) {
    checkmate::assert_list(
      group$aux,
      types = c("character", "logical", "null"),
      min.len = 1,
      .var.name = "group$aux"
    )
    group$aux <- lapply(group$aux, function(x) {
      if (is.null(x) || is.logical(x)) character(0) else unique(x[!is.na(x)])
    })
  } else {
    group$aux <- replicate(nrow(group), character(0), simplify = FALSE)
  }

  # Normalize parameters - always present as a list of lists
  if ("parameters" %in% names(group)) {
    checkmate::assert_list(
      group$parameters,
      types = c("list", "null"),
      .var.name = "group$parameters"
    )
    group$parameters <- lapply(group$parameters, function(p) {
      if (is.null(p)) list() else as.list(p)
    })
  } else {
    group$parameters <- replicate(nrow(group), list(), simplify = FALSE)
  }

  # Step 2: Set validation ----
  A <- obj_cn
  all_feats <- unlist(group$feature)
  B <- unique(c(all_feats, unlist(group$aux)))

  # A must be a subset of B
  not_in_B <- setdiff(A, B)
  if (length(not_in_B) > 0) {
    if (!allow_unmapped) {
      cli::cli_abort(c(
        "{length(not_in_B)} column{?s} {?has/have} no matching entry in {.arg group}: {fmt_trunc(not_in_B)}",
        "i" = "{cli::qty(length(not_in_B))}Add {?it/them} to {.arg group}, or set {.code allow_unmapped = TRUE} to leave {?it/them} untouched."
      ))
    } else {
      cli::cli_inform(c(
        "!" = "{length(not_in_B)} column{?s} {?has/have} no matching entry in {.arg group} and will be left untouched.",
        "i" = "{fmt_trunc(not_in_B)}"
      ))
    }
  }

  # Features must not appear in more than one group
  if (anyDuplicated(all_feats) > 0) {
    dups <- unique(all_feats[duplicated(all_feats)])
    cli::cli_abort(c(
      "{length(dups)} feature{?s} appear{?s/} in more than one group:",
      "x" = "{fmt_trunc(dups)}",
      "i" = "Assign each feature to one group."
    ))
  }

  # Prune B \ A (QC-dropped probes not in obj)
  iter <- seq_len(nrow(group))
  group$feature <- prune_to_A(group$feature, A = A, iter = iter)
  group$aux <- prune_to_A(group$aux, A = A, iter = iter)
  group$aux <- remove_feat_from_aux(group$aux, group$feature, iter = iter)

  # Step 3: Prune empty groups ----
  empty <- lengths(group$feature) == 0L
  if (any(empty)) {
    group <- drop_groups(group, empty, "no feature matches a column name.")
  }
  if (nrow(group) == 0) {
    cli::cli_abort(c(
      "No group has a feature that matches a column name.",
      "i" = "Check that {.code group$feature} holds the column names of the data."
    ))
  }

  # Step 3b: Apply subset, demote non-subset features to aux ----
  if (!is.null(subset)) {
    bad_cols <- setdiff(subset, A)
    if (length(bad_cols) > 0) {
      cli::cli_abort(c(
        "{length(bad_cols)} name{?s} in {.arg subset} {?does/do} not match a column name:",
        "x" = "{fmt_trunc(bad_cols)}"
      ))
    }
    all_feats_now <- unlist(group$feature)
    bad_feats <- setdiff(subset, all_feats_now)
    if (length(bad_feats) > 0) {
      cli::cli_abort(c(
        "{length(bad_feats)} feature{?s} in {.arg subset} {?is/are} not assigned to any group:",
        "x" = "{fmt_trunc(bad_feats)}",
        "i" = "{cli::qty(length(bad_feats))}Add {?it/them} to a group, or remove {?it/them} from {.arg subset}."
      ))
    }
    # Features in subset stay; the rest are demoted to aux. One hashed match
    # over all groups: a per-group intersect() rehashes `subset` every time.
    iter <- seq_len(nrow(group))
    feat_kept <- prune_to_A(group$feature, A = subset, iter = iter)
    demoted <- remove_feat_from_aux(group$feature, feat_kept, iter = iter)
    group$aux <- Map(c, group$aux, demoted)
    group$feature <- feat_kept
    keep <- lengths(group$feature) > 0L
    if (any(!keep)) {
      group <- drop_groups(
        group,
        !keep,
        "no features remaining after applying {.arg subset}."
      )
    }
  }

  # Step 4: Pad groups if min_group_size > 0 ----
  if (min_group_size > 0) {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    group_size <- lengths(group$feature) + lengths(group$aux)
    need <- pmax(min_group_size - group_size, 0L)
    if (any(need > 0)) {
      # only columns `group` maps. A column it leaves out is left untouched,
      # and that includes not lending it to another group as an aux column.
      mapped <- A[A %in% B]
      group$aux <- Map(
        function(feat, aux, n) {
          if (n == 0L) {
            return(aux)
          }
          pool <- setdiff(mapped, c(feat, aux))
          if (length(pool) < n) {
            cli::cli_abort(c(
              "{.arg min_group_size} is too large: {.arg group} maps {length(mapped)} column{?s}.",
              "i" = "Set {.arg min_group_size} to at most {length(mapped)}."
            ))
          }
          c(aux, sample(pool, size = n))
        },
        group$feature,
        group$aux,
        need
      )
    }
  }

  # Clean up list names
  for (col in names(group)) {
    if (is.list(group[[col]])) names(group[[col]]) <- NULL
  }
  # Add group names for diagnostic
  if (!"group" %in% names(group)) {
    group$group <- seq_len(nrow(group))
  }
  class(group) <- c("slideimp_tbl", "data.frame")
  return(group)
}

#' Grouped K-NN or PCA Imputation
#'
#' Perform K-NN or PCA imputation independently within feature groups.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#'
#' @param group A data frame or a character scalar. The feature groups, in one
#'   of these forms:
#'
#'   * A character scalar naming a supported Illumina platform, such as
#'     `"EPICv2"`. This needs the `slideimp.extra` package.
#'   * A long-format data frame with columns `group` and `feature`.
#'   * A list-column data frame with a `feature` list-column. Optional
#'     list-columns are `aux`, for auxiliary feature names, and `parameters`,
#'     for group-specific parameter lists.
#'
#' @param subset Optional character vector. The features to impute and
#'   return, by name, or `NULL` to impute every grouped feature and return
#'   every column. A group's other features become its auxiliary columns, and
#'   a group left with no feature is dropped with a message.
#' @param allow_unmapped Logical. If `FALSE`, every column in `colnames(obj)`
#'   must appear in `group`. If `TRUE`, columns with no group assignment are
#'   left untouched and are not used as auxiliary columns.
#' @param k Integer or `NULL`. Number of nearest neighbors for K-NN imputation.
#'   If `NULL`, `k` may be supplied through `group$parameters`.
#' @param ncp Integer or `NULL`. Number of components for PCA imputation. If
#'   `NULL`, `ncp` may be supplied through `group$parameters`.
#' @param method Character or `NULL`. One of `"euclidean"` or `"manhattan"`
#'   for K-NN imputation, or `"regularized"` or `"EM"` for PCA imputation.
#'   `NULL` uses `"euclidean"` for K-NN and `"regularized"` for PCA.
#' @param cores Integer. Number of cores for K-NN imputation only, as
#'   described under Parallelization.
#' @param .progress Logical. If `TRUE`, show progress.
#' @param min_group_size Integer or `NULL`. Minimum total number of columns
#'   per group, counting both features and auxiliary columns. Groups smaller
#'   than this are padded with columns sampled at random from the other
#'   grouped columns of `obj`.
#' @param seed Integer, numeric, or `NULL`. Random seed used when padding
#'   small groups and, for PCA imputation, passed to [pca_imp()].
#' @param pin_blas Logical. If `TRUE`, pin BLAS threads to 1 on the `mirai`
#'   daemons in a parallel run. Defaults to `TRUE`.
#' @param on_infeasible Character. One of `"error"`, `"skip"`, or `"mean"`.
#'   The action for a group that is infeasible for imputation, for example
#'   because `k` or `ncp` is too large for the columns that `colmax` leaves.
#'   Defaults to `"error"`.
#'
#' @details
#' Imputing each group on its own can reduce the runtime for large matrices.
#'
#' Specify `k` and related arguments to use K-NN imputation, or `ncp` and
#' related arguments to use PCA imputation. If both `k` and `ncp` are `NULL`,
#' `group$parameters` must supply either `k` or `ncp` for every group. Passing
#' an argument of the other method, such as `scale` with `k`, is an error.
#'
#' Group-specific parameters in `group$parameters` take priority over global
#' arguments, and global arguments fill in the parameters a group does not
#' set. All groups must use the same imputation method. A global argument
#' shared with [knn_imp()] or [pca_imp()] that is `NULL` takes that function's
#' default.
#'
#' Per-group `k` is capped at the group's column count minus one, counting
#' features and auxiliary columns. Per-group `ncp` is capped at the smaller of
#' that cap and `nrow(obj) - 2`. `group_imp()` warns when it caps a value.
#'
#' A group is infeasible when it is too small for any valid `k` or `ncp`, or
#' when [knn_imp()] or [pca_imp()] cannot impute it, for example because
#' `colmax` leaves too few columns. With `on_infeasible = "skip"`, an
#' infeasible group keeps its missing values. With `"mean"`, column means fill
#' them, and with `"error"`, `group_imp()` stops with an error. A group whose
#' features have no missing values is not imputed and is never infeasible. An
#' infeasible group adds nothing to the `mean_filled_cols` and
#' `mean_filled_cells` attributes of the result.
#'
#' With `subset`, the result holds only the `subset` columns, in the order
#' they have in `obj`. The other features still take part in the imputation,
#' as auxiliary columns of their groups.
#'
#' The `mirai` daemons run R's `"L'Ecuyer-CMRG"` generator. So with
#' `nb.init > 1`, a parallel PCA run draws different random initializations
#' from a sequential run under R's default generator. The two agree when the
#' session also runs `RNGkind("L'Ecuyer-CMRG")`.
#'
#' @inheritSection slideimp-package Column mean fill
#'
#' @section Parallelization:
#' * K-NN: use the `cores` argument. If `mirai` daemons are active, `cores`
#'   is set to `1` to avoid nested parallelism.
#' * PCA: use `mirai::daemons()` instead of `cores`.
#' * With `mirai` daemons, each daemon runs BLAS on one thread unless
#'   `pin_blas = FALSE`.
#'
#' @inheritSection slideimp-package Missing values and non-finite input
#'
#' @inheritSection pca_imp PCA Performance tips
#'
#' @note
#' A platform name for `group`, such as `"EPICv2"` or `"EPICv2_deduped"`, needs
#' the optional `slideimp.extra` package. Its `slideimp_arrays` object lists
#' the supported platforms.
#'
#' @returns A `slideimp_results` object: a numeric matrix with the rows of
#' `obj` and the columns that `subset` selects, with the missing values of
#' the grouped features imputed.
#'
#' @seealso
#' * [prep_groups()] to inspect and edit the groups and their parameters.
#' * [knn_imp()] and [pca_imp()] to impute a single matrix.
#'
#' @examples
#' set.seed(1234)
#' to_test <- sim_mat(10, 20, perc_total_na = 0.05, perc_col_na = 1)
#' obj <- to_test$input
#' group <- to_test$col_group
#' head(group)
#'
#' # Simple grouped K-NN imputation
#' results <- group_imp(obj, group = group, k = 2, .progress = FALSE)
#' results
#'
#' # Impute only a subset of features
#' subset_features <- sample(to_test$col_group$feature, size = 10)
#' knn_subset <- group_imp(
#'   obj,
#'   group = group,
#'   subset = subset_features,
#'   k = 2,
#'   .progress = FALSE
#' )
#'
#' # Use prep_groups() to inspect and edit per-group parameters
#' prepped <- prep_groups(colnames(obj), group)
#' prepped$parameters <- lapply(seq_len(nrow(prepped)), function(i) list(k = 2))
#' prepped$parameters[[2]]$k <- 4
#' knn_grouped <- group_imp(obj, group = prepped, .progress = FALSE)
#'
#' @examplesIf interactive() && requireNamespace("mirai", quietly = TRUE)
#' # PCA imputation with mirai parallelism
#' mirai::daemons(2)
#' pca_grouped <- group_imp(obj, group = group, ncp = 2)
#' mirai::daemons(0)
#' pca_grouped
#'
#' @export
group_imp <- function(
  obj,
  group,
  subset = NULL,
  allow_unmapped = FALSE,
  k = NULL,
  ncp = NULL,
  method = NULL,
  # K-NN arguments
  cores = 1,
  .progress = TRUE,
  min_group_size = NULL,
  colmax = NULL,
  post_imp = NULL,
  dist_pow = NULL,
  # PCA arguments
  scale = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  row.w = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL,
  solver = NULL,
  lobpcg_control = NULL,
  clamp = NULL,
  pin_blas = TRUE,
  on_infeasible = c("error", "skip", "mean")
) {
  checkmate::assert_matrix(
    obj,
    mode = "double",
    col.names = "unique",
    null.ok = FALSE,
    .var.name = "obj"
  )
  # the variances serve only PCA, and a global `k` makes every group K-NN or
  # is refused below. Each group's scan is this one indexed (scan_cols()).
  scan <- col_scan(obj, var = is.null(k))
  refuse_nonfinite(scan)
  checkmate::assert_flag(allow_unmapped, .var.name = "allow_unmapped")
  checkmate::assert_int(k, lower = 1, null.ok = TRUE, .var.name = "k")
  checkmate::assert_int(ncp, lower = 1, null.ok = TRUE, .var.name = "ncp")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_flag(.progress, .var.name = ".progress")
  checkmate::assert_flag(pin_blas, null.ok = FALSE, .var.name = "pin_blas")
  on_infeasible <- match.arg(on_infeasible)
  if (!is.null(solver)) {
    solver <- match.arg(solver, c("auto", "exact", "lobpcg"))
  }
  cn <- colnames(obj)
  # per-column missing counts, used further down to leave out the groups with
  # nothing to impute. `obj` stays the caller's matrix, attributes and all,
  # until fill_parts() copies it once at the end: the groups read slices of
  # it, and nothing writes into it, so the scan holds for it throughout.
  cmiss <- scan$n_miss

  # Step 1: Build canonical groups via prep_groups() ----
  # After this call, group$feature, group$aux, and group$parameters are all
  # guaranteed to exist as properly typed list-columns.
  # subset is also handled here - non-subset features are demoted to aux.
  group <- prep_groups(
    obj_cn = cn,
    group = group,
    subset = subset,
    allow_unmapped = allow_unmapped,
    min_group_size = if (is.null(min_group_size)) 0L else min_group_size,
    seed = seed
  )
  feat_lengths <- lengths(group$feature)
  aux_lengths <- lengths(group$aux)
  # messages name a group by its label, the way the result and print() do
  group_lab <- as.character(group$group)
  # the columns the result holds, ascending: NULL for every column.
  # prep_groups() has refused a `subset` that is empty, names a column twice,
  # or names one that is not in `obj` or in no group, so every `subset`
  # column is a feature of one group, and none of resolve_narrow_subset()'s
  # refusals can fire.
  result_cols <- resolve_narrow_subset(subset, obj)

  # Step 2: Resolve parameters ----
  if (!is.null(k) && !is.null(ncp)) {
    cli::cli_abort(
      "Supply at most one of {.arg k} (K-NN imputation) or {.arg ncp} (PCA imputation)."
    )
  }

  # the names each group set for itself, taken before a global fills any gap,
  # so that an unknown name is reported against `group$parameters` only when
  # it came from there
  own_param_names <- unique(unlist(
    lapply(group$parameters, names),
    use.names = FALSE
  ))

  # global fills gaps (group-wise wins)
  global_params <- list(
    k = k,
    method = method,
    colmax = colmax,
    post_imp = post_imp,
    dist_pow = dist_pow,
    ncp = ncp,
    scale = scale,
    coeff.ridge = coeff.ridge,
    threshold = threshold,
    row.w = row.w,
    seed = seed,
    nb.init = nb.init,
    maxiter = maxiter,
    miniter = miniter,
    lobpcg_control = lobpcg_control,
    solver = solver,
    clamp = clamp
  )
  global_params <- global_params[!vapply(global_params, is.null, logical(1))]
  fill_params <- function(parameters, fill) {
    lapply(parameters, function(p) {
      for (nm in names(fill)) {
        if (is.null(p[[nm]])) p[[nm]] <- fill[[nm]]
      }
      p
    })
  }

  # `k` and `ncp` first: they decide the backend, and the backend decides
  # which of the other globals apply
  group$parameters <- fill_params(
    group$parameters,
    global_params[intersect(names(global_params), c("k", "ncp"))]
  )

  # validate: each group has exactly one of k or ncp
  has_k <- vapply(group$parameters, \(p) "k" %in% names(p), logical(1))
  has_ncp <- vapply(group$parameters, \(p) "ncp" %in% names(p), logical(1))

  if (any(has_k & has_ncp)) {
    bad <- group_lab[has_k & has_ncp]
    cli::cli_abort(c(
      "{length(bad)} group{?s} {?has/have} both {.arg k} and {.arg ncp}:",
      "x" = "{fmt_trunc(bad)}",
      "i" = "Set one of {.arg k} or {.arg ncp} for each group, as a global argument or in {.code group$parameters}."
    ))
  }
  if (any(!has_k & !has_ncp)) {
    bad <- group_lab[!has_k & !has_ncp]
    cli::cli_abort(c(
      "{length(bad)} group{?s} {?has/have} neither {.arg k} nor {.arg ncp}:",
      "x" = "{fmt_trunc(bad)}",
      "i" = "Supply {.arg k} or {.arg ncp}, or set one of them in {.code group$parameters}."
    ))
  }
  if (any(has_k) && any(has_ncp)) {
    cli::cli_abort(c(
      "Some groups use {.arg k} and others use {.arg ncp}.",
      "i" = "Give every group {.arg k} for K-NN imputation, or every group {.arg ncp} for PCA imputation."
    ))
  }

  imp_method <- if (all(has_k)) "knn" else "pca"
  is_knn_mode <- imp_method == "knn"
  reg <- .param_registry[[imp_method]]

  # a global argument that belongs to the other backend. `seed` is exempt: it
  # also seeds the padding in prep_groups(), which both backends go through,
  # and K-NN has no other use for it.
  wrong_backend <- setdiff(names(global_params), c(reg$allowed, "seed"))
  if (length(wrong_backend) > 0L) {
    other_param <- if (is_knn_mode) "ncp" else "k"
    cli::cli_abort(c(
      "{cli::qty(length(wrong_backend))}Argument{?s} {.arg {wrong_backend}} {cli::qty(length(wrong_backend))}{?does/do} not apply to {.strong {method_label(imp_method)}} imputation.",
      "i" = "{cli::qty(length(wrong_backend))}Remove {?it/them}, or use {.arg {other_param}} to select the other method."
    ))
  }
  group$parameters <- fill_params(
    group$parameters,
    global_params[intersect(names(global_params), reg$allowed)]
  )

  # validate method values
  valid_methods <- reg$methods
  bad_method <- vapply(
    group$parameters,
    function(p) {
      !is.null(p$method) && !(p$method %in% valid_methods)
    },
    logical(1)
  )

  if (any(bad_method)) {
    bad <- group_lab[bad_method]
    cli::cli_abort(c(
      "{length(bad)} group{?s} {?has/have} an invalid {.arg method} for {method_label(imp_method)} imputation:",
      "x" = "{fmt_trunc(bad)}",
      "i" = "Set {.arg method} to {.val {cli::cli_vec(valid_methods, list(\"vec-last\" = \" or \"))}}."
    ))
  }

  # validate parameter names
  check_unknown_params(
    param_names = own_param_names,
    mode = imp_method,
    arg = "group$parameters"
  )

  # cap per-group k/ncp
  group_size <- feat_lengths + aux_lengths
  required_param <- reg$required
  cap <- group_size - 1L
  if (!is_knn_mode) {
    cap <- pmin(cap, nrow(obj) - 2L)
  }
  # a group that admits no `k` or `ncp` at all: one column, or for PCA an
  # `obj` of fewer than three rows. Such a group is infeasible, not an input
  # error, and `on_infeasible` decides what happens to it further down.
  too_small <- cap < 1L

  for (i in seq_len(nrow(group))) {
    p <- group$parameters[[i]]
    # the value is compared next, and a comparison coerces: a character "3"
    # is greater than an integer 11
    checkmate::assert_int(
      p[[required_param]],
      lower = 1,
      .var.name = sprintf("group$parameters[[%d]]$%s", i, required_param)
    )
    if (!too_small[i] && p[[required_param]] > cap[i]) {
      cli::cli_warn(
        "Group {group_lab[[i]]}: {.arg {required_param}} capped from {p[[required_param]]} to {cap[i]} (group size = {group_size[i]})."
      )
      p[[required_param]] <- cap[i]
      group$parameters[[i]] <- p
    }
  }

  # Step 3: Imputation loop ----
  # column-index lookups
  iter <- seq_len(nrow(group))

  # each group's columns of `obj`, by position, from the names in `lst`
  split_pos <- function(lst, lens) {
    out <- collapse::gsplit(
      collapse::fmatch(unlist(lst), cn),
      factor(rep.int(iter, lens), levels = iter)
    )
    names(out) <- NULL
    out
  }
  feat_splits <- split_pos(group$feature, feat_lengths)
  # the columns of each group's slice, its features then its aux columns
  group_cols <- Map(c, feat_splits, split_pos(group$aux, aux_lengths))

  # a group whose feature columns hold no missing value needs no imputation,
  # so it is not sent to the backend at all, and a group too small to impute
  # is infeasible only when it needs one. The backends would return no parts
  # for such a group themselves, through their `subset` early exit.
  needs_imp <- vapply(feat_splits, function(f) any(cmiss[f] > 0L), logical(1))

  # a group too small for any `k` or `ncp` goes through `on_infeasible`, like
  # a group the backend refuses. The backends cannot report it themselves:
  # they reject a one-column matrix as invalid input, not as infeasible.
  small_now <- too_small & needs_imp
  if (on_infeasible == "error" && any(small_now)) {
    i <- which(small_now)[[1L]]
    cli::cli_abort(
      c(
        "Group {group_lab[[i]]} is infeasible for imputation.",
        "x" = if (is_knn_mode) {
          "K-NN imputation needs at least 2 columns. The group has {group_size[i]} column{?s}."
        } else {
          "PCA imputation needs at least 3 rows and 2 columns. {.arg obj} has {nrow(obj)} row{?s}, and the group has {group_size[i]} column{?s}."
        },
        "i" = "Pad small groups with {.arg min_group_size}, or set {.arg on_infeasible} to {.val skip} or {.val mean}."
      ),
      class = "slideimp_infeasible"
    )
  }

  cli::cli_inform(
    "Imputing {nrow(group)} group{?s} using {.strong {method_label(imp_method)}}."
  )

  # parallelism resolution
  parallelize <- mirai_daemons_set()

  if (is_knn_mode) {
    if (cores > 1 && parallelize) {
      cli::cli_inform(c(
        "{.arg cores} = {cores} is ignored, to avoid nested parallelism with the {.pkg mirai} daemons.",
        "i" = "To use {cores} threads in place of the daemons, call {.code mirai::daemons(0)} first."
      ))
      cores <- 1
    }
  } else if (cores > 1 && !parallelize) {
    cli::cli_inform(c(
      "{.arg cores} = {cores} is ignored for PCA imputation.",
      "i" = "Use {.code mirai::daemons()} to impute the groups in parallel."
    ))
    cores <- 1
  }

  # build per-group call parameters. Both backends take the group's features
  # as `subset`, the first columns of the group's slice (`group_cols`), and
  # impute_group() keeps their parts over those columns only. `subset` is
  # not in the parameter registry, so no group can have set it.
  params <- lapply(iter, function(i) {
    p <- group$parameters[[i]]
    p$.progress <- FALSE
    p$subset <- seq_len(feat_lengths[[i]])

    if (is_knn_mode) {
      p$cores <- cores
    }

    p
  })

  imp_fn <- if (is_knn_mode) knn_imp_parts else pca_imp_parts

  if (parallelize) {
    cli::cli_inform("Running mode: {.strong mirai}")
  } else if (is_knn_mode && cores > 1) {
    cli::cli_inform("Running mode: {.strong threaded} ({cores} cores)")
  } else {
    cli::cli_inform("Running mode: {.strong sequential}")
  }

  # imputation. Each group hands back parts over its features, in `obj`'s
  # numbering, and they are filled into one copy of `obj` at the end; no
  # group's result is a matrix. Groups too small to impute are settled here,
  # in either running mode. `on_infeasible` is "skip" or "mean" at this
  # point: "error" aborted above. Under "skip" a group has no part, so its
  # holes stay NA.
  fallback_flags <- small_now
  # the parts of each group, a list per group, NULL where there are none
  group_parts <- vector("list", nrow(group))
  # what `post_imp` set to a column mean inside each group that ran, numbered
  # as the result's columns, and the `post_imp` the group ran with. NULL for a
  # group that never reached the backend.
  mean_filled <- vector("list", nrow(group))
  # the result's column of each group's features: their `obj` columns, or
  # with `subset` their positions in `result_cols`. Every feature is a
  # `subset` column, since prep_groups() demotes the others to aux.
  feat_out <- if (is.null(result_cols)) {
    feat_splits
  } else {
    lapply(feat_splits, match, table = result_cols)
  }
  if (on_infeasible == "mean") {
    for (i in which(small_now)) {
      f <- feat_splits[[i]]
      group_parts[[i]] <- list(fallback_part(
        obj,
        sort(f[cmiss[f] > 0L]),
        "mean"
      ))
    }
  }
  # the groups that go to the backend
  run <- which(needs_imp & !too_small)

  if (parallelize && length(run) > 0L) {
    # one task per group, each holding its own slice, sent by value, with no
    # shared memory; each hands back its parts, sized by its holes.
    # mirai sends the function with every task, so anything in the crate
    # travels once per group: with all groups' indices in it the payload
    # grows as groups^2.
    tasks <- lapply(run, function(i) {
      list(
        sub_mat = group_slice(obj, group_cols[[i]]),
        sub_scan = scan_cols(scan, group_cols[[i]]),
        col_idx = group_cols[[i]],
        feat_out = feat_out[[i]],
        params = params[[i]]
      )
    })

    crated_fn <- .crate(
      function(task) {
        if (pin_blas) {
          RhpcBLASctl::blas_set_num_threads(1)
          RhpcBLASctl::omp_set_num_threads(1)
        }
        out <- impute_group(
          task$sub_mat,
          task$sub_scan,
          task$col_idx,
          task$feat_out,
          task$params,
          imp_fn,
          on_infeasible
        )
        # handed back as a value. An error raised here reaches the host as
        # text, without the class a caller matches on, so the host raises it.
        # Only the message travels: the condition's call holds the slice.
        if (inherits(out, "slideimp_infeasible")) {
          out <- structure(
            list(message = conditionMessage(out)),
            class = "slideimp_infeasible_group"
          )
        }
        out
      },
      impute_group = impute_group,
      imp_fn = imp_fn,
      pin_blas = pin_blas,
      on_infeasible = on_infeasible
    )
    mapped <- mirai::mirai_map(tasks, crated_fn)
    # the slices have been sent; the host need not hold them while the
    # daemons work
    rm(tasks)
    res <- collect_mirai_map(mapped, .progress)
    rm(mapped)

    # A daemon that fails returns an `errorValue` in place of the group's
    # list. Abort before the write-back: the sequential path lets a
    # non-`slideimp_infeasible` error propagate, and collecting silently here
    # would hand back a matrix whose holes for that group are all NA.
    failed <- which(vapply(res, mirai::is_error_value, logical(1)))
    if (length(failed)) {
      first_msg <- mirai_error_text(res[[failed[[1L]]]])
      # impute_group() reaches slideimp's internals through the namespace
      # each daemon loaded, so a daemon on another build of slideimp fails
      # on an internal it lacks. R words that error in the session's
      # language, so the hint is shown only where it reads in English.
      build_hint <- if (
        grepl("could not find function", first_msg, fixed = TRUE)
      ) {
        "The {.pkg mirai} daemons may load another build of {.pkg slideimp} than this session. Restart them after installing this one."
      }
      cli::cli_abort(c(
        "Imputation failed in {length(failed)} group{?s}: {fmt_trunc(group_lab[run[failed]])}.",
        "x" = "Group {group_lab[[run[[failed[[1L]]]]]]}: {first_msg}",
        "i" = build_hint,
        "i" = "Call {.code mirai::daemons(0)} and run again to see the original error."
      ))
    }
    infeasible <- which(vapply(
      res,
      inherits,
      logical(1),
      what = "slideimp_infeasible_group"
    ))
    if (length(infeasible)) {
      first_msg <- res[[infeasible[[1L]]]]$message
      cli::cli_abort(
        c(
          "Group {group_lab[[run[[infeasible[[1L]]]]]]} is infeasible for imputation.",
          "x" = "{first_msg}",
          "i" = "Set {.arg on_infeasible} to {.val skip} or {.val mean}."
        ),
        class = "slideimp_infeasible"
      )
    }

    fallback_flags[run] <- vapply(res, `[[`, logical(1), "fallback")
    mean_filled[run] <- lapply(res, `[[`, "mean_filled")
    group_parts[run] <- lapply(res, `[[`, "parts")
    rm(res)
  } else {
    if (.progress && length(run) > 0L) {
      pb <- cli::cli_progress_bar(name = "Imputing", total = length(run))
    }
    for (i in run) {
      col_idx <- group_cols[[i]]
      out <- impute_group(
        group_slice(obj, col_idx),
        scan_cols(scan, col_idx),
        col_idx,
        feat_out[[i]],
        params[[i]],
        imp_fn,
        on_infeasible
      )
      if (inherits(out, "slideimp_infeasible")) {
        cli::cli_abort(
          c(
            "Group {group_lab[[i]]} is infeasible for imputation.",
            "i" = "Set {.arg on_infeasible} to {.val skip} or {.val mean}."
          ),
          parent = out,
          class = "slideimp_infeasible"
        )
      }
      group_parts[[i]] <- out$parts
      fallback_flags[i] <- out$fallback
      mean_filled[i] <- list(out$mean_filled)
      if (.progress) cli::cli_progress_update(id = pb)
    }
    if (.progress && length(run) > 0L) cli::cli_progress_done(id = pb)
  }

  # one write of every group's parts into one copy of the columns the result
  # holds. With no part, nothing was imputed or mean-filled, and those
  # columns are returned as they are. The parts cover disjoint columns, since
  # no feature is in two groups (prep_groups()), and columns no group
  # imputed come back unchanged.
  parts <- unlist(group_parts, recursive = FALSE)
  rm(group_parts)
  if (length(parts) > 0L) {
    obj <- fill_parts(obj, parts, keep = result_cols, n_miss = cmiss)
  } else if (!is.null(result_cols)) {
    obj <- obj[, result_cols, drop = FALSE]
  }
  rm(parts)

  # `post_imp` as the backend applied it, read off each group's result, not as
  # it was passed: the default here is NULL, which hands each group the
  # backend's own default, and `group$parameters` can set it group by group.
  # One value when the groups agree. NA when they do not, or when no group
  # reached the backend and none was passed, so that nothing was applied.
  post_imp_run <- unique(unlist(
    lapply(mean_filled, `[[`, "post_imp"),
    use.names = FALSE
  ))
  mean_filled <- combine_mean_filled(mean_filled, nrow(obj))
  obj <- new_slideimp_results(
    obj,
    imp_method,
    fallback = group$group[fallback_flags],
    post_imp = if (length(post_imp_run) == 1L) {
      post_imp_run
    } else if (length(post_imp_run) > 1L || is.null(post_imp)) {
      NA
    } else {
      post_imp
    },
    mean_filled_cols = mean_filled$cols,
    mean_filled_cells = mean_filled$cells
  )
  attr(obj, "metacaller") <- "group_imp"
  attr(obj, "fallback_action") <- on_infeasible
  obj
}
