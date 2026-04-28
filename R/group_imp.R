#' Resolve a group specification to a data.frame
#'
#' S3 generic for converting various group specifications into the
#' canonical data.frame form expected by [prep_groups()]. This generic
#' exists only to allow `slideimp.extra` to register the character method.
#'
#' @param x A group specification. The base package provides a method for
#' `data.frame`. `slideimp.extra` provides `character` for chip-name lookup.
#' @return A data.frame with at least a `feature` column, suitable for
#' passing to [prep_groups()].
#'
#' @examples
#' df <- data.frame(feature = c("cg1", "cg2"), group = c(1, 1))
#' slideimp_resolve_group(df)
#'
#' @keywords internal
#' @export
slideimp_resolve_group <- function(x) {
  UseMethod("slideimp_resolve_group")
}

#' @rdname slideimp_resolve_group
#' @export
slideimp_resolve_group.data.frame <- function(x) {
  unique(x)
}

#' @rdname slideimp_resolve_group
#' @export
slideimp_resolve_group.default <- function(x) {
  cli::cli_abort(c(
    "No {.fun slideimp_resolve_group} method for objects of class {.cls {class(x)[1]}}.",
    "i" = "Install and load {.pkg slideimp.extra} to enable character-based group lookup.",
    ">" = "See the {.pkg slideimp} README for installation instructions."
  ))
}

#' Remove Features from Auxiliary Columns
#'
#' Performs a fast, grouped set-difference using `collapse`. For each group
#' (element) in the lists, it removes any values in `aux` that are already
#' present in the corresponding `feat` element.
#'
#' @param aux A list of character vectors representing auxiliary column names
#' per group.
#' @param feat A list of character vectors representing feature column names
#' per group,
#'   must be the same length as `aux`.
#' @param iter Levels to be used with [collapse::gsplit()].
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

  keep <- is.na(collapse::fmatch(list(aux_gid, aux_flat), list(feat_gid, feat_flat)))
  gid_f <- factor(aux_gid[keep], levels = iter)
  out <- collapse::gsplit(aux_flat[keep], gid_f)
  names(out) <- NULL
  out
}

#' Prune List Elements to a Global Reference Set
#'
#' Helper function to filter elements within a list of vectors, keeping only
#' those present in a reference character vector `A`.
#'
#' @inheritParams group_imp
#'
#' @param lst A list of character vectors (e.g., group features or auxiliary
#' columns) to prune.
#' @param A Character vector to prune from.
#'
#' @returns A list of character vectors the same length as `lst`, keeping only
#'   elements that exist in `A`.
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

#' Compute Column Index Mappings for a Single Group
#'
#' Maps a group's feature and auxiliary column positions (relative to the
#' original matrix) into the structures needed by the imputation loop:
#' local feature indices within the submatrix, the combined column index
#' for submatrix extraction, and the feature names.
#'
#' @param g Integer scalar. The group index (position in `feat_splits` /
#'   `aux_splits` / `prep_groups`).
#' @param feat_splits A list of integer vectors, one per group, giving each
#'   group's feature column positions in the original matrix.
#' @param aux_splits A list of integer vectors (or `NULL`s), one per group,
#'   giving each group's auxiliary column positions in the original matrix.
#' @param prep_groups A list of character vectors, one per group, containing
#'   the feature names (i.e., `group$feature` from the prepped table).
#'
#' @returns A named list with three elements:
#'
#' * `features_idx_local`: Integer vector of feature positions within the
#'   extracted submatrix (`1:n_features`), used to slice imputed columns
#'   back out of the result.
#' * `col_idx`: Integer vector of all column positions (features then aux)
#'   in the original matrix, used to extract the submatrix via
#'   `obj[, col_idx]`.
#' * `features_names`: Character vector of feature names for this group.
#'
#' @noRd
#' @keywords internal
group_indices <- function(g, feat_splits, aux_splits, prep_groups) {
  f_idx <- feat_splits[[g]]
  col_idx <- c(f_idx, aux_splits[[g]])
  list(
    features_idx_local = seq_along(f_idx),
    col_idx = col_idx,
    features_names = prep_groups[[g]]
  )
}

#' Prepare Groups for Imputation
#'
#' Normalize and validate a grouping specification for use with [group_imp()].
#'
#' `prep_groups()` converts long-format or list-column group specifications
#' into a validated `slideimp_tbl`, enforces feature and auxiliary-column set
#' relationships, prunes dropped columns, and optionally pads small groups.
#'
#' @param obj_cn Character vector of column names from the data matrix, usually
#'   `colnames(obj)`.
#' @inheritParams group_imp
#' @param min_group_size Integer or `NULL`. Minimum total number of columns
#'   per group, counting both features and auxiliary columns. Groups smaller
#'   than this are padded with randomly sampled columns from `obj_cn`. If
#'   `NULL` or `0`, no padding is performed.
#' @param seed Integer, numeric, or `NULL`. Random seed used when padding small
#'   groups.
#'
#' @details
#' Let \eqn{A} be `obj_cn` and \eqn{B} be the union of all feature and
#' auxiliary names in `group`. When `allow_unmapped = FALSE`, the function
#' enforces \eqn{A \subseteq B}: every matrix column must appear somewhere in
#' the grouping specification.
#'
#' Elements in \eqn{B} but not in \eqn{A}, such as previously dropped probes,
#' are pruned from each group. Groups left with zero features after pruning or
#' subset demotion are removed with a diagnostic message.
#'
#' @returns A `data.frame` of class `slideimp_tbl` containing:
#' * `group`: original group labels, if provided, or sequential group labels.
#' * `feature`: a list-column of character vectors containing feature names.
#' * `aux`: a list-column of character vectors containing auxiliary names.
#' * `parameters`: a list-column of per-group configuration lists.
#'
#' @seealso [group_imp()]
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
    min.len = 1, unique = TRUE,
    any.missing = FALSE, .var.name = "obj_cn"
  )
  group <- slideimp_resolve_group(group)
  checkmate::assert_data_frame(group, min.rows = 1, .var.name = "group")
  checkmate::assert_names(
    colnames(group),
    must.include = "feature", .var.name = "group"
  )
  checkmate::assert_int(min_group_size, lower = 0, .var.name = "min_group_size")
  checkmate::assert_number(seed, null.ok = TRUE, .var.name = "seed")
  checkmate::assert_character(
    subset,
    min.len = 1, any.missing = FALSE, unique = TRUE,
    null.ok = TRUE, .var.name = "subset"
  )

  # Step 1: Normalize input format ----
  if ("group" %in% names(group) && inherits(group$feature, "character")) {
    if (anyNA(group$group)) {
      cli::cli_abort("{.code NA} is not allowed in {.code group$group}.")
    }

    checkmate::assert_character(
      group$feature,
      any.missing = FALSE, unique = TRUE,
      .var.name = "group$feature"
    )

    if (is.numeric(group$group) && !checkmate::test_integerish(group$group)) {
      cli::cli_warn(c(
        "Non-integer numeric {.code group$group} values detected.",
        "i" = "The {.field group} column is used to group the features in the {.field feature} column.",
        "i" = "Did you mean to use integer or character labels?"
      ))
    }

    group <- collapse::fsummarize(
      collapse::fgroup_by(group, group),
      feature = list(feature)
    )
  } else if (inherits(group$feature, "character")) {
    cli::cli_abort(c(
      "{.arg group} has a character {.field feature} column but no {.field group} column.",
      "i" = "Either add a {.field group} column to define groups, or use list-columns.",
      "i" = "See {.help prep_groups} for details."
    ))
  }

  checkmate::assert_list(
    group$feature,
    types = "character", min.len = 1,
    unique = TRUE, .var.name = "group$feature"
  )

  # Normalize aux — always present as a list of character vectors
  if ("aux" %in% names(group)) {
    checkmate::assert_list(
      group$aux,
      types = c("character", "logical", "null"),
      min.len = 1, .var.name = "group$aux"
    )
    group$aux <- lapply(group$aux, function(x) {
      if (is.null(x) || is.logical(x)) character(0) else unique(x[!is.na(x)])
    })
  } else {
    group$aux <- replicate(nrow(group), character(0), simplify = FALSE)
  }

  # Normalize parameters — always present as a list of lists
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
        "{length(not_in_B)} column{?s} in {.arg obj} {?has/have} no matching entry in {.arg group}: {fmt_trunc(not_in_B)}",
        "i" = "Run {.code setdiff(colnames(obj), group$feature)} to examine the unmapped probes.",
        "i" = "If you are sure {.arg group} is correct and want to leave these columns untouched, set {.code allow_unmapped = TRUE}."
      ))
    } else {
      cli::cli_inform(c(
        "!" = "{length(not_in_B)} column{?s} in {.arg obj} {?has/have} no matching entry in {.arg group} and will be left untouched.",
        "i" = "{fmt_trunc(not_in_B)}"
      ))
    }
  }

  # Features must not appear in more than one group
  if (anyDuplicated(all_feats) > 0) {
    dups <- unique(all_feats[duplicated(all_feats)])
    cli::cli_abort(c(
      "Features appear in more than one group:",
      "x" = "{fmt_trunc(dups)}"
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
    dropped <- which(empty)

    cli::cli_inform(c(
      "!" = "{length(dropped)} group{?s} dropped: no features remaining after matching {.arg obj_cn}.",
      "i" = "Dropped group indices: {paste(dropped, collapse = ', ')}"
    ))

    group <- group[!empty, , drop = FALSE]
    rownames(group) <- NULL
  }
  if (nrow(group) == 0) {
    cli::cli_abort(c(
      "No groups remain after pruning.",
      "i" = "Check that {.code group$feature} matches {.code colnames(obj)}."
    ))
  }

  # Step 3b: Apply subset, demote non-subset features to aux ----
  if (!is.null(subset)) {
    bad_cols <- setdiff(subset, A)
    if (length(bad_cols) > 0) {
      cli::cli_abort(c(
        "{cli::qty(length(bad_cols))}{.arg subset} contains column{?s} not in {.arg obj_cn}:",
        "x" = "{fmt_trunc(bad_cols)}"
      ))
    }
    all_feats_now <- unlist(group$feature)
    bad_feats <- setdiff(subset, all_feats_now)
    if (length(bad_feats) > 0) {
      cli::cli_abort(c(
        "{cli::qty(length(bad_feats))}{.arg subset} contains feature{?s} not assigned to any group:",
        "x" = "{fmt_trunc(bad_feats)}",
        "i" = "Add them to a group or remove them from {.arg subset}."
      ))
    }
    # Features in subset stay; the rest are demoted to aux
    for (i in seq_len(nrow(group))) {
      feat_g <- intersect(group$feature[[i]], subset)
      demoted <- setdiff(group$feature[[i]], feat_g)
      group$aux[[i]] <- c(group$aux[[i]], demoted)
      group$feature[[i]] <- feat_g
    }
    keep <- lengths(group$feature) > 0L
    if (!any(keep)) {
      cli::cli_abort("No groups have features to impute after applying {.arg subset}.")
    }
    if (any(!keep)) {
      dropped <- which(!keep)

      cli::cli_inform(c(
        "!" = "{length(dropped)} group{?s} dropped: no features remaining after applying {.arg subset}.",
        "i" = "Dropped group indices: {paste(dropped, collapse = ', ')}"
      ))

      group <- group[keep, , drop = FALSE]
      rownames(group) <- NULL
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
      group$aux <- Map(function(feat, aux, n) {
        if (n == 0L) {
          return(aux)
        }
        pool <- setdiff(A, c(feat, aux))
        if (length(pool) < n) {
          cli::cli_abort(c(
            "{.arg min_group_size} is too large.",
            "x" = "Not enough columns available to pad."
          ))
        }
        c(aux, sample(pool, size = n))
      }, group$feature, group$aux, need)
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
#' @param obj A numeric matrix with samples in rows and features in columns.
#' @param group Specification of how features should be grouped for imputation.
#'   Accepted formats are:
#'
#'   * A character scalar naming a supported Illumina platform; see Note.
#'   * A long-format `data.frame` with columns `group` and `feature`.
#'   * A list-column `data.frame` with a `feature` list-column. Optional
#'     list-columns are `aux`, for auxiliary feature names, and `parameters`,
#'     for group-specific parameter lists.
#'
#' @param subset Optional character vector of feature names to impute. If
#'   `NULL`, all grouped features are imputed. Features in a group but not in
#'   `subset` are demoted to auxiliary columns for that group. Groups left with
#'   zero features after demotion are dropped with a message.
#' @param allow_unmapped Logical. If `FALSE`, every column in `colnames(obj)`
#'   must appear in `group`. If `TRUE`, columns with no group assignment are
#'   left untouched and are not used as auxiliary columns.
#' @param k Integer or `NULL`. Number of nearest neighbors for K-NN imputation.
#'   If `NULL`, `k` may be supplied through `group$parameters`.
#' @param ncp Integer or `NULL`. Number of components for PCA imputation. If
#'   `NULL`, `ncp` may be supplied through `group$parameters`.
#' @param method Character or `NULL`. For K-NN imputation, one of
#'   `"euclidean"` or `"manhattan"`. For PCA imputation, one of
#'   `"regularized"` or `"EM"`. If `NULL`, the corresponding backend default
#'   is used unless overridden by `group$parameters`.
#' @param cores Integer. Number of cores for K-NN imputation only. For PCA
#'   imputation, use `mirai::daemons()` to parallelize across groups.
#' @param .progress Logical. If `TRUE`, show progress.
#' @param min_group_size Integer or `NULL`. Minimum total number of columns
#'   per group, counting both features and auxiliary columns. Groups smaller
#'   than this are padded with randomly sampled columns from `obj`.
#' @param post_imp Logical or `NULL`. If `TRUE`, replace any remaining missing
#'   values with column means after imputation. If `NULL`, the backend default
#'   is used unless overridden by `group$parameters`.
#' @param seed Integer, numeric, or `NULL`. Random seed for reproducibility.
#' @param pin_blas Logical. If `TRUE`, pin BLAS threads to 1 to reduce
#'   contention when using parallel PCA on systems linked with multithreaded
#'   BLAS.
#' @param on_infeasible Character. One of `"error"`, `"skip"`, or `"mean"`.
#'   Controls behavior when a group is infeasible for imputation, for example
#'   when `k` or `ncp` exceeds the number of usable columns after applying
#'   `colmax`.
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#'
#' @details
#' `group_imp()` performs K-NN or PCA imputation on feature groups
#' independently, which can substantially reduce run time for large matrices.
#'
#' Specify `k` and related arguments to use K-NN imputation, or `ncp` and
#' related arguments to use PCA imputation. If both `k` and `ncp` are `NULL`,
#' `group$parameters` must supply either `k` or `ncp` for every group.
#'
#' Group-specific parameters in `group$parameters` take priority over global
#' arguments. Global arguments fill in any missing values. All groups must use
#' the same imputation method.
#'
#' For method-specific arguments inherited from [knn_imp()] or [pca_imp()],
#' `NULL` means the backend default is used unless overridden by
#' `group$parameters`.
#'
#' Per-group `k` is capped at the number of usable columns in the group minus
#' one. Per-group `ncp` is capped at the maximum feasible number of PCA
#' components for that group's submatrix. A warning is issued when capping
#' occurs.
#'
#' @section Parallelization:
#' * K-NN: use the `cores` argument. If `mirai::daemons()` are active, `cores`
#'   is automatically set to `1` to avoid nested parallelism.
#' * PCA: use `mirai::daemons()` instead of `cores`.
#'
#' When running PCA imputation in parallel with `mirai`, set `pin_blas = TRUE`
#' in [tune_imp()] or [group_imp()] to prevent BLAS threads from
#' oversubscribing CPU cores. This relies on `RhpcBLASctl` and works with
#' OpenBLAS and MKL (typical on Linux, and on Windows after an OpenBLAS swap).
#' `pin_blas = TRUE` may have no effect on macOS.
#'
#' @inheritSection pca_imp Performance tips
#'
#' @note
#' A character scalar can be passed to `group` to name a supported Illumina
#' platform, such as `"EPICv2"` or `"EPICv2_deduped"`. This requires the
#' optional `slideimp.extra` package to be installed. Supported platforms are
#' listed in `slideimp.extra::slideimp_arrays`.
#'
#' @returns A numeric matrix of the same dimensions as `obj`, with missing
#' values imputed. The returned object has class `slideimp_results`.
#'
#' @seealso [prep_groups()], [knn_imp()], [pca_imp()]
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
  tree = NULL,
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
  pin_blas = FALSE,
  na_check = TRUE,
  on_infeasible = c("error", "skip", "mean")
) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric", col.names = "unique",
    null.ok = FALSE, .var.name = "obj"
  )
  check_finite(obj)
  checkmate::assert_flag(pin_blas, null.ok = FALSE, .var.name = "pin_blas")
  on_infeasible <- match.arg(on_infeasible)
  if (!is.null(solver)) {
    solver <- match.arg(solver, c("auto", "dsyevr", "lobpcg"))
  }
  cn <- colnames(obj)
  rn <- rownames(obj)
  attributes(obj) <- list(dim = dim(obj))

  # Step 1: Build canonical groups via prep_groups() ----
  # After this call, group$feature, group$aux, and group$parameters are all
  # guaranteed to exist as properly typed list-columns.
  # subset is also handled here — non-subset features are demoted to aux.
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

  # Step 2: Resolve parameters ----
  if (!is.null(k) && !is.null(ncp)) {
    cli::cli_abort("Cannot specify both {.arg k} and {.arg ncp} as global parameters.")
  }

  # Global fills gaps (group-wise wins)
  global_params <- list(
    k = k, method = method, colmax = colmax, post_imp = post_imp,
    dist_pow = dist_pow, tree = tree, ncp = ncp, scale = scale,
    coeff.ridge = coeff.ridge, threshold = threshold, row.w = row.w,
    seed = seed, nb.init = nb.init, maxiter = maxiter, miniter = miniter,
    lobpcg_control = lobpcg_control, solver = solver
  )
  global_params <- global_params[!vapply(global_params, is.null, logical(1))]

  group$parameters <- lapply(group$parameters, function(p) {
    for (nm in names(global_params)) {
      if (is.null(p[[nm]])) p[[nm]] <- global_params[[nm]]
    }
    p
  })

  # Validate: each group has exactly one of k or ncp
  has_k <- vapply(group$parameters, \(p) "k" %in% names(p), logical(1))
  has_ncp <- vapply(group$parameters, \(p) "ncp" %in% names(p), logical(1))

  if (any(has_k & has_ncp)) {
    bad <- which(has_k & has_ncp)
    cli::cli_abort(c(
      "{cli::qty(length(bad))}Group{?s} {fmt_trunc(bad)} {cli::qty(length(bad))}{?has/have} both {.arg k} and {.arg ncp} in parameters.",
      "i" = "Specify only one imputation method per group."
    ))
  }
  if (any(!has_k & !has_ncp)) {
    bad <- which(!has_k & !has_ncp)
    cli::cli_abort(c(
      "{cli::qty(length(bad))}Group{?s} {fmt_trunc(bad)} {cli::qty(length(bad))}{?has/have} neither {.arg k} nor {.arg ncp}.",
      "i" = "Specify global {.arg k}/{.arg ncp} or set them in {.code group$parameters}."
    ))
  }
  if (any(has_k) && any(has_ncp)) {
    cli::cli_abort(c(
      "Inconsistent imputation methods across groups.",
      "i" = "All groups must use either {.arg k} (K-NN) or {.arg ncp} (PCA)."
    ))
  }

  imp_method <- if (all(has_k)) "knn" else "pca"
  is_knn_mode <- imp_method == "knn"
  reg <- .param_registry[[imp_method]]

  # Validate method values
  valid_methods <- reg$methods
  bad_method <- vapply(group$parameters, function(p) {
    !is.null(p$method) && !(p$method %in% valid_methods)
  }, logical(1))

  if (any(bad_method)) {
    bad <- which(bad_method)
    cli::cli_abort(c(
      "{cli::qty(length(bad))}Invalid {.arg method} for {toupper(imp_method)} in group{?s} {fmt_trunc(bad)}.",
      "i" = "Must be one of: {.val {valid_methods}}."
    ))
  }

  # Validate parameter names
  all_param_names <- unique(unlist(lapply(group$parameters, names), use.names = FALSE))

  check_unknown_params(
    param_names = all_param_names,
    mode = imp_method,
    arg = "group$parameters"
  )

  # Cap per-group k/ncp
  group_size <- feat_lengths + aux_lengths
  required_param <- reg$required

  for (i in seq_len(nrow(group))) {
    p <- group$parameters[[i]]
    cap <- if (is_knn_mode) {
      group_size[i] - 1L
    } else {
      min(group_size[i] - 1L, nrow(obj) - 2L)
    }
    if (p[[required_param]] > cap) {
      cli::cli_warn(
        "Group {i}: {.arg {required_param}} capped from {p[[required_param]]} to {cap} (group size = {group_size[i]})."
      )
      p[[required_param]] <- cap
    }
    if (p[[required_param]] < 1L) {
      cli::cli_abort(c(
        "Group {i}: {.arg {required_param}} must be {.code >= 1} after capping.",
        "x" = "Group size = {group_size[i]}."
      ))
    }
    group$parameters[[i]] <- p
  }

  cli::cli_inform(
    "Imputing {nrow(group)} group{?s} using {.strong {toupper(imp_method)}}."
  )

  # Step 3: Imputation loop ----
  # Column-index lookups
  iter <- seq_len(nrow(group))

  all_feats_vec <- unlist(group$feature)
  all_feats_pos <- collapse::fmatch(all_feats_vec, cn)
  gid_feat <- factor(rep.int(iter, feat_lengths), levels = iter)
  feat_splits <- collapse::gsplit(all_feats_pos, gid_feat)
  names(feat_splits) <- NULL

  all_aux_pos <- collapse::fmatch(unlist(group$aux), cn)
  gid_aux <- factor(rep.int(iter, aux_lengths), levels = iter)
  aux_splits <- collapse::gsplit(all_aux_pos, gid_aux)
  names(aux_splits) <- NULL

  indices <- lapply(
    iter, group_indices,
    feat_splits = feat_splits,
    aux_splits = aux_splits,
    prep_groups = group$feature
  )

  # Parallelism resolution
  parallelize <- tryCatch(mirai::require_daemons(), error = function(e) FALSE)

  if (is_knn_mode) {
    if (cores > 1 && parallelize) {
      cli::cli_inform(c(
        "!" = "Both {.arg cores} > 1 and active {.pkg mirai} daemons detected.",
        "i" = "Setting {.arg cores} = 1 to avoid nested parallelism. Distribution will be handled by {.pkg mirai}."
      ))
      cores <- 1
    }
  } else if (cores > 1 && !parallelize) {
    cli::cli_inform(c(
      "!" = "{.arg cores} = {cores} is ignored for PCA imputation.",
      "i" = "Parallelism comes from threaded BLAS or {.pkg mirai} daemons across groups."
    ))
    cores <- 1
  }

  # Build per-group call parameters
  params <- lapply(iter, function(i) {
    p <- group$parameters[[i]]
    p$na_check <- FALSE
    if (is_knn_mode) {
      p$cores <- cores
      p$subset <- indices[[i]]$features_idx_local
      p$.progress <- FALSE
    }
    p
  })

  imp_fn <- if (is_knn_mode) knn_imp else pca_imp
  group_ids <- group$group

  if (parallelize) {
    cli::cli_inform("Running mode: {.strong mirai}")
  } else if (is_knn_mode && cores > 1) {
    cli::cli_inform("Running mode: {.strong threaded} ({cores} cores)")
  } else {
    cli::cli_inform("Running mode: {.strong sequential}")
  }

  # Imputation
  if (parallelize) {
    feat_cumsum <- cumsum(c(0L, feat_lengths))
    out_ranges <- lapply(iter, function(i) {
      (feat_cumsum[i] + 1L):feat_cumsum[i + 1L]
    })

    check_pin_blas(pin_blas)

    big_obj <- bigmemory::as.big.matrix(obj, shared = TRUE)
    big_obj_desc <- bigmemory::describe(big_obj)
    big_out <- bigmemory::big.matrix(
      nrow = nrow(obj), ncol = length(all_feats_pos),
      type = "double", shared = TRUE
    )
    big_out_desc <- bigmemory::describe(big_out)
    on.exit(
      {
        rm(big_obj, big_out)
        gc()
      },
      add = TRUE
    )

    crated_fn <- carrier::crate(
      function(i) {
        if (pin_blas) {
          RhpcBLASctl::blas_set_num_threads(1)
          RhpcBLASctl::omp_set_num_threads(1)
        }

        src <- bigmemory::attach.big.matrix(big_obj_desc)
        dst <- bigmemory::attach.big.matrix(big_out_desc)
        sub_mat <- src[, indices[[i]]$col_idx, drop = FALSE]

        imputed <- tryCatch(
          suppressMessages(do.call(imp_fn, c(list(obj = sub_mat), params[[i]]))),
          slideimp_infeasible = function(e) {
            switch(on_infeasible,
              error = cli::cli_abort(c(
                "Group {group_ids[[i]]} is infeasible for imputation.",
                "x" = conditionMessage(e),
                "i" = "Set {.code on_infeasible = 'skip'} or {.code on_infeasible = 'mean'} to continue."
              ), parent = e),
              skip = structure(sub_mat, fallback = TRUE),
              mean = structure(
                mean_imp_col(sub_mat, subset = indices[[i]]$features_idx_local),
                fallback = TRUE
              )
            )
          }
        )

        dst[, out_ranges[[i]]] <- imputed[, indices[[i]]$features_idx_local, drop = FALSE]
        return(isTRUE(attr(imputed, "fallback")))
      },
      big_obj_desc = big_obj_desc,
      big_out_desc = big_out_desc,
      imp_fn = imp_fn,
      indices = indices,
      params = params,
      pin_blas = pin_blas,
      out_ranges = out_ranges,
      on_infeasible = on_infeasible,
      group_ids = group_ids
    )
    m <- mirai::mirai_map(iter, crated_fn)
    fallback_flags <- unlist(m[.progress = .progress])
    obj[, all_feats_pos] <- big_out[, ]
  } else {
    if (.progress) pb <- cli::cli_progress_bar(name = "Imputing", total = length(iter))
    fallback_flags <- logical(length(iter))
    for (i in iter) {
      sub_mat <- obj[, indices[[i]]$col_idx, drop = FALSE]

      imputed <- tryCatch(
        suppressMessages(do.call(imp_fn, c(list(obj = sub_mat), params[[i]]))),
        slideimp_infeasible = function(e) {
          switch(on_infeasible,
            error = cli::cli_abort(c(
              "Group {group_ids[[i]]} is infeasible for imputation.",
              "x" = conditionMessage(e),
              "i" = "Set {.code on_infeasible = 'skip'} or {.code on_infeasible = 'mean'} to continue."
            ), parent = e),
            skip = structure(sub_mat, fallback = TRUE),
            mean = structure(
              mean_imp_col(sub_mat, subset = indices[[i]]$features_idx_local),
              fallback = TRUE
            )
          )
        }
      )

      obj[, feat_splits[[i]]] <- imputed[, indices[[i]]$features_idx_local, drop = FALSE]
      fallback_flags[i] <- isTRUE(attr(imputed, "fallback"))
      if (.progress) cli::cli_progress_update(id = pb)
    }
    if (.progress) cli::cli_progress_done(id = pb)
  }

  fallback_groups <- group$group[fallback_flags]
  has_remaining_na <- if (na_check) anyNA(obj[, all_feats_pos]) else NULL

  colnames(obj) <- cn
  rownames(obj) <- rn
  class(obj) <- c("slideimp_results", class(obj))
  attr(obj, "imp_method") <- imp_method
  attr(obj, "metacaller") <- "group_imp"
  attr(obj, "fallback") <- fallback_groups
  attr(obj, "fallback_action") <- on_infeasible
  attr(obj, "has_remaining_na") <- has_remaining_na
  attr(obj, "post_imp") <- post_imp
  obj
}
