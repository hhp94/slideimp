#' Remove Features from Auxiliary Columns
#'
#' Performs a fast, grouped set-difference using `collapse`. For each group
#' (element) in the lists, it removes any values in `aux` that are already
#' present in the corresponding `feat` element.
#'
#' @param aux A list of character vectors representing auxiliary column names per group.
#' @param feat A list of character vectors representing feature column names per group,
#'   must be the same length as `aux`.
#' @param iter Levels to be used with [collapse::gsplit()]
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
#' @param lst A list of character vectors (e.g., group features or auxiliary columns) to prune.
#' @param A Character vector to prune from.
#' @param iter Levels to be used with [collapse::gsplit()]
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

#' Get Cleaned Illumina Manifests
#'
#' Downloads pre-cleaned Illumina manifests via the `slideimp.extra`
#' package for use with [group_imp()].
#'
#' @param group Character. Name of the requested manifest.
#'
#' @returns A data.frame containing the cleaned Illumina manifest for
#' the requested array type.
#'
#' @seealso [group_imp()]
#'
#' @noRd
#' @keywords internal
slideimp_extra_manifests <- function(group = NULL) {
  if (!requireNamespace("slideimp.extra", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg slideimp.extra} package is required for this functionality.",
      "i" = "Install it with:",
      ">" = "See the {.pkg slideimp} README or website for installation instructions."
    ))
  }
  checkmate::assert_choice(group, choices = slideimp.extra::slideimp_arrays)
  deduped <- group %in% c("EPICv2_deduped", "MSA_deduped")
  if (deduped) {
    group <- switch(group,
      EPICv2_deduped = "EPICv2",
      MSA_deduped = "MSA"
    )
  }
  return(slideimp.extra::ilmn_manifest(chip = group, deduped = deduped))
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
#' Normalizes and validates a grouping specification for use with [group_imp()].
#' Converts long-format or canonical list-column input into a validated
#' `slideimp_tbl`, enforcing set relationships, pruning dropped columns,
#' and optionally padding small groups.
#'
#' @inheritParams group_imp
#'
#' @param obj_cn Character vector of column names from the data matrix
#' (e.g., `colnames(obj)`). Every element must appear in `group$feature` unless `allowed_unmapped = TRUE`.
#'
#' @details
#' ### Set Validation
#' Let \eqn{A} = `obj_cn` and \eqn{B} = the union of all feature and auxiliary names
#' in `group`. The function enforces \eqn{A \subseteq B}: every column in
#' the matrix must appear somewhere in the manifest.
#'
#' * `Pruning:` Elements in \eqn{B} but not in \eqn{A} (e.g., QC-dropped probes)
#'   are silently pruned from each group.
#' * `Dropping:` Groups left with zero features after pruning are
#'   removed entirely with a diagnostic message.
#'
#' @returns A `data.frame` of class `slideimp_tbl` containing:
#' * `group`: Original group labels (if provided).
#' * `feature`: A list-column of character vectors (feature names).
#' * `aux`: A list-column of character vectors (auxiliary names).
#' * `parameters`: A list-column of per-group configuration lists.
#'
#' @seealso [group_imp()]
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
  if (is.character(group)) {
    group <- slideimp_extra_manifests(group)
  } else {
    checkmate::assert_data_frame(group, min.rows = 1, .var.name = "group")
    group <- unique(group)
  }
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
    stopifnot("`NA` is not allowed in `group$group`" = !anyNA(group$group))
    checkmate::assert_character(
      group$feature,
      any.missing = FALSE, unique = TRUE,
      .var.name = "group$feature"
    )
    if (is.numeric(group$group) && !checkmate::test_integerish(group$group)) {
      warning(
        "Non-integer numeric 'group$group' values detected. ",
        "The `group` column is used to group the features in the `feature` ",
        "column. Did you mean to use integer or character labels?"
      )
    }
    group <- collapse::fsummarize(
      collapse::fgroup_by(group, group),
      feature = list(feature)
    )
    # group$group <- NULL
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
      message(paste0(
        length(not_in_B), " column(s) in `obj` have no matching entry in `group` ",
        "and will be left untouched: ",
        fmt_trunc(not_in_B), "."
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
    message(
      "Groups ", paste(which(empty), collapse = ", "),
      " dropped: no features remaining after matching obj columns."
    )
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
      message(
        "Groups ", paste(which(!keep), collapse = ", "),
        " dropped: no features remaining after applying subset."
      )
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
  class(group) <- c("slideimp_tbl", "data.frame")
  return(group)
}

#' Grouped K-NN or PCA Imputation
#'
#' Perform K-NN or PCA imputation independently on feature groups
#' (e.g. by chromosomes, flanking probes, or clustering-based groups).
#'
#' @inheritParams slide_imp
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#' @param group Specification of how features should be grouped for imputation.
#' Accepts three formats:
#'
#' * `character`: (Requires the `{slideimp.extra}` package, available on GitHub). A single string
#'   naming a supported Illumina platform (e.g., `"EPICv2"`, `"EPICv2_deduped"`).
#'   The manifest is fetched automatically. For installation instructions, see the
#'   package README. Supported platforms can be viewed via
#'   `slideimp.extra::slideimp_arrays`.
#' * `data.frame` (Long format):
#'     * `group`: Column identifying the group for each feature.
#'     * `feature`: Character column of individual feature names.
#' * `data.frame` (List-column format):
#'     * `feature`: List-column of character vectors to impute (e.g., from [prep_groups()]).
#'     * `aux`: (Optional) List-column of auxiliary names used for context.
#'     * `parameters`: (Optional) List-column of group-specific parameter lists.
#'
#' @param subset Character vector of feature names to impute (default `NULL`
#' means impute all features). Must be a subset of `obj_cn` (`colnames(obj)`)
#' and must appear in at least one group's `feature`. Features in a group
#' but not in `subset` are demoted to auxiliary columns for that group. Groups left
#' with zero features after demotion are dropped with a message.
#'
#' @param allow_unmapped Logical. If `FALSE`, every column in
#' `colnames(obj)` *must* appear in `group`. If `TRUE`, columns that have no
#' group assignment are left untouched (neither imputed nor used as auxiliary
#' columns) and a message is issued instead of an error.
#'
#' @param min_group_size Integer or `NULL`. Minimum number of columns
#' (features + aux) per group. Groups smaller than this are padded with
#' randomly sampled columns from `obj`. Passed to [prep_groups()] internally.
#'
#' @param cores The number of OpenMP cores for K-NN imputation **only**. For PCA
#' or mirai-based parallelism, use `mirai::daemons()` instead.
#'
#' @param .progress Show imputation progress (default = `TRUE`).
#' @param seed Numeric or `NULL`. Random seed for reproducibility when padding
#' for `min_group_size` and passed to [pca_imp()].
#'
#' @param pin_blas Logical. If `TRUE`, pin BLAS thread to 1 to help
#' with parallel performance on systems linked with multi-threaded BLAS.
#'
#' @details
#' This function performs K-NN or PCA imputation on groups of features
#' independently, which significantly reduces imputation time for large
#' datasets.
#'
#' Specify `k` and related arguments to use K-NN, `ncp` and related arguments
#' for PCA imputation. If `k` and `ncp` are both `NULL`, then
#' `group$parameters` must contain either `k` or `ncp` for every group.
#'
#' ## Parameter resolution
#' Group-wise parameters (in `group$parameters`) take priority. Global
#' arguments (`k`, `ncp`, `method`, etc.) fill in any gaps where a group has
#' no value set. All groups must agree on the imputation method (all KNN or
#' all PCA). Per-group `k` is capped at `group_size - 1` and `ncp` at
#' `min(nrow(group) - 2L, ncol(group) - 1L)`, with a warning when capping occurs.
#'
#' ## Strategies for grouping
#' - Breaking down search space by chromosomes
#' - Grouping features with their flanking values/neighbors
#' - Using clusters identified by column clustering techniques
#'
#' @section Parallelization:
#' Parallelization behavior depends on the imputation method:
#'
#' - **KNN**: use the `cores` argument (if OpenMP is available).
#'   If `mirai::daemons()` are also active, `cores` is automatically
#'   set to 1 to avoid nested parallelism.
#' - **PCA**: use `mirai::daemons()` instead of `cores`.
#'
#' **Linux / OpenBLAS / MKL users:** If your machine uses a multi-threaded
#' BLAS (e.g., OpenBLAS or Intel MKL), set `pin_blas = TRUE` when tuning
#' PCA imputation in parallel. Without it, BLAS threads and `mirai` workers
#' compete for cores, which can cause slowdowns (CPU thrashing).
#'
#' **macOS users:** OpenMP is typically unavailable on macOS unless manually
#' configured. `cores` will fall back to 1 automatically; use
#' `mirai::daemons()` for parallelization instead.
#'
#' @inherit knn_imp return
#'
#' @export
#'
#' @seealso [prep_groups()]
#'
#' @examples
#' # Generate example data with missing values
#' set.seed(1234)
#' to_test <- sim_mat(50, 20, perc_total_na = 0.3, perc_col_na = 1)
#' obj <- to_test$input
#' group <- to_test$col_group # metadata that maps `colnames(obj)` to groups
#' head(group)
#'
#' # Simple grouped K-NN imputation
#' results <- group_imp(obj, group = group, k = 2)
#'
#' # Impute only a subset of features
#' subset_features <- sample(to_test$col_group$feature, size = 10)
#' knn_subset <- group_imp(obj, group = group, subset = subset_features, k = 5)
#'
#' # Use prep_groups() to inspect and tweak per-group parameters
#' prepped <- prep_groups(colnames(obj), group)
#' prepped$parameters <- lapply(seq_len(nrow(prepped)), \(i) list(k = 5))
#' prepped$parameters[[2]]$k <- 10
#' knn_grouped <- group_imp(obj, group = prepped, cores = 2)
#'
#' # PCA imputation with mirai parallelism
#' mirai::daemons(2)
#' pca_grouped <- group_imp(obj, group = group, ncp = 2)
#' mirai::daemons(0)
#' pca_grouped
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
  max_cache = NULL,
  # PCA arguments
  scale = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  row.w = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL,
  pin_blas = FALSE,
  na_check = TRUE
) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric", col.names = "unique",
    null.ok = FALSE, .var.name = "obj"
  )
  checkmate::assert_flag(pin_blas, null.ok = FALSE, .var.name = "pin_blas")
  cn <- colnames(obj)
  rn <- rownames(obj)
  # obj_attrs <- attributes(obj)
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
    dist_pow = dist_pow, tree = tree, max_cache = max_cache,
    ncp = ncp, scale = scale, coeff.ridge = coeff.ridge,
    threshold = threshold, row.w = row.w, seed = seed,
    nb.init = nb.init, maxiter = maxiter, miniter = miniter
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
      "i" = "All groups must use either {.arg k} (KNN) or {.arg ncp} (PCA)."
    ))
  }

  imp_method <- if (all(has_k)) "knn" else "pca"
  is_knn_mode <- imp_method == "knn"

  # Validate method values
  valid_methods <- if (is_knn_mode) c("euclidean", "manhattan") else c("regularized", "EM")
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
  allowed_params <- if (is_knn_mode) {
    c("k", "method", "colmax", "post_imp", "dist_pow", "tree", "max_cache")
  } else {
    c(
      "ncp", "scale", "method", "coeff.ridge", "row.w",
      "threshold", "seed", "nb.init", "maxiter", "miniter",
      "colmax", "post_imp"
    )
  }
  all_param_names <- unique(unlist(lapply(group$parameters, names)))
  unknown_params <- setdiff(all_param_names, allowed_params)
  if (length(unknown_params) > 0) {
    cli::cli_abort(c(
      "{cli::qty(length(unknown_params))}Unknown parameter{?s} for {imp_method} method:",
      "x" = "{fmt_trunc(unknown_params, 10)}"
    ))
  }

  # Cap per-group k/ncp
  group_size <- feat_lengths + aux_lengths
  required_param <- if (is_knn_mode) "k" else "ncp"

  for (i in seq_len(nrow(group))) {
    p <- group$parameters[[i]]
    cap <- if (is_knn_mode) {
      group_size[i] - 1L
    } else {
      min(group_size[i] - 1L, nrow(obj) - 2L)
    }
    if (p[[required_param]] > cap) {
      warning(sprintf(
        "Group %d: %s capped from %d to %d (group size = %d).",
        i, required_param, p[[required_param]], cap, group_size[i]
      ))
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

  message(sprintf("Imputing %d group(s) using %s.", nrow(group), toupper(imp_method)))

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
    if (cores > 1) {
      if (!has_openmp()) {
        message("OpenMP not available (common on macOS). KNN will run single-threaded. Use mirai::daemons() for parallelization.")
        cores <- 1
      } else if (parallelize) {
        message(
          "Both `cores > 1` and `mirai::daemons()` detected. ",
          "Setting `cores = 1` to avoid nested parallelism. ",
          "Parallelization will be handled by `mirai`."
        )
        cores <- 1
      }
    }
  } else if (cores > 1) {
    warning(
      "`cores` is ignored for PCA imputation; parallelism comes from ",
      "threaded BLAS or mirai daemons across groups. Setting `cores = 1`."
    )
    cores <- 1
  }

  # Build per-group call parameters
  params <- lapply(iter, function(i) {
    p <- group$parameters[[i]]
    p$na_check <- FALSE
    if (is_knn_mode) {
      p$cores <- cores
      p$subset <- indices[[i]]$features_idx_local
    }
    p
  })

  imp_fn <- if (is_knn_mode) knn_imp else pca_imp

  if (parallelize) {
    message("Running Mode: parallel (mirai across groups)...")
  } else if (is_knn_mode && cores > 1) {
    message("Running Mode: parallel (OpenMP within groups)...")
  } else {
    message("Running Mode: sequential ...")
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
        imputed <- do.call(imp_fn, c(list(obj = sub_mat), params[[i]]))
        dst[, out_ranges[[i]]] <- imputed[, indices[[i]]$features_idx_local, drop = FALSE]
        return(isTRUE(attr(imputed, "fallback")))
      },
      big_obj_desc = big_obj_desc,
      big_out_desc = big_out_desc,
      imp_fn = imp_fn,
      indices = indices,
      params = params,
      pin_blas = pin_blas,
      out_ranges = out_ranges
    )
    m <- mirai::mirai_map(iter, crated_fn)
    fallback_flags <- unlist(m[.progress = .progress])
    obj[, all_feats_pos] <- big_out[, ]
  } else {
    if (.progress) pb <- cli::cli_progress_bar(total = length(iter))
    fallback_flags <- logical(length(iter))
    for (i in iter) {
      sub_mat <- obj[, indices[[i]]$col_idx, drop = FALSE]
      imputed <- do.call(imp_fn, c(list(obj = sub_mat), params[[i]]))
      obj[, feat_splits[[i]]] <- imputed[, indices[[i]]$features_idx_local, drop = FALSE]
      if (.progress) cli::cli_progress_update(id = pb)
      fallback_flags[i] <- isTRUE(attr(imputed, "fallback"))
    }
    if (.progress) cli::cli_progress_done(id = pb)
  }

  fallback_groups <- which(fallback_flags)
  has_remaining_na <- if (na_check) anyNA(obj[, all_feats_pos]) else NULL

  colnames(obj) <- cn
  rownames(obj) <- rn
  class(obj) <- c("slideimp_results", class(obj))
  attr(obj, "imp_method") <- imp_method
  attr(obj, "metacaller") <- "group_imp"
  attr(obj, "fallback") <- fallback_groups
  attr(obj, "has_remaining_na") <- has_remaining_na
  obj
}
