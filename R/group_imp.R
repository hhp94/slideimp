fmt_trunc <- function(x, n = 5) {
  paste(paste(x[seq_len(min(n, length(x)))], collapse = ", "), "... (see ?prep_groups)")
}

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
#' (i.e., `colnames(obj)`). Every element must appear in the `group$feature`.
#'
#' @details
#'
#' **Set validation.** Let A = `obj_cn` and B = the union of all feature and
#' aux names in `group$feature`. The function enforces \eqn{A \subseteq B}: every
#' column in the matrix must appear somewhere in the manifest (`group`). Elements
#' in B but not in A (QC-dropped probes) are silently pruned from each group.
#' Groups left with zero features after pruning are dropped with a message.
#'
#' @return A `data.frame` of class `slideimp_tbl` with columns:
#' - **feature**: A list-column of character vectors of feature names
#' - **aux**: A list-column of auxiliary column names
#' - **parameters**: A list-column of parameter lists
#'
#' @seealso [group_imp()]
#'
#' @export
#'
#' @examples
#' # Simulate some data, 10 rows, 100 columns.
#' sim_obj <- sim_mat(n = 10, p = 100, perc_col_na = 1)
#' obj <- sim_obj$input
#' obj_meta <- sim_obj$col_group
#'
#' # Inspect grouping structure before imputation
#' prepped <- prep_groups(colnames(obj), obj_meta)
#' prepped
#'
#' # Pad small groups to a minimum size (say 60)
#' prepped_padded <- prep_groups(colnames(obj), obj_meta, min_group_size = 60)
#' prepped_padded
#'
#' # Impute only a subset of features (others become aux)
#' sub <- sample(colnames(obj), 50)
#' prepped_sub <- prep_groups(colnames(obj), obj_meta, subset = sub)
#' prepped_sub
#'
#' # Tweak per-group parameters, then pass to group_imp
#' prepped$parameters <- lapply(seq_len(nrow(prepped)), \(i) list(k = 5))
#' prepped$parameters[[2]]$k <- 10
#' prepped$parameters
#' imputed <- group_imp(obj, prepped)
#' imputed
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
  checkmate::assert_data_frame(group, min.rows = 1, .var.name = "group")
  group <- unique(group)
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
    stop(
      "`group` has a character 'feature' column but no 'group' column. ",
      "Either add a 'group' column to define groups, or use list-columns ",
      "(see ?prep_groups)."
    )
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
      stop(
        "These column(s) in `obj` have no group assignment: ",
        fmt_trunc(not_in_B),
        "\nIf you are sure `group` is correct ",
        "and want to leave these columns untouched, set ",
        "`allow_unmapped = TRUE`."
      )
    } else {
      warning(
        length(not_in_B),
        " column(s) in `obj` have no group assignment ",
        "and will be left completely untouched: ",
        fmt_trunc(not_in_B),
        "."
      )
    }
  }

  # Features must not appear in more than one group
  if (anyDuplicated(all_feats) > 0) {
    dups <- unique(all_feats[duplicated(all_feats)])
    stop(
      "Features appear in more than one group: ",
      fmt_trunc(dups)
    )
  }

  # Prune B \ A (QC-dropped probes not in obj)
  group$feature <- lapply(group$feature, intersect, A)
  group$aux <- lapply(group$aux, intersect, A)
  group$aux <- Map(setdiff, group$aux, group$feature)

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
    stop("No groups remain after pruning. Check that group features match obj column names.")
  }

  # Step 3b: Apply subset, demote non-subset features to aux ----
  if (!is.null(subset)) {
    bad_cols <- setdiff(subset, A)
    if (length(bad_cols) > 0) {
      stop(
        "subset contains columns not in obj_cn: ",
        fmt_trunc(bad_cols)
      )
    }
    all_feats_now <- unlist(group$feature)
    bad_feats <- setdiff(subset, all_feats_now)
    if (length(bad_feats) > 0) {
      stop(
        "subset contains features not assigned to any group: ",
        fmt_trunc(bad_feats),
        ". Add them to a group or remove from subset."
      )
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
      stop("No groups have features to impute after applying subset.")
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
          stop("`min_group_size` is too large; not enough columns to pad.")
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
#' K-NN or PCA imputation by groups, such as chromosomes, flanking columns,
#' or clusters identified by column clustering techniques.
#'
#' @inheritParams slide_imp
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#' @param group A data.frame describing how features should be grouped for
#' imputation. Accepts two formats:
#'
#'  **Long format**:
#'    - `group`: A column identifying which group each feature belongs to
#'    - `feature`: A character column of individual feature names
#'
#'  **List-column format** (e.g., output of [prep_groups()]):
#'    - `feature`: A list-column of character vectors of feature column names
#'    to impute
#'    - `aux`: (Optional) A list-column of character vectors of auxiliary
#'    column names used for imputation but not imputed themselves
#'    - `parameters`: (Optional) A list-column of group-specific parameter
#'    lists. Group-level values take priority; global arguments fill gaps.
#'
#' @param subset Character vector of feature names to impute (default `NULL`
#' means impute all features). Must be a subset of `obj_cn` and must appear
#' in at least one group's `feature` list. Features in a group but not in
#' `subset` are demoted to auxiliary columns for that group. Groups left
#' with zero features after demotion are dropped with a message.
#'
#' @param allow_unmapped Logical. If `FALSE` (default), every column in
#' `colnames(obj)` *must* appear in at least one group's `feature` or `aux`.
#' If `TRUE`, columns that have no group assignment are left completely
#' untouched (neither imputed nor used as auxiliary columns) and a warning
#' is issued instead of an error.
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
#' @inherit knn_imp note return
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
#' group <- to_test$col_group
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
  cores = 1,
  .progress = TRUE,
  min_group_size = NULL,
  colmax = NULL,
  post_imp = NULL,
  dist_pow = NULL,
  tree = NULL,
  max_cache = NULL,
  scale = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  row.w = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL
) {
  checkmate::assert_matrix(
    obj,
    mode = "numeric", col.names = "unique",
    null.ok = FALSE, .var.name = "obj"
  )
  cn <- colnames(obj)

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

  # Step 2: Resolve parameters ----
  if (!is.null(k) && !is.null(ncp)) {
    stop("Cannot specify both 'k' and 'ncp' as global parameters.")
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
    stop(sprintf(
      "Group(s) %s have both 'k' and 'ncp' in parameters.",
      paste(which(has_k & has_ncp), collapse = ", ")
    ))
  }
  if (any(!has_k & !has_ncp)) {
    stop(sprintf(
      "Group(s) %s have neither 'k' nor 'ncp'. Specify global 'k'/'ncp' or set them in group$parameters.",
      paste(which(!has_k & !has_ncp), collapse = ", ")
    ))
  }
  if (any(has_k) && any(has_ncp)) {
    stop("Inconsistent imputation methods across groups; all groups must use either k (KNN) or ncp (PCA).")
  }

  imp_method <- if (all(has_k)) "knn" else "pca"
  is_knn_mode <- imp_method == "knn"

  # Validate method values
  valid_methods <- if (is_knn_mode) c("euclidean", "manhattan") else c("regularized", "EM")
  bad_method <- vapply(group$parameters, function(p) {
    !is.null(p$method) && !(p$method %in% valid_methods)
  }, logical(1))
  if (any(bad_method)) {
    stop(sprintf(
      "Invalid `method` for %s in group(s) %s. Must be one of: %s",
      toupper(imp_method),
      paste(which(bad_method), collapse = ", "),
      paste(valid_methods, collapse = ", ")
    ))
  }

  # Validate parameter names
  allowed_params <- if (is_knn_mode) {
    c("k", "method", "colmax", "post_imp", "dist_pow", "tree", "max_cache")
  } else {
    c(
      "ncp", "scale", "method", "coeff.ridge", "row.w",
      "threshold", "seed", "nb.init", "maxiter", "miniter"
    )
  }
  all_param_names <- unique(unlist(lapply(group$parameters, names)))
  unknown_params <- setdiff(all_param_names, allowed_params)
  if (length(unknown_params) > 0) {
    stop(
      "Unknown parameters for ", imp_method, " method: ",
      paste(unknown_params, collapse = ", ")
    )
  }

  # Cap per-group k/ncp
  feat_lengths <- lengths(group$feature)
  aux_lengths <- lengths(group$aux)
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
      stop(sprintf(
        "Group %d: %s must be >= 1 after capping (group size = %d).",
        i, required_param, group_size[i]
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
  feat_splits <- collapse::gsplit(all_feats_pos, rep(iter, feat_lengths))

  if (any(aux_lengths > 0)) {
    all_aux_pos <- collapse::fmatch(unlist(group$aux), cn)
    aux_splits <- collapse::gsplit(all_aux_pos, rep(iter, aux_lengths))
  } else {
    aux_splits <- vector("list", nrow(group))
  }

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
        message("OpenMP not available. Imputation will run single-threaded.")
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
    if (is_knn_mode) {
      p$cores <- cores
      p$subset <- group$feature[[i]]
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
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
        src <- bigmemory::attach.big.matrix(big_obj_desc)
        dst <- bigmemory::attach.big.matrix(big_out_desc)
        sub_mat <- src[, indices[[i]]$col_idx, drop = FALSE]
        imputed <- do.call(imp_fn, c(list(obj = sub_mat), params[[i]]))
        dst[, out_ranges[[i]]] <- imputed[, indices[[i]]$features_idx_local, drop = FALSE]
        invisible(NULL)
      },
      big_obj_desc = big_obj_desc,
      big_out_desc = big_out_desc,
      imp_fn = imp_fn,
      indices = indices,
      params = params,
      out_ranges = out_ranges
    )
    m <- mirai::mirai_map(iter, crated_fn)
    m[.progress = .progress]
    obj[, all_feats_pos] <- big_out[, ]
  } else {
    if (.progress) pb <- cli::cli_progress_bar(total = length(iter))
    for (i in iter) {
      sub_mat <- obj[, indices[[i]]$col_idx, drop = FALSE]
      imputed <- do.call(imp_fn, c(list(obj = sub_mat), params[[i]]))
      obj[, feat_splits[[i]]] <- imputed[, indices[[i]]$features_idx_local, drop = FALSE]
      if (.progress) cli::cli_progress_update(id = pb)
    }
    if (.progress) cli::cli_progress_done(id = pb)
  }

  class(obj) <- c("slideimp_results", class(obj))
  attr(obj, "imp_method") <- imp_method
  return(obj)
}
