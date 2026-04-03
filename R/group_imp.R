group_indices <- function(g, feat_splits, aux_splits, group_features) {
  f_idx <- feat_splits[[g]]
  col_idx <- c(f_idx, aux_splits[[g]])
  list(
    features_idx_local = seq_along(f_idx),
    col_idx = col_idx,
    features_names = group_features[[g]]
  )
}

process_group_params <- function(g, base_params, group, is_knn_mode, has_parameters) {
  group_params <- base_params
  if (is_knn_mode) {
    group_params$subset <- group$features[[g]]
  }
  if (has_parameters && !is.null(group$parameters[[g]])) {
    group_specific <- group$parameters[[g]]
    group_specific[["cores"]] <- NULL
    group_params[names(group_specific)] <- group_specific
  }
  group_params
}

normalize_list_column <- function(x) {
  if (is.null(x) || is.logical(x)) {
    return(character(0))
  }
  x <- x[!is.na(x)]
  return(unique(x))
}

#' Grouped K-NN or PCA Imputation
#'
#' K-NN or PCA imputation by groups, such as chromosomes, flanking columns, or clusters
#' identified by column clustering techniques.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#' @inheritParams slide_imp
#' @param group A data.frame/[tibble::tibble()] describing how features should be grouped for
#'   imputation. Accepts two formats:
#'
#'   **Long format**:
#'   - **group**: A column identifying which group each feature belongs to
#'   - **features**: A character column of individual feature names
#'
#'   **List-column format** (preferably created by [group_features()]):
#'   - **features**: A list-column of character vectors of feature column names
#'     to impute
#'   - **aux**: (Optional) A list-column of character vectors of auxiliary column
#'     names used for imputation but not imputed themselves
#'   - **parameters**: (Optional) A list-column of group-specific parameter lists
#'
#' @param cores Controls the number of cores to parallelize over for K-NN imputation with OpenMP only.
#' To setup parallelization for K-NN without OpenMP and PCA imputation, use `mirai::daemons()`.
#' @param .progress Show imputation progress (default = FALSE)
#'
#' @details
#' This function performs K-NN or PCA imputation on groups of features independently,
#' which significantly reduce imputation time for large datasets.
#'
#' Specify `k` and related arguments to use K-NN, `ncp` and related arguments for PCA imputation.
#' If `k` and `ncp` are both `NULL`, then the group-wise parameters column i.e., `group$parameters`
#' must be specified and must contains either `k` or `ncp` for all groups of group-wise parameters.
#'
#' Strategies for grouping may include:
#' - Breaking down search space by chromosomes
#' - Grouping features with their flanking values/neighbors (e.g., 1000 bp down/up stream of a CpG)
#' - Using clusters identified by column clustering techniques
#'
#' Only features in each group will be imputed, using the search space defined as the
#' union of the features and optional aux columns of that group. Columns that are in
#' aux or in the object but not in any features will be left unchanged.
#'
#' @inherit knn_imp note return
#'
#' @export
#'
#' @seealso [group_features()]
#'
#' @examples
#' # Generate example data with missing values
#' set.seed(1234)
#' to_test <- sim_mat(
#'   m = 20,
#'   n = 50,
#'   perc_NA = 0.3,
#'   perc_col_NA = 1,
#'   nchr = 2
#' )
#' # t() to put features in columns
#' obj <- t(to_test$input)
#' head(to_test$group_feature) # which group each feature belongs to
#'
#' # Long format API: pass a data.frame with 'group' and 'features' columns.
#' long_fmt <- to_test$group_feature
#' head(long_fmt)
#' knn_long <- group_imp(obj, group = long_fmt, k = 5, cores = 2)
#' knn_long
#'
#' # List-column format API: Use group_features() to create the group tibble.
#' # By setting `k = 5` in group_features(), we are doing K-NN imputation in group_imp().
#' # To make use of the `subset` argument in knn_imp(), we specify subset in group_features().
#' # For demonstration of different group-wise parameters we set `k = 10` for the second group.
#'
#' subset_features <- sample(to_test$group_feature$features, size = 10)
#' head(subset_features)
#' knn_df <- group_features(obj, to_test$group_feature, k = 5, subset = subset_features)
#' knn_df
#' knn_df$parameters[[2]]$k <- 10
#' knn_df$parameters
#'
#' # Run grouped imputation. `k` for K-NN has been specified in `knn_df`.
#' knn_grouped <- group_imp(obj, group = knn_df, cores = 2)
#' knn_grouped # only features in subset are imputed
#' @examplesIf requireNamespace("carrier", quietly = TRUE)
#' # Specify `ncp` for PCA directly in the group_imp() function (instead of in
#' # group_features()). We run in parallel with `mirai::daemons(2)`.
#' mirai::daemons(2) # Set up 2 cores for parallelization
#' pca_df <- group_features(obj, to_test$group_feature)
#' pca_grouped <- group_imp(obj, group = pca_df, ncp = 2)
#' mirai::daemons(0)
#' pca_grouped
#'
group_imp <- function(
  obj,
  group,
  # KNN-specific parameters
  k = NULL,
  colmax = NULL,
  knn_method = NULL,
  post_imp = NULL,
  dist_pow = NULL,
  tree = NULL,
  cores = 1,
  max_cache = NULL,
  # PCA-specific parameters
  ncp = NULL,
  scale = NULL,
  pca_method = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  row.w = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL,
  .progress = TRUE
) {
  # pre-conditioning  ----
  checkmate::assert_matrix(obj, mode = "numeric", col.names = "unique", null.ok = FALSE, .var.name = "obj")
  checkmate::assert_data_frame(group, min.rows = 1, .var.name = "group")
  checkmate::assert_names(colnames(group), must.include = c("features"), .var.name = "group")
  if (!tibble::is_tibble(group)) {
    group <- tibble::as_tibble(group)
  }
  if ("group" %in% names(group) && inherits(group$features, "character")) {
    stopifnot("NA is not allowed in group$group" = !anyNA(group$group))
    checkmate::assert_character(group$features, any.missing = FALSE, .var.name = "group$features")
    if (is.numeric(group$group) && !checkmate::test_integerish(group$group)) {
      warning(
        "Non-integer numeric 'group$group' values detected. ",
        "`group_imp` reserves the name `group` in the `group` parameter ",
        "to group the features in the `features` column. ",
        "Did you mean to use integer or character labels for `group$group`?"
      )
    }
    group <- collapse::fsummarize(collapse::fgroup_by(group, group), features = list(features))
    group$group <- NULL
  } else if (inherits(group$features, "character")) {
    stop(
      "`group` has a character 'features' column but no 'group' column. ",
      "Either add a 'group' column to define groups, or use list-columns (see ?group_features)."
    )
  }
  checkmate::assert_list(group$features,
    types = "character", min.len = 1,
    .var.name = "group$features"
  )
  group$features <- lapply(group$features, normalize_list_column)
  feat_lengths <- lengths(group$features)
  if (any(feat_lengths == 0L)) {
    stop(
      sprintf("Group(s) %s have no features.", paste(which(feat_lengths == 0L), collapse = ", "))
    )
  }
  if ("aux" %in% names(group)) {
    checkmate::assert_list(
      group$aux,
      types = c("character", "logical", "null"),
      min.len = 1, .var.name = "group$aux"
    )
    group$aux <- lapply(group$aux, normalize_list_column)
  } else {
    group$aux <- list(character(0))
  }
  aux_only <- purrr::map2(group$aux, group$features, setdiff)
  aux_lengths <- lengths(aux_only)
  ## k, ncp, and group$parameters logic ----
  global_k <- !is.null(k)
  global_ncp <- !is.null(ncp)
  if (global_k && global_ncp) {
    stop("Cannot specify both 'k' and 'ncp' as global parameters")
  }
  has_parameters <- "parameters" %in% names(group)
  use_global <- global_k || global_ncp
  if (use_global) {
    imp_method <- if (global_k) "knn" else "pca"
    if (use_global && has_parameters) {
      stop("Cannot specify both global 'k'/'ncp' and group$parameters. Please use one or the other")
    }
  } else {
    if (!has_parameters) {
      stop("Must specify either global 'k' for K-NN imputation, 'ncp' for PCA imputation, or provide group$parameters")
    }
    checkmate::assert_list(group$parameters, types = c("list", "null"), min.len = 1, .var.name = "group$parameters")
    has_k <- vapply(group$parameters, \(p) "k" %in% names(p), logical(1))
    has_ncp <- vapply(group$parameters, \(p) "ncp" %in% names(p), logical(1))
    both <- has_k & has_ncp
    neither <- !has_k & !has_ncp
    if (any(both)) {
      stop(sprintf("Group(s) %s have both 'k' and 'ncp' in parameters", paste(which(both), collapse = ", ")))
    }
    if (any(neither)) {
      stop(sprintf("Group(s) %s have neither 'k' nor 'ncp' in parameters", paste(which(neither), collapse = ", ")))
    }
    if (any(has_k) && any(has_ncp)) {
      stop("Inconsistent imputation methods across groups; all group-wise parameters must specify either just k for KNN imputation or ncp for PCA imputation")
    }
    imp_method <- if (all(has_k)) "knn" else "pca"
    if (imp_method == "knn") {
      allowed_params <- c(
        "k", "method", "colmax", "cores", "post_imp", "dist_pow", "tree",
        "max_cache"
      )
      required_param <- "k"
    } else {
      allowed_params <- c(
        "ncp", "scale", "method", "coeff.ridge", "row.w",
        "threshold", "seed", "nb.init", "maxiter", "miniter"
      )
      required_param <- "ncp"
    }
    all_param_names <- unique(unlist(lapply(group$parameters, names)))
    unknown_params <- setdiff(all_param_names, allowed_params)
    if (length(unknown_params) > 0) {
      stop(
        "Unknown parameters in group$parameters for ", imp_method, " method: ",
        paste(unknown_params, collapse = ", ")
      )
    }
  }
  message(sprintf(
    "Imputing %d group(s) using %s with %s parameters.",
    nrow(group), toupper(imp_method), if (use_global) "global" else "group-wise"
  ))
  group_size <- feat_lengths + aux_lengths
  if (has_parameters) {
    param_values <- vapply(group$parameters, function(p) p[[required_param]], numeric(1))
  } else {
    param_values <- rep(if (imp_method == "knn") k else ncp, length(group_size))
  }
  if (imp_method == "knn") {
    invalid_groups <- which(param_values >= group_size | param_values < 1)
  } else {
    invalid_groups <- which(param_values > pmin(group_size, nrow(obj) - 1) | param_values < 1)
  }
  if (length(invalid_groups) > 0) {
    stop(sprintf(
      "%s is either too large or < 1 for group(s): %s",
      if (imp_method == "knn") "k" else "ncp",
      paste(invalid_groups, collapse = ", ")
    ))
  }
  cn <- colnames(obj)
  # `all_feats` doesn't have to cover all cn. Uncovered columns are un-imputed
  all_feats <- unlist(group$features)
  all_aux <- unlist(group$aux)
  if (!all(unique(c(all_feats, all_aux)) %in% cn)) {
    stop("Some features or aux columns not found in `obj` column names")
  }
  if (anyDuplicated(all_feats) > 0) {
    stop("Same features can't be in more than 1 groups")
  }
  iter <- seq_len(nrow(group))

  # Column-index lookups for each group ----
  # Step 1: Map feature names to column positions in `obj`, in one batch.
  # Instead of looping per-group with match(), we unlist everything,
  # do a single fmatch(), then split the result back by group.
  all_feats_pos <- collapse::fmatch(all_feats, cn) # integer positions in obj
  # the recreate per-group splits: group 1 gets first feat_lengths[1] positions, etc.
  feat_splits <- collapse::gsplit(all_feats_pos, rep(iter, feat_lengths))
  # Step 2: Check if any aux-only columns exist (aux_only computed earlier).
  # Step 3: Same batch lookup for aux columns, but only if any exist.
  if (any(aux_lengths > 0)) {
    all_aux_pos <- collapse::fmatch(unlist(aux_only), cn)
    aux_splits <- collapse::gsplit(all_aux_pos, rep(iter, aux_lengths))
  } else {
    aux_splits <- vector("list", nrow(group))
  }
  # Step 4: Assemble per-group index bundles.
  # Each element contains:
  # - col_idx: combined column positions (features first, then aux)
  # used to subset obj into a sub-matrix for imputation
  # - features_idx_local: 1:window_sizeures where the feature columns sit inside
  # the sub-matrix (always the first columns, since features are pre-pended before
  # aux in col_idx)
  # - features_names: original column names, used to write results back
  indices <- lapply(
    iter,
    group_indices,
    feat_splits = feat_splits,
    aux_splits = aux_splits,
    group_features = group$features
  )

  # parallel logic ----
  is_knn_mode <- imp_method == "knn"
  parallelize <- tryCatch(mirai::require_daemons(), error = function(e) FALSE)

  if (cores > 1 & is_knn_mode) {
    if (!has_openmp()) {
      message("OpenMP not available. KNN will run single-threaded.")
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
  if (cores > 1 & !is_knn_mode) {
    message(
      sprintf(
        "cores = %d but is ignored for PCA imputation. Call `mirai::daemons(%d)` to set up parallelization.",
        cores,
        cores
      )
    )
  }

  # group level parameters ----
  if (imp_method == "knn") {
    base_params <- list(
      k = k, colmax = colmax, method = knn_method,
      cores = cores, post_imp = post_imp,
      dist_pow = dist_pow, tree = tree
    )
  } else {
    base_params <- list(
      ncp = ncp, scale = scale, method = pca_method,
      coeff.ridge = coeff.ridge, threshold = threshold,
      seed = seed, row.w = row.w, nb.init = nb.init,
      maxiter = maxiter, miniter = miniter
    )
  }
  # Strip NULLs
  base_params <- base_params[!vapply(base_params, is.null, logical(1))]

  params <- lapply(
    iter,
    process_group_params,
    base_params = base_params,
    group = group,
    is_knn_mode = is_knn_mode,
    has_parameters = has_parameters
  )

  # imputation ----
  imp_fn <- if (is_knn_mode) knn_imp else pca_imp

  if (parallelize || (is_knn_mode && cores > 1)) {
    message("Running Mode: parallel...")
  } else {
    message("Running Mode: sequential...")
  }

  if (parallelize) {
    # overall scheme: convert the beta matrix to a big.matrix. Also, create a
    # nrow(obj)*imputed features upfront called big_out. This is so that we don't have to
    # accumulate object in the purrr loop. The purrr loop will just write the
    # imputed CpGs to the big_out matrix in parallel. No race conds will be there
    # because we already enforced uniqueness of `group$features`. This means that
    # we have to map the incontiguous group$features index into contiguous
    # big_out matrix.
    feat_cumsum <- cumsum(c(0L, feat_lengths))
    out_ranges <- lapply(iter, function(i) {
      (feat_cumsum[i] + 1L):feat_cumsum[i + 1L]
    })

    big_obj <- bigmemory::as.big.matrix(obj)
    big_obj_desc <- bigmemory::describe(big_obj)
    big_out <- bigmemory::big.matrix(
      nrow = nrow(obj), ncol = length(all_feats_pos), type = "double"
    )
    big_out_desc <- bigmemory::describe(big_out)

    crated_fn <- purrr::in_parallel(
      function(i) {
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
    purrr::walk(iter, crated_fn, .progress = .progress)
    obj[, all_feats_pos] <- big_out[, ]
  } else {
    if (.progress) pb <- cli::cli_progress_bar(total = length(iter))
    for (i in iter) {
      sub_mat <- obj[, indices[[i]]$col_idx, drop = FALSE]
      imputed <- do.call(imp_fn, c(list(obj = sub_mat), params[[i]]))
      obj[, feat_splits[[i]]] <- imputed[, indices[[i]]$features_idx_local, drop = FALSE]
      if (.progress) {
        cli::cli_progress_update(id = pb)
      }
    }
    if (.progress) {
      cli::cli_progress_done(id = pb)
    }
  }

  class(obj) <- c("ImputedMatrix", class(obj))
  attr(obj, "imp_method") <- imp_method
  return(obj)
}

#' Group Features for Imputation
#'
#' Groups matrix columns (features) based on a provided grouping data.frame,
#' optionally preparing parameters for K-NN or PCA imputation. This function
#' organizes features into groups, handles imputation of only a subset of features,
#' and can pad groups to meet a minimum size.
#'
#' @inheritParams knn_imp
#' @param features_df A data.frame with exactly two columns: `features` and
#' `group`. Maps feature identifiers to their respective groups. No missing
#' values or duplicate `features` values are allowed.
#' @param k Integer or `NULL`. If specified, prepares parameters for K-NN
#' imputation with `k` neighbors. Cannot be used together with `ncp`.
#' @param ncp Integer or `NULL`. If specified, prepares parameters for PCA
#' imputation with `ncp` principal components. Cannot be used together with
#' `k`.
#' @param min_group_size Integer (default = `0`). Minimum number of features per
#' group. If a group has fewer features, additional features are randomly
#' sampled from remaining columns to meet this threshold.
#' @param seed Numeric or `NULL`. Random seed for reproducibility when sampling
#' for `min_group_size` padding.
#'
#' @return A `tibble::tibble()` with columns:
#'
#' - **features**: A list-column containing character vectors of feature column
#' names to impute
#' - **aux**: A list-column containing character vectors of auxiliary
#' column names used for imputation but not imputed themselves. Omitted if all
#' elements are `NULL`
#' - **parameters**: A list-column containing group-specific parameters if `k`
#' or `ncp` are specified
#'
#' @examples
#' sim_obj <- sim_mat(perc_col_NA = 1)
#'
#' obj <- t(sim_obj$input)
#' obj_meta <- sim_obj$group_feature
#'
#' # group `obj` based on the metadata
#' head(obj_meta)
#'
#' # create `group_df` which can then be used for `group_imp`. We can specify
#' # `k` for K-NN imputation and subset here as well.
#' group_df <- group_features(
#'   obj,
#'   obj_meta,
#'   subset = sample(obj_meta$features, size = 10),
#'   k = 10
#' )
#' group_df
#'
#' imputed_obj <- group_imp(obj, group_df)
#' imputed_obj
#'
#' @seealso [group_imp()]
#'
#' @export
group_features <- function(
  obj,
  features_df,
  k = NULL,
  ncp = NULL,
  subset = NULL,
  min_group_size = 0,
  seed = NULL
) {
  # pre-conditioning
  checkmate::assert_matrix(
    obj,
    mode = "numeric",
    row.names = "named",
    col.names = "unique",
    null.ok = FALSE,
    .var.name = "obj"
  )
  cn <- colnames(obj)
  checkmate::assert_data_frame(
    features_df,
    any.missing = FALSE,
    min.rows = 1,
    ncols = 2,
    null.ok = FALSE,
    .var.name = "features_df"
  )
  checkmate::assert_subset(
    colnames(features_df),
    c("features", "group"),
    .var.name = "features_df"
  )
  stopifnot(
    "`features_df$features` cannot contain duplicated values" =
      anyDuplicated(features_df$features) == 0L
  )
  if (!is.null(k) && !is.null(ncp)) {
    stop("Specify either `k` for K-NN imputation or `ncp` for PCA imputation or neither, not both")
  }
  if (!is.null(k)) {
    checkmate::assert_int(k, lower = 1, .var.name = "k")
  }
  if (!is.null(ncp)) {
    checkmate::assert_int(ncp, lower = 1, .var.name = "ncp")
  }
  checkmate::assert_int(min_group_size, lower = 0, .var.name = "min_group_size")
  checkmate::assert_character(
    subset,
    min.len = 1,
    any.missing = FALSE,
    unique = TRUE,
    null.ok = TRUE,
    .var.name = "subset"
  )
  if (!is.null(subset)) {
    if (length(setdiff(subset, cn)) > 0) {
      warning("Some elements of `subset` not found in colnames(obj) and are excluded")
    }
    subset <- intersect(subset, cn)
    if (length(subset) == 0) {
      stop("No element in `subset` found in colnames(obj)")
    }
  }
  checkmate::assert_number(seed, null.ok = TRUE, .var.name = "seed")
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # cn_df is used for joining
  cn_df <- data.frame(features = cn)
  if (!is.null(subset)) {
    # If there are subset, move all non subsets to `aux`
    cn_df$subset <- as.logical(
      collapse::fmatch(cn_df$features, subset, nomatch = 0)
    )
  } else {
    # else everything is in features
    cn_df$subset <- TRUE
  }

  matched <- collapse::join(
    cn_df,
    features_df,
    on = "features",
    how = "left",
    verbose = FALSE
  )
  if (sum(!is.na(matched$group)) == 0) {
    stop("`colnames(obj)` matched with no group in `features_df$group`")
  }

  # `features` list-column, to be imputed
  features <- collapse::fsubset(matched, subset)
  features <- tibble::tibble(
    features = collapse::gsplit(features$features, features$group, use.g.names = TRUE),
    group = names(features)
  )

  # `aux` list=column, first part is to handle the subset function
  aux <- collapse::fsubset(matched, !subset)
  aux <- tibble::tibble(
    aux = collapse::gsplit(aux$features, aux$group, use.g.names = TRUE),
    group = names(aux)
  )

  # combine
  group <- collapse::join(features, aux, on = "group", how = "left", verbose = FALSE)
  group$length <- lengths(group$features) + lengths(group$aux)

  # pad groups to `min_group_size` if needed
  if (min_group_size > 0) {
    group$need <- pmax(min_group_size - group$length, 0)
    group$min_group_size <- purrr::map2(group$features, group$need, \(x, y) {
      if (y == 0) {
        return(NULL)
      }
      pool <- setdiff(matched$features, x)
      if (length(pool) < y) {
        stop("`min_group_size` is too large")
      }
      sample(pool, size = y)
    })
    group$aux <- purrr::map2(group$aux, group$min_group_size, \(x, y) c(x, y))
    group[, c("need", "min_group_size")] <- NULL
    group$length <- lengths(group$features) + lengths(group$aux)
  }

  # add group specific `parameters` list column
  if (!is.null(k)) {
    group$parameters <- lapply(
      group$length,
      \(x) tibble::tibble(k = min(x - 1, k))
    )
  } else if (!is.null(ncp)) {
    group$parameters <- lapply(
      group$length,
      \(x) tibble::tibble(ncp = min(ncol(obj), x - 1, ncp))
    )
  }
  # do nothing if null k and ncp
  group[, "length"] <- NULL

  # clean up
  if (all(lengths(group$aux) == 0)) {
    group$aux <- NULL
  }

  for (i in names(group)) {
    if (is.list(group[[i]])) {
      names(group[[i]]) <- NULL
    }
  }

  return(group)
}
