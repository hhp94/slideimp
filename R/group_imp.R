#' Grouped K-NN or PCA Imputation
#'
#' K-NN or PCA imputation by groups, such as chromosomes, flanking columns, or clusters
#' identified by column clustering techniques.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#' @inheritParams slide_imp
#' @param group A data.frame/[tibble::tibble()] with columns:
#' \describe{
#' \item{features}{A list column containing character vectors of feature column names to impute}
#' \item{aux}{(Optional) A list column containing character vectors of auxiliary
#' column names used for imputation but not imputed themselves}
#' \item{parameters}{(Optional) A list column containing group-specific parameters}
#' }
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
#' Typical strategies for grouping may include:
#' \itemize{
#' \item Breaking down search space by chromosomes
#' \item Grouping features with their flanking values/neighbors (e.g., 1000 bp down/up stream of a CpG)
#' \item Using clusters identified by column clustering techniques
#' }
#'
#' Only features in each group (each row of the data.frame) will be imputed, using
#' the search space defined as the union of the features and optional aux columns
#' of that group. Columns that are in aux or in the object but not in any features
#' will be left unchanged.
#'
#' @inherit knn_imp note return
#'
#' @export
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
#'
#' # `group_1` will be all the CpGs on Chr1. Same for `group_2`
#' group_1 <- subset(to_test$group_feature, group == "chr1")$feature_id
#' group_2 <- subset(to_test$group_feature, group == "chr2")$feature_id
#'
#' # Impute only first 3 values of group 1, the rest are aux. Group 2 does 4 features.
#' # Also optionally vary the parameters by group
#' knn_df <- tibble::tibble(
#'   features = list(group_1[1:3], group_2[1:4]),
#'   aux = list(group_1, group_2),
#'   parameters = list(
#'     list(k = 3, dist_pow = 1),
#'     list(k = 4, method = "manhattan")
#'   )
#' )
#' knn_df
#'
#' # Run grouped imputation. t() to put features on the columns. `k` for K-NN has
#' # been specified in `knn_df`.
#' obj <- t(to_test$input)
#' knn_grouped <- group_imp(obj, group = knn_df)
#' knn_grouped
#'
#' # Specify `ncp` for PCA in the `group_imp` function since no group-wise parameters are
#' # specified.
#' pca_df <- tibble::tibble(
#'   features = list(group_1[1:3], group_2[1:4])
#' )
#' pca_grouped <- group_imp(obj, group = pca_df, ncp = 2)
#' pca_grouped
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
  # PCA-specific parameters
  ncp = NULL,
  scale = NULL,
  pca_method = NULL,
  coeff.ridge = NULL,
  threshold = NULL,
  seed = NULL,
  nb.init = NULL,
  maxiter = NULL,
  miniter = NULL,
  # Others
  cores = 1,
  .progress = TRUE
) {
  # pre-conditioning  ----
  checkmate::assert_matrix(obj, mode = "numeric", row.names = "named", col.names = "unique", null.ok = FALSE, .var.name = "obj")
  checkmate::assert_data_frame(group, min.rows = 1, .var.name = "group")
  checkmate::assert_names(colnames(group), must.include = c("features"), .var.name = "group")
  if (!tibble::is_tibble(group)) {
    group <- tibble::as_tibble(group)
  }
  checkmate::assert_list(group$features, types = "character", min.len = 1, .var.name = "group$features")
  for (i in seq_along(group$features)) {
    if (anyDuplicated(group$features[[i]]) > 0) {
      warning(paste("Group", i, "contains repeated features. Removing duplicates ..."))
      group$features[[i]] <- unique(group$features[[i]])
    }
  }
  features_length <- vapply(group$features, length, numeric(1))
  empty_features <- features_length == 0
  if (any(empty_features)) {
    stop(sprintf("Group(s) %s have no features.", paste(which(empty_features), collapse = ", ")))
  }
  if ("aux" %in% names(group)) {
    checkmate::assert_list(group$aux, types = c("character", "null"), min.len = 1, .var.name = "group$aux")
    for (i in seq_along(group$aux)) {
      if (anyDuplicated(group$aux[[i]]) > 0) {
        group$aux[[i]] <- unique(group$aux[[i]])
      }
    }
  } else {
    group$aux <- list(NULL)
  }
  aux_length <- vapply(group$aux, length, numeric(1))
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
    message(paste("Performing group-wise", toupper(imp_method), "imputation with the same parameters for all groups.", sep = " "))
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
      allowed_params <- c("k", "method", "colmax", "cores", "post_imp", "dist_pow", "tree")
      required_param <- "k"
    } else {
      allowed_params <- c(
        "ncp", "scale", "method", "row.w", "ind.sup", "quanti.sup", "coeff.ridge",
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
    message(paste("Performing group-wise", toupper(imp_method), "imputation with group-wise parameters.", sep = " "))
  }
  group_size <- features_length + aux_length
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
  if (is.null(cn)) {
    stop("`obj` must have column names for grouping")
  }
  # `all_feats` doesn't have to cover all cn. Uncovered columns are un-imputed
  all_feats <- do.call(c, group$features)
  all_aux <- do.call(c, group$aux)
  if (!all(unique(c(all_feats, all_aux)) %in% cn)) {
    stop("Some features or aux columns not found in `obj` column names")
  }
  if (anyDuplicated(all_feats) > 0) {
    stop("Same features can't be in more than 1 groups")
  }
  # imputation ----
  iter <- seq_len(nrow(group))
  indices <- lapply(iter, function(g) {
    # Have to convert to index to be able to use bigmemory
    columns <- unique(c(group$features[[g]], group$aux[[g]]))
    non_features <- setdiff(columns, group$features[[g]])
    features_idx <- match(group$features[[g]], cn)
    non_features_idx <- match(non_features, cn)
    return(
      list(
        features_idx_local = seq_along(features_idx),
        col_idx = c(features_idx, non_features_idx),
        features_names = cn[features_idx]
      )
    )
  })

  params <- lapply(iter, function(g) {
    if (imp_method == "knn") {
      group_params <- list(
        k = k,
        colmax = colmax,
        method = knn_method,
        cores = cores,
        post_imp = post_imp,
        subset = group$features[[g]], # only columns in features need to be imputed
        dist_pow = dist_pow,
        tree = tree
      )
    } else {
      # pca has to impute features + aux, but only features will be in the results
      group_params <- list(
        ncp = ncp,
        scale = scale,
        method = pca_method,
        coeff.ridge = coeff.ridge,
        threshold = threshold,
        seed = seed,
        nb.init = nb.init,
        maxiter = maxiter,
        miniter = miniter
      )
    }

    # use group-specific parameters if provided
    if (has_parameters && !is.null(group$parameters[[g]])) {
      group_specific <- group$parameters[[g]]
      for (param_name in setdiff(names(group_specific), "cores")) {
        group_params[[param_name]] <- group_specific[[param_name]]
      }
    }
    group_params <- group_params[!vapply(group_params, is.null, logical(1))]
    return(group_params)
  })

  if (imp_method == "knn") {
    if (.progress) {
      if (cores > 1) {
        message("Running in parallel...")
      } else {
        message("Running in sequential...")
      }
    }
    results <- purrr::map(
      iter,
      function(i) {
        do.call(
          knn_imp,
          c(list(obj = obj[, indices[[i]]$col_idx]), params[[i]])
        )[, indices[[i]]$features_idx_local]
      },
      .progress = .progress
    )
  } else {
    parallelize <- tryCatch(mirai::require_daemons(), error = function(e) FALSE)
    if (parallelize) {
      if (.progress) {
        message("Running in parallel...")
      }
      # have to use bigmemory big matrix. This will create only one copy in RAM.
      big_obj <- bigmemory::as.big.matrix(obj)
      big_obj_desc <- bigmemory::describe(big_obj)
      crated_fn <- purrr::in_parallel(
        function(i) {
          attached_matrix <- bigmemory::attach.big.matrix(big_obj_desc)
          imputed <- do.call(
            pca_imp,
            c(list(obj = attached_matrix[, indices[[i]]$col_idx]), params[[i]])
          )[, indices[[i]]$features_idx_local]
          colnames(imputed) <- indices[[i]]$features_names
          return(imputed)
        },
        big_obj_desc = big_obj_desc,
        pca_imp = pca_imp,
        indices = indices,
        params = params
      )
      results <- purrr::map(iter, crated_fn, .progress = .progress)
    } else {
      if (cores > 1) {
        warning(
          sprintf(
            "cores = %d but running in **sequential**. Call `mirai::daemons(%d)` to set up the parallelization",
            cores,
            cores
          )
        )
      }
      results <- purrr::map(
        iter,
        function(i) {
          do.call(pca_imp, c(list(obj = obj[, indices[[i]]$col_idx]), params[[i]]))[, indices[[i]]$features_idx_local]
        },
        .progress = .progress
      )
    }
  }

  for (res in results) {
    obj[, colnames(res)] <- res
  }

  class(obj) <- c("ImputedMatrix", class(obj))
  attr(obj, "imp_method") <- imp_method
  return(obj)
}
