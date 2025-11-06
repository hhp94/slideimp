#' Grouped K-NN Imputation
#'
#' K-NN imputation by groups, such as chromosomes, flanking columns, or clusters
#' identified by column clustering techniques.
#'
#' @inheritParams knn_imp
#' @inheritParams pca_imp
#' @inheritParams slide_imp
#' @param group A data.frame with columns:
#' \describe{
#' \item{features}{A list column containing character vectors of feature column names to impute}
#' \item{aux}{(Optional) A list column containing character vectors of auxiliary
#' column names used for imputation but not imputed themselves}
#' \item{parameters}{(Optional) A list column containing group-specific parameters}
#' }
#'
#' @details
#' This function performs K-NN or PCA imputation on groups of features independently,
#' which significantly reduce imputation time for large datasets. Typical strategies
#' for grouping may include:
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
#' group_df <- tibble::tibble(
#'   features = list(group_1[1:3], group_2[1:4]),
#'   aux = list(group_1, group_2),
#'   parameters = list(
#'     list(k = 3, dist_pow = 1),
#'     list(k = 4, method = "manhattan")
#'   )
#' )
#' group_df
#'
#' # Run grouped imputation. t() to put features on the columns
#' obj <- t(to_test$input)
#' grouped_results <- group_imp(obj, group = group_df, k = 5)
#' grouped_results
group_imp <- function(
  obj,
  group,
  # KNN-specific parameters
  k = NULL,
  colmax = 0.9,
  knn_method = c("euclidean", "manhattan"),
  cores = 1,
  post_imp = TRUE,
  dist_pow = 0,
  tree = NULL,
  # PCA-specific parameters
  ncp = NULL,
  scale = TRUE,
  pca_method = c("Regularized", "EM"),
  coeff.ridge = 1,
  ind.sup = NULL,
  threshold = 1e-6,
  seed = NULL,
  nb.init = 1,
  maxiter = 1000,
  miniter = 5,
  # Others
  .progress = TRUE
) {
  # pre-conditioning
  if (is.null(k) && is.null(ncp)) {
    stop("Specify either 'k' for K-NN imputation or 'ncp' for PCA imputation")
  }
  imp_method <- if (!is.null(k)) {
    "knn"
  } else {
    "pca"
  }
  checkmate::assert_data_frame(group, min.rows = 1, .var.name = "group")
  checkmate::assert_names(colnames(group), must.include = c("features"), .var.name = "group")
  checkmate::assert_list(group$features, types = "character", min.len = 1, .var.name = "group$features")
  if ("aux" %in% names(group)) {
    checkmate::assert_list(group$aux, types = c("character", "null"), min.len = 1, .var.name = "group$aux")
  } else {
    group$aux <- list(NULL)
  }
  # handling group wise parameters
  has_parameters <- "parameters" %in% names(group)
  if (has_parameters) {
    checkmate::assert_list(group$parameters, types = c("list", "null"), min.len = 1, .var.name = "group$parameters")
    # allowed parameters based on method
    if (imp_method == "knn") {
      allowed_params <- c(
        "k", "method", "colmax", "cores", "post_imp", "dist_pow", "tree"
      )
    } else {
      allowed_params <- c(
        "ncp", "scale", "method", "row.w", "ind.sup", "quanti.sup", "coeff.ridge",
        "threshold", "seed", "nb.init", "maxiter", "miniter"
      )
    }
    all_param_names <- unique(unlist(lapply(group$parameters, names)))
    unknown_params <- setdiff(all_param_names, allowed_params)
    if (length(unknown_params) > 0) {
      stop(
        "Unknown parameters in group$parameters for ", imp_method, " method: ",
        paste(unknown_params, collapse = ", ")
      )
    }
    message("Running with group-wise parameters")
  } else {
    message("Running with the same parameters for all groups")
  }
  # feats and aux pre-conditioning
  cn <- colnames(obj)
  if (is.null(cn)) {
    stop("`obj` must have column names for grouping")
  }
  # `all_feats` doesn't have to cover all cn. Uncovered columns are un-imputed
  all_feats <- purrr::list_c(group$features)
  all_aux <- purrr::list_c(group$aux)
  if (!all(unique(c(all_feats, all_aux)) %in% cn)) {
    stop("Some features or aux columns not found in `obj` column names")
  }
  if (any(duplicated(all_feats))) {
    stop("Same features can't be in more than 1 groups")
  }
  for (g in seq_len(nrow(group))) {
    if (.progress) {
      message("Imputing group ", g, "/", nrow(group))
    }
    # behavior: all unique columns in feat will be imputed using unique columns in
    # feats and aux. Probably have to docs this
    feats <- unique(group$features[[g]])
    aux <- unique(group$aux[[g]])
    columns <- unique(c(feats, aux))

    if (imp_method == "knn") {
      f_imp <- "knn_imp"
      group_params <- list(
        obj = obj[, columns, drop = FALSE],
        k = k,
        colmax = colmax,
        method = knn_method,
        cores = cores,
        post_imp = post_imp,
        subset = feats, # character vector of column names to impute
        dist_pow = dist_pow,
        tree = tree
      )
    } else {
      f_imp <- "pca_imp"
      group_params <- list(
        obj = obj[, columns, drop = FALSE],
        ncp = ncp,
        scale = scale,
        method = pca_method,
        # row.w = row.w,
        ind.sup = ind.sup,
        quanti.sup = NULL,
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
      # only allow overriding specific imputation related parameters
      for (param_name in names(group_specific)) {
        group_params[[param_name]] <- group_specific[[param_name]]
      }
    }
    # Call the imputation function
    obj[, feats] <- do.call(f_imp, group_params)[, feats]
  }

  return(obj)
}
