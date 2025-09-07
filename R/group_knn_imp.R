#' Grouped K-NN Imputation
#'
#' K-NN imputation by groups, such as chromosomes, flanking columns, or clusters identified by column clustering techniques.
#'
#' @inheritParams knn_imp
#' @param group A data.frame with columns:
#' \describe{
#' \item{features}{A list column containing character vectors of feature column names to impute}
#' \item{aux}{A list column containing character vectors of auxiliary column names used for imputation but not imputed themselves}
#' \item{parameters}{(Optional) A list column containing group-specific parameters that override global settings. Allowed parameters: `k`, `weighted`, `method`, `tree`, `dist_pow`}
#' }
#' @param .progress Shows progress
#'
#' @details
#' This function performs K-NN imputation on groups of features independently, which will significantly
#' reduce imputation time for large datasets. Typical strategies for grouping may include:
#' \itemize{
#' \item Breaking down search space by chromosomes
#' \item Grouping features with their flanking values/neighbors (e.g., 1000 bp down/up stream of a CpG)
#' \item Using clusters identified by column clustering techniques
#' }
#'
#' Only features in each group (each row of the data.frame) will be imputed, using the search space
#' defined as the union of the features and aux columns of that group. Columns that are in aux or in the object
#' but not in any features will be left unchanged.
#'
#' @returns A list of length `n_imp` containing numeric matrices or [bigmemory::big.matrix()]
#' objects (if `output` is specified) where only imputed features are returned. Missing values
#' are imputed using k-NN for columns with missingness below `colmax`, and mean
#' imputation for remaining missing values if `post_imp = TRUE`. Only the features
#' specified in the groups are imputed; other columns remain unchanged.
#'
#' Each list element represents an independent imputation. The only element of the
#' list when `n_pmm == -1` is a single imputed matrix.
#'
#' @inherit knn_imp note
#'
#' @seealso [SlideKnn()], [mean_impute_col()], [bigmemory::big.matrix()], [restore_dimnames()], [knn_imp()]
#'
#' @export
#'
#' @examples
#' # Generate example data with missing values. This simulates a 20x50 matrix with missing values
#' # and groups by chromosome. Here we are simulating 2 chromosomes.
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
#'   parameters = list(list(k = 3, weighted = TRUE), list(k = 4, method = "manhattan"))
#' )
#' group_df
#'
#' # Run grouped imputation. t() to put features on the columns
#' obj <- t(to_test$input)
#' grouped_results <- group_knn_imp(obj, group = group_df, k = 5)
#' grouped_results
group_knn_imp <- function(
    obj,
    group,
    k,
    colmax = 0.9,
    rowmax = 0.9,
    method = c("euclidean", "manhattan"),
    cores = 1,
    post_imp = TRUE,
    weighted = FALSE,
    dist_pow = 1,
    tree = NULL,
    n_imp = 1,
    n_pmm = -1,
    seed = 42,
    output = NULL,
    overwrite = FALSE,
    .progress = TRUE) {
  # most pre-conditioning handled by knn_imp
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, null.ok = FALSE, .var.name = "obj")
  checkmate::assert_data_frame(group, min.rows = 1)
  checkmate::assert_names(colnames(group), must.include = c("features", "aux"))
  checkmate::assert_list(group$features, types = "character", min.len = 1)
  checkmate::assert_list(group$aux, types = c("character", "null"), min.len = 1)
  # check for optional parameters column in group
  has_parameters <- "parameters" %in% colnames(group)
  if (has_parameters) {
    checkmate::assert_list(group$parameters, types = c("list", "null"), min.len = 1)
    # only these params will be group wise different-able
    allowed_params <- c("k", "weighted", "method", "tree", "dist_pow")
    all_param_names <- unique(unlist(lapply(group$parameters, names)))
    unknown_params <- setdiff(all_param_names, allowed_params)
    if (length(unknown_params) > 0) {
      stop("Unknown parameters in group$parameters: ", paste(unknown_params, collapse = ", "))
    }
    message("Running with group-wise parameters")
  } else {
    message("Running with the same parameters for all groups")
  }
  checkmate::assert_int(seed, lower = 0, .var.name = "seed", null.ok = FALSE)
  set.seed(seed)
  # file backed logic
  if (!is.null(output)) {
    checkmate::assert_path_for_output(output, .var.name = "output")
  }
  checkmate::assert_flag(overwrite, null.ok = FALSE, .var.name = "overwrite")
  # feats and aux pre-conditioning
  rn <- rownames(obj)
  cn <- colnames(obj)
  if (is.null(cn)) {
    stop("`obj` must have column names for grouping")
  }
  # `all_feats` doesn't have to cover all cn. Uncovered columns are unimputed
  all_feats <- purrr::list_c(group$features)
  all_aux <- purrr::list_c(group$aux)
  if (!all(unique(c(all_feats, all_aux)) %in% cn)) {
    stop("Some features or aux columns not found in `obj` column names")
  }
  if (any(duplicated(all_feats))) {
    stop("Same features can't be in more than 1 groups")
  }
  file_backed <- !is.null(output)
  # copy pasted logic from knn_imp
  paths <- check_result_list(output = output, n_imp = n_imp, overwrite = overwrite)
  result_list <- create_result_list(
    data_to_copy = obj,
    file_backed = file_backed,
    n_imp = n_imp,
    backfiles = paths$backfiles,
    descfiles = paths$descfiles,
    backingpath = paths$backingpath
  )
  # dimnames handling
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = TRUE)
  can_store_dimnames <- isTRUE(getOption("bigmemory.allow.dimnames")) || !file_backed
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)
  # main for loop over groups. Parallel at `knn_imp` level only for minimal overhead
  for (g in seq_len(nrow(group))) {
    if (.progress) {
      message("Imputing group ", g, "/", nrow(group))
    }
    # behavior: all unique columns in feat will be imputed using unique columns in
    # feats and aux. Probably have to docs this
    feats <- unique(group$features[[g]])
    aux <- unique(group$aux[[g]])
    columns <- unique(c(feats, aux))
    # use the do.call method. First get all the needed params.
    group_params <- list(
      obj = obj[, columns, drop = FALSE],
      k = k,
      colmax = colmax,
      rowmax = rowmax,
      method = method,
      cores = cores,
      post_imp = post_imp,
      subset = feats, # character vector of column names to impute
      weighted = weighted,
      dist_pow = dist_pow,
      tree = tree,
      n_imp = n_imp,
      n_pmm = n_pmm,
      seed = seed,
      output = NULL, # Force in-memory to avoid nested file-backing. Won't take
      # too much memory since groups are assumed to be smaller than obj
      overwrite = FALSE
    )
    # override with group-specific parameters if provided
    if (has_parameters && !is.null(group$parameters[[g]])) {
      group_specific <- group$parameters[[g]]
      # only allow overriding specific imputation related parameters
      for (param_name in names(group_specific)) {
        group_params[[param_name]] <- group_specific[[param_name]]
      }
    }
    # call `knn_imp` with the prepared parameters
    sub_imp <- do.call(knn_imp, group_params)
    # Update the imputed features in the main result_list
    for (j in seq_len(n_imp)) {
      result_list[[j]][, feats] <- sub_imp[[j]][, feats]
    }
  }
  if (file_backed) {
    for (i in seq_len(n_imp)) {
      bigmemory::flush(result_list[[i]])
    }
  }
  # S3 OOPs
  if (can_store_dimnames) {
    for (i in seq_len(n_imp)) {
      rownames(result_list[[i]]) <- rn
      colnames(result_list[[i]]) <- cn
    }
  }
  attr(result_list, "rownames") <- rn
  attr(result_list, "colnames") <- cn
  attr(result_list, "subset") <- match(all_feats, cn)
  attr(result_list, "ncol") <- ncol(obj)
  class(result_list) <- c("KnnImpList", class(result_list))
  return(result_list)
}
