#' K-Nearest Neighbor Imputation for Missing Values
#'
#' @description
#' Imputes missing values in numeric matrices using the k-nearest neighbors algorithm
#' with a two-stage approach: k-NN imputation for columns with missingness below a
#' threshold, followed by optional mean imputation for remaining missing values.
#'
#' @details
#' This function performs **column-wise** nearest neighbor calculations.
#'
#' When `weighted = TRUE`, imputed values are computed as distance-weighted averages
#' where weights are inverse distances raised to the power of `dist_pow`. This gives
#' closer neighbors greater influence in the imputation process, which can increase
#' predictive performance.
#'
#' The `tree` parameter enables faster neighbor search using spatial data structures
#' but requires pre-filling missing values with column means, which may introduce bias
#' in high-missingness data. Tree construction overhead may reduce performance for
#' low-dimensional data or small k values.
#'
#' @section Performance Optimization:
#' - **Tree methods**: Only use when imputation runtime becomes prohibitive and missingness is low (<20% missing)
#' - **Subset imputation**: Use `subset` parameter for efficiency when only specific columns need imputation
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @inheritParams slide_imp
#' @param k Number of nearest neighbors for imputation. 10 is a good starting point.
#' @param subset Character vector of column names or integer vector of column
#'   indices specifying which columns to impute.
#'
#' @return A matrix with `dim(obj)` with missing values imputed.
#'
#' @references
#' Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu (2002).
#' Diagnosis of multiple cancer types by shrunken centroids of gene expression
#' PNAS 99: 6567-6572. Available at www.pnas.org
#'
#' @examples
#' data(khanmiss1)
#' sum(is.na(khanmiss1))
#'
#' # Basic k-NN imputation (khanmiss1 has genes in rows, so transpose)
#' t_khanmiss1 <- t(khanmiss1)
#' result <- knn_imp(t_khanmiss1, k = 5)
#' result
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
  tree = NULL
) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, null.ok = FALSE, .var.name = "obj")
  checkmate::assert_int(k, lower = 1, upper = ncol(obj) - 1, .var.name = "k")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  checkmate::assert(
    checkmate::check_character(subset, min.len = 0, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
    checkmate::check_integerish(subset, lower = 1, upper = ncol(obj), min.len = 0, any.missing = FALSE, null.ok = TRUE, unique = TRUE),
    combine = "or",
    .var.name = "subset"
  )
  stopifnot(length(dist_pow) == 1, dist_pow >= 0, !is.infinite(dist_pow))
  checkmate::assert_choice(tree, choices = c("kd", "ball"), null.ok = TRUE, .var.name = "tree")
  # subset logic
  subset <- if (is.null(subset)) {
    seq_len(ncol(obj))
  } else if (length(subset) == 0) {
    integer(0)
  } else if (is.character(subset)) {
    if (is.null(colnames(obj))) {
      stop("`subset` contains characters but `obj` doesn't have column names")
    }
    matched_indices <- match(subset, colnames(obj), nomatch = NA)
    matched_indices[!is.na(matched_indices)]
  } else {
    subset
  }
  # Early return for empty subset or no missing data
  if (!anyNA(obj[, subset, drop = FALSE])) {
    return(obj)
  }
  # Calculate complement for further processing
  complement <- setdiff(seq_len(ncol(obj)), subset)
  miss <- is.na(obj)
  cmiss <- colSums(miss)
  if (any(cmiss / nrow(obj) == 1)) {
    stop("Col(s) with all missing detected. Remove before proceed")
  }
  # Partition
  knn_imp_cols <- (cmiss / nrow(obj)) < colmax
  pre_imp_cols <- obj[, knn_imp_cols, drop = FALSE]
  pre_imp_miss <- miss[, knn_imp_cols, drop = FALSE]
  pre_imp_cmiss <- cmiss[knn_imp_cols]
  knn_indices <- which(knn_imp_cols)
  complement_knn <- intersect(complement, knn_indices)
  pos_complement <- match(complement_knn, knn_indices)
  # Set all values outside of subset to be zero cmiss. This will make
  # `impute_knn_brute`/`impute_knn_mlpack` skip these columns
  pre_imp_cmiss[pos_complement] <- 0L
  # Impute
  if (is.null(tree)) {
    imputed_values <- impute_knn_brute(
      obj = pre_imp_cols,
      miss = pre_imp_miss,
      k = k,
      n_col_miss = pre_imp_cmiss,
      method = switch(method,
        "euclidean" = 0L,
        "manhattan" = 1L
      ),
      dist_pow = dist_pow,
      cores = cores
    )
  } else {
    imputed_values <- impute_knn_mlpack(
      # Has to pre-fill with col means
      obj = mean_imp_col(pre_imp_cols),
      miss = pre_imp_miss,
      k = k,
      n_col_miss = pre_imp_cmiss,
      method = switch(method,
        "euclidean" = 0L,
        "manhattan" = 1L
      ),
      tree = tree,
      dist_pow = dist_pow,
      cores = cores
    )
  }
  # Convert NaN values back to NA
  imputed_values[is.nan(imputed_values)] <- NA_real_
  # Map column indices from pre_imp_cols to original matrix columns and create
  # the index matrix for direct assignment
  imp_indices <- cbind(imputed_values[, 1], knn_indices[imputed_values[, 2]])
  # Impute the object
  obj[imp_indices] <- imputed_values[, 3] # 3rd column is the first imputation

  if (post_imp) {
    if (anyNA(obj[, subset, drop = FALSE])) {
      na_indices <- which(is.na(obj[, subset, drop = FALSE]), arr.ind = TRUE)
      sub_means <- colMeans(obj[, subset, drop = FALSE], na.rm = TRUE)
      i_vec <- na_indices[, 1]
      jj_vec <- na_indices[, 2]
      j_vec <- subset[jj_vec]
      obj[cbind(i_vec, j_vec)] <- sub_means[jj_vec]
    }
  }

  return(obj)
}
