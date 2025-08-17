#' Column Mean Imputation
#'
#' @description
#' Imputes missing values in a matrix by replacing them with the mean of their
#' respective columns.
#'
#' @details
#' This function calculates the mean for each column excluding missing values and
#' replaces all missing values in that column with the computed mean. If a column
#' consists entirely of missing values, those values will remain as `NA` since
#' the mean cannot be calculated.
#'
#' The `subset` parameter allows for selective imputation of only specific columns.
#'
#' @param obj A numeric matrix with samples in rows and features in columns.
#' @param subset Character vector of column names or integer vector of column indices
#'   specifying the subset of columns to perform imputation. If `NULL` (default),
#'   all columns are processed.
#'
#' @return A numeric matrix of the same dimensions as `obj` with missing values
#'   in the specified columns replaced by column means.
#'
#' @examples
#' # Create example matrix with missing values
#' mat <- matrix(c(1, 2, NA, 4, 5, 6, NA, 8, 9), nrow = 3)
#' colnames(mat) <- c("A", "B", "C")
#' mat
#'
#' # Impute missing values with column means
#' imputed_mat <- mean_impute_col(mat)
#' imputed_mat
#'
#' # Impute only specific columns by name
#' imputed_subset <- mean_impute_col(mat, subset = c("A", "C"))
#' imputed_subset
#'
#' # Impute only specific columns by index
#' imputed_idx <- mean_impute_col(mat, subset = c(1, 3))
#' imputed_idx
#'
#' # Example with real data
#' data(khanmiss1)
#'
#' # Transpose since khanmiss1 has genes in rows
#' khanmiss_t <- t(khanmiss1)
#' imputed_khan <- mean_impute_col(khanmiss_t)
#' sum(is.na(khanmiss_t)) # Original missing values
#' sum(is.na(imputed_khan)) # After imputation
#'
#' @export
mean_impute_col <- function(obj, subset = NULL) {
  checkmate::assert_matrix(obj, mode = "numeric", .var.name = "obj")
  checkmate::assert(
    checkmate::check_character(subset, min.len = 0, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
    checkmate::check_integerish(subset, lower = 1, upper = ncol(obj), min.len = 0, any.missing = FALSE, null.ok = TRUE, unique = TRUE),
    combine = "or",
    .var.name = "subset"
  )

  if (!is.null(subset)) {
    if (length(subset) == 0) {
      return(obj)
    }
    if (is.character(subset)) {
      if (is.null(colnames(obj))) {
        stop("`subset` contains characters but `obj` doesn't have column names")
      }
      matched <- match(subset, colnames(obj), nomatch = NA)
      if (any(is.na(matched))) {
        warning("Some subset names not found in column names and will be ignored")
      }
      subset <- matched[!is.na(matched)]
      if (length(subset) == 0) {
        warning("No valid column names found in subset")
        return(obj)
      }
    }
    obj_subset <- obj[, subset, drop = FALSE]
    na_indices <- which(is.na(obj_subset), arr.ind = TRUE)
    if (nrow(na_indices) > 0) {
      column_means <- colMeans(obj_subset, na.rm = TRUE)
      obj_subset[na_indices] <- column_means[na_indices[, 2]]
      obj[, subset] <- obj_subset
    }
  } else {
    na_indices <- which(is.na(obj), arr.ind = TRUE)
    if (nrow(na_indices) > 0) {
      column_means <- colMeans(obj, na.rm = TRUE)
      obj[na_indices] <- column_means[na_indices[, 2]]
    }
  }
  return(obj)
}

# Handle the creation of output files and handle the overwriting
check_result_list <- function(output, n_imp, overwrite) {
  # Return NULL if no output specified. Same as !file_backed
  if (is.null(output) || length(output) == 0) {
    return(list(backingpath = NULL, backfiles = NULL, descfiles = NULL))
  }
  # Prepare file paths
  output <- fs::as_fs_path(output)
  backingpath <- fs::path_dir(output)
  # Validate path
  if (!all(fs::file_access(backingpath, mode = c("exists", "read", "write")))) {
    stop("Provided path to `output` doesn't exist or is not readable/writable")
  }
  base <- fs::path_file(output)
  # Generate file names for all imputation iterations
  imp_indicies <- seq_len(n_imp)
  suffix <- paste0("_imp", imp_indicies)
  backfiles <- paste0(base, suffix, ".bin")
  descfiles <- paste0(base, suffix, ".desc")
  # Check existing files and handle overwrite
  for (i in imp_indicies) {
    files_to_check <- fs::path(backingpath, c(backfiles[i], descfiles[i]))
    if (any(fs::file_exists(files_to_check))) {
      if (!overwrite) {
        stop(
          "Output files already exist. Set `overwrite = TRUE` to overwrite them."
        )
      }
      # Delete existing files
      unlink(files_to_check, force = TRUE)
      # if (.Platform$OS.type == "windows") {
      #   Sys.sleep(0.1)
      # }
      # invisible(gc(verbose = FALSE, full = TRUE))
      # Verify deletion
      if (any(fs::file_exists(files_to_check))) {
        stop("Failed to delete output files.")
      }
    }
  }
  return(list(
    backingpath = backingpath,
    backfiles = backfiles,
    descfiles = descfiles
  ))
}

# Helper function for knn_imp to create bigmatrix or copy for a single imputation.
# No need to pass to SlideKnn because knn_imp is always in memory there.
create_result_list <- function(
    data_to_copy,
    file_backed,
    n_imp,
    backfiles = NULL,
    descfiles = NULL,
    backingpath = NULL) {
  result_list <- vector("list", n_imp)
  if (file_backed) {
    for (i in seq_len(n_imp)) {
      bigmat <- bigmemory::big.matrix(
        nrow = nrow(data_to_copy),
        ncol = ncol(data_to_copy),
        type = "double",
        init = NA,
        backingpath = backingpath,
        backingfile = backfiles[i],
        descriptorfile = descfiles[i]
      )
      bigmat[, ] <- data_to_copy
      bigmemory::flush(bigmat)
      result_list[[i]] <- bigmat
    }
  } else {
    # For in-memory, create n_imp copies
    for (i in seq_len(n_imp)) {
      result_list[[i]] <- data_to_copy[, , drop = FALSE]
    }
  }
  return(result_list)
}

#' K-Nearest Neighbor Imputation for Missing Values
#'
#' @description
#' Imputes missing values in numeric matrices using the k-nearest neighbors algorithm
#' with a two-stage approach: k-NN imputation for columns with missingness below a
#' threshold, followed by optional mean imputation for remaining missing values.
#'
#' @details
#' This function performs **column-wise** distance calculations, which is particularly
#' important for understanding its application domain. Unlike typical data where
#' samples are in columns and features in rows, this function is optimized for
#' intensive longitudinal data and epi-genomics where:
#' - **Samples are in rows** and **features are in columns**
#' - Spatially or temporally adjacent features carry mutual information
#' - Imputation leverages nearby features within the same sample
#'
#' ## Weighting and Tree Methods
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
#' - **Tree methods**: Only use when imputation runtime becomes prohibitive and missingness is low
#' - **KDTree** (`tree = "kd"`): Suitable for >5000 columns with <20% missingness
#' - **BallTree** (`tree = "ball"`): For high-dimensional data or Manhattan distance
#' - **File backing**: Use `output` parameter for `n_imp > 1` to avoid memory issues
#' - **Subset imputation**: Consider `subset` parameter for efficiency when only specific columns need imputation
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#'   Must have at least 1 row and 2 columns.
#' @inheritParams SlideKnn
#' @param subset Character vector of column names or integer vector of column
#'   indices specifying which columns to impute.
#' @param tree Character. k-NN search method: `NULL` (brute-force), `"kd"` (KDTree),
#'   or `"ball"` (BallTree). Tree methods may introduce bias with high missingness.
#'   Default: `NULL`.
#' @param ... Currently not implemented.
#'
#' @return A list of length `n_imp` containing numeric matrices or [bigmemory::big.matrix()]
#'   objects (if `output` is specified) with the same dimensions as `obj`. Missing values
#'   are imputed using k-NN for columns with missingness below `colmax`, and mean
#'   imputation for remaining missing values if `post_imp = TRUE`.
#'
#'   The list has class `"KnnImpList"` with attributes:
#'   \itemize{
#'     \item `rownames`: Original row names from input matrix
#'     \item `colnames`: Original column names from input matrix
#'     \item `subset`: Column indices that were processed for imputation
#'     \item `ncol`: Number of columns in original matrix
#'   }
#'
#'   Each list element represents an independent imputation for uncertainty quantification.
#'
#' @note
#' **File-backed Storage**: For file-backed results using `output`, set
#' `options(bigmemory.allow.dimnames = TRUE)` before calling this function to preserve
#' dimnames, otherwise they can be manually restored from the original matrix with
#' [restore_dimnames()].
#'
#' **Memory Considerations**: When `n_imp > 1` for large data, use `subset` to
#' specify only required columns and provide `output` for file-backed storage to
#' avoid memory constraints.
#'
#' @references
#' Robert Tibshirani, Trevor Hastie, Balasubramanian Narasimhan, and Gilbert Chu (2002).
#' Diagnosis of multiple cancer types by shrunken centroids of gene expression
#' PNAS 99: 6567-6572. Available at www.pnas.org
#'
#' @seealso [SlideKnn()], [mean_impute_col()], [bigmemory::big.matrix()], [restore_dimnames()]
#'
#' @examples
#' # Quick start example
#' data(khanmiss1)
#' sum(is.na(khanmiss1))
#'
#' # Basic k-NN imputation (khanmiss1 has genes in rows, so transpose)
#' t_khanmiss1 <- t(khanmiss1)
#' result <- knn_imp(t_khanmiss1, k = 5)[[1]]
#' sum(is.na(result)) # Should be 0
#'
#' # Using weighted imputation with custom distance power multiple cores
#' imputed_weighted <- knn_imp(
#'   obj = t_khanmiss1,
#'   k = 5,
#'   cores = 4,
#'   weighted = TRUE,
#'   dist_pow = 2,
#'   method = "euclidean"
#' )
#'
#' # Preview
#' imputed_weighted
#'
#' # Access the imputations just like a list
#' imputed_weighted[[1]][1:5, 1:5]
#'
#' # Using only a subset of columns. Massive time saver
#' imputed_subset <- knn_imp(
#'   obj = t_khanmiss1,
#'   k = 5,
#'   cores = 4,
#'   weighted = TRUE,
#'   dist_pow = 2,
#'   subset = c("g189", "g299", "g361"),
#'   method = "euclidean"
#' )
#'
#' # PMM imputation with 5 imputations and 5 donors. Increase `n_pmm` in real
#' # data to ensure realistic uncertainty quantification.
#' imputed_pmm <- knn_imp(
#'   obj = t_khanmiss1,
#'   k = 5,
#'   n_imp = 5,
#'   n_pmm = 5,
#'   output = withr::local_tempfile(),
#'   # Enables results as list of bigmatrix
#'   seed = 123
#' )
#' # The default of the bigmatrix package is to not allow dimnames so the output
#' # may lost the dimnames
#' imputed_pmm[[1]][, ][1:5, 1:5]
#'
#' if (interactive()) {
#'   # But this can be reassigned after setting the options to be TRUE like so
#'   options(bigmemory.allow.dimnames = TRUE)
#'   restore_dimnames(imputed_pmm)
#' }
#'
#' length(imputed_pmm) # Returns 5 imputed data
#'
#' # Bootstrap imputation for uncertainty injection
#' imputed_boot <- knn_imp(
#'   obj = t_khanmiss1,
#'   k = 3,
#'   n_imp = 5,
#'   n_pmm = 0,
#'   # n_pmm = 0 enables bootstrapping nearest neighbors
#'   seed = 123
#' )
#'
#' length(imputed_boot) # Returns 5 imputed data
#'
#' @export
knn_imp <- function(
    obj,
    k,
    colmax = 0.9,
    rowmax = 0.9,
    method = c("euclidean", "manhattan"),
    cores = 1,
    post_imp = TRUE,
    subset = NULL,
    weighted = FALSE,
    dist_pow = 1,
    tree = NULL,
    n_imp = 1,
    n_pmm = -1,
    seed = 42,
    output = NULL,
    overwrite = FALSE,
    ...) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, null.ok = FALSE, .var.name = "obj")
  checkmate::assert_int(k, lower = 1, upper = ncol(obj) - 1, .var.name = "k")
  checkmate::assert_int(cores, lower = 1, .var.name = "cores")
  checkmate::assert_number(colmax, lower = 0, upper = 1, .var.name = "colmax")
  checkmate::assert_number(rowmax, lower = 0, upper = 1, .var.name = "rowmax")
  checkmate::assert_flag(post_imp, null.ok = FALSE, .var.name = "post_imp")
  checkmate::assert_int(n_imp, lower = 1, .var.name = "n_imp")
  checkmate::assert_int(n_pmm, lower = -1, upper = nrow(obj), .var.name = "n_pmm")
  if (n_imp > 1 && n_pmm == -1) {
    warning("n_imp > 1 requires n_pmm >= 0. Setting n_imp = 1.")
    n_imp <- 1
  }
  checkmate::assert_int(seed, lower = 0, .var.name = "seed", null.ok = FALSE)
  checkmate::assert(
    checkmate::check_character(subset, min.len = 0, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
    checkmate::check_integerish(subset, lower = 1, upper = ncol(obj), min.len = 0, any.missing = FALSE, null.ok = TRUE, unique = TRUE),
    combine = "or",
    .var.name = "subset"
  )
  checkmate::assert_flag(weighted)
  if (weighted) {
    stopifnot(length(dist_pow) == 1, dist_pow > 0, !is.infinite(dist_pow))
  }
  checkmate::assert_choice(tree, choices = c("kd", "ball"), null.ok = TRUE, .var.name = "tree")
  if (!is.null(output)) {
    checkmate::assert_path_for_output(output, .var.name = "output")
  }
  checkmate::assert_flag(overwrite, null.ok = FALSE, .var.name = "overwrite")
  set.seed(seed)
  # Determine if output should be file-backed
  file_backed <- !is.null(output)
  # Handle file-backed output preparation
  paths <- check_result_list(output = output, n_imp = n_imp, overwrite = overwrite)
  backingpath <- paths$backingpath
  backfiles <- paths$backfiles
  descfiles <- paths$descfiles
  # Handle subset conversion and validation
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
  # Calculate subset names for metadata
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = TRUE)
  can_store_dimnames <- isTRUE(getOption("bigmemory.allow.dimnames")) || !file_backed
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)
  rn <- rownames(obj)
  cn <- colnames(obj)
  # Early return for empty subset or no missing data
  if (length(subset) == 0 || !anyNA(obj[, subset, drop = FALSE])) {
    if (length(subset) > 0) {
      message("No missing data in subset columns")
    }
    # Create result with proper dimnames handling
    result_list <- create_result_list(
      data_to_copy = obj,
      file_backed = file_backed,
      n_imp = n_imp,
      backfiles = backfiles,
      descfiles = descfiles,
      backingpath = backingpath
    )
    # Add metadata
    attr(result_list, "rownames") <- rn
    attr(result_list, "colnames") <- cn
    attr(result_list, "subset") <- subset
    attr(result_list, "ncol") <- ncol(obj)
    class(result_list) <- c("KnnImpList", class(result_list))
    return(result_list)
  }
  # Calculate complement for further processing
  complement <- setdiff(seq_len(ncol(obj)), subset)
  miss <- is.na(obj)
  rmiss <- rowSums(miss) / ncol(obj)
  if (any(rmiss >= rowmax)) {
    stop("Row(s) missing exceeded rowmax. Remove row(s) with too high NA %")
  }
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
  if (any(rowSums(pre_imp_miss) / ncol(pre_imp_cols) == 1)) {
    stop("Row(s) missing exceeded rowmax. Remove row(s) with too high NA %")
  }
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
      weighted = weighted,
      dist_pow = dist_pow,
      n_imp = n_imp,
      n_pmm = n_pmm,
      seed = seed,
      cores = cores
    )
  } else {
    imputed_values <- impute_knn_mlpack(
      # Has to pre-fill with colMeans
      obj = mean_impute_col(pre_imp_cols),
      miss = pre_imp_miss,
      k = k,
      n_col_miss = pre_imp_cmiss,
      method = switch(method,
        "euclidean" = 0L,
        "manhattan" = 1L
      ),
      tree = tree,
      weighted = weighted,
      dist_pow = dist_pow,
      n_imp = n_imp,
      n_pmm = n_pmm,
      seed = seed,
      cores = cores
    )
  }
  # Convert NaN values back to NA
  imputed_values[is.nan(imputed_values)] <- NA_real_
  # Map column indices from pre_imp_cols to original matrix columns and create
  # the index matrix for direct assignment
  imp_indices <- cbind(imputed_values[, 1], knn_indices[imputed_values[, 2]])
  # Initialize the result list
  result_list <- create_result_list(
    data_to_copy = obj,
    file_backed = file_backed,
    n_imp = n_imp,
    backfiles = backfiles,
    descfiles = descfiles,
    backingpath = backingpath
  )
  # Fill in imputed values
  for (i in seq_len(n_imp)) {
    result_list[[i]][imp_indices] <- imputed_values[, i + 2]
  }
  # Flush file-backed matrices to disk and post-impute
  if (file_backed) {
    for (i in seq_len(n_imp)) {
      bigmemory::flush(result_list[[i]])
    }
    if (post_imp) {
      for (i in seq_len(n_imp)) {
        bigmem_impute_colmeans(result_list[[i]]@address, subset, cores = cores)
      }
    }
  } else {
    if (post_imp) {
      for (i in seq_len(n_imp)) {
        subset_data <- result_list[[i]][, subset, drop = FALSE]
        if (anyNA(subset_data)) {
          na_indices <- which(is.na(subset_data), arr.ind = TRUE)
          sub_means <- colMeans(subset_data, na.rm = TRUE)
          i_vec <- na_indices[, 1]
          jj_vec <- na_indices[, 2]
          j_vec <- subset[jj_vec]
          result_list[[i]][cbind(i_vec, j_vec)] <- sub_means[jj_vec]
        }
      }
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
  attr(result_list, "subset") <- subset
  attr(result_list, "ncol") <- ncol(obj)
  class(result_list) <- c("KnnImpList", class(result_list))

  return(result_list)
}


#' k-NN Imputation Wrapper
#'
#' @description
#' A wrapper function for k-NN imputation that filters rows based on missing value
#' proportions before imputing. Rows with missing proportions exceeding `rowmax`
#' are not imputed in this step to avoid throwing exceptions in [knn_imp()].
#'
#' @details
#' This function serves as an intermediate layer between the sliding window algorithm
#' and the core k-NN imputation function. It ensures that only rows with acceptable
#' missing value proportions are processed, preventing errors and improving
#' computational efficiency.
#'
#' @inheritParams SlideKnn
#' @param knn_imp Function object for k-NN imputation. Required for function crating
#'   in parallel environments.
#' @param ... Additional arguments passed to `knn_imp`.
#'
#' @return A list of imputed matrices, with only qualifying rows imputed via k-NN;
#'   others remain unchanged. Length of list equals `n_imp`.
#'
#' @seealso [knn_imp()], [SlideKnn()]
#'
#' @keywords internal
#' @noRd
impute_knn <- function(
    obj,
    k,
    rowmax,
    colmax,
    cores,
    method,
    post_imp,
    subset,
    weighted,
    dist_pow,
    tree,
    knn_imp,
    n_imp,
    n_pmm,
    seed,
    ...) {
  na_mat <- is.na(obj)
  # Determine 'good_rows': rows where the proportion of NAs is less than 'rowmax'.
  # Rows with too many NAs (as defined by rowmax) are excluded from imputation.
  good_rows <- rowSums(na_mat) / ncol(na_mat) < rowmax
  # If no rows meet the criteria for imputation, return list with original object.
  if (sum(good_rows) == 0) {
    return(replicate(n_imp, obj, simplify = FALSE))
  }
  # Call knn_imp on good rows - returns list of imputed objects
  imputed_good_list <- knn_imp(
    obj = obj[good_rows, , drop = FALSE],
    k = k,
    rowmax = rowmax,
    colmax = colmax,
    method = method,
    cores = cores,
    post_imp = post_imp,
    subset = subset,
    weighted = weighted,
    dist_pow = dist_pow,
    tree = tree,
    n_imp = n_imp,
    n_pmm = n_pmm,
    seed = seed,
    ...
  )
  # Process each imputed object in the list
  result_list <- vector("list", length(imputed_good_list))

  for (i in seq_along(imputed_good_list)) {
    # Initialize the result matrix with the original object's values.
    result <- obj
    # Place the imputed values from this imputation iteration back into the
    # corresponding 'good_rows' in the result matrix.
    result[good_rows, ] <- imputed_good_list[[i]]
    result_list[[i]] <- result
  }

  return(result_list)
}

# Given a start and end vectors. Give the counts_vec that counts the number
# of times that position is covered and also the regions where counts are > 1
find_overlap_regions <- function(start, end) {
  max_pos <- max(end)
  delta <- rep(0, max_pos + 1)
  delta[start] <- delta[start] + 1
  valid_ends <- end + 1 <= max_pos + 1
  delta[end[valid_ends] + 1] <- delta[end[valid_ends] + 1] - 1
  counts_vec <- cumsum(delta)[seq_len(max_pos)]
  storage.mode(counts_vec) <- "integer"
  # extract regions where counts > 1
  overlaps <- counts_vec > 1
  rle_over <- rle(overlaps)
  # compute start and end positions of the TRUE runs
  ends <- cumsum(rle_over$lengths)
  starts <- ends - rle_over$lengths + 1
  # return the region and counts_vec
  return(
    list(
      region = cbind(
        start = starts[rle_over$values],
        end = ends[rle_over$values]
      ),
      counts_vec = counts_vec
    )
  )
}

#' Sliding Window k-NN Imputation for Large data
#'
#' @description
#' Performs k-nearest neighbor imputation on large numeric matrices using a sliding
#' window approach. The matrix is divided into overlapping windows to handle missing
#' values while maintaining computational efficiency and preserving local data structures.
#'
#' @details
#' The sliding window approach divides the input matrix into smaller, overlapping
#' segments and applies k-NN imputation to each window independently. This strategy
#' is particularly advantageous for large data where applying k-NN imputation
#' to the entire matrix would be computationally prohibitive or exceed memory limits.
#'
#' The algorithm maintains local data relationships by using overlapping regions
#' between consecutive windows. Values in overlapping areas are averaged across
#' windows to produce the final imputed result. This approach assumes that features
#' (columns) are ordered meaningfully (e.g., by genomic position, time series, etc.).
#'
#' For the underlying k-NN implementation details, see [knn_imp()].
#'
#' @note
#' **Multiple Imputation**: Setting `n_imp` > 1 requires `n_pmm` >= 0 and should
#' be used with the `subset` parameter and `output` for file-backed storage to
#' manage memory efficiently.
#'
#' **Big Matrix Compatibility**: When using [bigmemory::big.matrix()] objects,
#' dimnames may not be preserved in the output due to bigmemory package defaults.
#' To retain dimnames, set `options(bigmemory.allow.dimnames = TRUE)` and manually
#' reassign dimnames to the output using the original object's dimnames.
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#'   Features should be meaningfully ordered (e.g., by genomic position or time).
#'   Accepts: numeric matrix, [bigmemory::big.matrix()], or path to big.matrix
#'   description file.
#' @param n_feat Integer. Number of features (columns) per sliding window.
#'   Must be between 2 and `ncol(obj)`.
#' @param subset Character vector of column names or integer vector of column
#'   indices specifying which columns to impute. If `NULL` (default), all columns
#'   are processed. Required when `n_imp` > 1 and `n_pmm` >= 0.
#' @param n_overlap Integer. Number of overlapping features between consecutive
#'   windows. Must be between 0 and `n_feat - 1`. Default: 10.
#' @param k Integer. Number of nearest neighbors for imputation. Must be between
#'   1 and `n_feat - 1`. Default: 10.
#' @param rowmax Numeric. Maximum allowable proportion of missing values per row
#'   (0-1). Function stops with error if exceeded. Default: 0.9.
#' @param colmax Numeric. Threshold proportion of missing values per column (0-1).
#'   Columns exceeding this threshold are imputed using column means instead of
#'   k-NN when `post_imp = TRUE`. Default: 0.9.
#' @param cores Integer. Number of CPU cores for parallel distance computation.
#'   Requires [mirai::daemons()] setup for `cores` > 1. Default: 1.
#' @param method Character. Distance metric for k-NN: `"euclidean"` or `"manhattan"`.
#'   Default: `"euclidean"`.
#' @param tree Character. k-NN search method: `NULL` (brute-force), `"kd"` (KDTree),
#'   or `"ball"` (BallTree). Tree methods use mlpack implementation but may be
#'   biased with high missing value percentages. Default: `NULL`.
#' @param post_imp Logical. Whether to impute remaining missing values with
#'   column means after k-NN imputation. Default: `TRUE`.
#' @param weighted Logical. Whether to use distance-weighted mean for imputation.
#'   If `FALSE`, uses simple mean of k nearest neighbors. Default: `FALSE`.
#' @param dist_pow Numeric. Positive value controlling distance penalty in weighted
#'   imputation. Values < 1 apply softer penalty, 1 is linear, > 1 is harsher.
#'   Only used when `weighted = TRUE`. Default: 1.
#' @param n_imp Integer. Number of multiple imputations to perform. Automatically
#'   set to 1 if `n_pmm = -1`. Default: 1.
#' @param n_pmm Integer. Multiple imputation method control:
#'   \itemize{
#'     \item `-1`: Deterministic single imputation (default)
#'     \item `> 0`: Predictive Mean Matching using `n_pmm` closest donors
#'     (recommended for MI. `8` is a good starting point.).
#'     \item `0`: Bootstrap resampling from k-nearest neighbors
#'   }
#' @param seed Integer. Random seed for reproducible results in stochastic imputation
#'   methods. Default: 42.
#' @param .progress Logical. Whether to display progress messages during imputation.
#'   Default: `FALSE`.
#' @param output Character. File path stem for saving file-backed big.matrix results
#'   (format: `"path/stem"`). If `NULL`, results are stored in memory. Recommended
#'   for large data and multiple imputations.
#' @param overwrite Logical. Whether to overwrite existing files at `output` path.
#'   Default: `FALSE`.
#'
#' @return A list of [bigmemory::big.matrix()] objects with length `n_imp`. Each
#'   matrix has `nrow(obj)` rows and `length(subset)` columns with missing values
#'   imputed. The list has class `"SlideKnnList"` with additional attributes:
#'   \itemize{
#'     \item `rownames`: Original row names
#'     \item `colnames`: Column names for imputed subset
#'     \item `ncol`: Number of columns in original matrix
#'     \item `subset`: Column indices that were imputed
#'   }
#'
#' @seealso [knn_imp()], [mean_impute_col()], [bigmemory::big.matrix()], [restore_dimnames()]
#'
#' @examples
#' # Generate sample data with missing values with 20 samples and 100 columns
#' # where the column order is sorted (e.g., by genomic position or time)
#'
#' set.seed(1234)
#' beta_matrix <- t(sim_mat(100, 20)$input)
#'
#' # ========================================
#' # Basic Sliding k-NN Imputation
#' # ========================================
#'
#' # Simple imputation with default parameters
#' imputed_basic <- SlideKnn(
#'   beta_matrix,
#'   k = 5,
#'   n_feat = 50,
#'   n_overlap = 10
#' )
#'
#' # Preview
#' imputed_basic
#'
#' # Access the result like a list
#' imputed_basic[[1]][1:5, 1:5]
#'
#' # ========================================
#' # Big Matrix Usage (Memory Efficient)
#' # ========================================
#'
#' # Convert to big.matrix for better memory management
#' big_data <- bigmemory::as.big.matrix(beta_matrix, type = "double")
#'
#' # output is NULL, the backend is still in-memory
#' imputed_big_in_memory <- SlideKnn(
#'   big_data,
#'   k = 5,
#'   n_feat = 50,
#'   output = NULL
#' )
#'
#' # output is now set to a location to save. The backend is now filebacked matrix
#' # to handle massive data
#' imputed_big <- SlideKnn(
#'   big_data,
#'   k = 5,
#'   n_feat = 50,
#'   output = withr::local_tempfile()
#' )
#'
#' # Because strip_dimnames is `FALSE` and the default of the bigmatrix package
#' # is to not allow dimnames, the output may lost the dimnames. We strip it here to
#' # demonstrate
#'
#' # But this can be reassigned after setting the options to be TRUE like so
#' if (interactive()) {
#'   on.exit(options(bigmemory.allow.dimnames = getOption("bigmemory.allow.dimnames")), add = TRUE)
#'   options(bigmemory.allow.dimnames = TRUE)
#'   rownames(imputed_big[[1]]) <- NULL
#'   colnames(imputed_big[[1]]) <- NULL
#'   # Names are now strip
#'   print(imputed_big)
#'   # Restore dimnames
#'   restore_dimnames(imputed_big)
#'   # Now restored
#'   print(imputed_big)
#' }
#'
#' # ========================================
#' # Multiple Imputation
#' # ========================================
#'
#' # Predictive Mean Matching (PMM) - recommended
#' imputed_pmm <- SlideKnn(
#'   beta_matrix,
#'   k = 8,
#'   n_feat = 60,
#'   n_imp = 3, # 3 imputations
#'   n_pmm = 5, # 5 donors for PMM
#'   subset = c("feat1", "feat2"),
#'   output = withr::local_tempfile(),
#'   .progress = TRUE
#' )
#'
#' imputed_pmm
#'
#' # ========================================
#' # Parallel Processing
#' # ========================================
#' @examplesIf interactive()
#' # Set up parallel processing
#' library(mirai)
#' daemons(4) # Use 4 cores
#'
#' # Enable for multi-core efficiency
#' options(bigmemory.allow.dimnames = TRUE)
#'
#' # Parallel imputation
#' imputed_parallel <- SlideKnn(
#'   beta_matrix,
#'   k = 5,
#'   n_feat = 50,
#'   cores = 4
#' )
#'
#' imputed_parallel
#'
#' # Clean up
#' daemons(0)
#'
#' @export
SlideKnn <- function(
    obj,
    n_feat,
    subset = NULL,
    n_overlap = 10,
    k = 10,
    rowmax = 0.9,
    colmax = 0.9,
    cores = 1,
    method = c("euclidean", "manhattan"),
    tree = NULL,
    post_imp = TRUE,
    weighted = FALSE,
    dist_pow = 1,
    n_imp = 1,
    n_pmm = -1,
    seed = 42,
    .progress = FALSE,
    output = NULL,
    overwrite = FALSE) {
  # Pre-conditioning ----
  method <- match.arg(method)
  checkmate::assert_number(rowmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_int(cores, lower = 1, null.ok = FALSE, .var.name = "cores")
  checkmate::assert_int(n_imp, lower = -1, null.ok = FALSE, .var.name = "n_imp")
  checkmate::assert_int(seed, null.ok = FALSE, .var.name = "seed")
  if (cores > 1) {
    tryCatch(
      mirai::require_daemons(),
      error = function(e) {
        stop(sprintf("%d cores requested, but no mirai daemon is setup. Call mirai::daemons(%d) to set up the parallelization", cores, cores))
      }
    )
    fn <- purrr::in_parallel
  } else {
    fn <- carrier::crate
  }
  checkmate::assert_flag(post_imp, .var.name = "post_imp", null.ok = FALSE)
  checkmate::assert_flag(.progress, .var.name = ".progress", null.ok = FALSE)
  if (!is.null(output)) {
    checkmate::assert_path_for_output(output, .var.name = "output")
  }
  checkmate::assert(
    checkmate::check_character(subset, min.len = 1, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
    checkmate::check_integerish(
      subset,
      lower = 1, upper = ncol(obj), min.len = 1, any.missing = FALSE, null.ok = TRUE, unique = TRUE
    ),
    .var.name = "subset"
  )

  # Enforce subset requirement for multiple imputations
  if (n_imp > 1 && n_pmm >= 0 && is.null(subset)) {
    stop(
      "When n_imp > 1 and n_pmm >= 0, subset must be explicitly specified. ",
      "Please provide a subset of columns to impute."
    )
  }
  if (n_imp > 1 && n_pmm == 0 && weighted) {
    warning("If bootstrapping nearest neighbors, weighted will be forced to FALSE")
    weighted <- FALSE
  }

  # Handle input matrix
  if (!bigmemory::is.big.matrix(obj)) {
    if (is.character(obj)) {
      obj <- bigmemory::attach.big.matrix(obj)
    } else if (is.matrix(obj) && is.numeric(obj)) {
      obj <- bigmemory::as.big.matrix(obj, type = "double")
    } else {
      stop("`obj` has to be a numeric matrix or a (descriptor of a) double bigmemory::big.matrix.")
    }
  }
  checkmate::assert_int(n_feat, lower = 2, upper = ncol(obj), null.ok = FALSE, .var.name = "n_feat")
  checkmate::assert_int(n_overlap, lower = 0, upper = n_feat - 1, null.ok = FALSE, .var.name = "n_overlap")
  checkmate::assert_int(k, lower = 1, upper = n_feat - 1, null.ok = FALSE, .var.name = "k")

  if (n_imp > 1 && n_pmm == -1) {
    warning("n_imp > 1 requires n_pmm >= 0. Setting n_imp = 1.")
    n_imp <- 1
  }

  # Warning for multiple imputations without file backing
  file_backed <- !is.null(output)
  if (n_imp > 1 && !file_backed) {
    warning("n_imp > 1 should be used with output to file-backed big.matrix. Provide `output` parameter to enable file backing.")
  }

  rn <- rownames(obj)
  cn <- colnames(obj)
  ## subset ----
  if (!is.null(subset)) {
    if (is.character(subset)) {
      stopifnot("`subset` are characters but `obj` doesn't have colnames" = !is.null(cn))
      matched <- match(subset, cn, nomatch = NA)
      # subset converted to index
      subset <- matched[!is.na(matched)]
      # sort to avoid trouble
      subset <- sort(subset)
    }
    if (length(subset) == 0) {
      stop("`subset` is not found in `colnames(obj)`")
    }
  } else {
    # If subset is NULL, set it to all columns. Fine for n_pmm == -1
    subset <- seq_len(ncol(obj))
  }

  # Windowing Logic ----
  # [Keep all the windowing logic as is - lines for calculating start, end, etc.]
  idx <- 1 # R index at 1
  max_step <- ceiling((ncol(obj) - idx) / (n_feat - n_overlap))
  step <- 0:max_step
  start <- idx + (step * n_feat) - (step * n_overlap)
  end <- start + n_feat - 1
  # Edge case, the end might overshoot the number of rows of the obj
  n_overshoot <- sum(end > ncol(obj))
  # In which case trim off the runs that overshoot
  corrected_length <- length(end) - n_overshoot
  start <- start[1:corrected_length]
  end <- end[1:corrected_length]
  # And make the last window extra wide to cover the full end
  end[corrected_length] <- ncol(obj)
  # Calculate where the subset lies, offset by the start index
  subset_list <- lapply(seq_along(start), function(i) {
    first <- findInterval(start[i] - 1, subset) + 1
    last <- findInterval(end[i], subset)
    if (first <= last) {
      subset[first:last] - start[i] + 1
    } else {
      integer(0)
    }
  })
  # Then, we calculate the offsets needed to subset the intermediate matrix
  width <- end - start + 1
  offset_start <- c(1, cumsum(width)[-length(width)] + 1)
  offset_end <- cumsum(width)

  # Sliding Imputation ----
  ## Init ----
  nr <- nrow(obj)
  nc <- ncol(obj)
  # Getting the backingpath and backfiles/descfiles. NULL if output is NULL.
  temp_dir <- if (file_backed) {
    checkmate::assert_flag(overwrite, .var.name = "overwrite", null.ok = FALSE)
    withr::local_tempdir(pattern = paste0("SlideKnn_", Sys.getpid()))
  } else {
    NULL
  }
  # Create descriptors to pass around on workers
  obj_desc <- bigmemory::describe(obj)
  # Check if we can store dimnames
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = TRUE)
  can_store_dimnames <- isTRUE(getOption("bigmemory.allow.dimnames"))
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)
  if (can_store_dimnames) {
    # strip dimnames if can store dimnames to be restored to save memory
    rownames(obj) <- NULL
    colnames(obj) <- NULL
  }
  final_colnames <- if (!is.null(cn)) {
    cn[subset]
  } else {
    as.character(subset)
  }
  # Get final imputed info
  final_imputed_info <- check_result_list(output, n_imp, overwrite)

  # Create all final imputed big.matrices
  final_imputed_list <- lapply(seq_len(n_imp), function(imp_idx) {
    bigmemory::big.matrix(
      nrow = nr,
      ncol = length(subset),
      type = "double",
      init = 0.0,
      backingpath = final_imputed_info$backingpath,
      descriptorfile = final_imputed_info$descfiles[[imp_idx]],
      backingfile = final_imputed_info$backfiles[[imp_idx]]
    )
  })

  # Overlap regions to average over
  overlap <- find_overlap_regions(start, end)

  # Process each imputation sequentially for each n_imp
  for (imp_idx in seq_len(n_imp)) {
    if (.progress) {
      message(sprintf("Processing imputation %d/%d", imp_idx, n_imp))
    }
    # Seed for this imputation
    seed_imp <- seed + (imp_idx - 1)
    # Init/overwrite one intermediate matrix per imputation to avoid race conds
    intermediate_info <- check_result_list(
      fs::path(temp_dir, "intermediate"), 1,
      overwrite = TRUE
    )
    intermediate <- bigmemory::big.matrix(
      nrow = nr,
      ncol = sum(width),
      type = "double",
      init = 0.0,
      backingpath = intermediate_info$backingpath,
      descriptorfile = intermediate_info$descfiles[[1]],
      backingfile = intermediate_info$backfiles[[1]]
    )
    intermediate_desc <- bigmemory::describe(intermediate)

    ## Impute ----
    if (.progress) {
      message("Step 1/3: Imputing")
    }
    purrr::walk(
      seq_along(start),
      fn(
        function(i, ...) {
          window_cols <- start[i]:end[i]
          obj_big <- bigmemory::attach.big.matrix(obj_desc)
          # Get imputation results
          imp_list <- impute_knn(
            # Realize in memory just the window cols
            obj = obj_big[, window_cols, drop = FALSE],
            k = k,
            rowmax = rowmax,
            colmax = colmax,
            cores = 1L, # For Slide k-NN, fix cores = 1
            method = method,
            post_imp = post_imp,
            weighted = weighted,
            dist_pow = dist_pow,
            subset = subset_list[[i]],
            n_imp = 1L, # Fix n_imp to 1
            n_pmm = n_pmm, # Control multiple imputation only through `n_pmm`
            seed = seed_imp,
            knn_imp = knn_imp,
            tree = tree,
            impute_knn_brute = impute_knn_brute,
            impute_knn_mlpack = impute_knn_mlpack,
            mean_impute_col = mean_impute_col,
            check_result_list = check_result_list
          )
          # Fill intermediate matrix
          intermediate_big <- bigmemory::attach.big.matrix(intermediate_desc)
          intermediate_big[, offset_start[i]:offset_end[i]] <- imp_list[[1]]
        },
        impute_knn = impute_knn,
        subset_list = subset_list,
        knn_imp = knn_imp,
        weighted = weighted,
        dist_pow = dist_pow,
        impute_knn_brute = impute_knn_brute,
        mean_impute_col = mean_impute_col,
        check_result_list = check_result_list,
        impute_knn_mlpack = impute_knn_mlpack,
        start = start,
        end = end,
        obj_desc = obj_desc,
        tree = tree,
        intermediate_desc = intermediate_desc,
        k = k,
        rowmax = rowmax,
        colmax = colmax,
        method = method,
        post_imp = post_imp,
        offset_start = offset_start,
        offset_end = offset_end,
        n_pmm = n_pmm,
        seed_imp = seed_imp
      ),
      .progress = .progress
    )
    if (file_backed) {
      bigmemory::flush(intermediate)
    }

    # Init/overwrite the temporary result_imp matrix once per imputation
    result_imp_info <- check_result_list(
      fs::path(temp_dir, "result_imp"), 1,
      overwrite = TRUE
    )
    result_imp <- bigmemory::big.matrix(
      nrow = nr,
      ncol = nc,
      type = "double",
      init = 0.0,
      backingpath = result_imp_info$backingpath,
      descriptorfile = result_imp_info$descfiles[[1]],
      backingfile = result_imp_info$backfiles[[1]]
    )

    ## Averaging ----
    if (.progress) {
      message("Step 2/3: Overlapping")
    }
    # Add the windows from intermediate. Single threaded
    bigmem_add_windows(
      result_imp@address,
      intermediate@address, start, end, offset_start, offset_end
    )
    if (.progress) {
      message("Step 3/3: Averaging")
    }
    # Averaging the overlaps
    bigmem_avg(
      result_imp@address,
      start = overlap$region[, "start"],
      end = overlap$region[, "end"],
      counts_vec = overlap$counts_vec,
      cores = cores
    )
    ## post_imp ----
    if (post_imp) {
      if (.progress) {
        message("Post-imputation")
      }
      bigmem_impute_colmeans(
        result_imp@address,
        col_indices = subset,
        cores = cores
      )
    }
    # Store result
    bigmem_copy(
      final_imputed_list[[imp_idx]]@address,
      result_imp@address,
      col_idx_r = subset,
      cores = cores
    )

    # Set dimnames if possible
    if (can_store_dimnames) {
      rownames(final_imputed_list[[imp_idx]]) <- rn
      colnames(final_imputed_list[[imp_idx]]) <- final_colnames
    }
  }
  # S3 OOPs
  if (can_store_dimnames) {
    rownames(obj) <- rn
    colnames(obj) <- cn
  }
  attr(final_imputed_list, "rownames") <- rn
  attr(final_imputed_list, "colnames") <- final_colnames
  attr(final_imputed_list, "ncol") <- ncol(obj)
  attr(final_imputed_list, "subset") <- subset
  class(final_imputed_list) <- c("SlideKnnList", class(final_imputed_list))
  return(final_imputed_list)
}

#' Restore Dimnames of [SlideKnn()]/[knn_imp()] Output
#'
#' [SlideKnn()]/[knn_imp()] can return a list of [bigmemory::big.matrix()] objects
#' that will have dimnames stripped if `options(bigmemory.allow.dimnames = TRUE)`
#' is not set beforehand. This function attempts to restore the dimnames of each
#' object in the returned list using the "rownames" and "colnames" attributes
#' stored on the list.
#'
#' @param obj Output of [SlideKnn()] or [knn_imp()]; a list with class
#'   `SlideKnnList` or `KnnImpList`.
#'
#' @return The invisible input `obj` with dimnames restored if possible.
#'
#' @export
#'
#' @examples
#' to_test <- t(
#'   sim_mat(
#'     n = 280,
#'     m = 100,
#'     perc_NA = 0.5,
#'     perc_col_NA = 1
#'   )$input
#' )
#' # manually strip dimnames from results
#' results <- knn_imp(to_test, k = 10, output = withr::local_tempfile())
#' on.exit(options(bigmemory.allow.dimnames = getOption("bigmemory.allow.dimnames")), add = TRUE)
#' options(bigmemory.allow.dimnames = TRUE)
#' rownames(results[[1]]) <- NULL
#' colnames(results[[1]]) <- NULL
#' # We see that dimnames have been stripped
#' results
#' # `restore_dimnames` stored in the object
#' restore_dimnames(results)
#' results
restore_dimnames <- function(obj) {
  UseMethod("restore_dimnames")
}

#' @rdname restore_dimnames
#' @export
restore_dimnames.SlideKnnList <- function(obj) {
  old_opt <- getOption("bigmemory.allow.dimnames")
  options(bigmemory.allow.dimnames = TRUE)
  can_store_dimnames <- isTRUE(getOption("bigmemory.allow.dimnames")) || is.numeric(obj[[1]])
  on.exit(options(bigmemory.allow.dimnames = old_opt), add = TRUE)

  if (!can_store_dimnames) {
    stop("Failed to set `options(bigmemory.allow.dimnames = TRUE)` to restore dimnames")
  }

  # Get attributes and current dimnames
  attr_rn <- attr(obj, "rownames")
  attr_cn <- attr(obj, "colnames")
  curr_rn <- lapply(obj, rownames)
  curr_cn <- lapply(obj, colnames)

  # Helper to check if all elements are identical
  all_same <- function(x) length(unique(x)) <= 1

  # Early return if nothing to do
  if (is.null(attr_rn) && is.null(attr_cn) && all(sapply(curr_rn, is.null)) && all(sapply(curr_cn, is.null))) {
    return(invisible(obj))
  }

  # Validate consistency
  if (!all_same(curr_rn)) stop("All objects must have the same rownames or all NULL")
  if (!all_same(curr_cn)) stop("All objects must have the same colnames or all NULL")

  # Check for missing attributes when dimnames exist
  if (!is.null(curr_rn[[1]]) && is.null(attr_rn)) stop("Objects have rownames but no rownames attribute found")
  if (!is.null(curr_cn[[1]]) && is.null(attr_cn)) stop("Objects have colnames but no colnames attribute found")

  # Early return if already matching
  if (identical(curr_rn[[1]], attr_rn) && identical(curr_cn[[1]], attr_cn)) {
    return(invisible(obj))
  }

  # Restore dimnames from attributes
  if (is.null(curr_rn[[1]]) && !is.null(attr_rn)) {
    for (i in seq_along(obj)) rownames(obj[[i]]) <- attr_rn
  }
  if (is.null(curr_cn[[1]]) && !is.null(attr_cn)) {
    for (i in seq_along(obj)) colnames(obj[[i]]) <- attr_cn
  }

  return(invisible(obj))
}

#' @rdname restore_dimnames
#' @export
restore_dimnames.KnnImpList <- function(obj) {
  restore_dimnames.SlideKnnList(obj)
}
