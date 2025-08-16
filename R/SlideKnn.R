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
      if (.Platform$OS.type == "windows") {
        Sys.sleep(0.1)
      }
      invisible(gc(verbose = FALSE, full = TRUE))
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
    backingpath = NULL,
    dimnames = NULL) {
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
        descriptorfile = descfiles[i],
        dimnames = dimnames
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

#' K-Nearest Neighbor (k-NN) Imputation
#'
#' @description
#' This function imputes missing values in a numeric matrix using the k-Nearest
#' Neighbors algorithm. It follows a two-stage process: first, it imputes
#' columns with a proportion of missing values below `colmax` using k-NN;
#' second, if requested, any remaining missing values are imputed using the column mean.
#'
#' @details
#' This function calculates the distances for neighbors column-wise. This is an
#' **extremely** important detail. Outside of epi-genomics data, most datasets have
#' people in columns and features (e.g., weight, height, etc.) in rows for imputation.
#' However, in epi-genomics data or intensive longitudinal data, features for the same
#' sample that are spatially closer together may carry mutual information, so we have
#' features in columns and samples in rows; the algorithm will then impute values based
#' on nearby features for the same sample.
#'
#' The distance calculation between columns for identifying nearest neighbors is
#' scaled based on the number of non-missing value pairs. Specifically, the
#' raw distance is penalized by scaling it up for columns that have fewer
#' overlapping observations. This penalizes distances for columns with very few
#' shared observations used for distance calculations. See also [stats::dist()].
#'
#' When `weighted = TRUE`, imputed values are computed as weighted averages where
#' weights are the inverse of distances raised to the power of `dist_pow`. This
#' gives closer neighbors more influence in the imputation.
#'
#' The `tree` parameter allows for faster neighbor search using spatial data structures,
#' but requires pre-filling missing values with column means, which may introduce bias
#' in datasets with high missingness. The overhead of building the tree may make the function
#' slower in lower dimensions and smaller k values. Use "kd" trees for low-dimensional
#' data and "ball" trees for high-dimensional data or when using Manhattan distance.
#'
#' @section Parameter Guidance:
#' For most applications, start with `k = 5`, `colmax = 0.8`, `rowmax = 0.9`.
#' - `colmax = 0.5` means columns with > 50% missing values will use mean imputation
#' - `rowmax = 0.9` means rows with > 90% missing values will cause an error
#' - Set `n_pmm` to at least 10 in real applications for adequate uncertainty quantification
#' - Use `n_pmm > 0` for PMM-based uncertainty (recommended), and `n_pmm = 0` for
#' bootstrap-based uncertainty.
#'
#' @section Performance:
#' - Don't use `tree` unless imputation run-time is getting too slow and there's low missing.
#' - Use `tree = "kd"` for datasets with > 5000 columns and low missingness (< 20%). Benchmark your data.
#' - Use `tree = "ball"` for extremely high-dimensional data or Manhattan distance
#' - Use `output` parameter when `n_imp > 1` to avoid memory issues with large datasets
#' - Consider `subset` to impute only specific columns of interest for efficiency
#'
#' @inheritParams SlideKnn
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param ... Currently not implemented.
#'
#' @return A list of length `n_imp` containing numeric matrices (or [bigmemory::big.matrix()]
#'   objects if `output` is specified) with the same dimensions as `obj`. Missing values are
#'   imputed using k-NN for columns with missingness below `colmax`, and mean imputation
#'   for remaining missing values if `post_imp = TRUE`. Each list element represents an
#'   independent imputation for uncertainty quantification.
#'
#' @note
#' For file-backed storage using `output`, ensure `options(bigmemory.allow.dimnames = TRUE)`
#' is set before calling this function. Else rownames and colnames will have to be
#' restored (same as `obj`)
#'
#' Compared to `impute::impute.knn`, for columns with very high missingness, the
#' mean imputation uses any imputed values with original values for the mean calculation
#' instead of just the original values.
#'
#' When `n_imp` > 1 for large data, it is recommended to use `subset` and provide
#' `output` to use [bigmemory::big.matrix()] to save memory.
#'
#' @references
#' Troyanskaya, Olga, Michael Cantor, Gavin Sherlock, et al.
#' "Missing Value Estimation Methods for DNA Microarrays." Bioinformatics 17,
#' no. 6 (2001): 520–25. https://doi.org/10.1093/bioinformatics/17.6.520.
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
#' )[[1]]
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
#' )[[1]]
#'
#' # PMM imputation with 5 imputations and 5 donors. Increase `n_pmm` in real
#' # data to ensure realistic uncertainty quantification.
#' imputed_pmm <- knn_imp(
#'   obj = t_khanmiss1,
#'   k = 5,
#'   n_imp = 5,
#'   n_pmm = 5,
#'   output = withr::local_tempfile(), # Enables results as list of bigmatrix
#'   seed = 123
#' )
#' # The default of the bigmatrix package is to not allow dimnames so the output
#' # lost the dimnames
#' imputed_pmm[[1]][, ][1:5, 1:5]
#'
#' if (interactive()) {
#'   # But this can be reassigned after setting the options to be TRUE like so
#'   options(bigmemory.allow.dimnames = TRUE)
#'   for (i in length(imputed_pmm)) {
#'     rownames(imputed_pmm[[i]]) <- rownames(t_khanmiss1)
#'     colnames(imputed_pmm[[i]]) <- colnames(t_khanmiss1)
#'   }
#'   imputed_pmm[[1]][, ][1:5, 1:5]
#' }
#'
#' length(imputed_pmm) # Returns 5 imputed datasets
#'
#' # Bootstrap imputation for uncertainty injection
#' imputed_boot <- knn_imp(
#'   obj = t_khanmiss1,
#'   k = 3,
#'   n_imp = 5,
#'   n_pmm = 0, # n_pmm = 0 enables bootstrapping nearest neighbors
#'   seed = 123
#' )
#'
#' length(imputed_boot) # Returns 5 imputed datasets
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
  checkmate::assert_integerish(k, lower = 1, upper = ncol(obj) - 1, len = 1, .var.name = "k")
  checkmate::assert_integerish(cores, lower = 1, len = 1, .var.name = "cores")
  checkmate::assert_numeric(colmax, lower = 0, upper = 1, len = 1, .var.name = "colmax")
  checkmate::assert_numeric(rowmax, lower = 0, upper = 1, len = 1, .var.name = "rowmax")
  checkmate::assert_logical(post_imp, len = 1, null.ok = FALSE, any.missing = FALSE, .var.name = "post_imp")
  checkmate::assert_integerish(n_imp, lower = 1, len = 1, .var.name = "n_imp")
  checkmate::assert_integerish(n_pmm, lower = 0, upper = nrow(obj), len = 1, .var.name = "n_pmm")
  checkmate::assert_integerish(seed, lower = 0, len = 1, .var.name = "seed", null.ok = FALSE)
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
  checkmate::assert_logical(overwrite, len = 1, null.ok = FALSE, any.missing = FALSE, .var.name = "overwrite")

  set.seed(seed)

  # Store original dimnames
  rn <- rownames(obj)
  cn <- colnames(obj)
  dn <- list(rn, cn)
  # Determine if output should be file-backed
  file_backed <- !is.null(output)
  # Handle file-backed output preparation
  paths <- check_result_list(output = output, n_imp = n_imp, overwrite = overwrite)
  if (file_backed) {
    if (!isTRUE(getOption("bigmemory.allow.dimnames"))) {
      dn <- NULL
    }
  }
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
  # Early return for empty subset or no missing data
  if (length(subset) == 0 || !anyNA(obj[, subset, drop = FALSE])) {
    if (length(subset) > 0) {
      message("No missing data in subset columns")
    }
    return(create_result_list(
      data_to_copy = obj,
      file_backed = file_backed,
      n_imp = n_imp,
      backfiles = backfiles,
      descfiles = descfiles,
      backingpath = backingpath,
      dimnames = dn
    ))
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
    backingpath = backingpath,
    dimnames = dn
  )
  # Fill in imputed values
  for (i in seq_len(n_imp)) {
    result_list[[i]][imp_indices] <- imputed_values[, i + 2]
  }

  # Flush file-backed matrices to disk
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
  # Extract regions where counts > 1
  overlaps <- counts_vec > 1
  # Run Length Encoding
  rle_over <- rle(overlaps)
  # Compute start and end positions of the TRUE runs
  ends <- cumsum(rle_over$lengths)
  starts <- ends - rle_over$lengths + 1
  # Return the region and counts_vec
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

#' Sliding k-NN Imputation
#'
#' @description
#' Performs sliding window k-NN imputation on a numeric matrix to handle missing values.
#' The matrix is divided into overlapping windows, and imputation is applied to each window.
#' Overlapping regions are averaged to produce the final imputed matrix.
#'
#' @details
#' The sliding window approach is particularly useful for large datasets where applying
#' k-NN imputation to the entire matrix would be computationally prohibitive. By breaking
#' the data into smaller, overlapping windows, the algorithm maintains local structures
#' while keeping memory usage manageable.
#'
#' See [knn_imp()] for details about the underlying k-NN implementation.
#'
#' @note
#' Setting `n_imp` > 1 (i.e., multiple imputation) is intended to be used with `subset`.
#' and big.matrix.
#'
#' If your `obj` is a [bigmemory::big.matrix()] or description file, you must set
#' `strip_dimnames` to `TRUE` for the output big.matrix to have the same dimnames as
#' the `obj`. Alternatively, you can re-add the dimnames to the output using the
#' dimnames of the original object after setting `options(bigmemory.allow.dimnames = TRUE)`.
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#'   Ensure that the features in the columns are sorted (e.g., by genomic position).
#'   Can also be a path to the description file of, or a [bigmemory::big.matrix()].
#' @param n_feat Integer specifying the number of features (columns) in each window.
#'   Must be between 2 and the number of columns in `obj`.
#' @param subset Character vector of column names or integer vector of column indices
#'   specifying the subset of columns on which to perform imputation. If `NULL` (default),
#'   all columns are included.
#' @param n_overlap Integer specifying the number of features to overlap between
#'   consecutive windows. Default is 10. Must be between 0 and `n_feat - 1`.
#' @param k Integer specifying the number of nearest neighbors to use for imputation.
#'   Must be between 1 and (`n_feat` - 1).
#' @param rowmax Numeric between 0 and 1 specifying the maximum allowable proportion
#'   of missing values in any row. If exceeded, the function stops with an error.
#' @param colmax Numeric between 0 and 1 specifying the threshold proportion of
#'   missing values in a column above which the column is imputed using the mean
#'   instead of k-NN if `post_imp` is `TRUE`.
#' @param cores Integer specifying the number of cores to use for parallel computation
#'   of distances. Default is 1.
#' @param method Character string specifying the distance metric for k-NN. One of
#'   `"euclidean"` or `"manhattan"`. Defaults to `"euclidean"`.
#' @param tree Character string specifying the k-NN method. `NULL` (default) uses
#'   brute-force search. `"kd"` uses KDTree and `"ball"` uses BallTree as implemented by
#'   the mlpack package, where missing values are first filled with column means
#'   (biased at high percentage of missing values and less efficient at lower dimensions).
#' @param post_imp Logical; if `TRUE` (default), any missing values remaining after
#'   k-NN imputation will be imputed with [mean_impute_col()].
#' @param weighted Logical; controls whether the imputed value should be a simple
#'   mean or weighted mean by inverse distance. Default is `FALSE`.
#' @param dist_pow A positive double that controls the penalty for larger distances
#'   in the weighted mean imputation. Must be greater than zero: values between 0
#'   and 1 apply a softer penalty, 1 is linear (default), and values greater than
#'   1 apply a harsher penalty.
#' @param n_imp Integer specifying the number of imputations to perform. Forced
#'   to 1 if `n_pmm` = -1.
#' @param n_pmm Integer controlling the multiple imputation method.
#'   If `n_pmm` = -1 (default): Deterministic single imputation.
#'   If `n_pmm` > 0: (recommended for MI) PMM multiple imputation using `n_pmm` closest
#'   donors (will not exceed available non-missing values row-wise).
#'   If `n_pmm` = 0: Bootstrap multiple imputation via resampling from k-Nearest Neighbors.
#' @param seed Integer; random seed for reproducible bootstrap sampling. Default is 42.
#' @param .progress Logical; if `TRUE`, show a progress bar. Default is `FALSE`.
#' @param output Character; path to save the output big.matrix. Format should be
#'   `path/stem` (e.g. "./my_results"). If provided, the back end will switch
#'   to [bigmemory::filebacked.big.matrix()] instead of in-memory matrix.
#' @param overwrite Logical; if `TRUE`, overwrite existing files at `output`. Default to `FALSE`.
#' @param strip_dimnames Logical; if `FALSE` (default), dimnames will not be removed,
#'   which will increase memory usage. Should be set to `TRUE` to save memory from
#'   overhead, especially when `cores` > 1. See details.
#'
#' @return A list of numeric matrices/big.matrices of the same dimensions as `obj`
#'   with missing values imputed. Length of list equals `n_imp`. If `obj` is file-backed,
#'   returns a list of big.matrix objects; otherwise returns regular matrices.
#'
#' @seealso [knn_imp()], [mean_impute_col()], [bigmemory::big.matrix()]
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
#' # Access the result (returns a list, even for single imputation)
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
#'   # strip_dimnames = TRUE  # Recommended to set to TRUE for efficiency
#' )
#'
#' # Because strip_dimnames is `FALSE` and the default of the bigmatrix package
#' # is to not allow dimnames, the output lost the dimnames
#' imputed_big[[1]][, ][1:5, 1:5]
#'
#' # But this can be reassigned after setting the options to be TRUE like so
#' if (interactive()) {
#'   options(bigmemory.allow.dimnames = TRUE)
#'   rownames(imputed_big[[1]]) <- rownames(big_data)
#'   colnames(imputed_big[[1]]) <- colnames(big_data)
#'   imputed_big[[1]][, ][1:5, 1:5]
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
#'   output = withr::local_tempfile(),
#'   .progress = TRUE
#' )
#'
#' # ========================================
#' # Parallel Processing
#' # ========================================
#'
#' @examplesIf interactive()
#'
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
#'   cores = 4,
#'   strip_dimnames = TRUE # Important for parallel efficiency
#' )
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
    overwrite = FALSE,
    strip_dimnames = FALSE) {
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
  if (n_imp > 1 && n_pmm == 0 && weighted) {
    warning("If bootstrapping nearest neighbors, weighted will be forced to FALSE")
    weighted <- FALSE
  }
  # Default is options(bigmemory.allow.dimnames = NULL or FALSE) which prevents
  # strip names and assigning names to outputs so we have to guard
  if (strip_dimnames) {
    on.exit(options(bigmemory.allow.dimnames = getOption("bigmemory.allow.dimnames")))
    options(bigmemory.allow.dimnames = TRUE)
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
  checkmate::assert_integerish(n_feat, lower = 2, upper = ncol(obj), len = 1, null.ok = FALSE, .var.name = "n_feat")
  checkmate::assert_int(n_overlap, lower = 0, upper = n_feat - 1, null.ok = FALSE, .var.name = "n_overlap")
  checkmate::assert_int(k, lower = 1, upper = n_feat - 1, null.ok = FALSE, .var.name = "k")

  # Determine if output should be file-backed
  file_backed <- !is.null(output)
  # Warning for multiple imputations without file backing
  if (n_imp > 1 && !file_backed) {
    warning("n_imp > 1 should be used with output to file-backed big.matrix. Provide `output` parameter to enable file backing.")
  }
  if (n_imp > 1 && n_pmm == -1) {
    n_imp <- 1
  }
  # Remove obj names temporarily to reduce size of pointers being passed to workers
  rn <- rownames(obj)
  cn <- colnames(obj)
  dn <- list(rn, cn)
  if (strip_dimnames) {
    rownames(obj) <- NULL
    colnames(obj) <- NULL
    dn <- NULL
  }
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
    # If subset is NULL, set it to all columns
    subset <- seq_len(ncol(obj))
  }
  # Windowing Logic ----
  # Calculate the total number of steps/windows needed.
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
  # Offset the subset_list to make indices relative to each window's start
  subset_list <- lapply(seq_along(start), function(i) {
    # Find range of subset that falls in window [start[i], end[i]]
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
  # Get final imputed info (multiple if n_imp > 1)
  final_imputed_info <- check_result_list(output, n_imp, overwrite)
  # Create all final imputed big.matrices upfront
  final_imputed_list <- lapply(seq_len(n_imp), function(imp_idx) {
    bigmemory::big.matrix(
      nrow = nr,
      ncol = nc,
      type = "double",
      init = 0.0,
      backingpath = final_imputed_info$backingpath,
      descriptorfile = final_imputed_info$descfiles[[imp_idx]],
      backingfile = final_imputed_info$backfiles[[imp_idx]]
    )
  })
  # Overlap regions (computed once, shared)
  overlap <- find_overlap_regions(start, end)
  # Process each imputation sequentially for each n_imp
  for (imp_idx in seq_len(n_imp)) {
    if (.progress) {
      message(sprintf("Processing imputation %d/%d", imp_idx, n_imp))
    }
    # Seed for this imputation
    seed_imp <- seed + (imp_idx - 1)
    # Create single intermediate for this imp (overwritten if file-backed)
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
            # Realize in memory
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
          # Fill single intermediate matrix
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
    # Get the current final imputed
    final_imputed <- final_imputed_list[[imp_idx]]
    ## Averaging ----
    if (.progress) {
      message("Step 2/3: Overlapping")
    }
    # Add the windows from intermediate
    purrr::walk(
      seq_along(start),
      function(i) {
        window_cols <- start[i]:end[i]
        final_imputed[, window_cols] <- {
          final_imputed[, window_cols] + intermediate[, offset_start[i]:offset_end[i]]
        }
      },
      .progress = FALSE
    )
    if (file_backed) {
      bigmemory::flush(final_imputed)
    }
    if (.progress) {
      message("Step 3/3: Averaging")
    }
    # Process averaging
    bigmem_avg(
      final_imputed@address,
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
        final_imputed@address,
        col_indices = subset,
        cores = cores
      )
    }
    # Restore names
    if (strip_dimnames) {
      rownames(final_imputed) <- rn
      colnames(final_imputed) <- cn
    }
    # If not file-backed, convert back to normal matrix
    if (!file_backed) {
      final_imputed_list[[imp_idx]] <- bigmemory::as.matrix(final_imputed)
      gc()
    }
  }
  # Return results
  return(final_imputed_list)
}
