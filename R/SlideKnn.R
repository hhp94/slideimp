#' Sliding KNN Imputation
#'
#' @description
#' Performs sliding window KNN imputation on a numeric matrix to handle missing values.
#' The matrix is divided into overlapping windows, and imputation is applied to each window.
#' Overlapping regions are averaged to produce the final imputed matrix.
#'
#' @details
#' The sliding window approach is particularly useful for large datasets where applying
#' KNN imputation to the entire matrix would be computationally prohibitive. By breaking
#' the data into smaller, overlapping windows, the algorithm maintains local correlation
#' structures while keeping memory usage manageable.
#'
#' See [knn_imp()] for details about the underlying KNN implementation.
#'
#' @note
#' Setting `nboot` > 1 is meant to be used only with `subset`.
#'
#' If your `obj` is a [bigmemory::big.matrix()] or description file, you have to set
#' `strip_dimnames` to `TRUE` for the output big.matrix to have the same dimnames as
#' the `obj`. Otherwise, you can also re-add the dimnames to the output using the
#' dimnames of the original object after setting `options(bigmemory.allow.dimnames = TRUE)`.
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#'   Make sure the features in the columns are sorted (e.g., by genomic position).
#'   Can also be a path to the description file of, or a [bigmemory::big.matrix()].
#' @param n_feat Integer specifying the number of features (columns) in each window.
#'   Must be between 2 and the number of columns in `obj`.
#' @param subset Character vector of column names or integer vector of column indices
#'   specifying the subset of columns to perform imputation. If `NULL` (default),
#'   all columns are included.
#' @param n_overlap Integer specifying the number of features to overlap between
#'   consecutive windows. Default is 10. Must be between 0 and `n_feat - 1`.
#' @param k Integer specifying the number of nearest neighbors to use for imputation.
#'   Must be between 1 and (`n_feat` - 1).
#' @param rowmax Numeric between 0 and 1 specifying the maximum allowable proportion
#'   of missing values in any row. If exceeded, the function stops with an error.
#' @param colmax Numeric between 0 and 1 specifying the threshold proportion of
#'   missing values in a column above which the column is imputed using the mean
#'   instead of k-NN if `post_imp` is true.
#' @param cores Integer specifying the number of cores to use for parallel computation
#'   of distances. Default is 1.
#' @param method Character string specifying the distance metric for k-NN. One of
#'   `"euclidean"`, `"manhattan"`, or `"impute.knn"`. Defaults to `"euclidean"`.
#' @param tree Character string specifying the k-NN method. `NULL` (default) uses
#'   brute-force search. `"kd"` uses KDTree and `"ball"` uses BallTree as implemented by
#'   the mlpack package where missing values are first filled with column means
#'   (biased at high percentage missing and less efficient at lower dimensions).
#' @param post_imp Logical; if `TRUE` (default), any missing values remaining after
#'   k-NN imputations will be imputed with [mean_impute_col()].
#' @param weighted Logical; controls whether the imputed value should be a simple
#'   mean or weighted mean by inverse distance. Default is `FALSE`.
#' @param dist_pow A positive double that controls the penalty for larger distances
#'   in the weighted mean imputation. Must be greater than zero: values between 0
#'   and 1 apply a softer penalty, 1 is linear (default), and values greater than
#'   1 apply a harsher penalty.
#' @param nboot Integer specifying the number of bootstrap imputations to perform.
#'   Default is 1.
#' @param seed Integer; random seed for reproducible bootstrap sampling. Default is 42.
#' @param .progress Logical; if `TRUE`, show a progress bar. Default is `FALSE`.
#' @param output Character; path to save the output big.matrix if `obj` is file-backed.
#'   Required when `obj` is a [bigmemory::big.matrix()] or path to big.matrix descriptor.
#' @param overwrite Logical; if `TRUE` (default), overwrite existing files at `output`.
#' @param block Integer; block size for processing large matrices. If `NULL` (default),
#'   calculated automatically based on the matrix size and number of cores.
#' @param strip_dimnames Logical; if `FALSE` (default), dimnames will not be removed
#'   which will increase memory usage. Should set to `TRUE` to save memory from
#'   overhead especially when `cores` > 1. See details.
#'
#' @return A list of numeric matrices/big.matrices of the same dimensions as `obj`
#'   with missing values imputed. Length of list equals `nboot`. If `obj` is file-backed,
#'   returns a list of big.matrix objects; otherwise returns regular matrices.
#'
#' @seealso [knn_imp()], [mean_impute_col()], [bigmemory::big.matrix()]
#'
#' @examples
#' data(khanmiss1)
#'
#' # Basic sliding KNN imputation with bootstrap
#' # Set `strip_dimnames` to `TRUE` with `options(bigmemory.allow.dimnames = TRUE)`
#' # when running multiple cores to minimize overhead.
#'
#' # Using weighted imputation
#' imputed_weighted <- SlideKnn(
#'   t(khanmiss1),
#'   k = 5,
#'   n_feat = 50,
#'   weighted = TRUE,
#'   dist_pow = 2
#' )
#' # Bootstrap version
#' # imputed <- SlideKnn(t(khanmiss1), k = 10, n_feat = 100, n_overlap = 10, nboot = 3)
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
    method = c("euclidean", "manhattan", "impute.knn"),
    tree = NULL,
    post_imp = TRUE,
    weighted = FALSE,
    dist_pow = 1,
    nboot = 1,
    seed = 42,
    .progress = FALSE,
    output = NULL,
    overwrite = TRUE,
    block = NULL,
    strip_dimnames = FALSE) {
  # Pre-conditioning ----
  method <- match.arg(method)
  checkmate::assert_integerish(n_feat, lower = 2, upper = ncol(obj), len = 1, null.ok = FALSE, .var.name = "n_feat")
  checkmate::assert_integerish(n_overlap, lower = 0, upper = n_feat - 1, len = 1, null.ok = FALSE, .var.name = "n_overlap")
  checkmate::assert_integerish(k, lower = 1, upper = n_feat - 1, len = 1, null.ok = FALSE, .var.name = "k")
  checkmate::assert_numeric(rowmax, lower = 0, upper = 1, len = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_numeric(colmax, lower = 0, upper = 1, len = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")
  checkmate::assert_integerish(nboot, lower = 1, len = 1, null.ok = FALSE, .var.name = "nboot")
  checkmate::assert_integerish(seed, lower = 0, len = 1, null.ok = FALSE, .var.name = "seed")
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
  checkmate::assert_character(output, len = 1, null.ok = TRUE, .var.name = "output")
  checkmate::assert(
    checkmate::check_character(subset, min.len = 1, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
    checkmate::check_integerish(
      subset,
      lower = 1, upper = ncol(obj), min.len = 1, any.missing = FALSE, null.ok = TRUE, unique = TRUE
    ),
    .var.name = "subset"
  )
  if (nboot > 1 && weighted) {
    warning("If nboot > 1, weighted will be forced to FALSE")
  }

  # set.seed() ----
  set.seed(seed)

  # Check big.matrix ----
  if (strip_dimnames) {
    on.exit(options(bigmemory.allow.dimnames = getOption("bigmemory.allow.dimnames")))
    options(bigmemory.allow.dimnames = TRUE)
  }

  ## obj ----
  file_backed <- if (is.character(obj)) {
    if (!fs::file_exists(obj)) {
      stop("The provided path to the big.matrix descriptor and/or filebacking can not be found.")
    }
    # Handle if path is to .bin instead of .desc
    if (fs::path_ext(obj) == "bin") {
      desc_path <- fs::path_ext_set(obj, "desc")
      if (fs::file_exists(desc_path)) {
        obj <- desc_path
      } else {
        stop("Descriptor file not found for the provided backing file.")
      }
    }
    obj <- bigmemory::attach.big.matrix(obj)
    TRUE
  } else if (bigmemory::is.big.matrix(obj)) {
    TRUE
  } else if (is.matrix(obj) && is.numeric(obj)) {
    obj <- bigmemory::as.big.matrix(obj, type = "double")
    FALSE
  } else {
    stop("`obj` has to be a numeric matrix or a (descriptor of a) double bigmemory::big.matrix.")
  }
  if (nboot > 1 && !file_backed) {
    warning("If nboot > 1 should be used with output to bigmemory::big.matrix. Change `output` to a path to enable this.")
  }
  # Remove obj names temporary to reduce size of pointers being passed to workers
  rn <- rownames(obj)
  cn <- colnames(obj)
  if (strip_dimnames) {
    rownames(obj) <- NULL
    colnames(obj) <- NULL
  }
  if (file_backed && is.null(output)) {
    stop("`output` must be provided if obj is a big.matrix or path to big.matrix descriptor (i.e., './output.bin').")
  }
  if (!is.null(output)) {
    output <- fs::as_fs_path(output)
    backingpath <- fs::path_dir(output)
    if (!all(fs::file_access(backingpath, mode = c("exists", "read", "write")))) {
      stop("Provided path to `output` doesn't exist or is not readable/writable")
    }
    file_name <- fs::path_file(output)
    ext <- fs::path_ext(file_name)
    base_no_ext <- fs::path_ext_remove(file_name)
    if (ext == "" || ext == "desc") {
      backfile <- fs::path_ext_set(base_no_ext, "bin")
      descfile <- fs::path_ext_set(base_no_ext, "desc")
    } else {
      backfile <- fs::path_ext_set(base_no_ext, paste0(ext, ".bin"))
      descfile <- fs::path_ext_set(base_no_ext, paste0(ext, ".desc"))
    }
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
  subset_list <- purrr::map2(
    start,
    end,
    function(x, y) {
      r <- intersect(subset, x:y)
      r <- r[!is.na(r)]
      if (length(r) == length(x:y)) {
        # Don't need this if line, But its fine. Makes the code more robust
        return(NULL)
      }
      if (length(r) > 0) {
        return(sort(r - x + 1))
      } else {
        return(r)
      }
    }
  )
  # Then, we calculate the offsets needed to subset the intermediate matrix
  width <- end - start + 1
  offset_start <- c(1, cumsum(width)[-length(width)] + 1)
  offset_end <- cumsum(width)

  # Sliding Imputation ----
  ## Init ----
  nr <- nrow(obj)
  nc <- ncol(obj)
  if (is.null(block)) {
    # whichever is smaller, ncol(obj) or
    # whichever is larger, 1, number of cores * 100 iterations each core, or 0.005 * ncol()
    block <- min(nc, max(1, nc %/% (cores * 100), floor(0.005 * nc)))
  }
  checkmate::assert_integerish(
    block,
    lower = 1,
    upper = ncol(obj),
    len = 1,
    null.ok = FALSE,
    .var.name = "block"
  )

  # Create lists for multiple bootstrap iterations
  intermediate_list <- vector("list", nboot)
  final_imputed_list <- vector("list", nboot)
  counts_list <- vector("list", nboot)

  ### file_backed ----
  if (file_backed) {
    checkmate::assert_flag(overwrite, .var.name = "overwrite", null.ok = FALSE)

    # Create temporary directory for intermediate files
    temp_dir <- withr::local_tempdir(pattern = paste0("SlideKnn_", Sys.getpid()))

    # For file-backed matrices, create separate files for each bootstrap
    for (boot_idx in seq_len(nboot)) {
      # Check and handle output files
      if (nboot > 1) {
        boot_backfile <- fs::path_ext_set(fs::path_ext_remove(backfile), paste0("_boot", boot_idx, ".", fs::path_ext(backfile)))
        boot_descfile <- fs::path_ext_set(fs::path_ext_remove(descfile), paste0("_boot", boot_idx, ".", fs::path_ext(descfile)))
      } else {
        boot_backfile <- backfile
        boot_descfile <- descfile
      }

      files_to_check <- file.path(backingpath, c(boot_backfile, boot_descfile))
      if (any(fs::file_exists(files_to_check))) {
        if (!overwrite) {
          stop("Provided `output` already exists but overwrite is `FALSE`")
        }
        # Proceed to delete if overwrite is TRUE
        unlink(files_to_check, force = TRUE)
        if (any(file.exists(files_to_check))) {
          stop("Failed to delete existing files for `output`")
        }
      }

      # Create intermediate matrix for this bootstrap
      intermediate_list[[boot_idx]] <- bigmemory::filebacked.big.matrix(
        nrow = nr,
        ncol = sum(width),
        type = "double",
        init = 0.0,
        backingfile = paste0("intermediate_boot", boot_idx, ".bin"),
        descriptorfile = paste0("intermediate_boot", boot_idx, ".desc"),
        backingpath = temp_dir
      )

      # Create counts matrix for this bootstrap
      counts_list[[boot_idx]] <- bigmemory::filebacked.big.matrix(
        nrow = nc,
        ncol = 1,
        type = "double",
        init = 0.0,
        backingfile = paste0("counts_boot", boot_idx, ".bin"),
        descriptorfile = paste0("counts_boot", boot_idx, ".desc"),
        backingpath = temp_dir
      )

      # Create final imputed matrix for this bootstrap
      final_imputed_list[[boot_idx]] <- bigmemory::filebacked.big.matrix(
        nrow = nr,
        ncol = nc,
        type = "double",
        init = 0.0,
        backingpath = backingpath,
        backingfile = boot_backfile,
        descriptorfile = boot_descfile
      )
    }
  } else {
    ### in-memory ----
    for (boot_idx in seq_len(nboot)) {
      intermediate_list[[boot_idx]] <- bigmemory::big.matrix(
        nrow = nr,
        ncol = sum(width),
        type = "double",
        init = 0.0
      )

      counts_list[[boot_idx]] <- bigmemory::big.matrix(
        nrow = nc,
        ncol = 1,
        type = "double",
        init = 0.0
      )

      final_imputed_list[[boot_idx]] <- bigmemory::big.matrix(
        nrow = nr,
        ncol = nc,
        type = "double",
        init = 0.0
      )
    }
  }

  # Create descriptors for all bootstrap matrices
  obj_desc <- bigmemory::describe(obj)
  intermediate_desc_list <- lapply(intermediate_list, bigmemory::describe)
  counts_desc_list <- lapply(counts_list, bigmemory::describe)
  final_imputed_desc_list <- lapply(final_imputed_list, bigmemory::describe)

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

        # Get imputation results for all bootstrap iterations
        imp_list <- impute_knn(
          obj = obj_big[, window_cols, drop = FALSE],
          k = k,
          rowmax = rowmax,
          colmax = colmax,
          cores = 1L, # For Slide KNN, fix cores = 1
          method = method,
          post_imp = post_imp,
          weighted = weighted,
          dist_pow = dist_pow,
          subset = subset_list[[i]],
          nboot = nboot,
          seed = seed,
          knn_imp = knn_imp,
          tree = tree,
          impute_knn_brute = impute_knn_brute,
          impute_knn_mlpack = impute_knn_mlpack,
          mean_impute_col = mean_impute_col
        )

        # Fill intermediate matrices for each bootstrap iteration
        for (boot_idx in seq_len(nboot)) {
          intermediate_big <- bigmemory::attach.big.matrix(intermediate_desc_list[[boot_idx]])
          intermediate_big[, offset_start[i]:offset_end[i]] <- imp_list[[boot_idx]]
        }
      },
      impute_knn = impute_knn,
      subset_list = subset_list,
      knn_imp = knn_imp,
      weighted = weighted,
      dist_pow = dist_pow,
      impute_knn_brute = impute_knn_brute,
      mean_impute_col = mean_impute_col,
      impute_knn_mlpack = impute_knn_mlpack,
      start = start,
      end = end,
      obj_desc = obj_desc,
      tree = tree,
      intermediate_desc_list = intermediate_desc_list,
      k = k,
      rowmax = rowmax,
      colmax = colmax,
      method = method,
      post_imp = post_imp,
      offset_start = offset_start,
      offset_end = offset_end,
      nboot = nboot,
      seed = seed
    ),
    .progress = .progress
  )

  ## Averaging ----
  if (.progress) {
    message("Step 2/3: Overlapping")
  }

  # Process each bootstrap iteration
  for (boot_idx in seq_len(nboot)) {
    purrr::walk(
      seq_along(start),
      function(i) {
        window_cols <- start[i]:end[i]
        final_imputed_list[[boot_idx]][, window_cols] <- {
          final_imputed_list[[boot_idx]][, window_cols] + intermediate_list[[boot_idx]][, offset_start[i]:offset_end[i]]
        }
        counts_list[[boot_idx]][window_cols] <- counts_list[[boot_idx]][window_cols] + 1
      },
      .progress = FALSE
    )
  }

  # post-processing, average out values.
  # Block is the size of the block of columns to process the data by. Have to
  # recalculate the subset relative to w_start
  w_start <- seq(1, nc, by = block)
  w_end <- c(w_start[-1] - 1, nc)
  w_subset_list <- purrr::map2(
    w_start,
    w_end,
    function(x, y) {
      r <- intersect(subset, x:y)
      r <- r[!is.na(r)]
      if (length(r) == length(x:y)) {
        # Don't need this if line, But its fine. Makes the code more robust
        return(NULL)
      }
      if (length(r) > 0) {
        return(sort(r - x + 1))
      } else {
        return(r)
      }
    }
  )

  if (.progress) {
    message("Step 3/3: Averaging")
  }

  # Process averaging for each bootstrap iteration
  for (boot_idx in seq_len(nboot)) {
    purrr::walk(
      seq_along(w_start),
      fn(
        function(i) {
          window_cols <- w_start[i]:w_end[i]
          counts_big <- bigmemory::attach.big.matrix(counts_desc_list[[boot_idx]])
          final_imputed_big <- bigmemory::attach.big.matrix(final_imputed_desc_list[[boot_idx]])
          average <- counts_big[window_cols] > 1

          if (length(window_cols[average]) > 0) {
            final_imputed_big[, window_cols[average]] <- sweep(
              final_imputed_big[, window_cols[average], drop = F],
              MARGIN = 2,
              STATS = counts_big[window_cols[average]],
              FUN = "/"
            )
          }
        },
        counts_desc_list = counts_desc_list,
        final_imputed_desc_list = final_imputed_desc_list,
        w_start = w_start,
        w_end = w_end,
        boot_idx = boot_idx
      ),
      .progress = FALSE
    )
  }

  # Remove all the intermediate files to save memory
  if (file_backed) {
    rm(intermediate_list, counts_list)
    gc()
  }

  ## post_imp ----
  if (post_imp) {
    if (.progress) {
      message("Post-imputation")
    }

    for (boot_idx in seq_len(nboot)) {
      purrr::walk(
        seq_along(w_start),
        fn(
          function(i) {
            window_cols <- w_start[i]:w_end[i]
            final_imputed_big <- bigmemory::attach.big.matrix(final_imputed_desc_list[[boot_idx]])
            if (anyNA(final_imputed_big[, window_cols])) {
              final_imputed_big[, window_cols] <- mean_impute_col(
                final_imputed_big[, window_cols, drop = FALSE],
                subset = w_subset_list[[i]]
              )
            }
          },
          mean_impute_col = mean_impute_col,
          final_imputed_desc_list = final_imputed_desc_list,
          w_subset_list = w_subset_list,
          w_start = w_start,
          w_end = w_end,
          boot_idx = boot_idx
        ),
        .progress = FALSE
      )
    }
  }

  # Restore names for all bootstrap iterations
  if (strip_dimnames) {
    rownames(obj) <- rn
    colnames(obj) <- cn
    for (boot_idx in 1:nboot) {
      rownames(final_imputed_list[[boot_idx]]) <- rn
      colnames(final_imputed_list[[boot_idx]]) <- cn
    }
  }

  # Return results
  if (file_backed) {
    return(final_imputed_list)
  } else {
    # Convert to regular matrices and restore names
    result_list <- vector("list", nboot)
    for (boot_idx in seq_len(nboot)) {
      out <- bigmemory::as.matrix(final_imputed_list[[boot_idx]])
      rownames(out) <- rn
      colnames(out) <- cn
      result_list[[boot_idx]] <- out
    }
    return(result_list)
  }
}

#' KNN Imputation Wrapper
#'
#' @description
#' A wrapper function for KNN imputation that filters rows based on missing value
#' proportions before imputing. Rows with missing proportions exceeding `rowmax`
#' are not imputed in this step to avoid throwing exceptions in [knn_imp()].
#'
#' @details
#' This function serves as an intermediate layer between the sliding window algorithm
#' and the core KNN imputation function. It ensures that only rows with acceptable
#' missing value proportions are processed, preventing errors and improving
#' computational efficiency.
#'
#' @inheritParams SlideKnn
#' @param knn_imp Function object for KNN imputation. Required for function crating
#'   in parallel environments.
#' @param ... Additional arguments passed to `knn_imp`.
#'
#' @return A list of imputed matrices, with only qualifying rows imputed via KNN;
#'   others remain unchanged. Length of list equals `nboot`.
#'
#' @seealso [knn_imp()], [SlideKnn()]
#'
#' @keywords internal
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
    nboot,
    seed,
    ...) {
  na_mat <- is.na(obj)
  # Determine 'good_rows': rows where the proportion of NAs is less than 'rowmax'.
  # Rows with too many NAs (as defined by rowmax) are excluded from imputation.
  good_rows <- rowSums(na_mat) / ncol(na_mat) < rowmax
  # If no rows meet the criteria for imputation, return list with original object.
  if (sum(good_rows) == 0) {
    return(replicate(nboot, obj, simplify = FALSE))
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
    nboot = nboot,
    seed = seed,
    ...
  )

  # Process each imputed object in the list
  result_list <- vector("list", length(imputed_good_list))

  for (i in seq_along(imputed_good_list)) {
    # Initialize the result matrix with the original object's values.
    result <- obj
    # Place the imputed values from this bootstrap iteration back into the
    # corresponding 'good_rows' in the result matrix.
    result[good_rows, ] <- imputed_good_list[[i]]
    result_list[[i]] <- result
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
#' This implementation calculates the distances for neighbors column-wise. This is an
#' **extremely** important detail. Outside of microarray data, most datasets have
#' people in columns and features (e.g., weight, height, etc.) in rows for imputation.
#' However, in microarray data, genes or CpG sites for the same sample that are
#' spatially closer together may carry mutual information, so you can place genes/CpGs in columns
#' and samples in rows; the algorithm will then impute values based on nearby genes/CpGs
#' for the same sample.
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
#' in datasets with high missingness. The overhead of building tree may makes the function
#' slower in lower dimensions and k.
#'
#' @note Compared to `impute::impute.knn`, for columns with very high missingness, the
#' mean imputation uses the imputed values and original values for the mean calculation
#' instead of just the original values.
#'
#' When `nboot` > 1, output should be specified to use [bigmemory::big.matrix()] to
#' save memory.
#'
#' @inheritParams SlideKnn
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param output Character; path to save the output as list of [bigmemory::big.matrix()]
#' to save memory. Highly recommended for `nboot` > 1.
#' @param ... Currently not implemented.
#'
#' @return A list of numeric matrices or big.matrix of the same dimensions as `obj`
#' with missing values imputed. Length of list equals `nboot`.
#'
#' @seealso [SlideKnn()], [mean_impute_col()], [mean_impute_row()]
#'
#' @references
#' Troyanskaya, O., Cantor, M., Sherlock, G., Brown, P., Hastie, T., Tibshirani, R.,
#' Botstein, D. and Altman, R.B. (2001) Missing value estimation methods for DNA
#' microarrays. *Bioinformatics*, **17**(6), 520-525.
#'
#' @examples
#' # Load example data
#' data(khanmiss1)
#' sum(is.na(khanmiss1))
#'
#' # Perform k-NN imputation. `khanmiss1` stores genes in rows so we transpose.
#' # Set method to "impute.knn" to mimic how distance is scaled in impute::impute.knn.
#' imputed <- knn_imp(
#'   obj = t(khanmiss1),
#'   k = 3,
#'   colmax = 0.5,
#'   rowmax = 0.8,
#'   method = "euclidean"
#' )[[1]]
#'
#' # Check results
#' imputed[1:5, 1:20]
#' sum(is.na(imputed))
#'
#' # Using weighted imputation with custom distance power
#' imputed_weighted <- knn_imp(
#'   obj = t(khanmiss1),
#'   k = 5,
#'   weighted = TRUE,
#'   dist_pow = 2,
#'   method = "euclidean"
#' )[[1]]
#'
#' # Bootstrap imputation for uncertainty quantification
#' imputed_boot <- knn_imp(
#'   obj = t(khanmiss1),
#'   k = 3,
#'   nboot = 5,
#'   seed = 123
#' )
#' length(imputed_boot) # Returns 5 imputed datasets
#'
#' @export
knn_imp <- function(
    obj,
    k,
    colmax = 0.9,
    rowmax = 0.9,
    method = c("euclidean", "manhattan", "impute.knn"),
    cores = 1,
    post_imp = TRUE,
    subset = NULL,
    weighted = FALSE,
    dist_pow = 1,
    tree = NULL,
    nboot = 1,
    seed = 42,
    output = NULL,
    overwrite = FALSE,
    ...) {
  # Pre-conditioning
  method <- match.arg(method)
  checkmate::assert_matrix(obj, mode = "numeric", min.rows = 1, min.cols = 2, null.ok = FALSE, .var.name = "obj")
  checkmate::assert_true(sum(is.infinite(obj)) == 0, .var.name = "obj")
  checkmate::assert_integerish(k, lower = 1, upper = ncol(obj) - 1, len = 1, .var.name = "k")
  checkmate::assert_integerish(cores, lower = 1, len = 1, .var.name = "cores")
  checkmate::assert_numeric(colmax, lower = 0, upper = 1, len = 1, .var.name = "colmax")
  checkmate::assert_numeric(rowmax, lower = 0, upper = 1, len = 1, .var.name = "rowmax")
  checkmate::assert_logical(post_imp, len = 1, null.ok = FALSE, any.missing = FALSE, .var.name = "post_imp")
  checkmate::assert_integerish(nboot, lower = 1, len = 1, .var.name = "nboot")
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
  checkmate::assert_string(output, null.ok = TRUE, .var.name = "output")
  checkmate::assert_logical(overwrite, len = 1, null.ok = FALSE, any.missing = FALSE, .var.name = "overwrite")
  if (!is.null(tree)) {
    if (method == "impute.knn") {
      method <- "euclidean"
    }
  }
  set.seed(seed)
  # Store original dimnames
  rn <- rownames(obj)
  cn <- colnames(obj)
  # Check bigmemory options if output is specified
  file_backed <- !is.null(output)
  # Pre-compute file paths for all bootstrap iterations
  if (file_backed) {
    if (!isTRUE(getOption("bigmemory.allow.dimnames"))) {
      stop("`bigmemory.allow.dimnames` must be TRUE. Set it with: options(bigmemory.allow.dimnames = TRUE)")
    }
    # Prepare file paths
    output <- fs::as_fs_path(output)
    backingpath <- fs::path_dir(output)
    if (!all(fs::file_access(backingpath, mode = c("exists", "read", "write")))) {
      stop("Provided path to `output` doesn't exist or is not readable/writable")
    }
    file_name <- fs::path_file(output)
    ext <- fs::path_ext(file_name)
    base_no_ext <- fs::path_ext_remove(file_name)
    # Pre-compute all file paths
    boot_indices <- seq_len(nboot)
    suffix <- if (nboot == 1) "" else paste0("_boot", boot_indices)
    if (ext == "" || ext == "desc") {
      backfiles <- paste0(base_no_ext, suffix, ".bin")
      descfiles <- paste0(base_no_ext, suffix, ".desc")
    } else {
      backfiles <- paste0(base_no_ext, suffix, ".", ext, ".bin")
      descfiles <- paste0(base_no_ext, suffix, ".", ext, ".desc")
    }
    # Check if files exist and overwrite is FALSE
    for (i in boot_indices) {
      backfile_path <- fs::path(backingpath, backfiles[i])
      descfile_path <- fs::path(backingpath, descfiles[i])
      if (overwrite) {
        # Proceed to delete if overwrite is TRUE
        unlink(backfile_path, force = TRUE)
        unlink(descfile_path, force = TRUE)
        if (fs::file_exists(backfile_path) || fs::file_exists(descfile_path)) {
          stop("Failed to delete existing files for `output`")
        }
      } else {
        if (fs::file_exists(backfile_path) || fs::file_exists(descfile_path)) {
          stop(
            paste0(
              "Output files already exist for bootstrap ",
              i,
              ". Set overwrite = TRUE to overwrite them."
            )
          )
        }
      }
    }
  } else {
    # Not file-backed - set to NULL
    backingpath <- NULL
    boot_indices <- seq_len(nboot)
    backfiles <- NULL
    descfiles <- NULL
  }
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
  if (length(subset) == 0 || !sum(is.na(obj[, subset, drop = FALSE])) > 0) {
    if (length(subset) > 0) {
      message("No missing data in subset columns")
    }
    return(create_result_list(
      data_to_copy = obj,
      file_backed = file_backed,
      nboot = nboot,
      backfiles = backfiles,
      descfiles = descfiles,
      backingpath = backingpath,
      rn = rn,
      cn = cn
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
  # impute_knn_brute skip these columns
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
        "manhattan" = 1L,
        "impute.knn" = 2L
      ),
      weighted = weighted,
      dist_pow = dist_pow,
      nboot = nboot,
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
      nboot = nboot,
      seed = seed,
      cores = cores
    )
  }
  # Handle NaN values
  imputed_values[is.nan(imputed_values)] <- NA
  # Map column indices from pre_imp_cols to original matrix columns and create
  # the index matrix for direct assignment
  imp_indices <- cbind(imputed_values[, 1], knn_indices[imputed_values[, 2]])
  # Fill in imputed values
  if (file_backed) {
    # Create bigmemory matrices and modify them directly
    result_list <- vector("list", nboot)
    for (i in seq_len(nboot)) {
      # Create bigmemory matrix initialized with original data
      bigmat <- bigmemory::big.matrix(
        nrow = nrow(obj),
        ncol = ncol(obj),
        type = "double",
        init = NA,
        backingpath = backingpath,
        backingfile = backfiles[i],
        descriptorfile = descfiles[i],
        dimnames = list(rn, cn)
      )
      # Copy original data to bigmatrix
      bigmat[, ] <- obj
      # Directly modify bigmatrix with imputed values using pre-computed indices
      bigmat[imp_indices] <- imputed_values[, i + 2]
      # Post-imputation directly on bigmatrix if needed
      if (post_imp && anyNA(bigmat[, subset, drop = FALSE])) {
        subset_data <- bigmat[, subset, drop = FALSE]
        na_indices <- which(is.na(subset_data), arr.ind = TRUE)
        sub_means <- colMeans(subset_data, na.rm = TRUE)
        i_vec <- na_indices[, 1]
        jj_vec <- na_indices[, 2]
        j_vec <- subset[jj_vec]
        # Update bigmatrix using vectorized indexing
        bigmat[cbind(i_vec, j_vec)] <- sub_means[jj_vec]
      }
      # Flush to ensure changes are written to disk
      bigmemory::flush(bigmat)
      result_list[[i]] <- bigmat
    }
    return(result_list)
  } else {
    # Create all copies upfront with replicate
    result_list <- replicate(nboot, obj, simplify = FALSE)
    for (i in seq_len(nboot)) {
      # Directly modify result_list[[i]] with imputed values using pre-computed indices
      result_list[[i]][imp_indices] <- imputed_values[, i + 2]
      # Post-imputation step if needed
      if (post_imp && anyNA(result_list[[i]][, subset, drop = FALSE])) {
        na_indices <- which(is.na(result_list[[i]][, subset, drop = FALSE]), arr.ind = TRUE)
        sub_means <- colMeans(result_list[[i]][, subset, drop = FALSE], na.rm = TRUE)
        i_vec <- na_indices[, 1]
        jj_vec <- na_indices[, 2]
        j_vec <- subset[jj_vec]
        result_list[[i]][cbind(i_vec, j_vec)] <- sub_means[jj_vec]
      }
    }
    return(result_list)
  }
}

# Helper function for knn_imp to create bigmatrix or copy for a single bootstrap.
# No need to pass to SlideKnn because knn_imp is always in memory there.
create_result_list <- function(
    data_to_copy,
    file_backed,
    nboot,
    backfiles = NULL,
    descfiles = NULL,
    backingpath = NULL,
    rn = NULL,
    cn = NULL) {
  result_list <- vector("list", nboot)

  if (file_backed) {
    for (i in seq_len(nboot)) {
      bigmat <- bigmemory::big.matrix(
        nrow = nrow(data_to_copy),
        ncol = ncol(data_to_copy),
        type = "double",
        init = NA,
        backingpath = backingpath,
        backingfile = backfiles[i],
        descriptorfile = descfiles[i],
        dimnames = list(rn, cn)
      )
      bigmat[, ] <- data_to_copy
      result_list[[i]] <- bigmat
    }
  } else {
    # For in-memory, create nboot copies
    for (i in seq_len(nboot)) {
      result_list[[i]] <- data_to_copy
    }
  }

  return(result_list)
}

#' Row Mean Imputation
#'
#' @description
#' Imputes missing values (NA) in a matrix by replacing them with the mean of their
#' respective rows. This is a simple imputation method that preserves row-wise
#' patterns in the data.
#'
#' @details
#' This function calculates the mean for each row excluding missing values and
#' replaces all missing values in that row with the computed mean. If a row
#' consists entirely of missing values, those values will remain as `NA` since
#' the mean cannot be calculated.
#'
#' This method is most appropriate when features (columns) are on similar scales
#' and when the assumption that missing values should be close to the row average
#' is reasonable for your data context.
#'
#' @param obj A numeric matrix with samples in rows and features in columns.
#'
#' @return A numeric matrix of the same dimensions as `obj` with missing values
#'   replaced by row means.
#'
#' @seealso [mean_impute_col()], [knn_imp()]
#'
#' @examples
#' # Create example matrix with missing values
#' mat <- matrix(c(1, 2, NA, 4, 5, 6, NA, 8, 9), nrow = 3)
#' print(mat)
#'
#' # Impute missing values with row means
#' imputed_mat <- mean_impute_row(mat)
#' print(imputed_mat)
#'
#' # Example with real data
#' data(khanmiss1)
#' # Note: khanmiss1 has genes in rows, so we use it directly
#' imputed_khan <- mean_impute_row(khanmiss1)
#' sum(is.na(khanmiss1)) # Original missing values
#' sum(is.na(imputed_khan)) # After imputation
#'
#' @export
mean_impute_row <- function(obj) {
  checkmate::assert_matrix(obj, mode = "numeric", .var.name = "obj")

  na_indices <- which(is.na(obj), arr.ind = TRUE)
  if (nrow(na_indices) == 0) {
    return(obj) # No missing values to impute
  }

  row_means <- rowMeans(obj, na.rm = TRUE)
  obj[na_indices] <- row_means[na_indices[, 1]]
  return(obj)
}

#' Column Mean Imputation
#'
#' @description
#' Imputes missing values (NA) in a matrix by replacing them with the mean of their
#' respective columns. This is a simple imputation method that preserves column-wise
#' patterns in the data.
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
#' @seealso [mean_impute_row()], [knn_imp()]
#'
#' @examples
#' # Create example matrix with missing values
#' mat <- matrix(c(1, 2, NA, 4, 5, 6, NA, 8, 9), nrow = 3)
#' colnames(mat) <- c("A", "B", "C")
#' print(mat)
#'
#' # Impute missing values with column means
#' imputed_mat <- mean_impute_col(mat)
#' print(imputed_mat)
#'
#' # Impute only specific columns by name
#' imputed_subset <- mean_impute_col(mat, subset = c("A", "C"))
#' print(imputed_subset)
#'
#' # Impute only specific columns by index
#' imputed_idx <- mean_impute_col(mat, subset = c(1, 3))
#' print(imputed_idx)
#'
#' # Example with real data
#' data(khanmiss1)
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

dummy <- function() {
  Rcpp::evalCpp()
}
