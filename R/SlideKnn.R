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
#' @param n_imp Integer specifying the number of imputations to perform.
#'   Default is 1 for single imputation.
#' @param n_pmm Integer controlling the multiple imputation method when `n_imp` > 1.
#'   If `n_pmm` > 0: PMM multiple imputation using `n_pmm` closest donors (will not exceed
#'   available non-missing values column-wise).
#'   If `n_pmm` = 0: Bootstrap multiple imputation via resampling from k nearest neighbors.
#'   Ignored when `n_imp` = 1 (single imputation).
#' @param seed Integer; random seed for reproducible bootstrap sampling. Default is 42.
#' @param .progress Logical; if `TRUE`, show a progress bar. Default is `FALSE`.
#' @param output Character; path to save the output big.matrix if `obj` is file-backed.
#'   Required when `obj` is a [bigmemory::big.matrix()] or path to big.matrix descriptor.
#' @param overwrite Logical; if `TRUE` (default), overwrite existing files at `output`.
#' @param block Integer; block size for processing large matrices. If `NULL` (default),
#'   calculated automatically based on the matrix size and number of cores.
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
#' # Generate sample data with missing values with 10 samples and 200 columns
#' # where the column order is sorted (e.g., by genomic position or time)
#'
#' set.seed(1234)
#' beta_matrix <- t(sim_mat(200, 10)$input)
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
#' imputed_big <- SlideKnn(
#'   big_data,
#'   k = 5,
#'   n_feat = 50,
#'   output = withr::local_tempdir()
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
#' # ========================================
#' # File-backed Matrices (Large Datasets)
#' # ========================================
#'
#' # Create file-backed matrix for very large datasets
#' temp_dir <- withr::local_tempdir()
#' file_backed_data <- bigmemory::filebacked.big.matrix(
#'   nrow = nrow(beta_matrix),
#'   ncol = ncol(beta_matrix),
#'   type = "double",
#'   backingfile = "large_data.bin",
#'   descriptorfile = "large_data.desc",
#'   backingpath = temp_dir
#' )
#'
#' # Copy data to file-backed matrix
#' file_backed_data[, ] <- beta_matrix
#'
#' # Impute with automatic file output
#' imputed_filebacked <- SlideKnn(
#'   file_backed_data,
#'   k = 5,
#'   n_feat = 50,
#'   output = file.path(temp_dir, "imputed_result.bin"),
#'   overwrite = TRUE
#' )
#'
#' # Alternative: Load directly from descriptor file
#' desc_path <- file.path(temp_dir, "large_data.desc")
#' imputed_from_file <- SlideKnn(
#'   desc_path,
#'   k = 5,
#'   n_feat = 50,
#'   output = file.path(temp_dir, "imputed_from_desc.bin"),
#'   overwrite = TRUE
#' )
#'
#' # ========================================
#' # Multiple Imputation
#' # ========================================
#'
#' # Predictive Mean Matching (PMM) - recommended
#' imputed_pmm <- SlideKnn(
#'   bigmemory::as.big.matrix(beta_matrix),
#'   k = 8,
#'   n_feat = 60,
#'   n_imp = 3, # 3 imputations
#'   n_pmm = 5, # 5 donors for PMM
#'   output = temp_dir,
#'   overwrite = TRUE,
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
    n_pmm = 10,
    seed = 42,
    .progress = FALSE,
    output = NULL,
    overwrite = TRUE,
    block = NULL,
    strip_dimnames = FALSE) {
  # Pre-conditioning ----
  method <- match.arg(method)
  if (!is.character(obj)) {
    checkmate::assert_integerish(n_feat, lower = 2, upper = ncol(obj), len = 1, null.ok = FALSE, .var.name = "n_feat")
  }
  checkmate::assert_integerish(n_overlap, lower = 0, upper = n_feat - 1, len = 1, null.ok = FALSE, .var.name = "n_overlap")
  checkmate::assert_integerish(k, lower = 1, upper = n_feat - 1, len = 1, null.ok = FALSE, .var.name = "k")
  checkmate::assert_numeric(rowmax, lower = 0, upper = 1, len = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_numeric(colmax, lower = 0, upper = 1, len = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")
  checkmate::assert_integerish(n_imp, lower = 1, len = 1, null.ok = FALSE, .var.name = "n_imp")
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
  if (n_imp > 1 && n_pmm == 0 && weighted) {
    warning("If bootstrapping nearest neighbors, weighted will be forced to FALSE")
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
    checkmate::assert_integerish(n_feat, lower = 2, upper = ncol(obj), len = 1, null.ok = FALSE, .var.name = "n_feat")
    TRUE
  } else if (bigmemory::is.big.matrix(obj)) {
    TRUE
  } else if (is.matrix(obj) && is.numeric(obj)) {
    obj <- bigmemory::as.big.matrix(obj, type = "double")
    FALSE
  } else {
    stop("`obj` has to be a numeric matrix or a (descriptor of a) double bigmemory::big.matrix.")
  }
  if (n_imp > 1 && !file_backed) {
    warning("If n_imp > 1 should be used with output to bigmemory::big.matrix. Change `obj` to a bigmemory::big.matrix/path to bigmemory::big.matrix .desc file to enable this.")
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
  intermediate_list <- vector("list", n_imp)
  final_imputed_list <- vector("list", n_imp)
  counts_list <- vector("list", n_imp)

  ### file_backed ----
  if (file_backed) {
    checkmate::assert_flag(overwrite, .var.name = "overwrite", null.ok = FALSE)

    # Create temporary directory for intermediate files
    temp_dir <- withr::local_tempdir(pattern = paste0("SlideKnn_", Sys.getpid()))

    # For file-backed matrices, create separate files for each bootstrap
    for (boot_idx in seq_len(n_imp)) {
      # Check and handle output files
      if (n_imp > 1) {
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

        # Force garbage collection to release any file handles
        # This is especially important for memory-mapped files
        if (.Platform$OS.type == "windows") {
          invisible(gc(verbose = FALSE, full = TRUE))
        }

        # Try to delete
        tryCatch(
          {
            unlink(files_to_check, force = TRUE)
            if (.Platform$OS.type == "windows") {
              Sys.sleep(0.1)
              invisible(gc(verbose = FALSE, full = TRUE))
            }
          },
          warning = function(w) { }
        )

        # If files still exist after attempt, warn but continue
        if (any(file.exists(files_to_check))) {
          stop(sprintf(
            "Could not delete existing files before overwrite. Files: %s",
            paste(files_to_check[file.exists(files_to_check)], collapse = ", ")
          ))
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
    for (boot_idx in seq_len(n_imp)) {
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
          cores = 1L, # For Slide k-NN, fix cores = 1
          method = method,
          post_imp = post_imp,
          weighted = weighted,
          dist_pow = dist_pow,
          subset = subset_list[[i]],
          n_imp = n_imp,
          n_pmm = n_pmm,
          seed = seed,
          knn_imp = knn_imp,
          tree = tree,
          impute_knn_brute = impute_knn_brute,
          impute_knn_mlpack = impute_knn_mlpack,
          mean_impute_col = mean_impute_col
        )

        # Fill intermediate matrices for each bootstrap iteration
        for (boot_idx in seq_len(n_imp)) {
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
      n_imp = n_imp,
      n_pmm = n_pmm,
      seed = seed
    ),
    .progress = .progress
  )

  ## Averaging ----
  if (.progress) {
    message("Step 2/3: Overlapping")
  }

  # Process each bootstrap iteration
  for (boot_idx in seq_len(n_imp)) {
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
  for (boot_idx in seq_len(n_imp)) {
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
    Sys.sleep(0.1)
    invisible(gc(verbose = FALSE, full = TRUE))
  }

  ## post_imp ----
  if (post_imp) {
    if (.progress) {
      message("Post-imputation")
    }

    for (boot_idx in seq_len(n_imp)) {
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
    for (boot_idx in 1:n_imp) {
      rownames(final_imputed_list[[boot_idx]]) <- rn
      colnames(final_imputed_list[[boot_idx]]) <- cn
    }
  }

  # Return results
  if (file_backed) {
    return(final_imputed_list)
  } else {
    # Convert to regular matrices and restore names
    result_list <- vector("list", n_imp)
    for (boot_idx in seq_len(n_imp)) {
      out <- bigmemory::as.matrix(final_imputed_list[[boot_idx]])
      rownames(out) <- rn
      colnames(out) <- cn
      result_list[[boot_idx]] <- out
    }
    return(result_list)
  }
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
#' result <- knn_imp(t(khanmiss1), k = 5)[[1]]
#' sum(is.na(result)) # Should be 0
#'
#' # Detailed k-NN imputation with custom parameters
#' imputed <- knn_imp(
#'   obj = t(khanmiss1),
#'   k = 3,
#'   colmax = 0.5,
#'   rowmax = 0.8,
#'   method = "euclidean"
#' )[[1]]
#'
#' # Check results
#' imputed[1:5, 1:10]
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
#' # PMM imputation with 5 imputations and 3 donors. Increase `n_pmm` in real
#' # data to ensure realistic uncertainty quantification.
#' @examplesIf interactive()
#' imputed_pmm <- knn_imp(
#'   obj = t(khanmiss1),
#'   k = 5,
#'   n_imp = 5,
#'   n_pmm = 10,
#'   output = withr::local_tempdir(), # Enables results as list of bigmatrix
#'   seed = 123
#' )
#' # The default of the bigmatrix package is to not allow dimnames so the output
#' # lost the dimnames
#' imputed_pmm[[1]][, ][1:5, 1:5]
#'
#' # But this can be reassigned after setting the options to be TRUE like so
#' options(bigmemory.allow.dimnames = TRUE)
#' rownames(imputed_pmm[[1]]) <- rownames(big_data)
#' colnames(imputed_pmm[[1]]) <- colnames(big_data)
#' imputed_pmm[[1]][, ][1:5, 1:5]
#'
#' length(imputed_pmm) # Returns 5 imputed datasets
#'
#' # Bootstrap imputation for uncertainty injection
#' imputed_boot <- knn_imp(
#'   obj = t(khanmiss1),
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
    n_pmm = 10,
    seed = 42,
    output = NULL,
    overwrite = TRUE,
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
  checkmate::assert_string(output, null.ok = TRUE, .var.name = "output")
  checkmate::assert_logical(overwrite, len = 1, null.ok = FALSE, any.missing = FALSE, .var.name = "overwrite")

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
    boot_indices <- seq_len(n_imp)
    suffix <- if (n_imp == 1) "" else paste0("_boot", boot_indices)
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
    boot_indices <- seq_len(n_imp)
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
      n_imp = n_imp,
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
  # Handle NaN values
  imputed_values[is.nan(imputed_values)] <- NA
  # Map column indices from pre_imp_cols to original matrix columns and create
  # the index matrix for direct assignment
  imp_indices <- cbind(imputed_values[, 1], knn_indices[imputed_values[, 2]])
  # Fill in imputed values
  if (file_backed) {
    # Create bigmemory matrices and modify them directly
    result_list <- vector("list", n_imp)
    for (i in seq_len(n_imp)) {
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
    result_list <- replicate(n_imp, obj, simplify = FALSE)
    for (i in seq_len(n_imp)) {
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
    n_imp,
    backfiles = NULL,
    descfiles = NULL,
    backingpath = NULL,
    rn = NULL,
    cn = NULL) {
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
        dimnames = list(rn, cn)
      )
      bigmat[, ] <- data_to_copy
      result_list[[i]] <- bigmat
    }
  } else {
    # For in-memory, create n_imp copies
    for (i in seq_len(n_imp)) {
      result_list[[i]] <- data_to_copy
    }
  }

  return(result_list)
}

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
