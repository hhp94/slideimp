#' Sliding KNN Imputation
#'
#' Performs sliding window KNN imputation on a numeric matrix to handle missing values.
#' The matrix is divided into overlapping windows, and imputation is applied to each window.
#' Overlapping regions are averaged to produce the final imputed matrix.
#'
#' @details
#' See [knn_imp()] for details about the implementation.
#'
#' If your `obj` is a [bigmemory::big.matrix()] or description file, you have to set
#' `strip_dimnames` to `TRUE` for the output big.matrix to have the same dimnames as
#' the `obj`. Otherwise, you can also re-add the dimnames to the output using the
#' dimnames of the original object after setting `options(bigmemory.allow.dimnames = TRUE)`.
#'
#' @param obj A numeric matrix with **samples in rows** and **features in columns**. See details.
#' Can also be a path to the description file of, or a [bigmemory::big.matrix()].
#' @param n_feat Integer specifying the number of features (columns) in each window. Must be between 2 and the number of columns in `obj`.
#' @param subset Character vector of column names or integer vector of column indices specifying the subset of columns to perform imputation.
#' @param n_overlap Integer specifying the number of features to overlap between consecutive windows. Default is 10. Must be between 0 and `n_feat - 1`.
#' @param k Integer specifying the number of nearest neighbors to use for imputation. Must be between 1 and (number of columns - 1).
#' @param rowmax Numeric between 0 and 1 specifying the maximum allowable proportion of missing values in any row. If exceeded, the function stops with an error.
#' @param colmax Numeric between 0 and 1 specifying the threshold proportion of missing values in a column above which the column is imputed using the mean instead of k-NN.
#' @param cores Integer specifying the number of cores to use for parallel computation of distances.
#' @param method Character string specifying the distance metric for k-NN. One of `"euclidean"`, `"manhattan"`, or `"impute.knn"`. Defaults to `"euclidean"`.
#' @param tree Character string specifying the k-NN method. `NULL` (default) uses naive search.
#' `"kd"` uses KDTree and `"ball"` uses BallTree as implemented by the mlpack package where missing values are first filled with column means (biased at high percentage missing).
#' @param post_imp Logical; if TRUE (default), any missing after k-NN imputations will be imputed with [mean_impute_col()].
#' @param weighted Logical; controls for the imputed value to be a simple mean or weighted mean by inverse distance.
#' @param dist_pow  A positive double that controls the penalty for larger distances in the weighted mean imputation.
#' Must be greater than zero: values between 0 and 1 apply a softer penalty, 1 is linear (default), and values greater than 1 apply a harsher penalty.
#' @param .progress Logical; if TRUE, show a progress bar. Default is FALSE.
#' @param output Character; path to save the output big.matrix if `obj` is file-backed. Required in that case.
#' @param overwrite Logical; if TRUE (default), overwrite existing files at `output`.
#' @param block Integer; block size for processing large matrices. If NULL (default), calculated automatically.
#' @param strip_dimnames Logical; if FALSE (default), dimnames will not be removed which will increase memory usage.
#' Should set to TRUE to save memory from overhead especially when cores > 1. See details.
#'
#' @return A numeric matrix/big.matrix of the same dimensions as `obj` with missing values imputed.
#'
#' @examples
#'
#' data(khanmiss1)
#'
#' # Set `strip_dimnames` to `TRUE` with `options(bigmemory.allow.dimnames = TRUE)`
#' # when running multiple cores to minimize overhead.
#' imputed <- SlideKnn(t(khanmiss1), k = 10, n_feat = 100, n_overlap = 10)
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
    method = c("euclidean", "manhattan", "impute.knn"),
    tree = NULL,
    post_imp = TRUE,
    weighted = FALSE,
    dist_pow = 1,
    .progress = FALSE,
    output = NULL,
    overwrite = TRUE,
    block = NULL,
    strip_dimnames = FALSE) {
  # Pre-conditioning ----
  method <- match.arg(method)
  checkmate::assert_integerish(n_feat, lower = 2, upper = ncol(obj), len = 1, null.ok = FALSE, .var.name = "n_feat")
  checkmate::assert_integerish(n_overlap, lower = 0, upper = n_feat - 1, len = 1, null.ok = FALSE, .var.name = "n_overlap")
  checkmate::assert_integerish(k, lower = 1, upper = ncol(obj) - 1, len = 1, null.ok = FALSE, .var.name = "k")
  checkmate::assert_numeric(rowmax, lower = 0, upper = 1, len = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_numeric(colmax, lower = 0, upper = 1, len = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_integerish(cores, lower = 1, len = 1, null.ok = FALSE, .var.name = "cores")
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
  # Strategy: use big.memory and pass the pointer around across cores to avoid
  # copying the matrix over the cores to minimize overhead
  ### file_backed ----
  if (file_backed) {
    checkmate::assert_flag(overwrite, .var.name = "overwrite", null.ok = FALSE)
    files_to_check <- file.path(backingpath, c(backfile, descfile))
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
    temp_dir <- withr::local_tempdir(pattern = paste0("SlideKnn_", Sys.getpid()))
    # Description of the main matrices:
    # 1) `intermediate` is the matrix that does not contain overlap. This is so that
    # we can parallelize the impute_knn function, which takes a lot of time.
    # Do this to avoid race conditions in parallelization
    # 2) `final_imputed` holds the numerator value of imputation after adding the
    # overlapping windows
    # 3) `counts` is the column vector (ncol * 1) that holds the denominator to
    # normalize the column of final_imputed by
    # 4) `final_imputed` holds the averaged results
    intermediate <- bigmemory::filebacked.big.matrix(
      nrow = nr,
      ncol = sum(width),
      type = "double",
      init = 0.0,
      backingfile = "intermediate.bin",
      descriptorfile = "intermediate.desc",
      backingpath = temp_dir
    )
    # final_imputed and counts should be fast enough to be ran sequentially
    counts <- bigmemory::filebacked.big.matrix(
      nrow = nc,
      ncol = 1,
      type = "double",
      init = 0.0,
      backingfile = "counts.bin",
      descriptorfile = "counts.desc",
      backingpath = temp_dir
    )
    final_imputed <- bigmemory::filebacked.big.matrix(
      nrow = nr,
      ncol = nc,
      type = "double",
      init = 0.0,
      backingpath = backingpath,
      backingfile = backfile,
      descriptorfile = descfile
    )
  } else {
    ### in-memory ----
    intermediate <- bigmemory::big.matrix(
      nrow = nr,
      ncol = sum(width),
      type = "double",
      init = 0.0
    )
    # final_imputed and counts should be fast enough to be ran sequentially
    counts <- bigmemory::big.matrix(
      nrow = nc,
      ncol = 1,
      type = "double",
      init = 0.0
    )
    final_imputed <- bigmemory::big.matrix(
      nrow = nr,
      ncol = nc,
      type = "double",
      init = 0.0
    )
  }
  # These are pointers to the big.matrices. We are passing these around in
  # the for loops. These also avoid the fragile <<- solution
  obj_desc <- bigmemory::describe(obj)
  intermediate_desc <- bigmemory::describe(intermediate)
  counts_desc <- bigmemory::describe(counts)
  final_imputed_desc <- bigmemory::describe(final_imputed)
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
        intermediate_big <- bigmemory::attach.big.matrix(intermediate_desc)
        imp <- impute_knn(
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
          knn_imp = knn_imp,
          tree = tree,
          impute_knn_naive = impute_knn_naive,
          impute_knn_mlpack = impute_knn_mlpack,
          mean_impute_col = mean_impute_col
        )
        intermediate_big[, offset_start[i]:offset_end[i]] <- imp
      },
      impute_knn = impute_knn,
      subset_list = subset_list,
      knn_imp = knn_imp,
      weighted = weighted,
      dist_pow = dist_pow,
      impute_knn_naive = impute_knn_naive,
      mean_impute_col = mean_impute_col,
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
      offset_end = offset_end
    ),
    .progress = .progress
  )
  # Then we sequentially fill in the value to avoid race condition. This step
  # is fast enough without parallel so we don't need to pass the pointers.
  ## Averaging ----
  if (.progress) {
    message("Step 2/3: Overlapping")
  }
  purrr::walk(
    seq_along(start),
    function(i) {
      window_cols <- start[i]:end[i]
      final_imputed[, window_cols] <- {
        final_imputed[, window_cols] + intermediate[, offset_start[i]:offset_end[i]]
      }
      counts[window_cols] <- counts[window_cols] + 1
    },
    .progress = FALSE
  )
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

  # A bit of extra calculations here because subset won't apply. But its fine
  # because its fast enough. Otherwise the codes would be too verbose
  purrr::walk(
    seq_along(w_start),
    fn(
      function(i) {
        window_cols <- w_start[i]:w_end[i]
        counts_big <- bigmemory::attach.big.matrix(counts_desc)
        final_imputed_big <- bigmemory::attach.big.matrix(final_imputed_desc)
        average <- counts_big[window_cols] > 1
        # part 1 if average == FALSE (average == 1), then just assign `final_imputed` to `final_imputed`
        # final_imputed_big[, window_cols[!average]] <- final_imputed_big[, window_cols[!average], drop = F]
        # part 2 if average == TRUE (average == 1), then just sweep `final_imputed` to `final_imputed`
        if (length(window_cols[average]) > 0) {
          final_imputed_big[, window_cols[average]] <- sweep(
            final_imputed_big[, window_cols[average], drop = F],
            MARGIN = 2,
            STATS = counts_big[window_cols[average]],
            FUN = "/"
          )
        }
      },
      counts_desc = counts_desc,
      final_imputed_desc = final_imputed_desc,
      w_start = w_start,
      w_end = w_end
    ),
    .progress = FALSE
  )
  ## post_imp ----
  if (post_imp) {
    if (.progress) {
      message("Post-imputation")
    }
    purrr::walk(
      seq_along(w_start),
      fn(
        function(i) {
          window_cols <- w_start[i]:w_end[i]
          final_imputed_big <- bigmemory::attach.big.matrix(final_imputed_desc)
          if (anyNA(final_imputed_big[, window_cols])) {
            final_imputed_big[, window_cols] <- mean_impute_col(
              final_imputed_big[, window_cols, drop = FALSE],
              subset = w_subset_list[[i]]
            )
          }
        },
        mean_impute_col = mean_impute_col,
        final_imputed_desc = final_imputed_desc,
        w_subset_list = w_subset_list,
        w_start = w_start,
        w_end = w_end
      ),
      .progress = FALSE
    )
  }
  # Restore names
  if (strip_dimnames) {
    rownames(obj) <- rn
    colnames(obj) <- cn
    rownames(final_imputed) <- rn
    colnames(final_imputed) <- cn
  }

  if (file_backed) {
    return(final_imputed)
  } else {
    out <- bigmemory::as.matrix(final_imputed)
    rownames(out) <- rn
    colnames(out) <- cn
    return(out)
  }
}

#' KNN Imputation Wrapper
#'
#' A wrapper function for KNN imputation that filters rows based on missing value proportions before imputing.
#' Rows with missing proportions exceeding `rowmax` are not imputed in this step to avoid
#' throwing exceptions in [knn_imp()].
#'
#' @inheritParams SlideKnn
#' @param knn_imp function. Needed for crating.
#' @param ... Pass through to knn_imp.
#'
#' @return The imputed matrix, with only qualifying rows imputed via KNN; others remain unchanged.
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
    ...) {
  na_mat <- is.na(obj)

  # Determine 'good_rows': rows where the proportion of NAs is less than 'rowmax'.
  # Rows with too many NAs (as defined by rowmax) are excluded from imputation.
  good_rows <- rowSums(na_mat) / ncol(na_mat) < rowmax

  # If no rows meet the criteria for imputation, return the original object as is.
  if (sum(good_rows) == 0) {
    return(obj)
  }

  imputed_good <- knn_imp(
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
    ...
  )

  # Initialize the result matrix with the original object's values.
  result <- obj
  # Place the imputed values from 'imputed_good' back into the corresponding 'good_rows'
  # in the result matrix.
  result[good_rows, ] <- imputed_good
  return(result)
}

#' K-Nearest Neighbor (k-NN) Imputation
#'
#' This function imputes missing values in a numeric matrix using the k-Nearest
#' Neighbors algorithm. It follows a two-stage process. First, it imputes
#' columns with a proportion of missing values below `colmax` using k-NN.
#' Second, if requested, any remaining missing values are imputed using the column mean.
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
#' @note Compared to `impute::impute.knn`, for columns with very high missing, the
#' mean imputation uses the imputed values and original values for the mean calculation
#' instead of just the original values.
#'
#' @inheritParams SlideKnn
#' @param obj A numeric matrix with **samples in rows** and **features in columns**.
#' @param ... Not Implemented.
#'
#' @return A numeric matrix of the same dimensions as `obj` with missing values imputed.
#'
#' @export
#'
#' @examples
#' # See ?khanmiss1
#' data(khanmiss1)
#' sum(is.na(khanmiss1))
#'
#' # Perform k-NN imputation. `khanmiss1` stores genes in the row so we have to t().
#' # set method to "impute.knn" to mimic how distant is scaled in impute::impute.knn.
#' imputed <- knn_imp(obj = t(khanmiss1), k = 3, colmax = 0.5, rowmax = 0.8, method = "euclidean")
#' imputed[1:5, 1:20]
#' sum(is.na(imputed))
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
  if (!is.null(tree)) {
    if (method == "impute.knn") {
      stop("When tree is not NULL, method cannot be 'impute.knn'")
    }
  }
  if (!is.null(subset)) {
    if (length(subset) == 0) {
      # message("non-NULL `subset` of length 0 detected. Returning object unimputed.")
      return(obj)
    }
    if (is.character(subset)) {
      stopifnot("`subset` are characters but `obj` doesn't have colnames" = !is.null(colnames(obj)))
      matched <- match(subset, colnames(obj), nomatch = NA)
      # subset converted to index
      subset <- matched[!is.na(matched)]
    }
  } else {
    # If subset is NULL, set it to all columns
    subset <- seq_len(ncol(obj))
  }
  complement <- setdiff(seq_len(ncol(obj)), subset)
  if (cores >= 16) {
    warning("Setting cores too high can slow down runtime. Benchmark your data first.")
  }
  miss <- is.na(obj)
  if (sum(miss[, subset, drop = FALSE]) == 0) {
    message("No missing data")
    return(obj)
  }
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
  # Set all values outside of subset to be zero cmiss. This will make impute_knn_naive skip these columns
  pre_imp_cmiss[pos_complement] <- 0L
  if (any(rowSums(pre_imp_miss) / ncol(pre_imp_cols) == 1)) {
    stop("Row(s) missing exceeded rowmax. Remove row(s) with too high NA %")
  }
  # Impute
  ## knn imp cols. Note: only pre_imp_cols is imputed if post_imp is FALSE.
  if (is.null(tree)) {
    imputed_values <- impute_knn_naive(
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
      cores = cores
    )
  }
  # Impute `pre_imp_cols`
  imputed_values[is.nan(imputed_values)] <- NA
  pre_imp_cols[imputed_values[, 1]] <- imputed_values[, 2]

  ## Reconstruct obj
  obj[, knn_imp_cols] <- pre_imp_cols

  if (post_imp && anyNA(obj)) {
    obj <- mean_impute_col(obj)
  }

  return(obj)
}

#' @title Row Mean Imputation
#'
#' @description
#' Imputes missing values (NA) in a matrix by replacing them with the mean of their
#' respective rows.
#'
#' @inheritParams SlideKnn
#'
#' @return The matrix with NA values replaced by row means.
#'
#' @export
mean_impute_row <- function(obj) {
  na_indices <- which(is.na(obj), arr.ind = TRUE)
  row_means <- rowMeans(obj, na.rm = TRUE)
  obj[na_indices] <- row_means[na_indices[, 1]]
  return(obj)
}

#' @title Col Mean Imputation
#'
#' @description
#' Imputes missing values (NA) in a matrix by replacing them with the mean of their
#' respective cols.
#'
#' @inheritParams SlideKnn
#'
#' @return The matrix with NA values replaced by col means.
#'
#' @export
mean_impute_col <- function(obj, subset = NULL) {
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
      stopifnot("`subset` are characters but `obj` doesn't have colnames" = !is.null(colnames(obj)))
      matched <- match(subset, colnames(obj), nomatch = NA)
      subset <- matched[!is.na(matched)]
    }
    obj_subset <- obj[, subset, drop = FALSE]
    na_indices <- which(is.na(obj_subset), arr.ind = TRUE)
    if (length(na_indices) > 0) {
      column_means <- colMeans(obj_subset, na.rm = TRUE)
      obj_subset[na_indices] <- column_means[na_indices[, 2]]
      obj[, subset] <- obj_subset
    }
  } else {
    na_indices <- which(is.na(obj), arr.ind = TRUE)
    if (length(na_indices) > 0) {
      column_means <- colMeans(obj, na.rm = TRUE)
      obj[na_indices] <- column_means[na_indices[, 2]]
    }
  }
  return(obj)
}

dummy <- function() {
  Rcpp::evalCpp()
  mlpack::adaboost
}
