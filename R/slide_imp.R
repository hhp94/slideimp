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

#' Sliding Window k-NN Imputation
#'
#' @description
#' Performs k-nearest neighbor imputation on large numeric matrices using a sliding
#' window approach column-wise. This method assumes that columns are meaningfully sorted.
#'
#' @details
#' The sliding window approach divides the input matrix into smaller, overlapping
#' segments and applies k-NN imputation to each window independently. This strategy
#' is particularly advantageous for large data where applying k-NN imputation
#' to the entire matrix would be computationally prohibitive or exceed memory limits.
#'
#' Values in overlapping areas are averaged across windows to produce the final imputed result.
#' This approach assumes that features (columns) are ordered meaningfully (e.g.,
#' by genomic position, time series, etc.).
#'
#' For the underlying k-NN implementation details, see [knn_imp()].
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
#' @param cores Integer. Number of CPU cores for parallel distance computation. Default: 1.
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
#' @param seed Integer. Random seed for multiple imputation. Default: 42.
#' @param .progress Logical. Whether to display progress messages during imputation.
#'   Default: `FALSE`.
#' @param output Character. File path stem for saving file-backed big.matrix results
#'   (format: `"path/stem"`). If `NULL`, results are stored in memory. Recommended
#'   for large data and multiple imputations.
#' @param overwrite Logical. Whether to overwrite existing files at `output` path.
#'   Default: `FALSE`. See Notes for Windows.
#'
#' @return A list of [bigmemory::big.matrix()] objects with length `n_imp`. Each
#'   matrix has `nrow(obj)` rows and `length(subset)` columns with missing values
#'   imputed.
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
#' imputed_basic <- slide_imp(
#'   beta_matrix,
#'   k = 5,
#'   n_feat = 50,
#'   n_overlap = 10
#' )
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
#' imputed_big_in_memory <- slide_imp(
#'   big_data,
#'   k = 5,
#'   n_feat = 50,
#'   output = NULL
#' )
#'
#' # output is now set to a location to save. The backend is now file-backed matrix
#' # to handle large data
#' imputed_big <- slide_imp(
#'   big_data,
#'   k = 5,
#'   n_feat = 50,
#'   output = withr::local_tempfile()
#' )
#'
#' # Results may lose dimnames because bigmemory.allow.dimnames can't be set to TRUE.
#' # We can restore this with `restore_dimnames()`.
#' if (interactive()) {
#'   on.exit(options(bigmemory.allow.dimnames = getOption("bigmemory.allow.dimnames")), add = TRUE)
#'   options(bigmemory.allow.dimnames = TRUE)
#'   rownames(imputed_big[[1]]) <- NULL
#'   colnames(imputed_big[[1]]) <- NULL
#'   # We strip the dimnames manually as an example
#'   print(imputed_big)
#'   # Restore dimnames
#'   restore_dimnames(imputed_big)
#'   print(imputed_big)
#' }
#'
#' # ========================================
#' # Multiple Imputation
#' # ========================================
#'
#' # Predictive Mean Matching (PMM)
#' imputed_pmm <- slide_imp(
#'   beta_matrix,
#'   k = 8,
#'   n_feat = 60,
#'   n_imp = 3, # 3 imputations
#'   n_pmm = 5, # 5 donors for PMM. n_pmm = 0 enables neighbor bootstrapping.
#'   subset = c("feat1", "feat2"),
#'   output = withr::local_tempfile(),
#'   .progress = TRUE
#' )
#'
#' imputed_pmm
#'
#' @export
slide_imp <- function(
  obj,
  n_feat,
  subset = NULL,
  n_overlap = 10,
  rowmax = 0.9,
  colmax = 0.9,
  cores = 1,
  method = c("euclidean", "manhattan"),
  post_imp = TRUE,
  k = 10,
  dist_pow = 1,
  .progress = FALSE
) {
  # Pre-conditioning ----
  method <- match.arg(method)
  checkmate::assert_number(rowmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "rowmax")
  checkmate::assert_number(colmax, lower = 0, upper = 1, null.ok = FALSE, .var.name = "colmax")
  checkmate::assert_int(cores, lower = 1, null.ok = FALSE, .var.name = "cores")
  checkmate::assert_flag(post_imp, .var.name = "post_imp", null.ok = FALSE)
  checkmate::assert_flag(.progress, .var.name = ".progress", null.ok = FALSE)
  checkmate::assert(
    checkmate::check_character(subset, min.len = 1, any.missing = FALSE, unique = TRUE, null.ok = TRUE),
    checkmate::check_integerish(
      subset,
      lower = 1, upper = ncol(obj), min.len = 1, any.missing = FALSE, null.ok = TRUE, unique = TRUE
    ),
    .var.name = "subset"
  )
  checkmate::assert_int(n_feat, lower = 2, upper = ncol(obj), null.ok = FALSE, .var.name = "n_feat")
  checkmate::assert_int(n_overlap, lower = 0, upper = n_feat - 1, null.ok = FALSE, .var.name = "n_overlap")
  checkmate::assert_int(k, lower = 1, upper = n_feat - 1, null.ok = FALSE, .var.name = "k")
  checkmate::assert_matrix(obj)
  rn <- rownames(obj)
  cn <- colnames(obj)

  ## subset ----
  if (!is.null(subset)) {
    if (is.character(subset)) {
      stopifnot("`subset` are characters but `obj` doesn't have colnames" = !is.null(cn))
      matched <- match(subset, cn, nomatch = NA)
      subset <- matched[!is.na(matched)]
      subset <- sort(subset)
    }
    if (length(subset) == 0) {
      stop("`subset` is not found in `colnames(obj)`")
    }
  } else {
    subset <- seq_len(ncol(obj))
  }
  # Windowing Logic ----
  idx <- 1
  max_step <- ceiling((ncol(obj) - idx) / (n_feat - n_overlap))
  step <- 0:max_step
  start <- idx + (step * n_feat) - (step * n_overlap)
  end <- start + n_feat - 1
  # Overshoot
  n_overshoot <- sum(end > ncol(obj))
  corrected_length <- length(end) - n_overshoot
  start <- start[1:corrected_length]
  end <- end[1:corrected_length]
  end[corrected_length] <- ncol(obj)
  # Calculate where the subset lies
  subset_list <- lapply(seq_along(start), function(i) {
    first <- findInterval(start[i] - 1, subset) + 1
    last <- findInterval(end[i], subset)
    if (first <= last) {
      subset[first:last] - start[i] + 1
    } else {
      integer(0)
    }
  })
  # Calculate offsets
  width <- end - start + 1
  offset_start <- c(1, cumsum(width)[-length(width)] + 1)
  offset_end <- cumsum(width)
  # Sliding Imputation ----
  nr <- nrow(obj)
  nc <- ncol(obj)
  # Setup temporary directory if file-backed
  # Overlap regions to average over
  overlap <- find_overlap_regions(start, end)
  # Process each imputation sequentially
  # init
  # seed for this imputation
  # intermediate matrix for this imputation. In memory if not file backed
  intermediate <- bigmemory::big.matrix(
    nrow = nr,
    ncol = sum(width),
    type = "double",
    init = 0.0
  )
  if (.progress) {
    message("Step 1/3: Imputing")
  }
  # Sequential processing of windows
  for (i in seq_along(start)) {
    window_cols <- start[i]:end[i]
    # Get imputation results for this window
    imp_list <- knn_imp(
      obj = obj[, window_cols, drop = FALSE],
      k = k,
      colmax = colmax,
      cores = cores,
      method = method,
      post_imp = post_imp,
      dist_pow = dist_pow,
      subset = subset_list[[i]]
    )
    # Fill intermediate matrix with this window's results
    intermediate[, offset_start[i]:offset_end[i]] <- imp_list
  }
  # Create result matrix for this imputation
  result_imp <- bigmemory::big.matrix(
    nrow = nr,
    ncol = nc,
    type = "double",
    init = 0.0
  )
  ## Averaging ----
  if (.progress) {
    message("Step 2/3: Overlapping")
  }
  # Add the windows from intermediate
  bigmem_add_windows(
    result_imp@address,
    intermediate@address,
    start,
    end,
    offset_start,
    offset_end
  )
  if (.progress) {
    message("Step 3/3: Averaging")
  }
  # Average the overlaps
  bigmem_avg(
    result_imp@address,
    start = overlap$region[, "start"],
    end = overlap$region[, "end"],
    counts_vec = overlap$counts_vec,
    cores = cores
  )
  ## Post-imputation ----
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
  obj <- result_imp[, ]
  row.names(obj) <- rn
  colnames(obj) <- cn

  return(obj)
}
