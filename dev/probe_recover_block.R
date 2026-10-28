# Is one dgemm over a contiguous column block of Xhat cheaper than the dgemm
# over all of Xhat? Before the block restriction, recover_factors
# (src/svd_triplet.h:117-135) formed V = Xhat' * Vn over every column each EM
# iteration, but the loop only reads the rows for the changing columns, which
# the permutation puts at the tail of Xhat (src/armaSVD.cpp:482-490).
# Restricting it is still a single dgemm, on A = Xhat.colptr(n_fixed) with
# lda = nrX - the same shape reconstruct_cols() already uses
# (src/armaSVD.cpp:276-303); gemm_tn_block() (src/pca_linalg_utils.h) is
# that call, and this probe is the measurement that justified it.
#
# A column-major tail block of columns is contiguous with leading dimension
# nrX, so a separate N_ROW x n_chg matrix has the same memory layout and the
# same dgemm call. What this does not reproduce: the kernel would write into
# rows of the taller V (ldc = n_elig), not a fresh matrix; that output is
# n_chg x K and small next to the input.
#
# options(matprod = "blas") skips R's NaN/Inf pre-scan so crossprod() is the
# bare dgemm('T', 'N') call. OPENBLAS_NUM_THREADS is set before R starts by
# dev/probe_recover_block.sh.
#
# Usage, from the repository root (BENCH_OUT names an output directory):
#   Rscript dev/probe_recover_block.R run <label>
#   Rscript dev/probe_recover_block.R compare
if (!file.exists("DESCRIPTION")) stop("run from the repository root")

args <- commandArgs(trailingOnly = TRUE)
part <- args[[1L]]
out_dir <- Sys.getenv("BENCH_OUT")
stopifnot(nzchar(out_dir), dir.exists(out_dir))

N_ROW <- 1000L
N_COL <- 40000L
K <- 11L # ncp + 1 under method = "regularized", ncp = 10
N_CHG <- 4000L # 10% changing, as in dev/bench_wide_phases.R
MIN_BATCH <- 0.3 # seconds of wall time per timed batch
N_BATCH <- 5L

if (part == "run") {
  label <- args[[2L]]
  cat(sprintf(
    "%s: OPENBLAS_NUM_THREADS=%s\n",
    label, Sys.getenv("OPENBLAS_NUM_THREADS", "-")
  ))
  options(matprod = "blas")

  set.seed(20261028)
  x <- matrix(rnorm(N_ROW * N_COL), N_ROW)
  vn <- matrix(rnorm(N_ROW * K), N_ROW)
  x_chg <- x[, (N_COL - N_CHG + 1L):N_COL] # the contiguous tail block

  # the block result must equal the tail rows of the full result
  v_full <- crossprod(x, vn)
  v_chg <- crossprod(x_chg, vn)
  ref <- v_full[(N_COL - N_CHG + 1L):N_COL, ]
  rel_keep <- abs(ref) >= 1e-3
  stopifnot(sum(rel_keep) > 0L)
  cat(sprintf(
    "  block vs tail of full: max abs %.3g, max rel %.3g over %d of %d entries, bitwise %s\n",
    max(abs(v_chg - ref)), max(abs(v_chg - ref)[rel_keep] / abs(ref[rel_keep])),
    sum(rel_keep), length(ref), identical(v_chg, ref)
  ))

  cases <- list(full = x, block = x_chg)
  res <- list()
  for (nm in names(cases)) {
    a <- cases[[nm]]
    # size the batch from a loop long enough to time: a single block call is
    # below the elapsed-timer resolution on Windows and reads as 0
    reps <- 1L
    repeat {
      t_cal <- system.time(for (r in seq_len(reps)) crossprod(a, vn))[["elapsed"]]
      if (t_cal >= 0.05) break
      reps <- reps * 2L
    }
    reps <- max(1L, as.integer(ceiling(reps * MIN_BATCH / t_cal)))
    for (b in seq_len(N_BATCH)) {
      t0 <- proc.time()
      for (r in seq_len(reps)) crossprod(a, vn)
      dt <- proc.time() - t0
      res[[length(res) + 1L]] <- data.frame(
        label = label, case = nm, cols = ncol(a), batch = b, reps = reps,
        wall = dt[["elapsed"]], cpu = dt[["user.self"]] + dt[["sys.self"]]
      )
    }
    cat(sprintf("  %-5s (%5d columns): %d batches x %d calls\n", nm, ncol(a), N_BATCH, reps))
  }
  write.csv(
    do.call(rbind, res),
    file.path(out_dir, sprintf("recover_%s.csv", label)),
    row.names = FALSE
  )
}

if (part == "compare") {
  d <- do.call(rbind, lapply(
    list.files(out_dir, pattern = "^recover_.*\\.csv$", full.names = TRUE),
    read.csv
  ))
  d$ms_per_call <- 1000 * d$wall / d$reps
  agg <- do.call(rbind, lapply(split(d, list(d$label, d$case), drop = TRUE), function(x) {
    data.frame(
      label = x$label[1], case = x$case[1], cols = x$cols[1],
      calls = sum(x$reps), ms_per_call = mean(x$ms_per_call),
      se_ms = sd(x$ms_per_call) / sqrt(nrow(x)),
      cpu_over_wall = sum(x$cpu) / sum(x$wall)
    )
  }))
  agg <- agg[order(agg$label, agg$case), ]
  print(format(agg, digits = 3), row.names = FALSE)

  for (lb in sort(unique(agg$label))) {
    f <- agg[agg$label == lb & agg$case == "full", ]
    b <- agg[agg$label == lb & agg$case == "block", ]
    cat(sprintf(
      "%s: block / full time %.3f (column ratio %.3f), full - block %.1f ms (%.1f SE)\n",
      lb, b$ms_per_call / f$ms_per_call, b$cols / f$cols,
      f$ms_per_call - b$ms_per_call,
      (f$ms_per_call - b$ms_per_call) / sqrt(f$se_ms^2 + b$se_ms^2)
    ))
  }
}
