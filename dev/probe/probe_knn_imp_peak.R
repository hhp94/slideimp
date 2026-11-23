# Probe: knn_imp()'s peak memory, for a source tree given as an argument, so
# the same stages can be read before and after a change (Phase 7 of
# dev/si_mi_conversion.md at d5a34e5 moved its write from `obj[loc] <- value`
# into the matrix door). Each stage runs in a FRESH R process (callr),
# because a process's peak working set only ever grows.
#
# The DLL of <root> is built once here with load_all1(timer = FALSE), from
# this checkout's dev/dev-utils.R; every child loads that same DLL with
# pkgload::load_all(compile = FALSE), so no stage is measured on another
# build. A stage the tree's knn_imp() cannot run (an argument it does not
# have) is reported as such.
#
# Peak is the process peak working set (ps, Windows), which counts the
# kernel's C++ allocations that R's gc() figures miss. Each stage reports
#   extra = (peak after the stage - working set before it) / size of x
# so 1.00 means the stage held one more copy of the input on top of it. The
# input is built block by block in place; `masked` says when the stage
# stayed under the build's own peak.
#
# Run from the package root, naming the tree to load:
#   Rscript dev/probe/probe_knn_imp_peak.R <root>
# e.g. `.` for this checkout, or a worktree of another commit.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, dir.exists(args[[1]]))
root <- normalizePath(args[[1]], winslash = "/")
here <- getwd()
source(file.path(here, "dev", "dev-utils.R"))
setwd(root)
load_all1(timer = FALSE)
setwd(here)
cat(
  "tree:",
  root,
  " HEAD:",
  system2("git", c("-C", root, "rev-parse", "--short", "HEAD"), stdout = TRUE),
  " dirty:",
  length(system2(
    "git",
    c("-C", root, "status", "--porcelain", "--", "R", "src"),
    stdout = TRUE
  )) >
    0L,
  "\n"
)

n <- 2000L
p <- 4000L
frac_mean <- 0.1 # columns over colmax, so post_imp mean-fills them
cores <- 4L

stage <- function(root, scenario, n, p, frac_mean, cores) {
  suppressMessages(pkgload::load_all(root, compile = FALSE, quiet = TRUE))
  mem <- function() ps::ps_memory_full_info(ps::ps_handle())
  set.seed(1)
  ## The input, low rank plus noise, written block by block into one matrix.
  r <- 3L
  U <- matrix(rnorm(n * r), n, r)
  V <- matrix(rnorm(p * r), p, r)
  x <- matrix(0, n, p)
  for (j1 in seq(1L, p, by = 1000L)) {
    j <- j1:min(j1 + 999L, p)
    x[, j] <- U %*%
      t(V[j, , drop = FALSE]) +
      matrix(rnorm(n * length(j), sd = 0.5), n, length(j))
  }
  dimnames(x) <- list(paste0("r", seq_len(n)), paste0("c", seq_len(p)))
  ## Holes: 5% of the rows in every other column, and 95% of the rows in the
  ## last `frac_mean` of the columns (over the default colmax of 0.9).
  inel <- seq.int(p - as.integer(p * frac_mean) + 1L, p)
  for (j in seq_len(p)) {
    k <- if (j %in% inel) {
      as.integer(0.95 * n)
    } else if (j %% 2L == 0L) {
      as.integer(0.05 * n)
    } else {
      0L
    }
    if (k > 0L) {
      x[sort(sample.int(n, k)), j] <- NA
    }
  }
  rm(U, V)
  invisible(gc())
  X <- as.numeric(object.size(x))
  ## A tenth of the columns, two in every twenty: an odd column (no hole)
  ## and an even one (holes), eligible below p - p * frac_mean and
  ## ineligible above, so the K-NN part and the mean part are both cut.
  sel <- sort(c(seq(1L, p, by = 20L), seq(2L, p, by = 20L)))
  call_args <- switch(
    scenario,
    "post_imp = TRUE" = list(),
    "post_imp = FALSE" = list(post_imp = FALSE),
    "subset = 10%" = list(subset = sel),
    "output = mipca" = list(output = "mipca"),
    "subset 10%, mipca" = list(subset = sel, output = "mipca")
  )
  if (!all(names(call_args) %in% names(formals(knn_imp)))) {
    return(list(X_mb = X / 2^20, extra = NA_real_, masked = NA, ran = FALSE))
  }
  before <- mem()
  t0 <- proc.time()[["elapsed"]]
  res <- do.call(knn_imp, c(list(x, k = 10, cores = cores), call_args))
  secs <- proc.time()[["elapsed"]] - t0
  after <- mem()
  ## ps returns a named numeric vector, not a list.
  list(
    X_mb = X / 2^20,
    extra = (after[["peak_wset"]] - before[["wset"]]) / X,
    masked = after[["peak_wset"]] <= before[["peak_wset"]],
    ran = TRUE,
    secs = secs
  )
}

scenarios <- c(
  "post_imp = TRUE",
  "post_imp = FALSE",
  "subset = 10%",
  "output = mipca",
  "subset 10%, mipca"
)
cat(sprintf(
  "input %d x %d, %.0f%% of columns over colmax, k = 10, cores = %d\n",
  n,
  p,
  100 * frac_mean,
  cores
))
for (s in scenarios) {
  out <- callr::r(
    stage,
    args = list(
      root = root,
      scenario = s,
      n = n,
      p = p,
      frac_mean = frac_mean,
      cores = cores
    )
  )
  if (!out$ran) {
    cat(sprintf("  %-20s not in this tree's knn_imp()\n", s))
  } else {
    cat(sprintf(
      "  %-20s extra %.2f x  (input %.0f MB, %.1f s)%s\n",
      s,
      out$extra,
      out$X_mb,
      out$secs,
      if (isTRUE(out$masked)) "  MASKED by the build's peak" else ""
    ))
  }
}
