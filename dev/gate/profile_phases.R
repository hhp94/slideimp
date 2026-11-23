# Where the time of one imputer call goes, by the phase timers of the
# kernels, on a release build installed by dev/gate/profile_install.R. Dev only:
#   Rscript dev/gate/profile_phases.R <lib> <chr21.qs2> [n_calls] [sim_n] [sim_p]
# <lib> must be first on the library path (R_LIBS) before R starts; this
# script checks that slideimp came from it and stops otherwise.
#
# Workloads, those of dev/gate/eigen_time_worker.R: pca_imp(ncp = 5) under the
# exact solver and under LOBPCG, on chr21 (without its all-NA probes) and on
# sim_mat(sim_n, sim_p) under set.seed(1), and knn_imp(k = 10) on both. One
# untimed call first, then n_calls calls.
#
# PCA: rcpptimer writes a summary named after the timer into the global
# environment when the kernel returns. Each phase's time is summed over the
# calls and printed as a share of "pca_imp_internal_cpp_total"; "untimed" is
# that total minus the phases of the kernel's own level, so it holds
# whatever no phase covers. "form_gram", "eig" and "recover_factors" are
# inside "svd" and are printed as parts of it, not added to it.
# K-NN: the kernel times its parallel loop only, so what is printed is the
# wall time of knn_imp() from R against that loop: the rest is the R side
# plus the kernel's setup, not told apart.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "usage: Rscript dev/gate/profile_phases.R <lib> <chr21.qs2> [n_calls] ",
    "[sim_n] [sim_p]",
    call. = FALSE
  )
}
lib <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
chr21_path <- args[[2]]
n_calls <- if (length(args) >= 3L) as.integer(args[[3]]) else 5L
sim_n <- if (length(args) >= 4L) as.integer(args[[4]]) else 200L
sim_p <- if (length(args) >= 5L) as.integer(args[[5]]) else 5000L

suppressPackageStartupMessages(library(slideimp))
from <- normalizePath(dirname(find.package("slideimp")), winslash = "/")
if (!identical(tolower(from), tolower(lib))) {
  stop("slideimp loaded from ", from, ", not ", lib, call. = FALSE)
}

chr21 <- t(qs2::qs_read(chr21_path))
chr21 <- chr21[, colSums(!is.na(chr21)) > 0L]
set.seed(1)
sim <- sim_mat(n = sim_n, p = sim_p)$input
data <- list(chr21 = chr21, sim = sim)
for (nm in names(data)) {
  x <- data[[nm]]
  cat(sprintf(
    "DATA %s: %d x %d, %.2f%% NA\n",
    nm,
    nrow(x),
    ncol(x),
    100 * mean(is.na(x))
  ))
}

# One summary as a named vector of total microseconds per phase.
phase_totals <- function(tm) {
  if (!is.data.frame(tm) || !all(c("Microseconds", "Count") %in% names(tm))) {
    utils::str(tm)
    stop("a timer summary of another shape than expected", call. = FALSE)
  }
  stats::setNames(tm$Microseconds * tm$Count, rownames(tm))
}

take_timer <- function(name) {
  if (!exists(name, envir = globalenv(), inherits = FALSE)) {
    stop("no timer summary named ", name, ": not a LOC_TIMER build?")
  }
  tm <- get(name, envir = globalenv())
  rm(list = name, envir = globalenv())
  tm
}

inner <- c("form_gram", "eig", "recover_factors")
total_tag <- "pca_imp_internal_cpp_total"

profile_pca <- function(label, x, solver) {
  run <- function() pca_imp(x, ncp = 5, solver = solver)
  set.seed(1)
  warm <- run()
  stopifnot(is.matrix(warm), !anyNA(warm))
  first <- take_timer("pca_kernel")
  sums <- NULL
  counts <- NULL
  for (i in seq_len(n_calls)) {
    run()
    tm <- take_timer("pca_kernel")
    tot <- phase_totals(tm)
    cnt <- stats::setNames(tm$Count, rownames(tm))
    if (is.null(sums)) {
      sums <- tot
      counts <- cnt
    } else {
      stopifnot(identical(names(tot), names(sums)))
      sums <- sums + tot
      counts <- counts + cnt
    }
  }
  stopifnot(total_tag %in% names(sums))
  total <- sums[[total_tag]]
  outer <- setdiff(names(sums), c(total_tag, inner))
  untimed <- total - sum(sums[outer])
  cat(sprintf(
    "\nPCA %s: %d calls, kernel total %.3f ms a call\n",
    label,
    n_calls,
    total / n_calls / 1000
  ))
  show <- c(sums[outer], untimed = untimed)
  show <- show[order(-show)]
  for (nm in names(show)) {
    cat(sprintf(
      "  %-22s %10.3f ms a call  %6.2f%%\n",
      nm,
      show[[nm]] / n_calls / 1000,
      100 * show[[nm]] / total
    ))
  }
  for (nm in intersect(inner, names(sums))) {
    cat(sprintf(
      "    in svd: %-14s %8.3f ms a call  %6.2f%%  (%d a call)\n",
      nm,
      sums[[nm]] / n_calls / 1000,
      100 * sums[[nm]] / total,
      as.integer(counts[[nm]] / n_calls)
    ))
  }
  invisible(first)
}

profile_knn <- function(label, x) {
  run <- function() knn_imp(x, k = 10)
  set.seed(1)
  warm <- run()
  stopifnot(is.matrix(warm), !anyNA(warm))
  wall <- 0
  loop <- 0
  for (i in seq_len(n_calls)) {
    if (exists("knn_tm_raw", envir = globalenv(), inherits = FALSE)) {
      rm(list = "knn_tm_raw", envir = globalenv())
    }
    w0 <- Sys.time()
    run()
    wall <- wall + as.numeric(difftime(Sys.time(), w0, units = "secs"))
    raw <- get("knn_tm_raw", envir = globalenv())
    stopifnot(nrow(raw) >= 1L, all(raw$tag == "parallel_loop"))
    loop <- loop + sum(raw$duration_ns) / 1e9
  }
  cat(sprintf(
    paste0(
      "\nKNN %s: %d calls, knn_imp() %.3f ms a call, parallel loop %.3f ms ",
      "(%.2f%%), rest %.3f ms (%.2f%%)\n"
    ),
    label,
    n_calls,
    1000 * wall / n_calls,
    1000 * loop / n_calls,
    100 * loop / wall,
    1000 * (wall - loop) / n_calls,
    100 * (wall - loop) / wall
  ))
}

first <- NULL
for (nm in names(data)) {
  for (solver in c("exact", "lobpcg")) {
    tm <- profile_pca(paste(nm, solver), data[[nm]], solver)
    if (is.null(first)) {
      first <- tm
    }
  }
}
cat("\nOne raw summary, as rcpptimer wrote it (warm-up call, chr21 exact):\n")
print(first)

for (nm in names(data)) {
  profile_knn(nm, data[[nm]])
}
