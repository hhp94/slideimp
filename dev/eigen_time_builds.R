# Timing driver (dev/eigen_conversion.md, Phase 0, "Timing"): two installed
# builds of slideimp, each in its own library, timed in alternating fresh
# Rscript processes; reports medians over the interleaved repetitions.
#   Rscript dev/eigen_time_builds.R --lib-a=<lib> --lib-b=<lib> [options]
# Install each build first with dev/eigen_time_install.R; a is the reference,
# b the converted build. Proving the driver: pass the reference twice, as two
# separate installs, and the ratios should sit near 1.
#
# Options:
#   --reps=<n>          repetitions per build and workload (default 3)
#   --workloads=<a,b>   from dev/eigen_time_worker.R (default: all six)
#   --chr21=<path>      default dev/GSE286313.EPICv1_chr21.qs2
#   --sim-n=<n> --sim-p=<p>   the larger sim_mat() (default 200 x 5000)
#   --threads=<n>       sets OPENBLAS_NUM_THREADS and OMP_NUM_THREADS for the
#                       workers; unset by default, so each runs with the
#                       machine's defaults, OpenMP included
#   --out=<file.csv>    also write every run there
#
# Each worker is a fresh process, so no build sees the other's DLL, and the
# library is chosen with R_LIBS set before the process starts. The worker
# stops if slideimp loaded from anywhere else. Within a repetition, each
# workload runs a then b on odd repetitions and b then a on even ones, so a
# slow drift in the machine falls on both builds alike.
#
# Only the imputer calls are timed: each worker runs its call 5 times in one
# replicate() and reports the total (dev/eigen_time_worker.R, n_calls), so
# every time below is for 5 calls. Elapsed is wall time; cpu/wall is
# (user + system) / elapsed of those calls, which shows how many threads
# they kept busy (dev/CPP.md: a pin is verified by CPU time over wall time).

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
known <- c(
  "lib-a",
  "lib-b",
  "reps",
  "workloads",
  "chr21",
  "sim-n",
  "sim-p",
  "threads",
  "out"
)
bad <- setdiff(sub("=.*$", "", sub("^--", "", args)), known)
if (length(bad) || is.null(opt("lib-a")) || is.null(opt("lib-b"))) {
  stop(
    "usage: Rscript dev/eigen_time_builds.R --lib-a=<lib> --lib-b=<lib> ",
    "[--reps=n] [--workloads=a,b] [--chr21=path] [--sim-n=n] [--sim-p=p] ",
    "[--threads=n] [--out=file.csv]",
    call. = FALSE
  )
}

libs <- c(
  a = normalizePath(opt("lib-a"), winslash = "/", mustWork = TRUE),
  b = normalizePath(opt("lib-b"), winslash = "/", mustWork = TRUE)
)
reps <- as.integer(opt("reps", "3"))
workloads <- strsplit(
  opt(
    "workloads",
    paste(
      c(
        "chr21_pca_exact",
        "chr21_pca_lobpcg",
        "chr21_knn",
        "sim_pca_exact",
        "sim_pca_lobpcg",
        "sim_knn"
      ),
      collapse = ","
    )
  ),
  ","
)[[1]]
chr21 <- normalizePath(
  opt("chr21", "dev/GSE286313.EPICv1_chr21.qs2"),
  winslash = "/",
  mustWork = TRUE
)
sim_n <- opt("sim-n", "200")
sim_p <- opt("sim-p", "5000")
threads <- opt("threads")
out <- opt("out")
worker <- normalizePath("dev/eigen_time_worker.R", mustWork = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")

if (!is.null(threads)) {
  Sys.setenv(OPENBLAS_NUM_THREADS = threads, OMP_NUM_THREADS = threads)
}
cat(sprintf("a: %s\nb: %s\n", libs[["a"]], libs[["b"]]))
for (k in names(libs)) {
  f <- file.path(libs[[k]], "_build.txt")
  if (file.exists(f)) {
    cat(sprintf(
      "  %s built from: %s\n",
      k,
      paste(readLines(f), collapse = " ")
    ))
  }
}
cat(sprintf(
  "reps %d, sim %s x %s, OPENBLAS_NUM_THREADS '%s', OMP_NUM_THREADS '%s'\n",
  reps,
  sim_n,
  sim_p,
  Sys.getenv("OPENBLAS_NUM_THREADS"),
  Sys.getenv("OMP_NUM_THREADS")
))
cat(sprintf("workloads: %s\n\n", paste(workloads, collapse = ", ")))

run_one <- function(build, workload) {
  Sys.setenv(R_LIBS = libs[[build]])
  res <- suppressWarnings(system2(
    rscript,
    c(
      shQuote(worker),
      workload,
      shQuote(libs[[build]]),
      shQuote(chr21),
      sim_n,
      sim_p
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  line <- grep("^RESULT ", res, value = TRUE)
  if (!is.null(attr(res, "status")) || length(line) != 1L) {
    cat(res, sep = "\n")
    stop("worker failed: build ", build, ", ", workload, call. = FALSE)
  }
  f <- strsplit(line, " ", fixed = TRUE)[[1]]
  data.frame(
    build = build,
    workload = workload,
    elapsed = as.numeric(f[[3]]),
    cpu = as.numeric(f[[4]]) + as.numeric(f[[5]]),
    slideimp = f[[6]],
    rcppeigen = f[[7]]
  )
}

t_start <- proc.time()[["elapsed"]]
runs <- list()
for (r in seq_len(reps)) {
  order <- if (r %% 2L == 1L) c("a", "b") else c("b", "a")
  for (w in workloads) {
    for (b in order) {
      x <- run_one(b, w)
      x$rep <- r
      runs[[length(runs) + 1L]] <- x
      cat(sprintf(
        "  rep %d  %s  %-18s %8.3fs  cpu/wall %5.2f  slideimp %s, RcppEigen %s\n",
        r,
        b,
        w,
        x$elapsed,
        x$cpu / x$elapsed,
        x$slideimp,
        x$rcppeigen
      ))
    }
  }
}
runs <- do.call(rbind, runs)
if (!is.null(out)) {
  utils::write.csv(runs, out, row.names = FALSE)
}

cat(
  "\nmedian elapsed (s, 5 calls per run) over",
  reps,
  "rep(s); ratio is b / a:\n"
)
cat(sprintf(
  "  %-18s %9s %9s %7s   %s\n",
  "workload",
  "a",
  "b",
  "ratio",
  "cpu/wall a, b"
))
for (w in workloads) {
  ea <- runs$elapsed[runs$workload == w & runs$build == "a"]
  eb <- runs$elapsed[runs$workload == w & runs$build == "b"]
  ca <- runs$cpu[runs$workload == w & runs$build == "a"]
  cb <- runs$cpu[runs$workload == w & runs$build == "b"]
  cat(sprintf(
    "  %-18s %9.3f %9.3f %7.3f   %.2f, %.2f   (a %.3f-%.3f, b %.3f-%.3f)\n",
    w,
    median(ea),
    median(eb),
    median(eb) / median(ea),
    median(ca / ea),
    median(cb / eb),
    min(ea),
    max(ea),
    min(eb),
    max(eb)
  ))
}
cat(sprintf(
  "\n%d runs, %.1fs wall in total\n",
  nrow(runs),
  proc.time()[["elapsed"]] - t_start
))
