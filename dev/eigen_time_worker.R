# One timed workload, in a fresh process, for dev/eigen_time_builds.R. Not meant
# to be run by hand except to try a workload:
#   Rscript dev/eigen_time_worker.R <workload> <lib> <chr21.qs2> <sim_n> <sim_p>
# <lib> is the library the slideimp under test must load from. The driver
# puts it first through R_LIBS before this process starts; this script only
# checks that it did, and stops if slideimp came from anywhere else.
#
# Workloads: pca_exact, pca_lobpcg and knn, each on chr21 (dev/GSE286313...
# without its all-NA probes) or on sim_mat(sim_n, sim_p) under set.seed(1):
#   chr21_pca_exact chr21_pca_lobpcg chr21_knn sim_pca_exact sim_pca_lobpcg
#   sim_knn
# and pca_lobpcg1, the LOBPCG call with one component, and knn4, the K-NN
# call with cores = 4, on either.
# Only the imputer calls are timed, not loading the data. One untimed call
# runs first: the first call in a process pays for lazy loading and first
# BLAS use, 1.4 to 2.5 times a later one on the chr21 PCA workloads, and
# that cost is common to both builds, so it would pull every ratio towards 1.
# Then the call runs n_calls times in one replicate(), timed as a whole: a
# single PCA call on chr21 takes under 0.1s. The times reported are the
# total over the n_calls calls. Elapsed is read from Sys.time(), which steps
# in about 2 us here, not proc.time(), which steps in 10 ms on Windows, a
# few percent of a chr21 PCA total; user and system time come from
# proc.time() and keep its 10 ms steps, so cpu/wall is coarse on short runs.
#
# The last line printed is the result, for the driver to parse:
#   RESULT <workload> <elapsed> <user> <system> <slideimp> <RcppEigen>

n_calls <- 5L

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    "usage: Rscript dev/eigen_time_worker.R <workload> <lib> <chr21.qs2> ",
    "<sim_n> <sim_p>",
    call. = FALSE
  )
}
workload <- args[[1]]
lib <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
chr21_path <- args[[3]]
sim_n <- as.integer(args[[4]])
sim_p <- as.integer(args[[5]])

suppressPackageStartupMessages(library(slideimp))
from <- normalizePath(dirname(find.package("slideimp")), winslash = "/")
if (!identical(tolower(from), tolower(lib))) {
  stop("slideimp loaded from ", from, ", not ", lib, call. = FALSE)
}
## What the build compiled against, written by dev/eigen_time_install.R.
## packageVersion("RcppEigen") below is what this process has on its path,
## which is the same only when R_LIBS is the same as at install time.
build <- file.path(lib, "_build.txt")
if (file.exists(build)) {
  cat("BUILD", paste(readLines(build), collapse = " "), "\n")
}

parts <- regmatches(workload, regexec("^(chr21|sim)_(.+)$", workload))[[1]]
if (length(parts) != 3L) {
  stop("unknown workload: ", workload, call. = FALSE)
}
data_name <- parts[[2]]
call_name <- parts[[3]]

if (data_name == "chr21") {
  x <- t(qs2::qs_read(chr21_path))
  x <- x[, colSums(!is.na(x)) > 0L]
} else {
  set.seed(1)
  x <- sim_mat(n = sim_n, p = sim_p)$input
}

run <- switch(
  call_name,
  pca_exact = function() pca_imp(x, ncp = 5, solver = "exact"),
  pca_lobpcg = function() pca_imp(x, ncp = 5, solver = "lobpcg"),
  ## One component: every block of the LOBPCG iteration has one column, so
  ## every product is matrix times vector (item 9 of Phase 4).
  pca_lobpcg1 = function() pca_imp(x, ncp = 1, solver = "lobpcg"),
  knn = function() knn_imp(x, k = 10),
  ## The same call over four worker threads.
  knn4 = function() knn_imp(x, k = 10, cores = 4),
  stop("unknown workload: ", workload, call. = FALSE)
)

set.seed(1)
warm <- run()
stopifnot(is.matrix(warm), !anyNA(warm))
t0 <- proc.time()
w0 <- Sys.time()
res <- replicate(n_calls, run(), simplify = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), w0, units = "secs"))
dt <- proc.time() - t0
stopifnot(
  length(res) == n_calls,
  all(vapply(res, function(r) is.matrix(r) && !anyNA(r), logical(1)))
)

cat(sprintf(
  "RESULT %s %.6f %.4f %.4f %s %s\n",
  workload,
  elapsed,
  dt[["user.self"]],
  dt[["sys.self"]],
  as.character(packageVersion("slideimp")),
  as.character(packageVersion("RcppEigen"))
))
