# Install one build of slideimp into its own library, for the timing driver
# dev/eigen_time_builds.R (dev/eigen_conversion.md at 4f2b669, Phase 0,
# "Timing").
#   Rscript dev/eigen_time_install.R <source dir> <lib>
# <source dir> is a checkout of the build to time, outside this tree, since
# other sessions may be building here, and used for nothing else while this
# runs: R CMD INSTALL builds in place, and --preclean removes the objects a
# load_all1() build in the same folder relies on. <lib> is created if
# missing and should hold nothing but this build.
#
# Every build is installed the same way, plain R CMD INSTALL, so the
# reference and a converted build differ only in their source. The flags are
# R's plus whatever the user Makevars (~/.R/Makevars.win on Windows) adds,
# which applies to both alike; install both from the same machine and
# account. Never time a load_all1() build against an installed one.
#
# The build compiles against whichever RcppEigen and RcppArmadillo come first
# on the library path. Set R_LIBS before R starts to choose them; this script
# does not. It writes <lib>/_build.txt with the source commit and the
# versions the build compiled against, which the worker prints with every
# timing.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "usage: Rscript dev/eigen_time_install.R <source dir> <lib>",
    call. = FALSE
  )
}
src <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
lib <- args[[2]]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
lib <- normalizePath(lib, winslash = "/", mustWork = TRUE)

git <- function(...) {
  system2("git", c("-C", shQuote(src), ...), stdout = TRUE)
}
src_paths <- c("R", "src", "DESCRIPTION", "NAMESPACE")
info <- c(
  source = src,
  head = git("rev-parse", "--short", "HEAD"),
  source_commit = git("log", "-1", "--format=%h", "--", src_paths),
  source_dirty = length(git("status", "--porcelain", "--", src_paths)) > 0L,
  rcppeigen = as.character(packageVersion("RcppEigen")),
  rcppeigen_lib = dirname(find.package("RcppEigen")),
  rcpparmadillo = as.character(packageVersion("RcppArmadillo")),
  r_version = R.version.string,
  time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
cat(sprintf("%-14s %s\n", names(info), info), sep = "")

r_exe <- file.path(R.home("bin"), "R")
status <- system2(
  r_exe,
  c(
    "CMD",
    "INSTALL",
    "--preclean",
    paste0("--library=", shQuote(lib)),
    shQuote(src)
  )
)
if (!identical(status, 0L)) {
  stop("R CMD INSTALL failed with status ", status, call. = FALSE)
}
writeLines(sprintf("%s=%s", names(info), info), file.path(lib, "_build.txt"))
cat(
  "installed slideimp",
  as.character(packageVersion("slideimp", lib.loc = lib)),
  "into",
  lib,
  "\n"
)
