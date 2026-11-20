# Review of Phase 2 (dev/review_p2.md): the suite gate of the plan, rerun on
# a clone. Run from the ROOT OF A CLONE of the commit under test, never from
# the working tree, with R_LIBS naming an empty library, set before R starts
# so that the mirai daemons inherit it (dev/review_p2_suite.sh at aea6898 does
# both):
#   R_LIBS=<empty dir> Rscript dev/review_p2_suite.R --in-scratch-clone
# It installs the clone into that library, checks that a daemon loads that
# copy, then runs the whole suite, manual-gated tests included, against
# load_all1(timer = FALSE): the mipca tests refuse more than one chain under
# the timer build.
args <- commandArgs(trailingOnly = TRUE)
if (!identical(args, "--in-scratch-clone")) {
  stop("run from the root of a scratch clone; never in the working tree")
}
lib <- Sys.getenv("R_LIBS")
stopifnot(nzchar(lib), dir.exists(lib))
lib <- normalizePath(lib, winslash = "/")
cat("libPaths:", .libPaths(), sep = "\n  ")

tarball <- pkgbuild::build(
  dest_path = tempdir(),
  vignettes = FALSE,
  manual = FALSE,
  quiet = TRUE
)
install.packages(
  tarball,
  lib = lib,
  repos = NULL,
  type = "source",
  quiet = TRUE
)

mirai::daemons(1)
seen <- mirai::mirai(list(
  version = format(utils::packageVersion("slideimp")),
  path = normalizePath(find.package("slideimp"), winslash = "/")
))[]
mirai::daemons(0)
cat("\nWhat a daemon sees:\n")
str(seen)
stopifnot(identical(dirname(seen$path), lib))

source("dev/dev-utils.R")
load_all1(timer = FALSE)
res <- check_manual(reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf(
  "\nfiles: %d  tests: %d  failed: %d  errors: %d  skipped: %d  warnings: %d  passed expectations: %d\n",
  length(unique(df$file)),
  nrow(df),
  sum(df$failed),
  sum(df$error),
  sum(df$skipped),
  sum(df$warning),
  sum(df$passed)
))
bad <- df[df$failed > 0 | df$error | df$skipped, c("file", "test")]
if (nrow(bad)) {
  cat("\nFailing or skipped tests:\n")
  print(bad, row.names = FALSE)
}
