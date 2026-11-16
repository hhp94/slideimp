# The mirai tests run their workers in fresh R processes, which load slideimp
# from the library, not from load_all(). After a change to package internals
# the installed copy is stale, and those tests then fail on the version
# mismatch rather than on the change itself.
#
# This installs the working tree into the library named by R_LIBS (set before
# R starts, so the daemons inherit it), checks that a daemon actually sees that
# copy, then runs the test files that start mirai daemons, or with `all`,
# every test file.
#
# Run from the package root, with R_LIBS pointing at an empty directory:
#   R_LIBS=<dir> Rscript dev/test_mirai_with_dev_install.R [all]

run_all <- identical(commandArgs(trailingOnly = TRUE), "all")
lib <- Sys.getenv("R_LIBS")
stopifnot(nzchar(lib), dir.exists(lib))
cat("libPaths:", .libPaths(), sep = "\n  ")
cat("\n")

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
  path = find.package("slideimp"),
  formals = names(formals(slideimp:::new_slideimp_results)),
  has_fallback_part = exists(
    "fallback_part",
    envir = asNamespace("slideimp"),
    inherits = FALSE
  )
))[]
mirai::daemons(0)
cat("\nWhat a daemon sees:\n")
str(seen)
stopifnot(isTRUE(seen$has_fallback_part))

source("dev/dev-utils.R")
load_all1()
res <- check_manual(
  filter = if (run_all) NULL else "group_imp|tune_imp",
  reporter = "summary",
  stop_on_failure = FALSE
)
df <- as.data.frame(res)
cat(sprintf(
  "\nfiles: %d  tests: %d  failed: %d  errors: %d  skipped: %d  passed expectations: %d\n",
  length(unique(df$file)),
  nrow(df),
  sum(df$failed),
  sum(df$error),
  sum(df$skipped),
  sum(df$passed)
))
bad <- df[df$failed > 0 | df$error, c("file", "test")]
if (nrow(bad)) {
  cat("\nFailing tests:\n")
  print(bad, row.names = FALSE)
}
