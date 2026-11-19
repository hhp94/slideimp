# Fixes after the review of Phase 2 (dev/review_p2.md): the two PCA test
# files, whole, manual-gated tests included, against
# load_all1(timer = FALSE). Run from the ROOT OF A
# CLONE that holds the tests under check, never from the working tree
# (dev/review_p2_fix_check.sh makes the clone and copies the tests in):
#   Rscript dev/review_p2_fix_suite.R --in-scratch-clone
# test-pca_imp-mipca.R is left out: its mirai daemons load the installed
# package, which dev/review_p2_suite.R handles and this script does not.
args <- commandArgs(trailingOnly = TRUE)
if (!identical(args, "--in-scratch-clone")) {
  stop("run from the root of a scratch clone; never in the working tree")
}
source("dev/dev-utils.R")
load_all1(timer = FALSE)
res <- check_manual(
  filter = "^pca_imp(-subset)?$",
  reporter = "summary",
  stop_on_failure = FALSE
)
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
warned <- df[df$warning > 0, c("file", "test")]
if (nrow(warned)) {
  cat("\nTests with a warning:\n")
  print(warned, row.names = FALSE)
}
