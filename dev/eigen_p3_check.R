# Phase 3 of dev/eigen_conversion.md: devtools::document() or
# devtools::check() on a scratch clone, never on the working tree
# (dev/eigen_p3_gate.sh makes the clone). Run from the ROOT OF THAT CLONE:
#   Rscript dev/eigen_p3_check.R document --in-scratch-clone
#   Rscript dev/eigen_p3_check.R check <check dir> --in-scratch-clone
# Which RcppEigen the build uses is chosen by R_LIBS, set before R starts;
# the script prints the version and the library it was found in, so the log
# says which one it tested.
args <- commandArgs(trailingOnly = TRUE)
if (!"--in-scratch-clone" %in% args) {
  stop("run from the root of a scratch clone; never in the working tree")
}
args <- setdiff(args, "--in-scratch-clone")
mode <- args[[1]]
cat("libPaths:", .libPaths(), sep = "\n  ")
cat(sprintf(
  "\nRcppEigen %s from %s\n",
  format(utils::packageVersion("RcppEigen")),
  dirname(find.package("RcppEigen"))
))

if (identical(mode, "document")) {
  devtools::document()
} else if (identical(mode, "check")) {
  stopifnot(length(args) == 2L)
  check_dir <- args[[2]]
  dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)
  res <- devtools::check(
    document = FALSE,
    check_dir = check_dir,
    error_on = "never"
  )
  cat(sprintf(
    "\ncheck: %d errors, %d warnings, %d notes\n",
    length(res$errors),
    length(res$warnings),
    length(res$notes)
  ))
  for (kind in c("errors", "warnings", "notes")) {
    for (x in res[[kind]]) {
      cat("\n[", kind, "]\n", x, "\n", sep = "")
    }
  }
} else {
  stop("mode is `document` or `check`")
}
