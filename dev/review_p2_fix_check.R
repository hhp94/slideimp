# Fixes after the review of Phase 2 (dev/review_p2.md, findings 1 and 3):
# run ONE test of tests/testthat/test-pca_imp.R, chosen by its description,
# against an INSTALLED build, to show the test passes on the head and fails
# on a build that has the defect. Dev only, from the package root:
#   Rscript dev/review_p2_fix_check.R <lib> <which>
# <which> is `rank` (the regression test for the dsyevr eigenvalue buffer,
# finding 1) or `overflow` (the stop on a Gram that is not finite, finding
# 3). Both tests call exported functions only, so an installed build can run
# them. A build with the buffer defect ends this process: no summary line
# is printed and the exit status is not 0, which dev/review_p2_fix_check.sh
# reports. It also prints every warning the test's own fits raise, and the
# class and message of the overflow stop.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("usage: Rscript dev/review_p2_fix_check.R <lib> <rank|overflow>")
}
lib <- normalizePath(args[[1]], winslash = "/")
which <- match.arg(args[[2]], c("rank", "overflow"))
desc <- switch(
  which,
  rank = "pca_imp returns on replicate samples with ncp above their rank",
  overflow = "finite input whose squares overflow stops in the eigensolver"
)
library(slideimp, lib.loc = lib)
stopifnot(identical(
  normalizePath(dirname(find.package("slideimp")), winslash = "/"),
  lib
))
cat(
  "slideimp",
  format(utils::packageVersion("slideimp")),
  "from",
  lib,
  "\n"
)

# the warnings of the test's fits, on the same inputs
if (which == "rank") {
  set.seed(21)
  base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  x <- base6[rep(seq_len(6), times = 5), ]
  show_warnings <- function(label, expr) {
    cat("warnings of ", label, ":\n", sep = "")
    flush(stdout())
    n <- 0L
    withCallingHandlers(
      invisible(expr),
      warning = function(w) {
        n <<- n + 1L
        cat("  warning:", gsub("\n", " | ", conditionMessage(w)), "\n")
        invokeRestart("muffleWarning")
      }
    )
    cat("  count:", n, "\n")
    flush(stdout())
  }
  for (k in c(5L, 7L, 12L)) {
    set.seed(1)
    show_warnings(sprintf("pca_imp(x, ncp = %d)", k), pca_imp(x, ncp = k))
  }
  set.seed(1)
  show_warnings(
    "pca_imp(x, ncp = 7, method = \"EM\")",
    pca_imp(x, ncp = 7, method = "EM")
  )
} else {
  set.seed(1234)
  x0 <- sim_mat(30, 10, perc_total_na = 0.1, perc_col_na = 1)$input
  x <- x0 * 1e200
  n_warn <- 0L
  withCallingHandlers(
    invisible(pca_imp(x0, ncp = 2, scale = FALSE)),
    warning = function(w) {
      n_warn <<- n_warn + 1L
      cat("  warning:", gsub("\n", " | ", conditionMessage(w)), "\n")
      invokeRestart("muffleWarning")
    }
  )
  cat("warnings of the fit at the data's own scale, count:", n_warn, "\n")
  for (solver in c("exact", "lobpcg")) {
    e <- tryCatch(
      pca_imp(x, ncp = 2, scale = FALSE, solver = solver),
      error = function(e) e
    )
    cat(
      solver,
      ": [",
      paste(class(e), collapse = "/"),
      "] ",
      if (inherits(e, "error")) conditionMessage(e) else "a value",
      "\n",
      sep = ""
    )
  }
}
flush(stdout())

res <- testthat::test_file(
  "tests/testthat/test-pca_imp.R",
  desc = desc,
  reporter = "summary",
  stop_on_failure = FALSE,
  package = NULL
)
df <- as.data.frame(res)
cat(sprintf(
  "\n%s: tests %d, expectations passed %d, failed %d, errors %d, skipped %d, warnings %d\n",
  which,
  nrow(df),
  sum(df$passed),
  sum(df$failed),
  sum(df$error),
  sum(df$skipped),
  sum(df$warning)
))
stopifnot(nrow(df) == 1L)
