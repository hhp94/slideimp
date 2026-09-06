# Does knn_imp() reproduce impute::impute.knn()?
#
# This lives here rather than under tests/ because `impute` is a Bioconductor
# package and is not in Suggests: any `impute::` call under tests/ makes
# R CMD check report an unstated dependency, even behind skip_if_not_installed()
# and the MANUAL_TESTS gate. dev/ is build-ignored, so nothing here reaches the
# check. Declaring `impute` in Suggests would need an Additional_repositories
# entry and is the maintainer's call, not a side effect of adding a test.
#
# The claim used to sit in test-knn_imp.R as a block titled "Exactly replicate
# impute.knn" with the comparison itself commented out, so it had never been
# checked. It has now: the numbers are printed AND asserted.
#
# impute.knn searches over ROWS, knn_imp over COLUMNS, so its input and output
# are transposed. It averages its k neighbors unweighted, which is
# dist_pow = 0 here. maxp = the number of rows of its input disables its
# recursive two-means split, leaving a brute-force search on both sides.
# post_imp fills what K-NN could not, which impute.knn handles differently, so
# it is off.
#
# Usage: Rscript dev/probe_impute_knn.R
# Exits non-zero if any comparison exceeds the tolerance.

setwd("C:/Users/amser/Projects/slideimp")
source("R/dev-utils.R")
suppressMessages(load_all1(timer = FALSE))

if (!requireNamespace("impute", quietly = TRUE)) {
  stop("package 'impute' is not installed; nothing to compare against.")
}

# Tolerance: measured at 2.2e-16 absolute and 4.5e-16 relative, so 1e-12 sits
# four orders above the observed disagreement and far below anything that
# would change an answer. Bitwise equality is reported, never asserted - the
# two are different compilations doing the same arithmetic in different orders.
TOL <- 1e-12
failures <- 0L

report <- function(label, got, want) {
  stopifnot(!anyNA(got), !anyNA(want), identical(dim(got), dim(want)))
  abs_err <- max(abs(got - want))
  big <- abs(want) > 1e-8
  stopifnot(any(big))
  rel_err <- max(abs(got[big] - want[big]) / abs(want[big]))
  ok <- abs_err < TOL && rel_err < TOL
  if (!ok) {
    failures <<- failures + 1L
  }
  cat(sprintf(
    "%-22s max abs %10.3e  max rel %10.3e  over %d of %d cells  bitwise %-5s  %s\n",
    label,
    abs_err,
    rel_err,
    sum(big),
    length(want),
    identical(got, want),
    if (ok) "PASS" else "FAIL"
  ))
}

for (k in c(3L, 10L)) {
  set.seed(1234)
  obj <- sim_mat(100, 100, perc_total_na = 0.05)$input

  r1 <- knn_imp(
    obj,
    k = k,
    method = "euclidean",
    dist_pow = 0,
    post_imp = FALSE,
    na_check = FALSE
  )
  r2 <- t(impute::impute.knn(t(obj), k = k, maxp = ncol(obj))$data)

  na_idx <- which(is.na(obj))
  cat(sprintf("k = %d, %d missing cells\n", k, length(na_idx)))
  report("  imputed cells", r1[,][na_idx], r2[na_idx])
  report("  whole matrix", r1[,], r2)
}

if (failures > 0L) {
  stop(sprintf("%d comparison(s) exceeded the tolerance.", failures))
}
cat("\nall comparisons within", format(TOL), "\n")
