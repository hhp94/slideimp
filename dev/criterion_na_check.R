# Throwaway: print criterion_final on the inputs of the test "criterion_final
# is NA when the objective alone stopped the fit" and run that one test.
# Run from the ROOT OF A SCRATCH CLONE, never the working tree:
#   Rscript dev/criterion_na_check.R --in-scratch-clone
# Used to show the test fails on the build before the NA rule (894826d) and
# passes after it.
args <- commandArgs(trailingOnly = TRUE)
if (!identical(args, "--in-scratch-clone")) {
  stop("run from the root of a scratch clone; never in the working tree")
}
source("dev/dev-utils.R")
load_all1(timer = FALSE)
cat("commit:", system("git rev-parse --short HEAD", intern = TRUE), "\n")

set.seed(21)
base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
x <- base6[rep(seq_len(6), times = 5), ]
rownames(x) <- sprintf("s%02d", seq_len(nrow(x)))
for (ncp in c(5L, 6L, 7L, 12L)) {
  set.seed(1)
  res <- pca_imp(x, ncp = ncp)
  cat(sprintf(
    "ncp = %2d  converged %s  n_iter %d  criterion_final %s\n",
    ncp,
    attr(res, "converged"),
    attr(res, "n_iter"),
    format(attr(res, "criterion_final"), digits = 4)
  ))
}

res <- testthat::test_file(
  "tests/testthat/test-pca_imp.R",
  desc = "criterion_final is NA when the objective alone stopped the fit",
  reporter = "summary",
  stop_on_failure = FALSE
)
df <- as.data.frame(res)
cat(sprintf(
  "\ntests: %d  failed: %d  errors: %d  skipped: %d  passed expectations: %d\n",
  nrow(df),
  sum(df$failed),
  sum(df$error),
  sum(df$skipped),
  sum(df$passed)
))
