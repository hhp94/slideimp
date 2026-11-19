# Audit of dev/review_eigen_wrappers.md, findings 1 and 2: what the chain's
# own dsyevr call does when pca_imp() reaches it.
# dev/audit_eigen_wrappers_chain.sh runs this from the root of a scratch
# clone patched by dev/audit_eigen_wrappers_patch.R:
#   Rscript dev/audit_eigen_wrappers_chain.R <log file> [<real data file>]
# Prints, per pca_imp() call: how many dsyevr calls were made, and the most
# cells any of them wrote in W past the topk cells the chain allocates.
args <- commandArgs(trailingOnly = TRUE)
log <- args[[1]]
real <- if (length(args) > 1L) args[[2]] else NA_character_
Sys.setenv(SLIDEIMP_AUDIT_LOG = log)
source("dev/dev-utils.R")
suppressMessages(load_all1(timer = FALSE))
cat("R:", R.version.string, " LAPACK", La_version(), "\n")

run <- function(tag, obj, ...) {
  unlink(log)
  err <- tryCatch(
    {
      res <- pca_imp(obj, ...)
      sprintf("ok, %d NA left", sum(is.na(res)))
    },
    error = function(e) paste("ERROR:", conditionMessage(e))
  )
  calls <- if (file.exists(log)) {
    utils::read.table(
      log,
      col.names = c("n", "topk", "info", "m_out", "past")
    )
  } else {
    data.frame(n = 0L, topk = 0L, info = 0L, m_out = 0L, past = 0L)[0, ]
  }
  cat(sprintf(
    "%-52s %3d x %-5d %s\n    dsyevr calls %3d, Gram order %s, pairs asked %s, calls that wrote past topk %d, most cells past %d, info != 0 in %d, short count in %d\n",
    tag,
    nrow(obj),
    ncol(obj),
    err,
    nrow(calls),
    paste(unique(calls$n), collapse = ","),
    paste(unique(calls$topk), collapse = ","),
    sum(calls$past > 0),
    if (nrow(calls)) max(calls$past) else 0L,
    sum(calls$info != 0),
    sum(calls$m_out != calls$topk)
  ))
  invisible(calls)
}

set.seed(20261118)
rmat <- function(m, n) matrix(rnorm(m * n), m, n)
holes <- function(x, frac) {
  x[sample(length(x), round(frac * length(x)))] <- NA
  x
}

cat("\n-- controls: full-rank simulated data\n")
run("random 30 x 100, ncp 2", holes(rmat(30, 100), 0.05), ncp = 2)
run("random 30 x 100, ncp 7", holes(rmat(30, 100), 0.05), ncp = 7)
run("random 200 x 20, ncp 5", holes(rmat(200, 20), 0.05), ncp = 5)

cat("\n-- duplicated samples: 6 distinct rows, each 5 times, holes included\n")
base <- holes(rmat(6, 100), 0.05)
dup <- base[rep(1:6, each = 5), ]
for (ncp in c(2, 4, 5, 7, 12, 20)) {
  for (solver in c("exact", "auto")) {
    run(
      sprintf("6 distinct of 30, ncp %d, solver %s", ncp, solver),
      dup,
      ncp = ncp,
      solver = solver
    )
  }
}
run("6 distinct of 30, ncp 7, method EM", dup, ncp = 7, method = "EM")

cat(
  "\n-- duplicated columns: 10 columns of which 3 identical, holes included\n"
)
tall <- holes(rmat(40, 8), 0.05)
tall <- tall[, c(1:8, 8, 8)]
for (ncp in c(2, 7, 8, 9)) {
  run(
    sprintf("3 identical of 10 columns, ncp %d", ncp),
    tall,
    ncp = ncp,
    solver = "exact"
  )
}

cat(
  "\n-- the full range: ncp = eligible columns - 1, regularized (k = Gram order)\n"
)
run("random 40 x 6, ncp 5", holes(rmat(40, 6), 0.05), ncp = 5, solver = "exact")

if (!is.na(real) && requireNamespace("qs2", quietly = TRUE)) {
  cat("\n-- real data: GSE286313 EPIC v1 chr21, samples x probes\n")
  x <- t(qs2::qs_read(real))
  x <- x[, colSums(!is.na(x)) > 0, drop = FALSE]
  for (ncp in c(2, 10, 40)) {
    run(
      sprintf("chr21, ncp %d, solver exact", ncp),
      x,
      ncp = ncp,
      solver = "exact"
    )
  }
}
