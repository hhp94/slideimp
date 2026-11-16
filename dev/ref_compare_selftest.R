# Self-test for dev/ref_compare.R: prove it still fails where it should.
# Copies a capture to a temporary folder with one planted change per case
# below, runs the compare twice (tolerance 0, then the defaults), and checks
# that exactly the expected cases fail each time. Rerun it after changing the
# compare's walk or its default tolerances.
#   Rscript dev/ref_compare_selftest.R e5b3cf7
# Exit status 1 if either run fails a different set of cases than expected.

label <- commandArgs(trailingOnly = TRUE)[[1]]
src <- file.path("dev", "ref", label)
dst <- file.path(tempdir(), paste0(label, "-perturbed"))
dir.create(dst)
invisible(file.copy(list.files(src, full.names = TRUE), dst))
meta <- readRDS(file.path(dst, "_meta.rds"))
meta$label <- paste0(label, "-perturbed")
saveRDS(meta, file.path(dst, "_meta.rds"))

edit <- function(name, f) {
  p <- file.path(dst, paste0(name, ".rds"))
  saveRDS(f(readRDS(p)), p)
}

## A change of a few ulp in one cell: fails at tolerance 0 only.
edit("pca_x1_exact", function(r) {
  r$value[1, 1] <- r$value[1, 1] * (1 + 4 * .Machine$double.eps)
  r
})
## A change of 0.1 in one cell of the largest result.
edit("knn_chr21", function(r) {
  r$value[nrow(r$value), ncol(r$value)] <- r$value[
    nrow(r$value),
    ncol(r$value)
  ] +
    0.1
  r
})
## An NA where a value was.
edit("pca_x2_post", function(r) {
  r$value[2, 2] <- NA
  r
})
## A dropped attribute.
edit("pca_x1_lobpcg", function(r) {
  attr(r$value, "converged") <- NULL
  r
})
## A changed integer attribute.
edit("knn_x2_post", function(r) {
  attr(r$value, "mean_filled_cols") <- 1:4
  r
})
## A changed double attribute.
edit("pca_x1_em", function(r) {
  attr(r$value, "criterion_final") <- attr(r$value, "criterion_final") * 2
  r
})
## A changed warning message.
edit("slide_xi_knn_mean", function(r) {
  r$warnings[[1]]$message <- "changed"
  r
})
## A changed error message.
edit("pca_xa_allna_error", function(r) {
  r$value$message <- "changed"
  r
})
## A changed dimname.
edit("mean_xa_all", function(r) {
  colnames(r$value)[1] <- "renamed"
  r
})
## A changed NA placement.
edit("na_loc_x1", function(r) {
  r$value[[1]][1] <- r$value[[1]][1] + 1L
  r
})
## A case missing from the new capture.
invisible(file.remove(file.path(dst, "int_find_windows.rds")))

expect_default <- c(
  "int_find_windows",
  "knn_chr21",
  "knn_x2_post",
  "mean_xa_all",
  "na_loc_x1",
  "pca_x1_em",
  "pca_x1_lobpcg",
  "pca_x2_post",
  "pca_xa_allna_error",
  "slide_xi_knn_mean"
)
expect_zero <- sort(c(expect_default, "pca_x1_exact"))

run <- function(extra, expect) {
  out <- suppressWarnings(system2(
    "Rscript",
    c("dev/ref_compare.R", src, dst, extra),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  ## The per-case table: two spaces, the case name, then its status.
  hit <- regmatches(out, regexec("^  ([A-Za-z0-9_]+) +FAIL ", out))
  failed <- sort(vapply(Filter(length, hit), `[[`, "", 2L))
  ok <- identical(failed, sort(expect)) && identical(status, 1L)
  cat(sprintf(
    "%-24s %s: %d failed, %d expected, exit %s\n",
    if (length(extra)) paste(extra, collapse = " ") else "(defaults)",
    if (ok) "ok" else "WRONG",
    length(failed),
    length(expect),
    if (is.null(status)) "0" else status
  ))
  if (!ok) {
    cat("  unexpected:", setdiff(failed, expect), "\n")
    cat("  not caught:", setdiff(expect, failed), "\n")
  }
  ok
}

ok <- c(
  run(c("--tol-abs=0", "--tol-rel=0"), expect_zero),
  run(character(), expect_default)
)
unlink(dst, recursive = TRUE)
quit(status = if (all(ok)) 0L else 1L)
