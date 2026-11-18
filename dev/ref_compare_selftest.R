# Self-test for dev/ref_compare.R: prove it still fails where it should.
# Copies a capture twice to a temporary folder, a ref side and a new side,
# plants one change per case below (most on the new side; the pair that
# pins the relative figure's denominator on both), runs the compare twice
# (tolerance 0, then the defaults), and checks that exactly the expected
# cases fail each time. It also checks that a 0 turned into -0 passes but is
# reported as not bitwise. Rerun it after changing the compare's walk or its
# default tolerances.
#   Rscript dev/ref_compare_selftest.R e5b3cf7
# Exit status 1 if either run fails a different set of cases than expected,
# or the -0 is reported bitwise.

label <- commandArgs(trailingOnly = TRUE)[[1]]
src <- file.path("dev", "ref", label)
copy_side <- function(side) {
  d <- file.path(tempdir(), paste0(label, "-", side))
  dir.create(d)
  invisible(file.copy(list.files(src, full.names = TRUE), d))
  meta <- readRDS(file.path(d, "_meta.rds"))
  meta$label <- paste0(label, "-", side)
  saveRDS(meta, file.path(d, "_meta.rds"))
  d
}
ref <- copy_side("ref")
dst <- copy_side("perturbed")

edit <- function(name, f, dir = dst) {
  p <- file.path(dir, paste0(name, ".rds"))
  saveRDS(f(readRDS(p)), p)
}
set_cell <- function(v) {
  function(r) {
    r$value[1, 1] <- v
    r
  }
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
## The denominator pair: a value under the abs tolerance that one side holds
## as 0, once each way. Both must fail, whichever side the 0 is on.
edit("pca_x1_noscale", set_cell(5e-13), dir = ref)
edit("pca_x1_noscale", set_cell(0))
edit("pca_x1_roww", set_cell(0), dir = ref)
edit("pca_x1_roww", set_cell(5e-13))
## A 0 turned into -0: equal values, different bytes. Passes, reported as
## not bitwise.
edit("col_vars_x2", function(r) {
  z <- which(r$value == 0)
  stopifnot(length(z) > 0L)
  r$value[z[[1]]] <- -0
  r
})

expect_default <- c(
  "int_find_windows",
  "knn_chr21",
  "knn_x2_post",
  "mean_xa_all",
  "na_loc_x1",
  "pca_x1_em",
  "pca_x1_lobpcg",
  "pca_x1_noscale",
  "pca_x1_roww",
  "pca_x2_post",
  "pca_xa_allna_error",
  "slide_xi_knn_mean"
)
expect_zero <- sort(c(expect_default, "pca_x1_exact"))

run <- function(extra, expect) {
  out <- suppressWarnings(system2(
    "Rscript",
    c("dev/ref_compare.R", ref, dst, extra),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  ## The per-case table: two spaces, the case name, then its status.
  hit <- regmatches(out, regexec("^  ([A-Za-z0-9_]+) +FAIL ", out))
  failed <- sort(vapply(Filter(length, hit), `[[`, "", 2L))
  neg_zero <- grep("^  col_vars_x2 +ok .*NOT bitwise", out, value = TRUE)
  ok <- identical(failed, sort(expect)) &&
    identical(status, 1L) &&
    length(neg_zero) == 1L
  cat(sprintf(
    "%-24s %s: %d failed, %d expected, exit %s; -0 %s\n",
    if (length(extra)) paste(extra, collapse = " ") else "(defaults)",
    if (ok) "ok" else "WRONG",
    length(failed),
    length(expect),
    if (is.null(status)) "0" else status,
    if (length(neg_zero)) "passed, not bitwise" else "NOT SEEN"
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
unlink(c(ref, dst), recursive = TRUE)
quit(status = if (all(ok)) 0L else 1L)
