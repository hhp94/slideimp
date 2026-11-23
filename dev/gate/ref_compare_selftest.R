# Self-test for dev/gate/ref_compare.R: prove it still fails where it should.
# Copies a capture twice to a temporary folder, a ref side and a new side,
# plants one change per case below (most on the new side; the pair that
# pins the relative figure's denominator on both), runs the compare twice
# (tolerance 0, then the defaults), and checks that exactly the expected
# cases fail each time. It also checks that a 0 turned into -0 passes but is
# reported as not bitwise, and that a change to an excluded attribute
# (criterion_final) alone passes, is reported as excluded, and has its
# difference printed under --verbose, while the same change next to a value
# change, or a change to a double attribute that is not excluded, fails.
# Rerun it after changing the compare's walk, its exclusions or its default
# tolerances. The label names a capture taken by the current
# dev/gate/ref_capture.R, since the planted changes name its cases.
#   Rscript dev/gate/ref_compare_selftest.R e5b3cf7
# Exit status 1 if any run fails a different set of cases than expected, the
# -0 is reported bitwise, or an exclusion is not reported as described.

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
## An excluded attribute changed alone: passes under both runs, and is
## reported as excluded. Once doubled, once turned into NA.
edit("pca_x1_em", function(r) {
  stopifnot(isTRUE(attr(r$value, "criterion_final") > 0))
  attr(r$value, "criterion_final") <- attr(r$value, "criterion_final") * 2
  r
})
edit("pca_x1_auto", function(r) {
  stopifnot(is.finite(attr(r$value, "criterion_final")))
  attr(r$value, "criterion_final") <- NA_real_
  r
})
## An excluded attribute changed with a value cell: the exclusion does not
## hide the rest of the case.
edit("pca_x1_clamp", function(r) {
  attr(r$value, "criterion_final") <- attr(r$value, "criterion_final") * 2
  r$value[1, 1] <- r$value[1, 1] + 1e-6
  r
})
## A changed double attribute that is not excluded, planted on both sides
## under a name that extends the excluded one, so the exclusion is shown to
## match the whole name only.
edit(
  "pca_x1_nbinit",
  function(r) {
    attr(r$value, "criterion_final_x") <- 1e-7
    r
  },
  dir = ref
)
edit("pca_x1_nbinit", function(r) {
  attr(r$value, "criterion_final_x") <- 2e-7
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
edit("knn_x1_euc", function(r) {
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
  "knn_x1_euc",
  "knn_x2_post",
  "na_loc_x1",
  "pca_x1_clamp",
  "pca_x1_lobpcg",
  "pca_x1_nbinit",
  "pca_x1_noscale",
  "pca_x1_roww",
  "pca_x2_post",
  "pca_xa_allna_error",
  "slide_xi_knn_mean"
)
expect_zero <- sort(c(expect_default, "pca_x1_exact"))
## Cases whose only change is to an excluded attribute: they pass, and the
## per-case table reports the exclusion.
expect_excluded_ok <- c("pca_x1_auto", "pca_x1_em")

compare <- function(extra) {
  suppressWarnings(system2(
    "Rscript",
    c("dev/gate/ref_compare.R", ref, dst, extra),
    stdout = TRUE,
    stderr = TRUE
  ))
}

run <- function(extra, expect) {
  out <- compare(extra)
  status <- attr(out, "status")
  ## The per-case table: two spaces, the case name, then its status.
  hit <- regmatches(out, regexec("^  ([A-Za-z0-9_]+) +FAIL ", out))
  failed <- sort(vapply(Filter(length, hit), `[[`, "", 2L))
  neg_zero <- grep("^  col_vars_x2 +ok .*NOT bitwise", out, value = TRUE)
  excl_ok <- vapply(
    expect_excluded_ok,
    function(nm) {
      length(grep(
        sprintf("^  %s +ok .*\\(1 excluded\\)$", nm),
        out
      )) ==
        1L
    },
    TRUE
  )
  ok <- identical(failed, sort(expect)) &&
    identical(status, 1L) &&
    length(neg_zero) == 1L &&
    all(excl_ok)
  cat(sprintf(
    "%-24s %s: %d failed, %d expected, exit %s; -0 %s; excluded-only %d/%d ok\n",
    if (length(extra)) paste(extra, collapse = " ") else "(defaults)",
    if (ok) "ok" else "WRONG",
    length(failed),
    length(expect),
    if (is.null(status)) "0" else status,
    if (length(neg_zero)) "passed, not bitwise" else "NOT SEEN",
    sum(excl_ok),
    length(excl_ok)
  ))
  if (!ok) {
    cat("  unexpected:", setdiff(failed, expect), "\n")
    cat("  not caught:", setdiff(expect, failed), "\n")
    cat("  excluded-only not reported ok:", names(excl_ok)[!excl_ok], "\n")
  }
  ok
}

## --verbose on the doubled criterion_final alone, at tolerance 0: the case
## passes, and its excluded row carries the abs difference, which is the
## planted one.
run_verbose <- function() {
  a <- readRDS(file.path(ref, "pca_x1_em.rds"))
  want <- abs(attr(a$value, "criterion_final"))
  out <- compare(c(
    "--tol-abs=0",
    "--tol-rel=0",
    "--only=^pca_x1_em$",
    "--verbose"
  ))
  status <- attr(out, "status")
  row <- grep(
    "^  excl +value@criterion_final +excluded +abs ",
    out,
    value = TRUE
  )
  got <- if (length(row) == 1L) {
    as.numeric(sub("^.* abs ([^ ]+) .*$", "\\1", row))
  } else {
    NA_real_
  }
  ok <- is.null(status) &&
    length(row) == 1L &&
    isTRUE(abs(got - want) <= 0.01 * want)
  cat(sprintf(
    "%-24s %s: exit %s; excluded row %s, abs %s (planted %s)\n",
    "--verbose, one case",
    if (ok) "ok" else "WRONG",
    if (is.null(status)) "0" else status,
    if (length(row) == 1L) "printed" else "NOT SEEN",
    formatC(got, format = "e", digits = 2),
    formatC(want, format = "e", digits = 2)
  ))
  ok
}

ok <- c(
  run(c("--tol-abs=0", "--tol-rel=0"), expect_zero),
  run(character(), expect_default),
  run_verbose()
)
unlink(c(ref, dst), recursive = TRUE)
quit(status = if (all(ok)) 0L else 1L)
