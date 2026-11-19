# Review of Phase 2 (dev/review_p2.md): small readings of what
# dev/review_p2_cases.R wrote, which the one-line compare does not print.
#   Rscript dev/review_p2_show.R <old dir> <new dir>
# 1. the type of every `n_iter` and `iter` a mipca fit holds, on each side;
# 2. over the tinysweep group, how many fits changed their iteration count,
#    their `converged` flag or their solver counters.
#   Rscript dev/review_p2_show.R --diag <capture dir of dev/eigen_p2s6_diag.R>
# prints the name of each element of one diagnostics case, by position,
# since dev/ref_compare.R names a differing element by position only.
args <- commandArgs(trailingOnly = TRUE)
if (identical(args[[1]], "--diag")) {
  v <- readRDS(file.path(args[[2]], "xb_ncp40_lobpcg.rds"))$value
  print(
    data.frame(
      position = seq_along(v),
      name = names(v),
      type = vapply(v, typeof, ""),
      length = lengths(v)
    ),
    row.names = FALSE
  )
  quit(status = 0L)
}
old_dir <- args[[1]]
new_dir <- args[[2]]
rd <- function(d, nm) readRDS(file.path(d, paste0(nm, ".rds")))$value

cat("== mipca fit: type of each n_iter and iter, old then new\n")
o <- rd(old_dir, "general__tall_mipca_types")
n <- rd(new_dir, "general__tall_mipca_types")
print(data.frame(path = names(o), old = unname(o), new = unname(n[names(o)])))
cat("paths only in new:", setdiff(names(n), names(o)), "\n")

cat("== mipca fit: all.equal(old, new) at tolerance 0\n")
print(all.equal(
  rd(old_dir, "general__tall_mipca"),
  rd(new_dir, "general__tall_mipca"),
  tolerance = 0
))

cat("== tinysweep: attributes that changed, by shape\n")
files <- list.files(old_dir, "^tinysweep__.*\\.rds$")
shape <- sub("^tinysweep__(.*)_s[0-9]+\\.rds$", "\\1", files)
chg <- function(f, what) {
  a <- attr(readRDS(file.path(old_dir, f))$value, what)
  b <- attr(readRDS(file.path(new_dir, f))$value, what)
  !isTRUE(all.equal(as.double(a), as.double(b), tolerance = 0))
}
res <- data.frame(
  shape = shape,
  n_iter = vapply(files, chg, NA, what = "n_iter"),
  converged = vapply(files, chg, NA, what = "converged"),
  n_exact = vapply(files, chg, NA, what = "n_exact"),
  n_lobpcg_ok = vapply(files, chg, NA, what = "n_lobpcg_ok"),
  n_lobpcg_bad = vapply(files, chg, NA, what = "n_lobpcg_bad")
)
print(aggregate(. ~ shape, data = res, FUN = sum))
ran <- vapply(
  files,
  function(f) {
    v <- readRDS(file.path(new_dir, f))$value
    attr(v, "n_lobpcg_ok") + attr(v, "n_lobpcg_bad")
  },
  0
)
cat(
  "fits of the new build in which LOBPCG was tried at least once, by shape:\n"
)
print(tapply(ran > 0, shape, sum))
