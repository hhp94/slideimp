# Step 6 of Phase 2 (dev/eigen_conversion.md) made the `n_iter` attribute of
# pca_imp()'s result an integer, where it was a double holding a whole
# number. dev/ref_compare.R fails every case that carries the attribute, on
# its type alone, and then says no more about that attribute. This reads two
# captures and shows the type is all that changed: for every case that
# returned a value in both, it reports the type of `n_iter` on each side and
# whether the counts are equal, then whether the two results are identical()
# once the count is read as a double on both sides.
#   Rscript dev/probe_p2s6_niter.R <capture before> <capture after>
# Arguments are labels under dev/ref/, or paths to capture folders.
# identical() is reported, not asserted: the two captures are bytes from two
# builds. The script stops only if it found no case to compare.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "usage: Rscript dev/probe_p2s6_niter.R <before> <after>",
    call. = FALSE
  )
}
resolve <- function(x) if (dir.exists(x)) x else file.path("dev", "ref", x)
dir_a <- resolve(args[[1]])
dir_b <- resolve(args[[2]])
cases <- function(d) {
  setdiff(sub("\\.rds$", "", list.files(d, pattern = "\\.rds$")), "_meta")
}
both <- intersect(cases(dir_a), cases(dir_b))

n_cases <- 0L
n_with <- 0L
n_type <- 0L
n_equal <- 0L
n_identical <- 0L
for (nm in both) {
  a <- readRDS(file.path(dir_a, paste0(nm, ".rds")))
  b <- readRDS(file.path(dir_b, paste0(nm, ".rds")))
  ## a capture folder also holds files that are not cases (its timings)
  if (!is.list(a) || !is.list(b) || is.null(a$value) || is.null(b$value)) {
    next
  }
  n_cases <- n_cases + 1L
  a <- a$value
  b <- b$value
  ia <- attr(a, "n_iter")
  ib <- attr(b, "n_iter")
  if (is.null(ia) || is.null(ib)) {
    next
  }
  n_with <- n_with + 1L
  type_ok <- is.double(ia) && is.integer(ib)
  equal <- length(ia) == length(ib) && all(ia == ib)
  attr(a, "n_iter") <- as.numeric(ia)
  attr(b, "n_iter") <- as.numeric(ib)
  same <- identical(a, b)
  n_type <- n_type + type_ok
  n_equal <- n_equal + equal
  n_identical <- n_identical + same
  if (!type_ok || !equal || !same) {
    cat(sprintf(
      "  %-30s %s -> %s, counts %s, rest %s\n",
      nm,
      typeof(ia),
      typeof(ib),
      if (equal) "equal" else "DIFFER",
      if (same) "identical" else "NOT identical"
    ))
  }
}
cat(sprintf(
  paste0(
    "%d cases with a value in both captures, %d carry n_iter: double ",
    "before and integer after in %d, equal counts in %d, identical with ",
    "the count read as a double in %d\n"
  ),
  n_cases,
  n_with,
  n_type,
  n_equal,
  n_identical
))
if (n_with == 0L) {
  stop("no case carries n_iter in both captures", call. = FALSE)
}
