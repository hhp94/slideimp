# Probe: does a change to how src/mipca.cpp prints its progress line change
# the bytes printed?
#
# Run once on the old code with tag `before`, once on the new with tag
# `after`, then compare the two files byte for byte. Identical means the
# change printed nothing different on this input.
#
# Same input and calls as hhp94/test's dev/prototypes/probe-progress-bytes.R:
# a 40 x 12 uniform matrix with 30 holes, two chains of 60 warmup and 60
# sampling iterations, at refresh 1, 7 and 50. The "finished in" lines carry
# a timing and are dropped. Two chains need the timer off (the sampler
# refuses more than one chain under -DLOC_TIMER), hence timer = FALSE.
#   Rscript dev/probe_mipca_progress_bytes.R before
# Writes <out_dir>/progress-<tag>.txt, out_dir the second argument (default
# tempdir(), which does not survive the run: pass a directory).
args <- commandArgs(trailingOnly = TRUE)
tag <- args[1L]
stopifnot(!is.na(tag))
out_dir <- if (length(args) >= 2L) args[2L] else tempdir()

source("dev/dev-utils.R")
load_all1(timer = FALSE)

set.seed(3)
X <- matrix(
  runif(40L * 12L),
  40L,
  12L,
  dimnames = list(paste0("s", 1:40), paste0("c", 1:12))
)
X[cbind(sample.int(40L, 30L, TRUE), sample.int(12L, 30L, TRUE))] <- NA

out <- character(0)
for (refresh in c(1L, 7L, 50L)) {
  out <- c(
    out,
    sprintf("== refresh %d ==", refresh),
    capture.output(
      invisible(mipca(
        X,
        rank = 2L,
        transform = "identity",
        chains = 2L,
        iter_warmup = 60L,
        iter_sampling = 60L,
        refresh = refresh,
        seed = 1L
      ))
    )
  )
}
out <- out[!grepl("finished in", out, fixed = TRUE)]
path <- file.path(out_dir, sprintf("progress-%s.txt", tag))
writeLines(out, path, useBytes = TRUE)
cat(sprintf(
  "%d lines to %s; %d progress lines\n",
  length(out),
  path,
  sum(grepl("Iteration:", out, fixed = TRUE))
))
cat(head(out[grepl("Iteration:", out, fixed = TRUE)], 3), sep = "\n")
