# What dev/GSE286313.EPICv1_chr21.qs2 holds, so CLAUDE.md can describe it
# from a measurement: class, storage, dimensions, dimnames, value range and
# missingness, as stored. Lines are prefixed "RUN|".
# Usage, from the repository root: Rscript dev/probe_epic_chr21_shape.R
if (!file.exists("DESCRIPTION")) stop("run from the repository root")
f <- "dev/GSE286313.EPICv1_chr21.qs2"
t0 <- proc.time()
d <- qs2::qs_read(f)
cat(sprintf("RUN| read: %.2f s, file %.1f MB\n",
            (proc.time() - t0)[["elapsed"]], file.size(f) / 2^20))
cat("RUN| class:", paste(class(d), collapse = ", "), " storage:", storage.mode(d),
    " is.matrix:", is.matrix(d), "\n")
if (!is.matrix(d)) {
  cat("RUN| names:", paste(utils::head(names(d), 10), collapse = ", "), "\n")
  str(d, max.level = 1)
  quit(status = 0)
}
cat(sprintf("RUN| dim: %d x %d, %.1f MB in memory\n",
            nrow(d), ncol(d), as.numeric(object.size(d)) / 2^20))
cat("RUN| rownames:", paste(utils::head(rownames(d), 4), collapse = ", "), "\n")
cat("RUN| colnames:", paste(utils::head(colnames(d), 4), collapse = ", "), "\n")
cat("RUN| attributes:", paste(names(attributes(d)), collapse = ", "), "\n")
rng <- range(d, na.rm = TRUE)
cat(sprintf("RUN| value range: %.4g to %.4g, non-finite non-NA cells: %d\n",
            rng[1], rng[2], sum(is.infinite(d))))
na_cells <- sum(is.na(d))
cat(sprintf("RUN| NA cells: %d (%.3f%%)\n", na_cells, 100 * na_cells / length(d)))
rm_ <- rowSums(is.na(d))
cm_ <- colSums(is.na(d))
cat(sprintf("RUN| rows with any NA: %d of %d, max row NA fraction %.3f\n",
            sum(rm_ > 0), nrow(d), max(rm_) / ncol(d)))
cat(sprintf("RUN| cols with any NA: %d of %d, max col NA fraction %.3f\n",
            sum(cm_ > 0), ncol(d), max(cm_) / nrow(d)))
cat(sprintf("RUN| rows over 90%% NA: %d, cols over 90%% NA: %d\n",
            sum(rm_ / ncol(d) > 0.9), sum(cm_ / nrow(d) > 0.9)))
