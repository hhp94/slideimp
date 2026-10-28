# Compare two dev/bench_wide_phases.R output directories, BEFORE and AFTER a
# change: the whole-call wall time, the per-phase EM timings at each thread
# count and missingness level, and the imputed values each run saved, on max
# abs and max rel error (relative over cells with |before| >= 1e-3).
#
# Usage, from the repository root:
#   Rscript dev/compare_bench_dirs.R <before_dir> <after_dir>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, dir.exists(args[[1L]]), dir.exists(args[[2L]]))
dirs <- c(before = args[[1L]], after = args[[2L]])

# the same phase sets as dev/bench_wide_phases.R
LOOP_PHASES <- c(
  "restandardize", "svd", "form_gram", "eig", "recover_factors",
  "post_svd", "reconstruct", "objective"
)
TOP_PHASES <- c("restandardize", "svd", "post_svd", "reconstruct", "objective")
ONCE_PHASES <- c("gram_cache_init", "postprocessing", "pca_imp_internal_cpp_total")

rd <- function(dir, pat) {
  do.call(rbind, lapply(
    list.files(dir, pattern = pat, full.names = TRUE), read.csv
  ))
}
both <- function(pat) {
  do.call(rbind, lapply(names(dirs), function(w) {
    d <- rd(dirs[[w]], pat)
    d$when <- w
    d
  }))
}

w <- both("^wall_.*\\.csv$")
wall <- merge(
  w[w$when == "before", c("label", "chg", "n_chg", "wall", "cpu")],
  w[w$when == "after", c("label", "chg", "wall", "cpu")],
  by = c("label", "chg"), suffixes = c("_before", "_after")
)
wall$saved_s <- wall$wall_before - wall$wall_after
cat("== whole pca_imp() call (s) ==\n")
print(format(wall[order(wall$chg, wall$label), ], digits = 3), row.names = FALSE)

ph <- both("^phases_.*\\.csv$")
ph$ms <- ph$Microseconds / 1000
ph$se_ms <- ph$SD / sqrt(ph$Count) / 1000
keep <- c("label", "chg", "Name", "Count", "ms", "se_ms")
cmp <- merge(
  ph[ph$when == "before", keep], ph[ph$when == "after", keep],
  by = c("label", "chg", "Name", "Count"), suffixes = c("_before", "_after")
)
cmp$saved_ms <- cmp$ms_before - cmp$ms_after
# saving in units of its own standard error
cmp$diff_in_se <- cmp$saved_ms / sqrt(cmp$se_ms_before^2 + cmp$se_ms_after^2)

for (tag in sort(unique(cmp$chg))) {
  for (lb in sort(unique(cmp$label))) {
    d <- cmp[cmp$label == lb & cmp$chg == tag, ]
    lp <- d[d$Name %in% LOOP_PHASES, ]
    lp <- lp[order(match(lp$Name, LOOP_PHASES)), ]
    top <- lp$Name %in% TOP_PHASES
    cat(sprintf(
      "\n== %s %s: per EM iteration (ms), %d samples each; iteration total %.1f -> %.1f ms ==\n",
      tag, lb, lp$Count[1], sum(lp$ms_before[top]), sum(lp$ms_after[top])
    ))
    print(
      format(
        lp[, c("Name", "ms_before", "ms_after", "saved_ms", "diff_in_se")],
        digits = 3
      ),
      row.names = FALSE
    )
    op <- d[d$Name %in% ONCE_PHASES, ]
    cat("-- once per call (ms) --\n")
    print(
      format(op[, c("Name", "ms_before", "ms_after", "saved_ms")], digits = 3),
      row.names = FALSE
    )

    b <- readRDS(file.path(dirs[["before"]], sprintf("imp_%s_%s.rds", lb, tag)))
    a <- readRDS(file.path(dirs[["after"]], sprintf("imp_%s_%s.rds", lb, tag)))
    stopifnot(length(b) > 0L, length(a) == length(b))
    rel_keep <- abs(b) >= 1e-3
    stopifnot(sum(rel_keep) > 0L)
    cat(sprintf(
      "imputed cells after vs before: max abs %.3g, max rel %.3g over %d of %d cells, bitwise %s\n",
      max(abs(a - b)), max(abs(a - b)[rel_keep] / abs(b[rel_keep])),
      sum(rel_keep), length(b), identical(a, b)
    ))
  }
}
