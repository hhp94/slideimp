# Fixed per-call overhead of the K-NN chain at cores = 1.
#
# Two costs are under measurement here:
#   1. RcppThread::parallelFor() resizing the global thread pool twice per
#      .Call (down to `cores`, then back), which quickpool implements by
#      joining every worker and spawning new ones.
#   2. mean_imp_col() running in full after the kernel even when the kernel
#      left no missing value behind.
# Both are fixed costs, so they are measured on window-shaped inputs where the
# real arithmetic is small enough for them to matter.
#
# Run it the same way before and after a change and compare the two tables.
# Usage, from the repository root: Rscript dev/bench_knn_overhead.R [label]

if (!file.exists("DESCRIPTION") || !file.exists("R/dev-utils.R")) {
  stop("run this from the repository root: Rscript dev/bench_knn_overhead.R")
}
source("R/dev-utils.R")
suppressMessages(load_all1(timer = FALSE))

label <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "unlabelled"

# Median wall time per call, in milliseconds. The costs under measurement are
# around a millisecond, which is the resolution of the clock, so each timed
# block runs `inner` calls and divides.
time_ms <- function(expr, reps = 11L, inner = 20L) {
  e <- substitute(expr)
  env <- parent.frame()
  for (i in seq_len(3L)) {
    eval(e, env)
  }
  t <- numeric(reps)
  for (i in seq_len(reps)) {
    a <- Sys.time()
    for (j in seq_len(inner)) {
      eval(e, env)
    }
    t[i] <- as.numeric(difftime(Sys.time(), a, units = "secs")) * 1000 / inner
  }
  stats::median(t)
}

mk <- function(n, p, seed = 1L, perc = 0.5) {
  set.seed(seed)
  sim_mat(n, p, perc_col_na = perc)$input
}

rows <- list()
add <- function(what, ms) {
  rows[[length(rows) + 1L]] <<- data.frame(what = what, ms = round(ms, 3))
}

# --- fixed overhead of the parallel .Calls, with as little work as possible ---
m_small <- mk(100, 40)
grp <- which(colSums(is.na(m_small)) > 0L) - 1L
cpl <- which(colSums(is.na(m_small)) == 0L) - 1L
stopifnot(length(grp) > 0L, length(cpl) > 0L)

add(
  "impute_knn_brute .Call, 100x40, cores=1",
  time_ms(
    impute_knn_brute(
      obj = m_small,
      k = 5L,
      grp_impute = as.integer(grp),
      grp_miss_no_imp = integer(0),
      grp_complete = as.integer(cpl),
      method = 0L,
      dist_pow = 0,
      cores = 1L,
      pb = FALSE
    )
  )
)

set.seed(7)
m_full <- matrix(stats::rnorm(100 * 40), 100, 40)
stopifnot(!anyNA(m_full))
add(
  "mean_imp_col, 100x40, nothing to do",
  time_ms(
    mean_imp_col(m_full, subset = seq_len(ncol(m_full)), cores = 1L)
  )
)
add("col_vars, 100x40", time_ms(col_vars(m_full, cores = 1L)))
add(
  "col_miss (serial .Call), 100x40",
  time_ms(mat_miss(m_full, col = TRUE, prop = FALSE))
)
add("check_finite (serial .Call), 100x40", time_ms(check_finite(m_full)))

# --- whole-call timings on window-shaped inputs ---
for (dims in list(c(100, 40), c(500, 60), c(200, 300))) {
  m <- mk(dims[1], dims[2], seed = dims[2])
  nm <- sprintf("knn_imp %dx%d, k=5, cores=1", dims[1], dims[2])
  add(
    paste0(nm, ", post_imp=TRUE"),
    time_ms(
      knn_imp(m, k = 5L, cores = 1L, post_imp = TRUE, .progress = FALSE),
      inner = 5L
    )
  )
  add(
    paste0(nm, ", post_imp=FALSE"),
    time_ms(
      knn_imp(m, k = 5L, cores = 1L, post_imp = FALSE, .progress = FALSE),
      inner = 5L
    )
  )
}

out <- do.call(rbind, rows)
cat("\n==== knn overhead, label:", label, "====\n")
print(out, row.names = FALSE, right = FALSE)
cat("\nR", as.character(getRversion()), "|", Sys.info()[["nodename"]], "\n")
