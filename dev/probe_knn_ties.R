# Does the order of equal-distance neighbors change what knn_imp() returns?
#
# distance_vector_impl() (src/impute_knn_brute.cpp) appends the first k
# candidates in scan order, sorts them by distance, and from then on lets a
# strictly closer candidate evict top_k.back(). When several fill entries tie
# at the worst distance, the sort decides which of them sits at the back, so
# the sort's handling of ties decides which tied neighbor is dropped.
#
# The tie rule dev/knn_imp_structure.md states is "first-encountered wins"
# under the kernel's scan order: complete columns (group 3) in ascending
# column order, then the columns with missing values in ascending column
# order, minus the target. This script checks the kernel against that rule.
#
# Part 1 builds cases where the chosen neighbor SET can be read back out of
# the result. Each fill candidate is the target plus a +-1 pattern on the
# target's observed rows, so its distance is exactly 1 under both metrics,
# and candidate j holds 2^(j - 1) at the target's one missing row. Candidates
# placed after the fill in scan order use a +-1/2 pattern (distance 1/4
# euclidean, 1/2 manhattan), so each one evicts a tied fill entry. With
# dist_pow = 0 the imputed value times k is the sum of 2^(j - 1) over the
# chosen neighbors, exactly, and its binary digits are the chosen set.
#
# Part 2 sweeps random small-integer matrices, where exact ties are common,
# and compares every imputed cell against a pure-R reference applying the
# scan-order rule. With integer data and dist_pow = 0 every sum on both sides
# is exact, so a cell matches bitwise when both sides chose neighbors with the
# same values. The same cells are also compared against a "leftmost column
# wins" rule - what knn_ref() in tests/testthat/test-knn_imp.R gets from a
# stable order() over column indices - to show the two are different rules.
#
# Part 3 times knn_imp() on tie-free data where the per-column sort is as
# large a share of the work as it gets, so a change to the sort can be costed.
#
# Measured with GCC 14.2 (libstdc++) on Windows. With std::sort in the fill
# phase the kernel followed the scan-order rule at every k <= 16 tried (5, 10,
# 16) and broke it at every k >= 17 tried (17, 20, 30, 32), in all three
# Part 1 layouts, by evicting tied candidates from the front of the run; in
# Part 2 it changed 102 to 611 of 2400 cells at k >= 17. With std::stable_sort
# every row passes, bitwise, and an A/B/A run of Part 3 put the two sorts
# within run-to-run noise of each other.
#
# Usage, from the repository root: Rscript dev/probe_knn_ties.R [label]
# Exits non-zero if the kernel departs from the scan-order rule anywhere.

if (!file.exists("DESCRIPTION") || !file.exists("R/dev-utils.R")) {
  stop("run this from the repository root: Rscript dev/probe_knn_ties.R")
}
source("R/dev-utils.R")
suppressMessages(load_all1(timer = FALSE))

label <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "unlabelled"
TOL <- 1e-12
failures <- 0L

# ---------------------------------------------------------------------------
# Reference
# ---------------------------------------------------------------------------

# The kernel's scan order when knn_imp() is called with no `subset` and
# colmax = 1: every column with a missing value is in grp_impute, every other
# column in grp_complete, and group 3 is scanned first.
scan_order <- function(obj, target) {
  has_miss <- colSums(is.na(obj)) > 0L
  c(which(!has_miss), setdiff(which(has_miss), target))
}

leftmost_order <- function(obj, target) {
  setdiff(seq_len(ncol(obj)), target)
}

acc_euclidean <- function(x) sum(x^2)
acc_manhattan <- function(x) sum(abs(x))

# Distance from `target` to each column of `cand`: the mean per-row
# contribution over rows observed in both, +Inf when there are none.
ref_dist <- function(obj, target, cand, acc) {
  x <- obj[, target]
  obs <- !is.na(x)
  vapply(
    cand,
    function(j) {
      both <- obs & !is.na(obj[, j])
      n <- sum(both)
      if (n == 0L) Inf else acc(x[both] - obj[both, j]) / n
    },
    numeric(1)
  )
}

# order() leaves ties in their input order, so the order of `cand` IS the tie
# rule. Returns original column indices of the chosen neighbors, nearest first.
ref_neighbors <- function(obj, target, k, acc, cand) {
  d <- ref_dist(obj, target, cand, acc)
  ord <- order(d)[seq_len(k)]
  ord <- ord[is.finite(d[ord])]
  cand[ord]
}

# A column's neighbor set is decided by distances alone unless the k-th and
# (k + 1)-th smallest finite distances are equal.
tie_sensitive <- function(obj, target, k, acc) {
  d <- sort(ref_dist(obj, target, leftmost_order(obj, target), acc))
  d <- d[is.finite(d)]
  length(d) > k && d[k] == d[k + 1L]
}

# dist_pow = 0 imputation of every missing cell, with neighbors chosen under
# the tie rule `order_fun`. A neighbor missing at a row drops out of that row.
ref_impute <- function(obj, k, acc, order_fun) {
  out <- obj
  for (t in which(colSums(is.na(obj)) > 0L)) {
    nn <- ref_neighbors(obj, t, k, acc, order_fun(obj, t))
    for (r in which(is.na(obj[, t]))) {
      v <- obj[r, nn]
      ok <- !is.na(v)
      out[r, t] <- if (any(ok)) sum(v[ok]) / sum(ok) else NA_real_
    }
  }
  out
}

run_kernel <- function(obj, k, method) {
  r <- knn_imp(
    obj,
    k = k,
    method = method,
    colmax = 1,
    dist_pow = 0,
    post_imp = FALSE,
    na_check = FALSE
  )
  unclass(r)[, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Part 1 - read the chosen set back out of the imputed value
# ---------------------------------------------------------------------------

# Column 1 is the target, missing at row 1 only. Columns 2.. are candidates,
# candidate j in column j + 1. `kind` gives each candidate's group ("c" for
# complete, "m" for one missing cell, in a row the target observes) and
# `closer` marks the ones built at the smaller distance.
build_case <- function(kind, closer, seed) {
  set.seed(seed)
  n_obs <- 8L
  m <- length(kind)
  t_obs <- sample(-5:5, n_obs, replace = TRUE)
  sgn <- matrix(sample(c(-1, 1), n_obs * m, replace = TRUE), n_obs, m)
  step <- ifelse(closer, 0.5, 1)
  cand <- t_obs + sweep(sgn, 2, step, "*")
  for (j in which(kind == "m")) {
    cand[1L + (j %% n_obs), j] <- NA
  }
  obj <- rbind(c(NA, 2^(seq_len(m) - 1)), cbind(t_obs, cand))
  dimnames(obj) <- NULL
  obj
}

decode <- function(value, k, m) {
  s <- round(value * k)
  which(floor(s / 2^(seq_len(m) - 1)) %% 2 == 1) + 1L # back to columns
}

# Three layouts, all with exactly k tied candidates first in scan order and
# n_closer closer candidates after them:
#   complete - every candidate in group 3, scan order is column order
#   masked   - every candidate in group 1, scan order is column order
#   mixed    - tied candidates alternate complete / masked by column, so the
#              scan order (complete first) is NOT the column order, and the
#              closer ones are masked, at the end
layout_case <- function(layout, k, n_closer) {
  m <- k + n_closer
  closer <- c(rep(FALSE, k), rep(TRUE, n_closer))
  kind <- switch(
    layout,
    complete = rep("c", m),
    masked = rep("m", m),
    mixed = c(rep_len(c("c", "m"), k), rep("m", n_closer))
  )
  list(kind = kind, closer = closer)
}

fmt <- function(x) if (length(x)) paste(x, collapse = ",") else "-"

cat("\n==== Part 1: chosen neighbor set, read back from the imputed value ====\n")
p1 <- list()
for (layout in c("complete", "masked", "mixed")) {
  for (method in c("euclidean", "manhattan")) {
    acc <- if (method == "euclidean") acc_euclidean else acc_manhattan
    for (k in c(5L, 16L, 17L, 20L, 32L)) {
      n_closer <- 4L
      lc <- layout_case(layout, k, n_closer)
      obj <- build_case(lc$kind, lc$closer, seed = k)
      m <- length(lc$kind)

      # the case has to be the case it claims to be: exact ties at the fill
      # boundary, and closer candidates that the fill does not reach
      cand_scan <- scan_order(obj, 1L)
      d_scan <- ref_dist(obj, 1L, cand_scan, acc)
      stopifnot(
        length(unique(d_scan[seq_len(k)])) == 1L,
        all(d_scan[-seq_len(k)] < d_scan[1L]),
        sum(is.na(obj[, 1L])) == 1L
      )

      got <- run_kernel(obj, k, method)[1L, 1L]
      chosen <- decode(got, k, m)
      tied_cols <- cand_scan[seq_len(k)]
      evicted <- sort(setdiff(tied_cols, chosen))
      want_scan <- sort(setdiff(
        tied_cols,
        ref_neighbors(obj, 1L, k, acc, cand_scan)
      ))
      want_left <- sort(setdiff(
        tied_cols,
        ref_neighbors(obj, 1L, k, acc, leftmost_order(obj, 1L))
      ))
      ok <- identical(evicted, want_scan) && length(chosen) == k
      if (!ok) {
        failures <- failures + 1L
      }
      p1[[length(p1) + 1L]] <- data.frame(
        layout = layout,
        method = substr(method, 1, 3),
        k = k,
        kernel_evicted = fmt(evicted),
        scan_rule = fmt(want_scan),
        leftmost_rule = fmt(want_left),
        result = if (ok) "PASS" else "FAIL"
      )
    }
  }
}
print(do.call(rbind, p1), row.names = FALSE, right = FALSE)

# ---------------------------------------------------------------------------
# Part 2 - random small-integer matrices
# ---------------------------------------------------------------------------

sweep_case <- function(seed, n_row = 10L, n_col = 60L) {
  set.seed(seed)
  obj <- matrix(
    as.double(sample(0:3, n_row * n_col, replace = TRUE)),
    n_row,
    n_col
  )
  for (j in sample(n_col, n_col / 2)) {
    obj[sample(n_row, 2L), j] <- NA
  }
  obj
}

differs <- function(a, b) {
  xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b)
}

cat("\n==== Part 2: random 10 x 60 matrices, values in 0:3, 40 seeds ====\n")
seeds <- 1:40
p2 <- list()
for (method in c("euclidean", "manhattan")) {
  acc <- if (method == "euclidean") acc_euclidean else acc_manhattan
  for (k in c(5L, 10L, 16L, 17L, 20L, 30L)) {
    n_cells <- 0L
    n_cols <- 0L
    n_sens <- 0L
    bad_scan <- 0L
    bad_scan_cols <- 0L
    bad_left <- 0L
    abs_err <- 0
    rel_err <- 0
    bitwise <- TRUE
    for (s in seeds) {
      obj <- sweep_case(s)
      na <- is.na(obj)
      targets <- which(colSums(na) > 0L)
      got <- run_kernel(obj, k, method)
      want <- ref_impute(obj, k, acc, scan_order)
      left <- ref_impute(obj, k, acc, leftmost_order)
      both <- na & !is.na(got) & !is.na(want)
      d_scan <- differs(got, want) & na
      n_cells <- n_cells + sum(na)
      n_cols <- n_cols + length(targets)
      n_sens <- n_sens + sum(vapply(
        targets,
        function(t) tie_sensitive(obj, t, k, acc),
        logical(1)
      ))
      bad_scan <- bad_scan + sum(d_scan)
      bad_scan_cols <- bad_scan_cols + sum(colSums(d_scan) > 0L)
      bad_left <- bad_left + sum(differs(got, left) & na)
      if (any(both)) {
        e <- abs(got[both] - want[both])
        abs_err <- max(abs_err, e)
        big <- abs(want[both]) > 1e-8
        if (any(big)) {
          rel_err <- max(rel_err, e[big] / abs(want[both][big]))
        }
      }
      bitwise <- bitwise && identical(got[na], want[na])
    }
    ok <- bad_scan == 0L && abs_err < TOL && rel_err < TOL
    if (!ok) {
      failures <- failures + 1L
    }
    p2[[length(p2) + 1L]] <- data.frame(
      method = substr(method, 1, 3),
      k = k,
      cells = n_cells,
      target_cols = n_cols,
      tie_sensitive_cols = n_sens,
      cols_off_scan = bad_scan_cols,
      cells_off_scan = bad_scan,
      cells_off_leftmost = bad_left,
      max_abs = signif(abs_err, 3),
      max_rel = signif(rel_err, 3),
      bitwise = bitwise,
      result = if (ok) "PASS" else "FAIL"
    )
  }
}
print(do.call(rbind, p2), row.names = FALSE, right = FALSE)
cat(
  "max_abs / max_rel are over cells both sides filled; max_rel over those",
  "with |reference| > 1e-8.\n"
)

# ---------------------------------------------------------------------------
# Part 3 - cost
# ---------------------------------------------------------------------------

# Median wall time per call over `reps` timed blocks of `inner` calls each.
# Sys.time(), not proc.time(): the latter's elapsed clock ticks in 10 ms steps
# on Windows, which at 2 calls per block hides anything under 5 ms a call.
time_ms <- function(expr, reps = 15L, inner = 5L) {
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

# Every column is a target, so the sort runs once per column. On 10 x 60 a
# target has only 59 candidates of 10 rows each, so at k = 50 the sort is as
# large a share of the work as it gets; 20 x 2000 is the ordinary case where
# the distance loop dominates.
cat("\n==== Part 3: knn_imp() wall time, tie-free data, cores = 1 ====\n")
p3 <- list()
for (sh in list(c(10L, 60L, 1000L), c(20L, 2000L, 3L))) {
  set.seed(sh[2])
  obj <- matrix(stats::rnorm(sh[1] * sh[2]), sh[1], sh[2])
  for (j in seq_len(sh[2])) {
    obj[sample(sh[1], 2L), j] <- NA
  }
  for (k in c(5L, 20L, 50L)) {
    p3[[length(p3) + 1L]] <- data.frame(
      shape = sprintf("%dx%d", sh[1], sh[2]),
      k = k,
      ms = round(
        time_ms(run_kernel(obj, k, "euclidean"), reps = 11L, inner = sh[3]),
        4
      )
    )
  }
}
print(do.call(rbind, p3), row.names = FALSE, right = FALSE)

cat("\nlabel:", label, "|", R.version.string, "|", Sys.info()[["nodename"]], "\n")
cxx <- tryCatch(
  pkgbuild::with_build_tools(
    system2("g++", "--version", stdout = TRUE)[1L],
    required = FALSE
  ),
  error = function(e) "unknown"
)
cat("compiler:", cxx, "\n")
cat(sprintf("\n%d failing row(s)\n", failures))
quit(status = if (failures > 0L) 1L else 0L)
