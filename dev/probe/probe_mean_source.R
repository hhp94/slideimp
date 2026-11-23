# Probe: how mean_imp_col() on parts should take its column means
# (dev/si_mi_conversion.md at d5a34e5, Phase 2, "Open").
#
# 1. Difference. For every column with at least one hole and one observed
#    value, the mean today's C++ kernel writes (read back from its output at a
#    hole) against three means taken in R, one column at a time:
#      sum      sum(x, na.rm = TRUE) / n_obs   long double sum, double divide
#      colmeans .colMeans(x, n, 1, TRUE)      long double sum and divide
#      mean     mean(x[-holes])               long double, two passes
#    Reported: columns compared, columns not bitwise equal, max abs and max rel
#    error (rel over columns whose kernel mean is nonzero), and max distance in
#    ulps of the kernel mean. Then each mean, the kernel's included, against
#    an exact one: each observed value is split into three doubles on fixed
#    power-of-two grids, each piece is summed exactly in double (the script
#    asserts the bounds that make it exact), and the mean is formed from the
#    three sums, so it is within about one ulp of the true mean.
# 2. Speed and allocation. bench::mark of the C++ kernel that mean_imp_col()
#    called up to 1.2.0 (mean_imp_col_internal(), dimnames attached as it
#    did) at 1, 4 and 8 cores, against mean_imp_col() on parts as it ships
#    and the prototype with each of the three R means. mem_alloc is R's own
#    allocations; the kernel's Armadillo output is allocated outside R and
#    then copied into R memory by Rcpp, so only the copy shows.
#
# Inputs: chr21 (dev/GSE286313.EPICv1_chr21.qs2, all-NA probes kept, as
# mean_imp_col() takes them), and a 2000 x 20000 normal matrix with holes in
# 5% of the rows of every other column and in 95% of the rows of the last 10%
# of the columns, timed over every column and over that last 10%.
#   Rscript dev/probe/probe_mean_source.R
source("dev/dev-utils.R")
load_all1(timer = FALSE)

chr21 <- t(qs2::qs_read("dev/GSE286313.EPICv1_chr21.qs2"))

set.seed(1)
n <- 2000L
p <- 20000L
big <- matrix(0, n, p, dimnames = list(paste0("r", 1:n), paste0("c", 1:p)))
offset <- rnorm(p, sd = 100)
for (j1 in seq(1L, p, by = 1000L)) {
  j <- j1:min(j1 + 999L, p)
  big[, j] <- matrix(rnorm(n * length(j)), n) +
    rep(offset[j], each = n)
}
last <- seq.int(p - p %/% 10L + 1L, p)
for (j in seq_len(p)) {
  k <- if (j %in% last) {
    as.integer(0.95 * n)
  } else if (j %% 2L == 0L) {
    as.integer(0.05 * n)
  } else {
    0L
  }
  if (k > 0L) big[sample.int(n, k), j] <- NA
}

# The prototype: one mean part over `cols`, built from the input one column at
# a time. Columns with no observed value are left out of the part, as the
# kernel leaves them.
proto_part <- function(obj, cols, how) {
  n_row <- nrow(obj)
  cols <- sort(as.integer(cols))
  n_miss <- col_scan(obj, cols = cols)$n_miss
  keep <- cols[n_miss < n_row]
  fill <- cols[n_miss > 0 & n_miss < n_row]
  rows <- vector("list", length(fill))
  means <- numeric(length(fill))
  for (k in seq_along(fill)) {
    x <- obj[, fill[[k]]]
    h <- which(is.na(x))
    rows[[k]] <- h
    means[[k]] <- switch(
      how,
      sum = sum(x, na.rm = TRUE) / (n_row - length(h)),
      colmeans = .colMeans(x, n_row, 1L, na.rm = TRUE),
      mean = mean(x[-h])
    )
  }
  len <- lengths(rows)
  map <- cbind(
    row = unlist(rows, use.names = FALSE),
    col = rep.int(fill, len)
  )
  list(
    part = new_part(keep, map, rep.int(means, len), "mean"),
    fill = fill,
    means = means
  )
}
proto_imp <- function(obj, cols, how) {
  fill_parts(obj, list(proto_part(obj, cols, how)$part))
}

## 1. Difference ----

# The mean of `x` to within about one ulp. x = a + b + c exactly, with a on a
# 2^-20 grid, b on a 2^-40 grid and |c| < 2^-40. With |x| < 2^10, fewer than
# 2^11 values and no nonzero |x| below 2^-30, each piece's sum is an integer
# multiple of its grid below 2^53 of it, so sum() returns it exactly.
exact_mean <- function(x) {
  nz <- x[x != 0]
  stopifnot(
    length(x) < 2^11,
    max(abs(x)) < 2^10,
    length(nz) == 0L || min(abs(nz)) >= 2^-30
  )
  a <- trunc(x * 2^20) / 2^20
  r <- x - a
  b <- trunc(r * 2^40) / 2^40
  cc <- r - b
  n <- length(x)
  sum(a) / n + (sum(b) + sum(cc)) / n
}

diff_report <- function(label, obj) {
  cols <- seq_len(ncol(obj))
  ker <- mean_imp_col_internal(obj, cols - 1L, 1L)
  pp <- proto_part(obj, cols, "sum")
  fill <- pp$fill
  first_hole <- vapply(
    fill,
    function(j) which(is.na(obj[, j]))[[1]],
    integer(1)
  )
  k_mean <- ker[cbind(first_hole, fill)]
  stopifnot(length(fill) > 0L, !anyNA(k_mean))
  nz <- k_mean != 0
  stopifnot(any(nz))
  ulp <- 2^(floor(log2(abs(k_mean[nz]))) - 52)
  cat(sprintf(
    "\n%s: %d x %d, %d columns with a hole and an observed value\n",
    label,
    nrow(obj),
    ncol(obj),
    length(fill)
  ))
  cat(sprintf(
    "  %-9s %12s %12s %12s %9s\n",
    "R mean",
    "not bitwise",
    "max abs",
    "max rel",
    "max ulps"
  ))
  r_means <- list()
  for (how in c("sum", "colmeans", "mean")) {
    r_mean <- proto_part(obj, cols, how)$means
    r_means[[how]] <- r_mean
    d <- abs(r_mean - k_mean)
    cat(sprintf(
      "  %-9s %12d %12.3e %12.3e %9.1f\n",
      how,
      sum(r_mean != k_mean),
      max(d),
      max(d[nz] / abs(k_mean[nz])),
      max(d[nz] / ulp)
    ))
  }

  # against the exact mean, in ulps of it
  e_mean <- vapply(
    fill,
    function(j) {
      x <- obj[, j]
      exact_mean(x[!is.na(x)])
    },
    numeric(1)
  )
  enz <- e_mean != 0
  stopifnot(any(enz))
  e_ulp <- 2^(floor(log2(abs(e_mean[enz]))) - 52)
  cat(sprintf(
    "  against the exact mean, over %d columns with a nonzero mean:\n",
    sum(enz)
  ))
  cat(sprintf(
    "  %-9s %12s %12s %12s\n",
    "mean",
    "max ulps",
    "median ulps",
    "cols > 1 ulp"
  ))
  for (how in c("kernel", names(r_means))) {
    m <- if (how == "kernel") k_mean else r_means[[how]]
    u <- abs(m[enz] - e_mean[enz]) / e_ulp
    cat(sprintf(
      "  %-9s %12.1f %12.1f %12d\n",
      how,
      max(u),
      stats::median(u),
      sum(u > 1)
    ))
  }
  # the whole result, mean_imp_col() as it ships against the kernel's
  today <- ker
  dimnames(today) <- dimnames(obj)
  new <- mean_imp_col(obj)
  stopifnot(
    identical(is.na(today), is.na(new)),
    identical(dimnames(today), dimnames(new))
  )
  d <- abs(new - today)
  rel_over <- !is.na(today) & today != 0
  cat(sprintf(
    "  mean_imp_col() against the kernel: max abs %.3e, max rel %.3e over %d nonzero cells, bitwise %s\n",
    max(d, na.rm = TRUE),
    max(d[rel_over] / abs(today[rel_over])),
    sum(rel_over),
    identical(unname(today), unname(new))
  ))
}
diff_report("chr21", chr21)
diff_report("big", big)

## 2. Speed and allocation ----
kernel <- function(obj, cols, cores) {
  out <- mean_imp_col_internal(obj, cols - 1L, cores)
  dimnames(out) <- dimnames(obj)
  out
}
timing <- function(label, obj, cols, iterations) {
  cat(sprintf("\n%s, %d of %d columns:\n", label, length(cols), ncol(obj)))
  res <- bench::mark(
    `kernel, cores 1` = kernel(obj, cols, 1L),
    `kernel, cores 4` = kernel(obj, cols, 4L),
    `kernel, cores 8` = kernel(obj, cols, 8L),
    `mean_imp_col()` = mean_imp_col(obj, subset = cols),
    `parts, sum` = proto_imp(obj, cols, "sum"),
    `parts, colmeans` = proto_imp(obj, cols, "colmeans"),
    `parts, mean` = proto_imp(obj, cols, "mean"),
    check = FALSE,
    iterations = iterations,
    filter_gc = FALSE
  )
  print(res[, c("expression", "min", "median", "mem_alloc", "n_gc")])
}
timing("chr21", chr21, seq_len(ncol(chr21)), 20)
timing("big", big, seq_len(ncol(big)), 5)
timing("big", big, last, 5)

cat(sprintf(
  "\ninput sizes: chr21 %.1f MB, big %.1f MB\n",
  object.size(chr21) / 2^20,
  object.size(big) / 2^20
))
