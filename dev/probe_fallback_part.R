# Probe: fallback_part(), vectorized over its columns, against what it
# replaced and against base R.
#
# 1. Mean. fallback_part(method = "mean") against the per-column mean part of
#    Phase 2 (kept below as `percol_part()`, .colMeans() of one column at a
#    time): map and values compared with identical(), since both are built in
#    this process from the same input.
# 2. Median. fallback_part(method = "median"), which takes
#    collapse::fmedian(), against stats::median() of each column's observed
#    values: columns not bitwise equal, max abs and max rel error (rel over
#    columns whose median is nonzero).
# 3. Speed and R allocations, bench::mark, of the two mean parts and the
#    median part, and of stats::median() over the columns by apply().
#
# Inputs: chr21 (dev/GSE286313.EPICv1_chr21.qs2, all-NA probes dropped from
# the columns taken, as mean_imp_col() drops them), and a 2000 x 20000 normal
# matrix with column offsets of sd 100, holes in 5% of the rows of every
# other column and in 95% of the rows of the last 10% of the columns, taken
# over every column and over that last 10%. Each input has even and odd
# observed counts, so both branches of the median are reached; the script
# prints how many columns of each.
#   Rscript dev/probe_fallback_part.R
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
  big[, j] <- matrix(rnorm(n * length(j)), n) + rep(offset[j], each = n)
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

# Phase 2's mean part: one column at a time, sized from the hole counts.
percol_part <- function(obj, cols) {
  n_miss <- col_miss_internal(obj)
  n_row <- nrow(obj)
  n_hole <- as.vector(n_miss[cols])
  fill <- cols[n_hole > 0]
  n_hole <- n_hole[n_hole > 0]
  end <- cumsum(n_hole)
  map <- matrix(0L, sum(n_hole), 2L, dimnames = list(NULL, c("row", "col")))
  map[, 2L] <- rep.int(fill, n_hole)
  values <- numeric(nrow(map))
  for (k in seq_along(fill)) {
    first <- (fill[[k]] - 1) * as.double(n_row)
    x <- obj[seq.int(first + 1, first + n_row)]
    at <- seq.int(end[[k]] - n_hole[[k]] + 1, end[[k]])
    map[at, 1L] <- which(is.na(x))
    values[at] <- .colMeans(x, n_row, 1L, na.rm = TRUE)
  }
  new_part(cols = cols, map = map, values = values, kind = "mean")
}

usable <- function(obj, cols) {
  n_miss <- col_miss_internal(obj)
  cols[n_miss[cols] < nrow(obj)]
}

cases <- list(
  chr21 = list(x = chr21, cols = usable(chr21, seq_len(ncol(chr21)))),
  `big, every column` = list(x = big, cols = usable(big, seq_len(p))),
  `big, last 10%` = list(x = big, cols = usable(big, last))
)

for (nm in names(cases)) {
  x <- cases[[nm]]$x
  cols <- cases[[nm]]$cols
  cat("\n==", nm, "==", length(cols), "columns\n")

  ## 1. mean
  new <- fallback_part(x, cols, "mean")
  old <- percol_part(x, cols)
  cat(sprintf(
    "mean: %d holes; map identical %s; values identical %s\n",
    nrow(new$map),
    identical(new$map, old$map),
    identical(new$values, old$values)
  ))
  invisible(validate_part(new, nrow(x), ncol(x)))

  ## 2. median
  med <- fallback_part(x, cols, "median")
  invisible(validate_part(med, nrow(x), ncol(x)))
  stopifnot(identical(med$map, new$map))
  fm <- as.vector(collapse::fmedian(
    x[, cols, drop = FALSE],
    na.rm = TRUE,
    use.g.names = FALSE
  ))
  sm <- apply(x[, cols, drop = FALSE], 2L, stats::median, na.rm = TRUE)
  sm <- unname(sm)
  n_obs <- nrow(x) - as.vector(col_miss_internal(x)[cols])
  nz <- sm != 0
  stopifnot(any(nz))
  cat(sprintf(
    paste0(
      "median: observed counts even in %d columns, odd in %d; ",
      "not bitwise in %d; max abs %.3g; max rel %.3g (over %d nonzero)\n"
    ),
    sum(n_obs %% 2L == 0L),
    sum(n_obs %% 2L == 1L),
    sum(fm != sm),
    max(abs(fm - sm)),
    max(abs(fm[nz] - sm[nz]) / abs(sm[nz])),
    sum(nz)
  ))
  # the part's values are the medians at the holes
  stopifnot(identical(
    med$values,
    fm[match(med$map[, "col"], cols)]
  ))

  ## 3. speed and allocation
  res <- bench::mark(
    `fallback_part(mean)` = fallback_part(x, cols, "mean"),
    `per-column mean part` = percol_part(x, cols),
    `fallback_part(median)` = fallback_part(x, cols, "median"),
    `apply(stats::median)` = apply(
      x[, cols, drop = FALSE],
      2L,
      stats::median,
      na.rm = TRUE
    ),
    check = FALSE,
    iterations = 5,
    filter_gc = FALSE
  )
  res$x_mult <- as.numeric(res$mem_alloc) / as.numeric(object.size(x))
  print(res[, c("expression", "min", "median", "mem_alloc", "x_mult")])
}
