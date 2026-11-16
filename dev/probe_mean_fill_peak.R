# Probe: the memory of the mean fill, piece by piece, on the input of
# dev/probe_pca_imp_peak.R (2000 x 20000, holes in 5% of the rows of every
# other column and in 95% of the rows of the last 10% of the columns, which
# are the ones filled).
#
# Two measures:
# - R's own allocations summed over one call (bench's mem_alloc), in this
#   process. Deterministic, but blind to the kernel's Armadillo output, which
#   is allocated outside R (only Rcpp's copy of it into R shows).
# - The peak: each stage in a fresh R process (callr) against the DLL built
#   here, reporting
#     extra = (peak working set after the stage - working set before) / size of x
#   as dev/probe_pca_imp_peak.R does. It counts every allocation, but it also
#   counts garbage R has not yet collected, so it moves with the timing of
#   R's collections. Anything built before `before` (the input, a prebuilt
#   part) is not counted.
#   Rscript dev/probe_mean_fill_peak.R
source("dev/dev-utils.R")
load_all1(timer = FALSE)

n <- 2000L
p <- 20000L
frac_mean <- 0.1

build_input <- function(n, p, frac_mean) {
  set.seed(1)
  r <- 3L
  U <- matrix(rnorm(n * r), n, r)
  V <- matrix(rnorm(p * r), p, r)
  x <- matrix(0, n, p)
  for (j1 in seq(1L, p, by = 1000L)) {
    j <- j1:min(j1 + 999L, p)
    x[, j] <- U %*%
      t(V[j, , drop = FALSE]) +
      matrix(rnorm(n * length(j), sd = 0.5), n, length(j))
  }
  dimnames(x) <- list(paste0("r", seq_len(n)), paste0("c", seq_len(p)))
  inel <- seq.int(p - as.integer(p * frac_mean) + 1L, p)
  hole <- vector("list", p)
  for (j in seq_len(p)) {
    k <- if (j %in% inel) {
      as.integer(0.95 * n)
    } else if (j %% 2L == 0L) {
      as.integer(0.05 * n)
    } else {
      0L
    }
    if (k > 0L) {
      rows <- sort(sample.int(n, k))
      x[rows, j] <- NA
      hole[[j]] <- rows + (j - 1) * n
    }
  }
  list(x = x, inel = inel, hole_lin = unlist(hole, use.names = FALSE))
}

stage <- function(root, scenario, n, p, frac_mean, build_input) {
  suppressMessages(pkgload::load_all(root, compile = FALSE, quiet = TRUE))
  mem <- function() ps::ps_memory_full_info(ps::ps_handle())
  input <- build_input(n, p, frac_mean)
  x <- input$x
  inel <- input$inel
  hole_lin <- input$hole_lin
  rm(input)
  part <- NULL
  if (scenario == "fill_parts() alone, part prebuilt") {
    part <- fallback_part(x, inel)
  }
  invisible(gc())
  X <- as.numeric(object.size(x))
  before <- mem()

  if (scenario == "kernel, as 1.2.0 did") {
    out <- mean_imp_col_internal(x, inel - 1L, 1L)
    dimnames(out) <- dimnames(x)
  } else if (scenario == "mean_imp_col()") {
    out <- mean_imp_col(x, subset = inel)
  } else if (scenario == "fallback_part() alone") {
    part <- fallback_part(x, inel)
  } else if (scenario == "fill_parts() alone, part prebuilt") {
    out <- fill_parts(x, list(part))
  } else if (scenario == "write-back alone") {
    obj <- x
    obj[hole_lin] <- 0.5
    obj[, inel] <- x[, inel]
  } else if (scenario == "write-back, mean_imp_col()") {
    obj <- x
    obj[hole_lin] <- 0.5
    obj[, inel] <- x[, inel]
    obj <- mean_imp_col(obj, subset = inel)
  }
  after <- mem()
  list(
    X_mb = X / 2^20,
    extra = (after[["peak_wset"]] - before[["wset"]]) / X,
    masked = after[["peak_wset"]] <= before[["peak_wset"]]
  )
}

## R's allocations ----
local({
  input <- build_input(n, p, frac_mean)
  x <- input$x
  inel <- input$inel
  part <- fallback_part(x, inel)
  cat(sprintf(
    "x: %d x %d, %.0f MB; filled: the last %d columns, %d holes\n",
    n,
    p,
    as.numeric(object.size(x)) / 2^20,
    length(inel),
    nrow(part$map)
  ))
  res <- bench::mark(
    `kernel, as 1.2.0 did` = {
      out <- mean_imp_col_internal(x, inel - 1L, 1L)
      dimnames(out) <- dimnames(x)
      out
    },
    `mean_imp_col()` = mean_imp_col(x, subset = inel),
    `fallback_part() alone` = fallback_part(x, inel),
    `fill_parts() alone, part prebuilt` = fill_parts(x, list(part)),
    check = FALSE,
    iterations = 3,
    filter_gc = FALSE
  )
  res$x_mult <- as.numeric(res$mem_alloc) / as.numeric(object.size(x))
  cat("\nallocated by R over one call:\n")
  print(res[, c("expression", "median", "mem_alloc", "x_mult")])
})

## The peak ----
scenarios <- c(
  "kernel, as 1.2.0 did",
  "mean_imp_col()",
  "fallback_part() alone",
  "fill_parts() alone, part prebuilt",
  "write-back alone",
  "write-back, mean_imp_col()"
)
cat("\npeak working set, fresh process per stage:\n")
for (s in scenarios) {
  r <- callr::r(
    stage,
    args = list(
      root = getwd(),
      scenario = s,
      n = n,
      p = p,
      frac_mean = frac_mean,
      build_input = build_input
    )
  )
  cat(sprintf(
    "%-36s extra peak = %.2f x input%s\n",
    s,
    r$extra,
    if (r$masked) "   (MASKED by the build)" else ""
  ))
}
