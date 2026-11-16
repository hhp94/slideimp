# Probe: where pca_imp()'s peak memory sits today, and what the two si_mi
# output designs would cost after the kernel. Each stage runs in a FRESH R
# process (callr), because a process's peak working set only ever grows.
#
# The DLL is built once here with load_all1(timer = FALSE); every child loads
# that same DLL with pkgload::load_all(compile = FALSE), so no stage is
# measured on a different build.
#
# Peak is the process peak working set (ps, Windows), which counts C++
# (Armadillo) allocations that R's gc() figures miss. Each stage reports
#   extra = (peak after the stage - working set before it) / size of x
# so 1.00 means the stage held one more copy of the input on top of it.
# The input is built block by block in place, so building it does not set
# the peak; `masked` says when the stage stayed under the build's own peak.
#
# Run from the package root:
#   Rscript dev/probe_pca_imp_peak.R

source("dev/dev-utils.R")
load_all1(timer = FALSE)

n <- 2000L
p <- 20000L
frac_mean <- 0.1 # columns over colmax, so post_imp mean-fills them

stage <- function(root, scenario, n, p, frac_mean) {
  suppressMessages(pkgload::load_all(root, compile = FALSE, quiet = TRUE))
  mem <- function() ps::ps_memory_full_info(ps::ps_handle())
  set.seed(1)
  ## The input, low rank plus noise, written block by block into one matrix.
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
  ## Holes: 5% of the rows in every other column, and 95% of the rows in the
  ## last `frac_mean` of the columns (over the default colmax of 0.9).
  inel <- if (frac_mean > 0) {
    seq.int(p - as.integer(p * frac_mean) + 1L, p)
  } else {
    integer(0)
  }
  elig <- setdiff(seq_len(p), inel)
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
  hole_lin <- unlist(hole, use.names = FALSE)
  rm(hole, U, V)
  invisible(gc())
  X <- as.numeric(object.size(x))
  before <- mem()

  if (scenario == "pca_imp, post_imp = TRUE") {
    res <- suppressWarnings(pca_imp(
      x,
      ncp = 2,
      solver = "exact",
      maxiter = 10,
      miniter = 1,
      threshold = 1e-4
    ))
  } else if (scenario == "pca_imp, post_imp = FALSE") {
    res <- suppressWarnings(pca_imp(
      x,
      ncp = 2,
      solver = "exact",
      maxiter = 10,
      miniter = 1,
      threshold = 1e-4,
      post_imp = FALSE
    ))
  } else if (scenario == "today's post step alone") {
    ## pca_imp()'s write-back, then mean_imp_col() on the ineligible columns.
    obj <- x
    obj[hole_lin] <- 0.5
    obj[, inel] <- x[, inel] # put the ineligible holes back
    obj <- mean_imp_col(obj, subset = inel)
  } else if (scenario == "parts, matrix door") {
    ## Means read from x one column at a time, then one copy of x filled.
    means <- vapply(inel, function(j) mean(x[, j], na.rm = TRUE), 0)
    vals <- rep(0.5, length(hole_lin))
    col_of <- (hole_lin - 1) %/% n + 1
    m <- match(col_of, inel)
    vals[!is.na(m)] <- means[m[!is.na(m)]]
    out <- x
    out[hole_lin] <- vals
  } else if (scenario %in% c("slices, si_mi door", "slices + as.matrix")) {
    part <- function(cols) {
      sub <- x[, cols, drop = FALSE]
      lin <- hole_lin[((hole_lin - 1) %/% n + 1) %in% cols]
      col_g <- (lin - 1) %/% n + 1
      map <- cbind(
        row = as.integer((lin - 1) %% n + 1),
        col = match(col_g, cols)
      )
      d <- .diag_frame(
        map = map[integer(0), , drop = FALSE],
        draws = list(matrix(NA_real_, 0L, 1L)),
        p = ncol(sub),
        nms = colnames(sub),
        prop_diag = 0.5
      )
      .si_mi_new(
        chunks = list(list(
          obj = sub,
          map = map,
          draws = matrix(0.5, nrow(map), 1L),
          transform = list(name = "identity", backward = function(z, ...) z),
          diag = d,
          prov = list(
            chains = 1L,
            thin = 1L,
            ndraws = 1L,
            thin_final = 1L,
            thin_effective = 1L,
            seam_state = "never"
          ),
          call = quote(probe())
        )),
        rownames = rownames(x),
        nkeep = 1L,
        call = quote(probe())
      )
    }
    fit <- .si_mi_join(
      list(part(elig), part(inel)),
      call = quote(probe()),
      placement = c(elig, inel)
    )
    if (scenario == "slices + as.matrix") {
      out <- as.matrix(fit)
    }
  }
  after <- mem()
  ## ps returns a named numeric vector, not a list.
  list(
    X_mb = X / 2^20,
    extra = (after[["peak_wset"]] - before[["wset"]]) / X,
    masked = after[["peak_wset"]] <= before[["peak_wset"]]
  )
}

scenarios <- c(
  "pca_imp, post_imp = TRUE",
  "pca_imp, post_imp = FALSE",
  "today's post step alone",
  "parts, matrix door",
  "slices, si_mi door",
  "slices + as.matrix"
)
cat(sprintf(
  "x: %d x %d, %.0f%% of columns over colmax\n",
  n,
  p,
  100 * frac_mean
))
for (s in scenarios) {
  r <- callr::r(
    stage,
    args = list(
      root = getwd(),
      scenario = s,
      n = n,
      p = p,
      frac_mean = frac_mean
    )
  )
  cat(sprintf(
    "%-28s x = %4.0f MB   extra peak = %.2f x input%s\n",
    s,
    r$X_mb,
    r$extra,
    if (r$masked) "   (MASKED by the build)" else ""
  ))
}
