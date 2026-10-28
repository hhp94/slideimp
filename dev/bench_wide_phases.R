# Per-phase EM timings for one wide group at an EPICv2-like shape - 1000
# samples x 40000 probes, about a mid-size chromosome - under 1 vs 4 OpenBLAS
# threads, read from the kernel's LOC_TIMER phase timers.
#
# Two missingness levels bracket the Gram cache (src/gram_ops.h:63): probes
# with no NA are folded into the Gram once before the loop, so each
# iteration's dsyrk covers only the "changing" probes. CHG_FRACS sets the
# share of probes that carry NAs; each such probe gets NA_PER_COL missing
# rows. Every run does exactly ITERS EM iterations (threshold = 0 can never be
# met), so each loop phase has ITERS samples.
#
# Builds with load_all1() (timer on, diagnostics off), which builds at -O0.
# That slows the package's own loops, not OpenBLAS, a prebuilt DLL, and both
# thread counts build from the same source with the same flags.
#
# OPENBLAS_NUM_THREADS is set before R starts by dev/bench_wide_phases.sh;
# RhpcBLASctl cannot reach Rlapack.dll (dev/probe_blas_leak.R).
#
# Usage, from the repository root (BENCH_OUT names an output directory):
#   Rscript dev/bench_wide_phases.R run <label>
#   Rscript dev/bench_wide_phases.R compare
if (!file.exists("DESCRIPTION")) stop("run from the repository root")

args <- commandArgs(trailingOnly = TRUE)
part <- args[[1L]]
out_dir <- Sys.getenv("BENCH_OUT")
stopifnot(nzchar(out_dir), dir.exists(out_dir))

N_ROW <- 1000L
N_COL <- 40000L
RANK <- 10L
CHG_FRACS <- c(0.1, 1.0)
NA_PER_COL <- 10L
ITERS <- 10L
NCP <- 10L
SOLVER <- "exact"
LOOP_PHASES <- c(
  "restandardize", "svd", "form_gram", "eig", "recover_factors",
  "post_svd", "reconstruct", "objective"
)
# the loop's top-level phases; svd encloses form_gram, eig and recover_factors
TOP_PHASES <- c("restandardize", "svd", "post_svd", "reconstruct", "objective")
ONCE_PHASES <- c("gram_cache_init", "postprocessing", "pca_imp_internal_cpp_total")

if (part == "run") {
  label <- args[[2L]]
  cat(sprintf(
    "%s: OPENBLAS_NUM_THREADS=%s\n",
    label, Sys.getenv("OPENBLAS_NUM_THREADS", "-")
  ))

  source("R/dev-utils.R")
  t_build <- system.time(load_all1())[["elapsed"]]
  cat(sprintf("load_all1() build: %.1f s\n", t_build))

  # low-rank signal plus noise, the same matrix for every condition
  set.seed(20261028)
  x0 <- matrix(rnorm(N_ROW * RANK), N_ROW) %*% matrix(rnorm(RANK * N_COL), RANK) +
    matrix(rnorm(N_ROW * N_COL, sd = 0.5), N_ROW)

  phases <- list()
  walls <- list()
  for (f in CHG_FRACS) {
    tag <- sprintf("chg%03d", round(100 * f))
    set.seed(20261029)
    x <- x0
    chg <- sort(sample.int(N_COL, round(f * N_COL)))
    for (j in chg) x[sample.int(N_ROW, NA_PER_COL), j] <- NA
    n_chg <- sum(colSums(is.na(x)) > 0)
    stopifnot(n_chg == length(chg))

    if (exists("pca_imp_gram", envir = globalenv())) {
      rm("pca_imp_gram", envir = globalenv())
    }
    invisible(gc(verbose = FALSE))
    t0 <- proc.time()
    res <- suppressMessages(suppressWarnings(pca_imp(
      x,
      ncp = NCP, scale = TRUE, threshold = 0,
      maxiter = ITERS, miniter = ITERS, solver = SOLVER, seed = 1L
    )))
    dt <- proc.time() - t0
    wall <- dt[["elapsed"]]
    cpu <- dt[["user.self"]] + dt[["sys.self"]]

    tm <- get("pca_imp_gram", envir = globalenv())
    # rcpptimer puts the phase tags in the row names, not a column
    d <- data.frame(Name = rownames(tm), unclass(tm), row.names = NULL)
    # every loop phase must have one sample per EM iteration, or the means
    # below are over something other than what they claim
    loop <- d[d$Name %in% LOOP_PHASES, ]
    stopifnot(setequal(unique(loop$Name), LOOP_PHASES), all(loop$Count == ITERS))
    d$label <- label
    d$chg <- tag
    phases[[tag]] <- d
    walls[[tag]] <- data.frame(
      label = label, chg = tag, n_chg = n_chg, wall = wall, cpu = cpu,
      build = t_build
    )
    saveRDS(
      as.numeric(res[is.na(x)]),
      file.path(out_dir, sprintf("imp_%s_%s.rds", label, tag))
    )
    cat(sprintf(
      "  %s (%d changing probes): wall %.2f s | cpu %.2f s | cpu/wall %.2f\n",
      tag, n_chg, wall, cpu, cpu / wall
    ))
    rm(x, res)
  }
  write.csv(
    do.call(rbind, phases),
    file.path(out_dir, sprintf("phases_%s.csv", label)),
    row.names = FALSE
  )
  write.csv(
    do.call(rbind, walls),
    file.path(out_dir, sprintf("wall_%s.csv", label)),
    row.names = FALSE
  )
}

if (part == "compare") {
  rd <- function(pat) {
    do.call(rbind, lapply(
      list.files(out_dir, pattern = pat, full.names = TRUE), read.csv
    ))
  }
  ph <- rd("^phases_.*\\.csv$")
  w <- rd("^wall_.*\\.csv$")
  stopifnot(all(c("blas1", "blas4") %in% ph$label))

  cat("== whole pca_imp() call ==\n")
  w$cpu_over_wall <- w$cpu / w$wall
  w <- w[order(w$chg, w$label), ]
  print(format(w, digits = 3), row.names = FALSE)

  ph$ms <- ph$Microseconds / 1000
  ph$se_ms <- ph$SD / sqrt(ph$Count) / 1000
  keep <- c("chg", "Name", "Count", "ms", "se_ms")
  cmp <- merge(
    ph[ph$label == "blas1", keep], ph[ph$label == "blas4", keep],
    by = c("chg", "Name", "Count"), suffixes = c("_1", "_4")
  )
  cmp$speedup_4_over_1 <- cmp$ms_1 / cmp$ms_4
  # difference in units of its own standard error
  cmp$diff_in_se <- (cmp$ms_1 - cmp$ms_4) / sqrt(cmp$se_ms_1^2 + cmp$se_ms_4^2)

  for (tag in sort(unique(cmp$chg))) {
    d <- cmp[cmp$chg == tag, ]
    lp <- d[d$Name %in% LOOP_PHASES, ]
    lp <- lp[order(match(lp$Name, LOOP_PHASES)), ]
    # iteration total from the top-level phases, and each phase's share of
    # it at 1 thread
    top <- lp$Name %in% TOP_PHASES
    it_1 <- sum(lp$ms_1[top])
    it_4 <- sum(lp$ms_4[top])
    lp$share_1 <- lp$ms_1 / it_1
    cat(sprintf(
      "\n== %s: per EM iteration (ms), %d samples each; iteration total %.1f -> %.1f ms ==\n",
      tag, ITERS, it_1, it_4
    ))
    print(format(lp, digits = 3), row.names = FALSE)
    op <- d[d$Name %in% ONCE_PHASES, ]
    cat("-- once per call (ms) --\n")
    print(format(op[, c("Name", "ms_1", "ms_4", "speedup_4_over_1")], digits = 3), row.names = FALSE)

    # worst-case difference over the imputed cells; relative figure over
    # imputed cells with |blas1 value| >= 1e-3
    a <- readRDS(file.path(out_dir, sprintf("imp_blas4_%s.rds", tag)))
    b <- readRDS(file.path(out_dir, sprintf("imp_blas1_%s.rds", tag)))
    stopifnot(length(b) > 0L, length(a) == length(b))
    rel_keep <- abs(b) >= 1e-3
    stopifnot(sum(rel_keep) > 0L)
    cat(sprintf(
      "imputed cells blas4 vs blas1: max abs %.3g, max rel %.3g over %d of %d cells, bitwise %s\n",
      max(abs(a - b)), max(abs(a - b)[rel_keep] / abs(b[rel_keep])),
      sum(rel_keep), length(b), identical(a, b)
    ))
  }
}
