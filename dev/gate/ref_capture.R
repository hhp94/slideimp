# Reference capture: run a fixed set of calls on one build and save every
# result, with its type, dims and attributes, so a later build can be checked
# against it by dev/gate/ref_compare.R. Shared by the two plans on the si-mi
# branch, dev/si_mi_conversion.md at d5a34e5 and dev/eigen_conversion.md at
# 4f2b669 (Phase 0 in each), so it covers what both list.
#
# Run from the package root, with a label naming the capture:
#   Rscript dev/gate/ref_capture.R e5b3cf7
# It writes dev/ref/<label>/, one .rds per case plus _meta.rds. dev/ref/ is
# ignored (dev/ is opt-in) and regenerable: check out the reference commit and
# rerun. _meta.rds records the commit that last touched the package source,
# so a capture says which source it ran on whatever HEAD is.
#
# Each case runs under its own set.seed(1), so adding, removing or reordering
# cases leaves the others unchanged. A case that errors is stored as its
# condition class and message, and every case stores the warnings it raised,
# class and message, so the error surface is captured with the values.
# Messages are not stored. The sampler (mipca() and friends) is left out:
# its draws are compared on their own (dev/to-do.md, "Release condition:
# pass on the RcppEigen candidate").
#
# Inputs: sim_mat() output under fixed seeds, edited to hold columns over
# colmax, constant columns, an all-NA column or one Inf; an input with no
# holes; inputs sized or holed so pca_imp() takes an eigensystem too large
# for R to demote solver = "auto" to exact, and the no-Gram-cache branches,
# wide and tall; and the chr21 EPIC data (dev/GSE286313.EPICv1_chr21.qs2,
# which needs qs2). It has no genomic coordinates, so slide_imp() takes the
# column index as location, and group_imp() takes blocks of consecutive
# columns as groups.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1]])) {
  stop("usage: Rscript dev/gate/ref_capture.R <label>", call. = FALSE)
}
label <- args[[1]]
out_dir <- file.path("dev", "ref", label)
if (dir.exists(out_dir)) {
  stop(out_dir, " exists; pick another label or delete it", call. = FALSE)
}

source("dev/dev-utils.R")
load_all1(timer = FALSE)
dir.create(out_dir, recursive = TRUE)

git <- function(...) {
  system2("git", c(...), stdout = TRUE)
}
src_paths <- c("R", "src", "DESCRIPTION", "NAMESPACE")
meta <- list(
  label = label,
  head = git("rev-parse", "--short", "HEAD"),
  source_commit = git("log", "-1", "--format=%h", "--", src_paths),
  source_dirty = length(git("status", "--porcelain", "--", src_paths)) > 0L,
  r_version = R.version.string,
  la_library = La_library(),
  blas = extSoftVersion()[["BLAS"]],
  ## The build links to whichever RcppEigen and RcppArmadillo come first on
  ## the library path, which R_LIBS set before R starts can change
  ## (dev/eigen_conversion.md at 4f2b669, "Both RcppEigen versions").
  rcppeigen = as.character(packageVersion("RcppEigen")),
  rcppeigen_lib = dirname(find.package("RcppEigen")),
  rcpparmadillo = as.character(packageVersion("RcppArmadillo")),
  time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
saveRDS(meta, file.path(out_dir, "_meta.rds"))
cat(sprintf(
  "capture %s: source %s%s, HEAD %s\n",
  label,
  meta$source_commit,
  if (meta$source_dirty) " (DIRTY)" else "",
  meta$head
))
cat(sprintf(
  "RcppEigen %s from %s; RcppArmadillo %s\n",
  meta$rcppeigen,
  meta$rcppeigen_lib,
  meta$rcpparmadillo
))

## ---- inputs ---------------------------------------------------------------

set.seed(1)
sim1 <- sim_mat(n = 60, p = 300, perc_total_na = 0.1, perc_col_na = 0.5)
x1 <- sim1$input

## x2: x1 with the awkward columns. Columns 1 to 5 are 95% NA (over the
## default colmax of 0.9); 6 is constant with holes; 7 is constant with no
## holes; 8 is exactly 90% NA, on the colmax boundary.
set.seed(2)
x2 <- x1
for (j in 1:5) {
  x2[sample.int(nrow(x2), 57L), j] <- NA
}
x2[, 6] <- 0.5
x2[sample.int(nrow(x2), 6L), 6] <- NA
x2[, 7] <- 0.25
x2[, 8] <- x1[, 8]
x2[sample.int(nrow(x2), 54L), 8] <- NA

## xa: x2 with an all-NA column, which the imputers refuse. mean_imp_col(),
## which left it NA, had four cases here until it was dropped (Phase 9 of
## dev/si_mi_conversion.md at d5a34e5).
xa <- x2
xa[, 9] <- NA

## x1z: x1 with column 1 a copy of column 2, holed on other rows, so the two
## are at distance zero and the K-NN weights take their zero-distance side.
x1z <- x1
z_obs <- which(!is.na(x1[, 2]))
stopifnot(length(z_obs) >= 20L)
x1z[, 1] <- x1[, 2]
x1z[z_obs[seq(1L, length(z_obs), by = 4L)], 1] <- NA
stopifnot(
  anyNA(x1z[, 1]),
  sum(!is.na(x1z[, 1]) & !is.na(x1z[, 2])) >= 10L,
  identical(x1z[!is.na(x1z[, 1]), 1], x1z[!is.na(x1z[, 1]), 2])
)

set.seed(3)
x0 <- sim_mat(n = 40, p = 100, perc_total_na = 0)$input
stopifnot(!anyNA(x0))

## x2i: x2 with one Inf, for the messages that name the offending cell. It
## sits off the first row and column, so an off-by-one in the 1-based row
## and column those messages print shows. It sits in a column with holes for
## mean_imp_col_internal(), deleted in Phase 1 of dev/eigen_conversion.md at
## 4f2b669; the placement stays so the input matches the reference capture's.
x2i <- x2
inf_col <- which(colSums(is.na(x2)) > 0L & seq_len(ncol(x2)) >= 10L)[[1]]
inf_row <- which(!is.na(x2[, inf_col]) & seq_len(nrow(x2)) >= 3L)[[1]]
x2i[inf_row, inf_col] <- Inf

## Inputs that take the PCA paths x1 does not (dev/eigen_conversion.md at
## 4f2b669, Phase 0 review). xb: min(nrow, eligible columns) >= 250, the
## eigensystem size from which R stops demoting solver = "auto" to exact; it is
## captured under each forced solver only, since the kernel's own choice
## between them is timed and so differs between machines. xw and xt: every
## column has a hole, so the kernel has no Gram cache, wide and tall. xt2:
## tall, with complete columns, for LOBPCG.
set.seed(11)
xb <- sim_mat(n = 300, p = 400, perc_total_na = 0.1, perc_col_na = 0.5)$input
set.seed(12)
xw <- sim_mat(n = 60, p = 300, perc_total_na = 0.1, perc_col_na = 1)$input
set.seed(13)
xt <- sim_mat(n = 100, p = 50, perc_total_na = 0.1, perc_col_na = 1)$input
set.seed(14)
xt2 <- sim_mat(n = 100, p = 50, perc_total_na = 0.1, perc_col_na = 0.5)$input
stopifnot(all(colSums(is.na(xw)) > 0L), all(colSums(is.na(xt)) > 0L))

## xi: a window that K-NN cannot impute (one complete row), for
## on_infeasible = "mean"; the pattern the slide_imp() tests use.
set.seed(4)
xi <- sim_mat(20, 100, perc_total_na = 0.2)$input
xi[1:19, 1:10] <- NA
xi[20, 1:10] <- rnorm(10)

grp1 <- data.frame(
  group = sim1$col_group$group,
  feature = sim1$col_group$feature
)
## xi's groups: its first 10 columns alone, so that group is infeasible.
grp_i <- data.frame(
  group = ifelse(seq_len(ncol(xi)) <= 10L, "g_bad", "g_ok"),
  feature = colnames(xi)
)

## chr21 without its all-NA probes, which the imputers refuse; chr21_all
## keeps them.
chr21_all <- t(qs2::qs_read("dev/GSE286313.EPICv1_chr21.qs2"))
chr21 <- chr21_all[, colSums(!is.na(chr21_all)) > 0L]
grp_chr21 <- data.frame(
  group = sprintf("b%03d", (seq_len(ncol(chr21)) - 1L) %/% 500L),
  feature = colnames(chr21)
)
cat(sprintf(
  "chr21: %d x %d, %.2f%% NA (%d all-NA probes dropped), %d over colmax\n",
  nrow(chr21),
  ncol(chr21),
  100 * mean(is.na(chr21)),
  ncol(chr21_all) - ncol(chr21),
  sum(colMeans(is.na(chr21)) > 0.9)
))

## ---- the runner -----------------------------------------------------------

timings <- numeric()
case <- function(name, expr) {
  expr <- substitute(expr)
  warns <- list()
  set.seed(1)
  t0 <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      suppressMessages(eval(expr, envir = globalenv())),
      warning = function(w) {
        warns[[length(warns) + 1L]] <<- list(
          class = class(w),
          message = conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      structure(
        list(class = class(e), message = conditionMessage(e)),
        class = "ref_error"
      )
    }
  )
  dt <- proc.time()[["elapsed"]] - t0
  timings[[name]] <<- dt
  saveRDS(
    list(value = value, warnings = warns),
    file.path(out_dir, paste0(name, ".rds"))
  )
  cat(sprintf(
    "  %-40s %7.2fs%s%s\n",
    name,
    dt,
    if (inherits(value, "ref_error")) {
      paste0("  ERROR: ", strsplit(value$message, "\n")[[1]][[1]])
    } else {
      ""
    },
    if (length(warns)) sprintf("  %d warning(s)", length(warns)) else ""
  ))
  invisible()
}

## ---- cases ----------------------------------------------------------------

cat("cases:\n")

## Inputs themselves: sim_mat() and the edits above.
case("sim_mat_x1", sim1)
case("sim_mat_no_beta", sim_mat(n = 30, p = 40, beta = FALSE))
case("input_x2", x2)

## pca_imp()
##
## A fit left to converge stops at the first iteration from miniter whose
## criterion is under threshold. Near that point the criterion barely moves
## per iteration, so rounding, which differs by BLAS, can decide between two
## neighbouring iterations and so move the imputed values. Every case below
## that would run to convergence runs a fixed count instead, maxiter = miniter
## = n, with n per input and method chosen far enough past convergence that
## the criterion at n is under threshold / 100. So converged is TRUE and no
## non-convergence warning is raised whatever the BLAS. The cases that are
## stopped early on purpose, the error cases and the one with no holes keep
## the defaults.
n_x1 <- 100L
n_x1_em <- 300L
n_x2 <- 300L
n_chr21 <- 300L
n_xb <- 300L
n_xw <- 300L
n_xt <- 100L
n_xt2 <- 200L
case(
  "pca_x1_exact",
  pca_imp(x1, ncp = 3, solver = "exact", maxiter = n_x1, miniter = n_x1)
)
case(
  "pca_x1_lobpcg",
  pca_imp(x1, ncp = 3, solver = "lobpcg", maxiter = n_x1, miniter = n_x1)
)
## The default solver: R demotes it to exact, the eigensystem being small.
case("pca_x1_auto", pca_imp(x1, ncp = 3, maxiter = n_x1, miniter = n_x1))
case(
  "pca_x1_em",
  pca_imp(
    x1,
    ncp = 3,
    method = "EM",
    solver = "exact",
    maxiter = n_x1_em,
    miniter = n_x1_em
  )
)
case(
  "pca_x1_noscale",
  pca_imp(
    x1,
    ncp = 3,
    scale = FALSE,
    solver = "exact",
    maxiter = n_x1,
    miniter = n_x1
  )
)
case(
  "pca_x1_roww",
  pca_imp(
    x1,
    ncp = 3,
    row.w = seq_len(nrow(x1)),
    solver = "exact",
    maxiter = n_x1,
    miniter = n_x1
  )
)
case(
  "pca_x1_clamp",
  pca_imp(
    x1,
    ncp = 3,
    clamp = c(0, 1),
    solver = "exact",
    maxiter = n_x1,
    miniter = n_x1
  )
)
case(
  "pca_x1_nbinit",
  pca_imp(
    x1,
    ncp = 3,
    nb.init = 3,
    seed = 7,
    solver = "exact",
    maxiter = n_x1,
    miniter = n_x1
  )
)
case(
  "pca_x2_post",
  pca_imp(x2, ncp = 3, solver = "exact", maxiter = n_x2, miniter = n_x2)
)
case(
  "pca_x2_nopost",
  pca_imp(
    x2,
    ncp = 3,
    solver = "exact",
    post_imp = FALSE,
    maxiter = n_x2,
    miniter = n_x2
  )
)
case("pca_x0_noholes", pca_imp(x0, ncp = 2))
case("pca_xa_allna_error", pca_imp(xa, ncp = 3))
case(
  "pca_chr21_exact",
  pca_imp(
    chr21,
    ncp = 5,
    solver = "exact",
    maxiter = n_chr21,
    miniter = n_chr21
  )
)
case(
  "pca_chr21_lobpcg",
  pca_imp(
    chr21,
    ncp = 5,
    solver = "lobpcg",
    maxiter = n_chr21,
    miniter = n_chr21
  )
)
case(
  "pca_chr21_nopost",
  pca_imp(
    chr21,
    ncp = 5,
    solver = "exact",
    post_imp = FALSE,
    maxiter = n_chr21,
    miniter = n_chr21
  )
)
## xb under each forced solver; its default solver is left out (see xb).
case(
  "pca_xb_exact",
  pca_imp(xb, ncp = 3, solver = "exact", maxiter = n_xb, miniter = n_xb)
)
case(
  "pca_xb_lobpcg",
  pca_imp(xb, ncp = 3, solver = "lobpcg", maxiter = n_xb, miniter = n_xb)
)
## LOBPCG failing every iteration and falling back to exact.
case(
  "pca_x1_lobpcg_fallback",
  pca_imp(
    x1,
    ncp = 3,
    solver = "lobpcg",
    lobpcg_control = lobpcg_control(maxiter = 1, tol = 1e-14),
    maxiter = n_x1,
    miniter = n_x1
  )
)
## EM asks the eigensolver for ncp vectors, regularized for ncp + 1.
case(
  "pca_x1_em_lobpcg",
  pca_imp(
    x1,
    ncp = 3,
    method = "EM",
    solver = "lobpcg",
    maxiter = n_x1_em,
    miniter = n_x1_em
  )
)
## Stopped before converging: the warning that prints criterion_final.
case(
  "pca_x1_unconverged",
  pca_imp(x1, ncp = 3, solver = "exact", maxiter = 3, miniter = 1)
)
## No Gram cache, wide and tall; unequal row weights on the wide one.
case(
  "pca_xw_nocache",
  pca_imp(xw, ncp = 3, solver = "exact", maxiter = n_xw, miniter = n_xw)
)
case(
  "pca_xw_nocache_roww",
  pca_imp(
    xw,
    ncp = 3,
    row.w = seq_len(nrow(xw)),
    solver = "exact",
    maxiter = n_xw,
    miniter = n_xw
  )
)
case(
  "pca_xt_nocache",
  pca_imp(xt, ncp = 3, solver = "exact", maxiter = n_xt, miniter = n_xt)
)
case(
  "pca_xt2_lobpcg",
  pca_imp(xt2, ncp = 3, solver = "lobpcg", maxiter = n_xt2, miniter = n_xt2)
)
## LOBPCG twins of the cases above that force the exact solver, so each
## kernel path is captured under both. The warm-up is 1 iteration, not the
## default 10, so that the short fits reach LOBPCG too. pca_lo() takes the
## fixed count n, as maxiter = miniter = n.
lo <- lobpcg_control(warmup_iters = 1L)
pca_lo <- function(n, ...) {
  pca_imp(..., solver = "lobpcg", lobpcg_control = lo, maxiter = n, miniter = n)
}
case("pca_x1_noscale_lobpcg", pca_lo(n_x1, x1, ncp = 3, scale = FALSE))
case(
  "pca_x1_roww_lobpcg",
  pca_lo(n_x1, x1, ncp = 3, row.w = seq_len(nrow(x1)))
)
case("pca_x1_clamp_lobpcg", pca_lo(n_x1, x1, ncp = 3, clamp = c(0, 1)))
case(
  "pca_x1_nbinit_lobpcg",
  pca_lo(n_x1, x1, ncp = 3, nb.init = 3, seed = 7)
)
case("pca_x2_post_lobpcg", pca_lo(n_x2, x2, ncp = 3))
case("pca_x2_nopost_lobpcg", pca_lo(n_x2, x2, ncp = 3, post_imp = FALSE))
case(
  "pca_chr21_nopost_lobpcg",
  pca_lo(n_chr21, chr21, ncp = 5, post_imp = FALSE)
)
case(
  "pca_x1_unconverged_lobpcg",
  pca_imp(
    x1,
    ncp = 3,
    solver = "lobpcg",
    lobpcg_control = lo,
    maxiter = 3,
    miniter = 1
  )
)
case("pca_xw_nocache_lobpcg", pca_lo(n_xw, xw, ncp = 3))
case(
  "pca_xw_nocache_roww_lobpcg",
  pca_lo(n_xw, xw, ncp = 3, row.w = seq_len(nrow(xw)))
)
case("pca_xt_nocache_lobpcg", pca_lo(n_xt, xt, ncp = 3))
case("pca_x2i_inf_error", pca_imp(x2i, ncp = 3))

## knn_imp()
case("knn_x1_euc", knn_imp(x1, k = 5))
case("knn_x1_man", knn_imp(x1, k = 5, method = "manhattan"))
case("knn_x1_distpow", knn_imp(x1, k = 5, dist_pow = 1))
case("knn_x1z_zero_distpow", knn_imp(x1z, k = 5, dist_pow = 1))
case("knn_x1z_zero_distpow0", knn_imp(x1z, k = 5))
case("knn_x1_subset", knn_imp(x1, k = 5, subset = colnames(x1)[1:50]))
case("knn_x1_cores2", knn_imp(x1, k = 5, cores = 2))
case("knn_x2_post", knn_imp(x2, k = 5))
case("knn_x2_nopost", knn_imp(x2, k = 5, post_imp = FALSE))
case("knn_x0_noholes", knn_imp(x0, k = 5))
case("knn_chr21", knn_imp(chr21, k = 10))
case("knn_x2i_inf_error", knn_imp(x2i, k = 5))

## slide_imp()
case(
  "slide_x1_knn",
  slide_imp(
    x1,
    location = seq_len(ncol(x1)),
    window_size = 50,
    overlap_size = 10,
    min_window_n = 20,
    k = 5,
    .progress = FALSE
  )
)
case(
  "slide_x1_pca",
  slide_imp(
    x1,
    location = seq_len(ncol(x1)),
    window_size = 50,
    overlap_size = 10,
    min_window_n = 20,
    ncp = 2,
    .progress = FALSE
  )
)
case(
  "slide_x1_flank",
  slide_imp(
    x1,
    location = seq_len(ncol(x1)),
    window_size = 20,
    flank = TRUE,
    min_window_n = 10,
    subset = colnames(x1)[seq(5, 300, by = 25)],
    k = 5,
    .progress = FALSE
  )
)
case(
  "slide_xi_knn_mean",
  slide_imp(
    xi,
    location = seq_len(ncol(xi)),
    window_size = 10,
    min_window_n = 5,
    k = 3,
    on_infeasible = "mean",
    .progress = FALSE
  )
)
case(
  "slide_xi_pca_mean",
  slide_imp(
    xi,
    location = seq_len(ncol(xi)),
    window_size = 10,
    min_window_n = 5,
    ncp = 2,
    on_infeasible = "mean",
    .progress = FALSE
  )
)
case(
  "slide_chr21_knn",
  slide_imp(
    chr21,
    location = seq_len(ncol(chr21)),
    window_size = 1000,
    overlap_size = 100,
    min_window_n = 50,
    k = 10,
    .progress = FALSE
  )
)

## group_imp()
case("group_x1_knn", group_imp(x1, group = grp1, k = 5, .progress = FALSE))
case("group_x1_pca", group_imp(x1, group = grp1, ncp = 2, .progress = FALSE))
case(
  "group_xi_knn_mean",
  group_imp(xi, group = grp_i, k = 3, on_infeasible = "mean", .progress = FALSE)
)
case(
  "group_chr21_pca",
  group_imp(chr21, group = grp_chr21, ncp = 3, .progress = FALSE)
)

## Small exports and the tune_imp() NA placement
case("col_vars_x2", col_vars(x2))
case("col_vars_chr21", col_vars(chr21))
case("mat_miss_col", mat_miss(x2))
case("mat_miss_row_prop", mat_miss(x2, col = FALSE, prop = TRUE))
case("na_loc_x1", sample_na_loc(x1, n_cols = 20, n_reps = 3))
case("na_loc_x1_num", sample_na_loc(x1, num_na = 100, n_reps = 2))
case(
  "tune_x1_knn",
  tune_imp(
    x1,
    parameters = data.frame(k = c(3, 5)),
    .f = "knn_imp",
    n_reps = 2,
    num_na = 100,
    .progress = FALSE
  )
)
case(
  "tune_x1_pca",
  tune_imp(
    x1,
    parameters = data.frame(ncp = c(1, 3)),
    .f = "pca_imp",
    n_reps = 2,
    num_na = 100,
    .progress = FALSE
  )
)

## Internals called directly. The scan replaced col_vars_internal(),
## check_inf(), check_finite(), col_miss_internal() and row_miss_internal(),
## whose cases a capture taken before it holds instead.
case("int_scan", col_scan(x2))
case("int_scan_var_rows", col_scan(x2, var = TRUE, rows = TRUE))
case("int_scan_inf", col_scan(x2i, var = TRUE, rows = TRUE))
case("int_scan_allna", col_scan(xa, var = TRUE, rows = TRUE))
case("int_scan_cols", col_scan(x2, cols = c(9L, 1L, 4L), var = TRUE))
case("int_find_windows", find_windows(seq_len(300), window = 50, overlap = 10))

saveRDS(timings, file.path(out_dir, "_timings.rds"))
cat(sprintf(
  "%d cases, %.1fs, written to %s\n",
  length(timings),
  sum(timings),
  out_dir
))
