# Reference capture: run a fixed set of calls on one build and save every
# result, with its type, dims and attributes, so a later build can be checked
# against it by dev/ref_compare.R. Shared by the si-mi branch
# (dev/si_mi_conversion.md, Phase 0) and the eigen branch
# (dev/eigen_conversion.md, Phase 0), so it covers what both list.
#
# Run from the package root, with a label naming the capture:
#   Rscript dev/ref_capture.R e5b3cf7
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
# it is vendored and gated upstream (dev/vendored.md).
#
# Inputs: sim_mat() output under fixed seeds, edited to hold columns over
# colmax, constant columns and an all-NA column; an input with no holes; and
# the chr21 EPIC data (dev/GSE286313.EPICv1_chr21.qs2, which needs qs2). It
# has no genomic coordinates, so slide_imp() takes the column index as
# location, and group_imp() takes blocks of consecutive columns as groups.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1]])) {
  stop("usage: Rscript dev/ref_capture.R <label>", call. = FALSE)
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

## xa: x2 with an all-NA column, which mean_imp_col() leaves NA and the
## imputers refuse.
xa <- x2
xa[, 9] <- NA

set.seed(3)
x0 <- sim_mat(n = 40, p = 100, perc_total_na = 0)$input
stopifnot(!anyNA(x0))

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
## keeps them, for mean_imp_col().
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
case("pca_x1_exact", pca_imp(x1, ncp = 3, solver = "exact"))
case("pca_x1_lobpcg", pca_imp(x1, ncp = 3, solver = "lobpcg"))
case("pca_x1_auto", pca_imp(x1, ncp = 3))
case("pca_x1_em", pca_imp(x1, ncp = 3, method = "EM", solver = "exact"))
case("pca_x1_noscale", pca_imp(x1, ncp = 3, scale = FALSE, solver = "exact"))
case(
  "pca_x1_roww",
  pca_imp(x1, ncp = 3, row.w = seq_len(nrow(x1)), solver = "exact")
)
case("pca_x1_clamp", pca_imp(x1, ncp = 3, clamp = c(0, 1), solver = "exact"))
case(
  "pca_x1_nbinit",
  pca_imp(x1, ncp = 3, nb.init = 3, seed = 7, solver = "exact")
)
case("pca_x2_post", pca_imp(x2, ncp = 3, solver = "exact"))
case("pca_x2_nopost", pca_imp(x2, ncp = 3, solver = "exact", post_imp = FALSE))
case("pca_x0_noholes", pca_imp(x0, ncp = 2))
case("pca_xa_allna_error", pca_imp(xa, ncp = 3))
case("pca_chr21_exact", pca_imp(chr21, ncp = 5, solver = "exact"))
case("pca_chr21_lobpcg", pca_imp(chr21, ncp = 5, solver = "lobpcg"))
case(
  "pca_chr21_nopost",
  pca_imp(chr21, ncp = 5, solver = "exact", post_imp = FALSE)
)

## mean_imp_col()
case("mean_xa_all", mean_imp_col(xa))
case("mean_xa_subset", mean_imp_col(xa, subset = colnames(xa)[c(1, 6, 9, 20)]))
case("mean_xa_cores2", mean_imp_col(xa, cores = 2))
case("mean_chr21_all", mean_imp_col(chr21_all))

## knn_imp()
case("knn_x1_euc", knn_imp(x1, k = 5))
case("knn_x1_man", knn_imp(x1, k = 5, method = "manhattan"))
case("knn_x1_distpow", knn_imp(x1, k = 5, dist_pow = 1))
case("knn_x1_subset", knn_imp(x1, k = 5, subset = colnames(x1)[1:50]))
case("knn_x1_cores2", knn_imp(x1, k = 5, cores = 2))
case("knn_x2_post", knn_imp(x2, k = 5))
case("knn_x2_nopost", knn_imp(x2, k = 5, post_imp = FALSE))
case("knn_x0_noholes", knn_imp(x0, k = 5))
case("knn_chr21", knn_imp(chr21, k = 10))

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

## Armadillo-facing internals called directly (dev/eigen_conversion.md)
case("int_col_min_max", col_min_max(x2))
case("int_col_vars", col_vars_internal(x2))
case("int_mean_imp_col", mean_imp_col_internal(xa, c(0L, 5L, 8L)))
case("int_check_inf", check_inf(x2))
case("int_check_finite", check_finite(x2))
case("int_check_finite_allna_error", check_finite(xa))
case("int_col_miss", col_miss_internal(x2))
case("int_row_miss", row_miss_internal(x2))
case("int_find_windows", find_windows(seq_len(300), window = 50, overlap = 10))

saveRDS(timings, file.path(out_dir, "_timings.rds"))
cat(sprintf(
  "%d cases, %.1fs, written to %s\n",
  length(timings),
  sum(timings),
  out_dir
))
