# Inputs for step 6 of Phase 2 (dev/eigen_conversion.md), which took
# src/armaSVD.cpp off Armadillo. The Phase 0 capture runs the kernel through
# pca_imp(), in a build without diagnostics. This calls
# pca_imp_internal_cpp() directly in the diagnostics build
# (-DPCA_IMP_DIAGNOSTICS=1), where each iteration takes the singular values
# of V_prev' * V_new, the one call of the file that can fail, and where the
# kernel returns its histories and the solver's logs:
#
# - wide and tall inputs, with complete columns and with a hole in every
#   column, row weights (some whose sum depends on its order), no scaling,
#   the EM method, a random start;
# - 1 component, where the product is 1 x 1 and goes to dgemv, and up to 40,
#   past the 32 values up to which Armadillo summed the leading singular
#   values itself;
# - the SVD made to fail. Those cases set SLIDEIMP_P2S6_FAIL, which only a
#   build instrumented by dev/probe_p2s6_patch.R reads: on any other build
#   they are plain repeats of the same call.
#
# Each result is reduced to a form both versions of the kernel can give:
# the imputed cells as an n x 3 double matrix (row, col, value), counts and
# histories as double vectors. The two timings of the "auto" probe are
# dropped, and "auto" is not run: it chooses by wall-clock time.
#
# Run from the package root of a build, with a label naming the capture:
#   Rscript dev/eigen_p2s6_diag.R <label>
# It writes dev/ref/<label>/ in dev/ref_capture.R's format, so two builds'
# captures are compared with dev/ref_compare.R. dev/probe_p2s6_run.sh runs
# it on an instrumented clone of each build.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1]])) {
  stop("usage: Rscript dev/eigen_p2s6_diag.R <label>", call. = FALSE)
}
label <- args[[1]]
out_dir <- file.path("dev", "ref", label)
if (dir.exists(out_dir)) {
  stop(out_dir, " exists; pick another label or delete it", call. = FALSE)
}

source("dev/dev-utils.R")
load_all1(timer = FALSE, debug = TRUE)
dir.create(out_dir, recursive = TRUE)

git <- function(...) {
  system2("git", c(...), stdout = TRUE)
}
src_paths <- c("R", "src", "DESCRIPTION", "NAMESPACE")
saveRDS(
  list(
    label = label,
    head = git("rev-parse", "--short", "HEAD"),
    source_commit = git("log", "-1", "--format=%h", "--", src_paths),
    source_dirty = length(git("status", "--porcelain", "--", src_paths)) > 0L,
    r_version = R.version.string,
    la_library = La_library(),
    blas = extSoftVersion()[["BLAS"]],
    rcppeigen = as.character(packageVersion("RcppEigen")),
    rcppeigen_lib = dirname(find.package("RcppEigen")),
    rcpparmadillo = as.character(packageVersion("RcppArmadillo")),
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  file.path(out_dir, "_meta.rds")
)

## ---- inputs ---------------------------------------------------------------

## x1: wide, some columns complete. xw: wide, a hole in every column.
## xt: tall, a hole in every column. xt2: tall, some columns complete.
## xb: 300 x 400, for many components.
set.seed(1)
x1 <- sim_mat(n = 60, p = 300, perc_total_na = 0.1, perc_col_na = 0.5)$input
set.seed(2)
xw <- sim_mat(n = 40, p = 120, perc_total_na = 0.15, perc_col_na = 1)$input
set.seed(3)
xt <- sim_mat(n = 120, p = 40, perc_total_na = 0.15, perc_col_na = 1)$input
set.seed(4)
xt2 <- sim_mat(n = 150, p = 50, perc_total_na = 0.1, perc_col_na = 0.5)$input
set.seed(11)
xb <- sim_mat(n = 300, p = 400, perc_total_na = 0.1, perc_col_na = 0.5)$input

n_complete <- function(x) sum(colSums(is.na(x)) == 0L)
stopifnot(
  ncol(x1) > nrow(x1),
  n_complete(x1) > 0L,
  ncol(xw) > nrow(xw),
  n_complete(xw) == 0L,
  nrow(xt) > ncol(xt),
  n_complete(xt) == 0L,
  nrow(xt2) > ncol(xt2),
  n_complete(xt2) > 0L,
  all(vapply(list(x1, xw, xt, xt2, xb), anyNA, TRUE))
)

## ---- the kernel, called as both versions take it --------------------------

kernel <- function(
  x,
  ncp,
  solver = 0L,
  scale = TRUE,
  regularized = TRUE,
  init = 0L,
  iters = 12L,
  threshold = 0,
  row_w = rep(1, nrow(x)),
  warmup_iters = 3L,
  lobpcg_maxiter = 20L
) {
  res <- pca_imp_internal_cpp(
    obj = x,
    eligible_idx = seq_len(ncol(x)) - 1L,
    ncp = ncp,
    scale = scale,
    regularized = regularized,
    threshold = threshold,
    init = init,
    maxiter = iters,
    miniter = min(iters, 5L),
    row_w = row_w,
    coeff_ridge = 1,
    solver = solver,
    warmup_iters = warmup_iters,
    lobpcg_tol = 1e-9,
    lobpcg_maxiter = lobpcg_maxiter
  )
  cells <- if (is.null(res$imputed_values)) {
    cbind(res$loc, res$value)
  } else {
    res$imputed_values
  }
  storage.mode(cells) <- "double"
  dimnames(cells) <- NULL
  stopifnot(nrow(cells) == sum(is.na(x)), ncol(cells) == 3L)
  hist <- c("eigval_hist", "obj_hist", "subspace_cos_hist")
  stopifnot(hist %in% names(res))
  drop <- c(
    "imputed_values",
    "loc",
    "value",
    "auto_exact_mean_sec",
    "auto_lobpcg_mean_sec"
  )
  out <- res[setdiff(names(res), drop)]
  out$cells <- cells
  out$n_iter <- as.numeric(out$n_iter)
  out$obj_hist <- as.numeric(out$obj_hist)
  out$subspace_cos_hist <- as.numeric(out$subspace_cos_hist)
  ## the first iteration has no subspace before it; every later one has
  stopifnot(
    is.nan(out$subspace_cos_hist[[1]]),
    out$n_iter < 2 || all(is.finite(out$subspace_cos_hist[-1]))
  )
  out
}

## ---- the runner, as in dev/ref_capture.R ----------------------------------

n_cases <- 0L
case <- function(name, expr) {
  expr <- substitute(expr)
  warns <- list()
  set.seed(1)
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
  n_cases <<- n_cases + 1L
  saveRDS(
    list(value = value, warnings = warns),
    file.path(out_dir, paste0(name, ".rds"))
  )
  cat(sprintf(
    "  %-28s%s\n",
    name,
    if (inherits(value, "ref_error")) {
      paste0("  ERROR: ", strsplit(value$message, "\n")[[1]][[1]])
    } else {
      sprintf(
        "  %d iterations, smallest cosine %.6f",
        value$n_iter,
        min(value$subspace_cos_hist, na.rm = TRUE)
      )
    }
  ))
  invisible()
}

inj_id <- 0L
inj <- function(spec, expr) {
  inj_id <<- inj_id + 1L
  Sys.setenv(SLIDEIMP_P2S6_FAIL = sprintf("%s#%d", spec, inj_id))
  on.exit(Sys.unsetenv("SLIDEIMP_P2S6_FAIL"))
  expr
}

## ---- cases ----------------------------------------------------------------

cat("cases:\n")

inputs <- list(x1 = x1, xw = xw, xt = xt, xt2 = xt2)
for (nm in names(inputs)) {
  for (ncp in c(1L, 2L, 5L)) {
    eval(bquote(case(
      .(sprintf("%s_ncp%d_exact", nm, ncp)),
      kernel(inputs[[.(nm)]], ncp = .(ncp))
    )))
    eval(bquote(case(
      .(sprintf("%s_ncp%d_lobpcg", nm, ncp)),
      kernel(inputs[[.(nm)]], ncp = .(ncp), solver = 1L)
    )))
  }
}
case("x1_ncp3_noscale", kernel(x1, ncp = 3L, scale = FALSE))
case("xt_ncp3_noscale", kernel(xt, ncp = 3L, scale = FALSE))
case("x1_ncp3_em", kernel(x1, ncp = 3L, regularized = FALSE))
case(
  "xt2_ncp3_em_lobpcg",
  kernel(xt2, ncp = 3L, regularized = FALSE, solver = 1L)
)
case("x1_ncp3_init2", kernel(x1, ncp = 3L, init = 2L))
case("xt_ncp3_init3", kernel(xt, ncp = 3L, init = 3L))
case(
  "x1_ncp3_roww",
  kernel(x1, ncp = 3L, row_w = seq(0.5, 2, length.out = nrow(x1)))
)
case(
  "xt2_ncp3_roww_lobpcg",
  kernel(
    xt2,
    ncp = 3L,
    solver = 1L,
    row_w = seq(2, 0.25, length.out = nrow(xt2))
  )
)
## Weights whose sum depends on the order it is taken in (case() sets the
## seed before each).
case(
  "x1_ncp3_roww_runif",
  kernel(x1, ncp = 3L, row_w = runif(nrow(x1), 0.1, 3))
)
case(
  "xw_ncp2_roww_runif_lobpcg",
  kernel(xw, ncp = 2L, solver = 1L, row_w = runif(nrow(xw), 0.1, 3))
)
case(
  "xt_ncp2_roww_nmiss",
  kernel(xt, ncp = 2L, row_w = 1 - rowMeans(is.na(xt)))
)
case("x1_ncp3_converge", kernel(x1, ncp = 3L, iters = 200L, threshold = 1e-6))
case("x1_ncp2_one_iter", kernel(x1, ncp = 2L, iters = 1L))
for (ncp in c(31L, 32L, 33L, 40L)) {
  eval(bquote(case(
    .(sprintf("xb_ncp%d_exact", ncp)),
    kernel(xb, ncp = .(ncp), iters = 6L)
  )))
}
case("xb_ncp40_lobpcg", kernel(xb, ncp = 40L, solver = 1L, iters = 8L))
case("xb_ncp33_em", kernel(xb, ncp = 33L, regularized = FALSE, iters = 6L))

## The k-th SVD of the call fails. The first iteration takes none, so the
## k-th is taken in iteration k + 1.
for (k in c(1L, 4L)) {
  eval(bquote(case(
    .(sprintf("inj_svdnan_%d", k)),
    inj(.(sprintf("svdnan:%d", k)), kernel(x1, ncp = 3L))
  )))
  eval(bquote(case(
    .(sprintf("inj_svdinfo_%d", k)),
    inj(.(sprintf("svdinfo:%d", k)), kernel(x1, ncp = 3L))
  )))
}
case("inj_svdnan_ncp1", inj("svdnan:2", kernel(xt, ncp = 1L)))
## Past the last SVD of the call: nothing fails.
case("inj_svdnan_never", inj("svdnan:50", kernel(x1, ncp = 3L)))

stopifnot(n_cases == 49L)
cat(sprintf("%d cases, written to %s\n", n_cases, out_dir))
