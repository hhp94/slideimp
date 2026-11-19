# Inputs for step 4 of Phase 2 (dev/eigen_conversion.md), which moved the QR,
# the triangular solve and the dsyevd of src/lobpcg_warm.h onto si_la. The
# Phase 0 capture asks LOBPCG for at most 6 eigenpairs and never makes one of
# those three calls fail; this runs what it leaves out:
#
# - LOBPCG with 12 to 41 eigenpairs, where the Rayleigh-Ritz problems reach
#   orders 34 to 79 (dev/probe_p2s4_syevd_lwork.R says why that range);
# - samples and columns that are copies of each other, with more components
#   than the data has rank;
# - LOBPCG settings that stall it or stop it early;
# - each of the three calls made to fail on its k-th call, a NaN put into
#   the small Gram or the Gram itself, and a triangle made singular before
#   the solve. These cases set SLIDEIMP_P2S4_FAIL,
#   which only a build instrumented by dev/probe_p2s4_patch.R reads: on any
#   other build they are plain repeats of the same call.
#
# Run from the package root of a build, with a label naming the capture:
#   Rscript dev/eigen_p2s4_edges.R <label>
# It writes dev/ref/<label>/ in dev/ref_capture.R's format, so two builds'
# captures are compared with dev/ref_compare.R. dev/probe_p2s4_run.sh runs
# it on an instrumented clone of each build.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1]])) {
  stop("usage: Rscript dev/eigen_p2s4_edges.R <label>", call. = FALSE)
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
    fresh = nzchar(Sys.getenv("SLIDEIMP_P2S4_FRESH")),
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  file.path(out_dir, "_meta.rds")
)

## ---- inputs ---------------------------------------------------------------

set.seed(1)
x1 <- sim_mat(n = 60, p = 300, perc_total_na = 0.1, perc_col_na = 0.5)$input
set.seed(11)
xb <- sim_mat(n = 300, p = 400, perc_total_na = 0.1, perc_col_na = 0.5)$input

chr21_all <- t(qs2::qs_read("dev/GSE286313.EPICv1_chr21.qs2"))
chr21 <- chr21_all[, colSums(!is.na(chr21_all)) > 0L]

## xd: 30 samples of which 6 are distinct, holes included, so the completed
## data has rank 6 at most. xdc: 40 x 12 with columns 4 to 12 copies of
## columns 1 to 3.
set.seed(21)
base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
xd <- base6[rep(seq_len(6), times = 5), ]
rownames(xd) <- sprintf("s%02d", seq_len(nrow(xd)))
set.seed(22)
base3 <- matrix(rnorm(40 * 3), 40, 3)
xdc <- base3[, rep(1:3, times = 4)]
xdc[cbind(c(2, 9, 17, 30), c(1, 5, 8, 12))] <- NA
dimnames(xdc) <- list(sprintf("s%02d", 1:40), sprintf("c%02d", 1:12))
stopifnot(
  anyNA(xd),
  nrow(unique(xd)) == 6L,
  qr(base3)$rank == 3L
)

## ---- the runner, as in dev/ref_capture.R ----------------------------------

n_cases <- 0L
case <- function(name, expr) {
  expr <- substitute(expr)
  warns <- list()
  Sys.setenv(SLIDEIMP_P2S4_CASE = name)
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
    "  %-32s%s%s\n",
    name,
    if (inherits(value, "ref_error")) {
      paste0("  ERROR: ", strsplit(value$message, "\n")[[1]][[1]])
    } else {
      sprintf("  NA left %d", sum(is.na(value)))
    },
    if (length(warns)) sprintf("  %d warning(s)", length(warns)) else ""
  ))
  invisible()
}

## The k-th call at `site` fails, in an instrumented build. The tag makes
## every value of the variable new, which is what resets the count.
inj_id <- 0L
inj <- function(spec, expr) {
  inj_id <<- inj_id + 1L
  Sys.setenv(SLIDEIMP_P2S4_FAIL = sprintf("%s#%d", spec, inj_id))
  on.exit(Sys.unsetenv("SLIDEIMP_P2S4_FAIL"))
  expr
}

## ---- cases ----------------------------------------------------------------

cat("cases:\n")

## Many eigenpairs. The regularized method asks for ncp + 1, EM for ncp.
for (ncp in c(11L, 16L, 25L, 40L)) {
  eval(bquote(case(
    .(sprintf("big_xb_ncp%02d", ncp)),
    pca_imp(xb, ncp = .(ncp), solver = "lobpcg")
  )))
}
case(
  "big_xb_ncp17_em",
  pca_imp(xb, ncp = 17, method = "EM", solver = "lobpcg")
)
case("big_xb_ncp25_exact", pca_imp(xb, ncp = 25, solver = "exact"))
for (ncp in c(11L, 20L, 30L)) {
  eval(bquote(case(
    .(sprintf("big_chr21_ncp%02d", ncp)),
    pca_imp(chr21, ncp = .(ncp), solver = "lobpcg")
  )))
}

## Copies: more components than rank. These fits converge within the default
## warm-up of 10 exact solves, so LOBPCG would never run: it starts after one
## exact solve here, and the fit is held to 15 iterations.
early <- lobpcg_control(warmup_iters = 1)
for (ncp in c(5L, 7L, 12L)) {
  eval(bquote(case(
    .(sprintf("dup_xd_ncp%02d", ncp)),
    pca_imp(
      xd,
      ncp = .(ncp),
      solver = "lobpcg",
      lobpcg_control = early,
      miniter = 15
    )
  )))
}
case(
  "dup_xd_ncp07_em",
  pca_imp(
    xd,
    ncp = 7,
    method = "EM",
    solver = "lobpcg",
    lobpcg_control = early,
    miniter = 15
  )
)
case(
  "dup_xdc_ncp05",
  pca_imp(xdc, ncp = 5, solver = "lobpcg", lobpcg_control = early, miniter = 15)
)
case(
  "dup_xdc_ncp08",
  pca_imp(xdc, ncp = 8, solver = "lobpcg", lobpcg_control = early, miniter = 15)
)

## Settings that stall LOBPCG, stop it early, or start it at once.
case(
  "ctl_x1_tight",
  pca_imp(
    x1,
    ncp = 3,
    solver = "lobpcg",
    lobpcg_control = lobpcg_control(tol = 1e-14, maxiter = 40)
  )
)
case(
  "ctl_x1_maxiter2",
  pca_imp(
    x1,
    ncp = 3,
    solver = "lobpcg",
    lobpcg_control = lobpcg_control(maxiter = 2)
  )
)
case(
  "ctl_x1_warmup1",
  pca_imp(
    x1,
    ncp = 3,
    solver = "lobpcg",
    lobpcg_control = lobpcg_control(warmup_iters = 1)
  )
)
case(
  "ctl_x1_noscale",
  pca_imp(x1, ncp = 3, scale = FALSE, solver = "lobpcg")
)
case("ctl_x1_ncp1", pca_imp(x1, ncp = 1, method = "EM", solver = "lobpcg"))

## The k-th call of one kind fails. With the default warm-up of 10 exact
## solves, QRs 1 to 18 seed the state after exact solves 2 to 10, 19 and 20
## clean P at the entry of the first lobpcg_solve(), and the ones after are
## its residual and P blocks, then the refresh after it.
for (k in c(1:2, 17:30)) {
  eval(bquote(case(
    .(sprintf("inj_qr_%02d", k)),
    inj(.(sprintf("qr:%d", k)), pca_imp(x1, ncp = 3, solver = "lobpcg"))
  )))
}
for (k in 1:8) {
  eval(bquote(case(
    .(sprintf("inj_syevd_%02d", k)),
    inj(.(sprintf("syevd:%d", k)), pca_imp(x1, ncp = 3, solver = "lobpcg"))
  )))
}
for (k in c(1L, 3L, 6L)) {
  eval(bquote(case(
    .(sprintf("inj_syevdnan_%02d", k)),
    inj(.(sprintf("syevdnan:%d", k)), pca_imp(x1, ncp = 3, solver = "lobpcg"))
  )))
}
for (k in 1:2) {
  eval(bquote(case(
    .(sprintf("inj_solve_%02d", k)),
    inj(.(sprintf("solve:%d", k)), pca_imp(x1, ncp = 3, solver = "lobpcg"))
  )))
}
for (k in 1:2) {
  eval(bquote(case(
    .(sprintf("inj_gramnan_%02d", k)),
    inj(.(sprintf("gramnan:%d", k)), pca_imp(x1, ncp = 3, solver = "lobpcg"))
  )))
}
case(
  "inj_xb_qr_05",
  inj("qr:5", pca_imp(xb, ncp = 11, solver = "lobpcg"))
)
case(
  "inj_xb_syevd_04",
  inj("syevd:4", pca_imp(xb, ncp = 11, solver = "lobpcg"))
)
## A NaN in a small Gram of order 34 and more, where dsyevd takes its divide
## and conquer path (x1's are of order 8 to 12).
for (k in 1:2) {
  eval(bquote(case(
    .(sprintf("inj_xb16_syevdnan_%02d", k)),
    inj(.(sprintf("syevdnan:%d", k)), pca_imp(xb, ncp = 16, solver = "lobpcg"))
  )))
}
case(
  "inj_xb16_gramnan_01",
  inj("gramnan:1", pca_imp(xb, ncp = 16, solver = "lobpcg"))
)
## A triangle singular to working precision in the k-th solve.
for (k in c(1L, 2L, 5L)) {
  eval(bquote(case(
    .(sprintf("inj_singular_%02d", k)),
    inj(.(sprintf("singular:%d", k)), pca_imp(x1, ncp = 3, solver = "lobpcg"))
  )))
}
case(
  "inj_xb_singular_03",
  inj("singular:3", pca_imp(xb, ncp = 11, solver = "lobpcg"))
)

stopifnot(n_cases == 60L)
cat(sprintf("%d cases, written to %s\n", n_cases, out_dir))
