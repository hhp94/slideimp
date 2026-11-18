# Edge inputs for the exports Phase 1 of dev/eigen_conversion.md converts.
# Step 1 (src/mat_stats.cpp): empty dims, integer and logical input, a
# vector that is not a matrix, odd col_idx and cores values, and a zero held
# as both 0 and -0 at a column's extreme. The Phase 0 capture runs these
# exports on ordinary input only; this checks that a conversion accepts,
# coerces and refuses what the Armadillo signatures did. A later step adds
# its own cases here.
#
# Run from the package root of a build, with a label naming the capture:
#   Rscript dev/eigen_p1_edges.R <label>
# It writes dev/ref/<label>/ in dev/ref_capture.R's format, so two builds'
# captures are compared with dev/ref_compare.R; compare at tolerance 0. Run
# it on the start commit in a clone, as the reference capture was, and on
# the converted build. Every call is positional, since step 1 renamed the
# first argument from `mat` to `obj`. Character input is left out: under a
# load_all1() build it can kill the R process instead of raising.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !nzchar(args[[1]])) {
  stop("usage: Rscript dev/eigen_p1_edges.R <label>", call. = FALSE)
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
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  file.path(out_dir, "_meta.rds")
)
cat("formals:", names(formals(col_vars_internal)), "\n")

## The runner, as in dev/ref_capture.R.
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
  saveRDS(
    list(value = value, warnings = warns),
    file.path(out_dir, paste0(name, ".rds"))
  )
  cat(sprintf(
    "  %-34s %s%s\n",
    name,
    if (inherits(value, "ref_error")) {
      paste0("ERROR [", value$class[[1]], "] ", value$message)
    } else {
      paste(class(value), collapse = "/")
    },
    if (length(warns)) sprintf("  %d warning(s)", length(warns)) else ""
  ))
  invisible()
}

## ---- inputs ---------------------------------------------------------------

m0r <- matrix(numeric(0), 0, 2)
m0c <- matrix(numeric(0), 3, 0)
m00 <- matrix(numeric(0), 0, 0)
one <- matrix(c(1, 2), 1, 2)
mi <- matrix(c(1L, NA, 3L, 4L, 5L, NA), 3, 2)
ml <- matrix(c(TRUE, NA, FALSE, TRUE), 2, 2)
vec <- c(1, 2, 3)
## zeros of both signs at the minimum, in each order and row parity
zmin <- cbind(c(0, -0, 1), c(-0, 0, 1), c(1, 0, -0), c(1, -0, 0))
## and at the maximum
zmax <- -zmin
mnan <- cbind(
  c(NaN, 1, 2),
  c(NA, NA, NA),
  c(NaN, NaN, NaN),
  c(Inf, 1, 2),
  c(1, -Inf, NA)
)
mm <- cbind(c(1, NA, 3), c(NA, NA, NA), c(1, 2, 3), c(Inf, NA, 1), c(2, NaN, 4))
inf_late <- cbind(c(1, Inf), c(NA, NA))

## ---- cases ----------------------------------------------------------------

cat("cases:\n")

case("cmm_0rows", col_min_max(m0r))
case("cmm_0cols", col_min_max(m0c))
case("cmm_0x0", col_min_max(m00))
case("cmm_int", col_min_max(mi))
case("cmm_lgl", col_min_max(ml))
case("cmm_vec", col_min_max(vec))
case("cmm_zero_min", col_min_max(zmin))
case("cmm_zero_max", col_min_max(zmax))
case("cmm_nan", col_min_max(mnan))

case("cv_one_row", col_vars_internal(one))
case("cv_0rows", col_vars_internal(m0r))
case("cv_0cols", col_vars_internal(m0c))
case("cv_nan", col_vars_internal(mnan))
case("cv_int", col_vars_internal(mi))
case("cv_lgl", col_vars_internal(ml))
case("cv_vec", col_vars_internal(vec))
case("cv_cores_0", col_vars_internal(mnan, 0L))
case("cv_cores_neg", col_vars_internal(mnan, -3L))
case("cv_cores_na", col_vars_internal(mnan, NA_integer_))
case("cv_cores_frac", col_vars_internal(mnan, 2.5))
case("cv_cores_len2", col_vars_internal(mnan, c(1L, 2L)))
case("cv_cores_2", col_vars_internal(mnan, 2L))

case("mic_idx_empty", mean_imp_col_internal(mm, integer(0)))
case("mic_idx_mean", mean_imp_col_internal(mm, c(0L, 4L)))
case("mic_idx_allna", mean_imp_col_internal(mm, 1L))
case("mic_idx_complete", mean_imp_col_internal(mm, 2L))
case("mic_idx_dup", mean_imp_col_internal(mm, c(0L, 0L)))
case("mic_idx_double", mean_imp_col_internal(mm, c(0, 4)))
case("mic_idx_frac", mean_imp_col_internal(mm, 4.7))
case("mic_idx_neg_frac", mean_imp_col_internal(mm, -0.5))
case("mic_idx_neg1_frac", mean_imp_col_internal(mm, -1.5))
case("mic_idx_neg_int", mean_imp_col_internal(mm, -1L))
case("mic_idx_na_int", mean_imp_col_internal(mm, NA_integer_))
case("mic_idx_na_real", mean_imp_col_internal(mm, NA_real_))
case("mic_idx_nan", mean_imp_col_internal(mm, NaN))
case("mic_idx_ncol", mean_imp_col_internal(mm, 5L))
case("mic_idx_ncol_frac", mean_imp_col_internal(mm, 4.999))
case("mic_idx_big", mean_imp_col_internal(mm, 1e10))
case("mic_idx_huge", mean_imp_col_internal(mm, 1e20))
case("mic_idx_neg_huge", mean_imp_col_internal(mm, -1e20))
case("mic_idx_inf", mean_imp_col_internal(mm, Inf))
case("mic_idx_neg_inf", mean_imp_col_internal(mm, -Inf))
case("mic_idx_lgl", mean_imp_col_internal(mm, TRUE))
case("mic_inf_col", mean_imp_col_internal(mm, 3L))
case("mic_inf_unselected", mean_imp_col_internal(mm, c(0L, 2L)))
case("mic_cores_2", mean_imp_col_internal(mm, c(0L, 4L), 2L))
case("mic_0rows", mean_imp_col_internal(m0r, 0L))
case("mic_0cols", mean_imp_col_internal(m0c, integer(0)))
case("mic_int", mean_imp_col_internal(mi, c(0L, 1L)))
case("mic_vec", mean_imp_col_internal(vec, 0L))

case("ci_inf", check_inf(mnan))
case("ci_neg_inf", check_inf(matrix(c(1, NA, -Inf), 3, 1)))
case("ci_0x0", check_inf(m00))
case("ci_int", check_inf(mi))
case("ci_vec", check_inf(vec))

case("cf_allna_before_inf", check_finite(mnan))
case("cf_inf_before_allna", check_finite(inf_late))
case("cf_0rows", check_finite(m0r))
case("cf_0x0", check_finite(m00))
case("cf_nan_col", check_finite(cbind(c(NaN, NaN))))
case("cf_int", check_finite(mi))
case("cf_vec", check_finite(vec))

case("colm_nan", col_miss_internal(mnan))
case("colm_0rows", col_miss_internal(m0r))
case("colm_0cols", col_miss_internal(m0c))
case("colm_0x0", col_miss_internal(m00))
case("colm_int", col_miss_internal(mi))
case("colm_lgl", col_miss_internal(ml))
case("colm_vec", col_miss_internal(vec))

case("rowm_nan", row_miss_internal(mnan))
case("rowm_0rows", row_miss_internal(m0r))
case("rowm_0cols", row_miss_internal(m0c))
case("rowm_0x0", row_miss_internal(m00))
case("rowm_int", row_miss_internal(mi))
case("rowm_lgl", row_miss_internal(ml))
case("rowm_vec", row_miss_internal(vec))

cat(sprintf(
  "%d cases\n",
  length(list.files(out_dir, pattern = "^[^_].*\\.rds$"))
))
