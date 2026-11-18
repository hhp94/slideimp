# Edge inputs for the exports Phase 1 of dev/eigen_conversion.md converts:
# empty dims, integer and logical input, a vector that is not a matrix, NA,
# NaN and Inf, a zero held as both 0 and -0 at a column's extreme, and the
# sampler's odd arguments. The Phase 0 capture runs these exports on
# ordinary input only; this checks that the rewritten kernels compute what
# the Armadillo ones did on the inputs both accept.
#
# Run from the package root of a build, with a label naming the capture:
#   Rscript dev/eigen_p1_edges.R <label>
# It writes dev/ref/<label>/ in dev/ref_capture.R's format, so two builds'
# captures are compared with dev/ref_compare.R; compare at tolerance 0. Run
# it on the start commit in a clone, as the reference capture was, and on
# the converted build.
#
# Phase 1 changed the containers on purpose: vectors where Armadillo gave
# 1 x p or n x 1 matrices, integer counts and positions where it gave
# doubles. vals() reduces a result to its values as doubles, so the two
# builds compare on what they compute. Calls are positional and pass only
# arguments both builds take. Character input is left out: under a
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
cat("col_vars_internal formals:", names(formals(col_vars_internal)), "\n")

vals <- function(x) as.double(x)

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

case("cv_one_row", vals(col_vars_internal(one)))
case("cv_0rows", vals(col_vars_internal(m0r)))
case("cv_0cols", vals(col_vars_internal(m0c)))
case("cv_nan", vals(col_vars_internal(mnan)))
case("cv_int", vals(col_vars_internal(mi)))
case("cv_lgl", vals(col_vars_internal(ml)))
case("cv_vec", vals(col_vars_internal(vec)))

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

case("colm_nan", vals(col_miss_internal(mnan)))
case("colm_0rows", vals(col_miss_internal(m0r)))
case("colm_0cols", vals(col_miss_internal(m0c)))
case("colm_0x0", vals(col_miss_internal(m00)))
case("colm_int", vals(col_miss_internal(mi)))
case("colm_lgl", vals(col_miss_internal(ml)))
case("colm_vec", vals(col_miss_internal(vec)))

case("rowm_nan", vals(row_miss_internal(mnan)))
case("rowm_0rows", vals(row_miss_internal(m0r)))
case("rowm_0cols", vals(row_miss_internal(m0c)))
case("rowm_0x0", vals(row_miss_internal(m00)))
case("rowm_int", vals(row_miss_internal(mi)))
case("rowm_lgl", vals(row_miss_internal(ml)))
case("rowm_vec", vals(row_miss_internal(vec)))

## The sampler (src/sample_each_rep_cpp.cpp). Each case runs under
## set.seed(1), so the draws from R's generator must give identical
## positions.
set.seed(42)
sx <- matrix(round(rnorm(12 * 8), 1), 12, 8)
sx[c(2, 15, 40, 77)] <- NA
sx[, 8] <- 5 # constant: never usable
sargs <- list(
  pool = 0:7,
  need = c(2L, 3L, 3L),
  room = rep(2L, 12),
  cap = c(rep(4L, 7), 0L)
)
case(
  "ser_int",
  vals(sample_each_rep_cpp(
    sx,
    sargs$pool,
    sargs$need,
    sargs$room,
    sargs$cap,
    50L
  ))
)
case(
  "ser_double",
  vals(sample_each_rep_cpp(
    sx,
    as.numeric(sargs$pool),
    as.numeric(sargs$need),
    as.numeric(sargs$room),
    as.numeric(sargs$cap),
    50
  ))
)
case(
  "ser_tight_rows",
  vals(sample_each_rep_cpp(
    sx,
    sargs$pool,
    sargs$need,
    rep(1L, 12),
    sargs$cap,
    50L
  ))
)
case(
  "ser_infeasible",
  vals(sample_each_rep_cpp(
    sx,
    0:1,
    c(5L, 5L, 5L),
    rep(1L, 12),
    c(9L, 9L),
    20L
  ))
)
case(
  "ser_zero_attempts",
  vals(sample_each_rep_cpp(
    sx,
    sargs$pool,
    sargs$need,
    sargs$room,
    sargs$cap,
    0L
  ))
)
case(
  "ser_no_need",
  vals(sample_each_rep_cpp(
    sx,
    sargs$pool,
    integer(0),
    sargs$room,
    sargs$cap,
    5L
  ))
)
case(
  "ser_frac_attempts",
  vals(sample_each_rep_cpp(
    sx,
    sargs$pool,
    sargs$need,
    sargs$room,
    sargs$cap,
    2.9
  ))
)
case(
  "ser_int_matrix",
  vals(sample_each_rep_cpp(
    matrix(1:40, 10, 4),
    0:3,
    c(2L, 2L),
    rep(1L, 10),
    rep(3L, 4),
    10L
  ))
)

cat(sprintf(
  "%d cases\n",
  length(list.files(out_dir, pattern = "^[^_].*\\.rds$"))
))
