# Install the RcppEigen candidate (Eigen 5) into its own private library
# (dev/eigen_conversion.md at 4f2b669, "Both RcppEigen versions"). Never over
# the CRAN one: a run selects the candidate by setting R_LIBS=<lib> before R
# starts.
#   Rscript dev/eigen_install_candidate.R <lib> [--check]
# The candidate comes from the RcppCore r-universe, the source
# https://github.com/RcppCore/RcppEigen/issues/151 names. It installs only
# when <lib> has no RcppEigen or an older one than the r-universe's;
# --check reports both versions and installs or creates nothing. Either way
# it reads the r-universe's package index, which is a download. Dependencies
# are not installed into <lib>; Rcpp comes from the usual library.
#
# It prints packageVersion("RcppEigen") and the Eigen version the installed
# headers define, so the log says which candidate was tested. <lib> is
# regenerable and stays out of the repository.

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
pos <- setdiff(args, "--check")
if (length(pos) != 1L) {
  stop(
    "usage: Rscript dev/eigen_install_candidate.R <lib> [--check]",
    call. = FALSE
  )
}
lib <- pos[[1]]
## --check creates nothing, so a mistyped <lib> is not left behind.
if (!check_only) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
}
lib <- normalizePath(lib, winslash = "/", mustWork = !check_only)
repo <- "https://rcppcore.r-universe.dev"
pkg <- "RcppEigen"

installed_version <- function() {
  d <- file.path(lib, pkg)
  if (!dir.exists(d)) NULL else packageVersion(pkg, lib.loc = lib)
}
avail_version <- function(type) {
  a <- available.packages(repos = repo, type = type)
  if (pkg %in% rownames(a)) package_version(a[pkg, "Version"]) else NULL
}

have <- installed_version()
src_v <- avail_version("source")
bin_v <- if (.Platform$OS.type == "windows") avail_version("win.binary")
cat(sprintf(
  "%s in %s: %s\nr-universe source: %s%s\n",
  pkg,
  lib,
  if (is.null(have)) "none" else as.character(have),
  as.character(src_v),
  if (is.null(bin_v)) "" else paste0(", Windows binary: ", bin_v)
))
if (is.null(src_v)) {
  stop(pkg, " not found at ", repo, call. = FALSE)
}

if (!is.null(have) && have >= src_v) {
  cat("up to date, nothing installed\n")
} else if (check_only) {
  cat("would install", as.character(src_v), "(--check: nothing installed)\n")
} else {
  type <- if (!is.null(bin_v) && bin_v == src_v) "win.binary" else "source"
  cat("installing", as.character(src_v), "as", type, "\n")
  install.packages(
    pkg,
    lib = lib,
    repos = repo,
    type = type,
    dependencies = FALSE
  )
  have <- installed_version()
  if (is.null(have) || have < src_v) {
    stop("install did not leave ", pkg, " ", src_v, " in ", lib, call. = FALSE)
  }
}

if (!is.null(installed_version())) {
  cat(sprintf("%s %s in %s\n", pkg, as.character(installed_version()), lib))
  ## Eigen 5 defines its version in Eigen/Version; Eigen 3.4 in Macros.h.
  inc <- file.path(lib, pkg, "include", "Eigen")
  files <- c(
    file.path(inc, "Version"),
    file.path(inc, "src", "Core", "util", "Macros.h")
  )
  v <- character()
  for (f in files[file.exists(files)]) {
    v <- c(
      v,
      grep(
        "^#define EIGEN_((WORLD|MAJOR|MINOR|PATCH)_VERSION|VERSION_STRING)\\b",
        readLines(f),
        value = TRUE
      )
    )
  }
  cat(
    "Eigen headers:",
    if (length(v)) v else "(no version macros found)",
    sep = "\n  "
  )
  cat("\n")
}
