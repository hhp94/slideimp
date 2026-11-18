# Install the latest RcppArmadillo from CRAN into a private library, for the
# Armadillo oracle (dev/eigen_conversion.md, "The Armadillo oracle").
#   Rscript dev/oracle_install_armadillo.R <lib> [--check]
# It installs only when <lib> has no RcppArmadillo or an older one than
# CRAN's; otherwise it says so and leaves the library alone. --check reports
# the two versions and installs or creates nothing. Either way it reads
# CRAN's package index, which is a download.
#
# CRAN's version is the source package's. The Windows binary is used when it
# is the same version; when the binary lags, the source package is built,
# which needs Rtools. Dependencies are not installed into <lib>: Rcpp comes
# from the usual library. dev/oracle_run.R takes <lib> as its argument, puts
# it first on the library path and stops unless RcppArmadillo then resolves
# inside it, so R_LIBS_USER stays as it is.
#
# <lib> is regenerable and stays out of the repository. The oracle prints
# packageVersion("RcppArmadillo") on every run; this script prints the
# version it leaves in <lib>.

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
pos <- setdiff(args, "--check")
if (length(pos) != 1L) {
  stop(
    "usage: Rscript dev/oracle_install_armadillo.R <lib> [--check]",
    call. = FALSE
  )
}
lib <- pos[[1]]
## --check creates nothing, so a mistyped <lib> is not left behind.
if (!check_only) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
}
lib <- normalizePath(lib, winslash = "/", mustWork = !check_only)
cran <- "https://cloud.r-project.org"
pkg <- "RcppArmadillo"

installed_version <- function() {
  d <- file.path(lib, pkg)
  if (!dir.exists(d)) NULL else packageVersion(pkg, lib.loc = lib)
}
avail_version <- function(type) {
  a <- available.packages(repos = cran, type = type)
  if (pkg %in% rownames(a)) package_version(a[pkg, "Version"]) else NULL
}

have <- installed_version()
src_v <- avail_version("source")
bin_v <- if (.Platform$OS.type == "windows") avail_version("win.binary")
cat(sprintf(
  "%s in %s: %s\nCRAN source: %s%s\n",
  pkg,
  lib,
  if (is.null(have)) "none" else as.character(have),
  as.character(src_v),
  if (is.null(bin_v)) "" else paste0(", Windows binary: ", bin_v)
))
if (is.null(src_v)) {
  stop(pkg, " not found on CRAN", call. = FALSE)
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
    repos = cran,
    type = type,
    dependencies = FALSE
  )
  have <- installed_version()
  if (is.null(have) || have < src_v) {
    stop("install did not leave ", pkg, " ", src_v, " in ", lib, call. = FALSE)
  }
}
cat(sprintf("%s %s in %s\n", pkg, as.character(installed_version()), lib))
