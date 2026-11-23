# Install a release build of slideimp WITH the phase timers into its own
# library, for dev/gate/profile_phases.R. Dev only, run by hand:
#   Rscript dev/gate/profile_install.R <source dir> <lib>
# <source dir> is a clone outside this tree, used for nothing else while this
# runs: R CMD INSTALL builds in place. <lib> is created if missing.
#
# load_all1() also turns the timers on, but it builds without optimization,
# so its phase shares are not those of an installed package. This is plain
# R CMD INSTALL, with -DLOC_TIMER and rcpptimer's include directory appended
# to PKG_CPPFLAGS through the user Makevars, as load_all1() appends them:
# whatever the user Makevars already sets (the optimization level) stays.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "usage: Rscript dev/gate/profile_install.R <source dir> <lib>",
    call. = FALSE
  )
}
src <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
lib <- args[[2]]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
lib <- normalizePath(lib, winslash = "/", mustWork = TRUE)

inc <- system.file("include", package = "rcpptimer")
if (!nzchar(inc)) {
  stop("rcpptimer is not installed", call. = FALSE)
}
flags <- paste("-DLOC_TIMER", sprintf('-I"%s"', inc))

r_exe <- file.path(R.home("bin"), "R")
status <- withr::with_makevars(
  c(PKG_CPPFLAGS = flags),
  system2(
    r_exe,
    c(
      "CMD",
      "INSTALL",
      "--preclean",
      paste0("--library=", shQuote(lib)),
      shQuote(src)
    )
  ),
  assignment = "+="
)
if (!identical(status, 0L)) {
  stop("R CMD INSTALL failed with status ", status, call. = FALSE)
}
cat(
  "installed slideimp",
  as.character(packageVersion("slideimp", lib.loc = lib)),
  "with timers into",
  lib,
  "\n"
)
