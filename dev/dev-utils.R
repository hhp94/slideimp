# Dev-only helpers. They live in dev/, which is build-ignored, so none of
# this ships, and they are not part of the package namespace: load them with
# source("dev/dev-utils.R") from the package root.
#
# load_all1() is the way to load the package in development, not bare
# load_all(), because it decides which compile-time macros the DLL is built
# with:
#
# - `timer = TRUE` is the default and adds -DLOC_TIMER plus rcpptimer's
#   include directory. src/loc_timer.h compiles the whole LOC_TIMER_* family
#   away without it, so a bare load_all() silently produces a package with no
#   phase timers at all.
# - `debug = TRUE` sets -DPCA_IMP_DIAGNOSTICS=1, which guards the diagnostic
#   blocks in src/armaSVD.cpp. The source defaults it to 0, so it is off in
#   every build that does not ask.
#
# It calls clean_dll() first, every time. That is deliberate: a package DLL
# keeps whatever flags built it until a source edit forces a recompile, so
# switching macros without a clean leaves the next run measuring something it
# did not choose, with no source file recording the difference. Never work
# around the clean to save a rebuild, and never mix a load_all1() build with a
# bare load_all() one inside a single comparison.
#
# check_manual() runs devtools::test() with MANUAL_TESTS=true, the gate that
# tests/testthat/helper.R checks, unlocking tests that are skipped by default.
# dump_roxygen2() prints an `rg` command for surveying roxygen blocks and
# signatures; it prints, it does not run. scratch() opens dev/scratch.R.

load_all1 <- function(timer = TRUE, debug = FALSE) {
  checkmate::assert_flag(timer)
  checkmate::assert_flag(debug)

  flags <- character()
  if (timer) {
    # rcpptimer is a dev-only dependency: src/loc_timer.h includes its header
    # only under -DLOC_TIMER. It is deliberately absent from LinkingTo so that
    # ordinary builds do not require it, which means its include directory has
    # to be added by hand here.
    inc <- system.file("include", package = "rcpptimer")
    if (!nzchar(inc)) {
      cli::cli_abort(c(
        "{.pkg rcpptimer} is required by {.code timer = TRUE}.",
        i = 'Install it, or call {.code load_all1(timer = FALSE)}.'
      ))
    }
    flags <- c(flags, "-DLOC_TIMER", sprintf('-I"%s"', inc))
  }
  flags <- c(flags, sprintf("-DPCA_IMP_DIAGNOSTICS=%d", as.integer(debug)))

  devtools::clean_dll()
  # src/Makevars assigns PKG_CPPFLAGS, and a makefile assignment overrides an
  # environment variable of the same name, so the macros cannot go in through
  # the environment. They go in through the user Makevars instead, which R
  # reads after the package's, as an append on top of -DARMA_64BIT_WORD=1.
  withr::with_makevars(
    c(PKG_CPPFLAGS = paste(flags, collapse = " ")),
    devtools::load_all(),
    assignment = "+="
  )
}

scratch <- function() {
  file.edit("dev/scratch.R")
}

dump_roxygen2 <- function(output_file = "roxygen2.txt", dir = "R/") {
  cat(sprintf(
    "rg -nU --multiline-dotall \"^#'|^[a-zA-Z0-9_\\.]+ <- function\\(.*?\\) ?\\{\" %s > %s\n",
    dir,
    output_file
  ))
}

check_manual <- function(...) {
  withr::with_envvar(
    c(MANUAL_TESTS = "true"),
    devtools::test(...)
  )
}
