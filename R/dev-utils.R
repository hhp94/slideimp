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
