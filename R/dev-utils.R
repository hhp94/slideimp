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

  old <- Sys.getenv("PKG_CPPFLAGS", unset = "")
  new <- paste(c(old, flags), collapse = " ")

  devtools::clean_dll()
  withr::with_envvar(c(PKG_CPPFLAGS = new), devtools::load_all())
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
