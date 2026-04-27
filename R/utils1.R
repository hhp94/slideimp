load_all1 <- function(timer = TRUE, debug = FALSE) {
  checkmate::assert_flag(timer)
  checkmate::assert_flag(debug)

  flags <- character()
  if (timer) flags <- c(flags, "-DLOC_TIMER")
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
    dir, output_file
  ))
}
