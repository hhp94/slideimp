load_all1 <- function(recompile = c("eig_sym_sel.h", "fastSVD.cpp")) {
  if (is.null(recompile)) {
    devtools::clean_dll()
  } else {
    for (f in recompile) {
      path <- file.path("src", f)
      Sys.setFileTime(path, Sys.time())
    }
  }
  withr::with_envvar(c(PKG_CPPFLAGS = "-DLOC_TIMER"), devtools::load_all())
}

scratch <- function() {
  file.edit("dev/scratch.R")
}

dump_roxygen2 <- function(output_file = "roxygen2.txt", dir = "R/") {
  cmd <- paste0('rg -n "^#\'|^[a-zA-Z0-9_\\.]+ <- function" ', dir, " > ", output_file)
  cat(cmd)
}
