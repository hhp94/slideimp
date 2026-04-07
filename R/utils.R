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
