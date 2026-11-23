# Fails when a compiled object of the package holds Rcpp's as.list()
# fallback: a symbol naming Rcpp::Armor, convert_using_rfunction or
# r_true_cast<19> (the cast to VECSXP). An Rcpp::List built from a SEXP
# compiles it in, and rchk reports its Armor<SEXP> as errors (dev/CPP.md,
# "No Rcpp::List built from a SEXP").
#
# Usage, from the package root:
#   Rscript dev/check_rcpp_list.R [package dir] [--no-build]
#
# It cleans the package's objects and builds them first, with
# pkgbuild::compile_dll(debug = TRUE): -O0, no extra macros, the build of
# load_all1(timer = FALSE). At -O0 every inline function a translation unit
# uses is emitted as a symbol of its object. An optimized build inlines
# Armor<SEXP> and r_true_cast<19> away and keeps convert_using_rfunction only
# while the compiler declines to inline it, so it can read clean with the cast
# compiled in. --no-build
# reads the objects already in src/, for a build made by other means; they
# must come from that -O0 build. -DLOC_TIMER is left out because code under
# it is exempt from the rule. Run load_all1() again afterwards to get the
# timer build back.
#
# nm comes from PATH, else from Rtools (RTOOLS45_HOME, default C:/rtools45).

args <- commandArgs(trailingOnly = TRUE)
build <- !("--no-build" %in% args)
args <- setdiff(args, "--no-build")
pkg <- if (length(args)) args[[1]] else "."

nm <- Sys.which("nm")
if (!nzchar(nm)) {
  rtools <- Sys.getenv("RTOOLS45_HOME", "C:/rtools45")
  nm <- file.path(rtools, "x86_64-w64-mingw32.static.posix", "bin", "nm.exe")
  if (!file.exists(nm)) {
    stop("nm not found on PATH or in Rtools", call. = FALSE)
  }
}

if (build) {
  pkgbuild::clean_dll(pkg)
  pkgbuild::compile_dll(pkg, debug = TRUE, force = TRUE)
}

objs <- Sys.glob(file.path(pkg, "src", "*.o"))
if (!length(objs)) {
  stop("no object files in ", file.path(pkg, "src"), call. = FALSE)
}

pattern <- "Rcpp::Armor<|convert_using_rfunction|r_true_cast<19>"
hits <- 0L
for (o in objs) {
  syms <- system2(nm, c("-C", shQuote(o)), stdout = TRUE)
  bad <- unique(sub("^[0-9a-fA-F ]*[A-Za-z] ", "", grep(pattern, syms, value = TRUE)))
  if (length(bad)) {
    cat(basename(o), ":\n", paste0("  ", bad, "\n"), sep = "")
    hits <- hits + length(bad)
  }
}
cat(sprintf("%d object(s) read, %d symbol(s) of the as.list() fallback\n", length(objs), hits))
if (hits > 0L) {
  quit(save = "no", status = 1L)
}
