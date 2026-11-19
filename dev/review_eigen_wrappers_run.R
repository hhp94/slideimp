# Review of the wrapper layer and the oracle (dev/review_eigen_wrappers.md):
# run the author's oracle scripts as they are, each in a child R process,
# and keep the whole output. Dev only, run by hand, from the package root:
#   Rscript dev/review_eigen_wrappers_run.R <lib> <out dir> <what>...
# <lib>      the private RcppArmadillo library (dev/lib/arma on this machine).
# <out dir>  where logs and builds go; nothing is built in the working tree.
# <what>     any of: check, oracle, selftest
#   check     dev/oracle_install_armadillo.R <lib> --check (reads CRAN's
#             index, installs nothing)
#   oracle    dev/oracle_run.R <lib> --verbose
#   selftest  dev/oracle_selftest.R <lib>
# Rcpp comes from the usual library: set R_LIBS_USER before R starts when
# the shell's HOME is not the one R was set up under.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "usage: Rscript dev/review_eigen_wrappers_run.R <lib> <out dir> <what>...",
    call. = FALSE
  )
}
lib <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
what <- args[-(1:2)]
if (!file.exists("dev/oracle_run.R")) {
  stop("run from the package root", call. = FALSE)
}
rscript <- file.path(R.home("bin"), "Rscript")

run <- function(label, script_args, tail_n) {
  log <- file.path(out_dir, paste0(label, ".log"))
  t0 <- Sys.time()
  out <- suppressWarnings(system2(
    rscript,
    script_args,
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  writeLines(c(out, sprintf("exit: %d", status)), log)
  cat(sprintf(
    "== %s: exit %d, %.0f s, %d lines, log %s\n",
    label,
    status,
    secs,
    length(out),
    log
  ))
  cat(utils::tail(out, tail_n), sep = "\n")
  cat("\n")
  invisible(status)
}

for (w in what) {
  if (w == "check") {
    run(
      "arma-check",
      c("dev/oracle_install_armadillo.R", shQuote(lib), "--check"),
      10L
    )
  } else if (w == "oracle") {
    run(
      "oracle",
      c(
        "dev/oracle_run.R",
        shQuote(lib),
        shQuote(paste0("--build=", file.path(out_dir, "oracle-build"))),
        "--verbose"
      ),
      110L
    )
  } else if (w == "selftest") {
    run(
      "selftest",
      c(
        "dev/oracle_selftest.R",
        shQuote(lib),
        shQuote(paste0("--build=", file.path(out_dir, "selftest-build")))
      ),
      70L
    )
  } else {
    stop("unknown <what>: ", w, call. = FALSE)
  }
}
