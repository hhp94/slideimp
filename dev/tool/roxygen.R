# Driver for dev/tool/roxygen_utils.R. Run from the package root:
#
#   Rscript dev/tool/roxygen.R lint [--fill] [--dir=<path>]
#   Rscript dev/tool/roxygen.R wrap
#   Rscript dev/tool/roxygen.R md <path>
#
# lint     lint_roxygen() over R/ against dev/WRITING.md, sections 2, 4 and
#          6, and lint_seealso() over man/. It does not run document(), so
#          man/ is whatever the last devtools::document() wrote. --fill also
#          reports each block not filled to 80 columns.
#          --dir=<path> lints another R/ directory and skips the @seealso
#          check.
# wrap     refills every #' paragraph in R/ to 80 columns, in place, and
#          names the files it changed. It rewrites slideimp's own roxygen
#          layout, so run it only when that is wanted, then document() and
#          compare the rendered pages.
# md       writes every roxygen block in R/, with its signature, to <path>.
#
# R/RcppExports.R, which is generated, is skipped by all three.

if (
  !file.exists("DESCRIPTION") ||
    !identical(
      unname(read.dcf("DESCRIPTION", fields = "Package")[1, 1]),
      "slideimp"
    )
) {
  stop("run from the slideimp package root", call. = FALSE)
}
source(file.path("dev", "tool", "roxygen_utils.R"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("give a command: lint, wrap or md <path>", call. = FALSE)
}
cmd <- args[[1L]]
cat("skipped:", paste(ROXYGEN_SKIP, collapse = ", "), "\n\n")

if (cmd == "lint") {
  dir_arg <- grep("^--dir=", args[-1L], value = TRUE)
  fill <- "--fill" %in% args[-1L]
  dir <- if (length(dir_arg) == 1L) sub("^--dir=", "", dir_arg) else "R"
  rox <- lint_roxygen(dir, fill = fill)
  cat("lint_roxygen:", nrow(rox), "problem(s) in", dir, "\n")
  if (nrow(rox) > 0L) {
    print(rox, right = FALSE)
    cat("\nby problem:\n")
    print(table(sub(":.*", "", rox$problem)))
  }
  if (length(dir_arg) == 0L) {
    see <- lint_seealso("man")
    cat("\nlint_seealso:", nrow(see), "finding(s)\n")
    if (nrow(see) > 0L) {
      print(see, right = FALSE)
    }
  }
} else if (cmd == "wrap") {
  changed <- wrap_roxygen("R")
  cat("rewrapped", length(changed), "file(s)\n")
  if (length(changed) > 0L) {
    cat(paste0("  ", basename(changed)), sep = "\n")
  }
} else if (cmd == "md") {
  if (length(args) != 2L) {
    stop("md needs an output path", call. = FALSE)
  }
  path <- roxygen_md(args[[2L]], "R")
  cat("wrote", path, "\n")
} else {
  stop("unknown command: ", cmd, call. = FALSE)
}
