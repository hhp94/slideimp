# Run some test files against load_all1(timer = FALSE), with a per-test
# table and a TOTAL line.
# Usage: Rscript dev/tool/run_test_files.R [manual] all | <name> [<name> ...]
#   <name>: a file under tests/testthat/, with or without the "test-" prefix
#           and ".R" suffix, e.g. knn_imp or test-knn_imp.R
#   all:    every test file
#   manual: sets MANUAL_TESTS=true, unlocking skip_if_not_manual()
# timer = FALSE because the sampler refuses more than one chain under
# -DLOC_TIMER, which the mipca tests need. The tests that start mirai daemons
# run the INSTALLED slideimp there: set R_LIBS, before R starts, to a library
# holding an install of this checkout (dev/tool/test_mirai_with_dev_install.R
# makes one).

args <- commandArgs(trailingOnly = TRUE)
if ("manual" %in% args) {
  Sys.setenv(MANUAL_TESTS = "true")
  args <- setdiff(args, "manual")
}
if (length(args) == 0L) {
  stop(
    "usage: Rscript dev/tool/run_test_files.R [manual] all | <name> ...",
    call. = FALSE
  )
}
files <- if (identical(args, "all")) {
  list.files(
    file.path("tests", "testthat"),
    pattern = "^test-.*\\.R$",
    full.names = TRUE
  )
} else {
  vapply(
    args,
    function(a) {
      a <- sub("\\.R$", "", sub("^test-", "", basename(a)))
      file.path("tests", "testthat", paste0("test-", a, ".R"))
    },
    character(1)
  )
}
missing <- files[!file.exists(files)]
if (length(missing)) {
  stop("no such test file: ", paste(missing, collapse = ", "), call. = FALSE)
}
cat(
  "HEAD:",
  system("git rev-parse --short HEAD", intern = TRUE),
  " MANUAL_TESTS:",
  Sys.getenv("MANUAL_TESTS"),
  " R_LIBS:",
  Sys.getenv("R_LIBS"),
  "\n"
)

source("dev/dev-utils.R")
load_all1(timer = FALSE)

res <- lapply(files, function(f) {
  as.data.frame(testthat::test_file(
    f,
    reporter = "summary",
    stop_on_failure = FALSE,
    load_package = "none"
  ))
})
df <- do.call(rbind, res)
cat("\n=== per test\n")
print(
  df[, c("file", "test", "nb", "failed", "skipped", "error")],
  right = FALSE
)
cat(
  "\nTOTAL tests:",
  nrow(df),
  " failed:",
  sum(df$failed > 0),
  " errored:",
  sum(df$error),
  " skipped:",
  sum(df$skipped),
  "\n"
)
