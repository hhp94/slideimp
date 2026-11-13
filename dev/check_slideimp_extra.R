# Does prep_groups() accept the manifests slideimp.extra supplies?
#
# This lives here rather than under tests/ because `slideimp.extra` is a GitHub
# package that is not in Imports or Suggests, so CRAN and CI cannot install
# it, and any `slideimp.extra::` call under tests/ makes R CMD check report an
# unstated dependency. Its fixture is a 338 MB Illumina download, which has no
# place in the package source either. dev/ is build-ignored, so nothing here
# reaches the build or the check.
#
# It needs, on the machine that runs it:
# - slideimp.extra, from pak::pkg_install("hhp94/slideimp.extra").
# - The MSA manifest CSV (MSA-48v1-0_20102838_A1.csv, an Illumina download) in
#   the manifest directory. slideimp.extra keeps the manifests it derives from
#   it in the same directory. That directory is dev/manifests unless given as
#   an argument. Like everything else in dev/ it is not tracked.
#
# Usage, from the repository root:
#   Rscript dev/check_slideimp_extra.R [MANIFEST_DIR]
# Exits non-zero if a check fails.

if (!file.exists("DESCRIPTION") || !file.exists("dev/dev-utils.R")) {
  stop("run this from the repository root: Rscript dev/check_slideimp_extra.R")
}

args <- commandArgs(trailingOnly = TRUE)
manifest_dir <- if (length(args) > 0L) args[[1L]] else "dev/manifests"
if (!dir.exists(manifest_dir)) {
  stop(
    "manifest directory '",
    manifest_dir,
    "' does not exist; see the header of this script."
  )
}

# check that slideimp.extra is installed WITHOUT loading it. It imports
# slideimp, so loading it now would load the installed slideimp, and
# load_all1() would then replace the namespace slideimp.extra registered its
# methods into.
if (!nzchar(system.file(package = "slideimp.extra"))) {
  stop(
    "package 'slideimp.extra' is not installed; ",
    "install it with pak::pkg_install(\"hhp94/slideimp.extra\")."
  )
}

source("dev/dev-utils.R")
suppressMessages(load_all1(timer = FALSE))
invisible(loadNamespace("slideimp.extra"))

testthat::test_that("prep_groups() accepts slideimp.extra manifests", {
  slideimp.extra::set_slideimp_path(manifest_dir)
  withr::defer(slideimp.extra::set_slideimp_path(NULL))

  msa <- slideimp.extra::ilmn_manifest(
    "MSA",
    deduped = TRUE,
    rawdir = manifest_dir
  )
  testthat::expect_gt(length(msa$feature), 0L)

  # the manifest object itself, and the name slideimp.extra registers for it
  testthat::expect_no_error(prep_groups(msa$feature, group = msa))
  testthat::expect_no_error(prep_groups(msa$feature, group = "MSA_deduped"))
})
