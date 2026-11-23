# Setup for the MSA example at the top of README.qmd. A hidden chunk there
# sources this file, so the README shows only the calls a reader makes and
# their live output, while the setup that depends on this machine stays here.
#
# slideimp.extra is not in Imports or Suggests. README.qmd is build-ignored,
# so neither it nor this file reaches the tarball or R CMD check, and the
# README may call slideimp.extra.
#
# It needs, on the machine that renders the README:
# - slideimp.extra, from pak::pkg_install("hhp94/slideimp.extra").
# - The MSA manifest in dev/manifests, as for dev/check_slideimp_extra.R: the
#   Illumina CSV (MSA-48v1-0_20102838_A1.csv), or the cleaned copy that
#   slideimp.extra derives from it. Like everything else in dev/manifests it
#   is not tracked.
#
# It leaves in the calling environment:
# - MSA_beta_matrix: 20 simulated samples x every MSA probe outside chrM,
#   with beta values and a few NA. Dropping chrM gives the README's
#   group_imp() call a manifest group with no column to match.
#
# It also points slideimp.extra at dev/manifests for the rest of the session,
# so the README's `group = "MSA"` finds the manifest.

if (!file.exists("DESCRIPTION") || !file.exists("dev/readme_msa.R")) {
  stop("source this from the repository root, as README.qmd does")
}
if (!nzchar(system.file(package = "slideimp.extra"))) {
  stop(
    "package 'slideimp.extra' is not installed; ",
    "install it with pak::pkg_install(\"hhp94/slideimp.extra\")."
  )
}

slideimp.extra::set_slideimp_path("dev/manifests")
MSA_beta_matrix <- local({
  msa <- slideimp.extra::ilmn_manifest("MSA", verbose = FALSE)
  msa <- msa[msa$group != "M", ]

  set.seed(1234)
  m <- slideimp::sim_mat(n = 20, p = nrow(msa), rho = 0.75, beta = TRUE)$input
  colnames(m) <- msa$feature
  m[1, 1] <- NA
  m[3, 2] <- NA
  m[1, 4] <- NA
  m
})
