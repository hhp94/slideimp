# Step 6 of Phase 2 (dev/eigen_conversion.md): instrument the SVD that the
# diagnostics build of src/armaSVD.cpp takes each iteration, in a scratch
# clone, never in the working tree. dev/probe_p2s6_run.sh copies this into a
# clone and runs it from the clone's root:
#   Rscript dev/probe_p2s6_patch.R --in-scratch-clone
# The patch is chosen by what the file holds: arma::svd() of the product, as
# up to step 5, or si_la::matmul and si_la::svd_values, as from step 6.
#
# SLIDEIMP_P2S6_FAIL = "<site>:<k>#<tag>" makes the k-th SVD since the
# variable last changed fail:
#   svdnan   a NaN is written into V_prev before the product, so the SVD is
#            handed a matrix holding NaN (both versions)
#   svdinfo  the status of si_la::svd_values is set to 1, as dgesdd reports a
#            failure to converge (step 6 only: arma::svd() offers no way in,
#            so in the older file the case is a plain repeat)
# Nothing changes unless the variable is set.
args <- commandArgs(trailingOnly = TRUE)
if (
  length(args) != 1L ||
    args[[1]] != "--in-scratch-clone" ||
    !file.exists("src/armaSVD.cpp")
) {
  stop(
    "run from the root of a scratch clone, with --in-scratch-clone; ",
    "never in the working tree"
  )
}
path <- "src/armaSVD.cpp"
txt <- readLines(path)

## Replace the one line that equals `anchor` by `lines`.
swap <- function(txt, anchor, lines) {
  at <- which(txt == anchor)
  if (length(at) != 1L) {
    stop(sprintf("%d lines equal to: %s", length(at), anchor))
  }
  append(txt[-at], lines, after = at - 1L)
}

helper <- c(
  "#include <string>",
  "#include <cstdlib>",
  "static bool p2s6_hit(const char *site)",
  "{",
  "  static std::string last;",
  "  static long count = 0;",
  "  const char *e = std::getenv(\"SLIDEIMP_P2S6_FAIL\");",
  "  const std::string cur = e ? e : \"\";",
  "  if (cur != last)",
  "  {",
  "    last = cur;",
  "    count = 0;",
  "  }",
  "  const std::size_t colon = cur.find(':');",
  "  if (colon == std::string::npos || cur.substr(0, colon) != site)",
  "  {",
  "    return false;",
  "  }",
  "  return ++count == std::atol(cur.c_str() + colon + 1);",
  "}"
)
txt <- swap(txt, "#include <string>", helper)

nan_line <- "std::numeric_limits<double>::quiet_NaN();"
old_anchor <- "        arma::vec s = arma::svd(V_prev.t() * V_new);"
new_anchor <-
  "      int st = si_la::matmul(true, 1.0, V_prev.data(), n_gram, n_cp,"
info_anchor <-
  "      st = si_la::svd_values(n_cp, n_cp, V_cross.data(), V_cross_sv.data());"
if (old_anchor %in% txt) {
  version <- "arma::svd"
  txt <- swap(
    txt,
    old_anchor,
    c(
      "        if (p2s6_hit(\"svdnan\"))",
      paste0("          V_prev(0, 0) = ", nan_line),
      old_anchor
    )
  )
} else if (new_anchor %in% txt) {
  version <- "si_la::svd_values"
  txt <- swap(
    txt,
    new_anchor,
    c(
      "      if (p2s6_hit(\"svdnan\"))",
      paste0("        V_prev[0] = ", nan_line),
      new_anchor
    )
  )
  txt <- swap(
    txt,
    info_anchor,
    c(info_anchor, "      if (p2s6_hit(\"svdinfo\"))", "        st = 1;")
  )
} else {
  stop("src/armaSVD.cpp holds neither form of the diagnostics SVD")
}
writeLines(txt, path)
cat("patched src/armaSVD.cpp, the version calling", version, "\n")
