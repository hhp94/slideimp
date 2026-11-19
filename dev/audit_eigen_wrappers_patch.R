# Audit of dev/review_eigen_wrappers.md, finding 1: instrument
# src/eig_sym_sel.h in a scratch clone, never in the working tree.
# dev/audit_eigen_wrappers_chain.sh copies this into the clone and runs it
# from the clone's root:
#   Rscript dev/audit_eigen_wrappers_patch.R --in-scratch-clone
# dsyevr is handed a W of n cells filled with a sentinel instead of the
# chain's eigvals (topk cells); the first topk cells are copied to eigvals
# afterwards, so the results are unchanged. For every call, one line goes to
# the file named by SLIDEIMP_AUDIT_LOG: n, topk, info, m_out, and how many
# cells past topk the call changed. In the unpatched chain those cells lie
# past the end of eigvals.
f <- "src/eig_sym_sel.h"
## It rewrites a source file, so it runs only when told it is in a clone.
if (
  !file.exists(f) ||
    !identical(commandArgs(trailingOnly = TRUE), "--in-scratch-clone")
) {
  stop(
    "run from the root of a scratch clone, with --in-scratch-clone; ",
    "never in the working tree"
  )
}
txt <- paste(readLines(f), collapse = "\n")
rep1 <- function(old, new) {
  hits <- gregexpr(old, txt, fixed = TRUE)[[1]]
  n <- if (identical(as.integer(hits), -1L)) 0L else length(hits)
  if (n != 1L) {
    stop(sprintf("expected 1 occurrence, found %d: %s", n, old))
  }
  txt <<- sub(old, new, txt, fixed = TRUE)
}
rep1(
  "#include <cstring>",
  "#include <cstring>\n#include <cstdio>\n#include <cstdlib>"
)
rep1(
  "  eigvecs.set_size(ws.n, ws.topk);",
  paste(
    "  eigvecs.set_size(ws.n, ws.topk);",
    "  std::vector<double> audit_w(static_cast<std::size_t>(ws.n), -777.25);",
    sep = "\n"
  )
)
rep1(
  "          eigvals.memptr(), eigvecs.memptr(), &ws.ldc,",
  "          audit_w.data(), eigvecs.memptr(), &ws.ldc,"
)
rep1(
  "  if (info != 0 || m_out != ws.topk)",
  paste(
    "  {",
    "    int audit_past = 0;",
    "    for (std::size_t i = static_cast<std::size_t>(ws.topk);",
    "         i < static_cast<std::size_t>(ws.n); ++i)",
    "    {",
    "      if (!(audit_w[i] == -777.25))",
    "      {",
    "        ++audit_past;",
    "      }",
    "    }",
    "    std::copy(audit_w.begin(), audit_w.begin() + ws.topk, eigvals.memptr());",
    "    if (const char *audit_f = std::getenv(\"SLIDEIMP_AUDIT_LOG\"))",
    "    {",
    "      if (std::FILE *fp = std::fopen(audit_f, \"a\"))",
    "      {",
    "        std::fprintf(fp, \"%d %d %d %d %d\\n\", static_cast<int>(ws.n),",
    "                     static_cast<int>(ws.topk), static_cast<int>(info),",
    "                     static_cast<int>(m_out), audit_past);",
    "        std::fclose(fp);",
    "      }",
    "    }",
    "  }",
    "  if (info != 0 || m_out != ws.topk)",
    sep = "\n"
  )
)
writeLines(txt, f)
cat("patched", f, "\n")
