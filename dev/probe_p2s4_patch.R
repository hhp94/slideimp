# Instrument the LOBPCG path in a scratch clone, never in the working tree.
# The `edges` stage of dev/eigen_p4_gate.sh copies this into a clone and runs
# it from the clone's root:
#   Rscript dev/probe_p2s4_patch.R --in-scratch-clone old
#   Rscript dev/probe_p2s4_patch.R --in-scratch-clone new
# The mode only names the clone. It was written for step 4 of Phase 2
# (dev/eigen_conversion.md at 4f2b669) and patched every form the sources had
# through the conversion. It now patches the sources as clang-format lays
# them out, which is every commit from 571d176 on, with either form of the
# update after a QR: dtrtrs on the transposed operands (si_la::solve_tri), or
# one dtrsm from the right (si_la::solve_tri_right). For an older file, the
# Armadillo ones included, use `git show 2eb46bb:dev/probe_p2s4_patch.R`.
# Every text it replaces must occur exactly once, or it stops. Every version
# gets the same things, driven by environment variables read at run time, so
# one build serves every run:
#
# - SLIDEIMP_P2S4_LOG names a file that gets one line per event: every QR,
#   triangular solve and dsyevd the path makes, each failure, and the status
#   of every lobpcg_solve() call. A line starts with SLIDEIMP_P2S4_CASE.
# - SLIDEIMP_P2S4_FAIL = "<site>:<k>#<tag>" makes the k-th call at <site>
#   since the variable last changed fail, the way that site can fail:
#     qr        si_la::qr_econ returns 1
#     syevd     dsyevd reports info = 1
#     solve     the triangular solve returns 1
#     syevdnan  a NaN is written into the small Gram before dsyevd
#     gramnan   a NaN is written into the Gram handed to lobpcg_solve()
#     singular  the last pivot of the triangle handed to the solve is set
#               to 0, so the solve reports the zero and the block counts as
#               collapsed
#     resnan    the first column of A * X is set to NaN before the residuals
#               of a LOBPCG iteration, so one residual norm is NaN
#     eiginfo   the exact solver reports LAPACK's info = 1
#
# None of it changes a result unless SLIDEIMP_P2S4_FAIL is set. It also makes
# dev/ref_capture.R name each case in SLIDEIMP_P2S4_CASE.
args <- commandArgs(trailingOnly = TRUE)
if (
  length(args) != 2L ||
    args[[1]] != "--in-scratch-clone" ||
    !args[[2]] %in% c("old", "new") ||
    !file.exists("src/lobpcg_warm.h")
) {
  stop(
    "run from the root of a scratch clone, with --in-scratch-clone old|new; ",
    "never in the working tree"
  )
}
mode <- args[[2]]
lobpcg_txt <- paste(readLines("src/lobpcg_warm.h"), collapse = "\n")
holds <- function(s) grepl(s, lobpcg_txt, fixed = TRUE)

patch_file <- function(f, edits) {
  txt <- paste(readLines(f), collapse = "\n")
  for (e in edits) {
    hits <- gregexpr(e$old, txt, fixed = TRUE)[[1]]
    n <- if (identical(as.integer(hits), -1L)) 0L else length(hits)
    if (n != 1L) {
      stop(sprintf("%s: expected 1 occurrence, found %d: %s", f, n, e$old))
    }
    txt <- gsub(e$old, e$new, txt, fixed = TRUE)
  }
  writeLines(txt, f)
}
ed <- function(old, new) list(old = old, new = new)
lines <- function(...) paste(..., sep = "\n")

harness <- lines(
  "#include <cstdio>",
  "#include <cstdlib>",
  "#include <cstring>",
  "#include <string>",
  "namespace p2s4 {",
  "inline std::FILE *log_file() {",
  "  static std::FILE *f = nullptr;",
  "  static bool tried = false;",
  "  if (!tried) {",
  "    tried = true;",
  "    const char *p = std::getenv(\"SLIDEIMP_P2S4_LOG\");",
  "    if (p && *p) {",
  "      f = std::fopen(p, \"a\");",
  "    }",
  "  }",
  "  return f;",
  "}",
  "template <typename... T> inline void log(const char *fmt, T... args) {",
  "  if (std::FILE *f = log_file()) {",
  "    const char *c = std::getenv(\"SLIDEIMP_P2S4_CASE\");",
  "    std::fprintf(f, \"%s \", (c && *c) ? c : \"-\");",
  "    std::fprintf(f, fmt, args...);",
  "    std::fflush(f);",
  "  }",
  "}",
  "inline bool inject(const char *site) {",
  "  static std::string spec;",
  "  static long count = 0;",
  "  const char *p = std::getenv(\"SLIDEIMP_P2S4_FAIL\");",
  "  const std::string now = p ? p : \"\";",
  "  if (now != spec) {",
  "    spec = now;",
  "    count = 0;",
  "  }",
  "  const std::size_t len = std::strlen(site);",
  "  if (spec.compare(0, len, site) != 0 || spec.size() <= len ||",
  "      spec[len] != ':') {",
  "    return false;",
  "  }",
  "  ++count;",
  "  const bool hit = (count == std::atol(spec.c_str() + len + 1));",
  "  if (hit) {",
  "    log(\"injected %s %ld\\n\", site, count);",
  "  }",
  "  return hit;",
  "}",
  "inline unsigned long long ull(const double v) {",
  "  return static_cast<unsigned long long>(v);",
  "}",
  "} // namespace p2s4"
)
nan_cell <- "std::numeric_limits<double>::quiet_NaN()"

## src/hybrid_topk_eig.h: the status of every lobpcg_solve() call, the NaN in
## the Gram, and the exact solver made to report LAPACK's info = 1.
solve_line <- paste0(
  "    lob_res = lobpcg_solve(A, a_rows, a_cols, ctx.lobpcg_state, ",
  "ctx.lobpcg_opt);"
)
exact_status <- paste0(
  "  const int exact_st = eig_sym_sel(eigvals, eigvecs, A, a_rows, a_cols, ",
  "rng);"
)
patch_file(
  "src/hybrid_topk_eig.h",
  list(
    ed(
      exact_status,
      lines(
        "  const int exact_st =",
        "      p2s4::inject(\"eiginfo\")",
        "          ? 1",
        "          : eig_sym_sel(eigvals, eigvecs, A, a_rows, a_cols, rng);"
      )
    ),
    ed(
      solve_line,
      lines(
        "    if (p2s4::inject(\"gramnan\")) {",
        sprintf("      A[0] = %s;", nan_cell),
        "    }",
        solve_line,
        "    p2s4::log(\"lobpcg status %d iterations %d fail_iter %d k %llu\\n\",",
        "              lob_res.status, lob_res.iterations, lob_res.fail_iter,",
        "              p2s4::ull(ctx.lobpcg_state.X.ncol));"
      )
    )
  )
)
patch_file(
  "dev/ref_capture.R",
  list(
    ed(
      "case <- function(name, expr) {",
      "case <- function(name, expr) {\n  Sys.setenv(SLIDEIMP_P2S4_CASE = name)"
    )
  )
)

## The update after a QR, by what the file holds. Each form logs the solve
## with the fields dev/probe_p2s4_summary.R reads, rcond -1 standing for "not
## estimated".
log_solve <- lines(
  "    p2s4::log(\"solve %llu %llu rcond -1 approx 0\\n\", p2s4::ull(r),",
  "              p2s4::ull(AQ.nrow));"
)
solve_right <- lines(
  "    if (la_failed(",
  "            si_la::solve_tri_right('U', AQ.nrow, AQ.ncol, R.data(), AQ.data()),",
  "            \"dtrsm\")) {"
)
solve_left <- lines(
  "    if (la_failed(si_la::solve_tri('L', Rt.nrow, Y.ncol, Rt.data(), Y.data()),",
  "                  \"dtrtrs\")) {"
)
solve_ed <- if (holds(solve_right)) {
  ed(
    solve_right,
    lines(
      "    const int p2s4_st =",
      "        p2s4::inject(\"solve\")",
      "            ? 1",
      "            : si_la::solve_tri_right('U', AQ.nrow, AQ.ncol, R.data(),",
      "                                     AQ.data());",
      log_solve,
      "    if (la_failed(p2s4_st, \"dtrsm\")) {"
    )
  )
} else if (holds(solve_left)) {
  ed(
    solve_left,
    lines(
      "    const int p2s4_st =",
      "        p2s4::inject(\"solve\")",
      "            ? 1",
      "            : si_la::solve_tri('L', Rt.nrow, Y.ncol, Rt.data(), Y.data());",
      log_solve,
      "    if (la_failed(p2s4_st, \"dtrtrs\")) {"
    )
  )
} else {
  stop("src/lobpcg_warm.h holds neither form of the triangular solve")
}
from_right <- holds(solve_right)

## Both forms check the columns of A*Q before either reads R: the zero goes
## into the last pivot there.
cols_guard <- lines(
  "      Rcpp::stop(\"internal error: A*Q and Q differ in columns\");",
  "    }"
)
aq_line <- "    AQ = mul(A, n, n, Q);"
qr_call <- lines(
  "  if (la_failed(si_la::qr_econ(v_rows, v_cols, V.data(), Q.data(), R.data()),",
  "                \"dgeqrf or dorgqr\")) {"
)
syevd_call <- paste0(
  "  if (la_failed(si_la::syevd('U', n, M.data(), w.data()), \"dsyevd\")) {"
)
status_ok <- lines("  if (status == si_la::ok) {", "    return false;", "  }")
xonly_line <- "    // no usable P this call; fall back to X-only initial RR."

patch_file(
  "src/lobpcg_warm.h",
  list(
    ed(
      "#include \"si_lapack.h\"",
      lines("#include \"si_lapack.h\"", harness)
    ),
    ## The first Rayleigh-Ritz taken over X alone, because P did not survive
    ## its QR.
    ed(xonly_line, lines(xonly_line, "    p2s4::log(\"xonly_rr\\n\");")),
    ed(
      "    residual(AX, X, lam, R);",
      lines(
        "    if (p2s4::inject(\"resnan\")) {",
        "      for (std::size_t p2s4_i = 0; p2s4_i < AX.nrow; ++p2s4_i) {",
        sprintf("        AX(p2s4_i, 0) = %s;", nan_cell),
        "      }",
        "    }",
        "    residual(AX, X, lam, R);"
      )
    ),
    ed(
      cols_guard,
      lines(
        cols_guard,
        "    if (p2s4::inject(\"singular\")) {",
        "      R(R.nrow - 1, R.ncol - 1) = 0.0;",
        "    }"
      )
    ),
    ed(
      status_ok,
      lines(
        status_ok,
        "  p2s4::log(\"status %d from %s\\n\", status, routine);"
      )
    ),
    ed(
      qr_call,
      lines(
        "  p2s4::log(\"qr %llu %llu\\n\", p2s4::ull(V.nrow), p2s4::ull(V.ncol));",
        "  if (la_failed(p2s4::inject(\"qr\")",
        "                    ? 1",
        "                    : si_la::qr_econ(v_rows, v_cols, V.data(), Q.data(),",
        "                                     R.data()),",
        "                \"dgeqrf or dorgqr\")) {",
        "    p2s4::log(\"qr_failed\\n\");"
      )
    ),
    solve_ed,
    ed(
      aq_line,
      lines(
        "    p2s4::log(\"qr_dropped %llu of %llu\\n\", p2s4::ull(r),",
        "              p2s4::ull(old_cols));",
        aq_line
      )
    ),
    ed(
      syevd_call,
      lines(
        "  if (p2s4::inject(\"syevdnan\")) {",
        sprintf("    M(0, 0) = %s;", nan_cell),
        "  }",
        "  p2s4::log(\"syevd %d\\n\", static_cast<int>(n));",
        "  if (la_failed(p2s4::inject(\"syevd\")",
        "                    ? 1",
        "                    : si_la::syevd('U', n, M.data(), w.data()),",
        "                \"dsyevd\")) {",
        "    p2s4::log(\"syevd_failed\\n\");"
      )
    )
  )
)
cat(
  "patched the",
  mode,
  "clone: the update after a QR is",
  if (from_right) "dtrsm from the right" else "dtrtrs on the transposes",
  "\n"
)
