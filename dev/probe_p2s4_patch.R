# Step 4 of Phase 2 (dev/eigen_conversion.md): instrument the LOBPCG path in
# a scratch clone, never in the working tree. dev/probe_p2s4_run.sh copies
# this into a clone and runs it from the clone's root:
#   Rscript dev/probe_p2s4_patch.R --in-scratch-clone old
#   Rscript dev/probe_p2s4_patch.R --in-scratch-clone new
# The patch is chosen by what src/lobpcg_warm.h holds, not by the mode, which
# only names the clone: the file as step 3 left it (arma::qr_econ,
# arma::solve, arma::dsyevd_, as at eade032), or a file whose three calls go
# through si_la (si_la::qr_econ, si_la::solve_tri, si_la::syevd): ff0e394,
# whose containers are still arma::mat, or the second half of step 4, whose
# are lobpcg_detail::Block. Two anchors differ between those two and are
# taken from whichever the file holds; Phase 4 renamed the members of Block,
# and nm() writes every anchor in the names the file holds. The call into it in
# src/hybrid_topk_eig.h is likewise found by its form: the Gram an arma::mat
# before step 5 of Phase 2, a raw buffer with its two dimensions from step 5
# on. Every version gets the same three
# things, driven by environment variables read at run time, so one build
# serves every run:
#
# - SLIDEIMP_P2S4_LOG names a file that gets one line per event: every QR,
#   triangular solve and dsyevd the path makes, each failure, and the status
#   of every lobpcg_solve() call. A line starts with SLIDEIMP_P2S4_CASE.
# - SLIDEIMP_P2S4_FAIL = "<site>:<k>#<tag>" makes the k-th call at <site>
#   since the variable last changed fail, the way that site can fail:
#     qr        arma::qr_econ() returns false / si_la::qr_econ returns 1
#     syevd     dsyevd reports info = 1
#     solve     arma::solve() throws / si_la::solve_tri returns 1
#     syevdnan  a NaN is written into the small Gram before dsyevd
#     gramnan   a NaN is written into the Gram handed to lobpcg_solve()
#     singular  the last pivot of the triangle handed to the solve is set
#               to 0, so the solve takes its least-squares fallback, or,
#               from Phase 4 on, dtrtrs reports the zero and the block
#               counts as collapsed
#     resnan    the first column of A * X is set to NaN before the residuals
#               of a LOBPCG iteration, so one residual norm is NaN (the
#               Block form of the file only)
#     eiginfo   the exact solver reports LAPACK's info = 1 (from the
#               raw-buffer Gram on)
# - SLIDEIMP_P2S4_FRESH, the file of step 3 only: dsyevd gets the workspace
#   length of its own query on every call, where step 3 kept the longest one
#   so far. si_la::syevd always does.
#
# None of it changes a result unless SLIDEIMP_P2S4_FAIL or
# SLIDEIMP_P2S4_FRESH is set. It also makes dev/ref_capture.R name each case
# in SLIDEIMP_P2S4_CASE.
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
## The first of the candidate lines that the file holds.
anchor <- function(...) {
  cand <- c(...)
  hit <- cand[vapply(cand, holds, TRUE)]
  if (!length(hit)) {
    stop("src/lobpcg_warm.h holds none of: ", paste(cand, collapse = " | "))
  }
  hit[[1]]
}
## Phase 4 renamed the members of lobpcg_detail::Block off Armadillo's
## names. Before, a Block and an arma::mat both had n_rows, n_cols and
## memptr(); a file that holds the new names gets every anchor and every
## line written below in nrow, ncol and data(), through ed().
block_new_names <- holds("    std::size_t nrow = 0;")
nm <- function(s) {
  if (!block_new_names) {
    return(s)
  }
  s <- gsub(".n_rows", ".nrow", s, fixed = TRUE)
  s <- gsub(".n_cols", ".ncol", s, fixed = TRUE)
  gsub(".memptr()", ".data()", s, fixed = TRUE)
}
arma_calls <- holds("arma::qr_econ(Q, R, V)")
## Phase 4 made si_la::solve_tri plain dtrtrs, without the two outputs.
plain_solve <- !arma_calls && !holds(nm("Y.memptr(), &rcond, &approx)"))

patch_file <- function(f, edits) {
  txt <- paste(readLines(f), collapse = "\n")
  for (e in edits) {
    hits <- gregexpr(e$old, txt, fixed = TRUE)[[1]]
    n <- if (identical(as.integer(hits), -1L)) 0L else length(hits)
    want <- if (is.null(e$n)) 1L else e$n
    if (n != want) {
      stop(sprintf(
        "%s: expected %d occurrence(s), found %d: %s",
        f,
        want,
        n,
        e$old
      ))
    }
    txt <- gsub(e$old, e$new, txt, fixed = TRUE)
  }
  writeLines(txt, f)
}
ed <- function(old, new, n = 1L) list(old = nm(old), new = nm(new), n = n)
lines <- function(...) paste(..., sep = "\n")

harness <- lines(
  "#include <cstdio>",
  "#include <cstdlib>",
  "#include <cstring>",
  "#include <stdexcept>",
  "#include <string>",
  "namespace p2s4",
  "{",
  "  inline std::FILE *log_file()",
  "  {",
  "    static std::FILE *f = nullptr;",
  "    static bool tried = false;",
  "    if (!tried)",
  "    {",
  "      tried = true;",
  "      const char *p = std::getenv(\"SLIDEIMP_P2S4_LOG\");",
  "      if (p && *p)",
  "      {",
  "        f = std::fopen(p, \"a\");",
  "      }",
  "    }",
  "    return f;",
  "  }",
  "  template <typename... T>",
  "  inline void log(const char *fmt, T... args)",
  "  {",
  "    if (std::FILE *f = log_file())",
  "    {",
  "      const char *c = std::getenv(\"SLIDEIMP_P2S4_CASE\");",
  "      std::fprintf(f, \"%s \", (c && *c) ? c : \"-\");",
  "      std::fprintf(f, fmt, args...);",
  "      std::fflush(f);",
  "    }",
  "  }",
  "  inline bool inject(const char *site)",
  "  {",
  "    static std::string spec;",
  "    static long count = 0;",
  "    const char *p = std::getenv(\"SLIDEIMP_P2S4_FAIL\");",
  "    const std::string now = p ? p : \"\";",
  "    if (now != spec)",
  "    {",
  "      spec = now;",
  "      count = 0;",
  "    }",
  "    const std::size_t len = std::strlen(site);",
  "    if (spec.compare(0, len, site) != 0 || spec.size() <= len ||",
  "        spec[len] != ':')",
  "    {",
  "      return false;",
  "    }",
  "    ++count;",
  "    const bool hit = (count == std::atol(spec.c_str() + len + 1));",
  "    if (hit)",
  "    {",
  "      log(\"injected %s %ld\\n\", site, count);",
  "    }",
  "    return hit;",
  "  }",
  "  inline unsigned long long ull(const double v)",
  "  {",
  "    return static_cast<unsigned long long>(v);",
  "  }",
  "} // namespace p2s4"
)
nan_cell <- "std::numeric_limits<double>::quiet_NaN()"

## Common to both: the status of every lobpcg_solve() call, and the NaN in
## the Gram. Until step 5 the Gram is an arma::mat; from step 5 it is a raw
## buffer and the call takes its two dimensions.
hybrid_txt <- paste(readLines("src/hybrid_topk_eig.h"), collapse = "\n")
solve_arma <- "    lob_res = lobpcg_solve(A, ctx.lobpcg_state, ctx.lobpcg_opt);"
solve_raw <- lines(
  "    lob_res = lobpcg_solve(A, a_rows, a_cols, ctx.lobpcg_state,",
  "                           ctx.lobpcg_opt);"
)
gram_arma <- grepl(solve_arma, hybrid_txt, fixed = TRUE)
if (!gram_arma && !grepl(solve_raw, hybrid_txt, fixed = TRUE)) {
  stop("src/hybrid_topk_eig.h holds neither form of the lobpcg_solve() call")
}
solve_line <- if (gram_arma) solve_arma else solve_raw
gram_cell <- if (gram_arma) "A(0, 0)" else "A[0]"
## The exact solver made to report LAPACK's info = 1. Phase 4 made the call
## return its status; before, it was tested in place. Files older than the
## raw-buffer Gram do not get this injection.
exact_status <- "  const int exact_st = eig_sym_sel(eigvals, eigvecs, A, a_rows, a_cols, rng);"
exact_flag <- "  if (eig_sym_sel(eigvals, eigvecs, A, a_rows, a_cols, rng) != si_la::ok)"
exact_ed <- if (grepl(exact_status, hybrid_txt, fixed = TRUE)) {
  list(ed(
    exact_status,
    lines(
      "  const int exact_st =",
      "      p2s4::inject(\"eiginfo\")",
      "          ? 1",
      "          : eig_sym_sel(eigvals, eigvecs, A, a_rows, a_cols, rng);"
    )
  ))
} else if (grepl(exact_flag, hybrid_txt, fixed = TRUE)) {
  list(ed(
    exact_flag,
    lines(
      "  if (p2s4::inject(\"eiginfo\") ||",
      "      eig_sym_sel(eigvals, eigvecs, A, a_rows, a_cols, rng) != si_la::ok)"
    )
  ))
} else {
  list()
}
patch_file(
  "src/hybrid_topk_eig.h",
  c(
    exact_ed,
    list(
      ed(
        solve_line,
        lines(
          "    if (p2s4::inject(\"gramnan\"))",
          "    {",
          sprintf("      %s = %s;", gram_cell, nan_cell),
          "    }",
          solve_line,
          "    p2s4::log(\"lobpcg status %d iterations %d fail_iter %d k %llu\\n\",",
          "              lob_res.status, lob_res.iterations, lob_res.fail_iter,",
          "              p2s4::ull(ctx.lobpcg_state.X.n_cols));"
        )
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

## Common to both versions of src/lobpcg_warm.h: the first Rayleigh-Ritz
## taken over X alone, because P did not survive its QR.
xonly <- ed(
  "    // no usable P this call; fall back to X-only initial RR.",
  lines(
    "    // no usable P this call; fall back to X-only initial RR.",
    "    p2s4::log(\"xonly_rr\\n\");"
  )
)
## A triangle singular to working precision, handed to the solve.
singular <- lines(
  "      if (p2s4::inject(\"singular\"))",
  "      {",
  "        R(R.n_rows - 1, R.n_cols - 1) = 0.0;",
  "      }"
)

log_qr <- lines(
  "    p2s4::log(\"qr %llu %llu\\n\", p2s4::ull(V.n_rows),",
  "              p2s4::ull(V.n_cols));"
)
## The matvec recomputed when the QR dropped columns, and the transpose of R
## the solve is handed: arma::mat until ff0e394, lobpcg_detail::Block after.
aq_line <- anchor("      AQ = A * Q;", "      AQ = mul(A, n, n, Q);")
rt_line <- if (arma_calls) {
  NULL
} else {
  anchor(
    "      const arma::mat Rt = R.t();",
    "      const Block Rt = transpose(R);"
  )
}
log_dropped <- lines(
  "      p2s4::log(\"qr_dropped %llu of %llu\\n\", p2s4::ull(r),",
  "                p2s4::ull(old_cols));",
  aq_line
)

if (arma_calls) {
  patch_file(
    "src/lobpcg_warm.h",
    list(
      ed(
        "#include \"pca_linalg_utils.h\"",
        lines("#include \"pca_linalg_utils.h\"", harness)
      ),
      xonly,
      ed(
        "    if (!arma::qr_econ(Q, R, V))\n    {",
        lines(
          log_qr,
          "    if (p2s4::inject(\"qr\") || !arma::qr_econ(Q, R, V))",
          "    {",
          "      p2s4::log(\"qr_failed\\n\");"
        )
      ),
      ed(
        "      AQ = arma::solve(arma::trimatl(R.t()), AQ.t()).t();",
        lines(
          "      p2s4::log(\"solve %llu %llu\\n\", p2s4::ull(r),",
          "                p2s4::ull(AQ.n_rows));",
          "      if (p2s4::inject(\"solve\"))",
          "      {",
          "        throw std::runtime_error(\"solve(): solution not found\");",
          "      }",
          singular,
          "      AQ = arma::solve(arma::trimatl(R.t()), AQ.t()).t();"
        )
      ),
      ed(aq_line, log_dropped),
      ed(
        "    w.set_size(n);\n    arma::blas_int info = 0, lwork = -1, liwork = -1;",
        lines(
          "    w.set_size(n);",
          "    if (p2s4::inject(\"syevdnan\"))",
          "    {",
          sprintf("      M(0, 0) = %s;", nan_cell),
          "    }",
          "    if (p2s4::inject(\"syevd\"))",
          "    {",
          "      p2s4::log(\"syevd_failed 1\\n\");",
          "      return false;",
          "    }",
          "    arma::blas_int info = 0, lwork = -1, liwork = -1;"
        )
      ),
      ed(
        "    if (info != 0)\n    {\n      return false;",
        lines(
          "    if (info != 0)",
          "    {",
          "      p2s4::log(\"syevd_failed %d\\n\", static_cast<int>(info));",
          "      return false;"
        ),
        n = 2L
      ),
      ed(
        "    if (work_buf.size() < static_cast<size_t>(work_q))",
        lines(
          "    p2s4::log(\"syevd %d query %llu held %llu\\n\", static_cast<int>(n),",
          "              p2s4::ull(work_q), p2s4::ull(work_buf.size()));",
          "    if (std::getenv(\"SLIDEIMP_P2S4_FRESH\"))",
          "    {",
          "      work_buf.clear();",
          "      iwork_buf.clear();",
          "    }",
          "    if (work_buf.size() < static_cast<size_t>(work_q))"
        )
      )
    )
  )
} else {
  patch_file(
    "src/lobpcg_warm.h",
    list(
      ed(
        "#include \"si_lapack.h\"",
        lines("#include \"si_lapack.h\"", harness)
      ),
      xonly,
      ed(
        "    residual(AX, X, lam, R);",
        lines(
          "    if (p2s4::inject(\"resnan\"))",
          "    {",
          "      for (std::size_t p2s4_i = 0; p2s4_i < AX.n_rows; ++p2s4_i)",
          "      {",
          sprintf("        AX(p2s4_i, 0) = %s;", nan_cell),
          "      }",
          "    }",
          "    residual(AX, X, lam, R);"
        )
      ),
      ed(rt_line, lines(singular, rt_line)),
      ed(
        "    if (status == si_la::ok)\n    {\n      return false;\n    }",
        lines(
          "    if (status == si_la::ok)",
          "    {",
          "      return false;",
          "    }",
          "    p2s4::log(\"status %d from %s\\n\", status, routine);"
        )
      ),
      ed(
        "    if (la_failed(si_la::qr_econ(v_rows, v_cols, V.memptr(), Q.memptr(),\n                                 R.memptr()),\n                  \"dgeqrf or dorgqr\"))\n    {",
        lines(
          log_qr,
          "    if (la_failed(p2s4::inject(\"qr\")",
          "                      ? 1",
          "                      : si_la::qr_econ(v_rows, v_cols, V.memptr(),",
          "                                       Q.memptr(), R.memptr()),",
          "                  \"dgeqrf or dorgqr\"))",
          "    {",
          "      p2s4::log(\"qr_failed\\n\");"
        )
      ),
      if (plain_solve) {
        ## Phase 4: plain dtrtrs, no estimate and no fallback. The line keeps
        ## the fields dev/probe_p2s4_summary.R reads, with rcond -1 standing
        ## for "not estimated".
        ed(
          "      if (la_failed(si_la::solve_tri('L', Rt.n_rows, Y.n_cols, Rt.memptr(),\n                                     Y.memptr()),\n                    \"dtrtrs\"))",
          lines(
            "      const int p2s4_st =",
            "          p2s4::inject(\"solve\")",
            "              ? 1",
            "              : si_la::solve_tri('L', Rt.n_rows, Y.n_cols, Rt.memptr(),",
            "                                 Y.memptr());",
            "      p2s4::log(\"solve %llu %llu rcond -1 approx 0\\n\", p2s4::ull(r),",
            "                p2s4::ull(AQ.n_rows));",
            "      if (la_failed(p2s4_st, \"dtrtrs\"))"
          )
        )
      } else {
        ed(
          "      if (la_failed(si_la::solve_tri('L', Rt.n_rows, Y.n_cols, Rt.memptr(),\n                                     Y.memptr(), &rcond, &approx),\n                    \"dtrtrs or dgelsd\"))",
          lines(
            "      const int p2s4_st =",
            "          p2s4::inject(\"solve\")",
            "              ? 1",
            "              : si_la::solve_tri('L', Rt.n_rows, Y.n_cols, Rt.memptr(),",
            "                                 Y.memptr(), &rcond, &approx);",
            "      p2s4::log(\"solve %llu %llu rcond %.6g approx %d\\n\", p2s4::ull(r),",
            "                p2s4::ull(AQ.n_rows), rcond, static_cast<int>(approx));",
            "      if (la_failed(p2s4_st, \"dtrtrs or dgelsd\"))"
          )
        )
      },
      ed(aq_line, log_dropped),
      ed(
        "    if (la_failed(si_la::syevd('U', n, M.memptr(), w.memptr()), \"dsyevd\"))\n    {",
        lines(
          "    if (p2s4::inject(\"syevdnan\"))",
          "    {",
          sprintf("      M(0, 0) = %s;", nan_cell),
          "    }",
          "    p2s4::log(\"syevd %d\\n\", static_cast<int>(n));",
          "    if (la_failed(p2s4::inject(\"syevd\")",
          "                      ? 1",
          "                      : si_la::syevd('U', n, M.memptr(), w.memptr()),",
          "                  \"dsyevd\"))",
          "    {",
          "      p2s4::log(\"syevd_failed\\n\");"
        )
      )
    )
  )
}
cat(
  "patched the",
  mode,
  "clone:",
  if (arma_calls) "calls through Armadillo" else "calls through si_la",
  if (gram_arma) "(Gram as arma::mat)" else "(Gram as a raw buffer)",
  "\n"
)
