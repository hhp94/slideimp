# Proves the Armadillo oracle can fail (dev/eigen_conversion.md, "The
# Armadillo oracle": "It must be shown to fail"). Dev only, run by hand, from
# the package root:
#   Rscript dev/oracle_selftest.R <lib> [--build=<dir>] [--only=<regex>]
# <lib>     the private RcppArmadillo library, as for dev/oracle_run.R.
# --build   where the broken copies are written and built (default: a folder
#           under tempdir()). Each run's output is kept there as
#           selftest-<name>.log.
# --only    run only the mutations whose name matches.
#
# For each mutation below it copies src/si_lapack.h and src/si_lapack.cpp to
# a folder of its own, breaks one wrapper in that copy with one textual
# replacement, and runs dev/oracle_run.R on the copy in a child R process.
# The oracle must then exit with status 1 and name every routine the
# mutation was aimed at. The working tree's src/ is never written to, so
# there is nothing to restore; the last run is the oracle on the unbroken
# source, which must pass.
#
# A replacement must change the text the stated number of times, or the
# mutation is reported as not applied: a self-test that silently broke
# nothing would pass for the wrong reason. Each mutation leaves a call that
# is still valid for BLAS and LAPACK, so the wrong answer comes back as a
# value and not as a crash; the two that drop a leading-dimension check are
# the exception: BLAS refuses the call through xerbla. Under one BLAS that
# is a printed line and a wrapper status of 0, under another an R error,
# which the oracle turns into every refusal failing; it fails either way.
#
# Wrong wrappers the oracle is known not to catch, and so not listed: one
# that asks dgesdd for its optimum below 1024 cells, one that gives dorgqr
# max(m, n) cells, and one without the eigenpair count check.
# dev/eigen_conversion.md, "The review of the wrapper layer", has the
# measurement.
#
# The last line is "SELFTEST PASSED" or "SELFTEST FAILED", and the exit
# status is 1 on failure.

flip <- "(uplo == 'U') ? \"L\" : \"U\""
mutations <- list(
  list(
    name = "gemm_beta",
    what = "gemm passes alpha where beta belongs",
    expect = "gemm",
    old = "b, &ldb_i, &beta, c, &ldc_i FCONE FCONE);",
    new = "b, &ldb_i, &alpha, c, &ldc_i FCONE FCONE);",
    times = 1L
  ),
  list(
    name = "gemv_beta",
    what = "gemv passes alpha where beta belongs",
    expect = "gemv",
    old = "&beta, y, &inc FCONE);",
    new = "&alpha, y, &inc FCONE);",
    times = 1L
  ),
  list(
    name = "syrk_triangle",
    what = "syrk writes the other triangle",
    expect = c("syrk", "syrk_block"),
    old = "F77_CALL(dsyrk)(&uplo, &trans,",
    new = paste0("F77_CALL(dsyrk)(", flip, ", &trans,"),
    times = 1L
  ),
  list(
    name = "matmul_one_row",
    what = "matmul negates a product that has one row",
    expect = "matmul",
    old = "return gemv('T', b_rows, b_cols, alpha, b, b_rows, a, beta, c);",
    new = "return gemv('T', b_rows, b_cols, -alpha, b, b_rows, a, beta, c);",
    times = 1L
  ),
  list(
    name = "syevr_triangle",
    what = "syevr_range reads the other triangle",
    expect = "syevr_range",
    old = "F77_CALL(dsyevr)(&jobz, &range, &uplo,",
    new = paste0("F77_CALL(dsyevr)(&jobz, &range, ", flip, ","),
    times = 2L
  ),
  list(
    name = "syevr_count",
    what = "syevr_range reports a failure LAPACK did not have",
    expect = "syevr_range",
    old = "if (m_out < 0 || static_cast<std::size_t>(m_out) != count)",
    new = "if (true)",
    times = 1L
  ),
  list(
    name = "syevd_triangle",
    what = "syevd reads the other triangle",
    expect = "syevd",
    old = "F77_CALL(dsyevd)(&jobz, &uplo,",
    new = paste0("F77_CALL(dsyevd)(&jobz, ", flip, ","),
    times = 2L
  ),
  list(
    name = "qr_r_below",
    what = "qr_econ leaves ones below the diagonal of r",
    expect = "qr_econ",
    old = "std::fill(r + j * n + j + 1, r + j * n + n, 0.0);",
    new = "std::fill(r + j * n + j + 1, r + j * n + n, 1.0);",
    times = 1L
  ),
  list(
    name = "solve_triangle",
    what = "solve_tri solves with the other triangle",
    expect = c("solve_tri", "solve_tri_chain"),
    old = "F77_CALL(dtrtrs)(&uplo, &trans, &diag,",
    new = paste0("F77_CALL(dtrtrs)(", flip, ", &trans, &diag,"),
    times = 1L
  ),
  list(
    name = "solve_approx_flag",
    what = "solve_tri does not say its answer is the approximate one",
    expect = "solve_tri",
    old = "*approx = true;",
    new = "*approx = false;",
    times = 1L
  ),
  list(
    name = "svd_sign",
    what = "svd_values flips the sign of the largest singular value",
    expect = "svd_values",
    old = paste(
      "iwork.data(), &info FCONE);",
      "    return info;",
      "  }",
      "",
      "} // namespace si_la",
      sep = "\n"
    ),
    new = paste(
      "iwork.data(), &info FCONE);",
      "    s[0] = -s[0];",
      "    return info;",
      "  }",
      "",
      "} // namespace si_la",
      sep = "\n"
    ),
    times = 1L
  ),
  list(
    name = "syevr_w_past",
    what = "syevr_range writes n eigenvalues into the caller's w",
    expect = "syevr_range",
    old = "std::copy(w_all.begin(), w_all.begin() + count, w);",
    new = "std::copy(w_all.begin(), w_all.end(), w);",
    times = 1L
  ),
  list(
    name = "syevr_abstol",
    what = "syevr_range passes abstol = 1e-12, not 0",
    expect = "syevr_chain",
    old = "const double abstol = 0.0;",
    new = "const double abstol = 1e-12;",
    times = 1L
  ),
  list(
    name = "work_query_ignored",
    what = "every workspace is the documented minimum, whatever LAPACK asked for",
    expect = "syevr_chain",
    old = "if (queried > len)",
    new = "if (false && queried > len)",
    times = 1L
  ),
  list(
    name = "syevr_no_finite",
    what = "syevr_range does not refuse NaN or Inf",
    expect = c("syevr_range", "refusals"),
    old = "NaN in the result.\n    if (!tri_finite(uplo, n, a))",
    new = "NaN in the result.\n    if (false)",
    times = 1L
  ),
  list(
    name = "syevd_no_finite",
    what = "syevd does not refuse NaN or Inf",
    expect = "refusals",
    old = "// the result.\n    if (!tri_finite(uplo, n, a))",
    new = "// the result.\n    if (false)",
    times = 1L
  ),
  list(
    name = "finite_whole",
    what = "the eigen wrappers refuse NaN in the triangle LAPACK does not read",
    expect = c("syevr_range", "syevd", "refusals"),
    old = "if (!all_finite(a + j * n + first, len))",
    new = "if (!all_finite(a + j * n, n))",
    times = 1L
  ),
  list(
    name = "matmul_order",
    what = "matmul tests one column before one row",
    expect = "matmul",
    old = "if (c_rows == 1)",
    new = "if (c_rows == 1 && b_cols != 1)",
    times = 1L
  ),
  list(
    name = "matmul_column_gemm",
    what = "matmul sends a product with one column to dgemm",
    expect = "matmul",
    old = "if (b_cols == 1)",
    new = "if (false && b_cols == 1)",
    times = 1L
  ),
  list(
    name = "solve_threshold_up",
    what = "solve_tri falls back below 1000 * DBL_EPSILON",
    expect = "solve_tri",
    old = "solved = !(*rcond < eps) && !std::isnan(*rcond);",
    new = "solved = !(*rcond < 1000.0 * eps) && !std::isnan(*rcond);",
    times = 1L
  ),
  list(
    name = "solve_threshold_down",
    what = "solve_tri falls back below DBL_EPSILON / 1000",
    expect = "solve_tri",
    old = "solved = !(*rcond < eps) && !std::isnan(*rcond);",
    new = "solved = !(*rcond < eps / 1000.0) && !std::isnan(*rcond);",
    times = 1L
  ),
  list(
    name = "solve_no_nan_test",
    what = "solve_tri does not test the condition estimate for NaN",
    expect = "solve_tri",
    old = "solved = !(*rcond < eps) && !std::isnan(*rcond);",
    new = "solved = !(*rcond < eps);",
    times = 1L
  ),
  list(
    name = "solve_sv_cut",
    what = "solve_tri's fallback cuts singular values at DBL_EPSILON",
    expect = "solve_tri",
    old = "double rcond_sv = nd * eps;",
    new = "double rcond_sv = eps;",
    times = 1L
  ),
  list(
    name = "solve_code",
    what = "solve_tri returns err_arg where it refuses NaN or Inf",
    expect = "solve_tri",
    old = paste(
      "!all_finite(rhs.data(), n * nrhs))",
      "    {",
      "      return err_nonfinite;",
      sep = "\n"
    ),
    new = paste(
      "!all_finite(rhs.data(), n * nrhs))",
      "    {",
      "      return err_arg;",
      sep = "\n"
    ),
    times = 1L
  ),
  list(
    name = "svd_code",
    what = "svd_values returns err_arg where it refuses NaN or Inf",
    expect = "svd_values",
    old = paste(
      "if (!all_finite(a, m * n))",
      "    {",
      "      return err_nonfinite;",
      sep = "\n"
    ),
    new = paste(
      "if (!all_finite(a, m * n))",
      "    {",
      "      return err_arg;",
      sep = "\n"
    ),
    times = 1L
  ),
  list(
    name = "gemm_lda_trans",
    what = "gemm does not check lda for a transposed operand",
    expect = "refusals",
    old = "const std::size_t a_rows = (trans_a == 'N') ? m : k;",
    new = "const std::size_t a_rows = (trans_a == 'N') ? m : 1;",
    times = 1L
  ),
  list(
    name = "syrk_lda_trans",
    what = "syrk does not check lda for trans 'T'",
    expect = "refusals",
    old = "const std::size_t a_rows = (trans == 'N') ? n : k;",
    new = "const std::size_t a_rows = (trans == 'N') ? n : 1;",
    times = 1L
  ),
  list(
    name = "trans_conj",
    what = "the wrappers take 'C' as a transpose option",
    expect = "refusals",
    old = "return c == 'N' || c == 'T';",
    new = "return c == 'N' || c == 'T' || c == 'C';",
    times = 1L
  )
)

args <- commandArgs(trailingOnly = TRUE)
opts <- grep("^--", args, value = TRUE)
pos <- setdiff(args, opts)
if (length(pos) != 1L) {
  stop(
    "usage: Rscript dev/oracle_selftest.R <lib> [--build=<dir>] ",
    "[--only=<regex>]",
    call. = FALSE
  )
}
opt_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), opts, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
lib <- normalizePath(pos[[1]], winslash = "/", mustWork = TRUE)
build <- opt_value("build", file.path(tempdir(), "oracle-selftest"))
only <- opt_value("only", NULL)
if (!file.exists("dev/oracle_run.R") || !file.exists("src/si_lapack.cpp")) {
  stop("run from the package root", call. = FALSE)
}
dir.create(build, recursive = TRUE, showWarnings = FALSE)
build <- normalizePath(build, winslash = "/", mustWork = TRUE)
if (!is.null(only)) {
  keep <- grepl(only, vapply(mutations, `[[`, "", "name"))
  mutations <- mutations[keep]
}
if (!length(mutations)) {
  stop("no mutation selected", call. = FALSE)
}

rscript <- file.path(R.home("bin"), "Rscript")
source_text <- paste(readLines("src/si_lapack.cpp"), collapse = "\n")

## Count the non-overlapping occurrences of a fixed string.
occurrences <- function(pattern, text) {
  hits <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (identical(as.integer(hits), -1L)) 0L else length(hits)
}

## Run the oracle on the wrapper source in `src_dir`; return its exit status
## and the routines its last line names.
run_oracle <- function(name, src_dir) {
  out <- suppressWarnings(system2(
    rscript,
    c(
      "dev/oracle_run.R",
      shQuote(lib),
      shQuote(paste0("--src=", src_dir)),
      shQuote(paste0("--build=", file.path(build, paste0("build-", name))))
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  writeLines(out, file.path(build, paste0("selftest-", name, ".log")))
  verdict <- grep("^ORACLE (PASSED|FAILED)", out, value = TRUE)
  failed <- character()
  if (length(verdict) == 1L && startsWith(verdict, "ORACLE FAILED:")) {
    failed <- trimws(strsplit(sub("^ORACLE FAILED:", "", verdict), ",")[[1]])
  }
  list(
    status = status,
    verdict = if (length(verdict) == 1L) verdict else "(no verdict line)",
    failed = failed
  )
}

all_ok <- TRUE
cat(sprintf("%d mutations, each built under %s\n\n", length(mutations), build))
for (m in mutations) {
  n_found <- occurrences(m$old, source_text)
  if (n_found != m$times) {
    cat(sprintf(
      "%-18s NOT APPLIED: the text to replace occurs %d times, expected %d\n",
      m$name,
      n_found,
      m$times
    ))
    all_ok <- FALSE
    next
  }
  broken <- gsub(m$old, m$new, source_text, fixed = TRUE)
  src_dir <- file.path(build, paste0("src-", m$name))
  unlink(src_dir, recursive = TRUE)
  dir.create(src_dir)
  file.copy("src/si_lapack.h", src_dir)
  writeLines(broken, file.path(src_dir, "si_lapack.cpp"))

  res <- run_oracle(m$name, src_dir)
  caught <- res$status == 1L && all(m$expect %in% res$failed)
  all_ok <- all_ok && caught
  cat(sprintf(
    "%-18s %s\n  broken: %s\n  aimed at: %s\n  oracle exit %d, %s\n",
    m$name,
    if (caught) "caught" else "MISSED",
    m$what,
    paste(m$expect, collapse = ", "),
    res$status,
    res$verdict
  ))
}

cat("\nthe unbroken source:\n")
res <- run_oracle("unbroken", normalizePath("src", winslash = "/"))
clean <- res$status == 0L && identical(res$verdict, "ORACLE PASSED")
all_ok <- all_ok && clean
cat(sprintf(
  "%-18s %s\n  oracle exit %d, %s\n",
  "unbroken",
  if (clean) "passes" else "DOES NOT PASS",
  res$status,
  res$verdict
))

cat("\n")
if (!all_ok) {
  cat("SELFTEST FAILED\n")
  quit(status = 1L)
}
cat("SELFTEST PASSED\n")
