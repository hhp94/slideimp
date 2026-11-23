# Proves the test of the wrapper layer can fail (dev/eigen_conversion.md at
# 4f2b669, "Rulings on the last items"). Dev only, run by hand, from the package
# root:
#   Rscript dev/tool/si_lapack_selftest.R [--build=<dir>] [--only=<regex>]
# --build   where the broken copies are written and built (default: a folder
#           under tempdir()). Each run's output is kept there as
#           selftest-<name>.log.
# --only    run only the mutations whose name matches.
#
# For each mutation below it copies src/si_lapack.h and src/si_lapack.cpp to
# a folder of its own, breaks one wrapper in that copy with textual
# replacements, and runs dev/tool/si_lapack_test.R on the copy in a child R
# process. The test must then exit with status 1 and name every routine the
# mutation is aimed at. Where a mutation gives `check`, a failure line of
# each routine aimed at must name that check; with `only`, every failure line
# of those routines must name it, so the check is shown to be what catches
# the mutation. The working tree's src/ is never written to, so there is
# nothing to restore; the last run is the test on the unbroken source, which
# must pass.
#
# A replacement must change the text the stated number of times, or the
# mutation is reported as not applied: a self-test that silently broke
# nothing would pass for the wrong reason. Most mutations leave a call that
# is still valid for BLAS and LAPACK, so the wrong answer comes back as a
# value. The four that drop an argument check (gemm_ldc, gemm_lda_trans,
# qr_wide, syevr_il_iu) pass BLAS or LAPACK an illegal argument instead.
# Measured with R's BLAS and LAPACK replaced by OpenBLAS 0.3.34: the routine
# returns at once, so the wrapper reports 0 after dgemm and LAPACK's
# negative info after dorgqr and dsyevr, and the case fails on its status.
# Under R's own BLAS and LAPACK xerbla raises an R error instead (not run
# here), which the test records as a failing case.
#
# The first mutation is the one the ruling names: syevr_range hands dsyevr
# the caller's eigenvalue buffer, as the chain did before step 3 of Phase 2.
# Only the guard cells after w can see it: the first iu - il + 1 cells hold
# the right eigenvalues.
#
# Wrong wrappers not listed, and why. A size check that lets a size past int
# through: the test's size cases hand the wrapper a buffer of one cell, so
# the broken copy reads or writes far past it and the child crashes rather
# than failing. A workspace that ignores LAPACK's answer, or a nonzero
# abstol: the eigenpairs stay right to the test's tolerance, and this test
# compares with R's results, not with the chain's bits.
#
# The last line is "SELFTEST PASSED" or "SELFTEST FAILED", and the exit
# status is 1 on failure.

flip <- "(uplo == 'U') ? \"L\" : \"U\""
edit <- function(old, new, times = 1L) {
  list(old = old, new = new, times = times)
}
mutations <- list(
  list(
    name = "syevr_caller_w",
    what = "syevr_range hands dsyevr the caller's w and copies nothing",
    expect = "syevr_range",
    check = "the cells after w are as they were",
    only = TRUE,
    edits = list(
      edit(
        "&abstol, &m_out, w_all.data(), z, &n_i,",
        "&abstol, &m_out, w, z, &n_i,",
        2L
      ),
      edit("  std::copy(w_all.begin(), w_all.begin() + count, w);\n", "")
    )
  ),
  list(
    name = "syevr_triangle",
    what = "syevr_range reads the other triangle",
    expect = "syevr_range",
    edits = list(edit(
      "F77_CALL(dsyevr)(&jobz, &range, &uplo,",
      paste0("F77_CALL(dsyevr)(&jobz, &range, ", flip, ","),
      2L
    ))
  ),
  list(
    name = "syevr_no_finite",
    what = "syevr_range does not refuse NaN or Inf",
    expect = "syevr_range",
    edits = list(edit(
      "with NaN in the result.\n  if (!tri_finite(uplo, n, a))",
      "with NaN in the result.\n  if (false)"
    ))
  ),
  list(
    name = "syevr_copy_short",
    what = "syevr_range copies one eigenvalue too few",
    expect = "syevr_range",
    edits = list(edit(
      "std::copy(w_all.begin(), w_all.begin() + count, w);",
      "std::copy(w_all.begin(), w_all.begin() + count - 1, w);"
    ))
  ),
  list(
    name = "syevr_il_iu",
    what = "syevr_range does not refuse il > iu",
    expect = "syevr_range",
    edits = list(edit(
      "if (il < 1 || il > iu || iu > n)",
      "if (il < 1 || iu > n)"
    ))
  ),
  list(
    name = "syevd_triangle",
    what = "syevd reads the other triangle",
    expect = "syevd",
    edits = list(edit(
      "F77_CALL(dsyevd)(&jobz, &uplo,",
      paste0("F77_CALL(dsyevd)(&jobz, ", flip, ","),
      2L
    ))
  ),
  list(
    name = "syevd_no_finite",
    what = "syevd does not refuse NaN or Inf",
    expect = "syevd",
    edits = list(edit(
      "// the result.\n  if (!tri_finite(uplo, n, a))",
      "// the result.\n  if (false)"
    ))
  ),
  list(
    name = "finite_whole",
    what = "the eigen wrappers refuse NaN in the triangle LAPACK does not read",
    expect = c("syevr_range", "syevd"),
    edits = list(edit(
      "if (!all_finite(a + j * n + first, len))",
      "if (!all_finite(a + j * n, n))"
    ))
  ),
  list(
    name = "gemm_beta",
    what = "gemm passes alpha where beta belongs",
    expect = "gemm",
    edits = list(edit(
      "&ldb_i, &beta, c, &ldc_i FCONE FCONE);",
      "&ldb_i, &alpha, c, &ldc_i FCONE FCONE);"
    ))
  ),
  list(
    name = "gemm_k0",
    what = "gemm leaves C alone when the inner dimension is 0",
    expect = "gemm",
    edits = list(edit("scale_range(c + j * ldc, m, beta);", "(void)beta;"))
  ),
  list(
    name = "gemm_ldc",
    what = "gemm does not check ldc",
    expect = "gemm",
    edits = list(edit("if (ldc < m)", "if (false)"))
  ),
  list(
    name = "gemm_lda_trans",
    what = "gemm does not check lda for a transposed operand",
    expect = "gemm",
    edits = list(edit(
      "const std::size_t a_rows = (trans_a == 'N') ? m : k;",
      "const std::size_t a_rows = (trans_a == 'N') ? m : 1;"
    ))
  ),
  list(
    name = "gemv_beta",
    what = "gemv passes alpha where beta belongs",
    expect = "gemv",
    edits = list(edit("x, &inc, &beta, y,", "x, &inc, &alpha, y,"))
  ),
  list(
    name = "gemv_x_empty",
    what = "gemv leaves y alone when x is empty",
    expect = "gemv",
    edits = list(edit("scale_range(y, len_y, beta);", "(void)beta;"))
  ),
  list(
    name = "syrk_triangle",
    what = "syrk writes the other triangle",
    expect = "syrk",
    edits = list(edit(
      "F77_CALL(dsyrk)(&uplo, &trans,",
      paste0("F77_CALL(dsyrk)(", flip, ", &trans,")
    ))
  ),
  list(
    name = "syrk_k0_column",
    what = "syrk scales whole columns of C when k is 0",
    expect = "syrk",
    check = "the other triangle of C is as it was",
    edits = list(edit(
      "scale_range(c + j * ldc, j + 1, beta);",
      "scale_range(c + j * ldc, n, beta);"
    ))
  ),
  list(
    name = "qr_r_below",
    what = "qr_econ leaves ones below the diagonal of r",
    expect = "qr_econ",
    edits = list(edit(
      "std::fill(r + j * n + j + 1, r + j * n + n, 0.0);",
      "std::fill(r + j * n + j + 1, r + j * n + n, 1.0);"
    ))
  ),
  list(
    name = "qr_wide",
    what = "qr_econ does not refuse m < n",
    expect = "qr_econ",
    edits = list(edit("if (m < n)", "if (false)"))
  ),
  list(
    name = "solve_triangle",
    what = "solve_tri_right solves with the other triangle",
    expect = "solve_tri_right",
    edits = list(edit(
      "F77_CALL(dtrsm)(&side, &uplo, &trans, &diag,",
      paste0("F77_CALL(dtrsm)(&side, ", flip, ", &trans, &diag,")
    ))
  ),
  list(
    name = "solve_trans",
    what = "solve_tri_right solves with the transpose",
    expect = "solve_tri_right",
    edits = list(edit("const char trans = 'N';", "const char trans = 'T';"))
  ),
  list(
    name = "solve_info",
    what = "solve_tri_right reports success on a zero pivot",
    expect = "solve_tri_right",
    edits = list(edit("return static_cast<int>(i + 1);", "return ok;"))
  ),
  list(
    name = "svd_sign",
    what = "svd_values flips the sign of the largest singular value",
    expect = "svd_values",
    edits = list(edit(
      "iwork.data(), &info FCONE);\n  return info;",
      "iwork.data(), &info FCONE);\n  s[0] = -s[0];\n  return info;"
    ))
  ),
  list(
    name = "svd_code",
    what = "svd_values returns err_arg where it refuses NaN or Inf",
    expect = "svd_values",
    edits = list(edit(
      "if (!all_finite(a, m * n)) {\n    return err_nonfinite;",
      "if (!all_finite(a, m * n)) {\n    return err_arg;"
    ))
  ),
  list(
    name = "trans_conj",
    what = "the wrappers take 'C' as a transpose option",
    expect = c("gemm", "gemv", "syrk"),
    edits = list(edit(
      "return c == 'N' || c == 'T';",
      "return c == 'N' || c == 'T' || c == 'C';"
    ))
  )
)

args <- commandArgs(trailingOnly = TRUE)
unknown <- args[!grepl("^--(build|only)=.+$", args)]
if (length(unknown)) {
  stop(
    "usage: Rscript dev/tool/si_lapack_selftest.R [--build=<dir>] [--only=<regex>]",
    call. = FALSE
  )
}
opt_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
build <- opt_value("build", file.path(tempdir(), "si-lapack-selftest"))
only <- opt_value("only", NULL)
if (
  !file.exists("dev/tool/si_lapack_test.R") || !file.exists("src/si_lapack.cpp")
) {
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

## Run the test on the wrapper source in `src_dir`; return its exit status,
## its verdict line, the routines that line names, and its failure lines.
run_test <- function(name, src_dir) {
  out <- suppressWarnings(system2(
    rscript,
    c(
      "dev/tool/si_lapack_test.R",
      shQuote(paste0("--src=", src_dir)),
      shQuote(paste0("--build=", file.path(build, paste0("build-", name))))
    ),
    stdout = TRUE,
    stderr = TRUE,
    timeout = 600
  ))
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  writeLines(out, file.path(build, paste0("selftest-", name, ".log")))
  verdict <- grep("^WRAPPER TEST (PASSED|FAILED)", out, value = TRUE)
  failed <- character()
  if (length(verdict) == 1L && startsWith(verdict, "WRAPPER TEST FAILED:")) {
    failed <- trimws(strsplit(
      sub("^WRAPPER TEST FAILED:", "", verdict),
      ","
    )[[1]])
  }
  list(
    status = status,
    verdict = if (length(verdict) == 1L) verdict else "(no verdict line)",
    failed = failed,
    fail_lines = grep("^FAIL ", out, value = TRUE)
  )
}

## Whether the failure lines of `routine` name `check`: some of them, or,
## with `only`, every one.
names_check <- function(fail_lines, routine, check, only) {
  mine <- fail_lines[startsWith(fail_lines, paste0("FAIL ", routine, " | "))]
  hit <- grepl(check, mine, fixed = TRUE)
  length(mine) > 0L && (if (isTRUE(only)) all(hit) else any(hit))
}

all_ok <- TRUE
cat(sprintf("%d mutations, each built under %s\n\n", length(mutations), build))
for (m in mutations) {
  broken <- source_text
  applied <- TRUE
  for (e in m$edits) {
    n_found <- occurrences(e$old, broken)
    if (n_found != e$times) {
      cat(sprintf(
        "%-18s NOT APPLIED: a text to replace occurs %d times, expected %d\n",
        m$name,
        n_found,
        e$times
      ))
      applied <- FALSE
      break
    }
    broken <- gsub(e$old, e$new, broken, fixed = TRUE)
  }
  if (!applied) {
    all_ok <- FALSE
    next
  }
  src_dir <- file.path(build, paste0("src-", m$name))
  unlink(src_dir, recursive = TRUE)
  dir.create(src_dir)
  file.copy("src/si_lapack.h", src_dir)
  writeLines(broken, file.path(src_dir, "si_lapack.cpp"))

  res <- run_test(m$name, src_dir)
  caught <- res$status == 1L && all(m$expect %in% res$failed)
  by_check <- TRUE
  if (!is.null(m$check)) {
    by_check <- all(vapply(
      m$expect,
      function(r) names_check(res$fail_lines, r, m$check, m$only),
      TRUE
    ))
  }
  ok <- caught && by_check
  all_ok <- all_ok && ok
  cat(sprintf(
    "%-18s %s\n  broken: %s\n  aimed at: %s\n  test exit %d, %s\n",
    m$name,
    if (ok) "caught" else "MISSED",
    m$what,
    paste(m$expect, collapse = ", "),
    res$status,
    res$verdict
  ))
  if (!is.null(m$check)) {
    cat(sprintf(
      "  %s failure line names \"%s\": %s\n",
      if (isTRUE(m$only)) "every" else "a",
      m$check,
      if (by_check) "yes" else "NO"
    ))
  }
  cat(sprintf("  %d failure lines\n", length(res$fail_lines)))
}

cat("\nthe unbroken source:\n")
res <- run_test("unbroken", normalizePath("src", winslash = "/"))
clean <- res$status == 0L && identical(res$verdict, "WRAPPER TEST PASSED")
all_ok <- all_ok && clean
cat(sprintf(
  "%-18s %s\n  test exit %d, %s\n",
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
