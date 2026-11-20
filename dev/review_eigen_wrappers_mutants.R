# Review of the wrapper layer and the oracle (dev/review_eigen_wrappers.md):
# can dev/oracle_run.R pass with a wrong wrapper that dev/oracle_selftest.R
# does not try? The oracle and dev/lib/arma were removed after 141a9d0; the
# commit under review, 79cce0f, has the oracle's files (git show
# 79cce0f:<path>). Dev only, run by hand, from the package root:
#   Rscript dev/review_eigen_wrappers_mutants.R <lib> <out dir> [--only=<regex>]
# <lib>      the private RcppArmadillo library, as for dev/oracle_run.R.
# <out dir>  where the broken copies are written and built, and each run's
#            output is kept as mutant-<name>.log.
#
# Same mechanism as dev/oracle_selftest.R: src/si_lapack.h and
# src/si_lapack.cpp are copied to a folder of their own, the copy is changed
# by textual replacement, and dev/oracle_run.R runs on the copy in a child R
# process. src/ is never written to. Each replacement must change the text
# the stated number of times, or the mutant is reported as not applied.
#
# The first run is the unbroken source. For every later run the script
# prints the oracle's exit status and verdict, and each routine whose count
# of bitwise-identical cases differs from the unbroken run: the oracle
# reports that count and never asserts on it, so a mutant can move it and
# still pass. Nothing here asserts; the last line counts the survivors.

flip <- "(uplo == 'U') ? \"L\" : \"U\""
rep1 <- function(old, new, times = 1L) list(old = old, new = new, times = times)
mutants <- list(
  list(
    name = "work_min_only",
    what = "every workspace is the documented minimum; LAPACK's answer to the query is ignored",
    reps = list(rep1("if (queried > len)", "if (false && queried > len)"))
  ),
  list(
    name = "syevr_abstol",
    what = "syevr_range passes abstol = 1e-12, not 0",
    reps = list(rep1(
      "const double abstol = 0.0;",
      "const double abstol = 1e-12;"
    ))
  ),
  list(
    name = "syevr_no_count",
    what = "syevr_range does not check how many eigenpairs came back",
    reps = list(rep1(
      "if (m_out < 0 || static_cast<std::size_t>(m_out) != count)",
      "if (false)"
    ))
  ),
  list(
    name = "syevr_isuppz_half",
    what = "syevr_range allocates isuppz at half the documented length",
    reps = list(rep1("alloc(isuppz, 2 * count)", "alloc(isuppz, count)"))
  ),
  list(
    name = "matmul_order",
    what = "matmul tests one column before one row (a 1 x 1 result goes to the other dgemv)",
    reps = list(rep1("if (c_rows == 1)", "if (c_rows == 1 && b_cols != 1)"))
  ),
  list(
    name = "matmul_never_gemv",
    what = "matmul sends every product to dgemm, never dgemv",
    reps = list(
      rep1("if (c_rows == 1)", "if (false && c_rows == 1)"),
      rep1("if (b_cols == 1)", "if (false && b_cols == 1)")
    )
  ),
  list(
    name = "svd_always_query",
    what = "svd_values asks dgesdd for its optimum below 1024 cells too",
    reps = list(rep1("if (m * n >= 1024)", "if (true)"))
  ),
  list(
    name = "svd_no_finite_check",
    what = "svd_values does not refuse NaN or Inf; dgesdd is entered with it",
    reps = list(rep1(
      "if (!all_finite(a, m * n))",
      "if (false && !all_finite(a, m * n))"
    ))
  ),
  list(
    name = "solve_code",
    what = "solve_tri returns err_arg, not err_nonfinite, from its fallback",
    reps = list(rep1(
      paste(
        "if (!all_finite(tri.data(), n * n) || !all_finite(rhs.data(), n * nrhs))",
        "    {",
        "      return err_nonfinite;",
        sep = "\n"
      ),
      paste(
        "if (!all_finite(tri.data(), n * n) || !all_finite(rhs.data(), n * nrhs))",
        "    {",
        "      return err_arg;",
        sep = "\n"
      )
    ))
  ),
  list(
    name = "solve_threshold_up",
    what = "solve_tri falls back below 1000 * DBL_EPSILON, not DBL_EPSILON",
    reps = list(rep1(
      "solved = !(*rcond < eps) && !std::isnan(*rcond);",
      "solved = !(*rcond < 1000.0 * eps) && !std::isnan(*rcond);"
    ))
  ),
  list(
    name = "solve_threshold_down",
    what = "solve_tri falls back below DBL_EPSILON / 1000, not DBL_EPSILON",
    reps = list(rep1(
      "solved = !(*rcond < eps) && !std::isnan(*rcond);",
      "solved = !(*rcond < eps / 1000.0) && !std::isnan(*rcond);"
    ))
  ),
  list(
    name = "solve_no_nan_test",
    what = "solve_tri does not test the estimate for NaN",
    reps = list(rep1(
      "solved = !(*rcond < eps) && !std::isnan(*rcond);",
      "solved = !(*rcond < eps);"
    ))
  ),
  list(
    name = "solve_sv_cut",
    what = "solve_tri's fallback cuts singular values at DBL_EPSILON, not n * DBL_EPSILON",
    reps = list(rep1("double rcond_sv = nd * eps;", "double rcond_sv = eps;"))
  ),
  list(
    name = "gemm_lda_trans",
    what = "gemm does not check lda or ldb when the operand is transposed",
    reps = list(
      rep1(
        "const std::size_t a_rows = (trans_a == 'N') ? m : k;",
        "const std::size_t a_rows = (trans_a == 'N') ? m : 1;"
      ),
      rep1(
        "const std::size_t b_rows = (trans_b == 'N') ? k : n;",
        "const std::size_t b_rows = (trans_b == 'N') ? k : 1;"
      )
    )
  ),
  list(
    name = "syrk_lda_trans",
    what = "syrk does not check lda when trans is 'T'",
    reps = list(rep1(
      "const std::size_t a_rows = (trans == 'N') ? n : k;",
      "const std::size_t a_rows = (trans == 'N') ? n : 1;"
    ))
  ),
  list(
    name = "qr_dorgqr_query",
    what = "qr_econ passes dorgqr max(m, n) as its workspace, not dgeqrf's optimum",
    reps = list(rep1(
      paste(
        "F77_CALL(dorgqr)(&m_i, &p_i, &p_i, q, &m_i, tau.data(), work.data(),",
        "                     &lwork, &info);",
        sep = "\n"
      ),
      paste(
        "int lwork_small = std::max(m_i, n_i);",
        "    F77_CALL(dorgqr)(&m_i, &p_i, &p_i, q, &m_i, tau.data(), work.data(),",
        "                     &lwork_small, &info);",
        sep = "\n"
      )
    ))
  )
)

args <- commandArgs(trailingOnly = TRUE)
opts <- grep("^--", args, value = TRUE)
pos <- setdiff(args, opts)
if (length(pos) != 2L) {
  stop(
    "usage: Rscript dev/review_eigen_wrappers_mutants.R <lib> <out dir> ",
    "[--only=<regex>]",
    call. = FALSE
  )
}
lib <- normalizePath(pos[[1]], winslash = "/", mustWork = TRUE)
build <- pos[[2]]
dir.create(build, recursive = TRUE, showWarnings = FALSE)
build <- normalizePath(build, winslash = "/", mustWork = TRUE)
only <- sub("^--only=", "", grep("^--only=", opts, value = TRUE))
if (length(only)) {
  keep <- grepl(only[[1]], vapply(mutants, `[[`, "", "name"))
  mutants <- mutants[keep]
}
if (!file.exists("dev/oracle_run.R") || !file.exists("src/si_lapack.cpp")) {
  stop("run from the package root", call. = FALSE)
}
rscript <- file.path(R.home("bin"), "Rscript")
source_text <- paste(readLines("src/si_lapack.cpp"), collapse = "\n")

occurrences <- function(pattern, text) {
  hits <- gregexpr(pattern, text, fixed = TRUE)[[1]]
  if (identical(as.integer(hits), -1L)) 0L else length(hits)
}

## Per routine, "x of y" from the oracle's "bitwise identical in" lines.
bitwise_counts <- function(out) {
  head_at <- grep("^[a-z_]+ +(ok  |FAIL)  class ", out)
  res <- character()
  for (i in head_at) {
    nm <- sub(" .*$", "", out[[i]])
    blk <- out[i:min(length(out), i + 5L)]
    hit <- grep("bitwise identical in", blk, value = TRUE)
    if (length(hit)) {
      res[[nm]] <- sub("^.*identical in ([0-9]+ of [0-9]+).*$", "\\1", hit[[1]])
    }
  }
  res
}

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
  writeLines(
    c(out, sprintf("exit: %d", status)),
    file.path(build, paste0("mutant-", name, ".log"))
  )
  verdict <- grep("^ORACLE (PASSED|FAILED)", out, value = TRUE)
  list(
    status = status,
    verdict = if (length(verdict) == 1L) verdict else "(no verdict line)",
    bitwise = bitwise_counts(out),
    fails = grep("^  FAIL \\[", out, value = TRUE)
  )
}

clean <- run_oracle("unbroken", normalizePath("src", winslash = "/"))
cat(sprintf(
  "unbroken source: oracle exit %d, %s\n  bitwise: %s\n\n",
  clean$status,
  clean$verdict,
  paste(names(clean$bitwise), clean$bitwise, sep = " ", collapse = "; ")
))

survivors <- character()
for (m in mutants) {
  broken <- source_text
  applied <- TRUE
  for (r in m$reps) {
    n_found <- occurrences(r$old, broken)
    if (n_found != r$times) {
      cat(sprintf(
        "%-22s NOT APPLIED: the text occurs %d times, expected %d\n",
        m$name,
        n_found,
        r$times
      ))
      applied <- FALSE
      break
    }
    broken <- gsub(r$old, r$new, broken, fixed = TRUE)
  }
  if (!applied) {
    next
  }
  src_dir <- file.path(build, paste0("src-", m$name))
  unlink(src_dir, recursive = TRUE)
  dir.create(src_dir)
  file.copy("src/si_lapack.h", src_dir)
  writeLines(broken, file.path(src_dir, "si_lapack.cpp"))

  res <- run_oracle(m$name, src_dir)
  survived <- res$status == 0L && identical(res$verdict, "ORACLE PASSED")
  if (survived) {
    survivors <- c(survivors, m$name)
  }
  moved <- setdiff(
    paste(names(res$bitwise), res$bitwise),
    paste(names(clean$bitwise), clean$bitwise)
  )
  cat(sprintf(
    "%-22s %s\n  broken: %s\n  oracle exit %d, %s\n  bitwise counts that moved: %s\n",
    m$name,
    if (survived) "SURVIVED" else "caught",
    m$what,
    res$status,
    res$verdict,
    if (length(moved)) paste(moved, collapse = "; ") else "none"
  ))
  for (f in utils::head(res$fails, 4L)) {
    cat("  ", f, "\n", sep = "")
  }
}

cat(sprintf(
  "\n%d of %d mutants survived: %s\n",
  length(survivors),
  length(mutants),
  paste(survivors, collapse = ", ")
))
