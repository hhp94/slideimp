# The test of the wrapper layer's own (dev/eigen_conversion.md at 4f2b669,
# "Rulings on the last items"): every function of src/si_lapack.h, called
# through dev/si_lapack_test.cpp, checked against what its comment in the header
# says it does. The references are R's own; no other library is loaded. Dev
# only, run by hand, from the package root:
#   Rscript dev/si_lapack_test.R [--src=<dir>] [--build=<dir>] [--verbose]
# --src      the folder holding si_lapack.h and si_lapack.cpp (default:
#            src). dev/si_lapack_selftest.R points this at a copy with one
#            wrapper broken.
# --build    where the sources are staged and compiled (default: a folder
#            under tempdir()). Nothing is built in the working tree.
# --verbose  print every case, not only the failures and the worst per
#            routine.
#
# Both units are compiled with -DPCA_IMP_DIAGNOSTICS=1, so svd_values, which
# only the diagnostics build has, is tested with the rest.
#
# The references: products from R's own loops (options(matprod =
# "internal"), which calls no BLAS); eigen(symmetric = TRUE), which asks
# dsyevr for every eigenpair where the wrapper asks for a range; svd(), for
# the values of dgesdd; qr(), which is LINPACK's; backsolve() and
# forwardsolve().
#
# Per case:
# - the status, which every case states: a code of the wrapper's, or the
#   info LAPACK returns;
# - on success, each value against its reference: over the finite cells,
#   max(abs(got - ref)) <= tol$abs * s, where s is max(1, the largest
#   magnitude either side holds), and max(abs(got - ref) / m) <= tol$rel,
#   where m is max(abs(got), abs(ref)) of the cell, over the cells with
#   m > tol$rel_floor * s; the report gives how many cells that was. NaN and
#   Inf must sit in the same cells on both sides. Bitwise identity is
#   reported, not asserted: the references are other code.
# - with identical(): 16 guard cells holding a sentinel after every buffer
#   the wrapper writes, and n more after the eigenvalues and eigenvectors of
#   syevr_range; the rows past the last of an output whose leading
#   dimension is larger; the triangle syrk does not name; and every output
#   of a call that is refused or writes nothing.
# - checks of the routine's own, each named in the report: eigenvectors by
#   their residual and orthonormality, zeros below the diagonal of r, a
#   right-hand side left as it was on a zero pivot, and others.
# - a case that ends in an R error (xerbla, say) fails, and the run goes on.
#
# The guard cells after syevr_range's eigenvalues test something only on a
# range that makes dsyevr write past iu - il + 1 cells, which it does when
# the range cuts a group of equal eigenvalues. probe_dsyevr_past() calls
# dsyevr on every such case with a W of n cells and counts what it wrote
# past; the run fails if no case makes it write any.
#
# Tolerances: one set for every routine. The errors this test exists to
# find, an argument in the wrong place, are of the order of the values.
# tol$abs and tol$vec are about a hundredfold above the worst reading.
# tol$rel is set by its floor, not by a reading: one unit in the last place
# of the scale is 2.2e-10 relative in a cell at the floor, so a bound near
# the reading would test the floor and the BLAS, not the wrapper. Readings,
# 2026-11-19, R 4.6.2 with R's BLAS and LAPACK replaced by OpenBLAS 0.3.34:
# worst abs / scale 5.8e-15 (syevr_range, order 300), worst rel 6.7e-13
# (gemm, 40 x 30 over 1200 cells), worst vector check 8.3e-15 (a residual
# of syevr_range).
#
# Not reached by any case: err_count, which no known input gets once NaN and
# Inf are refused first, and err_alloc, which every wrapper reaches only
# past a size it accepts and a scan it passes.
#
# The last line is "WRAPPER TEST PASSED" or "WRAPPER TEST FAILED:
# <routines>", and the exit status is 1 on failure.

tol <- list(abs = 1e-12, rel = 1e-7, rel_floor = 1e-6, vec = 1e-12)

args <- commandArgs(trailingOnly = TRUE)
unknown <- args[!grepl("^--(src|build)=.+$|^--verbose$", args)]
if (length(unknown)) {
  stop(
    "usage: Rscript dev/si_lapack_test.R [--src=<dir>] [--build=<dir>] ",
    "[--verbose]",
    call. = FALSE
  )
}
opt_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
verbose <- "--verbose" %in% args
src <- normalizePath(opt_value("src", "src"), winslash = "/", mustWork = TRUE)
build <- opt_value("build", file.path(tempdir(), "si-lapack-test"))
test_cpp <- "dev/si_lapack_test.cpp"
if (!file.exists(test_cpp)) {
  stop("run from the package root: ", test_cpp, " not found", call. = FALSE)
}

cat("R:", R.version.string, "\n")
cat(
  "Rcpp:",
  as.character(packageVersion("Rcpp")),
  "at",
  dirname(find.package("Rcpp")),
  "\n"
)
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
cat("LAPACK:", La_library(), "version", La_version(), "\n")
cat("wrapper source:", src, "\n")

## Stage the three files in one folder: sourceCpp() compiles a .cpp that sits
## beside an included header of the same name as its own translation unit.
unlink(build, recursive = TRUE)
dir.create(build, recursive = TRUE)
staged <- c(
  test_cpp,
  file.path(src, "si_lapack.h"),
  file.path(src, "si_lapack.cpp")
)
if (!all(file.copy(staged, build))) {
  stop("could not stage the sources in ", build, call. = FALSE)
}
## The flags of the package build that matter here: the link line of
## src/Makevars and src/Makevars.win, which make expands, and the diagnostics
## define. Both reach both units.
Sys.setenv(
  PKG_CPPFLAGS = "-DPCA_IMP_DIAGNOSTICS=1",
  PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"
)
env <- new.env()
Rcpp::sourceCpp(
  file.path(build, "si_lapack_test.cpp"),
  env = env,
  cacheDir = file.path(build, "cache"),
  rebuild = TRUE,
  showOutput = TRUE
)
if (!identical(env$diagnostics_built(), 1L)) {
  stop("the test unit was not built with PCA_IMP_DIAGNOSTICS=1", call. = FALSE)
}
codes <- env$status_codes()
cat("\n")

options(matprod = "internal")
set.seed(20261119)
sentinel <- -777.25
big <- 2^31 # one past the largest int

# -- recording -----------------------------------------------------------------

values_log <- list()
cases_log <- list()
failures <- list()

maxabs <- function(x) if (length(x)) max(abs(x)) else 0

## got against ref, as the header describes.
compare <- function(got, ref) {
  got <- as.vector(got)
  ref <- as.vector(ref)
  bad <- function(why) {
    list(
      ok = FALSE,
      why = why,
      abs = NA,
      abs_s = NA,
      rel = NA,
      n_rel = 0L,
      bitwise = FALSE
    )
  }
  if (length(got) != length(ref)) {
    return(bad(sprintf("%d cells, the reference %d", length(got), length(ref))))
  }
  nf_g <- !is.finite(got)
  nf_r <- !is.finite(ref)
  if (!identical(nf_g, nf_r) || !identical(got[nf_g], ref[nf_r])) {
    return(bad("NaN or Inf in other cells than the reference's"))
  }
  g <- got[!nf_g]
  r <- ref[!nf_r]
  d <- abs(g - r)
  s <- max(1, abs(g), abs(r))
  m <- pmax(abs(g), abs(r))
  sel <- m > tol$rel_floor * s
  worst_abs <- if (length(d)) max(d) else 0
  worst_rel <- if (any(sel)) max(d[sel] / m[sel]) else 0
  ok <- worst_abs <= tol$abs * s && worst_rel <= tol$rel
  list(
    ok = ok,
    why = sprintf(
      "abs %.2e (bound %.2e), rel %.2e over %d cells (bound %.0e)",
      worst_abs,
      tol$abs * s,
      worst_rel,
      sum(sel),
      tol$rel
    ),
    abs = worst_abs,
    abs_s = worst_abs / s,
    rel = worst_rel,
    n_rel = sum(sel),
    bitwise = identical(got, ref)
  )
}

fail <- function(routine, case, what) {
  failures[[length(failures) + 1L]] <<- list(
    routine = routine,
    case = case,
    what = what
  )
  cat(sprintf("FAIL %s | %s | %s\n", routine, case, what))
}

vec_log <- list()

## A check on vectors: the largest deviation, already divided by the scale
## where it has one, against tol$vec. Logged for the report.
vec_ok <- function(routine, case, what, x) {
  vec_log[[length(vec_log) + 1L]] <<- data.frame(
    routine = routine,
    case = case,
    what = what,
    x = x
  )
  x <= tol$vec
}

## One case. `expect` is the name of a status code, or LAPACK's info as a
## number; `values` a named list of compare() results; `checks` a named
## logical vector.
record <- function(
  routine,
  case,
  status,
  expect,
  values = list(),
  checks = logical()
) {
  want <- if (is.character(expect)) codes[[expect]] else as.integer(expect)
  n_before <- length(failures)
  if (!identical(as.integer(status), as.integer(want))) {
    fail(
      routine,
      case,
      sprintf("status %d, expected %s (%d)", status, expect, want)
    )
  }
  for (nm in names(values)) {
    v <- values[[nm]]
    values_log[[length(values_log) + 1L]] <<- data.frame(
      routine = routine,
      case = case,
      value = nm,
      abs = v$abs,
      abs_s = v$abs_s,
      rel = v$rel,
      n_rel = v$n_rel,
      bitwise = v$bitwise
    )
    if (!isTRUE(v$ok)) {
      fail(routine, case, paste0(nm, ": ", v$why))
    }
  }
  for (nm in names(checks)) {
    if (!isTRUE(checks[[nm]])) {
      fail(routine, case, nm)
    }
  }
  passed <- length(failures) == n_before
  cases_log[[length(cases_log) + 1L]] <<- data.frame(
    routine = routine,
    case = case,
    passed = passed
  )
  if (verbose && passed) {
    cat(sprintf("ok   %s | %s\n", routine, case))
  }
}

## Calls the wrapper, then judges what came back; an R error in either is
## the case failing.
run <- function(routine, case, call, judge) {
  out <- tryCatch(call, error = function(e) e)
  if (inherits(out, "error")) {
    fail(routine, case, paste("R error:", conditionMessage(out)))
    cases_log[[length(cases_log) + 1L]] <<- data.frame(
      routine = routine,
      case = case,
      passed = FALSE
    )
    return(invisible())
  }
  res <- tryCatch(judge(out), error = function(e) e)
  if (inherits(res, "error")) {
    fail(routine, case, paste("R error in the judging:", conditionMessage(res)))
    cases_log[[length(cases_log) + 1L]] <<- data.frame(
      routine = routine,
      case = case,
      passed = FALSE
    )
  }
  invisible()
}

# -- inputs --------------------------------------------------------------------

rmat <- function(m, n) matrix(rnorm(m * n), m, n)
gram <- function(r, n) crossprod(rmat(r, n))

## M in the first rows of a buffer of ld rows; the rows below hold `fill`.
pad_rows <- function(M, ld, fill) {
  out <- matrix(fill, ld, ncol(M))
  if (nrow(M)) {
    out[seq_len(nrow(M)), ] <- M
  }
  out
}

## The triangle LAPACK does not read: TRUE below the diagonal for 'U'.
other_tri <- function(uplo, n) {
  M <- matrix(0, n, n)
  if (uplo == "U") lower.tri(M) else upper.tri(M)
}

# -- gemm ----------------------------------------------------------------------

## Padding rows of the operands hold NaN, so a wrong leading dimension reads
## them; padding rows of C hold the sentinel and must be left alone.
gemm_case <- function(
  ta,
  tb,
  m,
  n,
  k,
  alpha = 1,
  beta = 0,
  pad = 0,
  c_nan = FALSE,
  tag = ""
) {
  A <- if (ta == "N") rmat(m, k) else rmat(k, m)
  B <- if (tb == "N") rmat(k, n) else rmat(n, k)
  C <- rmat(m, n)
  if (c_nan) {
    C[] <- NaN
  }
  op_a <- if (ta == "N") A else t(A)
  op_b <- if (tb == "N") B else t(B)
  lda <- max(1, nrow(A)) + pad
  ldb <- max(1, nrow(B)) + pad
  ldc <- max(1, m) + pad
  a_buf <- pad_rows(A, lda, NaN)
  b_buf <- pad_rows(B, ldb, NaN)
  c_buf <- pad_rows(C, ldc, sentinel)
  prod <- if (k == 0) matrix(0, m, n) else op_a %*% op_b
  ref <- if (beta == 0) alpha * prod else alpha * prod + beta * C
  case <- trimws(sprintf("%s%s %d x %d, inner %d %s", ta, tb, m, n, k, tag))
  run(
    "gemm",
    case,
    env$t_gemm(
      ta,
      tb,
      m,
      n,
      k,
      alpha,
      a_buf,
      lda,
      b_buf,
      ldb,
      beta,
      c_buf,
      ldc
    ),
    function(out) {
      got <- matrix(out$c, ldc, n)
      below <- seq_len(ldc) > m
      record(
        "gemm",
        case,
        out$status,
        "ok",
        values = list(c = compare(got[!below, , drop = FALSE], ref)),
        checks = c(
          "the rows of C past m are as they were" = identical(
            got[below, , drop = FALSE],
            c_buf[below, , drop = FALSE]
          ),
          "the cells after C are as they were" = out$guard_c
        )
      )
    }
  )
}

for (ta in c("N", "T")) {
  for (tb in c("N", "T")) {
    gemm_case(ta, tb, 5, 4, 3)
    gemm_case(
      ta,
      tb,
      6,
      7,
      5,
      alpha = -1.5,
      beta = 0.5,
      pad = 2,
      tag = "alpha -1.5, beta 0.5, leading dimensions + 2"
    )
    gemm_case(ta, tb, 1, 6, 4, tag = "one row")
    gemm_case(ta, tb, 6, 1, 4, tag = "one column")
    gemm_case(ta, tb, 1, 1, 1000, tag = "row times column")
    gemm_case(ta, tb, 3, 3, 1, beta = 1, tag = "beta 1")
  }
}
gemm_case("N", "N", 40, 30, 20, tag = "larger")
gemm_case(
  "N",
  "N",
  4,
  3,
  0,
  beta = 0,
  c_nan = TRUE,
  pad = 1,
  tag = "beta 0, NaN in C"
)
gemm_case("N", "T", 4, 3, 0, beta = 2, tag = "beta 2")
gemm_case("T", "N", 4, 3, 0, beta = 1, tag = "beta 1")
gemm_case("N", "N", 0, 3, 2, pad = 1, tag = "no rows")
gemm_case("N", "N", 4, 0, 2, tag = "no columns")

## Refused before dgemm is entered: C is left as it was. The buffers are a
## few cells whatever the sizes say.
gemm_refused <- function(case, expect, ta, tb, m, n, k, lda, ldb, ldc) {
  a <- rep(1, 12)
  b <- rep(2, 12)
  c0 <- rep(3, 12)
  run(
    "gemm",
    case,
    env$t_gemm(ta, tb, m, n, k, 1, a, lda, b, ldb, 0, c0, ldc),
    function(out) {
      record(
        "gemm",
        case,
        out$status,
        expect,
        checks = c(
          "C is as it was" = identical(out$c, c0),
          "the cells after C are as they were" = out$guard_c
        )
      )
    }
  )
}
gemm_refused("trans_a 'C'", "err_arg", "C", "N", 2, 2, 2, 2, 2, 2)
gemm_refused("trans_b 'X'", "err_arg", "N", "X", 2, 2, 2, 2, 2, 2)
gemm_refused("ldc < m", "err_arg", "N", "N", 3, 2, 2, 3, 2, 2)
gemm_refused("lda < m, 'N'", "err_arg", "N", "N", 3, 2, 2, 2, 2, 3)
gemm_refused("lda < k, 'T'", "err_arg", "T", "N", 2, 2, 3, 2, 3, 2)
gemm_refused("ldb < k, 'N'", "err_arg", "N", "N", 2, 2, 3, 2, 2, 2)
gemm_refused("ldb < n, 'T'", "err_arg", "N", "T", 2, 3, 2, 2, 2, 2)
gemm_refused("m past int", "err_size", "N", "N", big, 1, 1, big, 1, big)
gemm_refused("k past int", "err_size", "N", "N", 1, 1, big, 1, big, 1)
gemm_refused("lda past int", "err_size", "N", "N", 1, 1, 1, big, 1, 1)

# -- gemv ----------------------------------------------------------------------

gemv_case <- function(
  trans,
  m,
  n,
  alpha = 1,
  beta = 0,
  pad = 0,
  y_nan = FALSE,
  tag = ""
) {
  A <- rmat(m, n)
  lda <- max(1, m) + pad
  a_buf <- pad_rows(A, lda, NaN)
  len_x <- if (trans == "N") n else m
  len_y <- if (trans == "N") m else n
  x <- rnorm(len_x)
  y <- rnorm(len_y)
  if (y_nan) {
    y[] <- NaN
  }
  op_a <- if (trans == "N") A else t(A)
  prod <- if (len_x == 0) rep(0, len_y) else as.vector(op_a %*% x)
  ref <- if (beta == 0) alpha * prod else alpha * prod + beta * y
  case <- trimws(sprintf("%s %d x %d %s", trans, m, n, tag))
  run(
    "gemv",
    case,
    env$t_gemv(trans, m, n, alpha, a_buf, lda, x, beta, y),
    function(out) {
      record(
        "gemv",
        case,
        out$status,
        "ok",
        values = list(y = compare(out$y, ref)),
        checks = c("the cells after y are as they were" = out$guard_y)
      )
    }
  )
}

for (trans in c("N", "T")) {
  gemv_case(trans, 5, 3)
  gemv_case(
    trans,
    6,
    4,
    alpha = 2.5,
    beta = -0.75,
    pad = 3,
    tag = "alpha 2.5, beta -0.75, lda + 3"
  )
  gemv_case(trans, 1, 5, tag = "one row")
  gemv_case(trans, 5, 1, tag = "one column")
  gemv_case(trans, 1, 1, beta = 1, tag = "beta 1")
}
gemv_case("T", 1000, 1, tag = "a column times a vector")
gemv_case("N", 1, 1000, tag = "a row times a vector")
gemv_case("N", 4, 0, beta = 0, y_nan = TRUE, tag = "x empty, beta 0, NaN in y")
gemv_case("T", 0, 4, beta = 3, tag = "x empty, beta 3")
gemv_case("N", 0, 4, tag = "y empty")
gemv_case("T", 4, 0, tag = "y empty")

gemv_refused <- function(case, expect, trans, m, n, lda) {
  a <- rep(1, 12)
  x <- rep(2, 12)
  y0 <- rep(3, 12)
  run(
    "gemv",
    case,
    env$t_gemv(trans, m, n, 1, a, lda, x, 0, y0),
    function(out) {
      record(
        "gemv",
        case,
        out$status,
        expect,
        checks = c(
          "y is as it was" = identical(out$y, y0),
          "the cells after y are as they were" = out$guard_y
        )
      )
    }
  )
}
gemv_refused("trans 'C'", "err_arg", "C", 2, 2, 2)
gemv_refused("lda < m", "err_arg", "N", 3, 2, 2)
gemv_refused("m past int", "err_size", "N", big, 1, big)
gemv_refused("n past int", "err_size", "N", 1, big, 1)

# -- syrk ----------------------------------------------------------------------

## The triangle not named holds the sentinel and must be left alone, as must
## padding rows.
syrk_case <- function(
  uplo,
  trans,
  n,
  k,
  alpha = 1,
  beta = 0,
  pad = 0,
  c_nan = FALSE,
  tag = ""
) {
  A <- if (trans == "N") rmat(n, k) else rmat(k, n)
  lda <- max(1, nrow(A)) + pad
  a_buf <- pad_rows(A, lda, NaN)
  C <- rmat(n, n)
  if (c_nan) {
    C[] <- NaN
  }
  other <- other_tri(uplo, n)
  C[other] <- sentinel
  ldc <- max(1, n) + pad
  c_buf <- pad_rows(C, ldc, sentinel)
  prod <- if (k == 0) {
    matrix(0, n, n)
  } else if (trans == "N") {
    A %*% t(A)
  } else {
    t(A) %*% A
  }
  ref <- if (beta == 0) alpha * prod else alpha * prod + beta * C
  case <- trimws(sprintf("%s%s n %d, k %d %s", uplo, trans, n, k, tag))
  run(
    "syrk",
    case,
    env$t_syrk(uplo, trans, n, k, alpha, a_buf, lda, beta, c_buf, ldc),
    function(out) {
      got <- matrix(out$c, ldc, n)
      below <- seq_len(ldc) > n
      top <- got[!below, , drop = FALSE]
      record(
        "syrk",
        case,
        out$status,
        "ok",
        values = list(c = compare(top[!other], ref[!other])),
        checks = c(
          "the other triangle of C is as it was" = identical(
            top[other],
            C[other]
          ),
          "the rows of C past n are as they were" = identical(
            got[below, , drop = FALSE],
            c_buf[below, , drop = FALSE]
          ),
          "the cells after C are as they were" = out$guard_c
        )
      )
    }
  )
}

for (uplo in c("U", "L")) {
  for (trans in c("N", "T")) {
    syrk_case(uplo, trans, 5, 3)
    syrk_case(
      uplo,
      trans,
      3,
      5,
      alpha = 0.5,
      beta = 2,
      pad = 2,
      tag = "alpha 0.5, beta 2, leading dimensions + 2"
    )
    syrk_case(uplo, trans, 1, 4, tag = "n 1")
    syrk_case(uplo, trans, 4, 1, tag = "k 1")
    syrk_case(uplo, trans, 1, 1, beta = 1, tag = "beta 1")
    syrk_case(uplo, trans, 6, 1000, tag = "1000 terms")
    syrk_case(
      uplo,
      trans,
      4,
      0,
      beta = 0,
      c_nan = TRUE,
      pad = 1,
      tag = "beta 0, NaN in the triangle"
    )
    syrk_case(uplo, trans, 4, 0, beta = -2, tag = "beta -2")
  }
}
syrk_case("U", "N", 0, 3, pad = 1, tag = "n 0")

syrk_refused <- function(case, expect, uplo, trans, n, k, lda, ldc) {
  a <- rep(1, 12)
  c0 <- rep(3, 12)
  run(
    "syrk",
    case,
    env$t_syrk(uplo, trans, n, k, 1, a, lda, 0, c0, ldc),
    function(out) {
      record(
        "syrk",
        case,
        out$status,
        expect,
        checks = c(
          "C is as it was" = identical(out$c, c0),
          "the cells after C are as they were" = out$guard_c
        )
      )
    }
  )
}
syrk_refused("uplo 'X'", "err_arg", "X", "N", 2, 2, 2, 2)
syrk_refused("trans 'C'", "err_arg", "U", "C", 2, 2, 2, 2)
syrk_refused("ldc < n", "err_arg", "U", "N", 3, 2, 3, 2)
syrk_refused("lda < n, 'N'", "err_arg", "U", "N", 3, 2, 2, 3)
syrk_refused("lda < k, 'T'", "err_arg", "L", "T", 2, 3, 2, 2)
syrk_refused("n past int", "err_size", "U", "N", big, 1, big, big)
syrk_refused("k past int", "err_size", "U", "T", 1, big, big, 1)

# -- syevr_range ---------------------------------------------------------------

## The triangle not named holds unrelated values, or NaN, which the wrapper
## must not read. With `cut`, the range cuts a group of equal eigenvalues,
## and the probe counts how many cells past the count dsyevr writes.
probe_cases <- 0L
probe_cells <- 0L
probe_hit <- 0L
syevr_case <- function(
  uplo,
  S,
  il,
  iu,
  tag,
  expect = "ok",
  junk_nan = FALSE,
  cut = FALSE
) {
  n <- nrow(S)
  A <- S
  other <- other_tri(if (uplo == "U") "U" else "L", n)
  A[other] <- if (junk_nan) NaN else rnorm(sum(other))
  count <- max(0, iu - il + 1)
  w0 <- rep(sentinel, count)
  z0 <- rep(sentinel, n * count)
  case <- trimws(sprintf("%s n %d pairs %d..%d %s", uplo, n, il, iu, tag))
  run(
    "syevr_range",
    case,
    env$t_syevr_range(uplo, n, il, iu, A, w0, z0),
    function(out) {
      checks <- c(
        "the cells after a are as they were" = out$guard_a,
        "the cells after w are as they were" = out$guard_w,
        "the cells after z are as they were" = out$guard_z
      )
      values <- list()
      if (identical(expect, "ok") && out$status == codes$ok) {
        ref <- rev(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
        w <- out$w
        z <- matrix(out$z, n, count)
        scale <- max(1, abs(S))
        values$w <- compare(w, ref[il:iu])
        checks <- c(
          checks,
          "w is ascending" = !is.unsorted(w),
          "z holds eigenvectors of S (residual)" = vec_ok(
            "syevr_range",
            case,
            "residual",
            maxabs(S %*% z - z %*% diag(w, count)) / scale
          ),
          "z is orthonormal" = vec_ok(
            "syevr_range",
            case,
            "orthonormality",
            maxabs(t(z) %*% z - diag(count))
          )
        )
      } else {
        checks <- c(
          checks,
          "w and z are as they were" = identical(out$w, w0) &&
            identical(out$z, z0)
        )
      }
      record("syevr_range", case, out$status, expect, values, checks)
    }
  )
  if (cut) {
    p <- env$probe_dsyevr_past(uplo, n, il, iu, A)
    probe_cases <<- probe_cases + 1L
    if (identical(p$info, 0L) && p$past > 0L) {
      probe_hit <<- probe_hit + 1L
      probe_cells <<- probe_cells + p$past
    }
  }
}

d6 <- diag(c(1, 1, 2, 2, 5, 5))
for (ul in c("U", "L")) {
  g <- gram(40, 10)
  g4 <- gram(4, 10)
  syevr_case(ul, g, 8, 10, "top 3")
  syevr_case(ul, g, 10, 10, "top 1")
  syevr_case(ul, g, 1, 10, "all")
  syevr_case(ul, g, 1, 2, "bottom 2")
  syevr_case(ul, g, 4, 6, "middle 3")
  syevr_case(ul, matrix(2.5), 1, 1, "order 1")
  syevr_case(ul, gram(80, 60), 56, 60, "top 5")
  syevr_case(ul, gram(400, 300), 291, 300, "top 10")
  syevr_case(ul, g4, 7, 10, "rank 4, top 4")
  syevr_case(ul, g4, 1, 10, "rank 4, all")
  syevr_case(ul, diag(6), 1, 6, "identity, all")
  syevr_case(ul, matrix(0, 6, 6), 1, 6, "zero, all")
  syevr_case(ul, d6, 5, 6, "repeated, top pair")
  indef <- rmat(12, 12)
  syevr_case(ul, (indef + t(indef)) / 2, 10, 12, "indefinite, top 3")
  syevr_case(ul, 100 * gram(30, 8), 6, 8, "scale 100, top 3")
  syevr_case(ul, g, 8, 10, "top 3, NaN in the other triangle", junk_nan = TRUE)

  ## Ranges that cut a group of equal eigenvalues.
  syevr_case(ul, diag(6), 5, 6, "identity, top 2", cut = TRUE)
  syevr_case(ul, d6, 6, 6, "top 1 of a pair", cut = TRUE)
  syevr_case(ul, d6, 4, 6, "top 3, one of a pair", cut = TRUE)
  syevr_case(ul, g4, 5, 10, "rank 4, top 6", cut = TRUE)
  syevr_case(ul, g4, 6, 10, "rank 4, top 5", cut = TRUE)
  ## What pca_imp() hands it when asked for more components than the data
  ## has rank: 30 samples of which 6 are distinct, columns centred, so the
  ## 30 x 30 Gram has rank 5.
  few <- rmat(6, 100)[rep(1:6, each = 5), ]
  few <- sweep(few, 2, colMeans(few))
  syevr_case(
    ul,
    few %*% t(few),
    23,
    30,
    "6 distinct rows of 30, top 8",
    cut = TRUE
  )

  ## NaN and Inf in the named triangle. NaN is tried on a proper subset only:
  ## dsyevr asked for every eigenpair of such a matrix may not return, and a
  ## self-test of a wrapper without its scan would hang.
  bad <- gram(20, 6)
  bad[2, 2] <- NaN
  syevr_case(ul, bad, 5, 6, "NaN on the diagonal", expect = "err_nonfinite")
  bad <- gram(20, 6)
  bad[1, 3] <- bad[3, 1] <- NaN
  syevr_case(ul, bad, 5, 6, "NaN off the diagonal", expect = "err_nonfinite")
  bad <- gram(20, 6)
  bad[1, 3] <- bad[3, 1] <- -Inf
  syevr_case(ul, bad, 5, 6, "-Inf off the diagonal", expect = "err_nonfinite")
  bad <- gram(20, 6)
  bad[2, 2] <- Inf
  syevr_case(
    ul,
    bad,
    1,
    6,
    "Inf on the diagonal, all",
    expect = "err_nonfinite"
  )

  syevr_case(ul, g, 0, 3, "il 0", expect = "err_arg")
  syevr_case(ul, g, 4, 3, "il > iu", expect = "err_arg")
  syevr_case(ul, g, 9, 11, "iu > n", expect = "err_arg")
}
syevr_case("X", gram(40, 10), 8, 10, "uplo 'X'", expect = "err_arg")

syevr_refused <- function(case, n, il, iu) {
  run(
    "syevr_range",
    case,
    env$t_syevr_range("U", n, il, iu, 1, sentinel, sentinel),
    function(out) {
      record(
        "syevr_range",
        case,
        out$status,
        "err_size",
        checks = c(
          "w and z are as they were" = identical(out$w, sentinel) &&
            identical(out$z, sentinel),
          "the cells after w are as they were" = out$guard_w,
          "the cells after z are as they were" = out$guard_z
        )
      )
    }
  )
}
syevr_refused("n past int", big, 1, 1)
syevr_refused("26 n past int", 1e8, 1e8, 1e8)

if (probe_hit == 0L) {
  fail(
    "syevr_range",
    "(the cut cases)",
    "dsyevr wrote past iu - il + 1 in no case: the guard cells after w test nothing"
  )
}

# -- syevd ---------------------------------------------------------------------

syevd_case <- function(uplo, S, tag, expect = "ok", junk_nan = FALSE) {
  n <- nrow(S)
  A <- S
  other <- other_tri(if (uplo == "U") "U" else "L", n)
  A[other] <- if (junk_nan) NaN else rnorm(sum(other))
  w0 <- rep(sentinel, n)
  case <- trimws(sprintf("%s n %d %s", uplo, n, tag))
  run("syevd", case, env$t_syevd(uplo, n, A, w0), function(out) {
    checks <- c(
      "the cells after a are as they were" = out$guard_a,
      "the cells after w are as they were" = out$guard_w
    )
    values <- list()
    if (identical(expect, "ok") && out$status == codes$ok && n > 0) {
      ref <- rev(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
      w <- out$w
      z <- matrix(out$a, n, n)
      scale <- max(1, abs(S))
      values$w <- compare(w, ref)
      checks <- c(
        checks,
        "w is ascending" = !is.unsorted(w),
        "a holds eigenvectors of S (residual)" = vec_ok(
          "syevd",
          case,
          "residual",
          maxabs(S %*% z - z %*% diag(w, n)) / scale
        ),
        "a is orthonormal" = vec_ok(
          "syevd",
          case,
          "orthonormality",
          maxabs(t(z) %*% z - diag(n))
        )
      )
    } else {
      checks <- c(
        checks,
        "a and w are as they were" = identical(out$a, as.vector(A)) &&
          identical(out$w, w0)
      )
    }
    record("syevd", case, out$status, expect, values, checks)
  })
}

for (ul in c("U", "L")) {
  syevd_case(ul, gram(40, 10), "")
  syevd_case(ul, gram(80, 60), "")
  syevd_case(ul, matrix(2.5), "order 1")
  syevd_case(ul, diag(6), "identity")
  syevd_case(ul, matrix(0, 6, 6), "zero")
  syevd_case(ul, d6, "repeated")
  indef <- rmat(12, 12)
  syevd_case(ul, (indef + t(indef)) / 2, "indefinite")
  syevd_case(ul, 100 * gram(30, 8), "scale 100")
  syevd_case(ul, gram(4, 10), "rank 4")
  syevd_case(ul, gram(40, 10), "NaN in the other triangle", junk_nan = TRUE)
  syevd_case(ul, matrix(0, 0, 0), "order 0")
  bad <- gram(20, 6)
  bad[2, 2] <- NaN
  syevd_case(ul, bad, "NaN on the diagonal", expect = "err_nonfinite")
  bad <- gram(20, 6)
  bad[1, 3] <- bad[3, 1] <- Inf
  syevd_case(ul, bad, "Inf off the diagonal", expect = "err_nonfinite")
}
syevd_case("X", gram(40, 10), "uplo 'X'", expect = "err_arg")

syevd_refused <- function(case, n) {
  run("syevd", case, env$t_syevd("U", n, 1, sentinel), function(out) {
    record(
      "syevd",
      case,
      out$status,
      "err_size",
      checks = c(
        "a and w are as they were" = identical(out$a, 1) &&
          identical(out$w, sentinel),
        "the cells after w are as they were" = out$guard_w
      )
    )
  })
}
syevd_refused("n past int", big)
syevd_refused("2 n^2 past int", 40000)

# -- qr_econ -------------------------------------------------------------------

## q r is compared with x, and |r| with |R| of base::qr() where x has full
## rank: the factor is then unique up to the sign of each row.
qr_case <- function(X, tag, expect = "ok", full_rank = TRUE) {
  m <- nrow(X)
  n <- ncol(X)
  q0 <- rep(sentinel, m * n)
  r0 <- rep(sentinel, n * n)
  case <- trimws(sprintf("%d x %d %s", m, n, tag))
  run("qr_econ", case, env$t_qr_econ(m, n, X, q0, r0), function(out) {
    checks <- c(
      "x is as it was" = identical(out$x, as.vector(X)),
      "the cells after q are as they were" = out$guard_q,
      "the cells after r are as they were" = out$guard_r
    )
    values <- list()
    if (identical(expect, "ok") && out$status == codes$ok && n > 0) {
      Q <- matrix(out$q, m, n)
      R <- matrix(out$r, n, n)
      values$qr <- compare(Q %*% R, X)
      checks <- c(
        checks,
        "q is orthonormal" = vec_ok(
          "qr_econ",
          case,
          "orthonormality",
          maxabs(t(Q) %*% Q - diag(n))
        ),
        "r is zero below the diagonal" = all(R[lower.tri(R)] == 0)
      )
      if (full_rank) {
        ref <- qr(X)
        checks <- c(
          checks,
          "base::qr() did not pivot" = identical(ref$pivot, seq_len(n))
        )
        values$abs_r <- compare(abs(R), abs(qr.R(ref)))
      }
    } else {
      checks <- c(
        checks,
        "q and r are as they were" = identical(out$q, q0) &&
          identical(out$r, r0)
      )
    }
    record("qr_econ", case, out$status, expect, values, checks)
  })
}

qr_case(rmat(10, 4), "tall")
qr_case(rmat(5, 5), "square")
qr_case(rmat(7, 1), "one column")
qr_case(matrix(-3), "1 x 1")
qr_case(rmat(200, 12), "tall, larger")
qr_case(1000 * rmat(30, 6), "scale 1000")
dup <- rmat(9, 2)
qr_case(cbind(dup, dup[, 1]), "a column repeated", full_rank = FALSE)
qr_case(matrix(0, 4, 3), "zero", full_rank = FALSE)
qr_case(matrix(0, 5, 0), "no columns")
qr_case(rmat(3, 5), "wide", expect = "err_arg")
run(
  "qr_econ",
  "m past int",
  env$t_qr_econ(big, 1, 1, sentinel, sentinel),
  function(out) {
    record(
      "qr_econ",
      "m past int",
      out$status,
      "err_size",
      checks = c(
        "q and r are as they were" = identical(out$q, sentinel) &&
          identical(out$r, sentinel),
        "the cells after q are as they were" = out$guard_q
      )
    )
  }
)

# -- solve_tri_right -----------------------------------------------------------

## X with X * tri(A) = B, B of m rows and n columns. A triangle well away
## from singular: a diagonal of either sign and magnitude above n,
## unit-variance entries elsewhere. The triangle not named holds NaN, which
## the wrapper must not read. A zero pivot at `zero_at` makes the wrapper
## return its index and leave b alone. The reference solves the transposed
## system, t(tri) %*% t(X) = t(B), with forwardsolve() or backsolve().
solve_case <- function(uplo, m, n, tag, zero_at = 0L, b_nan = FALSE) {
  A <- rmat(n, n)
  diag(A) <- sample(c(-1, 1), n, replace = TRUE) * (n + abs(rnorm(n)))
  other <- other_tri(if (uplo == "U") "U" else "L", n)
  tri <- A
  tri[other] <- 0
  A[other] <- NaN
  if (zero_at > 0L) {
    A[zero_at, zero_at] <- 0
    tri[zero_at, zero_at] <- 0
  }
  B <- rmat(m, n)
  if (b_nan) {
    B[1, 1] <- NaN
  }
  expect <- if (zero_at > 0L) zero_at else "ok"
  case <- trimws(sprintf("%s m %d, n %d %s", uplo, m, n, tag))
  run(
    "solve_tri_right",
    case,
    env$t_solve_tri_right(uplo, m, n, A, B),
    function(out) {
      checks <- c("the cells after b are as they were" = out$guard_b)
      values <- list()
      if (zero_at > 0L || m == 0 || n == 0) {
        checks <- c(checks, "b is as it was" = identical(out$b, as.vector(B)))
      } else {
        X <- matrix(out$b, m, n)
        ref <- if (uplo == "U") {
          t(forwardsolve(t(tri), t(B)))
        } else {
          t(backsolve(t(tri), t(B)))
        }
        values$x <- compare(X, ref)
        if (!b_nan) {
          values$residual <- compare(X %*% tri, B)
        }
      }
      record("solve_tri_right", case, out$status, expect, values, checks)
    }
  )
}

for (ul in c("U", "L")) {
  solve_case(ul, 1, 6, "")
  solve_case(ul, 4, 6, "")
  solve_case(ul, 1, 1, "")
  solve_case(ul, 3, 40, "")
  solve_case(ul, 200, 6, "")
  solve_case(ul, 2, 6, "zero pivot 3", zero_at = 3L)
  solve_case(ul, 2, 6, "zero pivot 1", zero_at = 1L)
  solve_case(ul, 2, 6, "zero pivot 6", zero_at = 6L)
  solve_case(ul, 2, 6, "NaN in b, not refused", b_nan = TRUE)
  solve_case(ul, 2, 0, "")
  solve_case(ul, 0, 4, "")
}

solve_refused <- function(case, expect, uplo, m, n) {
  a <- rep(1, 12)
  b0 <- rep(3, 12)
  run(
    "solve_tri_right",
    case,
    env$t_solve_tri_right(uplo, m, n, a, b0),
    function(out) {
      record(
        "solve_tri_right",
        case,
        out$status,
        expect,
        checks = c(
          "b is as it was" = identical(out$b, b0),
          "the cells after b are as they were" = out$guard_b
        )
      )
    }
  )
}
solve_refused("uplo 'X'", "err_arg", "X", 2, 2)
solve_refused("m past int", "err_size", "U", big, 1)
solve_refused("n past int", "err_size", "U", 1, big)

# -- svd_values ----------------------------------------------------------------

svd_case <- function(X, tag, expect = "ok") {
  m <- nrow(X)
  n <- ncol(X)
  p <- min(m, n)
  s0 <- rep(sentinel, p)
  case <- trimws(sprintf("%d x %d %s", m, n, tag))
  run("svd_values", case, env$t_svd_values(m, n, X, s0), function(out) {
    checks <- c(
      "the cells after a are as they were" = out$guard_a,
      "the cells after s are as they were" = out$guard_s
    )
    values <- list()
    if (identical(expect, "ok") && out$status == codes$ok && p > 0) {
      values$s <- compare(out$s, svd(X, nu = 0, nv = 0)$d)
      checks <- c(checks, "s is descending" = !is.unsorted(rev(out$s)))
    } else {
      checks <- c(checks, "s is as it was" = identical(out$s, s0))
    }
    record("svd_values", case, out$status, expect, values, checks)
  })
}

svd_case(rmat(10, 4), "tall")
svd_case(rmat(4, 10), "wide")
svd_case(rmat(6, 6), "square")
svd_case(rmat(1, 5), "one row")
svd_case(rmat(5, 1), "one column")
svd_case(rmat(8, 2) %*% rmat(2, 5), "rank 2")
svd_case(matrix(0, 3, 3), "zero")
svd_case(1000 * rmat(20, 7), "scale 1000")
svd_case(rmat(60, 40), "larger")
svd_case(matrix(0, 0, 3), "no rows")
svd_case(matrix(0, 3, 0), "no columns")
bad <- rmat(5, 4)
bad[2, 3] <- NaN
svd_case(bad, "NaN", expect = "err_nonfinite")
bad <- rmat(5, 4)
bad[4, 1] <- -Inf
svd_case(bad, "-Inf", expect = "err_nonfinite")
run(
  "svd_values",
  "m past int",
  env$t_svd_values(big, 1, 1, sentinel),
  function(out) {
    record(
      "svd_values",
      "m past int",
      out$status,
      "err_size",
      checks = c(
        "s is as it was" = identical(out$s, sentinel),
        "the cells after s are as they were" = out$guard_s
      )
    )
  }
)

# -- report --------------------------------------------------------------------

## Per routine: the cases, the worst absolute difference as it is and divided
## by its scale (the figure asserted on), the worst relative difference, and
## the worst of each vector check. Each worst figure names its case and
## value.
cases <- do.call(rbind, cases_log)
vals <- do.call(rbind, values_log)
vecs <- do.call(rbind, vec_log)
at <- function(df, i, col = "value") sprintf("%s, %s", df$case[i], df[[col]][i])
cat("\n")
for (r in unique(cases$routine)) {
  cr <- cases[cases$routine == r, ]
  vr <- vals[vals$routine == r & is.finite(vals$abs), ]
  failed_cases <- unique(vapply(
    Filter(function(f) f$routine == r, failures),
    `[[`,
    "",
    "case"
  ))
  cat(sprintf(
    "%s: %d cases, %d failed; %d values compared, %d bitwise\n",
    r,
    nrow(cr),
    length(failed_cases),
    nrow(vr),
    sum(vr$bitwise)
  ))
  if (nrow(vr)) {
    ia <- which.max(vr$abs)
    is <- which.max(vr$abs_s)
    ir <- which.max(vr$rel)
    cat(sprintf("  worst abs        %.2e (%s)\n", vr$abs[ia], at(vr, ia)))
    cat(sprintf("  worst abs/scale  %.2e (%s)\n", vr$abs_s[is], at(vr, is)))
    cat(sprintf(
      "  worst rel        %.2e over %d cells (%s)\n",
      vr$rel[ir],
      vr$n_rel[ir],
      at(vr, ir)
    ))
  }
  vv <- vecs[vecs$routine == r, ]
  for (w in unique(vv$what)) {
    vw <- vv[vv$what == w, ]
    iw <- which.max(vw$x)
    cat(sprintf(
      "  worst %-16s %.2e (%s)\n",
      w,
      vw$x[iw],
      vw$case[iw]
    ))
  }
}
cat(sprintf(
  paste0(
    "\nsyevr_range cut cases: dsyevr wrote past iu - il + 1 in %d of %d, ",
    "%d cells in all\n"
  ),
  probe_hit,
  probe_cases,
  probe_cells
))
cat(sprintf(
  "tolerances: abs %.0e x scale, rel %.0e over cells above %.0e x scale, vectors %.0e\n\n",
  tol$abs,
  tol$rel,
  tol$rel_floor,
  tol$vec
))

failed_routines <- unique(vapply(failures, `[[`, "", "routine"))
if (length(failed_routines)) {
  cat("WRAPPER TEST FAILED:", paste(failed_routines, collapse = ", "), "\n")
  quit(status = 1L)
}
cat("WRAPPER TEST PASSED\n")
