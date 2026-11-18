# The Armadillo oracle (dev/eigen_conversion.md, "The Armadillo oracle"): run
# every wrapper of src/si_lapack.h beside the Armadillo call it stands in
# for, and fail when a pair disagrees. Dev only, run by hand, from the
# package root:
#   Rscript dev/oracle_run.R <lib> [--src=<dir>] [--build=<dir>] [--verbose]
# <lib>     the private library dev/oracle_install_armadillo.R installed
#           RcppArmadillo into. It is put first on the library path here, and
#           the run stops unless RcppArmadillo then resolves inside it.
# --src     the folder holding si_lapack.h and si_lapack.cpp (default: src).
#           dev/oracle_selftest.R points this at a copy with one wrapper
#           broken.
# --build   where the sources are staged and compiled (default: a folder
#           under tempdir()). Nothing is built in the working tree.
# --verbose print every case, not only the worst per routine.
# Rcpp comes from the usual library, so R_LIBS_USER stays as it is.
#
# dev/oracle_pairs.cpp holds the pairs. This script builds the inputs, calls
# each pair and decides. Per case it requires:
# - the two sides agree on whether the call failed; when both failed, no
#   value is compared;
# - the same shape, and NaN and Inf in the same cells;
# - over the finite cells, max(abs(w - a)) <= tol_abs * max(1, s), where s is
#   the largest magnitude either side holds, and max(abs(w - a) / m) <=
#   tol_rel, where m is max(abs(w), abs(a)) of the cell. The relative figure
#   is taken over the cells with m > rel_floor * s; with rel_floor = 0 that is
#   every cell either side holds as nonzero. The report says how many cells.
# Per routine it prints the worst absolute and the worst relative difference
# over its cases, with the case each came from, and how many cases were
# bitwise identical. Bitwise identity is reported, never asserted.
#
# Three tolerance classes. The readings are from 2026-11-18, RcppArmadillo
# 15.6.2-1, on one machine with one BLAS.
# - "same": the pair enters the same BLAS or LAPACK routine with the same
#   arguments, so any difference is a finding in itself. Reading: every case
#   bitwise identical. The tolerance is a few hundred units in the last
#   place of a cell near 1, far below any wrong answer, so a BLAS that does
#   not repeat itself to the bit shows as "not bitwise" without failing.
# The other two are set about a hundredfold above their readings.
# - "tiny": matmul_tiny, the products Armadillo runs through its own loops
#   (square operands up to 4 x 4) where the wrapper calls BLAS. Reading:
#   worst abs 8.9e-16, worst rel 3.8e-15, 18 of 42 cases bitwise.
# - "eig": syevr_range. Armadillo calls dsyevr nowhere, so the pair is
#   against eig_sym(), which is dsyevd: another algorithm on the same matrix.
#   Reading: worst abs 1.7e-13, in the case scaled by 100; worst rel
#   7.6e-10. The relative figure is that large by construction: an absolute
#   error near machine epsilon times the scale, divided by a cell as small
#   as 1e-6 of the scale.
# The last line is "ORACLE PASSED" or "ORACLE FAILED: <routines>", and the
# exit status is 1 on failure.

tol <- list(
  same = list(abs = 1e-13, rel = 1e-12, rel_floor = 0),
  tiny = list(abs = 1e-13, rel = 1e-12, rel_floor = 1e-6),
  eig = list(abs = 1e-11, rel = 1e-7, rel_floor = 1e-6)
)

args <- commandArgs(trailingOnly = TRUE)
opts <- grep("^--", args, value = TRUE)
pos <- setdiff(args, opts)
if (length(pos) != 1L) {
  stop(
    "usage: Rscript dev/oracle_run.R <lib> [--src=<dir>] [--build=<dir>] ",
    "[--verbose]",
    call. = FALSE
  )
}
opt_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), opts, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
verbose <- "--verbose" %in% opts
lib <- normalizePath(pos[[1]], winslash = "/", mustWork = TRUE)
src <- normalizePath(opt_value("src", "src"), winslash = "/", mustWork = TRUE)
build <- opt_value("build", file.path(tempdir(), "oracle-build"))
pairs_file <- "dev/oracle_pairs.cpp"
if (!file.exists(pairs_file)) {
  stop("run from the package root: ", pairs_file, " not found", call. = FALSE)
}

.libPaths(c(lib, .libPaths()))
arma_at <- normalizePath(
  dirname(find.package("RcppArmadillo")),
  winslash = "/"
)
if (!identical(arma_at, lib)) {
  stop(
    "RcppArmadillo resolves to ",
    arma_at,
    ", not to <lib> ",
    lib,
    call. = FALSE
  )
}
cat("R:", R.version.string, "\n")
cat(
  "RcppArmadillo:",
  as.character(packageVersion("RcppArmadillo")),
  "at",
  arma_at,
  "\n"
)
cat(
  "Rcpp:",
  as.character(packageVersion("Rcpp")),
  "at",
  dirname(find.package("Rcpp")),
  "\n"
)
si <- sessionInfo()
cat("BLAS:", si$BLAS, "\n")
cat("LAPACK:", si$LAPACK, "version", La_version(), "\n")
cat("wrapper source:", src, "\n")

## Stage the three files in one folder: sourceCpp() compiles a .cpp that sits
## beside an included header of the same name as its own translation unit.
unlink(build, recursive = TRUE)
dir.create(build, recursive = TRUE)
staged <- c(
  pairs_file,
  file.path(src, "si_lapack.h"),
  file.path(src, "si_lapack.cpp")
)
if (!all(file.copy(staged, build))) {
  stop("could not stage the sources in ", build, call. = FALSE)
}
pair_env <- new.env()
Rcpp::sourceCpp(
  file.path(build, "oracle_pairs.cpp"),
  env = pair_env,
  cacheDir = file.path(build, "cache"),
  rebuild = TRUE,
  showOutput = TRUE
)
info <- pair_env$oracle_info()
cat(
  "Armadillo headers:",
  info$armadillo,
  "| warn level",
  info$warn_level,
  "| uword bytes",
  info$uword_bytes,
  "| blas_int bytes",
  info$blas_int_bytes,
  "\n\n"
)
codes <- pair_env$wrapper_codes()

# -- comparing ---------------------------------------------------------------

rows <- list()

diff_stats <- function(w, a, rel_floor) {
  if (!identical(dim(w), dim(a)) || length(w) != length(a)) {
    return(list(shape = FALSE))
  }
  w <- as.vector(w)
  a <- as.vector(a)
  bad_w <- !is.finite(w)
  bad_a <- !is.finite(a)
  fin <- !bad_w & !bad_a
  d <- abs(w[fin] - a[fin])
  m <- pmax(abs(w[fin]), abs(a[fin]))
  s <- if (length(m)) max(m) else 0
  over <- m > rel_floor * s & m > 0
  list(
    shape = TRUE,
    nonfinite_same = identical(bad_w, bad_a) && identical(w[bad_w], a[bad_a]),
    n_cells = length(w),
    n_finite = sum(fin),
    scale = s,
    max_abs = if (length(d)) max(d) else 0,
    max_rel = if (any(over)) max(d[over] / m[over]) else 0,
    n_rel = sum(over),
    bitwise = identical(writeBin(w, raw()), writeBin(a, raw()))
  )
}

## One case of one routine. `w` and `a` are named lists of the arrays each
## side returned, compared name by name; `extra` holds further requirements
## as named logicals, all of which must be TRUE.
record <- function(
  routine,
  case,
  w_failed,
  a_failed,
  w = list(),
  a = list(),
  class = "same",
  extra = logical()
) {
  t <- tol[[class]]
  problems <- character()
  st <- list(
    n_cells = 0,
    n_finite = 0,
    n_rel = 0,
    max_abs = 0,
    max_rel = 0,
    bitwise = NA
  )
  if (!identical(w_failed, a_failed)) {
    problems <- c(
      problems,
      sprintf("failure: wrapper %s, Armadillo %s", w_failed, a_failed)
    )
  } else if (!w_failed) {
    st$bitwise <- TRUE
    for (nm in names(w)) {
      d <- diff_stats(w[[nm]], a[[nm]], t$rel_floor)
      if (!d$shape) {
        problems <- c(
          problems,
          sprintf(
            "%s: shape %s against %s",
            nm,
            paste(dim(w[[nm]]) %||% length(w[[nm]]), collapse = "x"),
            paste(dim(a[[nm]]) %||% length(a[[nm]]), collapse = "x")
          )
        )
        st$bitwise <- FALSE
        next
      }
      st$n_cells <- st$n_cells + d$n_cells
      st$n_finite <- st$n_finite + d$n_finite
      st$n_rel <- st$n_rel + d$n_rel
      st$max_abs <- max(st$max_abs, d$max_abs)
      st$max_rel <- max(st$max_rel, d$max_rel)
      st$bitwise <- st$bitwise && d$bitwise
      if (!d$nonfinite_same) {
        problems <- c(problems, sprintf("%s: NaN or Inf in other cells", nm))
      }
      if (d$max_abs > t$abs * max(1, d$scale)) {
        problems <- c(
          problems,
          sprintf("%s: abs %.3g at scale %.3g", nm, d$max_abs, d$scale)
        )
      }
      if (d$max_rel > t$rel) {
        problems <- c(
          problems,
          sprintf("%s: rel %.3g over %d cells", nm, d$max_rel, d$n_rel)
        )
      }
    }
  }
  for (nm in names(extra)) {
    if (!isTRUE(extra[[nm]])) {
      problems <- c(problems, nm)
    }
  }
  rows[[length(rows) + 1L]] <<- c(
    list(
      routine = routine,
      case = case,
      class = class,
      both_failed = w_failed && a_failed,
      ok = length(problems) == 0L,
      problems = paste(problems, collapse = "; ")
    ),
    st
  )
  invisible(NULL)
}

## Run one pair; an R error raised inside it is a failure of that case, not
## of the run.
run_pair <- function(routine, case, expr, handle) {
  out <- tryCatch(expr, error = function(e) e)
  if (inherits(out, "error")) {
    record(
      routine,
      case,
      FALSE,
      FALSE,
      extra = structure(FALSE, names = paste("R error:", conditionMessage(out)))
    )
    return(invisible(NULL))
  }
  handle(out)
}

simple <- function(routine, case, class = "same") {
  function(out) {
    record(
      routine,
      case,
      out$wrapper$status != codes[["ok"]],
      out$arma$failed,
      w = list(value = out$wrapper$value),
      a = list(value = out$arma$value),
      class = class
    )
  }
}

# -- inputs ------------------------------------------------------------------

set.seed(20261118)
rmat <- function(m, n) matrix(rnorm(m * n), m, n)
coefs <- list(c(1, 0), c(-1, 1), c(2.5, -0.5))
## How many cases each input set holds, checked non-empty at the end.
sets <- integer()
count_set <- function(name, n) {
  sets[[name]] <<- (if (name %in% names(sets)) sets[[name]] else 0L) +
    as.integer(n)
}

## gemm: op(A) is m x k, op(B) is k x n.
gemm_shapes <- list(
  c(7, 5, 4),
  c(1, 5, 4),
  c(7, 5, 1),
  c(7, 1, 4),
  c(1, 1, 1),
  c(4, 4, 4),
  c(60, 45, 30),
  c(200, 150, 120),
  c(7, 0, 4),
  c(0, 5, 4),
  c(7, 5, 0)
)
gemm_case <- function(ta, tb, sh, cf, A = NULL, tag = "") {
  m <- sh[[1]]
  k <- sh[[2]]
  n <- sh[[3]]
  if (is.null(A)) {
    A <- if (ta == "N") rmat(m, k) else rmat(k, m)
  }
  B <- if (tb == "N") rmat(k, n) else rmat(n, k)
  C0 <- rmat(m, n)
  case <- sprintf(
    "%s%s %dx%dx%d alpha %g beta %g%s",
    ta,
    tb,
    m,
    k,
    n,
    cf[[1]],
    cf[[2]],
    tag
  )
  run_pair(
    "gemm",
    case,
    pair_env$pair_gemm(ta, tb, A, B, cf[[1]], cf[[2]], C0),
    simple("gemm", case)
  )
  count_set("gemm", 1L)
}
for (ta in c("N", "T")) {
  for (tb in c("N", "T")) {
    for (sh in gemm_shapes) {
      for (cf in coefs) {
        gemm_case(ta, tb, sh, cf)
      }
    }
    ## rank-deficient, all zero, columns on very different scales, and NaN
    ## and Inf in an operand.
    sh <- c(9, 6, 5)
    dims <- if (ta == "N") c(9, 6) else c(6, 9)
    dup <- rmat(dims[[1]], dims[[2]])
    dup[, 2] <- dup[, 1]
    scaled <- sweep(
      rmat(dims[[1]], dims[[2]]),
      2,
      10^seq(-8, 8, length.out = dims[[2]]),
      "*"
    )
    holed <- rmat(dims[[1]], dims[[2]])
    holed[2, 3] <- NaN
    holed[1, 1] <- Inf
    gemm_case(ta, tb, sh, c(1, 0), dup, " rank-deficient A")
    gemm_case(ta, tb, sh, c(1, 0), dup * 0, " zero A")
    gemm_case(ta, tb, sh, c(1, 0), scaled, " scaled A")
    gemm_case(ta, tb, sh, c(-1, 1), holed, " NaN and Inf in A")
  }
}

## gemm on a block of a larger buffer, the two forms the chain writes.
gemm_block_case <- function(kind, A, B, j0, nj, C0) {
  case <- sprintf(
    "%s A %dx%d B %dx%d block %d+%d",
    kind,
    nrow(A),
    ncol(A),
    nrow(B),
    ncol(B),
    j0,
    nj
  )
  run_pair(
    "gemm_block",
    case,
    pair_env$pair_gemm_block(kind, A, B, j0, nj, C0),
    simple("gemm_block", case)
  )
  count_set("gemm_block", 1L)
}
for (spec in list(
  c(6, 9, 4, 3, 6),
  c(20, 50, 3, 12, 38),
  c(20, 50, 3, 0, 50)
)) {
  ## C.rows(j0 ..) = A.cols(j0 ..)' B : A is r x c, B is r x q, C is c x q.
  gemm_block_case(
    "tn_rows",
    rmat(spec[[1]], spec[[2]]),
    rmat(spec[[1]], spec[[3]]),
    spec[[4]],
    spec[[5]],
    rmat(spec[[2]], spec[[3]])
  )
  ## C.cols(j0 ..) = A B.rows(j0 ..)' : A is r x q, B is c x q, C is r x c.
  gemm_block_case(
    "nt_cols",
    rmat(spec[[1]], spec[[3]]),
    rmat(spec[[2]], spec[[3]]),
    spec[[4]],
    spec[[5]],
    rmat(spec[[1]], spec[[2]])
  )
}

## gemv: A is m x n.
gemv_case <- function(tr, m, n, cf) {
  A <- rmat(m, n)
  len_x <- if (tr == "N") n else m
  len_y <- if (tr == "N") m else n
  case <- sprintf("%s %dx%d alpha %g beta %g", tr, m, n, cf[[1]], cf[[2]])
  run_pair(
    "gemv",
    case,
    pair_env$pair_gemv(tr, A, rnorm(len_x), cf[[1]], cf[[2]], rnorm(len_y)),
    simple("gemv", case)
  )
  count_set("gemv", 1L)
}
for (tr in c("N", "T")) {
  for (sh in list(
    c(7, 5),
    c(1, 5),
    c(7, 1),
    c(1, 1),
    c(4, 4),
    c(5, 5),
    c(300, 200)
  )) {
    for (cf in coefs) {
      gemv_case(tr, sh[[1]], sh[[2]], cf)
    }
  }
  ## an empty x, compared with the Armadillo product; an empty y.
  if (tr == "N") {
    gemv_case(tr, 5, 0, c(1, 0))
    gemv_case(tr, 0, 5, c(1, 0))
  } else {
    gemv_case(tr, 0, 5, c(1, 0))
    gemv_case(tr, 5, 0, c(1, 0))
  }
}

## syrk: C is n x n; A is n x k for 'N' and k x n for 'T'.
syrk_case <- function(ul, tr, n, k, cf, A = NULL, tag = "") {
  if (is.null(A)) {
    A <- if (tr == "N") rmat(n, k) else rmat(k, n)
  }
  case <- sprintf(
    "%s%s n %d k %d alpha %g beta %g%s",
    ul,
    tr,
    n,
    k,
    cf[[1]],
    cf[[2]],
    tag
  )
  run_pair(
    "syrk",
    case,
    pair_env$pair_syrk(ul, tr, A, cf[[1]], cf[[2]], rmat(n, n)),
    simple("syrk", case)
  )
  count_set("syrk", 1L)
}
for (ul in c("U", "L")) {
  for (tr in c("N", "T")) {
    for (sh in list(
      c(6, 4),
      c(6, 1),
      c(1, 4),
      c(1, 1),
      c(5, 9),
      c(40, 25),
      c(150, 220),
      c(6, 0),
      c(0, 4)
    )) {
      for (cf in list(c(1, 0), c(1, 1), c(0.5, 2))) {
        syrk_case(ul, tr, sh[[1]], sh[[2]], cf)
      }
    }
    dims <- if (tr == "N") c(8, 5) else c(5, 8)
    dup <- rmat(dims[[1]], dims[[2]])
    if (tr == "N") {
      dup[2, ] <- dup[1, ]
    } else {
      dup[, 2] <- dup[, 1]
    }
    syrk_case(ul, tr, 8, 5, c(1, 0), dup, " rank-deficient A")
    syrk_case(ul, tr, 8, 5, c(1, 1), dup * 0, " zero A")
  }
}
for (spec in list(c(30, 12, 5), c(30, 12, 0), c(200, 40, 13))) {
  case <- sprintf("X %dx%d from column %d", spec[[1]], spec[[2]], spec[[3]])
  run_pair(
    "syrk_block",
    case,
    pair_env$pair_syrk_block(rmat(spec[[1]], spec[[2]]), spec[[3]]),
    simple("syrk_block", case)
  )
  count_set("syrk_block", 1L)
}

## matmul against operator*. Which cases Armadillo runs through its own
## loops, as READ in armadillo_bits (glue_times::apply, gemm, gemv): the
## routine is chosen in this order, and the tiny rule is tested on the matrix
## operand of that routine.
arma_tiny <- function(trans_a, dim_a, dim_b) {
  c_rows <- if (trans_a) dim_a[[2]] else dim_a[[1]]
  tiny_sq <- function(d) d[[1]] == d[[2]] && d[[1]] <= 4
  if (c_rows == 1) {
    return(tiny_sq(dim_b))
  }
  if (dim_b[[2]] == 1) {
    return(tiny_sq(dim_a))
  }
  tiny_sq(dim_a) && tiny_sq(dim_b) && dim_a[[1]] == dim_b[[1]]
}
matmul_case <- function(trans_a, mode, A, B, tag = "") {
  c_rows <- if (trans_a) ncol(A) else nrow(A)
  empty <- length(A) == 0L || length(B) == 0L
  routine <- if (!empty && arma_tiny(trans_a, dim(A), dim(B))) {
    "matmul_tiny"
  } else {
    "matmul"
  }
  case <- sprintf(
    "%s %s A %dx%d B %dx%d%s",
    if (trans_a) "A'B" else "AB",
    mode,
    nrow(A),
    ncol(A),
    nrow(B),
    ncol(B),
    tag
  )
  run_pair(
    routine,
    case,
    pair_env$pair_matmul(trans_a, mode, A, B, rmat(c_rows, ncol(B))),
    simple(routine, case, if (routine == "matmul_tiny") "tiny" else "same")
  )
  count_set(routine, 1L)
}
## (rows of A, cols of A, cols of B); B has as many rows as op(A) has cols.
matmul_shapes <- list(
  plain = list(
    c(8, 3, 5),
    c(8, 3, 1),
    c(1, 6, 4),
    c(1, 1, 1),
    c(1, 5, 1),
    c(2, 2, 2),
    c(3, 3, 3),
    c(4, 4, 4),
    c(5, 5, 5),
    c(4, 4, 2),
    c(3, 3, 1),
    c(4, 4, 1),
    c(5, 5, 1),
    c(1, 3, 3),
    c(1, 5, 5),
    c(50, 50, 1),
    c(50, 50, 3),
    c(50, 3, 3),
    c(120, 12, 12),
    c(200, 200, 8),
    c(0, 3, 2),
    c(4, 3, 0)
  ),
  trans = list(
    c(8, 3, 5),
    c(8, 1, 5),
    c(8, 3, 1),
    c(1, 1, 1),
    c(5, 1, 1),
    c(2, 2, 2),
    c(3, 3, 3),
    c(4, 4, 4),
    c(5, 5, 5),
    c(3, 3, 1),
    c(3, 1, 3),
    c(4, 1, 4),
    c(50, 3, 3),
    c(50, 1, 1),
    c(200, 8, 8),
    c(4, 0, 2),
    c(4, 3, 0)
  )
)
for (trans_a in c(FALSE, TRUE)) {
  for (sh in matmul_shapes[[if (trans_a) "trans" else "plain"]]) {
    inner <- if (trans_a) sh[[1]] else sh[[2]]
    for (mode in c("set", "minus", "plus")) {
      matmul_case(trans_a, mode, rmat(sh[[1]], sh[[2]]), rmat(inner, sh[[3]]))
    }
  }
  ## An empty inner dimension: only as an assignment. Armadillo's in-place
  ## forms pass it on to dgemm with a leading dimension of 0, which BLAS
  ## refuses through xerbla.
  if (trans_a) {
    matmul_case(trans_a, "set", rmat(0, 3), rmat(0, 2))
  } else {
    matmul_case(trans_a, "set", rmat(4, 0), rmat(0, 3))
  }
  dup <- rmat(9, 6)
  dup[, 2] <- dup[, 1]
  holed <- rmat(9, 6)
  holed[2, 3] <- NaN
  inner <- if (trans_a) 9 else 6
  matmul_case(trans_a, "set", dup, rmat(inner, 4), " rank-deficient A")
  matmul_case(trans_a, "minus", dup * 0, rmat(inner, 4), " zero A")
  matmul_case(trans_a, "set", holed, rmat(inner, 4), " NaN in A")
}
## Dimensions that do not conform: both sides must refuse.
run_pair(
  "matmul",
  "AB set A 4x3 B 2x2, not conformable",
  pair_env$pair_matmul(FALSE, "set", rmat(4, 3), rmat(2, 2), rmat(4, 2)),
  simple("matmul", "AB set A 4x3 B 2x2, not conformable")
)
count_set("matmul", 1L)

## syevr_range against eig_sym(). The triangle that is not named holds
## unrelated finite values, so reading it changes the answer. The eigenvalues
## are compared, and Z diag(w) Z', which does not depend on the basis chosen
## inside a repeated eigenvalue as long as the range does not cut one.
syevr_case <- function(ul, S, il, iu, tag) {
  n <- nrow(S)
  A <- S
  junk <- rmat(n, n)
  other <- if (ul == "U") lower.tri(A) else upper.tri(A)
  A[other] <- junk[other]
  case <- sprintf("%s n %d pairs %d..%d %s", ul, n, il, iu, tag)
  run_pair(
    "syevr_range",
    case,
    pair_env$pair_syevr(ul, A, il, iu),
    function(out) {
      recon <- function(z, w) z %*% (w * t(z))
      failed_w <- out$wrapper$status != codes[["ok"]]
      record(
        "syevr_range",
        case,
        failed_w,
        out$arma$failed,
        w = if (failed_w) {
          list()
        } else {
          list(
            values = as.vector(out$wrapper$values),
            rebuilt = recon(out$wrapper$vectors, as.vector(out$wrapper$values))
          )
        },
        a = if (out$arma$failed) {
          list()
        } else {
          list(
            values = as.vector(out$arma$values),
            rebuilt = recon(out$arma$vectors, as.vector(out$arma$values))
          )
        },
        class = "eig",
        extra = c(
          "the reference was built from the named triangle alone" = identical(
            out$arma$sym,
            S
          )
        )
      )
    }
  )
  count_set("syevr_range", 1L)
}
gram <- function(m, n) crossprod(rmat(m, n)) / m
for (ul in c("U", "L")) {
  syevr_case(ul, gram(5, 1), 1, 1, "one by one")
  syevr_case(ul, gram(9, 2), 1, 2, "all")
  syevr_case(ul, gram(9, 2), 2, 2, "largest")
  syevr_case(ul, gram(20, 3), 2, 3, "top 2")
  syevr_case(ul, gram(40, 10), 8, 10, "top 3")
  syevr_case(ul, gram(40, 10), 1, 10, "all")
  syevr_case(ul, gram(40, 10), 1, 2, "bottom 2")
  syevr_case(ul, gram(40, 10), 4, 6, "middle")
  syevr_case(ul, gram(200, 60), 56, 60, "top 5")
  syevr_case(ul, gram(500, 300), 289, 300, "top 12")
  ## rank 4 of 10: the six zero eigenvalues are below the range, or all in it.
  syevr_case(ul, gram(4, 10), 7, 10, "rank-deficient, top 4")
  syevr_case(ul, gram(4, 10), 1, 10, "rank-deficient, all")
  syevr_case(ul, diag(6), 1, 6, "identity, all")
  syevr_case(ul, matrix(0, 6, 6), 1, 6, "zero, all")
  syevr_case(ul, diag(c(1, 1, 2, 2, 5, 5)), 5, 6, "repeated, top pair")
  indef <- rmat(12, 12)
  syevr_case(ul, (indef + t(indef)) / 2, 10, 12, "indefinite, top 3")
  syevr_case(ul, 100 * gram(30, 8), 6, 8, "scale 100, top 3")
}

## syevd against arma::dsyevd_, with the other triangle holding unrelated
## values as above.
syevd_case <- function(ul, S, extra, tag) {
  n <- nrow(S)
  A <- S
  junk <- rmat(n, n)
  other <- if (ul == "U") lower.tri(A) else upper.tri(A)
  A[other] <- junk[other]
  case <- sprintf("%s n %d workspace +%d %s", ul, n, extra, tag)
  run_pair("syevd", case, pair_env$pair_syevd(ul, A, extra), function(out) {
    record(
      "syevd",
      case,
      out$wrapper$status != codes[["ok"]],
      out$arma$failed,
      w = list(
        values = as.vector(out$wrapper$values),
        vectors = out$wrapper$vectors
      ),
      a = list(
        values = as.vector(out$arma$values),
        vectors = out$arma$vectors
      )
    )
  })
  count_set("syevd", 1L)
}
for (ul in c("U", "L")) {
  for (extra in c(0L, 500L)) {
    for (n in c(1, 2, 3, 6, 30, 90, 200)) {
      sym <- rmat(n, n)
      syevd_case(ul, (sym + t(sym)) / 2, extra, "symmetric")
    }
    syevd_case(ul, gram(4, 10), extra, "rank-deficient")
    syevd_case(ul, diag(6), extra, "identity")
    syevd_case(ul, matrix(0, 6, 6), extra, "zero")
    syevd_case(ul, diag(c(1, 1, 2, 2, 5, 5)), extra, "repeated")
  }
  syevd_case(ul, matrix(0, 0, 0), 0L, "empty")
}

## qr_econ against arma::qr_econ().
qr_case <- function(X, tag) {
  case <- sprintf("%dx%d %s", nrow(X), ncol(X), tag)
  run_pair("qr_econ", case, pair_env$pair_qr_econ(X), function(out) {
    record(
      "qr_econ",
      case,
      out$wrapper$status != codes[["ok"]],
      out$arma$failed,
      w = list(q = out$wrapper$q, r = out$wrapper$r),
      a = list(q = out$arma$q, r = out$arma$r)
    )
  })
  count_set("qr_econ", 1L)
}
qr_case(rmat(50, 5), "tall")
qr_case(rmat(20, 1), "one column")
qr_case(rmat(5, 5), "square")
qr_case(rmat(1, 1), "one by one")
qr_case(rmat(3, 7), "wide")
qr_case(rmat(1, 5), "one row")
qr_case(rmat(400, 30), "tall, larger")
qr_case(rmat(6, 5), "tall by one row")
dup <- rmat(12, 6)
dup[, 4] <- dup[, 2]
qr_case(dup, "rank-deficient")
dup <- rmat(6, 6)
dup[, 6] <- dup[, 1]
qr_case(dup, "square, rank-deficient")
zero_first <- rmat(8, 3)
zero_first[, 1] <- 0
qr_case(zero_first, "first column zero")
qr_case(matrix(0, 8, 3), "all zero")
qr_case(
  sweep(rmat(30, 6), 2, 10^seq(-8, 8, length.out = 6), "*"),
  "scaled columns"
)
qr_case(matrix(0, 0, 3), "no rows")
qr_case(matrix(0, 3, 0), "no columns")
qr_case(matrix(0, 0, 0), "empty")
holed <- rmat(8, 3)
holed[5, 2] <- NaN
qr_case(holed, "NaN")

## solve_tri against solve(trimatu / trimatl). The other triangle holds
## unrelated finite values. `path` is what the case was built to make
## Armadillo do, and the case fails if it did something else: "exact",
## "approx" (its warning, then the least-squares answer), "failed", or
## "empty" (no rows, or no right-hand side).
tri_case <- function(ul, Tm, B, tag, path) {
  n <- nrow(Tm)
  A <- Tm
  if (n > 0) {
    junk <- rmat(n, n)
    other <- if (ul == "U") lower.tri(A) else upper.tri(A)
    A[other] <- junk[other]
  }
  case <- sprintf("%s n %d nrhs %d %s", ul, n, ncol(B), tag)
  run_pair("solve_tri", case, pair_env$pair_solve_tri(ul, A, B), function(out) {
    failed_w <- out$wrapper$status != codes[["ok"]]
    took <- if (out$arma$failed) {
      "failed"
    } else if (path == "empty") {
      "empty"
    } else if (out$arma$approx) {
      "approx"
    } else {
      "exact"
    }
    extra <- logical()
    extra[[sprintf("built for the %s path, Armadillo took %s", path, took)]] <-
      identical(took, path)
    w <- list(value = out$wrapper$value)
    a <- list(value = out$arma$value)
    if (path == "empty") {
      ## Armadillo's estimate stays 0 for an empty system, so it warns and
      ## goes through its fallback, which returns the empty result. The
      ## wrapper returns the empty result and does not call it approximate.
      extra[["empty system: Armadillo warns, the wrapper sets no flag"]] <-
        isTRUE(out$arma$approx) && !isTRUE(out$wrapper$approx)
    } else {
      if (!failed_w) {
        extra[["approx flag against Armadillo's warning"]] <-
          identical(out$wrapper$approx, out$arma$approx)
      }
      w$rcond <- out$wrapper$rcond
      a$rcond <- out$arma$rcond
    }
    record("solve_tri", case, failed_w, out$arma$failed, w, a, extra = extra)
    tri_paths[[length(tri_paths) + 1L]] <<- took
  })
  count_set("solve_tri", 1L)
}
tri_paths <- list()
## An n x n triangle from a QR factor. `d` replaces its diagonal, which with
## the off-diagonal left as it is makes the inverse grow far faster than
## 1 / min(d). `rows` scales whole rows instead, so the condition number
## grows only like max(rows) / min(rows).
make_tri <- function(ul, n, d = NULL, rows = NULL) {
  r <- qr.R(qr(rmat(n + 3, n)))
  if (!is.null(rows)) {
    r <- r * rows
  }
  if (!is.null(d)) {
    diag(r) <- d
  }
  if (ul == "U") r else t(r)
}
for (ul in c("U", "L")) {
  for (sh in list(c(1, 1), c(3, 3), c(10, 4), c(10, 1), c(40, 40), c(12, 30))) {
    tri_case(
      ul,
      make_tri(ul, sh[[1]]),
      rmat(sh[[1]], sh[[2]]),
      "ordinary",
      "exact"
    )
  }
  n <- 8
  tri_case(
    ul,
    make_tri(ul, n, rows = 10^seq(0, -8, length.out = n)),
    rmat(n, 3),
    "ill-conditioned, above eps",
    "exact"
  )
  tri_case(
    ul,
    make_tri(ul, n, 10^seq(0, -10, length.out = n)),
    rmat(n, 3),
    "ill-conditioned, below eps",
    "approx"
  )
  tri_case(
    ul,
    make_tri(ul, n, c(rep(1, n - 1), 1e-20)),
    rmat(n, 3),
    "near-singular",
    "approx"
  )
  for (z in c(1, 4, n)) {
    d <- rep(1, n)
    d[[z]] <- 0
    tri_case(
      ul,
      make_tri(ul, n, d),
      rmat(n, 3),
      sprintf("zero at diagonal %d", z),
      "approx"
    )
  }
  tri_case(ul, matrix(0, n, n), rmat(n, 2), "all zero", "approx")
  holed <- rmat(n, 3)
  holed[2, 2] <- NaN
  tri_case(ul, make_tri(ul, n), holed, "NaN in B, sound A", "exact")
  d <- rep(1, n)
  d[[3]] <- NaN
  tri_case(ul, make_tri(ul, n, d), rmat(n, 3), "NaN on the diagonal", "failed")
  d <- rep(1, n)
  d[[n]] <- 0
  holed <- rmat(n, 3)
  holed[1, 1] <- Inf
  tri_case(ul, make_tri(ul, n, d), holed, "singular A, Inf in B", "failed")
  tri_case(ul, matrix(0, 0, 0), matrix(0, 0, 3), "no rows", "empty")
  tri_case(ul, make_tri(ul, 3), matrix(0, 3, 0), "no right-hand side", "empty")
}

## The solve as src/lobpcg_warm.h writes it: R is k x k upper triangular, AQ
## is n x k.
chain_case <- function(n, k, d = NULL, tag = "ordinary") {
  R <- qr.R(qr(rmat(n, k)))
  if (!is.null(d)) {
    diag(R) <- d
  }
  case <- sprintf("AQ %dx%d %s", n, k, tag)
  run_pair(
    "solve_tri_chain",
    case,
    pair_env$pair_solve_chain(R, rmat(n, k)),
    function(out) {
      failed_w <- out$wrapper$status != codes[["ok"]]
      record(
        "solve_tri_chain",
        case,
        failed_w,
        out$arma$failed,
        w = list(value = out$wrapper$value),
        a = list(value = out$arma$value),
        extra = c(
          "approx flag against Armadillo's warning" = identical(
            out$wrapper$approx,
            out$arma$approx
          ) ||
            failed_w
        )
      )
    }
  )
  count_set("solve_tri_chain", 1L)
}
chain_case(50, 3)
chain_case(200, 12)
chain_case(20, 1)
chain_case(60, 6, c(1, 1, 1, 1, 1, 1e-20), "near-singular R")

## svd_values against svd().
svd_case <- function(X, tag, expect_failed = FALSE) {
  case <- sprintf("%dx%d %s", nrow(X), ncol(X), tag)
  run_pair("svd_values", case, pair_env$pair_svd(X), function(out) {
    failed_w <- out$wrapper$status != codes[["ok"]]
    record(
      "svd_values",
      case,
      failed_w,
      out$arma$failed,
      w = list(value = as.vector(out$wrapper$value)),
      a = list(value = as.vector(out$arma$value)),
      extra = c(
        "the path the case was built for (failure)" = identical(
          out$arma$failed,
          expect_failed
        )
      )
    )
    svd_failed[[length(svd_failed) + 1L]] <<- out$arma$failed
  })
  count_set("svd_values", 1L)
}
svd_failed <- list()
svd_case(rmat(30, 5), "tall")
svd_case(rmat(5, 30), "wide")
svd_case(rmat(6, 6), "square")
svd_case(rmat(20, 1), "one column")
svd_case(rmat(1, 20), "one row")
svd_case(rmat(1, 1), "one by one")
svd_case(rmat(40, 30), "1200 cells, workspace asked for")
svd_case(rmat(32, 32), "1024 cells, workspace asked for")
svd_case(rmat(31, 33), "1023 cells, minimum workspace")
dup <- rmat(12, 6)
dup[, 4] <- dup[, 2]
svd_case(dup, "rank-deficient")
svd_case(matrix(0, 7, 4), "all zero")
## as the diagnostics build calls it: V_prev' V_new of two orthonormal blocks.
v1 <- qr.Q(qr(rmat(60, 5)))
v2 <- qr.Q(qr(rmat(60, 5)))
svd_case(crossprod(v1, v2), "two orthonormal blocks")
holed <- rmat(9, 4)
holed[3, 2] <- NaN
svd_case(holed, "NaN", expect_failed = TRUE)
holed[3, 2] <- Inf
svd_case(holed, "Inf", expect_failed = TRUE)
svd_case(matrix(0, 0, 3), "no rows")
svd_case(matrix(0, 3, 0), "no columns")

## What the wrappers refuse on their own, before any Fortran routine.
expected_refusals <- c(
  gemm_trans = "err_arg",
  gemm_conj = "err_arg",
  gemm_lda = "err_arg",
  gemm_ldb = "err_arg",
  gemm_ldc = "err_arg",
  gemm_size = "err_size",
  gemv_trans = "err_arg",
  gemv_lda = "err_arg",
  gemv_size = "err_size",
  syrk_uplo = "err_arg",
  syrk_trans = "err_arg",
  syrk_lda = "err_arg",
  syrk_ldc = "err_arg",
  syrk_size = "err_size",
  matmul_inner = "err_arg",
  matmul_size = "err_size",
  syevr_uplo = "err_arg",
  syevr_il_zero = "err_arg",
  syevr_il_past_iu = "err_arg",
  syevr_iu_past_n = "err_arg",
  syevr_size = "err_size",
  syevr_work_size = "err_size",
  syevd_uplo = "err_arg",
  syevd_size = "err_size",
  syevd_work_size = "err_size",
  qr_size = "err_size",
  solve_uplo = "err_arg",
  solve_size = "err_size",
  svd_size = "err_size"
)
got <- pair_env$wrapper_refusals()
for (nm in union(names(expected_refusals), names(got))) {
  want <- unname(codes[expected_refusals[nm]])
  have <- unname(got[nm])
  record(
    "refusals",
    nm,
    FALSE,
    FALSE,
    extra = structure(
      isTRUE(have == want),
      names = sprintf(
        "returned %s, expected %s (%s)",
        have,
        want,
        expected_refusals[nm]
      )
    )
  )
  count_set("refusals", 1L)
}

# -- report ------------------------------------------------------------------

tab <- do.call(
  rbind,
  lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE))
)
routines <- unique(tab$routine)
if (verbose) {
  for (i in seq_len(nrow(tab))) {
    cat(sprintf(
      "%-16s %-58s abs %-9.3g rel %-9.3g %s%s\n",
      tab$routine[[i]],
      tab$case[[i]],
      tab$max_abs[[i]],
      tab$max_rel[[i]],
      if (tab$both_failed[[i]]) {
        "both failed"
      } else if (isTRUE(tab$bitwise[[i]])) {
        "bitwise"
      } else {
        "not bitwise"
      },
      if (tab$ok[[i]]) "" else paste("  FAIL:", tab$problems[[i]])
    ))
  }
  cat("\n")
}
for (cl in names(tol)) {
  cat(sprintf(
    "tolerance class %-4s: abs %g, rel %g over %s\n",
    cl,
    tol[[cl]]$abs,
    tol[[cl]]$rel,
    if (tol[[cl]]$rel_floor == 0) {
      "every cell either side holds as nonzero"
    } else {
      sprintf("the cells above %g of the largest", tol[[cl]]$rel_floor)
    }
  ))
}
cat(
  "abs is asserted against the tolerance times max(1, largest magnitude).\n\n"
)
failed_routines <- character()
for (rt in routines) {
  x <- tab[tab$routine == rt, , drop = FALSE]
  cmp <- x[!x$both_failed, , drop = FALSE]
  i_abs <- which.max(cmp$max_abs)
  i_rel <- which.max(cmp$max_rel)
  n_bit <- sum(cmp$bitwise %in% TRUE)
  cat(sprintf(
    "%-16s %s  class %s\n",
    rt,
    if (all(x$ok)) "ok  " else "FAIL",
    x$class[[1]]
  ))
  cat(sprintf(
    "  cases %d, of which both sides failed in %d; cells compared %d, finite %d\n",
    nrow(x),
    sum(x$both_failed),
    sum(cmp$n_cells),
    sum(cmp$n_finite)
  ))
  if (rt != "refusals" && nrow(cmp) > 0L) {
    cat(sprintf(
      "  worst abs %.3g  [%s]\n",
      cmp$max_abs[[i_abs]],
      cmp$case[[i_abs]]
    ))
    cat(sprintf(
      "  worst rel %.3g  [%s], relative taken over %d cells in all\n",
      cmp$max_rel[[i_rel]],
      cmp$case[[i_rel]],
      sum(cmp$n_rel)
    ))
    cat(sprintf(
      "  bitwise identical in %d of %d compared cases\n",
      n_bit,
      nrow(cmp)
    ))
  }
  bad <- x[!x$ok, , drop = FALSE]
  for (i in seq_len(nrow(bad))) {
    cat(sprintf("  FAIL [%s]: %s\n", bad$case[[i]], bad$problems[[i]]))
  }
  if (!all(x$ok)) {
    failed_routines <- c(failed_routines, rt)
  }
}

## A check that compared nothing passes silently, so each input set, and each
## path a set was built to take, must be there.
need <- function(what, n) {
  if (!isTRUE(n > 0)) {
    cat(sprintf("EMPTY: %s\n", what))
    failed_routines <<- c(failed_routines, paste0("empty:", what))
  }
}
cat("\ninput sets:\n")
for (nm in names(sets)) {
  x <- tab[tab$routine == nm, , drop = FALSE]
  cat(sprintf("  %-16s %d cases\n", nm, sets[[nm]]))
  need(paste(nm, "cases"), sets[[nm]])
  need(paste(nm, "cases recorded"), nrow(x) == sets[[nm]])
  if (nm != "refusals") {
    need(paste(nm, "finite cells compared"), sum(x$n_finite))
    need(paste(nm, "cells in the relative figure"), sum(x$n_rel))
  }
}
for (nm in c(
  "gemm",
  "gemm_block",
  "gemv",
  "syrk",
  "syrk_block",
  "matmul",
  "matmul_tiny",
  "syevr_range",
  "syevd",
  "qr_econ",
  "solve_tri",
  "solve_tri_chain",
  "svd_values",
  "refusals"
)) {
  need(paste(nm, "input set"), if (nm %in% names(sets)) sets[[nm]] else 0L)
}
took <- unlist(tri_paths)
cat(sprintf(
  "  solve_tri paths Armadillo took: exact %d, approximate fallback %d, failed %d, empty system %d\n",
  sum(took == "exact"),
  sum(took == "approx"),
  sum(took == "failed"),
  sum(took == "empty")
))
for (p in c("exact", "approx", "failed", "empty")) {
  need(paste("solve_tri cases on the", p, "path"), sum(took == p))
}
cat(sprintf(
  "  svd_values paths: Armadillo failed in %d cases\n",
  sum(unlist(svd_failed))
))
need("svd_values cases that fail", sum(unlist(svd_failed)))
need(
  "matmul cases both sides refuse",
  sum(tab$routine == "matmul" & tab$both_failed)
)

cat("\n")
if (length(failed_routines)) {
  cat("ORACLE FAILED:", paste(unique(failed_routines), collapse = ", "), "\n")
  quit(status = 1L)
}
cat("ORACLE PASSED\n")
