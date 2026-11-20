# Step 4 of Phase 2 (dev/eigen_conversion.md at 4f2b669): does dsyevd return the
# same bits whatever workspace length it is given above its own answer?
#
# Before step 4, small_eig_desc() (src/lobpcg_warm.h) kept its dsyevd
# workspace on the LOBPCGState and never shrank it, so a Rayleigh-Ritz
# problem of order n that followed a larger one, of order m, was solved with
# the workspace of order m. si_la::syevd asks for the workspace of order n on
# every call. This script measures whether that changes the eigenpairs.
#
# It compiles dev/oracle_pairs.cpp as dev/oracle_run.R does and calls its
# pair_syevd(), which runs si_la::syevd beside dsyevd with `extra` cells
# added to the queried workspace. Here `extra` is the difference between the
# documented workspace of order m, 1 + 6 m + 2 m^2, and that of order n, for
# m = n (the control: must be bitwise), n + 1, 1.5 n, 2 n and 3 n.
#
# Dev only, run by hand, from the package root:
#   Rscript dev/probe_p2s4_syevd_lwork.R <lib> [<n max>]
# <lib> is the private library dev/oracle_install_armadillo.R installed
# RcppArmadillo into (dev/lib/arma).
#
# The oracle and dev/lib/arma were removed after 141a9d0, and from 97525a2
# on the oracle does not compile: this runs in a clone of an earlier commit,
# which has the oracle's files (git show 141a9d0:<path> prints each).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    "usage: Rscript dev/probe_p2s4_syevd_lwork.R <lib> [<n max>]",
    call. = FALSE
  )
}
lib <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
n_max <- if (length(args) >= 2L) as.integer(args[[2]]) else 140L
.libPaths(c(lib, .libPaths()))
stopifnot(identical(
  normalizePath(dirname(find.package("RcppArmadillo")), winslash = "/"),
  lib
))
cat("R:", R.version.string, "\n")
cat("LAPACK:", sessionInfo()$LAPACK, "version", La_version(), "\n")
cat("BLAS:", sessionInfo()$BLAS, "\n")

build <- file.path(tempdir(), "probe-syevd-lwork")
unlink(build, recursive = TRUE)
dir.create(build, recursive = TRUE)
stopifnot(all(file.copy(
  c("dev/oracle_pairs.cpp", "src/si_lapack.h", "src/si_lapack.cpp"),
  build
)))
pair_env <- new.env()
Rcpp::sourceCpp(
  file.path(build, "oracle_pairs.cpp"),
  env = pair_env,
  cacheDir = file.path(build, "cache"),
  rebuild = TRUE,
  showOutput = FALSE
)

lwmin <- function(n) 1 + 6 * n + 2 * n^2
bitwise <- function(a, b) {
  identical(
    writeBin(as.vector(a), raw()),
    writeBin(as.vector(b), raw())
  )
}

set.seed(1)
factors <- c(same = 1, plus1 = NA, x1.5 = 1.5, x2 = 2, x3 = 3)
rows <- list()
for (n in 2:n_max) {
  ## Two kinds of matrix: a random symmetric one, and the kind the
  ## Rayleigh-Ritz step solves, a projection of a Gram onto nearly
  ## orthonormal columns, with a decaying spectrum.
  s <- matrix(rnorm(n * n), n, n)
  mats <- list(
    sym = (s + t(s)) / 2,
    gram = {
      q <- qr.Q(qr(matrix(rnorm(n * n), n, n)))
      q %*% (c(100 / seq_len(n)) * t(q))
    }
  )
  for (kind in names(mats)) {
    A <- mats[[kind]]
    A <- (A + t(A)) / 2
    for (f in names(factors)) {
      m <- if (f == "plus1") n + 1 else ceiling(factors[[f]] * n)
      extra <- as.integer(lwmin(m) - lwmin(n))
      out <- pair_env$pair_syevd("U", A, extra)
      stopifnot(out$wrapper$status == 0L, !out$arma$failed)
      dv <- abs(out$wrapper$vectors - out$arma$vectors)
      dw <- abs(as.vector(out$wrapper$values) - as.vector(out$arma$values))
      rows[[length(rows) + 1L]] <- data.frame(
        n = n,
        kind = kind,
        buffer = f,
        m = m,
        bitwise = bitwise(out$wrapper$vectors, out$arma$vectors) &&
          bitwise(out$wrapper$values, out$arma$values),
        abs_vec = max(dv),
        abs_val = max(dw),
        rel_val = max(dw / abs(as.vector(out$arma$values)))
      )
    }
  }
}
res <- do.call(rbind, rows)
stopifnot(nrow(res) == (n_max - 1L) * 2L * length(factors))

cat(sprintf("\n%d calls, orders 2 to %d\n", nrow(res), n_max))
ctl <- res[res$buffer == "same", ]
cat(sprintf(
  "control (workspace of the same order): %d of %d bitwise\n",
  sum(ctl$bitwise),
  nrow(ctl)
))
for (f in setdiff(names(factors), "same")) {
  r <- res[res$buffer == f, ]
  off <- r[!r$bitwise, ]
  cat(sprintf(
    "\nworkspace of order %s: %d of %d calls not bitwise\n",
    f,
    nrow(off),
    nrow(r)
  ))
  if (nrow(off)) {
    cat(sprintf(
      "  orders that differ: %d to %d (%d distinct)\n",
      min(off$n),
      max(off$n),
      length(unique(off$n))
    ))
    cat(sprintf(
      "  orders in that span that do not: %s\n",
      paste(
        setdiff(seq(min(off$n), max(off$n)), unique(off$n)),
        collapse = " "
      )
    ))
    w <- off[which.max(off$abs_vec), ]
    cat(sprintf(
      "  worst abs, vectors: %.3g (n %d, %s); values: abs %.3g, rel %.3g, the relative figure over every eigenvalue of each call\n",
      w$abs_vec,
      w$n,
      w$kind,
      max(off$abs_val),
      max(off$rel_val)
    ))
  }
}
ok <- all(ctl$bitwise)
cat(if (ok) "\nCONTROL OK\n" else "\nCONTROL FAILED\n")
quit(status = if (ok) 0L else 1L)
