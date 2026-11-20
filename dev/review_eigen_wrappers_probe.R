# Review of the wrapper layer and the oracle (dev/review_eigen_wrappers.md):
# direct probes of src/si_lapack.h that dev/oracle_run.R does not make. The
# oracle and dev/lib/arma were removed after 141a9d0; the commit under
# review, 79cce0f, has the oracle's files (git show 79cce0f:<path>). The
# probe calls si_la::matmul and si_la::solve_tri with its two outputs, both
# gone since Phase 4, so it builds only in a clone of 79cce0f, not against
# src/ as it is now. Dev only, run by hand, from the root of that clone:
#   Rscript dev/review_eigen_wrappers_probe.R <lib> <out dir> [--group=<name>]
# <lib>      the private RcppArmadillo library, as for dev/oracle_run.R.
# <out dir>  where the sources are staged and built; nothing is built in the
#            working tree, and src/ and the oracle are only read.
# Without --group every group below runs in a child R process of its own,
# with a time limit, because a group may hang or kill R (NaN handed to
# LAPACK); the child's whole output is printed. nonfinite_syevr runs one
# child per case, each limited to 60 s. --group=<name> runs that group alone.
# Groups:
#   w_length         how far past iu - il + 1 dsyevr writes in w
#   refusals         refusals the oracle does not try
#   fuzz             every wrapper over small dimensions, 0 included
#   order            whether a larger 1 x 1 product tells the dgemv routes apart
#   fallback_large   solve_tri's dgelsd fallback above n = 8, with Armadillo
#   nonfinite_syevd  NaN and Inf handed to si_la::syevd
#   nonfinite_syevr  NaN and Inf handed to si_la::syevr_range, with Armadillo
# Nothing here asserts: it prints what was measured.

args <- commandArgs(trailingOnly = TRUE)
opts <- grep("^--", args, value = TRUE)
pos <- setdiff(args, opts)
if (length(pos) != 2L) {
  stop(
    "usage: Rscript dev/review_eigen_wrappers_probe.R <lib> <out dir> ",
    "[--group=<name>]",
    call. = FALSE
  )
}
lib <- normalizePath(pos[[1]], winslash = "/", mustWork = TRUE)
out_dir <- pos[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
group <- sub("^--group=", "", grep("^--group=", opts, value = TRUE))
case_id <- as.integer(sub("^--case=", "", grep("^--case=", opts, value = TRUE)))
if (!file.exists("dev/review_eigen_wrappers_probe.cpp")) {
  stop("run from the package root", call. = FALSE)
}
groups <- c(
  "w_length",
  "refusals",
  "fuzz",
  "order",
  "fallback_large",
  "nonfinite_syevd",
  "nonfinite_syevr"
)

# -- parent: one child per group ---------------------------------------------

## nonfinite_syevr runs one case per child, 20 cases: a case that hangs must
## not take the others with it.
n_syevr_cases <- 20L
if (!length(group) || (group[[1]] == "nonfinite_syevr" && !length(case_id))) {
  rscript <- file.path(R.home("bin"), "Rscript")
  todo <- if (length(group)) group[[1]] else groups
  for (g in todo) {
    ids <- if (g == "nonfinite_syevr") seq_len(n_syevr_cases) else NA
    cat(sprintf("==== group %s\n", g))
    for (id in ids) {
      out <- suppressWarnings(system2(
        rscript,
        c(
          "dev/review_eigen_wrappers_probe.R",
          shQuote(lib),
          shQuote(out_dir),
          paste0("--group=", g),
          if (!is.na(id)) paste0("--case=", id)
        ),
        stdout = TRUE,
        stderr = TRUE,
        timeout = if (is.na(id)) 300 else 60
      ))
      status <- attr(out, "status")
      if (is.null(status)) {
        status <- 0L
      }
      cat(out, sep = "\n")
      cat(sprintf(
        "---- exit %d%s; lines that look like xerbla: %d\n",
        status,
        if (status == 124L) {
          " (the time limit: the call did not return)"
        } else {
          ""
        },
        length(grep("illegal value|On entry to|gave error code", out))
      ))
    }
    cat("\n")
  }
  quit(status = 0L)
}

# -- child: build, then one group --------------------------------------------

group <- group[[1]]
.libPaths(c(lib, .libPaths()))
## The link line of src/Makevars, for the unit that has no RcppArmadillo.
Rcpp::registerPlugin("review_lapack", function() {
  list(env = list(PKG_LIBS = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"))
})
stage <- function(cpp, tag) {
  build <- file.path(out_dir, paste0("probe-build-", tag))
  dir.create(build, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(
    c(cpp, "src/si_lapack.h", "src/si_lapack.cpp"),
    build,
    overwrite = TRUE
  )
  if (!all(ok)) {
    stop("could not stage ", cpp, call. = FALSE)
  }
  env <- new.env()
  Rcpp::sourceCpp(
    file.path(build, basename(cpp)),
    env = env,
    cacheDir = file.path(build, "cache"),
    rebuild = FALSE,
    showOutput = FALSE
  )
  env
}
say <- function(...) {
  cat(sprintf(...), "\n", sep = "")
  flush(stdout())
}
set.seed(20261118)
rmat <- function(m, n) matrix(rnorm(m * n), m, n)
gram <- function(m, n) crossprod(rmat(m, n)) / m
show_nonfinite <- function(x) {
  sprintf(
    "%d NaN, %d Inf of %d",
    sum(is.nan(x)),
    sum(is.infinite(x)),
    length(x)
  )
}

if (group == "w_length") {
  probe <- stage("dev/review_eigen_wrappers_probe.cpp", "plain")
  say("R %s, LAPACK %s", R.version.string, La_version())
  say(
    "w is handed n cells; the header promises only the first iu - il + 1 are written."
  )
  say(
    "'past' = how many of the other cells the call changed; err = max abs difference"
  )
  say(
    "of the returned values from eigen(); res = max abs of S z - z diag(w).\n"
  )
  w_case <- function(tag, S, il, iu, ul = "U") {
    n <- nrow(S)
    out <- probe$probe_syevr(ul, S, il, iu)
    ref <- rev(eigen(S, symmetric = TRUE, only.values = TRUE)$values)[il:iu]
    res <- if (out$status == 0L) {
      max(abs(
        S %*% out$vectors - out$vectors %*% diag(out$values, iu - il + 1)
      ))
    } else {
      NA_real_
    }
    say(
      "%-46s %s n %3d pairs %3d..%-3d status %d  past %2d of %2d  err %.2e  res %.2e",
      tag,
      ul,
      n,
      il,
      iu,
      out$status,
      out$changed_past_count,
      n - (iu - il + 1),
      max(abs(out$values - ref)),
      res
    )
    if (out$changed_past_count > 0) {
      hit <- out$past_count[out$past_count != -777.25]
      say(
        "    written past w[count - 1]: %s",
        paste(signif(hit, 6), collapse = " ")
      )
    }
    invisible(out$changed_past_count)
  }
  total <- 0L
  say("-- the range does not cut a group of equal eigenvalues")
  total <- total + w_case("random Gram, top 3", gram(40, 10), 8, 10)
  total <- total + w_case("random Gram, top 5", gram(200, 60), 56, 60)
  total <- total + w_case("random Gram, top 1", gram(200, 60), 60, 60)
  d6 <- diag(c(1, 1, 2, 2, 5, 5))
  total <- total +
    w_case("diag(1,1,2,2,5,5), top pair (the oracle's)", d6, 5, 6)
  g4 <- gram(4, 10)
  total <- total +
    w_case("rank 4 of 10, top 4 (the oracle's)", g4, 7, 10)
  say("\n-- the range cuts a group of equal eigenvalues")
  for (ul in c("U", "L")) {
    total <- total + w_case("identity, top 2", diag(6), 5, 6, ul)
    total <- total + w_case("diag(1,1,2,2,5,5), top 1", d6, 6, 6, ul)
    total <- total + w_case("diag(1,1,2,2,5,5), top 3", d6, 4, 6, ul)
    total <- total + w_case("rank 4 of 10, top 6", g4, 5, 10, ul)
    total <- total + w_case("rank 4 of 10, top 5", g4, 6, 10, ul)
  }
  q <- qr.Q(qr(rmat(6, 6)))
  dense <- q %*% diag(c(1, 1, 1, 2, 3, 3)) %*% t(q)
  dense <- (dense + t(dense)) / 2
  total <- total + w_case("dense, eigenvalues 1,1,1,2,3,3, top 1", dense, 6, 6)
  total <- total + w_case("dense, eigenvalues 1,1,1,2,3,3, top 5", dense, 2, 6)
  say("\n-- what pca_imp() can hand it: more components than the data has rank")
  ## 30 samples of which 6 are distinct, 100 columns, columns centred: the
  ## 30 x 30 Gram has rank 5 and 25 eigenvalues at rounding level.
  base <- rmat(6, 100)
  x <- base[rep(1:6, each = 5), ]
  x <- sweep(x, 2, colMeans(x))
  total <- total + w_case("6 distinct rows of 30, top 5", tcrossprod(x), 26, 30)
  total <- total + w_case("6 distinct rows of 30, top 6", tcrossprod(x), 25, 30)
  total <- total + w_case("6 distinct rows of 30, top 8", tcrossprod(x), 23, 30)
  ## an exactly rank-3 matrix, 20 x 50
  lr <- rmat(20, 3) %*% t(rmat(50, 3))
  total <- total + w_case("rank 3, 20 x 50, top 3", tcrossprod(lr), 18, 20)
  total <- total + w_case("rank 3, 20 x 50, top 5", tcrossprod(lr), 16, 20)
  say("\ncells written past w[count - 1], all cases: %d", total)
}

if (group == "refusals") {
  probe <- stage("dev/review_eigen_wrappers_probe.cpp", "plain")
  got <- probe$probe_refusals()
  codes <- got[c("err_arg", "err_size", "err_alloc", "err_nonfinite")]
  want <- c(
    gemm_T_lda = "err_arg",
    gemm_T_ldb = "err_arg",
    gemv_T_lda = "err_arg",
    syrk_T_lda = "err_arg",
    gemm_k0_ld0 = "ok",
    solve_alloc = "err_alloc",
    syevd_n_32767 = "err_size",
    syevr_n_82595525 = "err_size"
  )
  for (nm in names(want)) {
    expected <- if (want[[nm]] == "ok") 0L else codes[[want[[nm]]]]
    say(
      "%-18s returned %6d, expected %6d (%s)%s",
      nm,
      got[[nm]],
      expected,
      want[[nm]],
      if (identical(as.integer(got[[nm]]), as.integer(expected))) {
        ""
      } else {
        "   <-- differs"
      }
    )
  }
  say("\nstatus of the failure paths the oracle accepts as 'any failure':")
  n <- 8
  tri <- qr.R(qr(rmat(n + 3, n)))
  d <- rep(1, n)
  d[[3]] <- NaN
  bad <- tri
  diag(bad) <- d
  out <- probe$probe_solve_tri("U", bad, rmat(n, 3))
  say(
    "solve_tri, NaN on the diagonal:   status %d (err_nonfinite is %d)",
    out$status,
    codes[["err_nonfinite"]]
  )
  d <- rep(1, n)
  d[[n]] <- 0
  bad <- tri
  diag(bad) <- d
  b <- rmat(n, 3)
  b[1, 1] <- Inf
  out <- probe$probe_solve_tri("U", bad, b)
  say(
    "solve_tri, singular A, Inf in B:  status %d (err_nonfinite is %d)",
    out$status,
    codes[["err_nonfinite"]]
  )
  x <- rmat(9, 4)
  x[3, 2] <- NaN
  say(
    "svd_values, NaN:                  status %d (err_nonfinite is %d)",
    probe$probe_svd(x)$status,
    codes[["err_nonfinite"]]
  )
}

if (group == "fuzz") {
  probe <- stage("dev/review_eigen_wrappers_probe.cpp", "plain")
  got <- probe$probe_fuzz(4L)
  say(
    "dimensions 0..4, every option: %d calls, %d returned a status other than 0",
    got[["calls"]],
    got[["nonzero_status"]]
  )
}

if (group == "order") {
  probe <- stage("dev/review_eigen_wrappers_probe.cpp", "plain")
  say("A 1 x k row times a k x 1 column. matmul takes the one-row test first")
  say("(dgemv 'T' on the column); the other order is dgemv 'N' on the row.")
  say("The oracle's only such products have k = 1 and k = 5.\n")
  for (k in c(1, 5, 8, 16, 50, 200, 1000, 5000)) {
    differs <- 0L
    worst <- 0
    reps <- 200L
    for (i in seq_len(reps)) {
      r <- probe$probe_dot_routes(rnorm(k), rnorm(k))
      if (!identical(r[["row_first"]], r[["column_first"]])) {
        differs <- differs + 1L
        worst <- max(worst, abs(r[["row_first"]] - r[["column_first"]]))
      }
    }
    say(
      "k %5d: the two dgemv routes differ in %3d of %d random products, worst abs %.3g",
      k,
      differs,
      reps,
      worst
    )
  }
}

if (group == "fallback_large") {
  pairs <- stage("dev/oracle_pairs.cpp", "pairs")
  say("RcppArmadillo %s", as.character(packageVersion("RcppArmadillo")))
  say("solve_tri's dgelsd fallback against arma::solve(); the oracle's")
  say(
    "fallback cases all have n = 8, below dgelsd's divide-and-conquer size.\n"
  )
  tri_case <- function(ul, n, nrhs, kind) {
    r <- qr.R(qr(rmat(n + 3, n)))
    if (kind == "zero diagonal") {
      diag(r)[[n %/% 2]] <- 0
    } else if (kind == "tiny diagonal") {
      diag(r)[[n]] <- 1e-20
    } else if (kind == "rank n/2") {
      r[(n %/% 2 + 1):n, ] <- 0
    }
    a <- if (ul == "U") r else t(r)
    junk <- rmat(n, n)
    other <- if (ul == "U") lower.tri(a) else upper.tri(a)
    a[other] <- junk[other]
    out <- pairs$pair_solve_tri(ul, a, rmat(n, nrhs))
    same_shape <- identical(dim(out$wrapper$value), dim(out$arma$value))
    say(
      "%s n %3d nrhs %2d %-14s wrapper status %d approx %-5s | Armadillo failed %-5s approx %-5s | max abs diff %.3g, bitwise %s",
      ul,
      n,
      nrhs,
      kind,
      out$wrapper$status,
      out$wrapper$approx,
      out$arma$failed,
      out$arma$approx,
      if (same_shape) max(abs(out$wrapper$value - out$arma$value)) else NA,
      same_shape && identical(out$wrapper$value, out$arma$value)
    )
  }
  for (ul in c("U", "L")) {
    for (n in c(8, 26, 30, 60, 130, 400)) {
      for (kind in c("zero diagonal", "tiny diagonal", "rank n/2")) {
        tri_case(ul, n, 3, kind)
      }
    }
    tri_case(ul, 60, 60, "zero diagonal")
    tri_case(ul, 130, 1, "zero diagonal")
  }
}

sym_cases <- function(n) {
  s <- gram(3 * n, n)
  mk <- function(i, j, v) {
    x <- s
    x[i, j] <- v
    x[j, i] <- v
    x
  }
  list(
    "NaN on the diagonal" = mk(2, 2, NaN),
    "NaN off the diagonal" = mk(1, 3, NaN),
    "Inf on the diagonal" = mk(2, 2, Inf),
    "-Inf off the diagonal" = mk(1, 3, -Inf),
    "all NaN" = matrix(NaN, n, n)
  )
}

if (group == "nonfinite_syevd") {
  probe <- stage("dev/review_eigen_wrappers_probe.cpp", "plain")
  say("si_la::syevd with NaN or Inf in the matrix; the wrapper has no check.")
  for (n in c(6, 60)) {
    cases <- sym_cases(n)
    for (nm in names(cases)) {
      for (ul in c("U", "L")) {
        say("calling: n %d, %s, %s", n, nm, ul)
        out <- probe$probe_syevd(ul, cases[[nm]])
        say(
          "  status %d; values: %s; vectors: %s",
          out$status,
          show_nonfinite(out$values),
          show_nonfinite(out$vectors)
        )
      }
    }
  }
}

if (group == "nonfinite_syevr") {
  probe <- stage("dev/review_eigen_wrappers_probe.cpp", "plain")
  pairs <- stage("dev/oracle_pairs.cpp", "pairs")
  say("si_la::syevr_range with NaN or Inf in the matrix; the wrapper has no")
  say("check. 'Armadillo' is arma::eig_sym() on the same matrix, the call the")
  say("oracle compares syevr_range with.")
  ## The same matrices in every child: built in one order from one seed.
  todo <- list()
  for (n in c(6, 60)) {
    cases <- sym_cases(n)
    for (nm in names(cases)) {
      for (rng in list(c(n - 1, n), c(1, n))) {
        todo[[length(todo) + 1L]] <- list(
          n = n,
          name = nm,
          s = cases[[nm]],
          il = rng[[1]],
          iu = rng[[2]]
        )
      }
    }
  }
  stopifnot(length(todo) == n_syevr_cases, length(case_id) == 1L)
  k <- todo[[case_id]]
  say("case %d: n %d, %s, pairs %d..%d", case_id, k$n, k$name, k$il, k$iu)
  out <- probe$probe_syevr("U", k$s, k$il, k$iu)
  say(
    "  wrapper status %d; values: %s; vectors: %s",
    out$status,
    show_nonfinite(out$values),
    show_nonfinite(out$vectors)
  )
  ## The oracle's pair runs the wrapper before Armadillo, so it is reached
  ## only when the wrapper returned above.
  ref <- pairs$pair_syevr("U", k$s, k$il, k$iu)
  say(
    "  in the oracle's pair: wrapper status %d, Armadillo failed %s",
    ref$wrapper$status,
    ref$arma$failed
  )
}
