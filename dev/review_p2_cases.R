# Review of Phase 2 of dev/eigen_conversion.md at 4f2b669 (dev/review_p2.md):
# cases run through the public pca_imp() on one installed build of slideimp,
# each result or error written to its own .rds, so two builds can be compared by
# dev/review_p2_compare.R. Dev only, run by hand from the package root:
#   Rscript dev/review_p2_cases.R <lib> <out dir> <group> [<case>]
# <lib> holds the build under test (dev/eigen_time_install.R made it) and
# slideimp must load from there, or the script stops. <group> is one of
# tiny, tinysweep, general, big, validation, overflow, dup; with <case>,
# only that case runs, which is how the script gives each overflow and dup
# case a process of its own with a time limit.
#   Rscript dev/review_p2_cases.R <lib> <out dir> --list <group>
# prints the case names of a group and runs nothing.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "usage: Rscript dev/review_p2_cases.R <lib> <out dir> <group> [<case>]",
    call. = FALSE
  )
}
lib <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- args[[2]]
list_only <- identical(args[[3]], "--list")
group <- if (list_only) args[[4]] else args[[3]]
only_case <- if (!list_only && length(args) >= 4L) args[[4]] else NULL

suppressPackageStartupMessages(library(slideimp, lib.loc = lib))
from <- normalizePath(
  dirname(getNamespaceInfo("slideimp", "path")),
  winslash = "/"
)
if (!identical(tolower(from), tolower(lib))) {
  stop("slideimp loaded from ", from, ", not ", lib, call. = FALSE)
}

## ---- inputs ---------------------------------------------------------------

# n x p of rank `r` plus noise, with `n_na` holes placed so that every row
# keeps at least one observed cell and every column at least two.
mk <- function(n, p, seed, n_na, r = 2L, noise = 0.3) {
  set.seed(seed)
  r <- min(r, n, p)
  x <- matrix(rnorm(n * r), n, r) %*%
    matrix(rnorm(r * p), r, p) +
    noise * matrix(rnorm(n * p), n, p)
  holes <- 0L
  tries <- 0L
  while (holes < n_na && tries < 100000L) {
    tries <- tries + 1L
    i <- sample.int(n, 1L)
    j <- sample.int(p, 1L)
    if (is.na(x[i, j])) {
      next
    }
    if (sum(!is.na(x[i, ])) <= 1L || sum(!is.na(x[, j])) <= 2L) {
      next
    }
    x[i, j] <- NA
    holes <- holes + 1L
  }
  stopifnot(holes == n_na)
  x
}

lob <- function(warmup = 0L, maxiter = 20L, tol = 1e-9) {
  lobpcg_control(warmup_iters = warmup, tol = tol, maxiter = maxiter)
}

cases <- list()
add <- function(grp, name, f) {
  cases[[length(cases) + 1L]] <<- list(group = grp, name = name, f = f)
}

## tiny: a Gram of order 2 to 4 under LOBPCG, where Armadillo ran its own
## loops for a product and si_la::matmul calls BLAS
## (dev/pca_imp_structure.md, "Two things Armadillo does that si_la::matmul
## does not"). k is the eigenpair count, n the order of the Gram.
tiny_fit <- function(x, ncp, method, solver) {
  pca_imp(
    x,
    ncp = ncp,
    method = method,
    solver = solver,
    lobpcg_control = if (solver == "lobpcg") lob() else NULL,
    threshold = 1e-12,
    maxiter = 300
  )
}
for (solver in c("lobpcg", "exact")) {
  local({
    s <- solver
    # k == n: every eigenpair of the Gram
    add("tiny", paste0("t2_reg1_", s), function() {
      tiny_fit(mk(12, 2, 11, 3, r = 1L), 1, "regularized", s)
    })
    add("tiny", paste0("t3_reg2_", s), function() {
      tiny_fit(mk(12, 3, 12, 4), 2, "regularized", s)
    })
    add("tiny", paste0("t4_reg3_", s), function() {
      tiny_fit(mk(15, 4, 13, 6, r = 3L), 3, "regularized", s)
    })
    # k == 1: the Gram times one column
    add("tiny", paste0("t2_em1_", s), function() {
      tiny_fit(mk(12, 2, 14, 3, r = 1L), 1, "EM", s)
    })
    add("tiny", paste0("t3_em1_", s), function() {
      tiny_fit(mk(12, 3, 15, 4, r = 1L), 1, "EM", s)
    })
    add("tiny", paste0("t4_em1_", s), function() {
      tiny_fit(mk(15, 4, 16, 6, r = 1L), 1, "EM", s)
    })
    add("tiny", paste0("w3_em1_", s), function() {
      tiny_fit(mk(3, 9, 17, 2, r = 1L), 1, "EM", s)
    })
    add("tiny", paste0("w4_em1_", s), function() {
      tiny_fit(mk(4, 9, 18, 3, r = 1L), 1, "EM", s)
    })
    # just outside the rule: order 5, and order 4 with k of 2 and 3
    add("tiny", paste0("t5_reg4_", s), function() {
      tiny_fit(mk(16, 5, 19, 6, r = 4L), 4, "regularized", s)
    })
    add("tiny", paste0("t5_em1_", s), function() {
      tiny_fit(mk(16, 5, 20, 6, r = 1L), 1, "EM", s)
    })
    add("tiny", paste0("t4_reg2_", s), function() {
      tiny_fit(mk(15, 4, 21, 6), 2, "regularized", s)
    })
    add("tiny", paste0("t4_reg1_", s), function() {
      tiny_fit(mk(15, 4, 22, 6, r = 1L), 1, "regularized", s)
    })
  })
}

## tinysweep: the same shapes as a user would call them, solver = "lobpcg"
## with the default control and the default threshold, 40 inputs a shape,
## to count how often such a fit moves and by how much.
sweep_shapes <- list(
  t2_em1 = list(n = 12, p = 2, na = 3, r = 1L, ncp = 1, method = "EM"),
  t3_em1 = list(n = 12, p = 3, na = 4, r = 1L, ncp = 1, method = "EM"),
  t4_em1 = list(n = 15, p = 4, na = 6, r = 1L, ncp = 1, method = "EM"),
  t3_reg2 = list(
    n = 12,
    p = 3,
    na = 4,
    r = 2L,
    ncp = 2,
    method = "regularized"
  ),
  t4_reg3 = list(
    n = 15,
    p = 4,
    na = 6,
    r = 3L,
    ncp = 3,
    method = "regularized"
  ),
  w3_em1 = list(n = 3, p = 9, na = 2, r = 1L, ncp = 1, method = "EM"),
  w4_em1 = list(n = 4, p = 9, na = 3, r = 1L, ncp = 1, method = "EM"),
  # outside the rule
  t5_em1 = list(n = 16, p = 5, na = 6, r = 1L, ncp = 1, method = "EM"),
  t5_reg4 = list(n = 16, p = 5, na = 6, r = 4L, ncp = 4, method = "regularized")
)
for (shape in names(sweep_shapes)) {
  for (seed in seq_len(40L)) {
    local({
      sh <- sweep_shapes[[shape]]
      sd <- 1000L + seed
      add("tinysweep", sprintf("%s_s%02d", shape, seed), function() {
        pca_imp(
          mk(sh$n, sh$p, sd, sh$na, r = sh$r),
          ncp = sh$ncp,
          method = sh$method,
          solver = "lobpcg"
        )
      })
    })
  }
}

## general: every public switch once, at ordinary sizes
x_tall <- function() mk(60, 40, 31, 240)
x_wide <- function() mk(40, 60, 32, 240)
add("general", "tall_lobpcg5", function() {
  pca_imp(x_tall(), ncp = 5, solver = "lobpcg", lobpcg_control = lob(2L))
})
add("general", "tall_exact5", function() {
  pca_imp(x_tall(), ncp = 5, solver = "exact")
})
add("general", "wide_lobpcg3", function() {
  pca_imp(x_wide(), ncp = 3, solver = "lobpcg", lobpcg_control = lob(2L))
})
add("general", "wide_exact3_em", function() {
  pca_imp(x_wide(), ncp = 3, solver = "exact", method = "EM")
})
add("general", "tall_default", function() pca_imp(x_tall()))
add("general", "tall_nbinit3_seed7", function() {
  pca_imp(x_tall(), ncp = 2, nb.init = 3, seed = 7, solver = "exact")
})
add("general", "tall_nbinit3_seed7_lobpcg", function() {
  pca_imp(
    x_tall(),
    ncp = 2,
    nb.init = 3,
    seed = 7,
    solver = "lobpcg",
    lobpcg_control = lob(2L)
  )
})
add("general", "tall_roww_runif", function() {
  set.seed(5)
  w <- runif(60, 0.2, 3)
  pca_imp(x_tall(), ncp = 2, row.w = w, solver = "exact")
})
add("general", "wide_roww_nmiss", function() {
  pca_imp(x_wide(), ncp = 2, row.w = "n_miss", solver = "exact")
})
add("general", "tall_noscale", function() {
  pca_imp(x_tall(), ncp = 2, scale = FALSE, solver = "exact")
})
add("general", "tall_ridge_half", function() {
  pca_imp(x_tall(), ncp = 2, coeff.ridge = 0.5, solver = "exact")
})
add("general", "tall_clamp", function() {
  pca_imp(x_tall(), ncp = 2, clamp = c(-0.5, 0.5), solver = "exact")
})
add("general", "tall_subset", function() {
  pca_imp(x_tall(), ncp = 2, subset = c(3L, 7L, 21L), solver = "exact")
})
add("general", "tall_mipca", function() {
  x <- x_tall()
  dimnames(x) <- list(paste0("s", seq_len(60)), paste0("f", seq_len(40)))
  unclass(pca_imp(x, ncp = 2, output = "mipca", solver = "exact"))
})
add("general", "tall_mipca_types", function() {
  # the type of every `n_iter` and `iter` the fit holds, by its path
  x <- x_tall()
  dimnames(x) <- list(paste0("s", seq_len(60)), paste0("f", seq_len(40)))
  fit <- unclass(pca_imp(x, ncp = 2, output = "mipca", solver = "exact"))
  out <- character()
  walk <- function(v, path) {
    if (is.list(v) && !is.null(names(v))) {
      for (nm in names(v)) {
        if (nzchar(nm)) walk(v[[nm]], c(path, nm))
      }
    } else if (is.list(v)) {
      for (i in seq_along(v)) {
        walk(v[[i]], c(path, sprintf("[[%d]]", i)))
      }
    }
    leaf <- path[length(path)]
    if (length(path) && leaf %in% c("n_iter", "iter")) {
      out[[paste(path, collapse = "$")]] <<- typeof(v)
    }
  }
  walk(fit, character())
  out
})
add("general", "tall_no_holes", function() {
  x <- x_tall()
  x[is.na(x)] <- 0
  pca_imp(x, ncp = 2)
})
add("general", "tall_maxiter3", function() {
  pca_imp(x_tall(), ncp = 2, maxiter = 5, miniter = 5, threshold = 0)
})
add("general", "tall_progress", function() {
  txt <- utils::capture.output(
    res <- pca_imp(x_tall(), ncp = 2, solver = "exact", .progress = TRUE)
  )
  list(res = res, printed = txt)
})
add("general", "auto_probe_260x300", function() {
  x <- mk(260, 300, 33, 3000, r = 3L)
  res <- pca_imp(x, ncp = 3)
  # the auto verdict rests on wall time, so it may differ between two runs
  # of one build; the compare reports it and does not count it
  res
})
add("general", "integer_input", function() {
  set.seed(34)
  x <- matrix(sample.int(50L, 40L * 12L, replace = TRUE), 40L, 12L)
  x[sample(length(x), 30L)] <- NA
  pca_imp(x, ncp = 2, solver = "exact")
})
add("general", "dup_cols_exact", function() {
  x <- mk(30, 8, 35, 20)
  x <- cbind(x, x)
  pca_imp(x, ncp = 9, solver = "exact")
})
add("general", "dup_rows_exact", function() {
  x <- mk(8, 30, 36, 20)
  x <- rbind(x, x)
  pca_imp(x, ncp = 9, solver = "exact")
})

## big: the two deliberate moves of the phase, for scale
add("big", "exact_ncp33", function() {
  pca_imp(mk(80, 70, 41, 500, r = 5L), ncp = 33, solver = "exact")
})
add("big", "exact_ncp32", function() {
  pca_imp(mk(80, 70, 41, 500, r = 5L), ncp = 32, solver = "exact")
})
add("big", "lobpcg_ncp16", function() {
  pca_imp(
    mk(120, 100, 42, 900, r = 5L),
    ncp = 16,
    solver = "lobpcg",
    lobpcg_control = lob(2L)
  )
})
add("big", "lobpcg_ncp15", function() {
  pca_imp(
    mk(120, 100, 42, 900, r = 5L),
    ncp = 15,
    solver = "lobpcg",
    lobpcg_control = lob(2L)
  )
})

## validation: what a user can pass that the kernel's new int arguments
## could treat differently
v <- function(...) pca_imp(x_tall(), ...)
add("validation", "maxiter_3e9", function() v(maxiter = 3e9))
add("validation", "maxiter_intmax", function() {
  v(maxiter = .Machine$integer.max, solver = "exact")
})
add("validation", "maxiter_intmax_dbl", function() {
  v(maxiter = 2147483647, solver = "exact")
})
add("validation", "maxiter_0", function() v(maxiter = 0))
add("validation", "maxiter_na", function() v(maxiter = NA_integer_))
add("validation", "miniter_3e9", function() v(miniter = 3e9, maxiter = 3e9))
add("validation", "miniter_gt_maxiter", function() v(miniter = 9, maxiter = 3))
add("validation", "ncp_na", function() v(ncp = NA_integer_))
add("validation", "ncp_half", function() v(ncp = 2.5))
add("validation", "ncp_dbl", function() v(ncp = 2.0, solver = "exact"))
add("validation", "ncp_cap", function() v(ncp = 39, solver = "exact"))
add("validation", "ncp_over_cap", function() v(ncp = 40))
add("validation", "lobpcg_maxiter_3e9", function() {
  v(solver = "lobpcg", lobpcg_control = list(maxiter = 3e9))
})
add("validation", "lobpcg_maxiter_intmax", function() {
  v(
    solver = "lobpcg",
    lobpcg_control = list(maxiter = .Machine$integer.max, warmup_iters = 2L)
  )
})
add("validation", "lobpcg_warmup_intmax", function() {
  v(
    solver = "lobpcg",
    lobpcg_control = list(warmup_iters = .Machine$integer.max)
  )
})
add("validation", "roww_short", function() v(row.w = rep(1, 59)))
add("validation", "roww_zero", function() v(row.w = rep(0, 60)))
add("validation", "roww_floor", function() {
  v(row.w = rep(1e-10, 60), solver = "exact")
})
add("validation", "roww_one_heavy", function() {
  v(row.w = c(1e10, rep(1e-10, 59)), solver = "exact")
})
add("validation", "obj_inf", function() {
  x <- x_tall()
  x[1, 1] <- Inf
  pca_imp(x)
})
add("validation", "nbinit_0", function() v(nb.init = 0))
add("validation", "nbinit_big_seed", function() {
  v(nb.init = 3, seed = .Machine$integer.max)
})

## overflow: input whose arithmetic leaves the double range inside the fit,
## so a Gram can hold Inf or NaN though the input holds neither. Each runs in
## its own process under a time limit: the chain's own dsyevr call was
## measured not to return on some such matrices
## (dev/review_eigen_wrappers.md, finding 2).
ovf <- function(scale_by, scale, solver, n = 30, p = 10) {
  x <- mk(n, p, 51, 30) * scale_by
  pca_imp(
    x,
    ncp = 2,
    scale = scale,
    solver = solver,
    lobpcg_control = if (solver == "lobpcg") lob() else NULL
  )
}
for (solver in c("exact", "lobpcg")) {
  local({
    s <- solver
    add("overflow", paste0("e200_noscale_", s), function() ovf(1e200, FALSE, s))
    add("overflow", paste0("e200_scale_", s), function() ovf(1e200, TRUE, s))
    add("overflow", paste0("e160_noscale_", s), function() ovf(1e160, FALSE, s))
    add("overflow", paste0("e150_noscale_", s), function() ovf(1e150, FALSE, s))
    add("overflow", paste0("em200_noscale_", s), function() {
      ovf(1e-200, FALSE, s)
    })
    # ncp 9 of 10 columns, regularized: every eigenpair, the range dsyevr
    # was measured not to return on
    add("overflow", paste0("e200_noscale_full_", s), function() {
      x <- mk(30, 10, 52, 30) * 1e200
      pca_imp(
        x,
        ncp = 9,
        scale = FALSE,
        solver = s,
        lobpcg_control = if (s == "lobpcg") lob() else NULL
      )
    })
  })
}

## dup: the rank-deficient inputs of dev/eigen_p2s4_edges.R, built as that
## script builds them and called as it calls them (set.seed(1) first), each
## in its own process: the instrumented clone of 59920f4 died on one of them.
## xd is 30 samples of which 6 are distinct, holes included; xdc is 40 x 12
## with columns 4 to 12 copies of columns 1 to 3.
dup_inputs <- function() {
  set.seed(21)
  base6 <- sim_mat(n = 6, p = 80, perc_total_na = 0.1, perc_col_na = 0.5)$input
  xd <- base6[rep(seq_len(6), times = 5), ]
  rownames(xd) <- sprintf("s%02d", seq_len(nrow(xd)))
  set.seed(22)
  base3 <- matrix(rnorm(40 * 3), 40, 3)
  xdc <- base3[, rep(1:3, times = 4)]
  xdc[cbind(c(2, 9, 17, 30), c(1, 5, 8, 12))] <- NA
  dimnames(xdc) <- list(sprintf("s%02d", 1:40), sprintf("c%02d", 1:12))
  list(xd = xd, xdc = xdc)
}
dup_fit <- function(which, ncp, solver, method = "regularized") {
  x <- dup_inputs()[[which]]
  set.seed(1)
  if (solver == "lobpcg") {
    pca_imp(
      x,
      ncp = ncp,
      method = method,
      solver = "lobpcg",
      lobpcg_control = lobpcg_control(warmup_iters = 1),
      miniter = 15
    )
  } else {
    pca_imp(x, ncp = ncp, method = method, solver = solver)
  }
}
for (solver in c("lobpcg", "exact", "auto")) {
  local({
    s <- solver
    for (ncp in c(5L, 7L, 12L)) {
      local({
        k <- ncp
        add("dup", sprintf("xd_ncp%02d_%s", k, s), function() {
          dup_fit("xd", k, s)
        })
      })
    }
    add("dup", paste0("xd_ncp07_em_", s), function() {
      dup_fit("xd", 7L, s, "EM")
    })
    add("dup", paste0("xdc_ncp05_", s), function() dup_fit("xdc", 5L, s))
    add("dup", paste0("xdc_ncp08_", s), function() dup_fit("xdc", 8L, s))
  })
}

## ---- run ------------------------------------------------------------------

in_group <- Filter(function(cs) identical(cs$group, group), cases)
if (!length(in_group)) {
  stop("no case in group ", group, call. = FALSE)
}
if (list_only) {
  cat(vapply(in_group, function(cs) cs$name, ""), sep = "\n")
  quit(status = 0L)
}
if (!is.null(only_case)) {
  in_group <- Filter(function(cs) identical(cs$name, only_case), in_group)
  if (length(in_group) != 1L) {
    stop("no case ", only_case, " in group ", group, call. = FALSE)
  }
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# the overflow and dup groups: one child process per case, with a time
# limit, since a LAPACK call that does not return cannot be interrupted from
# R and a process that dies takes the cases after it along.
if (group %in% c("overflow", "dup") && is.null(only_case)) {
  script <- sub(
    "^--file=",
    "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]]
  )
  limit <- 120
  n_stuck <- 0L
  for (cs in in_group) {
    st <- system2(
      file.path(R.home("bin"), "Rscript"),
      c(shQuote(script), shQuote(lib), shQuote(out_dir), group, cs$name),
      timeout = limit
    )
    if (!identical(st, 0L)) {
      n_stuck <- n_stuck + 1L
      cat(sprintf(
        "%-32s NO RESULT: child status %s%s\n",
        cs$name,
        st,
        if (identical(st, 124L)) sprintf(", no return in %d s", limit) else ""
      ))
    }
  }
  cat(sprintf(
    "%d %s cases, each in its own process with a %d s limit; %d gave no result\n",
    length(in_group),
    group,
    limit,
    n_stuck
  ))
  quit(status = 0L)
}

run_case <- function(f) {
  warns <- character()
  msgs <- character()
  t0 <- Sys.time()
  val <- withCallingHandlers(
    tryCatch(
      f(),
      error = function(e) {
        structure(
          list(class = class(e), message = conditionMessage(e)),
          class = "review_error"
        )
      }
    ),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(
    value = val,
    warnings = warns,
    messages = msgs,
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
  )
}

n_err <- 0L
for (cs in in_group) {
  res <- run_case(cs$f)
  saveRDS(res, file.path(out_dir, paste0(group, "__", cs$name, ".rds")))
  is_err <- inherits(res$value, "review_error")
  n_err <- n_err + is_err
  cat(sprintf(
    "%-32s %s (%.2fs)%s\n",
    cs$name,
    if (is_err) "ERROR" else "value",
    res$seconds,
    if (is_err) paste0(": ", substr(res$value$message, 1L, 90L)) else ""
  ))
}
cat(sprintf(
  "%d cases, %d ended in an error; slideimp %s from %s\n",
  length(in_group),
  n_err,
  as.character(packageVersion("slideimp")),
  lib
))
