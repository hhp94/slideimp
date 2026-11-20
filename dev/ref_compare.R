# Reference compare: check one capture from dev/ref_capture.R against
# another. The usual pair is the reference build and a changed build, each
# captured under its own label:
#   Rscript dev/ref_capture.R e5b3cf7       # once, at the reference commit
#   Rscript dev/ref_capture.R wip           # on the changed build
#   Rscript dev/ref_compare.R e5b3cf7 wip
# Arguments are labels under dev/ref/, or paths to capture folders. Options:
#   --tol-abs=<x>  --tol-rel=<x>   tolerances for double values (defaults
#                                  below, set from the noise floor)
#   --only=<regex>                 compare only cases whose name matches
#   --verbose                      print every compared element, not only
#                                  the worst per case
#
# What it checks, walking each stored result and its attributes:
# - Type, length, attribute names, and every non-double value (integer,
#   logical, character, language, the NA placement, index maps): identical().
#   Dims, dimnames and class are attributes, so they are checked this way.
#   Attribute names are compared as a set, not in order: identical() ignores
#   their order too, and nothing in the package reads an attribute by
#   position.
# - Doubles: NA, NaN and Inf in the same cells, identical(); over the finite
#   cells it reports max(abs(a - b)) and max(abs(a - b) / m), where m is
#   max(abs(a), abs(b)) of that cell, the relative figure taken over the
#   cells where m != 0 (it says how many). Dividing by the larger of the two
#   makes the figure the same whichever capture is a: dividing by one side
#   only lets a cell that side holds as 0 pass while the other holds a value
#   under the abs tolerance. It asserts both maxima against the tolerances,
#   and reports whether the two are bitwise identical without asserting it,
#   comparing the bytes, so a 0 that became -0 shows as not bitwise.
# - A case that errored in either capture: both must have errored with the
#   same class and message. The warnings each case raised: same classes and
#   messages, in order.
# - Functions and environments are skipped and counted.
# Cases present in only one capture fail. Exit status is 1 on any failure.
#
# The noise floor (dev/eigen_conversion.md at 4f2b669, Phase 0): capture the
# reference build twice, in two processes, and compare the two. On the reference
# build that run is bitwise, so the defaults below are not tighter than it. They
# are loose enough for a change of summation order in the last bits (a few ulp)
# and tight enough to fail on anything larger. A step that needs looser passes
# --tol-abs and --tol-rel and says why in its commit. dev/ref_compare_selftest.R
# proves this script still fails where it should; rerun it after changing the
# walk or the defaults.

tol_abs <- 1e-12
tol_rel <- 1e-10

args <- commandArgs(trailingOnly = TRUE)
opts <- grep("^--", args, value = TRUE)
pos <- setdiff(args, opts)
if (length(pos) != 2L) {
  stop(
    "usage: Rscript dev/ref_compare.R <ref> <new> [--tol-abs=x] ",
    "[--tol-rel=x] [--only=regex] [--verbose]",
    call. = FALSE
  )
}
opt_value <- function(name) {
  hit <- grep(paste0("^--", name, "="), opts, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else NULL
}
if (!is.null(v <- opt_value("tol-abs"))) {
  tol_abs <- as.numeric(v)
}
if (!is.null(v <- opt_value("tol-rel"))) {
  tol_rel <- as.numeric(v)
}
only <- opt_value("only")
verbose <- "--verbose" %in% opts

resolve <- function(x) {
  if (dir.exists(x)) x else file.path("dev", "ref", x)
}
dir_a <- resolve(pos[[1]])
dir_b <- resolve(pos[[2]])
for (d in c(dir_a, dir_b)) {
  if (!file.exists(file.path(d, "_meta.rds"))) {
    stop("no capture at ", d, call. = FALSE)
  }
}

show_meta <- function(d, role) {
  m <- readRDS(file.path(d, "_meta.rds"))
  ## extSoftVersion() leaves BLAS empty when Rblas.dll is swapped for
  ## another build, so the LAPACK path is shown alongside it.
  cat(sprintf(
    "%s: %s, source %s%s, %s, BLAS '%s', LAPACK %s, %s\n",
    role,
    m$label,
    m$source_commit,
    if (m$source_dirty) " (DIRTY)" else "",
    m$r_version,
    m$blas,
    basename(m$la_library),
    m$time
  ))
  ## Captures older than these fields leave them NULL.
  if (!is.null(m$rcppeigen)) {
    cat(sprintf(
      "  RcppEigen %s from %s, RcppArmadillo %s\n",
      m$rcppeigen,
      m$rcppeigen_lib,
      m$rcpparmadillo
    ))
  }
}
show_meta(dir_a, "ref")
show_meta(dir_b, "new")
cat(sprintf("tolerance: abs %g, rel %g\n\n", tol_abs, tol_rel))

## ---- the walk -------------------------------------------------------------

## Each compared element appends one row here.
rows <- list()
add <- function(
  path,
  kind,
  ok,
  note = "",
  abs = NA,
  rel = NA,
  n_rel = NA,
  bitwise = NA
) {
  rows[[length(rows) + 1L]] <<- data.frame(
    path = path,
    kind = kind,
    ok = ok,
    abs = abs,
    rel = rel,
    n_rel = n_rel,
    bitwise = bitwise,
    note = note
  )
}

bare <- function(x) {
  attributes(x) <- NULL
  x
}

walk <- function(a, b, path) {
  if (!identical(typeof(a), typeof(b))) {
    add(path, "type", FALSE, sprintf("%s vs %s", typeof(a), typeof(b)))
    return(invisible())
  }
  type <- typeof(a)
  if (
    type %in% c("closure", "builtin", "special", "environment", "externalptr")
  ) {
    add(path, "skip", TRUE, type)
    return(invisible())
  }

  ## Attributes first: same names, then each one walked.
  at_a <- attributes(a)
  at_b <- attributes(b)
  if (!identical(sort(names(at_a)), sort(names(at_b)))) {
    add(
      path,
      "attr names",
      FALSE,
      sprintf(
        "{%s} vs {%s}",
        paste(sort(names(at_a)), collapse = ", "),
        paste(sort(names(at_b)), collapse = ", ")
      )
    )
  } else {
    for (nm in names(at_a)) {
      walk(at_a[[nm]], at_b[[nm]], sprintf("%s@%s", path, nm))
    }
  }

  if (length(a) != length(b)) {
    add(path, "length", FALSE, sprintf("%d vs %d", length(a), length(b)))
    return(invisible())
  }

  if (type == "list") {
    for (i in seq_along(a)) {
      walk(a[[i]], b[[i]], sprintf("%s[[%d]]", path, i))
    }
  } else if (type == "double") {
    va <- bare(a)
    vb <- bare(b)
    same_special <- identical(is.na(va), is.na(vb)) &&
      identical(is.nan(va), is.nan(vb)) &&
      identical(va[is.infinite(va)], vb[is.infinite(vb)]) &&
      identical(is.infinite(va), is.infinite(vb))
    if (!same_special) {
      add(path, "double", FALSE, "NA/NaN/Inf cells differ")
      return(invisible())
    }
    fin <- is.finite(va)
    d <- abs(va[fin] - vb[fin])
    m <- pmax(abs(va[fin]), abs(vb[fin]))
    nz <- m != 0
    max_abs <- if (length(d)) max(d) else 0
    max_rel <- if (any(nz)) max(d[nz] / m[nz]) else 0
    ok <- max_abs <= tol_abs && max_rel <= tol_rel
    add(
      path,
      "double",
      ok,
      note = sprintf("%d finite", sum(fin)),
      abs = max_abs,
      rel = max_rel,
      n_rel = sum(nz),
      ## identical() treats 0 and -0 as equal; the bytes do not.
      bitwise = identical(writeBin(va, raw()), writeBin(vb, raw()))
    )
  } else {
    add(path, type, identical(bare(a), bare(b)))
  }
  invisible()
}

compare_case <- function(name, ra, rb) {
  ea <- inherits(ra$value, "ref_error")
  eb <- inherits(rb$value, "ref_error")
  if (ea || eb) {
    ok <- ea && eb && identical(unclass(ra$value), unclass(rb$value))
    add(
      "value",
      "error",
      ok,
      sprintf(
        "%s | %s",
        if (ea) ra$value$message else "(no error)",
        if (eb) rb$value$message else "(no error)"
      )
    )
  } else {
    walk(ra$value, rb$value, "value")
  }
  add(
    "warnings",
    "warnings",
    identical(ra$warnings, rb$warnings),
    sprintf("%d vs %d", length(ra$warnings), length(rb$warnings))
  )
}

## ---- run ------------------------------------------------------------------

case_names <- function(d) {
  f <- list.files(d, pattern = "\\.rds$")
  sub("\\.rds$", "", f[!startsWith(f, "_")])
}
na <- case_names(dir_a)
nb <- case_names(dir_b)
all_names <- union(na, nb)
if (!is.null(only)) {
  all_names <- grep(only, all_names, value = TRUE)
}
if (!length(all_names)) {
  stop("no cases to compare", call. = FALSE)
}

fmt <- function(x) {
  ifelse(is.na(x), "", formatC(x, format = "e", digits = 2))
}

n_fail <- 0L
summary_rows <- list()
for (name in all_names) {
  rows <- list()
  if (!(name %in% na) || !(name %in% nb)) {
    add(
      "",
      "case",
      FALSE,
      if (name %in% na) {
        "missing in new"
      } else {
        "missing in ref"
      }
    )
  } else {
    compare_case(
      name,
      readRDS(file.path(dir_a, paste0(name, ".rds"))),
      readRDS(file.path(dir_b, paste0(name, ".rds")))
    )
  }
  r <- do.call(rbind, rows)
  dbl <- r[r$kind == "double" & !is.na(r$abs), , drop = FALSE]
  case_ok <- all(r$ok)
  n_fail <- n_fail + !case_ok
  summary_rows[[name]] <- data.frame(
    case = name,
    status = if (case_ok) "ok" else "FAIL",
    elements = nrow(r),
    doubles = nrow(dbl),
    max_abs = if (nrow(dbl)) max(dbl$abs) else NA,
    max_rel = if (nrow(dbl)) max(dbl$rel) else NA,
    bitwise = if (nrow(dbl)) all(dbl$bitwise) else NA,
    skipped = sum(r$kind == "skip")
  )
  bad <- r[!r$ok, , drop = FALSE]
  show <- if (verbose) r else bad
  if (nrow(show)) {
    cat(sprintf("%s:\n", name))
    for (i in seq_len(nrow(show))) {
      cat(sprintf(
        "  %-4s %-40s %-10s abs %-9s rel %-9s (over %s) %s %s\n",
        if (show$ok[i]) "ok" else "FAIL",
        show$path[i],
        show$kind[i],
        fmt(show$abs[i]),
        fmt(show$rel[i]),
        ifelse(is.na(show$n_rel[i]), "-", show$n_rel[i]),
        ifelse(
          is.na(show$bitwise[i]),
          "",
          ifelse(show$bitwise[i], "bitwise", "not bitwise")
        ),
        show$note[i]
      ))
    }
  }
}

s <- do.call(rbind, summary_rows)
cat(paste0(
  "\nper case (doubles: worst element; rel over cells where ",
  "max(|ref|, |new|) != 0):\n"
))
cat(
  sprintf(
    "  %-28s %-4s %5s %5s  abs %-9s rel %-9s %s%s\n",
    s$case,
    s$status,
    s$elements,
    s$doubles,
    fmt(s$max_abs),
    fmt(s$max_rel),
    ifelse(is.na(s$bitwise), "", ifelse(s$bitwise, "bitwise", "NOT bitwise")),
    ifelse(s$skipped > 0, sprintf("  (%d skipped)", s$skipped), "")
  ),
  sep = ""
)
cat(sprintf(
  "\n%d case(s), %d failed; worst abs %s, worst rel %s, %d not bitwise\n",
  nrow(s),
  n_fail,
  fmt(suppressWarnings(max(s$max_abs, na.rm = TRUE))),
  fmt(suppressWarnings(max(s$max_rel, na.rm = TRUE))),
  sum(s$bitwise %in% FALSE)
))
quit(status = if (n_fail > 0L) 1L else 0L)
