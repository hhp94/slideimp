# Review of Phase 2 of dev/eigen_conversion.md at 4f2b669 (dev/review_p2.md):
# compare the cases dev/review_p2_cases.R wrote for two builds. Dev only:
#   Rscript dev/review_p2_compare.R <old dir> <new dir> [<group regex>]
# Per case it prints one line: whether both sides ended the same way; for two
# values, whether they are identical(), and if not, whether they are once
# every integer is read as a double (the n_iter ruling of step 6), the worst
# absolute and the worst relative difference over the numeric cells (the
# relative one over the cells the old side holds as nonzero, whose count is
# printed), and the attributes that differ; for two errors, whether class and
# message agree. It asserts nothing: the note reads the lines.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "usage: Rscript dev/review_p2_compare.R <old dir> <new dir> [<regex>]",
    call. = FALSE
  )
}
old_dir <- args[[1]]
new_dir <- args[[2]]
pat <- if (length(args) >= 3L) args[[3]] else "."

files <- sort(union(
  list.files(old_dir, "\\.rds$"),
  list.files(new_dir, "\\.rds$")
))
files <- grep(pat, files, value = TRUE)
if (!length(files)) {
  stop("no case matches", call. = FALSE)
}

# every integer and logical vector as a double, attributes kept, recursively
as_dbl <- function(x) {
  if (is.list(x)) {
    at <- attributes(x)
    x <- lapply(x, as_dbl)
    attributes(x) <- at
    return(x)
  }
  at <- attributes(x)
  if (is.integer(x)) {
    x <- as.double(x)
    attributes(x) <- at
  }
  if (!is.null(at)) {
    for (nm in names(at)) {
      if (is.integer(at[[nm]]) && !nm %in% c("dim")) {
        attr(x, nm) <- as.double(at[[nm]])
      }
    }
  }
  x
}

# the numeric leaves of a value, in order, as one double vector
num_leaves <- function(x) {
  if (is.list(x)) {
    return(unlist(lapply(x, num_leaves), use.names = FALSE))
  }
  if (is.numeric(x)) as.double(x) else double(0)
}

worst_abs <- 0
worst_rel <- 0
worst_abs_case <- ""
worst_rel_case <- ""
n_over_default <- 0L
worst <- function(a, b, nm) {
  if (length(a) != length(b)) {
    return("numeric cells differ in count")
  }
  if (!length(a)) {
    return("no numeric cells")
  }
  both_na <- is.na(a) & is.na(b)
  one_na <- xor(is.na(a), is.na(b))
  d <- abs(a - b)
  d[both_na] <- 0
  nz <- !is.na(a) & a != 0 & !is.na(b)
  m_abs <- max(d, na.rm = TRUE)
  m_rel <- if (any(nz)) max(d[nz] / abs(a[nz])) else 0
  if (m_abs > worst_abs) {
    worst_abs <<- m_abs
    worst_abs_case <<- nm
  }
  if (m_rel > worst_rel) {
    worst_rel <<- m_rel
    worst_rel_case <<- nm
  }
  # dev/ref_compare.R's defaults: a cell fails past both 1e-12 and 1e-10
  if (any(d[nz] > 1e-12 & d[nz] / abs(a[nz]) > 1e-10)) {
    n_over_default <<- n_over_default + 1L
  }
  sprintf(
    "worst abs %.3g, worst rel %.3g over %d nonzero cells, %d cells differ of %d%s",
    max(d, na.rm = TRUE),
    if (any(nz)) max(d[nz] / abs(a[nz])) else NA_real_,
    sum(nz),
    sum(d > 0, na.rm = TRUE) + sum(one_na),
    length(a),
    if (any(one_na)) sprintf(", %d NA on one side only", sum(one_na)) else ""
  )
}

attr_diff <- function(a, b) {
  nms <- union(names(attributes(a)), names(attributes(b)))
  out <- character()
  for (nm in nms) {
    x <- attr(a, nm, exact = TRUE)
    y <- attr(b, nm, exact = TRUE)
    if (!identical(x, y)) {
      show <- function(z) {
        if (is.null(z)) {
          "absent"
        } else if (length(z) == 1L && is.atomic(z)) {
          sprintf("%s %s", typeof(z), format(z, digits = 17))
        } else {
          sprintf("%s[%d]", typeof(z), length(z))
        }
      }
      out <- c(out, sprintf("%s: %s -> %s", nm, show(x), show(y)))
    }
  }
  out
}

tally <- c(
  same = 0L,
  same_as_dbl = 0L,
  differ = 0L,
  error_same = 0L,
  error_differ = 0L,
  ended_differently = 0L,
  missing = 0L
)
for (f in files) {
  nm <- sub("\\.rds$", "", f)
  fo <- file.path(old_dir, f)
  fn <- file.path(new_dir, f)
  if (!file.exists(fo) || !file.exists(fn)) {
    tally[["missing"]] <- tally[["missing"]] + 1L
    cat(sprintf(
      "%-40s MISSING on the %s side (no result: timed out or crashed)\n",
      nm,
      if (file.exists(fo)) "new" else "old"
    ))
    if (file.exists(fo) || file.exists(fn)) {
      r <- readRDS(if (file.exists(fo)) fo else fn)
      cat(sprintf(
        "    the other side: %s\n",
        if (inherits(r$value, "review_error")) {
          paste0("ERROR ", r$value$message)
        } else {
          "a value"
        }
      ))
    }
    next
  }
  o <- readRDS(fo)
  n <- readRDS(fn)
  eo <- inherits(o$value, "review_error")
  en <- inherits(n$value, "review_error")
  if (eo && en) {
    same <- identical(o$value$class, n$value$class) &&
      identical(o$value$message, n$value$message)
    key <- if (same) "error_same" else "error_differ"
    tally[[key]] <- tally[[key]] + 1L
    cat(sprintf(
      "%-40s both ERROR, %s\n",
      nm,
      if (same) "same class and message" else "DIFFERENT"
    ))
    if (!same) {
      cat(
        "    old: [",
        paste(o$value$class, collapse = "/"),
        "] ",
        gsub("\n", " | ", o$value$message),
        "\n",
        sep = ""
      )
      cat(
        "    new: [",
        paste(n$value$class, collapse = "/"),
        "] ",
        gsub("\n", " | ", n$value$message),
        "\n",
        sep = ""
      )
    } else {
      cat(
        "    ",
        gsub("\n", " | ", substr(o$value$message, 1L, 150L)),
        "\n",
        sep = ""
      )
    }
  } else if (eo || en) {
    tally[["ended_differently"]] <- tally[["ended_differently"]] + 1L
    cat(sprintf("%-40s ENDED DIFFERENTLY\n", nm))
    for (side in list(list("old", o), list("new", n))) {
      r <- side[[2]]
      cat(sprintf(
        "    %s: %s\n",
        side[[1]],
        if (inherits(r$value, "review_error")) {
          paste0(
            "ERROR [",
            paste(r$value$class, collapse = "/"),
            "] ",
            gsub("\n", " | ", r$value$message)
          )
        } else {
          "a value"
        }
      ))
    }
  } else {
    if (identical(o$value, n$value)) {
      tally[["same"]] <- tally[["same"]] + 1L
      cat(sprintf("%-40s identical\n", nm))
    } else if (identical(as_dbl(o$value), as_dbl(n$value))) {
      tally[["same_as_dbl"]] <- tally[["same_as_dbl"]] + 1L
      cat(sprintf("%-40s identical once integers are read as doubles\n", nm))
      ad <- attr_diff(o$value, n$value)
      if (length(ad)) cat(paste0("    ", ad, "\n"), sep = "")
    } else {
      tally[["differ"]] <- tally[["differ"]] + 1L
      cat(sprintf(
        "%-40s DIFFER: %s\n",
        nm,
        worst(num_leaves(unclass(o$value)), num_leaves(unclass(n$value)), nm)
      ))
      ad <- attr_diff(o$value, n$value)
      if (length(ad)) cat(paste0("    ", ad, "\n"), sep = "")
    }
  }
  if (!identical(o$warnings, n$warnings)) {
    cat(
      "    warnings differ\n      old: ",
      paste(o$warnings, collapse = " || "),
      "\n      new: ",
      paste(n$warnings, collapse = " || "),
      "\n",
      sep = ""
    )
  }
  if (!identical(o$messages, n$messages)) {
    cat("    messages differ\n")
  }
}
if (tally[["differ"]] > 0L) {
  cat(sprintf(
    paste0(
      "over the %d that differ, the result's own cells only (attributes ",
      "are listed per case): worst abs %.3g (%s), worst rel %.3g (%s); ",
      "%d hold a cell past both 1e-12 absolute and 1e-10 relative\n"
    ),
    tally[["differ"]],
    worst_abs,
    worst_abs_case,
    worst_rel,
    worst_rel_case,
    n_over_default
  ))
}
cat(sprintf(
  "%d cases: %s\n",
  length(files),
  paste(sprintf("%s %d", names(tally), tally), collapse = ", ")
))
