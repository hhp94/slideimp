# Step 4 of Phase 2 (dev/eigen_conversion.md): read the event logs that the
# instrumented clones of dev/probe_p2s4_run.sh wrote, one of the code before
# the step and one of the code after it, over the same run.
#   Rscript dev/probe_p2s4_summary.R <old log> <new log>
# It prints, for each: how often each call was made and failed, how every
# lobpcg_solve() call ended, and the orders dsyevd was asked for; then the
# cases in which the two builds did not end their lobpcg_solve() calls alike.
#
# With two capture folders after the logs, it also prints every case that
# ended in an error in either, with the class and message each build raised:
#   Rscript dev/probe_p2s4_summary.R <old log> <new log> <old dir> <new dir>
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 4L)) {
  stop(
    "usage: Rscript dev/probe_p2s4_summary.R <old log> <new log> ",
    "[<old capture dir> <new capture dir>]",
    call. = FALSE
  )
}

status_names <- c(
  "-1" = "NOT_RUN",
  "0" = "CONVERGED",
  "1" = "INITIAL_RR_FAILED",
  "2" = "STALLED",
  "3" = "RESIDUAL_COLLAPSED",
  "4" = "P_COLLAPSED",
  "5" = "RR_FAILED",
  "6" = "MAXITER",
  "7" = "BAD_SEED"
)

read_log <- function(path) {
  x <- readLines(path)
  parts <- strsplit(x, " ", fixed = TRUE)
  list(
    path = path,
    lines = x,
    case = vapply(parts, `[[`, "", 1L),
    event = vapply(parts, `[[`, "", 2L),
    parts = parts
  )
}
num <- function(log, event, pos) {
  hit <- log$event == event
  as.numeric(vapply(log$parts[hit], `[[`, "", pos))
}
show_table <- function(tab, indent = "  ") {
  for (nm in names(tab)) {
    cat(sprintf("%s%-28s %d\n", indent, nm, tab[[nm]]))
  }
}

lobpcg_runs <- function(log) {
  hit <- log$event == "lobpcg"
  data.frame(
    case = log$case[hit],
    status = as.integer(vapply(log$parts[hit], `[[`, "", 4L)),
    iterations = as.integer(vapply(log$parts[hit], `[[`, "", 6L)),
    fail_iter = as.integer(vapply(log$parts[hit], `[[`, "", 8L)),
    k = as.integer(vapply(log$parts[hit], `[[`, "", 10L))
  )
}

summarize <- function(log, is_old) {
  cat(sprintf("\n== %s\n", log$path))
  cat(sprintf(
    "%d events in %d cases\n",
    length(log$lines),
    length(unique(log$case))
  ))
  stopifnot(length(log$lines) > 0L)
  cat("events:\n")
  show_table(table(log$event))

  runs <- lobpcg_runs(log)
  cat(sprintf("lobpcg_solve() calls: %d, ended as\n", nrow(runs)))
  show_table(table(status_names[as.character(runs$status)]))
  cat(sprintf(
    "  eigenpairs asked of LOBPCG: %s\n",
    paste(sort(unique(runs$k)), collapse = " ")
  ))

  ord <- num(log, "syevd", 3L)
  stopifnot(length(ord) > 0L)
  mid <- ord >= 34 & ord <= 79
  cat(sprintf(
    "dsyevd: %d calls, orders %d to %d; %d of order 34 to 79\n",
    length(ord),
    min(ord),
    max(ord),
    sum(mid)
  ))
  if (is_old) {
    query <- num(log, "syevd", 5L)
    held <- num(log, "syevd", 7L)
    longer <- held > query
    cat(sprintf(
      "  given a longer workspace than its own query: %d calls, %d of them of order 34 to 79, in cases: %s\n",
      sum(longer),
      sum(longer & mid),
      paste(
        unique(log$case[log$event == "syevd"][longer & mid]),
        collapse = " "
      )
    ))
    cat(sprintf(
      "triangular solve: %d calls\nfirst Rayleigh-Ritz over X alone: %d\n",
      sum(log$event == "solve"),
      sum(log$event == "xonly_rr")
    ))
    fl <- log$lines[log$event == "syevd_failed"]
    cat(sprintf("dsyevd failures: %d, by info\n", length(fl)))
    if (length(fl)) {
      show_table(table(sub("^[^ ]+ syevd_failed ", "info ", fl)))
    }
  } else {
    rc <- num(log, "solve", 6L)
    ap <- num(log, "solve", 8L)
    cs <- log$case[log$event == "solve"]
    if (length(rc)) {
      ## An injected "solve" failure never enters the wrapper, which leaves
      ## rcond at 0 with the flag unset; it is not a solve.
      real <- !(rc == 0 & ap == 0)
      exact <- real & ap == 0
      cat(sprintf(
        "triangular solve: %d calls, %d by the dgelsd fallback (cases: %s); smallest rcond of the others %.3g, in case %s\n",
        sum(real),
        sum(ap != 0),
        paste(unique(cs[ap != 0]), collapse = " "),
        min(rc[exact]),
        cs[exact][which.min(rc[exact])]
      ))
    }
    cat(sprintf(
      "first Rayleigh-Ritz over X alone: %d\n",
      sum(log$event == "xonly_rr")
    ))
    st <- log$lines[log$event == "status"]
    cat(sprintf("non-zero si_la statuses: %d\n", length(st)))
    if (length(st)) {
      show_table(table(sub("^[^ ]+ ", "", st)))
    }
  }
  inj <- log$lines[log$event == "injected"]
  cat(sprintf("injected failures: %d\n", length(inj)))
  runs
}

old <- read_log(args[[1]])
new <- read_log(args[[2]])
## A file that calls Armadillo for its dsyevd logs the workspace it held and
## its query; one that calls si_la logs the order alone. Each log is read by
## its own format, so either build may be on either side.
arma_log <- function(log) {
  any(log$event == "syevd" & lengths(log$parts) >= 7L)
}
runs_old <- summarize(old, is_old = arma_log(old))
runs_new <- summarize(new, is_old = arma_log(new))

## Per case, how each build's lobpcg_solve() calls ended, in order.
key <- function(runs) {
  vapply(
    split(runs, factor(runs$case, levels = unique(runs$case))),
    function(r) {
      paste(r$status, r$iterations, r$fail_iter, sep = "/", collapse = " ")
    },
    ""
  )
}
ko <- key(runs_old)
kn <- key(runs_new)
cases <- union(names(ko), names(kn))
differ <- cases[
  !(cases %in% names(ko)) |
    !(cases %in% names(kn)) |
    ko[cases] != kn[cases]
]
differ <- differ[!is.na(differ)]
cat(sprintf(
  "\n== cases whose lobpcg_solve() calls (status/iterations/fail_iter) differ between the two: %d of %d\n",
  length(differ),
  length(cases)
))
for (cs in differ) {
  o <- strsplit(ko[cs], " ", fixed = TRUE)[[1]]
  n <- strsplit(kn[cs], " ", fixed = TRUE)[[1]]
  len <- min(length(o), length(n))
  first <- which(o[seq_len(len)] != n[seq_len(len)])
  cat(sprintf(
    "  %-20s old %d calls, new %d; first difference at call %s: old %s, new %s\n",
    cs,
    length(o),
    length(n),
    if (length(first)) first[[1]] else len + 1L,
    if (length(first)) o[[first[[1]]]] else "-",
    if (length(first)) n[[first[[1]]]] else "-"
  ))
}

## What each injected failure led to, in each build.
ended <- function(runs, cs) {
  r <- runs[runs$case == cs, ]
  bad <- r[r$status != 0L, ]
  sprintf(
    "%d solves, not converged: %s",
    nrow(r),
    if (nrow(bad)) {
      paste(
        sprintf(
          "%s@%d",
          status_names[as.character(bad$status)],
          bad$fail_iter
        ),
        collapse = " "
      )
    } else {
      "none"
    }
  )
}
inj_cases <- unique(new$case[new$event == "injected"])
cat(sprintf("\n== injected cases: %d\n", length(inj_cases)))
for (cs in inj_cases) {
  what <- sub(
    "^[^ ]+ injected ",
    "",
    new$lines[new$case == cs & new$event == "injected"]
  )
  e_new <- ended(runs_new, cs)
  e_old <- ended(runs_old, cs)
  cat(sprintf(
    "  %-22s %-12s %s; X-only RR %d%s\n",
    cs,
    what[[1]],
    e_new,
    sum(new$case == cs & new$event == "xonly_rr"),
    if (identical(e_new, e_old)) "" else paste0("  [old: ", e_old, "]")
  ))
}

## Cases that ended in an error, with what each build raised.
if (length(args) == 4L) {
  rd <- function(dir, cs) {
    f <- file.path(dir, paste0(cs, ".rds"))
    if (!file.exists(f)) {
      return("missing")
    }
    v <- readRDS(f)$value
    if (inherits(v, "ref_error")) {
      sprintf("[%s] %s", v$class[[1]], strsplit(v$message, "\n")[[1]][[1]])
    } else {
      sprintf("a result, %d NA, %d NaN", sum(is.na(v)), sum(is.nan(v)))
    }
  }
  all_cases <- sub(
    "\\.rds$",
    "",
    setdiff(list.files(args[[4]], pattern = "\\.rds$"), "_meta.rds")
  )
  cat("\n== cases that ended in an error in either build\n")
  n_err <- 0L
  for (cs in all_cases) {
    o <- rd(args[[3]], cs)
    n <- rd(args[[4]], cs)
    if (startsWith(o, "[") || startsWith(n, "[")) {
      n_err <- n_err + 1L
      cat(sprintf("  %s\n    old: %s\n    new: %s\n", cs, o, n))
    }
  }
  cat(sprintf("%d of %d cases\n", n_err, length(all_cases)))
}
