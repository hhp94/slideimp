# Review of the wrapper layer (dev/review_eigen_wrappers.md): the probe
# dev/review_eigen_wrappers_probe.R runs cases that never return, each in a
# child R process with a time limit. On Windows the limit ends Rscript.exe
# and can leave the Rterm.exe it started spinning inside LAPACK. This lists
# the R processes on the machine and, with --kill, ends only those whose
# command line names that probe script with a --case= argument.
#   Rscript dev/review_eigen_wrappers_procs.R [--kill]
# Needs the ps package.

kill <- "--kill" %in% commandArgs(trailingOnly = TRUE)
procs <- ps::ps()
is_r <- grepl("^(Rterm|Rscript|R)(\\.exe)?$", procs$name, ignore.case = TRUE)
procs <- procs[is_r & procs$pid != Sys.getpid(), , drop = FALSE]
cat(sprintf("%d other R processes\n", nrow(procs)))
for (i in seq_len(nrow(procs))) {
  cmd <- tryCatch(
    paste(ps::ps_cmdline(procs$ps_handle[[i]]), collapse = " "),
    error = function(e) "(command line not readable)"
  )
  mine <- grepl("review_eigen_wrappers_probe.R", cmd, fixed = TRUE) &&
    grepl("--case=", cmd, fixed = TRUE)
  cat(sprintf(
    "pid %6d %-12s started %s cpu %6.0f s  %s\n  %s\n",
    procs$pid[[i]],
    procs$name[[i]],
    format(procs$created[[i]], "%H:%M:%S"),
    sum(procs$user[[i]], procs$system[[i]], na.rm = TRUE),
    if (mine) "LEFT BY THE PROBE" else "",
    substr(cmd, 1L, 240L)
  ))
  if (mine && kill) {
    ok <- tryCatch(
      {
        ps::ps_kill(procs$ps_handle[[i]])
        TRUE
      },
      error = function(e) FALSE
    )
    cat(sprintf("  %s\n", if (ok) "ended" else "could not be ended"))
  }
}
