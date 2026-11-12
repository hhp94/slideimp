# CLAUDE.md

slideimp is a released R package: K-NN and PCA imputation for numeric
matrices, with grouped and sliding-window strategies and hyperparameter
tuning. It targets high-dimensional numeric data, Illumina DNA methylation
arrays and whole-genome bisulfite data in particular.

It is a CRAN package. Its exported API, its error surface and its
documentation are public, and callers pin on all three - so a change to any
of them is a deliberate act, not a side effect of tidying. That is the fact
most arguments here come back to.

`slideimp.extra` is the one companion package. It supplies Illumina array
manifests through `ilmn_manifest()` and `set_slideimp_path()`, which
`prep_groups()` accepts as its `group` argument. It is not in `Imports` or
`Suggests` and slideimp never requires it: every path that uses it is
optional and must degrade cleanly when it is absent.

This file holds invariants only: what to do, stated so it stays true as the
code changes. A measurement, a date, a pass count, or an account of something
that was once broken does not belong here no matter how important it is.

## House rules

**ASCII only, everywhere except `dev/`.** This is an accessibility
requirement, not a style preference: dashes of differing lengths and curly
quotes are hard for the maintainer to read. Use `-` and `--` for dashes,
straight quotes, `x` for multiplication, `>=` and `<=`, `->` for arrows, `u`
for micro, `^2` for superscripts. Greek letters are acceptable in genuine
mathematical context and nowhere else. The test is WHO READS IT: code,
comments, roxygen, commit messages and chat replies are all covered.
`man/`, `docs/` and `README.md` are generated, so the rule is enforced on
their sources. `LICENSE.md` is exempt as verbatim license text; `dev/` is
exempt so prose is never blocked, but write ASCII there anyway. Check with:

```bash
git grep -IlP '[^\x00-\x7F]' -- . ':!dev' ':!man' ':!docs' ':!LICENSE.md'
```

It prints offending file names, so no output means clean. The one legitimate
hit is `README.md` when an evaluated chunk captures the glyphs cli prints,
which section 7 of [dev/WRITING.md](dev/WRITING.md) keeps; any other
character there is a leak to fix in `README.Rmd`.

**Never assume - measure when it can be measured.** Do not guess where time
goes, whether an optimization helped, how much of a matrix is missing, or
whether a dependency has the argument you are about to call. Profile it, time
it, count it, print it. This applies to claims about behaviour as much as
about speed: if an assertion can be turned into a script that prints a number,
write the script. Measure the real thing, not a paraphrase of it. When
something genuinely cannot be measured, say so explicitly and label it an
assumption rather than letting it pass as fact.

**Never run anything inline in PowerShell.** Write a script - `.ps1`, `.py`,
or `.R` - then call it. No `-Command` one-liners, no inline `Rscript -e`, no
multi-statement pipelines typed straight into the shell. An inline command
leaves nothing behind: it cannot be re-read, corrected, re-run, or reviewed
after the fact, and its quoting breaks in ways that are invisible until they
silently do the wrong thing. Throwaway inspection scripts go in `dev/`.

**No heredocs or here-strings for writing files.** No `cat > file <<EOF`, no
`Set-Content` with a here-string, no `printf` of a multi-line body. On this
Windows machine a heredoc crosses Git Bash's quoting, PowerShell's parser and
CRLF line endings, and it fails as a truncated file or an "unexpected EOF",
not as a clean error. Write file bodies with the Write and Edit tools;
otherwise use whatever tool fits the job, `sed -i` included.

**C++ under `src/`: read [dev/CPP.md](dev/CPP.md) first.** It states two bans
in full. No hand-rolled loop unrolling, strip-mining for vectorization, or
OpenMP pragmas: write the simple loop and let the compiler decide, and spend
parallelism at the chain level through RcppThread and mirai. No thread
pinning from C++: pin from R with RhpcBLASctl, best-effort, and verify a pin
by CPU time over wall time, never by asking the library.

**Reach for the solution a dependency publishes; no low-level hacks.** When a
dependency has a designed, documented answer to a problem - RcppThread for
thread-safe interruption and stop, R's own APIs for anything touching R,
collapse's primitives for grouped statistics - use that, even when a
hand-rolled shortcut looks smaller. A hack against a library's internals
compiles today, breaks on that library's next release, and is exactly what
CRAN review exists to be unhappy about. If no published solution exists, that
is a finding worth writing down - and usually a reason to redesign - not a
license to reach into internals.

**Errors and messages go through cli.** `cli::cli_abort()`, `cli_warn()`,
`cli_inform()`, with inline markup (`{.arg x}`, `{.code y}`, `{.fn z}`) and
pluralization. `cli` is a hard dependency and the error surface is public API:
callers pin on wording and condition class, so a message changes deliberately,
not in passing. Input validation goes through `checkmate` with an explicit
`.var.name` so the message names the user's argument rather than an internal
one. How the text of a message, a help page or a vignette is written lives in
[dev/WRITING.md](dev/WRITING.md); read it before writing any of them.

**Full roxygen, markdown on.** Every exported function carries a title, a
description, `@param`, `@returns`, and runnable `@examples`. Shared prose
is written once and pulled in with `@inheritSection` or `@inheritParams`
rather than duplicated - a claim that appears in five help pages has to stay
true in five places otherwise. Run `devtools::document()` after touching any
roxygen block; `man/` is generated and never hand-edited.

**Run `air format .` after every change to an R file, before calling the work
done.** `air` is a command-line program configured by `air.toml`; run it from
the repository root on the whole tree, not on single files, and in the same
change as the edit, or its reformat lands in an unrelated commit and buries
the real change. It also reformats every R script in `dev/`, tracked or
ignored: a tracked `dev/` script showing as modified afterwards was
reformatted, not un-ignored. It moves lines, so re-check the `file:line`
citations in the structural maps after it runs. A diff from air that changes
what the code does is a bug to report, not to commit. Do not fight its output
by hand: a wrong layout is fixed in `air.toml`, agreed with the maintainer.

**Never touch `NEWS.md` without approval.** Its style is the maintainer's and
is not reconstructable from the diff. Propose the entry in chat and wait.

**Numerical tests report absolute and relative error, and assert on the worst
case.** `max(abs(a - b))` and `max(abs(a - b) / abs(b))` over the elements -
not correlation, R-squared, RMSE or any mean, which hide the one
catastrophically wrong cell in a million right ones. Report both, since
absolute error misleads across column scales and relative error blows up near
zero; state which elements the relative figure was taken over, and assert on
the maximum.

**A comparison against a STORED reference asserts on a tolerance, never
`identical()`.** A stored baseline is bytes from whichever build wrote it, so
asserting bitwise equality pins the compiler rather than the code. Report
bitwise-ness without asserting on it, so the log shows the day it stops being
true. `identical()` stays correct for two objects built in the SAME process,
and for integer maps and structural quantities anywhere.

**A test that cannot fail is worse than no test.** Check that a new test's
inputs actually exercise the branch it claims to. If it asserts something
about a subset, assert the subset is non-empty in the same test - `max()` over
nothing returns `-Inf` and passes silently. After relaxing a comparator, prove
it can still fail. When you fix a bug, add the regression test pinned to that
bug.

**Anything a tracked file cites must itself be tracked.** The maintainer
works from more than one machine, and a session that finds a cited file
missing reasons about why and produces a confident wrong answer. In `dev/`,
add the file's `!` line to `dev/.gitignore` in the same commit that first
cites it. Elsewhere, do not add an ignore rule for something the maintainer
produced or fetched by hand: an ignore entry is invisible in every later diff,
and its cost lands on a future session. Large regenerable DATA is the
exception and stays out. If a file genuinely cannot be tracked, say so and
ask; never decide it silently.

## Internals

The R-to-C++ boundary of `pca_imp()` and `knn_imp()` is not obvious from
either side alone. Structural maps of both chains, with file-and-line
citations throughout:

- [dev/pca_imp_structure.md](dev/pca_imp_structure.md) - R entry point, Rcpp
  boundary, EM outer loop, eigensolver layer (LOBPCG warm start included),
  down to BLAS.
- [dev/knn_imp_structure.md](dev/knn_imp_structure.md) - R entry point, Rcpp
  boundary, distance over incomplete rows, neighbor selection, weighting, and
  the `slideimp_infeasible` abort sites.

Read the relevant one before changing anything in `src/`, or before changing
an argument in `R/pca_imp.R` or `R/knn_imp.R` that reaches C++. They are
descriptive maps, not specifications: when a map and the code disagree, the
code is right and the map needs fixing in the same change.

## Layout

- `R/`, `src/` - package source. `src/` is Rcpp with RcppArmadillo, plus
  RcppEigen in the multiple-imputation sampler, and RcppThread for
  parallelism; `R/RcppExports.R` and `src/RcppExports.cpp` are generated.
- `R/dev-utils.R` - dev-only helpers, unexported and build-ignored. See
  "Working commands" below.
- `tests/testthat/` - testthat edition 3, run against `load_all()` so
  internals are visible. `tests/testthat/dev/` is untracked fixture space for
  the optional `slideimp.extra` manifest tests, holding regenerable vendor
  downloads: testthat runs from `tests/testthat/`, which is why
  `set_slideimp_path("dev")` resolves there.
- `vignettes/`, `README.Rmd`, `_pkgdown.yml` - user-facing documentation.
  `README.md` and `docs/` are generated; edit the sources.
- `dev/` - design notes, structural maps, dev scripts. Build-ignored so it
  never ships, and tracked opt-in: a file is ignored unless it has a `!` line
  in `dev/.gitignore`. Code here is exploratory and nothing in the package may
  depend on it; anything that turns durable is promoted to `R/` or `src/`.
- `dev/GSE286313.EPICv1_chr21.qs2` - real data, for measuring a claim on
  something other than `sim_mat()` output: Illumina EPIC v1 beta values for
  the chromosome 21 probes of GEO series GSE286313, with their real
  missingness (all-`NA` probes, probes over the default `colmax`, scattered
  `NA` in samples). Read it with `qs2::qs_read()` and `t()` it, since it is
  stored probes x samples. It has probe IDs and no genomic coordinates, so it
  cannot supply a `location` vector for `slide_imp()`. Dev scripts only:
  `qs2` is not a declared dependency, so nothing under `R/`, `tests/` or
  `vignettes/` may read it.
- [dev/to-do.md](dev/to-do.md) - work that is agreed on but not yet done.
  Read it before starting anything in `R/` or `src/`: an entry there may
  already cover the change, or may constrain how it has to be made. Add an
  entry when a change is decided on but deferred, and DELETE the entry when
  the change lands rather than marking it done.

## Working commands

Per the PowerShell rule, these go in a script that you then call.

**R's home directory depends on the shell.** R on Windows takes `HOME` from
the environment when it is set, so under Git Bash - the Bash tool here - `~`
is the user profile, while under PowerShell it is the Documents folder. An
`.Renviron` or `R_LIBS_USER` setting under `~` is seen from one shell and not
the other, and it has to be in place BEFORE R starts: `.libPaths()` from inside
a running script is too late. If a package that should be installed reports
`there is no package called '<x>'`, print `Sys.getenv("HOME")` and
`.libPaths()` before concluding it is missing.

**Load the package with `load_all1()`, not bare `load_all()`.** It decides
which compile-time macros the DLL is built with, and cleans the DLL first so
a change of macros always rebuilds. The header of `R/dev-utils.R` explains
each flag and the other helpers. Never mix a `load_all1()` build with a bare
`load_all()` one inside a single comparison.

```r
load_all1()                            # timer on, diagnostics off
load_all1(timer = FALSE)               # no rcpptimer needed
load_all1(debug = TRUE)                # -DPCA_IMP_DIAGNOSTICS=1

devtools::document()                   # after touching any roxygen block
devtools::test()
check_manual()                         # test() with MANUAL_TESTS=true
devtools::check()                      # before proposing anything release-shaped
```
