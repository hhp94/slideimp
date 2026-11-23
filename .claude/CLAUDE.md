# CLAUDE.md

slideimp is an R package on CRAN: K-NN and PCA imputation for numeric
matrices, with grouped and sliding-window strategies and hyperparameter
tuning. It targets high-dimensional numeric data, Illumina DNA methylation
arrays and whole-genome bisulfite data in particular.

Its exported API, its error surface and its documentation are public, and
callers pin on all three, so a change to any of them is deliberate, never a
side effect of tidying.

`slideimp.extra` is the one companion package. It supplies Illumina array
manifests (`ilmn_manifest()`, `set_slideimp_path()`) and registers a
`slideimp_resolve_group()` method, so `prep_groups()` and `group_imp()` take a
manifest, or a platform name such as `"EPICv2"`, as `group`. It is not in
`Imports` or `Suggests` and slideimp never requires it: every path that uses
it is optional and must degrade cleanly when it is absent.

This file holds invariants only: what to do, stated so it stays true as the
code changes. A measurement, a date, a pass count or an account of something
that was once broken does not belong here, however important.

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

It prints offending files; no output means clean. `README.md` may hit
legitimately where an evaluated chunk captures the glyphs cli prints, which
section 7 of [dev/WRITING.md](../dev/WRITING.md) keeps; any other character
there is a leak to fix in `README.qmd`.

**Never assume - measure when it can be measured.** Do not guess where time
goes, whether an optimization helped, how much of a matrix is missing, or
whether a dependency has the argument you are about to call: profile, time,
count or print it. Claims about behaviour count as much as claims about
speed: if an assertion can be turned into a script that prints a number,
write the script, and measure the real thing, not a paraphrase of it. What
genuinely cannot be measured is labelled an assumption, never passed off as
fact.

**Never run anything inline in PowerShell.** Write a script - `.ps1`, `.py`
or `.R` - then call it: no `-Command` one-liners, no inline `Rscript -e`, no
multi-statement pipelines typed into the shell. An inline command cannot be
re-read, corrected, re-run or reviewed afterwards, and its quoting breaks
silently. Throwaway inspection scripts go in `scratch/` (see "Layout").

**No heredocs or here-strings for writing files.** No `cat > file <<EOF`, no
`Set-Content` with a here-string, no `printf` of a multi-line body. On this
Windows machine a heredoc crosses Git Bash's quoting, PowerShell's parser and
CRLF line endings, and fails as a truncated file or an "unexpected EOF", not
a clean error. Write file bodies with the Write and Edit tools; otherwise use
whatever tool fits, `sed -i` included.

**C++ under `src/`: read [dev/CPP.md](../dev/CPP.md) first.** It states four
bans in full. No hand-rolled loop unrolling, strip-mining for vectorization, or
OpenMP pragmas: write the simple loop and let the compiler decide, and spend
parallelism at the chain level through RcppThread and mirai. No thread
pinning from C++: pin from R with RhpcBLASctl, best-effort, and verify a pin
by CPU time over wall time, never by asking the library. No dimnames in C++:
a matrix returns with none, and the R caller names it right after the call.
No `Rcpp::List` built from a `SEXP`: a list from R arrives as `SEXP`, is
checked with `TYPEOF()` and read with `VECTOR_ELT()`, and
`dev/check_rcpp_list.R` checks the compiled objects.

**Run `clang-format -i src/*.cpp src/*.h` after every change to C++ under
`src/`, before calling the work done.** Run it from the repository root, in
the same change as the edit, with clang-format 23.1.5 (`uv tool install
clang-format==23.1.5`): the repository does not pin it, and another release
lays code out differently. `.clang-format` holds the layout and
`.clang-format-ignore` the files it leaves alone: `src/RcppExports.cpp`,
and `src/Spectra/`, whose layout is Spectra's. Like air,
it splits and joins lines, so run the citation checker after it. A diff from
it that changes more than whitespace is a bug to report, not to commit; a
wrong layout is fixed in `.clang-format`, agreed with the maintainer.

**`obj` is a double matrix, and nothing in the package converts it.** Every
exported function that takes the data matrix accepts a double matrix only,
and refuses any other through checkmate (`mode = "double"`, with `.var.name
= "obj"`). No code under `R/` or `src/` converts `obj` or a slice of it: no
`storage.mode(obj) <-`, `as.double()` or `as.numeric()` on it. A conversion
copies the whole matrix, so a caller who needs one makes it once, before the
call, as the package help page says. Rcpp converts silently: an
`Rcpp::NumericMatrix` built from an integer matrix holds a new double copy.
So the type is checked in R before `obj` reaches C++, and C++ that cannot
rely on that check takes a `SEXP` and refuses anything but `REALSXP`, as
`col_scan_cpp()` does.

**Reach for the solution a dependency publishes; no low-level hacks.** When a
dependency has a designed, documented answer - RcppThread for thread-safe
interruption and stop, R's own APIs for anything touching R, collapse's
primitives for grouped statistics - use it, even when a hand-rolled shortcut
looks smaller. A hack against a library's internals breaks on its next
release and is what CRAN review flags. If no published solution exists, write
that down as a finding, usually a reason to redesign, not a license to reach
into internals.

**New errors and messages are base R; the cli calls already here stay.** New
code raises with `stop()` and `warning()`, `call. = FALSE`, and reports with
`message()`, even in a function whose other messages are cli. A condition a
caller is expected to catch carries a `slideimp_*` class through
`errorCondition()` or `warningCondition()`. The existing `cli_abort()`,
`cli_warn()` and `cli_inform()` calls are not converted: moving one to base R
changes its condition class as well as its text, and an edit to one for a
real reason keeps it in cli. Callers pin on wording and condition class, so a
message changes deliberately, not in passing. Input validation goes through
`checkmate` with an explicit `.var.name`, so the message names the user's
argument rather than an internal one. How a message, a help page or a
vignette is written lives in [dev/WRITING.md](../dev/WRITING.md); read it
before writing any of them.

**Full roxygen, markdown on.** Every exported function carries a title, a
description, `@param`, `@returns` and runnable `@examples`. Shared prose is
written once and pulled in with `@inheritSection` or `@inheritParams`, since a
claim repeated in five help pages has to stay true in five places. Run
`devtools::document()` after touching any roxygen block; `man/` is generated
and never hand-edited.

**Run `air format .` after every change to an R file, before calling the work
done.** Run it from the repository root on the whole tree, not on single
files, and in the same change as the edit, or its reformat lands in an
unrelated commit and buries the real change. It is configured by `air.toml`
and also reformats the R scripts in `dev/`, ignored ones included. It can
split or join lines, which breaks code fragments the notes cite, so run the
citation checker (see "Internals") after it. A diff from air that changes
what the code does is a bug to report, not to commit. Do not fight its output
by hand: a wrong layout is fixed in `air.toml`, agreed with the maintainer.

**`NEWS.md` is written once per release, from the diff.** Do not write or
propose a NEWS entry as part of a change. Before a release, the whole section
is drafted in one pass from the diff since the last release tag, with the
maintainer, whose style it is; outside that pass, never touch `NEWS.md`. When
a change has a user-visible effect that would be easy to miss in that diff -
a removed argument, a changed default, a new error - add one line for it to
the NEWS entry in [dev/to-do.md](../dev/to-do.md), creating the entry if there
is none. That line is a reminder for the draft, not NEWS text.

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
true.

**Doubles computed through BLAS or LAPACK are compared on a tolerance, even
two built in the same process.** A BLAS need not return the same bits for the
same call twice, so `identical()` between two such results pins the library,
not the code: assert with `expect_worst_err()` at a tolerance sized to
rounding, and report bitwise-ness without asserting on it. This covers every
result that passes through the PCA chain or the sampler, `slide_imp()`,
`group_imp()` and `tune_imp()` running PCA included. `identical()` stays
correct for integer maps and structural quantities anywhere, for values
copied rather than computed (the observed cells of a result), and for two
doubles built in one process by code that reaches no BLAS or LAPACK.

**A test that cannot fail is worse than no test.** Check that a new test's
inputs actually exercise the branch it claims to. If it asserts something
about a subset, assert the subset is non-empty in the same test - `max()` over
nothing returns `-Inf` and passes silently. After relaxing a comparator, prove
it can still fail. When you fix a bug, add the regression test pinned to that
bug.

**Anything a tracked file cites must itself be tracked.** The maintainer
works from more than one machine, and when a session finds a cited file
missing, it reasons about why and produces a confident wrong answer. In `dev/`,
add the file's `!` line to `dev/.gitignore` in the same commit that first
cites it. Elsewhere, add no ignore rule for something the maintainer produced
or fetched by hand: an ignore entry is invisible in every later diff. Two
things stay out regardless: large regenerable data, and `local/` (see
"Layout"), which a tracked file may mention only as a folder on one machine.
If a file genuinely cannot be tracked, say so and ask; never decide it
silently.

**`src/Spectra/` is vendored and never edited by hand.** It is a
third-party library, re-vendored whole by the script that
[dev/vendored.md](../dev/vendored.md) names. Every other file in the
repository is slideimp's own and in scope for every pass.

## Internals

The R-to-C++ boundary of `pca_imp()`, `knn_imp()` and the
multiple-imputation sampler is not obvious from either side alone.
Structural maps of the three chains, citing the code throughout:

- [dev/pca_imp_structure.md](../dev/pca_imp_structure.md) - R entry point, Rcpp
  boundary, EM outer loop, eigensolver layer (LOBPCG warm start included),
  down to BLAS.
- [dev/knn_imp_structure.md](../dev/knn_imp_structure.md) - R entry point, Rcpp
  boundary, distance over incomplete rows, neighbor selection, weighting, and
  the `slideimp_infeasible` abort sites.
- [dev/mipca_structure.md](../dev/mipca_structure.md) - the sampler's doors
  and `mipca_core()`, the column build, Rcpp boundary, one iteration of the
  data-augmentation loop, finalize, resume, and the random-number streams.

Read the relevant one before changing anything in `src/`, or an argument in
`R/pca_imp.R`, `R/knn_imp.R` or `R/mi_pca_*.R` that reaches C++. They describe, they do not
specify: when a map and the code disagree, the code is right and the map is
fixed in the same change.

The maps, and every other note in `dev/`, cite code by file and by a fragment
of one line that occurs exactly once in that file (`R/knn_imp.R`:
`knn_imp <- function(`), never by line number. A line number moves with every
edit above it and still lands on a real line, so it goes wrong without
failing; a stale fragment matches nothing. `uv run dev/check_map_anchors.py`
checks every citation in the three maps and in this file; its header gives the
forms, and naming another note on its command line checks that note too. Run
it after changing code that a map cites, and fix what it reports in the same
change.

## Layout

- `R/`, `src/` - package source. `src/` is plain Rcpp over R's own types,
  and RcppThread for parallelism; the PCA chain reaches R's BLAS and LAPACK
  through the wrapper layer in `src/si_lapack.h`. The multiple-imputation
  sampler in `src/mipca.cpp` adds RcppEigen and the Spectra headers vendored
  in `src/Spectra/`. `R/RcppExports.R` and `src/RcppExports.cpp` are
  generated.
- `tests/testthat/` - testthat edition 3, run against `load_all()` so
  internals are visible. Nothing under `tests/` may call a package that
  `DESCRIPTION` does not declare, even behind `skip_if_not_installed()`: R CMD
  check reports it as an undeclared dependency. A check that needs one is a
  script in `dev/`, run by hand: `dev/check_impute_knn.R` for `impute`,
  `dev/check_slideimp_extra.R` for `slideimp.extra`, whose header says what
  it needs.
  A test's tier belongs to its file, set on the file's first line with
  nothing above it. CRAN runs `test-smoke.R` alone, which has no skip: one
  smoke per export and per compiled path, each on a small input, at most 2
  cores and no mirai daemons, asserting shape, attributes and one value
  with a tolerance, under 0.2 s. The local tier, every file that starts
  with `skip_on_cran()`, is what `devtools::test()` runs; a new test goes
  there unless it is a smoke. The on-request tier, `test-<function>-manual.R`
  files starting with `skip_if_not_manual()`, holds what is slow or needs
  more than the session (comparisons with another package, mirai daemons,
  the PCA kernel fuzz); `check_manual()` runs it, before a CRAN submission
  and after any change to the kernel. The error surface is pinned in the
  local tier, not in the smokes.
- `vignettes/`, `README.qmd`, `_pkgdown.yml` - user-facing documentation.
  `README.md` and `docs/` are generated; edit the sources.
- `dev/` - design notes, structural maps, and dev scripts kept for reuse or
  cited by a note. Build-ignored so it never ships, and tracked opt-in: a file
  is ignored unless it has a `!` line in `dev/.gitignore`. Code here is
  exploratory and nothing in the package may depend on it; anything that
  turns durable is promoted to `R/` or `src/`.
  Scripts may be R or Python; run Python ones with `uv run <script>`, not
  `python`: uv supplies the interpreter, and a script lists the packages it
  needs in an inline `# /// script` header.
- `dev/dev-utils.R` - the dev helpers, `load_all1()` and `check_manual()`
  among them. Not part of the package, so `source()` it first.
- `dev/GSE286313.EPICv1_chr21.qs2` - real data, for measuring a claim on
  something other than `sim_mat()` output: Illumina EPIC v1 beta values for
  the chromosome 21 probes of GEO series GSE286313, with their real
  missingness (all-`NA` probes, probes over the default `colmax`, scattered
  `NA` in samples). Stored probes x samples: read it with `qs2::qs_read()`
  and `t()` it. It has probe IDs but no genomic coordinates, so it cannot
  supply `location` for `slide_imp()`. Dev scripts only: `qs2` is not a
  declared dependency, so nothing under `R/`, `tests/` or `vignettes/` may
  read it.
- [dev/to-do.md](../dev/to-do.md) - work that is agreed on but not yet done.
  Read it before starting anything in `R/` or `src/`: an entry there may
  already cover the change, or constrain how it is made. Add an entry when a
  change is decided on but deferred, and DELETE it when the change lands
  rather than marking it done.
- [dev/DECISIONS.md](../dev/DECISIONS.md) - closed decisions, stated as facts
  about the code. Read it before proposing a change in `R/` or `src/`, and do
  not reopen an entry on your own initiative: not as a to-do entry, a review
  finding or a drive-by fix. Only the maintainer reopens one, and the entry
  changes in the same change as the code.
- `local/` - machine-local, git-ignored and build-ignored by design: it names
  machines, addresses and keys, and this repository is public. It exists on
  one machine only, so its absence anywhere else is expected, not a loss.
  Nothing in the package, `tests/` or `dev/` may depend on it.
- `scratch/` - the default place for everything an agent makes that is not
  asked to last: probe and wrapper scripts, saved CI logs, downloaded or
  unpacked packages, measurement output. Git-ignored and build-ignored, never
  tracked, and never sourced, read or cited from anything tracked, `dev/`
  included, so it can be emptied at any time without breaking anything. A
  script that turns out worth keeping moves to `dev/` with its `!` line, in
  the change that first cites it.

## Working commands

Per the PowerShell rule, these go in a script that you then call.

**R's home directory depends on the shell.** R on Windows takes `HOME` from
the environment when it is set, so under Git Bash (the Bash tool here) `~` is
the user profile, and under PowerShell it is the Documents folder. An
`.Renviron` or `R_LIBS_USER` setting under `~` is seen from one shell and not
the other, and it must be in place BEFORE R starts: `.libPaths()` inside a
running script is too late. If a package that should be installed reports
`there is no package called '<x>'`, print `Sys.getenv("HOME")` and
`.libPaths()` before concluding it is missing.

**Load the package with `load_all1()`, not bare `load_all()`.** It decides
which compile-time macros the DLL is built with, and cleans the DLL first so
a change of macros always rebuilds; the header of `dev/dev-utils.R` explains
each flag and the other helpers. Never mix a `load_all1()` build with a bare
`load_all()` one inside a single comparison.

```r
source("dev/dev-utils.R")              # load_all1(), check_manual(), ...
load_all1()                            # timer on, diagnostics off
load_all1(timer = FALSE)               # no rcpptimer needed
load_all1(debug = TRUE)                # -DPCA_IMP_DIAGNOSTICS=1

devtools::document()                   # after touching any roxygen block
devtools::test()
check_manual()                         # test() with MANUAL_TESTS=true
devtools::check()                      # before proposing anything release-shaped
```

**Local CI on rhub's Linux images** is described under `local/`, on the one
machine that has it; read it before a run. A run takes the shared VM, so
start one only when asked. Which CRAN check each platform stands for, and
the run tiers of a release, are in [dev/ci.md](../dev/ci.md).
