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
requirement, not a style preference. Dashes of differing lengths and curly
quotes are hard for the maintainer to read, so plain ASCII hyphens and `--`
are the only dash forms allowed. Use `-` and `--` for dashes, straight quotes,
`x` for multiplication, `>=` and `<=`, `->` for arrows, `u` for micro, `^2`
for superscripts. Greek letters are acceptable in genuine mathematical context
and nowhere else.

This covers code, comments, roxygen, commit messages, and anything written
back to the user in chat. The test is WHO READS IT. `man/` and `docs/` are
generated, so the rule is enforced at the roxygen source rather than on the
output. `LICENSE.md` is exempt because it is verbatim license text that must
not be altered. `dev/` is exempt so prose is never blocked - not because
the characters are readable there, so write ASCII in it anyway.

Check it with:

```bash
git grep -IlP '[^\x00-\x7F]' -- . ':!dev' ':!man' ':!docs' ':!LICENSE.md'
```

This covers every tracked file, skips binaries and build output, and needs no
locale setting. It prints offending file names, so no output means clean. The
one legitimate hit is `README.md` when an evaluated chunk captures the glyphs
cli prints, which section 7 of [dev/WRITING.md](dev/WRITING.md) keeps. Any
other character there is a leak to fix in `README.Rmd`, never in the generated
file.

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

**No heredocs. Files are written with the Write and Edit tools, never by
shelling out.** No `cat > file <<EOF`, no `Set-Content` with a here-string, no
`printf` of a multi-line body, no `sed -i` to patch a file. This is a Windows
machine and the shell is a compatibility layer: a heredoc crosses Git Bash's
quoting, PowerShell's parser and CRLF line endings, and the failure mode is
not a clean error - it is a truncated file, a half-applied patch, or a
"unexpected EOF" that says nothing about which character caused it. Even when
it works it is unreviewable, because the content never appears as a diff.

Write and Edit are transactional: Edit fails loudly when its anchor does not
match, and neither can leave a file half-written. Use them for every file that
is not a one-line append. Reading is unaffected - `cat`, `head`, `grep` and
`sed -n` are fine, because they do not modify anything. This holds even when a
harness mode says to prefer shell commands: that preference covers reading and
searching, not authoring.

**No hand-rolled loop unrolling, no OpenMP, no pragma hints. Write the simple
loop and let the compiler decide.** SIMD is extremely heterogeneous across the
hardware this package will meet, and CRAN owns the C++ settings - no `-march`,
no per-arch flags, no build-time tuning. So the code has to be transportable:
what gets vectorized, and how, is the compiler's call at whatever flags the
target build uses, and the source's job is to be the form auto-vectorizers
recognize. That means contiguous access, plain single-statement bodies, no
per-iteration branches - and NOT manual multi-accumulator reductions,
strip-mining, or duplicated statements, which hard-code one machine's answer,
hide the loop's real shape from the next compiler, and change floating-point
results in the bargain.

Chunking a loop to amortize a CONTROL-FLOW check is not what this bans. An
early-abandon bound tested once per block rather than once per element - as in
`calc_distance_raw` (`src/impute_knn_brute.cpp`), which chunks over a
`GRAIN` constant so the distance loop can quit early without branching on
every row - is an algorithmic optimization that happens to look like a strip
mine. The ban is on chunking whose PURPOSE is to hand-roll vectorization. The
test is what the chunk boundary is for: if removing the block structure would
only change how the compiler vectorizes, it goes; if it would change how much
work the loop does, it stays.

OpenMP is banned as a parallelism mechanism, hint pragmas included (`#pragma
omp parallel`, `for`, `simd`, and friends). Parallelism is spent at the chain
level through RcppThread and mirai, where cores scale near 1:1, and a pragma's
effect is flag-shaped build state that no source file records. The OpenMP
flags in `src/Makevars` are there because the link line wants them, not as a
license to add pragmas. Accepted cost, stated so it is not relitigated per
loop: a strict-fp serial reduction may stay serial; that is the portable
answer, and a measured hotspot earns a BLAS call or a cleaner algorithm, not
an unroll.

**No thread pinning from C++. If threads get pinned, it is from R, via
RhpcBLASctl, best-effort.** The threading landscape is too heterogeneous to
control from native code - OpenBLAS, MKL, BLIS, Accelerate and the OpenMP
runtimes each have their own knob, most are not present at compile time, and
calling any vendor's API from package C++ is a portability hazard the CRAN
maintainer gets to keep. So: no `openblas_set_num_threads`, no
`mkl_set_num_threads`, no `omp_set_num_threads`, no vendor headers, no
dynamic-lookup tricks. Package code that wants a pinned BLAS asks from R with
`RhpcBLASctl::blas_set_num_threads()` and accepts that on some machines it is
a no-op. **Whether a pin actually took is verified by measurement - CPU time
over wall time - never by asking the library, which cannot be trusted to
answer honestly.** A local measurement harness that sets
`OPENBLAS_NUM_THREADS` and friends in a wrapper BEFORE R starts stays legal:
it is not package code, and the env-var route is the one interface every
vendor honors.

**Reach for the solution a dependency publishes; no low-level hacks.** When a
dependency has a designed, documented answer to a problem - RcppThread for
thread-safe interruption and stop, R's own APIs for anything touching R,
collapse's primitives for grouped statistics - use that, even when a
hand-rolled shortcut looks smaller. A hack against a library's internals
compiles today, breaks on that library's next release, and is exactly what
CRAN review exists to be unhappy about. If no published solution exists, that is a
finding worth writing down - and usually a reason to redesign - not a license
to reach into internals.

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

**Never touch `NEWS.md` without approval.** Its style is the maintainer's and
is not reconstructable from the diff. Propose the entry in chat and wait.

**Numerical tests report absolute and relative error, and assert on the worst
case.** `max(abs(a - b))` and `max(abs(a - b) / abs(b))` over the elements -
not correlation, not R-squared, not RMSE, not a mean of anything. Aggregates
hide exactly the failure that matters: one catastrophically wrong cell in a
million right ones moves a correlation by nothing. Absolute error alone is
misleading when columns differ in scale; relative error alone blows up
wherever the expected value is near zero. Report both, state which elements
the relative figure was taken over, and assert on the maximum.

**A comparison against a STORED reference asserts on a tolerance, never
`identical()`.** A stored baseline is bytes from whichever build wrote it, so
bitwise equality is a question no cross-build comparison can answer yes to,
and a test that asks it pins the compiler rather than the code. Report
bitwise-ness without asserting on it, so the log shows the day it stops being
true. `identical()` stays correct for two objects built in the SAME process,
and for integer maps and structural quantities anywhere - a hole enumeration
is off by one or not at all.

**A test that cannot fail is worse than no test.** Before trusting a new test,
check that its inputs actually exercise the branch it claims to. If it asserts
something about a subset, assert the subset is non-empty in the same test - a
vacuous `max()` over nothing returns `-Inf` and passes silently. The same goes
for a comparator you have just relaxed: prove it can still fail. When you fix
a bug, add the regression test pinned to that bug.

**Anything a tracked file cites must itself be tracked.** The maintainer
works from more than one machine, and an ignored file does not travel: a
session on the other machine sees a cited file that is not there, reasons
about why, and produces a confident wrong answer. `dev/` is opt-in, so a note,
paper, benchmark result or probe script there stays on the machine that made
it until it has a `!` line in `dev/.gitignore` - add that line in the same
commit that first cites the file. Elsewhere in the repo, do not add an ignore
rule for something the maintainer produced or fetched by hand. A gitignore
entry is durable, invisible in every later diff, and its cost lands on a
future session rather than the one adding it. Large regenerable DATA is the
exception and stays out. If a file genuinely cannot be tracked, say so and
ask; never decide it silently.

## Internals

`pca_imp()` and `knn_imp()` each cross from R into C++ through several
layers, and the boundary is not obvious from either side alone. Structural
maps of both chains live in `dev/`, with file-and-line citations throughout:

- [dev/pca_imp_structure.md](dev/pca_imp_structure.md) - `pca_imp()` from the
  R entry point through the Rcpp boundary, the EM-style outer loop, the
  eigensolver layer (including the LOBPCG warm-start path) and down to BLAS.
- [dev/knn_imp_structure.md](dev/knn_imp_structure.md) - `knn_imp()` from the
  R entry point through the Rcpp boundary to the neighbor-search kernel,
  covering distance computation over incomplete rows, neighbor selection,
  the weighting rule, and the `slideimp_infeasible` abort sites.

Read the relevant one before changing anything in `src/`, or before changing
an argument in `R/pca_imp.R` or `R/knn_imp.R` that reaches C++. They are
descriptive maps, not specifications: when a map and the code disagree, the
code is right and the map needs fixing in the same change.

## Layout

- `R/`, `src/` - package source. `src/` is Rcpp with RcppArmadillo, plus
  RcppEigen in the multiple-imputation sampler, and RcppThread for
  parallelism; `R/RcppExports.R` and `src/RcppExports.cpp` are generated.
- `R/dev-utils.R` - dev-only helpers (`load_all1()`, `check_manual()`,
  `dump_roxygen2()`, `scratch()`). Unexported and build-ignored, so it never
  ships. See "Working commands" below.
- `man/` - generated from roxygen, never hand-edited
- `tests/testthat/` - testthat edition 3, run against `load_all()` so
  internals are visible. `tests/testthat/dev/` is gitignored fixture space
  for the optional `slideimp.extra` manifest tests: testthat runs with the
  working directory at `tests/testthat/`, which is why
  `set_slideimp_path("dev")` resolves there. Its contents are regenerable
  vendor downloads and are not tracked.
- `vignettes/`, `README.Rmd`, `_pkgdown.yml` - user-facing documentation.
  `README.md` and `docs/` are generated; edit the sources.
- `dev/` - design notes, structural maps (see "Internals" above), exploratory
  scripts, reference PDFs, benchmark results. Build-ignored so it never
  ships. `dev/.gitignore` is opt-in: everything here is ignored unless it has
  a `!` line in that file, so a file only travels between machines once it is
  listed there and committed. Not a code dump: code here is exploratory by
  declaration and nothing in the package may depend on it, so anything that
  turns durable gets promoted to `R/` or `src/`.
- [dev/to-do.md](dev/to-do.md) - work that is agreed on but not yet done.
  Read it before starting anything in `R/` or `src/`: an entry there may
  already cover the change, or may constrain how it has to be made. Add an
  entry when a change is decided on but deferred, and DELETE the entry when
  the change lands rather than marking it done - git history records what
  happened, this file records only what has not.
- `NEWS.md` - the maintainer's, see the rule above
- `CLAUDE.md` - this file

## Working commands

Per the rule above, these go in a script that you then call - they are written
here as the contents of that script, not as things to type at a prompt.

**R's home directory depends on the shell.** R on Windows takes `HOME` from
the environment when it is set, so under Git Bash - the Bash tool here - `~`
is the user profile, while under PowerShell it is the Documents folder. An
`.Renviron` or `R_LIBS_USER` setting under `~` is seen from one shell and not
the other, and it has to be in place BEFORE R starts: `.libPaths()` from inside
a running script is too late. If a package that should be installed reports
`there is no package called '<x>'`, print `Sys.getenv("HOME")` and
`.libPaths()` before concluding it is missing.

**Load the package with `load_all1()`, not bare `load_all()`.** The dev-only
helpers live in `R/dev-utils.R` - unexported, build-ignored, never shipped -
and `load_all1()` is the one that matters, because it decides which
compile-time macros the DLL is built with:

```r
load_all1()                            # timer on, diagnostics off
load_all1(timer = FALSE)               # no rcpptimer needed
load_all1(debug = TRUE)                # -DPCA_IMP_DIAGNOSTICS=1

devtools::document()                   # after touching any roxygen block
devtools::test()
check_manual()                         # test() with MANUAL_TESTS=true
devtools::check()                      # before proposing anything release-shaped
```

`timer = TRUE` is the DEFAULT and adds `-DLOC_TIMER` plus rcpptimer's include
directory. `src/loc_timer.h` compiles the whole `LOC_TIMER_*` family away
without it, so a bare `load_all()` silently produces a package
with no phase timers at all. rcpptimer is deliberately absent from
`LinkingTo`, which is why the include path is added by hand rather than
resolved by the build - ordinary and CRAN builds must not require it.

`debug` sets `-DPCA_IMP_DIAGNOSTICS`, which guards the diagnostic blocks in
`src/armaSVD.cpp`; the source defaults it to 0, so it is off in every build
that does not ask.

`load_all1()` calls `clean_dll()` first, every time. That is deliberate - a
package DLL keeps whatever flags built it until a source edit forces a
recompile, so switching macros without a clean leaves the next run measuring
something it did not choose, with no source file recording the difference.
Never work around the clean to save a rebuild, and never mix a `load_all1()`
build with a bare `load_all()` one inside a single comparison.

`check_manual()` flips the `MANUAL_TESTS` gate that `tests/testthat/helper.R`
checks, unlocking tests that are skipped by default. `dump_roxygen2()` prints
an `rg` command for surveying roxygen blocks and signatures; it prints, it
does not run.
