# CI for a CRAN release

Which CRAN check each platform stands for, and the run tiers of a release.
A local run is an R CMD check in one of rhub's Linux images, in a VM on one
machine, started by a script in local/ there (the local script, below). How
to start one is in local/, a folder on that machine only.

## Where CRAN reports

- The package page, https://cran.r-project.org/web/checks/check_results_slideimp.html:
  one row per platform, and an "Additional issues" line under the table when
  one of CRAN's extra checks has flagged something.
- One platform's full log:
  https://www.r-project.org/nosvn/R.check/<platform>/slideimp-00check.html,
  `<platform>` as the package page names it.
- The kinds of extra check: https://cran.r-project.org/web/checks/check_issue_kinds.html.
  Each links to its README, https://www.stats.ox.ac.uk/pub/bdr/<kind>/README.txt
  for most, github.com/kalibera/cran-checks for rchk, rcnst and rlibro. A
  kind's results are in its README's folder, one file per package.

## Where each check runs

"VM image" is an rhub container, run by the local script. "rhub <name>" is
a platform of `rhub::rhub_check()`, which runs on GitHub's runners. An empty
Elsewhere cell means the VM image is the only place.

CRAN's platforms:

| CRAN platform                     | VM image        | Elsewhere |
|-----------------------------------|-----------------|-----------|
| r-release-linux-x86_64            | ubuntu-release  | |
| r-patched-linux-x86_64            | ubuntu-next     | |
| r-devel-linux-x86_64-debian-clang | ubuntu-clang, or debian-clang | |
| r-devel-linux-x86_64-debian-gcc   | ubuntu-gcc16    | |
| r-devel-linux-x86_64-fedora-clang | none exact; nearest clang23 | |
| r-devel-linux-x86_64-fedora-gcc   | none exact; nearest gcc16 | |
| r-devel-windows-x86_64            | none            | `devtools::check_win_devel()` |
| r-release-windows-x86_64          | none            | `devtools::check()` on Windows; `devtools::check_win_release()` |
| r-oldrel-windows-x86_64           | none            | `devtools::check_win_oldrelease()` |
| r-release-macos-arm64             | none            | `devtools::check_mac_release()`; rhub macos-arm64 |
| r-oldrel-macos-arm64              | none            | none |
| r-release-macos-x86_64            | none            | rhub macos |
| r-oldrel-macos-x86_64             | none            | none |

CRAN's extra checks:

| Kind                     | VM image       | Elsewhere |
|--------------------------|----------------|-----------|
| clang-ASAN, clang-UBSAN  | clang-asan, clang-ubsan | |
| gcc-ASAN, gcc-UBSAN      | gcc-asan (one image for both) | |
| valgrind                 | valgrind       | |
| rchk                     | rchk           | |
| noLD                     | nold           | |
| noSuggests               | nosuggests     | |
| MKL, BLIS, ATLAS         | mkl, blis, atlas | |
| OpenBLAS                 | ubuntu-release | |
| LTO                      | lto            | |
| clang23                  | clang23        | |
| gcc                      | gcc16          | |
| vnu                      | vnu            | |
| M1mac                    | none           | `devtools::check_mac_release()`; rhub macos-arm64 |
| M1-SAN                   | none           | rhub m1-san |

Extra checks with no run anywhere: rcnst, musl, linux-arm64 and rlibro.
rcnst reproduces on any VM image with `R_COMPILE_PKGS=1 R_JIT_STRATEGY=4
R_CHECK_CONSTANTS=5`, which the local script does not set.

Extra checks that need no run: donttest, noOMP, BLAS, C23, Strict, noRemap,
Intel and 0len (dev/DECISIONS.md, "CI").

## Run tiers

A run tier is a set of platforms to check on. It is not a test tier, the
CRAN, local or manual set of tests that `.claude/CLAUDE.md` defines; a local
run picks its test tier by flag: none for the local tier (NOT_CRAN=true),
--cran for the CRAN tier, --manual for the manual tier.

Run tiers 1 and 2 run on the VM, on committed HEAD, and need no push; tier 2
includes tier 1. Tier 3 needs no push either; tier 4 does.

1. Minimal: after a change to `R/` or `src/`, before merging into the dev
   branch. ubuntu-release, clang-asan and nosuggests, in the local test tier.
2. Full: before a release, and after any change to the PCA kernel.
   - First rebuild the saved images from rhub's latest (the local script's
     --fresh).
   - Local test tier on these groups, run in this order:
     - CRAN's Linux platforms: ubuntu-release, ubuntu-next, ubuntu-clang,
       ubuntu-gcc16, gcc16, clang23.
     - Memory and undefined behavior: clang-asan, clang-ubsan, gcc-asan,
       valgrind.
     - Static analysis: rchk.
     - Numerics: nold, atlas, mkl, blis.
     - Packaging: nosuggests, lto, vnu, and an r-ver image (rocker's image
       of an old R) at DESCRIPTION's R floor.
   - CRAN test tier (--cran) on ubuntu-release, nold, atlas, mkl and blis.
   - Manual test tier (--manual) on atlas, mkl and blis.
3. CRAN's builders: `devtools::check_win_devel()`,
   `devtools::check_mac_release()`, and `devtools::check(remote = TRUE,
   manual = TRUE)` on Windows. `devtools::check_win_release()` and
   `devtools::check_win_oldrelease()` when a change touches Windows-only
   code.
4. rhub, after the push: `rhub::rhub_check()` on macos-arm64, macos and
   m1-san, through `.github/workflows/rhub.yaml`. It checks what was pushed
   to GitHub.

## Release issue checklist

Add these steps to the Prepare for release section of the checklist
use_release_issue() writes. Here `*` marks usethis's own lines and `+` the
added ones; in the issue all are `*`:

```
* [ ] `devtools::check(remote = TRUE, manual = TRUE)`
* [ ] `devtools::check_win_devel()`
+ [ ] Local CI, full tier: needs no push, and its results feed cran-comments.md
+ [ ] `devtools::check_mac_release()`
* [ ] Update `cran-comments.md`
* [ ] `git push`
+ [ ] `rhub::rhub_check()` on macos-arm64, macos and m1-san: needs the push
```

The issue is public: write it without naming the VM's machine, its address
or anything in local/.
