#!/usr/bin/env bash
# Phase 3 of dev/eigen_conversion.md at 4f2b669: every build and run of the
# phase, on clones in a scratch folder outside the repository. The working tree
# is read and never built in. Dev only, run by hand from the package root:
#   bash dev/eigen_p3_gate.sh <scratch dir> <stage>
# Stages, in order:
#   prep     clone HEAD twice: `base` as it is, `p3` with the working tree's
#            changed tracked files copied in. devtools::document() runs on
#            `p3` (CRAN's RcppEigen), and its src/RcppExports.cpp is copied
#            back into the working tree, the one file this stage writes there
#   capture  the Phase 0 capture on both clones, each against CRAN's
#            RcppEigen and the candidate; `p3` against `base` at tolerance 0
#            on each, and `p3` against the stored 76d4c60 at the defaults
#   suite    dev/review_p2_suite.R on `p3`, the whole suite with manual
#            tests, once per RcppEigen, each with its own private library
#            for the mirai daemons
#   check    devtools::check() on `p3`, once per RcppEigen; needs Pandoc on
#            PATH or RSTUDIO_PANDOC set to the folder that holds it
# The candidate is dev/lib/eigen5 (dev/eigen_install_candidate.R installs
# it); it goes first on R_LIBS, so it wins over the CRAN copy in the user
# library. Under Git Bash R does not read Documents/.Renviron, so the caller
# sets R_LIBS_USER to the user library; the script stops when it is unset.
set -uo pipefail
if [ -z "${R_LIBS_USER:-}" ]; then
  echo "set R_LIBS_USER to the user library before calling this"
  exit 2
fi
repo="$(pwd)"
mkdir -p "$1"
scratch="$(cd "$1" && pwd)"
stage="$2"
w() { cygpath -m "$1"; }
eigen5="$(w "$repo/dev/lib/eigen5")"
base="$scratch/base"
p3="$scratch/p3"

# R_LIBS for one side: `cran` is the user library alone, `eigen5` puts the
# candidate first; an extra library, when given, goes in front of both
libs_for() {
  local side="$1" extra="${2:-}"
  local out=""
  [ -n "$extra" ] && out="$extra"
  if [ "$side" = eigen5 ]; then
    out="${out:+$out;}$eigen5"
  fi
  echo "$out"
}

case "$stage" in
prep)
  rm -rf "$base" "$p3"
  git clone --quiet "$repo" "$base"
  git clone --quiet "$repo" "$p3"
  echo "== clones of $(git -C "$base" rev-parse --short HEAD)"
  git -C "$repo" diff --name-only HEAD | while read -r f; do
    cp "$repo/$f" "$p3/$f"
    echo "   copied $f"
  done
  for f in $(git -C "$repo" ls-files --others --exclude-standard); do
    mkdir -p "$(dirname "$p3/$f")"
    cp "$repo/$f" "$p3/$f"
    echo "   copied untracked $f"
  done
  (cd "$p3" && Rscript dev/eigen_p3_check.R document --in-scratch-clone \
    > "$scratch/document.log" 2>&1)
  echo "== document on p3: exit $?"
  tail -5 "$scratch/document.log"
  echo "== what document changed in p3, against the copied tree:"
  git -C "$p3" status --short
  cp "$p3/src/RcppExports.cpp" "$repo/src/RcppExports.cpp"
  ;;
capture)
  for clone in base p3; do
    for side in cran eigen5; do
      label="$clone-$side"
      rm -rf "$scratch/$clone/dev/ref/$label"
      (cd "$scratch/$clone" && R_LIBS="$(libs_for "$side")" \
        Rscript dev/ref_capture.R "$label" > "$scratch/capture-$label.log" 2>&1)
      echo "== capture $label: exit $?"
      grep -n "^RcppEigen \|cases, " "$scratch/capture-$label.log"
    done
  done
  cmp() {
    local title="$1"
    shift
    Rscript dev/ref_compare.R "$@" > "$scratch/compare-$title.out" 2>&1
    echo "== compare $title: exit $?"
    grep -c "NOT bitwise" "$scratch/compare-$title.out" |
      sed 's/^/   lines saying NOT bitwise: /'
    tail -1 "$scratch/compare-$title.out"
  }
  for side in cran eigen5; do
    cmp "$side-tol0" --tol-abs=0 --tol-rel=0 \
      "$(w "$base/dev/ref/base-$side")" "$(w "$p3/dev/ref/p3-$side")"
  done
  cmp 76d4c60-default dev/ref/76d4c60 "$(w "$p3/dev/ref/p3-cran")"
  ;;
suite)
  # dev/review_p2_suite.R takes R_LIBS as one library, so the candidate is
  # copied into the suite's own library rather than put second on R_LIBS
  for side in cran eigen5; do
    lib="$scratch/suite-lib-$side"
    rm -rf "$lib"
    mkdir -p "$lib"
    [ "$side" = eigen5 ] && cp -r "$repo/dev/lib/eigen5/RcppEigen" "$lib/"
    (cd "$p3" && R_LIBS="$(w "$lib")" \
      Rscript dev/review_p2_suite.R --in-scratch-clone \
      > "$scratch/suite-$side.log" 2>&1)
    echo "== suite $side: exit $?"
    grep -n -A4 "What a daemon sees" "$scratch/suite-$side.log"
    grep -n "^files: " "$scratch/suite-$side.log"
    grep -n -A20 "Failing or skipped" "$scratch/suite-$side.log" || true
  done
  ;;
check)
  # the vignettes need Pandoc, which R finds outside RStudio only through
  # RSTUDIO_PANDOC (the folder of RStudio's bundled pandoc.exe) or PATH
  if [ -z "${RSTUDIO_PANDOC:-}" ] && ! command -v pandoc > /dev/null; then
    echo "set RSTUDIO_PANDOC to the folder holding pandoc.exe"
    exit 2
  fi
  for side in cran eigen5; do
    out="$scratch/check-$side"
    rm -rf "$out"
    (cd "$p3" && R_LIBS="$(libs_for "$side")" \
      Rscript dev/eigen_p3_check.R check "$(w "$out")" --in-scratch-clone \
      > "$scratch/check-$side.log" 2>&1)
    echo "== check $side: exit $?"
    grep -n "^RcppEigen " "$scratch/check-$side.log"
    grep -o -m1 'RcppEigen/include' "$out"/slideimp.Rcheck/00install.out |
      sed 's/^/   include seen in 00install.out: /'
    grep -h -o -m1 '[^" ]*RcppEigen/include' "$out"/slideimp.Rcheck/00install.out
    grep -n -A40 "^check: " "$scratch/check-$side.log"
  done
  ;;
*)
  echo "unknown stage $stage"
  exit 2
  ;;
esac
