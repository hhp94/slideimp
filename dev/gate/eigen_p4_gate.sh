#!/usr/bin/env bash
# Phase 4 of dev/eigen_conversion.md at 4f2b669: the gate of one item, on clones
# in a scratch folder outside the repository. The working tree is read and never
# built in. Dev only, run by hand from the package root:
#   bash dev/gate/eigen_p4_gate.sh <scratch dir> <stage>
# Stages, in order:
#   prep     clone HEAD twice: `base` as it is, `p4` with the working tree's
#            changed and untracked files copied in
#   capture  the Phase 0 capture on both clones, each against CRAN's
#            RcppEigen and the candidate; `p4` against `base` at tolerance 0
#            and at the defaults on each. Both captures are of one run on
#            one machine, so tolerance 0 is a fair question where an item is
#            expected to be bitwise
#   time     the timing driver on installs of `base` and `p4`, PCA workloads
#            only; P4_TIME_REPS sets the repetitions (default 5)
#   suite    dev/gate/review_p2_suite.R on `p4`, the whole suite with manual
#            tests, once per RcppEigen, each with its own private library
#            for the mirai daemons
# Two stages are gone: `diag` and `edges` ran the kernel's edge inputs
# (dev/eigen_p2s6_diag.R and dev/eigen_p2s4_edges.R at d5a34e5) on two
# builds. Those inputs are tests now, tests/testthat/test-pca_kernel-manual.R,
# with the exact solver as the oracle, so `suite` runs them under both
# RcppEigen versions. What `edges` also did, make a factorization fail on an
# instrumented build, has no test.
# The candidate is dev/lib/eigen5 (dev/gate/eigen_install_candidate.R installs
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
p4="$scratch/p4"

libs_for() {
  if [ "$1" = eigen5 ]; then
    echo "$eigen5"
  fi
}

cmp() {
  # cmp <title> <args of dev/gate/ref_compare.R>
  local title="$1"
  shift
  Rscript dev/gate/ref_compare.R "$@" > "$scratch/compare-$title.out" 2>&1
  echo "== compare $title: exit $?"
  tail -1 "$scratch/compare-$title.out"
}

case "$stage" in
prep)
  rm -rf "$base" "$p4"
  git clone --quiet "$repo" "$base"
  git clone --quiet "$repo" "$p4"
  echo "== clones of $(git -C "$base" rev-parse --short HEAD)"
  # a file the working tree deleted is listed too: it is removed from `p4`,
  # not copied
  git -C "$repo" diff --no-renames --name-only HEAD | while read -r f; do
    if [ ! -e "$repo/$f" ]; then
      rm -f "$p4/$f"
      echo "   removed $f"
    elif mkdir -p "$(dirname "$p4/$f")" && cp "$repo/$f" "$p4/$f"; then
      echo "   copied $f"
    else
      echo "   FAILED to copy $f"
      exit 1
    fi
  done || exit 1
  for f in $(git -C "$repo" ls-files --others --exclude-standard); do
    mkdir -p "$(dirname "$p4/$f")"
    cp "$repo/$f" "$p4/$f"
    echo "   copied untracked $f"
  done
  git -C "$p4" status --short
  ;;
capture)
  for clone in base p4; do
    for side in cran eigen5; do
      label="$clone-$side"
      rm -rf "$scratch/$clone/dev/ref/$label"
      (cd "$scratch/$clone" && R_LIBS="$(libs_for "$side")" \
        Rscript dev/gate/ref_capture.R "$label" > "$scratch/capture-$label.log" 2>&1)
      echo "== capture $label: exit $?"
      grep -n "^RcppEigen \|cases, " "$scratch/capture-$label.log"
    done
  done
  for side in cran eigen5; do
    cmp "$side-tol0" --tol-abs=0 --tol-rel=0 \
      "$(w "$base/dev/ref/base-$side")" "$(w "$p4/dev/ref/p4-$side")"
    cmp "$side-default" \
      "$(w "$base/dev/ref/base-$side")" "$(w "$p4/dev/ref/p4-$side")"
  done
  ;;
time)
  # for an item that adds work: both clones installed by
  # dev/gate/eigen_time_install.R, and the four PCA workloads of the timing
  # driver, `base` as a and `p4` as b
  for clone in base p4; do
    rm -rf "$scratch/time-lib-$clone"
    Rscript dev/gate/eigen_time_install.R "$(w "$scratch/$clone")" \
      "$(w "$scratch/time-lib-$clone")" > "$scratch/time-install-$clone.log" 2>&1
    echo "== install $clone: exit $?"
  done
  Rscript dev/gate/eigen_time_builds.R \
    --lib-a="$(w "$scratch/time-lib-base")" --lib-b="$(w "$scratch/time-lib-p4")" \
    --reps="${P4_TIME_REPS:-5}" \
    --workloads=chr21_pca_exact,chr21_pca_lobpcg,sim_pca_exact,sim_pca_lobpcg \
    --out="$(w "$scratch/time.csv")" > "$scratch/time.out" 2>&1
  echo "== time: exit $?"
  tail -9 "$scratch/time.out"
  ;;
suite)
  # dev/gate/review_p2_suite.R takes R_LIBS as one library, so the candidate is
  # copied into the suite's own library rather than put second on R_LIBS
  for side in cran eigen5; do
    lib="$scratch/suite-lib-$side"
    rm -rf "$lib"
    mkdir -p "$lib"
    [ "$side" = eigen5 ] && cp -r "$repo/dev/lib/eigen5/RcppEigen" "$lib/"
    (cd "$p4" && R_LIBS="$(w "$lib")" \
      Rscript dev/gate/review_p2_suite.R --in-scratch-clone \
      > "$scratch/suite-$side.log" 2>&1)
    echo "== suite $side: exit $?"
    grep -n -A4 "What a daemon sees" "$scratch/suite-$side.log"
    grep -n "^files: " "$scratch/suite-$side.log"
    grep -n -A20 "Failing or skipped" "$scratch/suite-$side.log" || true
  done
  ;;
*)
  echo "unknown stage $stage"
  exit 2
  ;;
esac
