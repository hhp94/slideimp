#!/usr/bin/env bash
# Phase 4 of dev/eigen_conversion.md at 4f2b669: the gate of one item, on clones
# in a scratch folder outside the repository. The working tree is read and never
# built in. Dev only, run by hand from the package root:
#   bash dev/eigen_p4_gate.sh <scratch dir> <stage>
# Stages, in order:
#   prep     clone HEAD twice: `base` as it is, `p4` with the working tree's
#            changed and untracked files copied in
#   capture  the Phase 0 capture on both clones, each against CRAN's
#            RcppEigen and the candidate; `p4` against `base` at tolerance 0
#            and at the defaults on each. Both captures are of one run on
#            one machine, so tolerance 0 is a fair question where an item is
#            expected to be bitwise
#   diag     dev/eigen_p2s6_diag.R, the diagnostics build, on both clones
#            against CRAN's RcppEigen, compared the same way: the one build
#            that compiles and calls si_la::svd_values
#   edges    two more clones, `edge-base` and `edge-p4`, instrumented by
#            dev/probe_p2s4_patch.R, and dev/eigen_p2s4_edges.R on each with
#            the event log: LOBPCG past what the capture asks of it, and
#            each of its factorizations made to fail. Compared the same way,
#            then dev/probe_p2s4_summary.R on the two logs
#   time     the timing driver on installs of `base` and `p4`, PCA workloads
#            only; P4_TIME_REPS sets the repetitions (default 5)
#   suite    dev/review_p2_suite.R on `p4`, the whole suite with manual
#            tests, once per RcppEigen, each with its own private library
#            for the mirai daemons
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
p4="$scratch/p4"

libs_for() {
  if [ "$1" = eigen5 ]; then
    echo "$eigen5"
  fi
}

cmp() {
  # cmp <title> <args of dev/ref_compare.R>
  local title="$1"
  shift
  Rscript dev/ref_compare.R "$@" > "$scratch/compare-$title.out" 2>&1
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
  git -C "$repo" diff --name-only HEAD | while read -r f; do
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
        Rscript dev/ref_capture.R "$label" > "$scratch/capture-$label.log" 2>&1)
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
diag)
  for clone in base p4; do
    label="$clone-diag"
    rm -rf "$scratch/$clone/dev/ref/$label"
    (cd "$scratch/$clone" && R_LIBS="" \
      Rscript dev/eigen_p2s6_diag.R "$label" > "$scratch/diag-$clone.log" 2>&1)
    echo "== diag $clone: exit $?"
    grep -n "cases, " "$scratch/diag-$clone.log"
  done
  cmp diag-tol0 --tol-abs=0 --tol-rel=0 \
    "$(w "$base/dev/ref/base-diag")" "$(w "$p4/dev/ref/p4-diag")"
  cmp diag-default \
    "$(w "$base/dev/ref/base-diag")" "$(w "$p4/dev/ref/p4-diag")"
  ;;
edges)
  for clone in base p4; do
    dir="$scratch/edge-$clone"
    rm -rf "$dir" "$scratch/edge-$clone.log"
    cp -r "$scratch/$clone" "$dir"
    rm -rf "$dir/dev/ref/edges"
    cp "$repo/dev/probe_p2s4_patch.R" "$repo/dev/eigen_p2s4_edges.R" "$dir/dev/"
    (cd "$dir" && Rscript dev/probe_p2s4_patch.R --in-scratch-clone new &&
      SLIDEIMP_P2S4_LOG="$(w "$scratch/edge-$clone.log")" \
        Rscript dev/eigen_p2s4_edges.R edges > "$scratch/edges-$clone.out" 2>&1)
    echo "== edges $clone: exit $?"
    grep "cases, " "$scratch/edges-$clone.out"
  done
  cmp edges-tol0 --tol-abs=0 --tol-rel=0 \
    "$(w "$scratch/edge-base/dev/ref/edges")" "$(w "$scratch/edge-p4/dev/ref/edges")"
  cmp edges-default \
    "$(w "$scratch/edge-base/dev/ref/edges")" "$(w "$scratch/edge-p4/dev/ref/edges")"
  Rscript dev/probe_p2s4_summary.R \
    "$(w "$scratch/edge-base.log")" "$(w "$scratch/edge-p4.log")" \
    "$(w "$scratch/edge-base/dev/ref/edges")" "$(w "$scratch/edge-p4/dev/ref/edges")" \
    > "$scratch/summary-edges.out" 2>&1
  echo "== summary edges: exit $?"
  ;;
time)
  # for an item that adds work: both clones installed by
  # dev/eigen_time_install.R, and the four PCA workloads of the timing
  # driver, `base` as a and `p4` as b
  for clone in base p4; do
    rm -rf "$scratch/time-lib-$clone"
    Rscript dev/eigen_time_install.R "$(w "$scratch/$clone")" \
      "$(w "$scratch/time-lib-$clone")" > "$scratch/time-install-$clone.log" 2>&1
    echo "== install $clone: exit $?"
  done
  Rscript dev/eigen_time_builds.R \
    --lib-a="$(w "$scratch/time-lib-base")" --lib-b="$(w "$scratch/time-lib-p4")" \
    --reps="${P4_TIME_REPS:-5}" \
    --workloads=chr21_pca_exact,chr21_pca_lobpcg,sim_pca_exact,sim_pca_lobpcg \
    --out="$(w "$scratch/time.csv")" > "$scratch/time.out" 2>&1
  echo "== time: exit $?"
  tail -9 "$scratch/time.out"
  ;;
suite)
  # dev/review_p2_suite.R takes R_LIBS as one library, so the candidate is
  # copied into the suite's own library rather than put second on R_LIBS
  for side in cran eigen5; do
    lib="$scratch/suite-lib-$side"
    rm -rf "$lib"
    mkdir -p "$lib"
    [ "$side" = eigen5 ] && cp -r "$repo/dev/lib/eigen5/RcppEigen" "$lib/"
    (cd "$p4" && R_LIBS="$(w "$lib")" \
      Rscript dev/review_p2_suite.R --in-scratch-clone \
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
