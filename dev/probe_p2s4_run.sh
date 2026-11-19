#!/usr/bin/env bash
# Step 4 of Phase 2 (dev/eigen_conversion.md): run the Phase 0 capture and
# dev/eigen_p2s4_edges.R on two instrumented clones outside the working tree,
# one of src/lobpcg_warm.h before the step and one after, and compare them.
# Dev only, run by hand, from the package root:
#   bash dev/probe_p2s4_run.sh <scratch dir> [<old commit>]
# <old commit> defaults to HEAD. For the first half of the step it was the
# last commit whose src/lobpcg_warm.h calls Armadillo for its QR, triangular
# solve and dsyevd (eade032 is the one measured); for the second half, which
# replaces the containers and products, it is ff0e394, the first half
# committed; for step 5, which changes the files around it, it is the
# commit of step 4. The "new" clone is HEAD with every file under src/ that
# the working tree changes against HEAD copied over it.
# dev/probe_p2s4_patch.R tells the versions apart by their text.
# Under Git Bash, set R_LIBS_USER before calling this if R's library sits
# under another HOME.
#
# It leaves, in <scratch dir>: p2s4-old/ and p2s4-new/, the clones, each with
# its captures under dev/ref/ (cov, edges, and edges-fresh in the old one);
# p2s4-<mode>-<run>.log, the event logs dev/probe_p2s4_summary.R reads; and
# p2s4-<mode>-<run>.out, what each run printed. The edges-fresh run is made
# only when the old clone's file is the one that kept its dsyevd workspace.
set -uo pipefail
repo="$(pwd)"
scratch="$1"
old_commit="${2:-HEAD}"
mkdir -p "$scratch"
scratch="$(cd "$scratch" && pwd)"

run() {
  # run <mode> <run name> <script> <label>, with the environment as set
  local mode="$1" name="$2" script="$3" label="$4"
  rm -f "$scratch/p2s4-$mode-$name.log"
  SLIDEIMP_P2S4_LOG="$scratch/p2s4-$mode-$name.log" \
    Rscript "$script" "$label" > "$scratch/p2s4-$mode-$name.out" 2>&1
  echo "== $mode $name: exit $?"
  grep "cases, " "$scratch/p2s4-$mode-$name.out"
  grep -c "warning: solve()" "$scratch/p2s4-$mode-$name.out" |
    sed 's/^/   lines of Armadillo solve() warnings on the console: /'
}

for mode in old new; do
  clone="$scratch/p2s4-$mode"
  rm -rf "$clone"
  git clone --quiet "$repo" "$clone"
  if [ "$mode" = old ]; then
    git -C "$clone" checkout --quiet "$old_commit"
  else
    git -C "$repo" diff --name-only --diff-filter=d HEAD -- src |
      while read -r f; do
        cp "$repo/$f" "$clone/$f"
      done
  fi
  cp "$repo/dev/probe_p2s4_patch.R" "$repo/dev/eigen_p2s4_edges.R" "$clone/dev/"
  cd "$clone"
  echo "== $mode clone at $(git rev-parse --short HEAD)"
  Rscript dev/probe_p2s4_patch.R --in-scratch-clone "$mode" || exit 1
  git diff --stat -- src
  run "$mode" capture dev/ref_capture.R cov
  run "$mode" edges dev/eigen_p2s4_edges.R edges
  if [ "$mode" = old ] && grep -q SLIDEIMP_P2S4_FRESH src/lobpcg_warm.h; then
    export SLIDEIMP_P2S4_FRESH=1
    run "$mode" edges-fresh dev/eigen_p2s4_edges.R edges-fresh
    unset SLIDEIMP_P2S4_FRESH
  fi
  cd "$repo"
done

compare() {
  # compare <title> <args to dev/ref_compare.R ...>
  local title="$1"
  shift
  local out="$scratch/p2s4-compare-$title.out"
  Rscript dev/ref_compare.R "$@" > "$out" 2>&1
  echo "== compare $title: exit $?"
  grep "FAIL\|not bitwise" "$out" | grep -v "^  FAIL " | head -40
  tail -1 "$out"
}
old="$scratch/p2s4-old/dev/ref"
new="$scratch/p2s4-new/dev/ref"
# The instrumentation changes no result: each clone's capture against the
# capture of the build it stands for, when the working tree has one.
[ -d dev/ref/p2s3 ] && compare cov-old-p2s3 --tol-abs=0 --tol-rel=0 dev/ref/p2s3 "$old/cov"
[ -d dev/ref/p2s4 ] && compare cov-new-p2s4 --tol-abs=0 --tol-rel=0 dev/ref/p2s4 "$new/cov"
compare edges-old-new-tol0 --tol-abs=0 --tol-rel=0 "$old/edges" "$new/edges"
compare edges-old-new-default "$old/edges" "$new/edges"
[ -d "$old/edges-fresh" ] &&
  compare edges-fresh-new-tol0 --tol-abs=0 --tol-rel=0 "$old/edges-fresh" "$new/edges"
