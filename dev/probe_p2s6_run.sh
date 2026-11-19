#!/usr/bin/env bash
# Step 6 of Phase 2 (dev/eigen_conversion.md): run dev/eigen_p2s6_diag.R on
# two instrumented clones outside the working tree, one of the kernel before
# the step and one after, and compare them. Dev only, run by hand, from the
# package root:
#   bash dev/probe_p2s6_run.sh <scratch dir> [<old commit>]
# <old commit> defaults to HEAD; the one measured is 8c15416, step 5. The
# "new" clone is HEAD with every file under src/ and R/ that the working tree
# changes against HEAD copied over it. dev/probe_p2s6_patch.R tells the two
# versions of src/armaSVD.cpp apart by their text.
#
# It leaves, in <scratch dir>: p2s6-old/ and p2s6-new/, the clones, each with
# its capture under dev/ref/diag, and p2s6-<mode>.out, what each run printed.
set -uo pipefail
repo="$(pwd)"
scratch="$1"
old_commit="${2:-HEAD}"
mkdir -p "$scratch"
scratch="$(cd "$scratch" && pwd)"

for mode in old new; do
  clone="$scratch/p2s6-$mode"
  rm -rf "$clone"
  git clone --quiet "$repo" "$clone"
  if [ "$mode" = old ]; then
    git -C "$clone" checkout --quiet "$old_commit"
  else
    git -C "$repo" diff --name-only --diff-filter=d HEAD -- src R |
      while read -r f; do
        cp "$repo/$f" "$clone/$f"
      done
  fi
  cp "$repo/dev/probe_p2s6_patch.R" "$repo/dev/eigen_p2s6_diag.R" "$clone/dev/"
  cd "$clone"
  echo "== $mode clone at $(git rev-parse --short HEAD)"
  Rscript dev/probe_p2s6_patch.R --in-scratch-clone || exit 1
  Rscript dev/eigen_p2s6_diag.R diag > "$scratch/p2s6-$mode.out" 2>&1
  echo "== $mode: exit $?"
  grep "cases, \|ERROR" "$scratch/p2s6-$mode.out"
  cd "$repo"
done

out="$scratch/p2s6-compare.out"
Rscript dev/ref_compare.R --tol-abs=0 --tol-rel=0 \
  "$scratch/p2s6-old/dev/ref/diag" "$scratch/p2s6-new/dev/ref/diag" > "$out" 2>&1
echo "== compare at tolerance 0: exit $?"
grep "FAIL\|not bitwise" "$out" | head -40
tail -1 "$out"
