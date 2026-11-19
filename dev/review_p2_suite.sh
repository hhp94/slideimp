#!/usr/bin/env bash
# Review of Phase 2 (dev/review_p2.md): run dev/review_p2_suite.R on a fresh
# clone of 597d28e in the scratch folder, with a private library for the
# mirai daemons. Dev only, from the package root:
#   bash dev/review_p2_suite.sh <scratch dir>
set -uo pipefail
repo="$(pwd)"
mkdir -p "$1"
scratch="$(cd "$1" && pwd)"
clone="$scratch/suite-src"
lib="$scratch/suite-lib"
rm -rf "$clone" "$lib"
mkdir -p "$lib"
git clone --quiet "$repo" "$clone"
git -C "$clone" checkout --quiet 597d28e
cp "$repo/dev/review_p2_suite.R" "$clone/dev/"
(cd "$clone" && R_LIBS="$(cygpath -m "$lib")" \
  Rscript dev/review_p2_suite.R --in-scratch-clone > "$scratch/suite.log" 2>&1)
echo "== suite on $(git -C "$clone" rev-parse --short HEAD): exit $?"
grep -n -A4 "What a daemon sees" "$scratch/suite.log"
grep -n "^files: " "$scratch/suite.log"
grep -n -A20 "Failing or skipped" "$scratch/suite.log"
