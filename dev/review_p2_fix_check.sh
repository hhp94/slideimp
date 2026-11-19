#!/usr/bin/env bash
# Fixes after the review of Phase 2 (dev/review_p2.md, findings 1 and 3):
# run the two tests added to tests/testthat/test-pca_imp.R against the
# builds dev/review_p2_run.sh and dev/review_p2_tag.sh installed in the
# scratch folder, each test in its own process, and print how each ended.
# Then run the two PCA test files whole, against load_all1(timer = FALSE),
# on a fresh clone of HEAD that takes the working tree's tests/. Stage
# `suite`, run only when named, is the whole suite on such a clone. Nothing
# is built in the working tree. Dev only, from the package root:
#   bash dev/review_p2_fix_check.sh <scratch dir> [installed|clone|suite]
set -uo pipefail
repo="$(pwd)"
scratch="$(cd "$1" && pwd)"
stage="${2:-all}"
w() { cygpath -m "$1"; }

if [ "$stage" = all ] || [ "$stage" = installed ]; then
  for ref in v1.2.0 old 905e1ac eade032 new; do
    for which in rank overflow; do
      out="$scratch/fix-$ref-$which.out"
      Rscript dev/review_p2_fix_check.R "$(w "$scratch/lib-$ref")" "$which" \
        > "$out" 2>&1
      echo "== $ref $which: exit $?"
      grep -a "warning\|^exact\|^lobpcg\|: tests \|Segmentation\|Failure\|Error" "$out"
    done
  done
fi

if [ "$stage" = all ] || [ "$stage" = clone ]; then
  clone="$scratch/fix-src"
  rm -rf "$clone"
  git clone --quiet "$repo" "$clone"
  cp "$repo"/tests/testthat/*.R "$clone/tests/testthat/"
  cp "$repo/dev/review_p2_fix_suite.R" "$clone/dev/"
  (cd "$clone" && Rscript dev/review_p2_fix_suite.R --in-scratch-clone \
    > "$scratch/fix-suite.log" 2>&1)
  echo "== PCA test files on a clone of $(git -C "$clone" rev-parse --short HEAD): exit $?"
  git -C "$clone" status --short
  grep -n "^files: " "$scratch/fix-suite.log"
  grep -n -A20 "Failing or skipped\|with a warning" "$scratch/fix-suite.log" || true
fi

# the whole suite, as dev/review_p2_suite.sh runs it, on the same kind of
# clone: HEAD with the working tree's tests, a private library for the
# mirai daemons
if [ "$stage" = suite ]; then
  clone="$scratch/fix-suite-src"
  lib="$scratch/fix-suite-lib"
  rm -rf "$clone" "$lib"
  mkdir -p "$lib"
  git clone --quiet "$repo" "$clone"
  cp "$repo"/tests/testthat/*.R "$clone/tests/testthat/"
  cp "$repo/dev/review_p2_suite.R" "$clone/dev/"
  (cd "$clone" && R_LIBS="$(w "$lib")" \
    Rscript dev/review_p2_suite.R --in-scratch-clone \
    > "$scratch/fix-suite-all.log" 2>&1)
  echo "== suite on a clone of $(git -C "$clone" rev-parse --short HEAD): exit $?"
  grep -n -A4 "What a daemon sees" "$scratch/fix-suite-all.log"
  grep -n "^files: " "$scratch/fix-suite-all.log"
  grep -n -A20 "Failing or skipped" "$scratch/fix-suite-all.log" || true
fi
