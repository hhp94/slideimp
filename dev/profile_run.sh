#!/usr/bin/env bash
# The phase profile of the kernels, on a clone in a scratch folder outside the
# repository. Dev only, run by hand from the package root:
#   bash dev/profile_run.sh <scratch dir> <stage> [n_calls]
# Stages, in order:
#   prep     clone HEAD into <scratch>/prof and copy the working tree's
#            changed and untracked files under R/ and src/ over it
#   install  dev/profile_install.R: a release build with the phase timers
#            into dev/lib/prof
#   run      dev/profile_phases.R against that library
# RSCRIPT names the Rscript to use when none is on the path. Under Git Bash R
# does not read Documents/.Renviron, so the caller sets R_LIBS_USER to the
# user library; the script stops when it is unset.
set -uo pipefail
if [ -z "${R_LIBS_USER:-}" ]; then
  echo "set R_LIBS_USER to the user library before calling this"
  exit 2
fi
rscript="${RSCRIPT:-Rscript}"
repo="$(pwd)"
mkdir -p "$1"
scratch="$(cd "$1" && pwd)"
stage="$2"
n_calls="${3:-5}"
w() { cygpath -m "$1"; }
clone="$scratch/prof"
lib="$repo/dev/lib/prof"

case "$stage" in
prep)
  rm -rf "$clone"
  git clone -q "$repo" "$clone" || exit 1
  git -C "$repo" status --porcelain -- R src | while read -r st path; do
    echo "copying $st $path"
    cp "$repo/$path" "$clone/$path" || exit 1
  done
  git -C "$clone" status --short
  ;;
install)
  mkdir -p "$lib"
  "$rscript" dev/profile_install.R "$(w "$clone")" "$(w "$lib")"
  ;;
run)
  R_LIBS="$(w "$lib")" "$rscript" dev/profile_phases.R "$(w "$lib")" \
    "$(w "$repo/dev/GSE286313.EPICv1_chr21.qs2")" "$n_calls"
  ;;
*)
  echo "unknown stage: $stage"
  exit 2
  ;;
esac
