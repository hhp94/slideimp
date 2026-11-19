#!/usr/bin/env bash
# Audit of dev/review_eigen_wrappers.md, finding 1 (dev/eigen_conversion.md,
# "The review of the wrapper layer"): does the chain's own dsyevr call write
# past its eigenvalue buffer when pca_imp() reaches it? Clones a commit
# outside the working tree, instruments the clone's src/eig_sym_sel.h, and
# runs pca_imp() on it. Dev only, run by hand, from the package root:
#   bash dev/audit_eigen_wrappers_chain.sh <scratch dir> [<commit>]
# <commit> defaults to HEAD; it must be one in which src/eig_sym_sel.h still
# calls dsyevr itself (79cce0f is the one measured). Under Git Bash, set
# R_LIBS_USER before calling this if R's library sits under another HOME.
set -euo pipefail
repo="$(pwd)"
scratch="$1"
commit="${2:-HEAD}"
clone="$scratch/audit-chain"

rm -rf "$clone"
git clone --quiet "$repo" "$clone"
git -C "$clone" checkout --quiet "$commit"
cp "$repo/dev/audit_eigen_wrappers_patch.R" "$repo/dev/audit_eigen_wrappers_chain.R" "$clone/dev/"
cd "$clone"
echo "== clone at $(git rev-parse --short HEAD)"
Rscript dev/audit_eigen_wrappers_patch.R --in-scratch-clone
git diff --stat -- src
if ! Rscript dev/audit_eigen_wrappers_chain.R "$scratch/audit-chain-calls.log" "$repo/dev/GSE286313.EPICv1_chr21.qs2" > "$scratch/audit-chain.log" 2>&1; then
  tail -40 "$scratch/audit-chain.log"
  exit 1
fi
sed -n '/^R: /,$p' "$scratch/audit-chain.log"
