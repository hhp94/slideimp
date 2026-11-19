#!/usr/bin/env bash
# Review of Phase 2 of dev/eigen_conversion.md (dev/review_p2.md): install
# one more commit or tag into its own library in the scratch folder and run
# groups of dev/review_p2_cases.R on it, to compare with the two builds of
# dev/review_p2_run.sh. Used for v1.2.0, the release, and for 905e1ac and
# eade032, the commits either side of step 3. Dev only, from the package
# root:
#   bash dev/review_p2_tag.sh <scratch dir> <commit or tag> "<groups>"
set -uo pipefail
repo="$(pwd)"
mkdir -p "$1"
scratch="$(cd "$1" && pwd)"
ref="$2"
groups="$3"
w() { cygpath -m "$1"; }

clone="$scratch/src-$ref"
if [ ! -d "$scratch/lib-$ref" ]; then
  rm -rf "$clone"
  git clone --quiet "$repo" "$clone"
  git -C "$clone" checkout --quiet "$ref"
  Rscript dev/eigen_time_install.R "$(w "$clone")" "$(w "$scratch/lib-$ref")" \
    > "$scratch/install-$ref.log" 2>&1
  echo "== install $ref: exit $?"
  tail -1 "$scratch/install-$ref.log"
fi
for group in $groups; do
  rm -f "$scratch/cases-$ref/${group}__"*.rds
  Rscript dev/review_p2_cases.R "$(w "$scratch/lib-$ref")" \
    "$(w "$scratch/cases-$ref")" "$group" > "$scratch/cases-$ref-$group.log" 2>&1
  echo "== cases $ref $group: exit $?"
  grep "NO RESULT\|ERROR" "$scratch/cases-$ref-$group.log"
  tail -1 "$scratch/cases-$ref-$group.log"
done
