#!/usr/bin/env bash
# Review of Phase 2 of dev/eigen_conversion.md (dev/review_p2.md): every
# build and run of the review, each stage on clones in a scratch folder
# outside the repository. The working tree is read and never built in. Dev
# only, run by hand from the package root:
#   bash dev/review_p2_run.sh <scratch dir> <stage> [...]
# Stages, in the order the note uses them:
#   build    clone 59920f4 (old) and 597d28e (new) and install each into its
#            own library with dev/eigen_time_install.R
#   cases    dev/review_p2_cases.R on both libraries, then
#            dev/review_p2_compare.R, a group at a time; arguments: the
#            groups, quoted as one word (default: all seven)
#   capture  the Phase 0 capture on a fresh clone of 597d28e, compared with
#            the stored captures 59920f4, p2s6 and 76d4c60
#   edges    dev/probe_p2s4_run.sh against 597d28e; arguments: the old
#            commit (default 59920f4, whose instrumented clone dies on a
#            rank-deficient case; eade032, step 3, is the first that does not)
#   diag     dev/probe_p2s6_run.sh, 59920f4 against 597d28e
#   time     the timing driver, old as a and new as b; arguments: reps
# Under Git Bash, set R_LIBS_USER before calling this if R's library sits
# under another HOME.
set -uo pipefail
repo="$(pwd)"
mkdir -p "$1"
scratch="$(cd "$1" && pwd)"
stage="$2"
old_commit=59920f4
new_commit=597d28e
w() { cygpath -m "$1"; }

case "$stage" in
build)
  for side in old new; do
    commit="$old_commit"
    [ "$side" = new ] && commit="$new_commit"
    clone="$scratch/src-$side"
    rm -rf "$clone" "$scratch/lib-$side"
    git clone --quiet "$repo" "$clone"
    git -C "$clone" checkout --quiet "$commit"
    Rscript dev/eigen_time_install.R "$(w "$clone")" "$(w "$scratch/lib-$side")" \
      > "$scratch/install-$side.log" 2>&1
    echo "== install $side ($commit): exit $?"
    grep -c "Wignored-attributes" "$scratch/install-$side.log" |
      sed 's/^/   -Wignored-attributes lines: /'
    grep -n -i "warning\|error" "$scratch/install-$side.log" |
      grep -v "Wignored-attributes\|In file included\|required from" | head -20
    tail -1 "$scratch/install-$side.log"
  done
  ;;
cases)
  groups="${3:-tiny tinysweep general big validation overflow dup}"
  for side in old new; do
    lib="$(w "$scratch/lib-$side")"
    out="$(w "$scratch/cases-$side")"
    for group in $groups; do
      rm -f "$scratch/cases-$side/${group}__"*.rds
      Rscript dev/review_p2_cases.R "$lib" "$out" "$group" \
        > "$scratch/cases-$side-$group.log" 2>&1
      echo "== cases $side $group: exit $?"
      tail -1 "$scratch/cases-$side-$group.log"
    done
  done
  echo "== R processes running now (only another stage's, if one is running):"
  tasklist 2>/dev/null | grep -i "rterm\|rscript" || echo "   none"
  for group in $groups; do
    Rscript dev/review_p2_compare.R "$(w "$scratch/cases-old")" \
      "$(w "$scratch/cases-new")" "^${group}__" \
      > "$scratch/cases-compare-$group.out" 2>&1
    echo "== compare $group: exit $?"
    tail -2 "$scratch/cases-compare-$group.out"
  done
  ;;
capture)
  clone="$scratch/cap-new"
  rm -rf "$clone"
  git clone --quiet "$repo" "$clone"
  git -C "$clone" checkout --quiet "$new_commit"
  (cd "$clone" && Rscript dev/ref_capture.R review > "$scratch/capture.log" 2>&1)
  echo "== capture on $new_commit: exit $?"
  grep -n "^capture \|^RcppEigen \|cases, " "$scratch/capture.log"
  cmp() {
    local title="$1"
    shift
    Rscript dev/ref_compare.R "$@" > "$scratch/capture-compare-$title.out" 2>&1
    echo "== compare $title: exit $?"
    grep -c "NOT bitwise" "$scratch/capture-compare-$title.out" |
      sed 's/^/   lines saying NOT bitwise: /'
    tail -1 "$scratch/capture-compare-$title.out"
  }
  cmp 59920f4-tol0 --tol-abs=0 --tol-rel=0 dev/ref/59920f4 "$(w "$clone/dev/ref/review")"
  cmp p2s6-tol0 --tol-abs=0 --tol-rel=0 dev/ref/p2s6 "$(w "$clone/dev/ref/review")"
  cmp 76d4c60-default dev/ref/76d4c60 "$(w "$clone/dev/ref/review")"
  echo "== n_iter, stored 59920f4 against the fresh capture"
  Rscript dev/probe_p2s6_niter.R dev/ref/59920f4 "$(w "$clone/dev/ref/review")"
  ;;
edges)
  from="${3:-$old_commit}"
  bash dev/probe_p2s4_run.sh "$scratch/edges-$from" "$from" \
    > "$scratch/edges-$from-run.out" 2>&1
  echo "== edges probe, $from against $new_commit: exit $?"
  grep "^== \|cases, \|lines of\|rror" "$scratch/edges-$from-run.out"
  ;;
diag)
  bash dev/probe_p2s6_run.sh "$scratch/diag" "$old_commit" \
    > "$scratch/diag-run.out" 2>&1
  echo "== diag probe: exit $?"
  grep -v "^warning" "$scratch/diag-run.out" | tail -40
  ;;
time)
  reps="${3:-5}"
  Rscript dev/eigen_time_builds.R \
    --lib-a="$(w "$scratch/lib-old")" --lib-b="$(w "$scratch/lib-new")" \
    --reps="$reps" \
    --workloads=chr21_pca_exact,chr21_pca_lobpcg,sim_pca_exact,sim_pca_lobpcg \
    --out="$(w "$scratch/time-$reps.csv")" > "$scratch/time-$reps.log" 2>&1
  echo "== driver: exit $?"
  tail -30 "$scratch/time-$reps.log"
  ;;
*)
  echo "unknown stage $stage"
  exit 2
  ;;
esac
