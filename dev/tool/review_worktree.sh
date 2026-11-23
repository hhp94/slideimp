#!/usr/bin/env bash
# Set up a throwaway worktree for a code review of a commit range.
#
#   bash dev/tool/review_worktree.sh <start> <end> <path>
#
# Run from the slideimp root. The worktree's files and index are at <end>
# and its HEAD is at <start>, so `git diff HEAD` there is the whole range,
# added and deleted files included, and /code-review with no target reviews
# it. HEAD is detached: no branch is made. Remove it afterwards with
#
#   git worktree remove --force <path>
#
# The worktree's .claude/CLAUDE.md and dev/ notes are <end>'s, which is
# what the code under review was written against.

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "usage: bash dev/tool/review_worktree.sh <start> <end> <path>" >&2
  exit 2
fi
start=$1
end=$2
path=$3

if [ ! -f DESCRIPTION ] || ! grep -q '^Package: slideimp$' DESCRIPTION; then
  echo "run from the slideimp package root" >&2
  exit 2
fi
if [ -e "$path" ]; then
  echo "$path exists; pick another path or remove it first" >&2
  exit 2
fi

git rev-parse --verify --quiet "$start^{commit}" >/dev/null ||
  { echo "not a commit: $start" >&2; exit 2; }
git rev-parse --verify --quiet "$end^{commit}" >/dev/null ||
  { echo "not a commit: $end" >&2; exit 2; }
git merge-base --is-ancestor "$start" "$end" ||
  { echo "$start is not an ancestor of $end" >&2; exit 2; }

git worktree add --detach "$path" "$end"
git -C "$path" reset --soft "$start"

echo
echo "worktree: $path"
echo "HEAD:     $(git -C "$path" log -1 --format='%h %s' HEAD)"
echo "content:  $(git log -1 --format='%h %s' "$end")"
echo "commits:  $(git rev-list --count "$start..$end")"
git -C "$path" diff --cached --shortstat
