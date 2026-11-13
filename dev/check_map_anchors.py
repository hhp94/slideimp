# /// script
# requires-python = ">=3.10"
# dependencies = []
# ///
"""Check the code citations in the tracked notes against the code.

The structural maps (dev/knn_imp_structure.md, dev/pca_imp_structure.md) and
the other notes cite code by file and by a piece of the code itself, never by
line number. Line numbers move with every edit above them and with every
`air format .`, and a stale one still lands on a real line, so it misleads
without failing. A stale fragment matches nothing, which this script reports.

Usage, from the repository root:

    uv run dev/check_map_anchors.py             every tracked dev/*.md, and CLAUDE.md
    uv run dev/check_map_anchors.py FILE...     just these files
    uv run dev/check_map_anchors.py --list      also print where each fragment lands

Prints one line per problem as `file:line: message` and exits 1 if there is
any. What it checks:

1. Anchors. In prose, a citation is written in parentheses, as

       (`R/utils.R`: `resolve_subset <- function(`)
       (`R/utils.R`: `fragment one`, `fragment two`)       two places in one file
       (`R/utils.R`: `fragment`; `src/mat_stats.cpp`: `fragment`)   two files

   An anchor may wrap across lines of the note: a paragraph is read with its
   lines joined by one space, as markdown renders it. Each fragment is a
   piece of ONE line of code, matched as a fixed string, and must occur
   exactly once in the file. A fragment that occurs nowhere
   has gone stale; one that occurs more than once is ambiguous. Either way,
   choose a different fragment from the code being cited. Something that
   looks like an anchor but is not in this form is reported too.
2. Paths. Every path cited -- in an anchor, in backticks in prose, or in a
   call tree -- must be a file git tracks, or a new file it does not ignore,
   because anything a tracked file cites must itself be tracked. Directories
   are accepted in call trees. A path outside this repository (another
   package's header, say) is not checked unless it is written as an anchor.
   The one exception is a file the note itself says how to recover, by
   writing `git show <commit>:<path>` somewhere in it; that commit must
   actually have the file.
3. Call trees: fenced blocks with no language tag. Paths are written in full
   (`src/armaSVD.cpp`, not `armaSVD.cpp`), and every name written as a call,
   `name()` or `name<...>()`, on a line that names a file must occur in at
   least one of the files named on that line.
4. No line numbers: `R/knn_imp.R:96` and `:96-106` are errors anywhere.
"""

import argparse
import os
import re
import subprocess
import sys

PATH = r"(?:DESCRIPTION|NAMESPACE|[\w.-]*[/.][\w./-]*)"
FRAG = r"`[^`\n]+`"
GROUP = rf"`{PATH}`: {FRAG}(?:, {FRAG})*"
ANCHOR_RE = re.compile(rf"\({GROUP}(?:; {GROUP})*\)")
GROUP_RE = re.compile(rf"`({PATH})`: ({FRAG}(?:, {FRAG})*)")
CODE_SPAN_RE = re.compile(r"`([^`\n]+)`")
LOOSE_RE = re.compile(rf"`{PATH}`: `")

TOP_DIRS = r"(?:R|src|tests|dev|vignettes|inst|man|\.github)"
NUMS = r"\d+(?:-\d+)?(?:,\s?\d+(?:-\d+)?)*"
LINE_NUMBER_RE = re.compile(
    rf"(?:{TOP_DIRS}/[\w./-]*\w|[\w./-]+\.(?:R|cpp|h|hpp|c|Rmd|md)):{NUMS}"
    rf"|(?:^|(?<=[\s`(,])):{NUMS}\b"
)

PROSE_PATH_RE = re.compile(rf"`({TOP_DIRS}/[^`\s*<>]+\.\w+)`")
FENCE_PATH_RE = re.compile(rf"(?<![\w/.-])({TOP_DIRS}/[\w./-]*\w)")
BARE_FILE_RE = re.compile(r"(?<![\w/.-])([\w-]+\.(?:R|cpp|h|hpp))\b")
CALL_RE = re.compile(r"(?<![\w.])(\.?[A-Za-z_][\w.:]*)(?:<[^<>()]*>)?\(")
HISTORY_RE = re.compile(r"git show ([0-9a-f]{6,40}):([\w./-]+)")
NOT_NAMES = {"if", "for", "while", "switch", "function", "return", "sizeof"}


class Checker:
    def __init__(self, listing):
        self.listing = listing
        self.problems = []
        self.n_fragments = 0
        self.n_paths = 0
        self.texts = {}
        self.history = set()
        self.anchors = []
        # tracked, plus new files that .gitignore lets through
        out = subprocess.run(
            ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        self.tracked = set(p for p in out.split("\0") if p and os.path.isfile(p))

    def problem(self, doc, n, msg):
        self.problems.append(f"{doc}:{n}: {msg}")

    def text(self, path):
        if path not in self.texts:
            with open(path, encoding="utf-8", errors="replace") as fh:
                self.texts[path] = fh.read()
        return self.texts[path]

    def check_path(self, doc, n, path, allow_dir=False):
        """True when `path` is a tracked file whose text can be searched."""
        self.n_paths += 1
        path = path.removeprefix("./")
        if path in self.tracked:
            return True
        if path in self.history:
            return False
        if allow_dir and os.path.isdir(path):
            return False
        if os.path.exists(path):
            self.problem(doc, n, f"`{path}` exists but git ignores it")
        else:
            self.problem(doc, n, f"`{path}`: no such file")
        return False

    def check_doc(self, doc):
        with open(doc, encoding="utf-8", errors="replace") as fh:
            whole = fh.read()
        lines = whole.split("\n")
        self.history = set()
        for m in HISTORY_RE.finditer(whole):
            rev, path = m.groups()
            found = subprocess.run(
                ["git", "cat-file", "-e", f"{rev}:{path}"], capture_output=True
            )
            if found.returncode == 0:
                self.history.add(path)
            else:
                n = whole[: m.start()].count("\n") + 1
                self.problem(doc, n, f"`{rev}:{path}` is not in git history")
        in_fence = False
        tagged = False
        paragraph = []
        for n, line in enumerate(lines, 1):
            stripped = line.strip()
            if stripped.startswith("```"):
                self.check_prose(doc, paragraph)
                paragraph = []
                if in_fence:
                    in_fence = False
                else:
                    in_fence = True
                    tagged = bool(stripped[3:].strip())
                continue
            for m in LINE_NUMBER_RE.finditer(line):
                self.problem(
                    doc,
                    n,
                    f"line-number citation `{m.group(0)}`; cite a code fragment instead",
                )
            if in_fence:
                self.check_fence_line(doc, n, line, tagged)
            elif not stripped or stripped.startswith("#"):
                self.check_prose(doc, paragraph)
                paragraph = []
                if stripped:
                    self.check_prose(doc, [(n, stripped)])
            else:
                paragraph.append((n, stripped))
        self.check_prose(doc, paragraph)

    def check_prose(self, doc, paragraph):
        """Check one paragraph, its lines joined by a space as markdown
        renders them, so an anchor may wrap across lines."""
        if not paragraph:
            return
        text = ""
        starts = []
        for n, piece in paragraph:
            if text:
                text += " "
            starts.append((len(text), n))
            text += piece

        def line_at(pos):
            return [n for start, n in starts if start <= pos][-1]

        for am in ANCHOR_RE.finditer(text):
            n = line_at(am.start())
            for gm in GROUP_RE.finditer(am.group(0)):
                path = gm.group(1)
                frags = CODE_SPAN_RE.findall(gm.group(2))
                self.anchors.extend((path, f) for f in frags)
                if not self.check_path(doc, n, path):
                    continue
                code = self.text(path)
                for frag in frags:
                    self.n_fragments += 1
                    hits = code.count(frag)
                    if hits == 0:
                        self.problem(doc, n, f"`{frag}` not found in {path}")
                    elif hits > 1:
                        self.problem(
                            doc, n, f"`{frag}` occurs {hits} times in {path}"
                        )
                    elif self.listing:
                        at = code[: code.index(frag)].count("\n") + 1
                        print(f"{doc}:{n}: {path}:{at}  {frag}")
        masked = ANCHOR_RE.sub(lambda m: " " * len(m.group(0)), text)
        for m in LOOSE_RE.finditer(masked):
            self.problem(
                doc,
                line_at(m.start()),
                "looks like an anchor but is not in the (`path`: `fragment`) form: "
                + masked[m.start() :].strip()[:60],
            )
        for m in PROSE_PATH_RE.finditer(masked):
            self.check_path(doc, line_at(m.start()), m.group(1))

    def check_fence_line(self, doc, n, line, tagged):
        files = [
            p
            for p in FENCE_PATH_RE.findall(line)
            if self.check_path(doc, n, p, allow_dir=True)
        ]
        if tagged:
            return
        for m in BARE_FILE_RE.finditer(line):
            self.problem(doc, n, f"`{m.group(1)}`: write the path in full")
        if not files:
            return
        for name in CALL_RE.findall(line):
            if name in NOT_NAMES:
                continue
            pattern = re.compile(rf"(?<![\w.]){re.escape(name)}(?!\w)")
            if not any(pattern.search(self.text(p)) for p in files):
                self.problem(
                    doc, n, f"`{name}` is not in {', '.join(files)}"
                )


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if not os.path.exists("DESCRIPTION"):
        sys.exit("run this from the repository root")
    docs = args.files
    if not docs:
        out = subprocess.run(
            ["git", "ls-files", "--", "dev/*.md", "CLAUDE.md"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        docs = out.split()
    checker = Checker(args.list)
    for doc in docs:
        checker.check_doc(doc)
    for p in checker.problems:
        print(p)
    print(
        f"{len(docs)} files, {checker.n_fragments} fragments, "
        f"{checker.n_paths} paths: {len(checker.problems)} problems"
    )
    sys.exit(1 if checker.problems else 0)


if __name__ == "__main__":
    main()
