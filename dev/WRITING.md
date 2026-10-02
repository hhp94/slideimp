# Writing rules: roxygen prose and cli messages

**This is the single source for how user-facing text is written in this package.** It is
self-contained on purpose: an agent picking up one batch of topics should need this file and
`CLAUDE.md`, nothing else.

**`CLAUDE.md` owns the invariants this file builds on, and this file does not restate them:** ASCII
everywhere, errors through cli, input validation through checkmate with an explicit `.var.name`,
and the fact that the exported API, the error surface and the documentation are public and pinned
by callers. Everything about how the text itself is written -- R1 to R9, the cli mechanics, the
roxygen template -- lives here.

Ported from the methylCIPHERv2 rules of the same name. What carried over is what is about English
and about cli and roxygen mechanics. What did not is everything specific to that package's objects,
its assets, its `@seealso` groups and its linters, none of which exist here.

---

## 0. Existing text is public API

**These rules bind text that is being written or edited. They are not a licence to sweep.** A cli
message is part of the error surface, and a help page is part of the documentation, and callers pin
on both (`CLAUDE.md`). Rewording a message that already ships so that it satisfies R3 is still a
change to public API, and it is made deliberately, with the maintainer, not in passing while
touching the file for something else.

So: a new message or a new `@param` follows this file in full. An existing one that you are already
changing for a real reason follows it too. An existing one that merely breaches a rule is reported,
not fixed.

---

## 1. Which channel: cli or `stop()`

**The line is audience, not transport.**

- A message about **input the user chose** is cli, whatever function raises it. This includes a
  message raised deep inside an internal helper, if the value it complains about arrived through an
  exported argument.
- A message about a **package defect** -- a "cannot happen" state, a missing dispatch branch, an
  internal invariant that broke -- is a plain `stop(..., call. = FALSE)`.

A defect message is **hard-coded and greppable**: a fixed prefix first, values appended after it,
so a bug report can be located from pasted text with no stack trace.

**A condition the caller is expected to handle carries a class.** `slideimp_infeasible` is the
model: callers and tests match on the class, not the text. A new condition that a caller could
reasonably want to catch gets a `slideimp_*` class at the `cli_abort()` site.

**Progress and status output is user-facing too.** It goes through `cli_inform()` (or cli's
progress API), not base `message()`, for the same reason errors do.

The code does not fully follow this section yet. Base `stop()` and `message()` on user-input paths
are known divergences, not precedent. Find them with:

```bash
git grep -nE '(^|[^_.a-z])(stop|warning|message)\(' -- R ':!R/dev-utils.R' ':!R/RcppExports.R'
```

Converting one is a change to the error surface (section 0), so it is proposed, not done in passing.

---

## 2. The English: R1 to R9

These bind **every word a user can see**: cli message text, roxygen prose, `vignettes/*.Rmd`,
`vignettes/articles/*.Rmd`, and `README.Rmd`. They do **not** bind code comments, `stop()` text
aimed at developers, or `dev/` docs. In those, ASCII `--` is still required and these rules do not
apply.

- **R1. ASD-STE100 Simplified Technical English.** One instruction per sentence. About 20 words for
  an instruction, 25 for a description. The simple word over the elaborate one. One word for one
  meaning. Articles present. No noun cluster longer than three words. No ambiguous `-ing` form.
- **R2. No first person and no "please". No contractions.** Prefer the data or the object as the
  subject over a bare passive. A function name is allowed as a grammatical subject (`group_imp()`
  splits `obj` ...). Second person is allowed where the user really is the actor.
- **R3. No `--` and no `;` in user-facing text.** Period, comma, colon, and a single spaced hyphen.
  This is an accessibility requirement, not a style preference. A `;` is how a long sentence gets
  smuggled past a period: prefer a period and a short new sentence. Find the current breaches with:

  ```bash
  git grep -nE "^#'.*( -- |;)" -- R ':!R/dev-utils.R'
  ```

- **R4. Describe the problem and give an actionable next step.** Never state that imputation
  continues or that nothing was stopped. Name the function to call **next**, never the one
  currently running.
- **R5. No vector sized by user input reaches cli uncapped.** The door is `fmt_trunc()` in
  `R/utils.R`, interpolated as a value (`{fmt_trunc(bad, 6)}`), never pasted into the template. The
  true total belongs in the lead line (`{length(bad)} column{?s} ...`), so the truncated list does
  not have to carry it. Where the cost is upstream of the render, cap the input, not the output.
- **R6. Every R language object in prose carries markup.** cli `{.arg}` / `{.fn}` / `{.val}` /
  `{.code}` / `{.cls}` / `{.path}` / `{.field}`; roxygen backticks and `[fn()]` links. **Only R
  language objects.** Domain and platform terms -- K-NN, PCA, LOBPCG, EPIC, 450K, WGBS, CpG, "beta
  value" -- stay plain prose.
- **R7. Keep a marked span short.** Mark the identifier, not the surrounding phrase and not a whole
  call with long arguments. roxygen2 renders a long backticked span badly and it breaks the PDF
  manual.
- **R8. No internal vocabulary.** Name an object the way the reader's own code names it. The test
  is mechanical: if the word is not a function name, an argument name, a class name, a column name
  in a returned object, or a word already in a message the user sees, it is jargon. Words that live
  in `dev/*_structure.md` and in `src/` -- the kernel, the chain, the eigenblock, the hole
  enumeration -- are ours, not the reader's. A word in the package `Description` (LOBPCG, warm
  start) has been put in front of the reader on purpose and passes.
- **R9. State the fix. Do not editorialise about the reader.** Cut the words that assign effort or
  blame: "yourself", "simply", "just", "manually", "you forgot", "of course". The test is
  mechanical: **delete the word and see whether the instruction changed.** If it did not, the word
  was not part of the instruction. An imperative is still the right mood, and second person is
  still allowed.

### Settled word choices

- **"data frame" is the concept and takes no markup. `data.frame` is the function or the class and
  takes markup.** Write `data frame`, plain, for the kind of object: `A data frame.` Write
  `data.frame()` marked, with parentheses, for the function. Write `{.cls data.frame}` in cli and
  `` `data.frame` `` in roxygen for the class. This is R6 applied to the right token: the concept is
  an English noun and names no R object.
- **A cli message does not hand out a recipe for transforming the data values.** Name the problem
  and the value the data should have. Do not give the code that converts the user's own matrix: no
  M-value to beta conversion, no rescaling, no `log2()`. A user who copies a numeric conversion and
  misapplies it gets a plausible but wrong imputation that traces back to our message, and the
  scale of their input is theirs to fix in their own pipeline. Name the likely cause instead. **A
  trivial structural fix that fails loudly is the exception and stays**: `t(obj)` for a transposed
  matrix, `as.matrix()` for a data frame, setting `colnames(obj)`. None of those can silently
  corrupt a result.
- **`slideimp.extra` is named, never linked.** It is not in `Imports` or `Suggests`, so a roxygen
  link to it cannot resolve, and a cli `{.fn slideimp.extra::x}` points at a package the reader may
  not have. Write it as `` `slideimp.extra` `` in roxygen and `{.pkg slideimp.extra}` in cli.

### What no rule covers

Wordiness is not rule-shaped. A sentence can break nothing above and still be bad, and the fix is a
judgement. **The instrument is the rendered page, not the tags.** After `devtools::document()`, read
the topic through `tools::Rd2txt()` or `?fn` under `load_all1()`. Weak sentences that are invisible
while writing the source are usually obvious on one reading of the output.

---

## 3. cli mechanics

- **`sprintf` output must never become cli input.** cli parses every message and bullet as a
  template, so a `{` arriving inside a built string is read as syntax. Hand data in as interpolated
  values, which are never re-parsed, or build the line with `cli::format_inline()`. `sprintf` is
  still right for a plain `stop()`.
- **Bind every `{?}` plural marker with an explicit `cli::qty()`** unless the quantity is the
  interpolation immediately before it. cli binds a marker to the last interpolated value earlier in
  the same string, and the quantity does not carry across elements of a `c()` message vector. Get
  this wrong and the handler throws in place of the real diagnostic, or silently pluralizes against
  the wrong value. Safe form: `"{cli::qty(length(bad))}Group{?s} {fmt_trunc(bad)} ..."`.
- **cli reflows whitespace.** A pre-aligned block collapses onto one line. Use `cli_verbatim()`
  where alignment matters. Inside `cli_abort()` / `cli_inform()`, bullets carry no alignment and
  each emits one self-contained bullet per row.
- **Tests assert the condition, never its wording.** `expect_error(..., class = "slideimp_*")`
  where the condition has a class, otherwise `expect_error()` alone. The wording is public API, and
  a deliberate rewording should be one change in `R/`, not a second hunt through `tests/`.

### How many bullets, and which

**Two bullets is the shape to aim for**, because a cause and a fix are usually two different facts.
Two is not a budget to spend, and several messages are correct at one bullet or none.

Bullets that **stay**:

- **The fix.** An action the reader takes. Every message wants exactly one, and an "Or ..."
  alternative merges into it unless the two fixes answer two different absences.
- **The offending value or location.** The columns that do not match, the unknown option names,
  the groups that failed.
- **The cause**, where the lead cannot carry it. Where the lead is one short clause, fold the cause
  in and save the bullet.

Bullets that **go**:

- **A mechanism of our own code.** How the package works is not a next step.
- **A restatement of the lead** from another angle.
- **A consequence the reader cannot act on**, unless it is the *cost* of ignoring a warning. A
  warning whose only effect is "these cells stay `NA`" keeps that line, because without it the
  reader does not know what the warning is worth.

**Name the offending input, and nothing else in list form.** The duplicate column names, the
unknown parameters, the groups with neither `k` nor `ncp`: the reader cannot look these up, and
`fmt_trunc()` is right for them. A list the reader could get from an exported function is that
function's job, and the message points at it instead.

**Nothing may refer across a list.** A list is the widest block in any message, so a bullet cannot
point back over one. Put the value in the bullet that uses it, or in the lead.

**An `inform` takes a lead.** The first line is the lead and the rest are bullets. A message that
opens with an `"i"` renders an info bullet with nothing above it.

**Do not abbreviate an API noun in a message.** Say the whole word the argument or function uses.

---

## 4. Roxygen: the template

### Tag order, every exported topic

```
Title Case Noun Phrase
(blank)
Description: one or more sentences. What it does.
(blank)
@inheritParams / @inheritSection   (only if it takes shared text)
@param      one per remaining formal, in signature order
@details    only if needed
@returns    always
@seealso    only if needed
@examples   or @examplesIf, always
@export
```

- **The title and description are the first two paragraphs, untagged.** That is the form the
  package uses, and it is what `CLAUDE.md`'s "`@title`, `@description`" means: every exported topic
  carries both, not that both tags are written out. An explicit `@description` is right only where
  the description has to sit after something else.
- **Titles are noun phrases.** `K-Nearest Neighbor Imputation for Numeric Matrices`, `Column Mean
  Imputation`. Title case, no trailing period. A verb-first title (`Calculate Matrix Column
  Variances`) is a breach.
- **`@returns`, not `@return`.**
- **Markdown is on** (`Roxygen: list(markdown = TRUE)`), so backticks and `[fn()]` links work.

### `@param` form

```
@param name <Type>. <What it is, one sentence.> [Defaults to <value>.]
```

The type fragment the package already uses, by kind of argument:

| kind | type fragment |
|---|---|
| matrix | `A numeric matrix.` |
| count or size | `Integer.` |
| count that may be absent | ``Integer or `NULL`.`` |
| flag | `Logical.` |
| one of a fixed set | ``Character. One of `"a"`, `"b"`, or `"c"`.`` |
| proportion | ``Numeric scalar between `0` and `1`.`` |
| other number | `Numeric.` |
| seed | ``Integer, numeric, or `NULL`.`` |
| optional column selector | `Optional character or integer vector.` |
| package object | ``A `slideimp_results` object.`` (and the other `slideimp_*` classes) |
| data frame | `A data frame.` |
| `...` on a method that ignores it | `Not used.` |

- **Where a default is stated, the form is `Defaults to <value>.`** Most params do not state one
  and leave it to the `\usage` line, so the sentence is optional. It is never stated for an argument
  that has no default, and an argument is never described as "Required."
- **One sentence of description.** If it needs two, the second usually belongs in `@details`. A
  long enumerated contract (the `group` forms, the `na_loc` formats) is the exception and may run
  to a short list.

### Examples

- **No `\dontrun{}`, ever.** Use `@examplesIf`, which keeps the example real and visible and only
  skips execution where the condition is false.
- **Paths that spin up mirai daemons get `@examplesIf interactive() && requireNamespace("mirai",
  quietly = TRUE)`.** Daemons are not something `R CMD check` should start.
- **Everything else runs unconditionally**, with no network and no `slideimp.extra`. Build inputs
  with `sim_mat()`.
- **An example never needs `slideimp.extra`.** It is not a declared dependency, so an example that
  uses it is an undeclared one.
- Examples do not share state across topics, so each block builds what it needs.

---

## 5. Shared text: write once, inherit

There is no dedicated donor topic. Shared text lives on the topic that owns it and is pulled in by
name: `@inheritParams knn_imp`, `@inheritParams pca_imp`, `@inheritParams group_imp`, and
`@inheritSection` from `slideimp-package` (`Missing values and non-finite input`), `pca_imp`
(`PCA Performance tips`) and `group_imp` (`Parallelization`). Re-read the sources rather than
trusting this list.

- **Shared text stays general and does not enumerate.** An enumeration in shared text is a list
  that has to stay true at every topic that inherits it.
- **Inheritance matches on the argument name alone.** This is the live footgun. `method` means a
  distance in `knn_imp()` and an algorithm in `pca_imp()`. `cores` is K-NN-only in some topics and
  general in others. `x`, `n` and `p` mean different things on each `print` method. A topic that
  inherits from a donor whose same-named argument means something else gets confidently wrong text,
  not an error. Override with a local `@param`, which always wins.
- **The package topic must keep `@keywords internal`.** It donates a section, so it has to exist as
  a topic, and `@keywords internal` keeps it out of the index without breaking the donation. `@noRd`
  produces no topic and silently drops every section it donates.

---

## 6. `@seealso`

- **A dangling link is an `R CMD check` WARNING.** Check every target exists after
  `devtools::document()`.
- **Never link to `slideimp.extra`** (section 2).
- Form: one sentence for a single link, a plain list for two or more, each item `[fn()]` plus a
  short clause saying what the reader gets there.

---

## 7. The prose files: vignettes and README

R1 to R9 bind them in full. Everything below is in addition.

### The three files are not built the same way

| | `vignettes/*.Rmd` | `vignettes/articles/*.Rmd` | `README.Rmd` |
|---|---|---|---|
| built by `R CMD check` and CRAN | yes | no | no |
| rendered by | check, pkgdown | pkgdown only | the maintainer |
| ships in the tarball | yes | no | no, but `README.md` does |

- **A `vignettes/*.Rmd` chunk that evaluates must run offline**, without `slideimp.extra`, quickly,
  and without starting mirai daemons. Anything else is `eval = FALSE` with its real output pasted
  below the call.
- **An article and `README.Rmd` may evaluate heavier work**, because nothing on CRAN runs them.
  Prefer a chunk that evaluates: real output cannot rot.
- **Pasted output is a claim about behaviour and rots like a `@details` sentence.** Run the call,
  then paste what it actually printed. Never retype it and never adjust it.
- **Pin anything random.** `set.seed()` before the first chunk that simulates, or the rendered file
  churns on every build.
- **Say what the reader does, not what the package does internally.** The pull toward internal
  vocabulary is strongest here, because the mechanism is interesting. R8 still applies.

### ASCII

The rule binds what an author types, not what a run prints.

- **Never hand-write a non-ASCII character**, in prose or in a chunk you author.
- **Non-ASCII that cli generated in captured output is fine and stays.** Do not set
  `options(cli.unicode = FALSE)` to launder it: the rendered file would then disagree with what the
  reader's own console shows.
- **Pandoc manufactures non-ASCII from ASCII input if it is allowed to.** Its `smart` extension
  turns a straight apostrophe or quote into a curly one. Disable it in the YAML of any `.Rmd`
  rendered to markdown:

  ```
  output:
    github_document:
      md_extensions: -smart
  ```

  The `CLAUDE.md` ASCII check covers `README.md`, so a smart quote that arrives this way shows up
  there, and the fix is the YAML, not hand-editing the generated file.

---

## 8. Auditing the manual: known-good exceptions

Everything in this section is intended. Reporting one of these as a defect is a false positive.

1. **`@examplesIf interactive() && ...` looks like an unguarded example, and is not.** `Rd2txt`
   prints the body with the guard stripped. The mirai-backed examples use it on purpose (section 4).
2. **Internal topics carry no `@returns` or `@examples`.** A `@keywords internal` topic generates
   no `\usage` block, so `R CMD check` asks for neither. An **exported** topic missing either is a
   real finding.
3. **`--` and `;` are required in this file and banned in package messages.** The R3 ban is scoped
   to text a user can see.

### What an auditor should actually check

Read for what no tool here sees: a `@param` sentence that does not match the formal it names, a
`@details` paragraph that describes behaviour the code no longer has, a title that is not a noun
phrase, an example that would need the network, `slideimp.extra` or a daemon, an inherited
`@param` whose donor means something else by the same name, and any breach of R1 to R9 in text a
user can see. Report breaches in existing text. Do not fix them (section 0).
