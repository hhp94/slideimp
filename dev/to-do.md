# To do

Work that is agreed on but not yet done. An entry stays here until the change
lands; then it is deleted, not marked done - git history is the record of what
happened, this file is the record of what has not.

An entry says what to change and why, enough that a session with no memory of
the conversation can pick it up. If an entry turns out to be wrong, delete it
and say so in chat.

## Draft `NEWS.md` for the release after 1.2.0

`NEWS.md` is written once per release, from the whole diff since the last
release, not one entry per change (see `CLAUDE.md`). Draft the section from
`git diff v1.2.0..HEAD` and `git log v1.2.0..HEAD`, with the maintainer, who
owns its style.

The development-version section already in `NEWS.md` was written one change at
a time, before this rule. Treat it as material for the draft, not as finished
text.

Easy to miss in the diff, noted so the draft does not:

- **Breaking:** the `na_check` argument of `knn_imp()`, `pca_imp()`,
  `slide_imp()` and `group_imp()` is removed, along with the
  `has_remaining_na` attribute it controlled and the note `print()` showed
  from it. `anyNA()` on the result answers the same question.
- `pca_imp()` on a matrix with no missing values now returns a
  `slideimp_results` with the same attributes as any other call, where it
  returned the bare input. `print()` shows it with the result header, and the
  call now copies the input once to attach the class.
