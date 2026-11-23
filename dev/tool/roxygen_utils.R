# Roxygen tooling, sourced by dev/tool/roxygen.R. Base R and tools only; nothing
# here is package code.
#
# The @param type table and the R8 words are the ones dev/WRITING.md names.
# The fill-to-80 check is opt-in, since filling the roxygen to 80 columns is
# not one of slideimp's rules.
#
# Every function takes paths relative to the repo root, which is where the
# driver runs from.

# R/RcppExports.R is generated, so it is neither linted nor rewrapped.
ROXYGEN_SKIP <- "RcppExports.R"

# ---- reading blocks ---------------------------------------------------------

# The roxygen blocks of every R file in dir but the skipped ones, in file
# order.
package_blocks <- function(dir = "R", skip = ROXYGEN_SKIP) {
  files <- sort(list.files(dir, pattern = "[.][Rr]$", full.names = TRUE))
  files <- files[!basename(files) %in% skip]
  unlist(lapply(files, doc_blocks), recursive = FALSE)
}

# Each top-level expression of file with a #' line in the comments right
# above it: its name, first line, those comments, its code down to its
# body's brace, and the expression.
doc_blocks <- function(file) {
  lines <- readLines(file, warn = FALSE)
  exprs <- parse(file, keep.source = TRUE)
  refs <- attr(exprs, "srcref")
  blocks <- list()
  for (i in seq_along(exprs)) {
    start <- refs[[i]][1L]
    first <- start
    while (first > 1L && grepl("^\\s*#", lines[first - 1L])) {
      first <- first - 1L
    }
    above <- lines[seq_len(start - first) + first - 1L]
    if (any(grepl("^\\s*#'", above))) {
      end <- signature_end(exprs[[i]], refs[[i]][3L])
      blocks[[length(blocks) + 1L]] <- list(
        file = file,
        line = first,
        name = object_name(exprs[[i]]),
        comments = above,
        code = lines[start:end],
        expr = exprs[[i]]
      )
    }
  }
  blocks
}

# The line of a function body's opening brace, else the expression's last.
signature_end <- function(expr, last) {
  value <- assigned_function(expr)
  if (!is.null(value)) {
    brace <- attr(value[[3L]], "srcref")
    if (length(brace) > 0L) {
      return(brace[[1L]][1L])
    }
  }
  last
}

object_name <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-"))) {
    return(deparse(expr[[2L]]))
  }
  deparse(expr)[1L]
}

# The function call of `name <- function(...) body`, else NULL.
assigned_function <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], as.name("<-"))) {
    value <- expr[[3L]]
    if (is.call(value) && identical(value[[1L]], as.name("function"))) {
      return(value)
    }
  }
  NULL
}

# ---- the markdown dump ------------------------------------------------------

# Writes every roxygen block in dir to path as one markdown file, each with
# its object's code down to the end of its signature. For reading a whole
# pass at once.
roxygen_md <- function(path, dir = "R") {
  out <- c("# Roxygen in R/", "")
  file <- ""
  for (block in package_blocks(dir)) {
    if (block[["file"]] != file) {
      file <- block[["file"]]
      out <- c(out, paste0("## ", file), "")
    }
    out <- c(
      out,
      paste0("### ", block[["name"]], " (line ", block[["line"]], ")"),
      "",
      "```r",
      block[["comments"]],
      block[["code"]],
      "```",
      ""
    )
  }
  writeLines(out, path)
  invisible(path)
}

# ---- wrapping ---------------------------------------------------------------

# Rewraps the prose of every #' block in dir to width, in place, as air does
# code: a paragraph is joined and then filled line by line, a tag's later
# lines indented two spaces, a list item's under its text. A backtick span or
# a [link] is never split, so a line holding a long one can run past width.
# Examples, code fences, headings, tables, blank lines and a bare tag line
# (`@details`, `@section Title:`) stay as written. Each file keeps its own
# line endings: writeLines() alone writes CRLF on Windows, which turned an LF
# file into CRLF when this was measured. Returns the files it changed,
# invisibly.
wrap_roxygen <- function(dir = "R", width = 80L, skip = ROXYGEN_SKIP) {
  files <- sort(list.files(dir, pattern = "[.][Rr]$", full.names = TRUE))
  files <- files[!basename(files) %in% skip]
  changed <- character()
  for (file in files) {
    lines <- readLines(file, warn = FALSE)
    joined <- wrap_lines(lines, width)
    if (!identical(joined, lines)) {
      raw <- readBin(file, "raw", file.info(file)$size)
      eol <- if (any(raw == as.raw(13L))) "\r\n" else "\n"
      con <- file(file, open = "wb")
      writeLines(joined, con, sep = eol)
      close(con)
      changed <- c(changed, file)
    }
  }
  invisible(changed)
}

# The tags whose text is code, names or a title, kept as written to the
# next tag: a joined @inheritSection would merge its title into its text.
WRAP_VERBATIM <- c(
  "examples",
  "examplesIf",
  "usage",
  "eval",
  "evalRd",
  "evalNamespace",
  "rawRd",
  "rawNamespace",
  "import",
  "importFrom",
  "importClassesFrom",
  "importMethodsFrom",
  "useDynLib",
  "inherit",
  "inheritParams",
  "inheritSection",
  "inheritDotParams",
  "rdname",
  "name",
  "aliases",
  "keywords",
  "family",
  "concept",
  "export",
  "exportS3Method",
  "method",
  "docType",
  "encoding",
  "order",
  "template",
  "templateVar",
  "includeRmd"
)

# A file's lines with each run of #' lines filled to width. Inf joins.
wrap_lines <- function(lines, width) {
  rox <- grepl("^\\s*#'", lines)
  out <- character()
  i <- 1L
  while (i <= length(lines)) {
    if (!rox[i]) {
      out <- c(out, lines[i])
      i <- i + 1L
      next
    }
    end <- i
    while (end < length(lines) && rox[end + 1L]) {
      end <- end + 1L
    }
    out <- c(out, wrap_block(lines[i:end], width))
    i <- end + 1L
  }
  out
}

# One run of #' lines, filled. Each line is kept as written, opens a
# paragraph, or joins the open one; a paragraph is filled when the next
# line does not join it.
wrap_block <- function(lines, width) {
  out <- character()
  para <- NULL
  verbatim <- FALSE
  fence <- FALSE
  for (line in lines) {
    lead <- sub("#'.*", "", line)
    text <- sub("^\\s*#' ?", "", line)
    tag <- regmatches(text, regexpr("^@[A-Za-z]+", text))
    item <- regmatches(text, regexpr("^\\s*([-*+]|[0-9]+[.)])\\s+", text))
    is_fence <- grepl("^\\s*```", text)
    if (length(tag) == 1L && !fence) {
      verbatim <- sub("@", "", tag) %in% WRAP_VERBATIM
    }
    keep <- fence ||
      is_fence ||
      (verbatim && length(tag) == 0L) ||
      !grepl("\\S", text) ||
      grepl("^\\s*(#|\\||>)", text)
    joins <- !keep &&
      length(tag) == 0L &&
      length(item) == 0L &&
      !is.null(para)
    if (!joins && !is.null(para)) {
      out <- c(out, fill_tokens(para, width))
      para <- NULL
    }
    if (is_fence) {
      fence <- !fence
    }
    if (keep) {
      out <- c(out, line)
    } else if (joins) {
      # Markdown drops a joined line's indent, a code span's included.
      para[["text"]] <- paste(para[["text"]], trimws(text, "left"))
    } else if (length(tag) == 1L) {
      # A verbatim tag, or a bare one such as @details or @section Title:
      # alone, keeps its line.
      words <- switch(
        tag,
        "@param" = "^@param\\s+\\S+",
        "@section" = "^@section[^:]*:",
        "^@\\S+"
      )
      if (verbatim || !grepl("\\S", sub(words, "", text))) {
        out <- c(out, line)
      } else {
        para <- list(lead = lead, rest = "  ", text = text)
      }
    } else if (length(item) == 1L) {
      para <- list(lead = lead, rest = strrep(" ", nchar(item)), text = text)
    } else {
      para <- list(lead = lead, rest = "", text = text)
    }
  }
  if (!is.null(para)) {
    out <- c(out, fill_tokens(para, width))
  }
  out
}

# A paragraph's words filled to width, after its first line's own indent
# (a list item's) and its later lines' indent.
fill_tokens <- function(para, width) {
  text <- para[["text"]]
  indent <- regmatches(text, regexpr("^\\s*", text))
  tokens <- wrap_tokens(trimws(text))
  first <- paste0(para[["lead"]], "#' ", indent)
  rest <- paste0(para[["lead"]], "#' ", para[["rest"]])
  lines <- character()
  line <- paste0(first, tokens[1L])
  for (t in tokens[-1L]) {
    # A word that would open a line as a tag, a heading or a list item
    # stays on the line before.
    opens <- grepl("^(@|#|\\||>|```|([-*+]|[0-9]+[.)])$)", t)
    if (nchar(line) + 1L + nchar(t) > width && !opens) {
      lines <- c(lines, line)
      line <- paste0(rest, t)
    } else {
      line <- paste(line, t)
    }
  }
  c(lines, line)
}

# The words of text, a backtick span or a [link] kept as one.
wrap_tokens <- function(text) {
  out <- character()
  open <- NULL
  for (w in strsplit(text, " ", fixed = TRUE)[[1L]]) {
    if (is.null(open)) {
      if (!nzchar(w)) {
        next
      }
      open <- w
    } else {
      open <- paste(open, w)
    }
    if (span_closed(open)) {
      out <- c(out, open)
      open <- NULL
    }
  }
  c(out, open)
}

span_closed <- function(x) {
  count <- function(ch) lengths(regmatches(x, gregexpr(ch, x, fixed = TRUE)))
  count("`") %% 2L == 0L && count("[") <= count("]")
}

# ---- lint against slideimp's WRITING.md -------------------------------------

# The type fragments an @param opens with: dev/WRITING.md, section 4,
# "@param form", as regexes. The table lists what the package already uses,
# so a miss is a question for the audit, not an error.
DOC_TYPES <- c(
  "A double matrix\\.",
  "A numeric matrix\\.",
  "Integer\\.",
  "Integer or `NULL`\\.",
  "Logical\\.",
  "Character\\. One of ",
  "Numeric scalar between `0` and `1`\\.",
  "Numeric\\.",
  "Integer, numeric, or `NULL`\\.",
  "Optional character or integer vector\\.",
  # A package object: A `slideimp_results` object. and the other classes.
  "An? `[A-Za-z0-9_.]+` object\\.",
  "A data frame\\.",
  "Character\\.",
  "Optional character vector\\.",
  "A function\\.",
  "An environment\\.",
  "An object\\.",
  "An atomic vector\\.",
  "A list of (character|integer|numeric) vectors\\b",
  # A union: a kind from the table, then "or" before the last kind.
  paste0(
    "(Character|Numeric|Integer|Logical|A data frame|A list|A matrix|",
    "A function|An? `[A-Za-z0-9_.]+` object)[^.]* or [^.]+\\."
  ),
  "Not used\\.",
  "Further arguments passed to ",
  "Named objects\\."
)

# WRITING.md R9: words to delete. The test there is "delete the word and see
# whether the instruction changed"; this only finds the candidates.
DOC_FILLER <- "\\b(yourself|simply|just|manually|of course|you forgot)\\b"

# WRITING.md R8 names these as ours, not the reader's: "the kernel, the
# chain, the hole enumeration". Candidates only; R8's test is whether the
# word names something the reader's own code names.
DOC_INTERNAL <- "\\b(kernel|chain|hole enumeration)\\b"

# Every rule lint_block() checks, over the #' lines in dir: a data frame of
# topic, param and problem, with no rows when clean. With fill = TRUE it
# also reports a block that wrap_roxygen() would refill to 80 columns.
lint_roxygen <- function(dir = "R", fill = FALSE) {
  found <- lapply(package_blocks(dir), lint_block, fill = fill)
  out <- do.call(rbind, found)
  if (is.null(out)) {
    out <- data.frame(
      topic = character(),
      param = character(),
      problem = character()
    )
  }
  rownames(out) <- NULL
  out
}

lint_block <- function(block, fill = FALSE) {
  comments <- block[["comments"]]
  tags <- doc_tags(comments)
  defaults <- doc_defaults(block[["expr"]])
  found <- character()
  for (i in seq_len(nrow(tags))) {
    found <- c(
      found,
      tag_problems(tags[["tag"]][i], tags[["text"]][i], defaults)
    )
  }
  # A block that rewrapping to 80 columns would change. wrap_roxygen() fixes
  # it. Opt-in.
  rox <- comments[grepl("^\\s*#'", comments)]
  if (fill && !identical(wrap_lines(rox, 80L), rox)) {
    found <- c(found, "(block)" = "is not filled to 80 columns")
  }
  if (length(found) == 0L) {
    return(NULL)
  }
  data.frame(
    topic = block[["name"]],
    param = names(found),
    problem = unname(found)
  )
}

# The problems of one tag, each named by the formal or the tag it is about.
tag_problems <- function(tag, text, defaults) {
  about <- paste0("@", tag)
  out <- character()
  if (tag == "param") {
    about <- sub("\\s.*", "", text)
    text <- sub("^\\S+\\s*", "", text)
    if (!any(starts_with_regex(text, DOC_TYPES))) {
      out <- c(out, "type is not in WRITING.md's table")
    }
    if (grepl("Default is ", text, fixed = TRUE)) {
      out <- c(out, "says `Default is`, not `Defaults to`")
    }
    says <- grepl("Defaults to ", text, fixed = TRUE)
    if (about %in% names(defaults) && !defaults[[about]] && says) {
      out <- c(out, "states a default it does not have")
    }
    if (grepl("\\bRequired\\b", text)) {
      out <- c(out, "says Required")
    }
  }
  if (tag == "return") {
    out <- c(out, "is @return, not @returns")
  }
  if (!tag %in% c("examples", "examplesIf")) {
    for (pattern in c(filler = DOC_FILLER, internal = DOC_INTERNAL)) {
      words <- regmatches(
        text,
        gregexpr(pattern, text, ignore.case = TRUE)
      )[[1L]]
      if (length(words) > 0L) {
        kind <- if (identical(pattern, DOC_FILLER)) "R9 word" else "R8 word"
        out <- c(out, paste0(kind, ": ", paste(unique(words), collapse = ", ")))
      }
    }
  }
  names(out) <- rep(about, length(out))
  out
}

# A block's roxygen as a data frame of tag and text, its lines joined. The
# lines before the first tag take the tag "intro".
doc_tags <- function(comments) {
  text <- sub("^\\s*#' ?", "", comments[grepl("^\\s*#'", comments)])
  group <- cumsum(grepl("^@[A-Za-z]+", text))
  tags <- character()
  texts <- character()
  for (g in unique(group)) {
    part <- text[group == g]
    if (g == 0L) {
      tags <- c(tags, "intro")
    } else {
      tags <- c(tags, sub("^@([A-Za-z]+).*", "\\1", part[1L]))
      part[1L] <- sub("^@[A-Za-z]+\\s*", "", part[1L])
    }
    texts <- c(texts, gsub("\\s+", " ", trimws(paste(part, collapse = " "))))
  }
  data.frame(tag = tags, text = texts)
}

# Each formal of a block's function but `...`, TRUE where it has a default.
# Empty for a block that is not a function.
doc_defaults <- function(expr) {
  value <- assigned_function(expr)
  if (is.null(value)) {
    return(logical())
  }
  args <- as.list(value[[2L]])
  args <- args[names(args) != "..."]
  vapply(
    args,
    function(v) !(is.symbol(v) && !nzchar(as.character(v))),
    NA
  )
}

starts_with_regex <- function(text, patterns) {
  vapply(patterns, function(p) grepl(paste0("^", p), text), NA)
}

# ---- @seealso links ---------------------------------------------------------

# Every link of each @seealso in dir's Rd files: a data frame of topic,
# target and finding. "no topic" is a dangling link, an R CMD check WARNING
# (WRITING.md section 6). "one-way" means the topic linked to does not link
# back; slideimp does not require that, so it is a note for the audit. A
# link into another package ([stats::update()]) is not checked. Run it after
# devtools::document().
lint_seealso <- function(dir = "man") {
  rds <- list.files(dir, pattern = "\\.Rd$", full.names = TRUE)
  alias <- list()
  links <- list()
  for (rd in rds) {
    db <- tools::parse_Rd(rd)
    name <- trimws(paste(unlist(rd_nodes(db, "\\name")), collapse = ""))
    alias[[name]] <- trimws(vapply(
      rd_nodes(db, "\\alias"),
      function(e) paste(unlist(e), collapse = ""),
      ""
    ))
    see <- rd_nodes(db, "\\seealso")
    if (length(see) > 0L) {
      links[[name]] <- unique(rd_link_targets(see[[1L]]))
    }
  }
  # Each alias, named by the topic that owns it, looked up by name.
  owner_of <- stats::setNames(
    rep(names(alias), lengths(alias)),
    unlist(alias, use.names = FALSE)
  )
  owner <- function(a) unname(owner_of[sub("\\(\\)$", "", a)])
  found <- list()
  for (name in names(links)) {
    for (target in links[[name]]) {
      topic <- owner(target)
      if (is.na(topic)) {
        finding <- "no topic"
      } else {
        back <- vapply(links[[topic]], owner, "")
        finding <- if (name %in% back) "" else "one-way"
      }
      if (nzchar(finding)) {
        found[[length(found) + 1L]] <- data.frame(
          topic = name,
          target = if (is.na(topic)) target else topic,
          finding = finding
        )
      }
    }
  }
  out <- do.call(rbind, found)
  if (is.null(out)) {
    out <- data.frame(
      topic = character(),
      target = character(),
      finding = character()
    )
  }
  rownames(out) <- NULL
  out
}

# Every \link target under an Rd node that stays in this package: the option
# of \link[=x]{y}, else the text. \link[pkg:x]{y} goes to another package
# and is left out.
rd_link_targets <- function(node) {
  if (identical(attr(node, "Rd_tag"), "\\link")) {
    opt <- attr(node, "Rd_option")
    if (is.null(opt)) {
      return(paste(unlist(node), collapse = ""))
    }
    opt <- paste(unlist(opt), collapse = "")
    if (!startsWith(opt, "=")) {
      return(character())
    }
    return(sub("^=", "", opt))
  }
  if (is.list(node)) {
    return(unlist(lapply(node, rd_link_targets), use.names = FALSE))
  }
  character()
}

# The top-level nodes of a parsed Rd with tag.
rd_nodes <- function(db, tag) {
  db[vapply(db, function(e) identical(attr(e, "Rd_tag"), tag), NA)]
}
