# vendor_spectra.R -- copy the Spectra headers mipca.cpp needs into src/Spectra/.
#
# Source: a clone of https://github.com/yixuan/spectra at dev/spectra, at
# the pinned commit below. The clone is ignored (dev/ is tracked opt-in) and
# per machine, so the script runs on either machine: when
# dev/spectra is absent it clones UPSTREAM there and checks out the pin
# (detached). A clone that is present is never moved: the script refuses if
# its HEAD is anything other than the pin.
#
# What it does, in order, stopping at the first failure:
#   1. clones if absent, then checks the clone and the pin;
#   2. walks the quoted #include graph from the two roots mipca.cpp includes,
#      inside the pinned tree, and refuses unless it reaches exactly FILES --
#      no file missing, none extra -- and the only non-std includes outside
#      Spectra are <Eigen/Core> and <Eigen/Jacobi>;
#   3. deletes src/Spectra/ and writes each file from `git cat-file blob`
#      (raw blob bytes, no filters);
#   4. re-hashes every written file with `git hash-object --no-filters` and
#      refuses unless it equals the upstream blob hash;
#   5. writes src/Spectra/MANIFEST: upstream, commit, and one blob hash per
#      file.
#
# Why blob hashes and not file hashes: a Windows checkout may convert line
# endings (core.autocrlf), so working-tree bytes are never what is compared.
# The upstream blobs carry no CR (checked here), so after `git add` the blob
# in this repo is the upstream blob and `git ls-files -s src/Spectra` shows
# the same hashes as the manifest.
#
# Patches: none. If one is ever needed (a Linux compiler warning inside a
# vendored header), it is applied here after step 4, recorded in this header,
# and the manifest hashes the patched bytes with the patch named beside them.
#
# From the repository root:
#   Rscript dev/tool/vendor_spectra.R

PIN <- "db1d5cc3279752ca7ea3e33da44ba2a85e4e4a95" # v1.2.0-4-gdb1d5cc
UPSTREAM <- "https://github.com/yixuan/spectra"
ROOTS <- c("SymEigsSolver.h", "MatOp/DenseSymMatProd.h")
FILES <- c(
  "SymEigsSolver.h",
  "HermEigsBase.h",
  "Util/Version.h",
  "Util/TypeTraits.h",
  "Util/SelectionRule.h",
  "Util/CompInfo.h",
  "Util/SimpleRandom.h",
  "MatOp/DenseSymMatProd.h",
  "MatOp/internal/ArnoldiOp.h",
  "LinAlg/Lanczos.h",
  "LinAlg/Arnoldi.h",
  "LinAlg/TridiagEigen.h",
  "LinAlg/UpperHessenbergQR.h",
  "LinAlg/Givens.h",
  "LinAlg/DoubleShiftQR.h"
)
EIGEN_OK <- c("Eigen/Core", "Eigen/Jacobi")

root <- getwd() # run from the repository root
stopifnot(file.exists(file.path(root, "DESCRIPTION")))
src <- file.path(root, "dev", "spectra") # ignored: dev/ is tracked opt-in
dst <- file.path(root, "src", "Spectra")

git <- function(...) {
  out <- suppressWarnings(system2("git", c(...), stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0L) {
    stop(
      "git ",
      paste(c(...), collapse = " "),
      " failed:\n",
      paste(out, collapse = "\n"),
      call. = FALSE
    )
  }
  out
}

# ---- 1. the clone and the pin ----------------------------------------------
if (!dir.exists(file.path(src, ".git"))) {
  if (file.exists(src)) {
    stop(
      src,
      " exists but is not a git clone; move it aside and re-run",
      call. = FALSE
    )
  }
  cat("no Spectra clone at", src, "-- cloning", UPSTREAM, "\n")
  git("clone", "--no-checkout", UPSTREAM, src)
  git("-C", src, "checkout", "--detach", PIN)
}
head <- git("-C", src, "rev-parse", "HEAD")
if (!identical(head, PIN)) {
  stop("dev/spectra is at ", head, ", not the pin ", PIN, call. = FALSE)
}
cat("source:", src, "at", git("-C", src, "describe", "--tags", PIN), "\n")

up <- function(f) paste0(PIN, ":include/Spectra/", f)
blob <- function(f) git("-C", src, "rev-parse", up(f))
read_up <- function(f) git("-C", src, "cat-file", "blob", up(f))

# ---- 2. the include closure, inside the pinned tree ------------------------
# Only preprocessor lines count: the headers also carry `/// #include` in
# their doc comments.
norm <- function(p) {
  parts <- strsplit(p, "/", fixed = TRUE)[[1L]]
  keep <- character(0)
  for (s in parts) {
    if (s == "..") {
      if (!length(keep)) {
        stop("include climbs out of Spectra: ", p, call. = FALSE)
      }
      keep <- keep[-length(keep)]
    } else if (s != ".") {
      keep <- c(keep, s)
    }
  }
  paste(keep, collapse = "/")
}
seen <- character(0)
angle <- character(0)
queue <- ROOTS
while (length(queue)) {
  f <- queue[1L]
  queue <- queue[-1L]
  if (f %in% seen) {
    next
  }
  seen <- c(seen, f)
  inc <- grep("^\\s*#\\s*include\\s*[<\"]", read_up(f), value = TRUE)
  q <- sub(
    "^\\s*#\\s*include\\s*\"([^\"]+)\".*$",
    "\\1",
    grep("\"", inc, value = TRUE)
  )
  a <- sub(
    "^\\s*#\\s*include\\s*<([^>]+)>.*$",
    "\\1",
    grep("<", inc, value = TRUE)
  )
  dir <- dirname(f)
  q <- vapply(q, function(x) norm(if (dir == ".") x else file.path(dir, x)), "")
  queue <- c(queue, setdiff(q, seen))
  angle <- union(angle, a)
}
missing <- setdiff(seen, FILES)
extra <- setdiff(FILES, seen)
if (length(missing) || length(extra)) {
  stop(
    "include closure differs from FILES.\n  reached, not listed: ",
    paste(missing, collapse = " "),
    "\n  listed, not reached: ",
    paste(extra, collapse = " "),
    call. = FALSE
  )
}
eig <- grep("^Eigen/|^unsupported/", angle, value = TRUE)
if (length(setdiff(eig, EIGEN_OK))) {
  stop(
    "vendored headers include Eigen modules beyond ",
    paste(EIGEN_OK, collapse = ", "),
    ": ",
    paste(setdiff(eig, EIGEN_OK), collapse = " "),
    call. = FALSE
  )
}
if (any(grepl("/", setdiff(angle, eig)))) {
  stop(
    "unexpected non-std angle includes: ",
    paste(grep("/", setdiff(angle, eig), value = TRUE), collapse = " "),
    call. = FALSE
  )
}
cat(
  "closure:",
  length(seen),
  "files from",
  paste(ROOTS, collapse = " + "),
  "\noutside Spectra:",
  paste(sort(angle), collapse = " "),
  "\n"
)

# ---- 3. write, 4. verify ---------------------------------------------------
unlink(dst, recursive = TRUE)
hashes <- character(length(FILES))
bytes <- 0
for (i in seq_along(FILES)) {
  f <- FILES[i]
  out <- file.path(dst, f)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  st <- system2("git", c("-C", src, "cat-file", "blob", up(f)), stdout = out)
  if (st != 0L) {
    stop("git cat-file failed for ", f, call. = FALSE)
  }
  raw <- readBin(out, "raw", file.size(out))
  if (any(raw == as.raw(13L))) {
    stop(f, " carries a CR byte", call. = FALSE)
  }
  if (any(raw > as.raw(127L))) {
    stop(f, " carries a non-ASCII byte", call. = FALSE)
  }
  want <- blob(f)
  got <- git("hash-object", "--no-filters", out)
  if (!identical(got, want)) {
    stop(
      f,
      ": written bytes hash to ",
      got,
      ", upstream blob is ",
      want,
      call. = FALSE
    )
  }
  hashes[i] <- want
  bytes <- bytes + length(raw)
}
written <- sort(list.files(dst, recursive = TRUE))
if (!identical(written, sort(FILES))) {
  stop("src/Spectra holds something other than FILES", call. = FALSE)
}

# ---- 5. manifest -----------------------------------------------------------
# Binary connection so the manifest is LF on Windows too.
con <- file(file.path(dst, "MANIFEST"), "wb")
on.exit(close(con), add = TRUE)
writeLines(
  c(
    "Spectra headers vendored by dev/tool/vendor_spectra.R. Do not edit by hand;",
    "re-run the script. Hashes are git blob hashes of the upstream files, which",
    "are also the blob hashes in this repository (git ls-files -s src/Spectra).",
    "",
    paste("upstream:", UPSTREAM),
    paste("commit:  ", PIN, "(v1.2.0-4-gdb1d5cc, untagged)"),
    "source:   include/Spectra/<path> -> src/Spectra/<path>",
    "license:  MPL-2.0, https://mozilla.org/MPL/2.0/",
    "patches:  none",
    "",
    paste(hashes, FILES, sep = "  ")
  ),
  con,
  sep = "\n",
  useBytes = TRUE
)

cat(sprintf(
  "wrote %d files, %d lines, %.1f KB to src/Spectra/, every blob hash verified\n",
  length(FILES),
  sum(vapply(file.path(dst, FILES), function(p) length(readLines(p)), 0L)),
  bytes / 1024
))
