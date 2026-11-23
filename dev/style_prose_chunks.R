# Style the R chunks of the prose files with styler. air formats .R files
# only, so the chunks in the vignette, the articles and README.qmd are styled
# here. Prints the styler version and whether each file changed.
# Usage: Rscript dev/style_prose_chunks.R
cat("styler", format(packageVersion("styler")), "\n")
files <- c(
  Sys.glob("vignettes/*.qmd"),
  Sys.glob("vignettes/articles/*.qmd"),
  "README.qmd"
)
files <- files[file.exists(files)]
res <- styler::style_file(files)
print(res)
