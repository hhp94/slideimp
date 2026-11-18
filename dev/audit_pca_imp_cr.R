# Audit of a code review of e5b3cf7..3c84803, the si_mi conversion
# (dev/si_mi_conversion.md, "Review of Phases 1 to 6"): measures the
# findings that can be measured. Numbered as that section numbers them.
#   F1  output = "mipca": does as.matrix() carry the input's class/attributes?
#   F2  fallback_part(): is the map double when which() returns double?
#   F4  mean_imp_col(): R-heap cost of copying complete columns
#   F5  cost of one col_miss_internal() pass against a pca_imp() call, and
#       the integer input's extra copy
#
# Run from the package root:
#   Rscript dev/audit_pca_imp_cr.R

source("dev/dev-utils.R")
load_all1(timer = FALSE)

extra_attrs <- function(m) setdiff(names(attributes(m)), c("dim", "dimnames"))
show <- function(label, m) {
  cat(sprintf(
    "  %-44s class = %s | extra attrs = %s\n",
    label,
    paste(class(m), collapse = "/"),
    if (length(extra_attrs(m))) paste(extra_attrs(m), collapse = ", ") else "-"
  ))
}

cat("== F1: attributes through output = 'mipca' ==\n")
set.seed(1)
x <- sim_mat(60, 12, perc_total_na = 0.1, perc_col_na = 1)$input
dimnames(x) <- list(
  paste0("r", seq_len(nrow(x))),
  paste0("c", seq_len(ncol(x)))
)
r <- knn_imp(x, k = 5)
show("knn_imp() result (input to the cases)", r)
cat("  knn result imp_method:", attr(r, "imp_method"), "\n")

r2 <- r
r2[cbind(1:5, 1:5)] <- NA
show("knn result with 5 new holes", r2)
fit_a <- suppressWarnings(pca_imp(r2, ncp = 2, output = "mipca"))
cat("  (a) chunks in fit:", length(fit_a$chunks), "\n")
m_a <- as.matrix(fit_a)
show("(a) every column eligible: as.matrix(fit)", m_a)
cat(
  "      imp_method on as.matrix(fit):",
  format(attr(m_a, "imp_method")),
  "\n"
)
mc <- mipca_complete(fit_a, m = 1L)
show("(a) mipca_complete(fit, 1)[[1]]", if (is.list(mc)) mc[[1]] else mc)

fit_b <- suppressMessages(pca_imp(r, ncp = 2, output = "mipca"))
show("(b) no-holes exit: as.matrix(fit)", as.matrix(fit_b))

y <- x
attr(y, "foo") <- "bar"
show(
  "(c) plain matrix + attr foo: as.matrix(fit)",
  as.matrix(suppressWarnings(pca_imp(y, ncp = 2, output = "mipca")))
)

r3 <- r2
r3[1:58, 12] <- NA # 58 / 60 over colmax 0.9: one ineligible column
fit_d <- suppressWarnings(pca_imp(r3, ncp = 2, output = "mipca"))
cat("  (d) chunks in fit:", length(fit_d$chunks), "\n")
show("(d) one ineligible column: as.matrix(fit)", as.matrix(fit_d))
show(
  "(e) same as (a), output = 'matrix'",
  suppressWarnings(pca_imp(r2, ncp = 2))
)

cat("\n== F2: arrayInd() type with a double index ==\n")
cat("  R:", R.version.string, "\n")
b <- deparse(body(arrayInd))
cat("  body has storage.mode integer:", any(grepl("storage.mode", b)), "\n")
cat(paste0("    ", grep("storage.mode", b, value = TRUE)), sep = "\n")
ai_d <- arrayInd(c(1, 5, 6), c(2L, 3L))
ai_i <- arrayInd(c(1L, 5L, 6L), c(2L, 3L))
cat("  typeof(arrayInd(double ind)):", typeof(ai_d), "\n")
cat("  typeof(arrayInd(integer ind)):", typeof(ai_i), "\n")
cat("  identical:", identical(ai_d, ai_i), "\n")
# which() itself hands its index to arrayInd(); show the path
cat(
  "  which() calls arrayInd():",
  any(grepl("arrayInd", deparse(body(which)))),
  "\n"
)

cat("\n== F4: mean_imp_col() over every column vs only the holed ones ==\n")
set.seed(2)
n <- 2000L
p <- 5000L
z <- matrix(rnorm(n * p), n, p)
holed <- seq(1L, p, by = 10L) # 10% of the columns
for (j in holed) {
  z[sample.int(n, 100L), j] <- NA
}
X_mb <- as.numeric(object.size(z)) / 2^20
heap_extra <- function(f) {
  invisible(gc(reset = TRUE))
  before <- gc()["Vcells", 2]
  res <- f()
  after <- gc()["Vcells", 6]
  list(res = res, extra = (after - before) / X_mb)
}
# the part before the fix: every column with an observed value, complete
# ones included, which is what mean_imp_col() built
every <- heap_extra(function() {
  fill_parts(z, list(fallback_part(z, seq_len(p), "mean")))
})
now <- heap_extra(function() mean_imp_col(z))
cat(sprintf(
  "  input %.0f MB, %d of %d columns with holes\n",
  X_mb,
  length(holed),
  p
))
cat(sprintf(
  "  part over every column (before the fix): extra R-heap peak = %.2f x input\n",
  every$extra
))
cat(sprintf(
  "  mean_imp_col() as it is now:             extra R-heap peak = %.2f x input\n",
  now$extra
))
cat("  same result:", identical(every$res, now$res), "\n")

cat("\n== F5: one hole-count pass against a pca_imp() call ==\n")
t_pass <- system.time(
  for (i in 1:10) {
    col_miss_internal(z)
  }
)[["elapsed"]] /
  10
t_pca <- system.time(
  suppressWarnings(pca_imp(
    z,
    ncp = 2,
    solver = "exact",
    maxiter = 10,
    miniter = 1
  ))
)[["elapsed"]]
cat(sprintf("  col_miss_internal(): %.4f s per pass\n", t_pass))
cat(sprintf(
  "  pca_imp(maxiter = 10): %.2f s, so 2 passes = %.2f%% of it\n",
  t_pca,
  100 * 2 * t_pass / t_pca
))
small <- z[1:200, 1:100]
t_pass_s <- system.time(
  for (i in 1:1000) {
    col_miss_internal(small)
  }
)[["elapsed"]] /
  1000
t_mic_s <- system.time(
  for (i in 1:200) {
    mean_imp_col(small)
  }
)[["elapsed"]] /
  200
cat(sprintf(
  "  window-sized 200 x 100: pass %.6f s, mean_imp_col() %.6f s\n",
  t_pass_s,
  t_mic_s
))

cat("\n== F5b: integer input, R-heap peak of pca_imp() ==\n")
zi <- z[, 1:1000]
zi[] <- round(zi * 100)
zd <- zi
storage.mode(zi) <- "integer"
Xd_mb <- as.numeric(object.size(zd)) / 2^20
peak <- function(m) {
  invisible(gc(reset = TRUE))
  before <- gc()["Vcells", 2]
  invisible(suppressWarnings(pca_imp(
    m,
    ncp = 2,
    solver = "exact",
    maxiter = 5,
    miniter = 1
  )))
  (gc()["Vcells", 6] - before) / Xd_mb
}
cat(sprintf(
  "  double input:  extra R-heap peak = %.2f x the double size\n",
  peak(zd)
))
cat(sprintf(
  "  integer input: extra R-heap peak = %.2f x the double size\n",
  peak(zi)
))
