# Generate tests/testthat/fixtures/pca_recover_block_ref.rds: the kernel's
# imputed values and mse for three fixed-iteration problems, recorded from
# the code BEFORE SVD_triplet()'s wide branch was restricted to the changing
# columns. The test "wide recover_factors block agrees with the stored
# reference" in tests/testthat/test-pca_imp.R compares the current kernel
# against it on a tolerance: the fixture is a cross-build reference, never a
# bitwise one.
#
# Cases: wide with NA-free (fixed) columns, wide with an NA in every column
# (no fixed block, so the Gram cache is inactive), and tall with fixed
# columns, whose path the change does not touch. The inputs are stored in the
# fixture rather than regenerated so the test does not depend on the RNG.
# Row weights are non-uniform so the row-scaled Vn path is exercised.
#
# Re-running this against changed code records the changed code's answers,
# so regenerate it only on purpose and say so in the commit.
#
# Usage, from the repository root:
#   Rscript dev/gen_pca_recover_block_ref.R
if (!file.exists("DESCRIPTION")) stop("run from the repository root")

source("R/dev-utils.R")
load_all1(timer = FALSE)

N_ITER <- 10L
NCP <- 3L
NA_PER_COL <- 3L
RANK <- 3L
OUT <- "tests/testthat/fixtures/pca_recover_block_ref.rds"

make_input <- function(n, p, n_na_cols, seed) {
  set.seed(seed)
  x <- matrix(rnorm(n * RANK), n) %*% matrix(rnorm(RANK * p), RANK) +
    matrix(rnorm(n * p, sd = 0.5), n)
  na_cols <- sort(sample.int(p, n_na_cols))
  for (j in na_cols) x[sample.int(n, NA_PER_COL), j] <- NA
  row_w <- runif(n, 0.5, 1.5)
  list(x = x, row_w = row_w / sum(row_w))
}

specs <- list(
  wide_fixed = list(n = 24L, p = 60L, n_na_cols = 40L, seed = 20261028L),
  wide_all = list(n = 24L, p = 60L, n_na_cols = 60L, seed = 20261029L),
  tall_fixed = list(n = 60L, p = 24L, n_na_cols = 16L, seed = 20261030L)
)

cases <- lapply(names(specs), function(nm) {
  s <- specs[[nm]]
  d <- make_input(s$n, s$p, s$n_na_cols, s$seed)
  res <- pca_imp_internal_cpp(
    obj = d$x,
    eligible_idx = seq_len(ncol(d$x)) - 1L,
    ncp = NCP,
    scale = TRUE,
    regularized = TRUE,
    threshold = 0,
    init = 0L,
    maxiter = N_ITER,
    miniter = N_ITER,
    row_w = d$row_w,
    coeff_ridge = 1,
    solver = 0L,
    warmup_iters = 0L,
    lobpcg_tol = 1e-8,
    lobpcg_maxiter = 1L
  )
  n_fixed <- sum(colSums(is.na(d$x)) == 0L)
  cat(sprintf(
    "%-10s %3d x %3d, n_fixed %2d, %4d imputed cells, mse %.12g\n",
    nm, s$n, s$p, n_fixed, nrow(res$imputed_values), res$mse
  ))
  list(
    x = d$x,
    row_w = d$row_w,
    imputed_values = res$imputed_values,
    mse = res$mse
  )
})
names(cases) <- names(specs)

out <- list(
  provenance = list(
    commit = system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
    r = R.version.string,
    blas = extSoftVersion()[["BLAS"]],
    lapack = La_version()
  ),
  ncp = NCP,
  iters = N_ITER,
  cases = cases
)
dir.create(dirname(OUT), showWarnings = FALSE)
saveRDS(out, OUT, version = 2)
cat(sprintf("wrote %s (%d bytes)\n", OUT, file.size(OUT)))
