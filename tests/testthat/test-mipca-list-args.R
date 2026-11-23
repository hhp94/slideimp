skip_on_cran()

# The sampler's kernels take their per-chain lists as SEXP and refuse anything
# that is not a list, so nothing is coerced through as.list() on the way in.

# arguments of the right type for every parameter of .mipca_da_loop_multi(),
# with no chain in them: the call gets past the list checks and stops at the
# chain count.
loop_args <- function() {
  list(
    xd_list = list(),
    g_fixed = numeric(0),
    etd_list = list(),
    xtild_list = list(),
    eig_start_list = list(),
    sigma_in = numeric(0),
    disp = numeric(0),
    resuming = FALSE,
    rng_state_in = character(0),
    miss_idx_l = integer(0),
    mi_row = integer(0),
    mi_col_l = integer(0),
    store_idx = integer(0),
    mi_col_store = integer(0),
    draws_list = list(),
    streams = integer(0),
    S = 1L,
    K = 1L,
    p = 1L,
    dfP = 1,
    nu = 1,
    warmup = 0L,
    ndraws = 1L,
    thin = 1L,
    gram_max_cond = 1e8,
    spectra_tol = 1e-10,
    spectra_ncv = 2L,
    spectra_maxitr = 10L,
    refresh = 0L
  )
}

test_that(".mipca_da_loop_multi() refuses a per-chain argument that is not a list", {
  # with lists everywhere the list checks pass, and the next check answers
  expect_error(
    do.call(.mipca_da_loop_multi, loop_args()),
    "internal error: the sampler needs at least one chain",
    fixed = TRUE
  )

  lists <- c(
    "xd_list",
    "etd_list",
    "xtild_list",
    "eig_start_list",
    "draws_list"
  )
  expect_true(all(lists %in% names(loop_args())))
  for (nm in lists) {
    args <- loop_args()
    # a matrix where the list of matrices belongs
    args[[nm]] <- matrix(0, 2L, 2L)
    expect_error(
      do.call(.mipca_da_loop_multi, args),
      sprintf("internal error: '%s' must be a list", nm),
      fixed = TRUE
    )
  }
})

test_that("mipca_cell_diag() refuses draws that are not a list", {
  expect_error(
    mipca_cell_diag(matrix(0, 2L, 2L), 1L, 256L),
    "internal error: 'draws' must be a list",
    fixed = TRUE
  )
  # a list of one chain's matrix is read, and returns one row per cell
  set.seed(20261122L)
  out <- mipca_cell_diag(list(matrix(stats::rnorm(40), 2L, 20L)), 1:2, 256L)
  expect_identical(dim(out), c(2L, 2L))
})
