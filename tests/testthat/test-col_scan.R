skip_on_cran()

# The column scan (R/col_scan.R) and the policies the entry points apply to it.

scan_input <- function() {
  set.seed(11)
  X <- matrix(stats::rnorm(40 * 10), 40, 10)
  X[sample.int(length(X), 60)] <- NA_real_
  X[, 1] <- NA_real_ # no observed value
  X[-1, 2] <- NA_real_ # one observed value
  X[3, 3] <- Inf
  X[4, 4] <- -Inf
  X[5:9, 4] <- NaN
  X[, 5] <- 2 # constant
  X[7, 6] <- Inf
  X[8, 6] <- NA_real_
  X
}

test_that("col_scan() counts as base R does, in both of its loops", {
  X <- scan_input()
  with_var <- col_scan(X, var = TRUE, rows = TRUE)
  counts <- col_scan(X, rows = TRUE)

  expect_identical(with_var$n_row, nrow(X))
  expect_identical(with_var$n_miss, as.integer(colSums(is.na(X))))
  expect_identical(with_var$has_inf, unname(colSums(is.infinite(X)) > 0))
  expect_identical(with_var$row_miss, as.integer(rowSums(is.na(X))))
  # the loop without the variance is a second loop
  expect_identical(counts$n_miss, with_var$n_miss)
  expect_identical(counts$has_inf, with_var$has_inf)
  expect_identical(counts$row_miss, with_var$row_miss)
  expect_null(counts$var)
  expect_null(col_scan(X)$row_miss)
  expect_true(any(with_var$has_inf) && any(with_var$n_miss == nrow(X)))

  # the variance: NA below two observed values, NaN with an Inf, and
  # otherwise base R's to rounding, over the finite columns
  v <- with_var$var
  expect_true(all(is.na(v[1:2]) & !is.nan(v[1:2])))
  expect_true(all(is.nan(v[c(3, 4, 6)])))
  fin <- c(5, 7:10)
  ref <- apply(X[, fin], 2, stats::var, na.rm = TRUE)
  # the relative error is over the columns whose reference is not 0
  expect_true(any(ref != 0))
  abs_err <- max(abs(v[fin] - ref))
  rel_err <- max(abs(v[fin] - ref)[ref != 0] / abs(ref[ref != 0]))
  expect_lt(abs_err, 1e-12)
  expect_lt(rel_err, 1e-12)
  expect_identical(v[5], 0)
})

test_that("col_scan() counts the tails as base R does, in both of its loops", {
  X <- scan_input()
  X[10:11, 7] <- c(-0.5, 0.5) # one cell at each bound, in neither tail
  tl <- c(-0.5, 0.5)
  lo <- as.integer(colSums(X < tl[1], na.rm = TRUE))
  hi <- as.integer(colSums(X > tl[2], na.rm = TRUE))
  # -Inf is below, Inf above, and the constant 2 above in every row
  expect_gt(lo[4], 0L)
  expect_gt(hi[3], 0L)
  expect_identical(hi[5], nrow(X))
  for (v in c(FALSE, TRUE)) {
    scan <- col_scan(X, var = v, rows = TRUE, tails = tl)
    expect_identical(scan$n_below, lo, info = paste("var =", v))
    expect_identical(scan$n_above, hi, info = paste("var =", v))
    plain <- col_scan(X, var = v, rows = TRUE)
    rest <- setdiff(names(plain), c("n_below", "n_above"))
    expect_identical(scan[rest], plain[rest], info = paste("var =", v))
    expect_null(plain$n_below)
    expect_null(plain$n_above)
  }
  expect_error(
    col_scan(X, tails = 0.5),
    "internal error: col_scan() takes tails as two doubles",
    fixed = TRUE
  )
})

test_that("a scan of some columns is the full scan indexed", {
  X <- scan_input()
  j <- c(10L, 3L, 1L, 6L)
  full <- col_scan(X, var = TRUE, tails = c(-0.5, 0.5))
  fields <- c("n_miss", "has_inf", "var", "n_below", "n_above")
  expect_identical(
    col_scan(X, cols = j, var = TRUE, tails = c(-0.5, 0.5))[fields],
    scan_cols(full, j)[fields]
  )
  expect_null(scan_cols(col_scan(X, rows = TRUE), j)$row_miss)
  expect_length(col_scan(X, cols = integer(0))$n_miss, 0L)
})

test_that("col_scan() refuses a matrix that is not double", {
  Xi <- matrix(1:6, 2, 3)
  expect_error(
    col_scan(Xi),
    "internal error: col_scan() takes a double matrix",
    fixed = TRUE
  )
  expect_error(
    col_scan_cpp(Xi, 0:2, FALSE, FALSE, double(0)),
    "internal error: col_scan_cpp() takes a double matrix",
    fixed = TRUE
  )
  expect_error(
    col_scan_cpp(matrix(0, 2, 3), 3L, FALSE, FALSE, double(0)),
    "'cols' holds 3 at position 1",
    fixed = TRUE
  )
  expect_error(
    col_scan_cpp(matrix(0, 2, 3), 0L, FALSE, FALSE, 0.5),
    "internal error: col_scan_cpp() takes no tails or two",
    fixed = TRUE
  )
})

# A small input every entry point can run on: 30 x 12, named, with holes,
# proportions in (0, 1), and two groups.
door_input <- function() {
  set.seed(21)
  X <- matrix(stats::runif(30 * 12, 0.05, 0.95), 30, 12)
  dimnames(X) <- list(sprintf("r%02d", 1:30), sprintf("c%02d", 1:12))
  X[sample.int(length(X), 25)] <- NA_real_
  X
}
door_calls <- function() {
  grp <- data.frame(
    feature = sprintf("c%02d", 1:12),
    group = rep(c("a", "b"), each = 6)
  )
  mi <- function(Y) {
    mipca(
      Y,
      rank = 2L,
      transform = "identity",
      chains = 1L,
      iter_sampling = 2L,
      iter_warmup = 1L,
      seed = 1L
    )
  }
  list(
    knn_imp = function(Y) knn_imp(Y, k = 3),
    pca_imp = function(Y) pca_imp(Y, ncp = 2),
    group_imp = function(Y) group_imp(Y, group = grp, k = 2),
    slide_imp = function(Y) {
      slide_imp(
        Y,
        location = 1:12,
        window_size = 6,
        overlap_size = 2,
        min_window_n = 4,
        k = 2
      )
    },
    tune_imp = function(Y) {
      tune_imp(Y, parameters = data.frame(k = 2), .f = "knn_imp", n_reps = 1)
    },
    sample_na_loc = function(Y) sample_na_loc(Y, n_cols = 2, n_rows = 2),
    mat_miss = function(Y) mat_miss(Y),
    col_vars = function(Y) col_vars(Y),
    mipca = mi,
    mipca_tune = function(Y) {
      mipca_tune(
        Y,
        rank = 1:2,
        transform = "identity",
        n_cols = 3,
        n_rows = 2,
        chains = 1L,
        iter_sampling = 2L,
        iter_warmup = 1L,
        seed = 1L
      )
    },
    mipca_flank_plan = function(Y) {
      mipca_flank_plan(
        Y,
        location = 1:12,
        window_size = 4,
        subset = 1:3,
        min_window_n = 4
      )
    },
    abb_cols = function(Y) abb_cols(Y)
  )
}
quiet <- function(expr) suppressMessages(suppressWarnings(expr))
# The message of the error `expr` raises, or NULL when it raises none, so a
# loop over entry points can name the one that failed through `info`.
# expect_no_error() takes no label: its `message` is a pattern that selects
# which errors count.
error_of <- function(expr) {
  tryCatch(
    {
      expr
      NULL
    },
    error = function(e) conditionMessage(e)
  )
}

test_that("every entry point takes a double matrix and refuses any other", {
  X <- door_input()
  Xi <- round(X * 100)
  storage.mode(Xi) <- "integer"
  Xl <- is.na(X)
  calls <- door_calls()
  expect_length(calls, 12L)
  for (nm in names(calls)) {
    expect_null(error_of(quiet(calls[[nm]](X))), info = nm)
    expect_error(
      quiet(calls[[nm]](Xi)),
      "Assertion on 'obj' failed: Must store doubles",
      fixed = TRUE,
      info = nm
    )
    expect_error(
      quiet(calls[[nm]](Xl)),
      "Assertion on 'obj' failed: Must store doubles",
      fixed = TRUE,
      info = nm
    )
  }
})

test_that("the imputers refuse an Inf and an all-NA column with one message", {
  X <- door_input()
  Xinf <- X
  Xinf[4, 7] <- -Inf
  Xna <- X
  Xna[, 9] <- NA_real_
  calls <- door_calls()
  both <- c(
    "knn_imp",
    "pca_imp",
    "group_imp",
    "slide_imp",
    "tune_imp",
    "mipca",
    "mipca_tune"
  )
  for (nm in both) {
    for (Y in list(Xinf, Xna)) {
      expect_error(
        quiet(calls[[nm]](Y)),
        "'obj' holds infinite values or columns with no observed value",
        fixed = TRUE,
        info = nm
      )
    }
  }

  # the entry points that accept a column with no observed value refuse an
  # Inf with a message of its own
  for (nm in c("sample_na_loc", "mipca_flank_plan", "abb_cols")) {
    expect_error(
      quiet(calls[[nm]](Xinf)),
      "'obj' holds infinite values, which are not supported",
      fixed = TRUE,
      info = nm
    )
    expect_null(error_of(quiet(calls[[nm]](Xna))), info = nm)
  }

  # and the statistics refuse neither
  for (nm in c("mat_miss", "col_vars")) {
    expect_null(error_of(calls[[nm]](Xinf)), info = nm)
    expect_null(error_of(calls[[nm]](Xna)), info = nm)
  }
})

test_that("group_imp() and slide_imp() take one scan, not one per block", {
  X <- door_input()
  real_scan <- col_scan
  # the scans of the call that read an input, and those of the input
  # itself. A scan taken by fill_parts() reads what it wrote, to check the
  # write: each window's fill takes one, and the one write at the end.
  n_scan <- 0L
  n_input <- 0L
  local_mocked_bindings(
    col_scan = function(obj, ...) {
      caller <- as.character(sys.call(-1L)[[1L]])
      if (!identical(caller, "fill_parts")) {
        n_scan <<- n_scan + 1L
        n_input <<- n_input + identical(obj, X)
      }
      real_scan(obj, ...)
    }
  )
  calls <- door_calls()
  three_groups <- data.frame(
    feature = sprintf("c%02d", 1:12),
    group = rep(c("a", "b", "c"), each = 4)
  )
  runs <- list(
    group_knn = function() calls$group_imp(X),
    group_pca = function() group_imp(X, group = three_groups, ncp = 2),
    slide_knn = function() calls$slide_imp(X)
  )
  # each call reads its input once, at its door, whatever the number of
  # groups or windows. A scan of a group's slice or a window would be a
  # second.
  for (nm in names(runs)) {
    n_scan <- 0L
    n_input <- 0L
    quiet(runs[[nm]]())
    expect_identical(n_scan, 1L, info = nm)
    expect_identical(n_input, 1L, info = nm)
  }
})

test_that("mipca_flank() takes one scan of its windows' columns, not one per window", {
  X <- door_input()
  X[2:3, 1:3] <- NA_real_
  plan <- mipca_flank_plan(
    X,
    location = 1:12,
    window_size = 4,
    subset = 1:3,
    min_window_n = 4
  )
  expect_identical(sum(plan$targets$status == "impute"), 3L)
  held <- sort(unique(unlist(lapply(1:3, function(t) {
    match(colnames(mipca_flank_window(plan, t)), colnames(X))
  }))))
  # some columns lie outside every window, so a scan of every column
  # differs from one of the windows' columns
  expect_lt(length(held), ncol(X))
  real_scan <- col_scan
  n_scan <- 0L
  seen <- NULL
  local_mocked_bindings(
    col_scan = function(obj, cols = NULL, ...) {
      n_scan <<- n_scan + 1L
      seen <<- cols
      real_scan(obj, cols = cols, ...)
    }
  )
  fit <- quiet(mipca_flank(
    plan,
    rank = 1L,
    seed = 1L,
    transform = "identity",
    chains = 1L,
    iter_sampling = 2L,
    iter_warmup = 1L,
    .progress = FALSE
  ))
  expect_s3_class(fit, "si_mi")
  expect_identical(n_scan, 1L)
  expect_identical(seen, held)
})

sampler_input <- function() {
  set.seed(31)
  X <- matrix(stats::rnorm(30 * 6), 30, 6)
  dimnames(X) <- list(sprintf("r%02d", 1:30), sprintf("c%d", 1:6))
  X
}
mi_fit <- function(Y, transform = "identity", ...) {
  mipca(
    Y,
    rank = 2L,
    transform = transform,
    chains = 1L,
    iter_sampling = 2L,
    iter_warmup = 1L,
    seed = 1L,
    ...
  )
}
# the cells a finalized fit drew, over its chunks
n_drawn <- function(fit) {
  sum(vapply(fit$chunks, function(ch) NROW(ch$draws), 0L))
}

test_that("with nothing to impute, the sampler reads only the subset", {
  X <- sampler_input()
  X[5, 4] <- Inf
  # c1 has no hole, and c4's Inf is not read
  fit <- mi_fit(X, subset = "c1")
  expect_s3_class(fit, "si_mi")
  expect_identical(n_drawn(fit), 0L)
  # an Inf in a column read is refused
  expect_error(
    mi_fit(X, subset = "c4"),
    "'obj' holds infinite values or columns with no observed value",
    fixed = TRUE
  )
  # a fit that samples reads every column
  X[3, 1] <- NA_real_
  expect_error(
    mi_fit(X, subset = "c1"),
    "'obj' holds infinite values or columns with no observed value",
    fixed = TRUE
  )
})

test_that("a complete matrix whose total overflows has nothing to impute", {
  # the sampler used to test the grand total for a finite sum, sent this
  # matrix on to a fit, and stopped with "nothing to resume"
  X <- sampler_input()
  X[, 1] <- 1e307
  expect_false(is.finite(sum(X)))
  fit <- mi_fit(X)
  expect_s3_class(fit, "si_mi")
  expect_identical(n_drawn(fit), 0L)
})

test_that("a column listed in 'abb' is refused with the imputers' message", {
  X <- sampler_input()
  X[3:6, 2] <- NA_real_
  X[2, 5] <- Inf
  expect_error(
    mi_fit(X, abb = "c5"),
    "'obj' holds infinite values or columns with no observed value",
    fixed = TRUE
  )
  X[, 5] <- NA_real_
  expect_error(
    mi_fit(X, abb = "c5"),
    "'obj' holds infinite values or columns with no observed value",
    fixed = TRUE
  )
})

test_that("a transform whose forward() returns integer values is refused", {
  X <- sampler_input()
  X[3:6, 2] <- NA_real_
  tf <- list(
    name = "rounded",
    forward = function(y, ...) {
      z <- round(y * 10)
      storage.mode(z) <- "integer"
      z
    },
    backward = function(z, ...) z / 10
  )
  expect_error(
    mi_fit(X, transform = tf),
    "the transform's 'forward' returned values of type \"integer\"",
    fixed = TRUE
  )
})

test_that("a fit or a plan that holds a matrix that is not double is refused", {
  # Regression: such an object, made before every entry point took double
  # only, reached mipca_core() and col_scan() and stopped there with an
  # "internal error".
  X <- sampler_input()
  X[3:6, 2] <- NA_real_
  fit <- mi_fit(X, keep_state = TRUE)
  expect_null(error_of(update(fit, iter_sampling = 2L)))
  old <- mi_fit(X, keep_state = TRUE)
  old$obj <- round(old$obj)
  storage.mode(old$obj) <- "integer"
  expect_true(is.integer(old$obj))
  expect_error(
    update(old, iter_sampling = 2L),
    "this fit holds a matrix that is not double",
    fixed = TRUE
  )

  D <- door_input()
  plan <- mipca_flank_plan(
    D,
    location = 1:12,
    window_size = 4,
    subset = 1:3,
    min_window_n = 4
  )
  plan$obj <- round(plan$obj * 100)
  storage.mode(plan$obj) <- "integer"
  for (run in list(
    function() mipca_flank_window(plan, 1),
    function() {
      mipca_flank(plan, rank = 1L, seed = 1L, transform = "identity")
    }
  )) {
    expect_error(
      run(),
      "'plan' holds a matrix that is not double",
      fixed = TRUE
    )
  }
})
