# sample_na_loc ----
test_that("returns a list of length n_reps, each a 2-col row/col matrix", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 3, n_rows = 2, n_reps = 4)

  expect_type(out, "list")
  expect_length(out, 4)
  for (rep in out) {
    expect_equal(nrow(rep), 3 * 2)
  }
})

# n_cols / n_rows semantics
test_that("each rep uses exactly n_cols distinct columns, each with n_rows NAs", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 3, n_reps = 5)
  for (rep in out) {
    tab <- table(rep[, "col"])
    expect_length(tab, 4) # exactly n_cols distinct columns
    expect_true(all(tab == 3)) # each with exactly n_rows NAs
  }
})

test_that("no duplicated (row, col) pairs within a rep", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 3, n_reps = 5)
  for (rep in out) {
    key <- paste(rep[, "row"], rep[, "col"], sep = ":")
    expect_false(anyDuplicated(key) > 0)
  }
})

# num_na distribution
test_that("num_na divisible by n_rows distributes evenly", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, num_na = 12, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 12)
  tab <- table(rep[, "col"])
  expect_length(tab, 4) # 12 / 3 = 4 columns
  expect_true(all(tab == 3))
})

test_that("num_na not divisible by n_rows bumps last buckets by at least 1", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  # num_na = 13, n_rows = 3 -> n_cols = 4, na_per_col = c(4, 3, 3, 3)
  # (remainder = 1 column gets a +1)
  out <- sample_na_loc(m, num_na = 13, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 13)
  tab <- sort(as.integer(table(rep[, "col"])))
  expect_equal(tab, c(3L, 3L, 3L, 4L))
})

# na_col_subset handling
test_that("numeric na_col_subset restricts columns to the pool", {
  m <- sim_mat(p = 10, perc_total_na = 0)$input
  pool <- c(2L, 4L, 6L, 8L)
  out <- sample_na_loc(
    m,
    n_cols = 3,
    n_rows = 2,
    na_col_subset = pool,
    n_reps = 10
  )
  used <- unique(unlist(lapply(out, function(r) r[, "col"])))
  expect_true(all(used %in% pool))
})

test_that("character na_col_subset resolves via colnames", {
  m <- sim_mat(p = 6, perc_total_na = 0)$input # colnames feature1 ... feature6
  out <- sample_na_loc(
    m,
    n_cols = 2,
    n_rows = 1,
    na_col_subset = c("feature2", "feature5"),
    n_reps = 5
  )
  used <- unique(unlist(lapply(out, function(r) r[, "col"])))
  expect_setequal(used, c(2L, 5L))
})

# zero-variance cols protection
test_that("columns keep >= 2 distinct observed values after injection", {
  m <- sim_mat(20, p = 6, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 5, n_reps = 10)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    touched <- unique(rep[, "col"])
    for (j in touched) {
      obs <- m2[, j][!is.na(m2[, j])]
      expect_gte(length(unique(obs)), 2L)
    }
  }
})

test_that("columns without enough eligible rows (after keeping 2 uniques) are skipped", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  # Column 3 has only 4 observed values with exactly 2 uniques.
  # After keeping 2 (one of each), left_over = 2 < needed = 3 -> must be skipped.
  m[1:4, 3] <- c(1, 1, 2, 2)
  m[5:20, 3] <- NA
  needed <- 3L

  # colmax = 1, so the column cap is not what excludes column 3. At the
  # default 0.9 its room under colmax is 18 - 16 = 2 < 3, and the column would
  # be left out for that reason, before the count of observed rows is looked
  # at. Both halves are asserted so the test cannot drift back to that.
  colmax <- 1
  col_room <- max_miss_count(nrow(m), colmax) - sum(is.na(m[, 3]))
  expect_gte(col_room, needed)
  expect_lt(sum(!is.na(m[, 3])) - 2L, needed)

  out <- sample_na_loc(
    m,
    n_cols = 4,
    n_rows = needed,
    colmax = colmax,
    na_col_subset = 1:5,
    n_reps = 10
  )

  for (rep in out) {
    expect_false(
      3L %in% rep[, "col"],
      label = "column 3 should never be selected (not enough sacrificable rows)"
    )
  }
  # with column 3 the only one in the pool, the request cannot be met
  expect_error(
    sample_na_loc(
      m,
      n_cols = 1,
      n_rows = needed,
      colmax = colmax,
      na_col_subset = 3
    ),
    class = "slideimp_infeasible"
  )
})

test_that("a column that cannot take NAs is left out, in the pool or not", {
  m <- sim_mat(20, 6, perc_total_na = 0)$input
  m[, 3] <- 7 # constant
  m[1:19, 5] <- NA # over the default colmax of 0.9, one observed value

  # Regression: both columns used to abort the whole call ("Some columns
  # already have zero (or NA) variance", "Some columns already exceed
  # `colmax`"), even with a pool that does not contain them.
  set.seed(11)
  out <- sample_na_loc(
    m,
    n_cols = 2,
    n_rows = 1,
    na_col_subset = c(1, 2, 4, 6),
    n_reps = 5
  )
  expect_length(out, 5L)

  # in the pool they are never selected: 4 of the 6 columns can take an NA,
  # and all 4 are used
  out <- sample_na_loc(m, n_cols = 4, n_rows = 1, n_reps = 20)
  used <- lapply(out, function(r) sort(unname(r[, "col"])))
  expect_true(all(vapply(used, identical, logical(1), c(1L, 2L, 4L, 6L))))

  # one column more than can take an NA is refused up front, by count
  err <- expect_error(
    sample_na_loc(m, n_cols = 5, n_rows = 1),
    class = "slideimp_infeasible"
  )
  expect_match(conditionMessage(err), "only 4 columns can", fixed = TRUE)
})

test_that("an infeasible request stops before sampling, naming the counts", {
  m <- sim_mat(20, 8, perc_total_na = 0)$input

  # Regression: colmax = 0 leaves no column able to take an NA. This used to
  # run every attempt and stop with an unclassed "Failed to sample NA
  # locations after 100 attempts".
  err <- expect_error(
    sample_na_loc(m, n_cols = 2, n_rows = 1, colmax = 0),
    class = "slideimp_infeasible"
  )
  expect_match(conditionMessage(err), "only 0 columns can", fixed = TRUE)

  # unequal demands: num_na = 11, n_rows = 3 asks for c(3, 4, 4). Leave room
  # for 4 NAs in column 1 only, and for 3 in the rest.
  colmax <- 0.5 # cap = 10
  m[1:6, 1] <- NA # room 4
  m[1:7, 2:8] <- NA # room 3
  expect_identical(max_miss_count(nrow(m), colmax), 10L)
  # rowmax = 1: rows 1 to 6 are missing in every column
  err <- expect_error(
    sample_na_loc(m, num_na = 11, n_rows = 3, colmax = colmax, rowmax = 1),
    class = "slideimp_infeasible"
  )
  expect_match(
    conditionMessage(err),
    "2 columns must take at least 4 NAs each, but only 1 column can",
    fixed = TRUE
  )
  # one NA fewer asks for c(3, 3, 4), which fits
  set.seed(12)
  out <- sample_na_loc(m, num_na = 10, n_rows = 3, colmax = colmax, rowmax = 1)
  tab <- table(out[[1]][, "col"])
  expect_identical(sort(as.integer(tab)), c(3L, 3L, 4L))
  # the demand of 4 can only have gone to column 1
  expect_identical(names(tab)[tab == 4L], "1")
})

test_that("a non-finite cell is reported as Inf, not as zero variance", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  m[3, 4] <- Inf

  # col_vars() returns NaN for this column, which the zero-variance pre-check
  # used to claim as a variance problem. The Inf gate now runs first.
  expect_true(is.nan(col_vars(m)[4]))
  expect_error(
    sample_na_loc(m, n_cols = 4, n_rows = 1, na_col_subset = 1:5),
    "Infinite"
  )
})

test_that("aborts when requested n_cols exceeds available pool", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  err <- expect_error(
    sample_na_loc(m, n_cols = 10, n_rows = 1, na_col_subset = 1:5),
    class = "slideimp_infeasible"
  )
  expect_match(
    conditionMessage(err),
    "Cannot place NAs in 10 columns",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "na_col_subset", fixed = TRUE)

  # Regression: without `na_col_subset` the message still blamed it
  # ("`na_col_subset` provides only 5 columns"), and an `n_cols` request was
  # reported in the `num_na` wording.
  err <- expect_error(
    sample_na_loc(m, n_cols = 10, n_rows = 1),
    class = "slideimp_infeasible"
  )
  expect_match(conditionMessage(err), "`obj` has only 5 columns", fixed = TRUE)
  expect_no_match(conditionMessage(err), "na_col_subset", fixed = TRUE)

  # the `num_na` wording stays for a `num_na` request
  err <- expect_error(
    sample_na_loc(m, num_na = 12, n_rows = 2),
    class = "slideimp_infeasible"
  )
  expect_match(conditionMessage(err), "Cannot place 12 NAs", fixed = TRUE)
})

# row / col budget enforcement
test_that("resulting matrix respects colmax", {
  m <- sim_mat(20, 8, perc_total_na = 0)$input
  colmax <- 0.5
  out <- sample_na_loc(m, n_cols = 4, n_rows = 5, colmax = colmax, n_reps = 10)
  cap <- floor(nrow(m) * colmax)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    col_miss <- colSums(is.na(m2))
    expect_true(all(col_miss <= cap))
  }
})

test_that("resulting matrix respects rowmax", {
  m <- sim_mat(30, 10, perc_total_na = 0)$input
  rowmax <- 0.4
  out <- sample_na_loc(m, n_cols = 5, n_rows = 4, rowmax = rowmax, n_reps = 10)
  cap <- floor(ncol(m) * rowmax)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    row_miss <- rowSums(is.na(m2))
    expect_true(all(row_miss <= cap))
  }
})

# rep independence / reproducibility
test_that("reps are independently sampled (not identical)", {
  m <- sim_mat(30, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 4, n_rows = 3, n_reps = 5)
  # Hash each rep; at least 2 distinct with high probability.
  sigs <- vapply(
    out,
    function(r) {
      paste(r[, "row"], r[, "col"], collapse = ",")
    },
    character(1)
  )
  expect_gt(length(unique(sigs)), 1)
})

test_that("set.seed makes sample_na_loc reproducible", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  set.seed(42)
  a <- sample_na_loc(m, n_cols = 3, n_rows = 2, n_reps = 3)
  set.seed(42)
  b <- sample_na_loc(m, n_cols = 3, n_rows = 2, n_reps = 3)
  expect_equal(a, b)
})

# failure path
test_that("aborts when budgets make sampling infeasible", {
  # 10 rows allows n_rows=3 under the "keep >=2 values" hard bound,
  # but colmax=0.1 -> floor(10*0.1)=1 NA per column, so needing 3 is impossible.
  m <- sim_mat(10, 6, perc_total_na = 0)$input
  expect_error(
    sample_na_loc(m, n_cols = 2, n_rows = 3, colmax = 0.1, max_attempts = 3),
    class = "slideimp_infeasible"
  )
})

test_that("a search that runs out of attempts raises slideimp_infeasible", {
  # Every column can take an NA, so the request passes the up-front check,
  # but rowmax = 0 leaves no row able to take one. Only the search can find
  # that out. Regression: this was an Rcpp::exception from Rcpp::stop().
  m <- sim_mat(10, 6, perc_total_na = 0)$input
  err <- expect_error(
    sample_na_loc(
      m,
      n_cols = 2,
      n_rows = 1,
      rowmax = 0,
      max_attempts = 3,
      n_reps = 2
    ),
    class = "slideimp_infeasible"
  )
  expect_s3_class(err, "rlang_error")
  expect_match(
    conditionMessage(err),
    "repetition 1 in 3 attempts",
    fixed = TRUE
  )
})

# pre-existing NAs
test_that("sampled positions never collide with pre-existing NAs", {
  set.seed(1)
  m <- sim_mat(30, 10, perc_total_na = 0)$input
  # scatter some NAs, keeping columns healthy
  preset <- cbind(row = c(1, 2, 3, 4, 5), col = c(1, 2, 3, 4, 5))
  m[preset] <- NA

  out <- sample_na_loc(m, n_cols = 5, n_rows = 3, n_reps = 20)
  for (rep in out) {
    # none of the sampled (row, col) pairs should already be NA
    expect_true(all(!is.na(m[rep])))
  }
})

test_that("final colmax / rowmax account for pre-existing NAs", {
  set.seed(2)
  m <- sim_mat(20, 8, perc_total_na = 0)$input
  # preload column 1 with 4 NAs and row 1 with 2 NAs (avoiding overlap with col 1)
  m[1:4, 1] <- NA
  m[1, 3:4] <- NA

  colmax <- 0.5 # cap = floor(20 * 0.5) = 10
  rowmax <- 0.5 # cap = floor(8 * 0.5)  = 4

  out <- sample_na_loc(
    m,
    n_cols = 4,
    n_rows = 3,
    colmax = colmax,
    rowmax = rowmax,
    n_reps = 20
  )
  col_cap <- floor(nrow(m) * colmax)
  row_cap <- floor(ncol(m) * rowmax)
  for (rep in out) {
    m2 <- apply_na(m, rep)
    expect_true(all(colSums(is.na(m2)) <= col_cap))
    expect_true(all(rowSums(is.na(m2)) <= row_cap))
  }
})

test_that("columns with exhausted col_room are skipped even if individually healthy", {
  set.seed(3)
  m <- sim_mat(20, 6, perc_total_na = 0)$input
  colmax <- 0.5 # cap = 10
  # column 2 already has 9 NAs -> col_room = 1, so needed = 3 can't fit
  na_rows <- sample.int(20, 9)
  m[na_rows, 2] <- NA

  out <- sample_na_loc(
    m,
    n_cols = 4,
    n_rows = 3,
    colmax = colmax,
    n_reps = 20
  )
  for (rep in out) {
    expect_false(2L %in% rep[, "col"])
  }
})

test_that("num_na with remainder > n_cols distributes via larger bumps", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  # num_na = 5, n_rows = 3 -> n_cols = 1, one column takes all 5
  out <- sample_na_loc(m, num_na = 5, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 5)
  tab <- as.integer(table(rep[, "col"]))
  expect_equal(tab, 5L)
})

test_that("num_na = 11 with n_rows = 3 yields c(3, 4, 4)", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, num_na = 11, n_rows = 3, n_reps = 1)
  rep <- out[[1]]
  expect_equal(nrow(rep), 11)
  expect_equal(sort(as.integer(table(rep[, "col"]))), c(3L, 4L, 4L))
})

# default num_na
test_that("with num_na and n_cols both NULL, the default fits the columns", {
  set.seed(21)
  m <- matrix(rnorm(200 * 60), nrow = 200)

  # Regression: the default was min(500, 5% of the cells) whatever the column
  # count. 12000 cells gave 500 NAs, which at n_rows = 2 needs 250 columns:
  # "Cannot place 500 NAs (requires 250 columns ...)" on a 60-column matrix.
  expect_message(
    out <- sample_na_loc(m),
    "Using default `num_na` = 120",
    fixed = TRUE
  )
  expect_identical(nrow(out[[1]]), 120L)
  expect_true(all(table(out[[1]][, "col"]) == 2L))

  # a column that cannot take NAs does not count toward the default
  m[, 7] <- 1
  expect_message(out <- sample_na_loc(m), "= 118", fixed = TRUE)
  expect_false(7L %in% out[[1]][, "col"])

  # tune_imp() with every sampling argument at its default
  res <- suppressMessages(tune_imp(
    m,
    data.frame(k = 2),
    .f = "knn_imp",
    .progress = FALSE
  ))
  expect_identical(nrow(res$result[[1]]), 118L)
  expect_true(is.na(res$error))
})

test_that("default_num_na() works in double and in whole columns", {
  # 5% of the cells, capped at 500
  expect_identical(default_num_na(1e6, n_rows = 2, n_fit = 1000), 500L)
  expect_identical(default_num_na(1000, n_rows = 2, n_fit = 1000), 50L)
  # capped by the columns that can take `n_rows` NAs
  expect_identical(default_num_na(1e6, n_rows = 2, n_fit = 10), 20L)
  # whole columns: 500 %/% 3 = 166 columns of 3
  expect_identical(default_num_na(1e6, n_rows = 3, n_fit = 1000), 498L)
  # never fewer than one column
  expect_identical(default_num_na(10, n_rows = 2, n_fit = 5), 2L)
  expect_identical(default_num_na(1e6, n_rows = 600, n_fit = 5), 600L)

  # Regression: nrow(obj) * ncol(obj) is integer arithmetic and is NA past
  # 2^31 - 1 cells, which turned the default into NA. 2400 x 900000 is an
  # EPIC-sized cohort. length() of such a matrix is a double.
  n_total <- 2400 * 900000
  expect_gt(n_total, .Machine$integer.max)
  expect_true(is.na(suppressWarnings(2400L * 900000L)))
  expect_no_warning(big <- default_num_na(n_total, n_rows = 2, n_fit = 9e5))
  expect_identical(big, 500L)
})

# caps
test_that("max_miss_count() is the largest count the methods admit", {
  # Regression: floor(100 * 0.29) is 28, but 29 / 100 <= 0.29 is TRUE, which
  # is the test knn_imp() and pca_imp() apply to `colmax`.
  expect_identical(floor(100 * 0.29), 28)
  expect_true(29 / 100 <= 0.29)
  expect_identical(max_miss_count(100, 0.29), 29L)

  # against the definition, over a grid: the largest k with k / n <= prop
  worst <- 0L
  n_checked <- 0L
  for (n in c(1L, 2L, 3L, 7L, 10L, 20L, 71L, 100L, 1000L)) {
    for (prop in c(0, 0.1, 0.29, 0.3, 0.5, 0.57, 0.7, 0.9, 0.99, 1)) {
      want <- max(which((0:n) / n <= prop)) - 1L
      worst <- max(worst, abs(max_miss_count(n, prop) - want))
      n_checked <- n_checked + 1L
    }
  }
  expect_identical(n_checked, 90L)
  expect_identical(worst, 0L)
})

test_that("a column at the inclusive colmax boundary is usable", {
  set.seed(22)
  m <- matrix(rnorm(100 * 4), nrow = 100)
  m[1:29, 1] <- NA
  m[1:28, 2] <- NA

  # Regression: column 1 has 29 of 100 missing, which is not over
  # colmax = 0.29, but aborted the call with "Some columns already exceed
  # `colmax`".
  out <- sample_na_loc(
    m,
    n_cols = 1,
    n_rows = 1,
    colmax = 0.29,
    na_col_subset = 3:4
  )
  expect_identical(nrow(out[[1]]), 1L)

  # column 2 has room for exactly one more NA under the same cap
  out <- sample_na_loc(
    m,
    n_cols = 1,
    n_rows = 1,
    colmax = 0.29,
    na_col_subset = 2
  )
  expect_identical(unname(out[[1]][, "col"]), 2L)
  # and column 1 has none
  expect_error(
    sample_na_loc(m, n_cols = 1, n_rows = 1, colmax = 0.29, na_col_subset = 1),
    class = "slideimp_infeasible"
  )
})

# protected rows
test_that("rows at their rowmax budget are the ones kept observed", {
  # 6 x 10. Rows 1 to 3 sit at the rowmax = 0.5 cap (5 NAs each, in columns
  # 6 to 10), so only rows 4 to 6 can take an NA. Columns 1 to 3 must each
  # give up exactly rows 4 to 6 and keep rows 1 to 3, which hold 3 distinct
  # values. That design is valid.
  set.seed(21)
  m <- matrix(rnorm(60), nrow = 6)
  m[1:3, 6:10] <- NA
  hand <- cbind(row = rep(4:6, 3), col = rep(1:3, each = 3))
  chk <- apply_na(m, hand)
  expect_lte(max(rowSums(is.na(chk))), max_miss_count(10, 0.5))
  expect_true(all(
    apply(chk[, 1:3], 2, function(x) length(unique(x[!is.na(x)]))) >= 2L
  ))

  # Regression: the two protected rows were drawn uniformly over all six
  # observed rows, so a column only succeeded when both landed in rows 1 to 3
  # (probability 0.2), and 95 of 200 seeds failed all 100 attempts. Now the
  # rows that cannot take an NA are the protected ones, and every seed finds
  # the one valid design.
  n_ok <- 0L
  for (seed in 1:40) {
    set.seed(seed)
    out <- sample_na_loc(
      m,
      n_cols = 3,
      n_rows = 3,
      rowmax = 0.5,
      na_col_subset = 1:3
    )
    got <- out[[1]][order(out[[1]][, "col"], out[[1]][, "row"]), ]
    n_ok <- n_ok + identical(unname(got), unname(hand))
  }
  expect_identical(n_ok, 40L)
})

test_that("sampled locations satisfy every constraint on tight matrices", {
  # Random small matrices with pre-existing NAs, repeated values and tight
  # caps, so each branch of the protected-row logic is reached: no row at its
  # budget, such rows holding one value, and such rows holding two.
  n_done <- 0L
  n_refused <- 0L
  for (seed in 1:150) {
    set.seed(seed)
    m <- matrix(sample(c(0.1, 0.2, 0.3, 0.4), 8 * 6, replace = TRUE), nrow = 8)
    m[sample.int(length(m), 12)] <- NA
    rowmax <- 0.5
    colmax <- 0.75
    out <- tryCatch(
      sample_na_loc(
        m,
        n_cols = 3,
        n_rows = 2,
        rowmax = rowmax,
        colmax = colmax
      ),
      slideimp_infeasible = function(e) NULL,
      error = function(e) {
        # a matrix whose own rows are already over rowmax
        expect_match(conditionMessage(e), "already exceed", fixed = TRUE)
        NULL
      }
    )
    if (is.null(out)) {
      n_refused <- n_refused + 1L
      next
    }
    loc <- out[[1]]
    n_done <- n_done + 1L
    expect_false(anyNA(m[loc])) # observed cells only
    expect_identical(anyDuplicated(loc), 0L)
    expect_identical(as.integer(table(loc[, "col"])), c(2L, 2L, 2L))
    m2 <- apply_na(m, loc)
    expect_lte(max(rowSums(is.na(m2))), max_miss_count(ncol(m), rowmax))
    expect_lte(
      max(colSums(is.na(m2))[unique(loc[, "col"])]),
      max_miss_count(nrow(m), colmax)
    )
    for (j in unique(loc[, "col"])) {
      expect_gte(length(unique(m2[!is.na(m2[, j]), j])), 2L)
    }
  }
  # not vacuous: most fixtures are sampled, and some are refused
  expect_gt(n_done, 75L)
  expect_identical(n_done + n_refused, 150L)
})

test_that("columns are drawn uniformly from the pool", {
  # one column per draw from a pool of 10: a lazy shuffle that favoured a
  # position would show up as a skewed count
  set.seed(23)
  m <- matrix(rnorm(8 * 10), nrow = 8)
  out <- sample_na_loc(m, n_cols = 1, n_rows = 1, n_reps = 4000)
  counts <- tabulate(
    vapply(out, function(r) r[1, "col"], integer(1)),
    nbins = 10
  )
  # expected 400 each, binomial sd about 19: 5 sd either way
  expect_true(all(abs(counts - 400) < 95))

  # and without replacement inside a rep: 10 columns out of 10 is all of them
  out <- sample_na_loc(m, n_cols = 10, n_rows = 1, n_reps = 20)
  expect_true(all(vapply(
    out,
    function(r) identical(sort(unname(r[, "col"])), 1:10),
    logical(1)
  )))
  # in more than one order
  expect_gt(length(unique(lapply(out, function(r) unname(r[, "col"])))), 1L)
})

test_that("sample_na_loc() does not build an n x p logical mask", {
  n <- 200L
  p <- 2000L
  set.seed(24)
  m <- matrix(rnorm(n * p), n)
  m[sample.int(length(m), 500)] <- NA

  # Regression: `!is.na(obj)` allocated two n x p logicals (1.6 MB each
  # here) to count missing values per row and column.
  f <- tempfile(fileext = ".out")
  Rprofmem(f, threshold = 4 * n * p)
  out <- sample_na_loc(m, n_cols = 5, n_rows = 2)
  Rprofmem(NULL)
  big <- grep("^[0-9]+ :", readLines(f), value = TRUE)
  unlink(f)
  expect_length(big, 0L)
  expect_identical(nrow(out[[1]]), 10L)
})

# argument checks and messages
test_that("sample_na_loc() messages state the right reason", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input

  # Regression: "each selected column needs at least one NA"
  err <- expect_error(sample_na_loc(m, num_na = 1, n_rows = 3))
  expect_match(
    conditionMessage(err),
    "`num_na` (1) must be at least `n_rows` (3)",
    fixed = TRUE
  )
  expect_no_match(conditionMessage(err), "at least one NA", fixed = TRUE)

  # Regression: "Maximum of 5 NAs per column (from `n_rows` = 3)". The 5 is
  # `num_na` folded into one column.
  set.seed(25)
  err <- expect_error(
    sample_na_loc(matrix(rnorm(36), 6), num_na = 5, n_rows = 3),
    class = "slideimp_infeasible"
  )
  expect_match(conditionMessage(err), "puts 5 NAs in one column", fixed = TRUE)
  expect_match(conditionMessage(err), "has 6 rows", fixed = TRUE)

  # Regression: "Maximum of 1 NAs" and "only 1 rows"
  err <- expect_error(
    sample_na_loc(matrix(rnorm(5), 1), n_cols = 1, n_rows = 1),
    class = "slideimp_infeasible"
  )
  expect_match(conditionMessage(err), "puts 1 NA in one column", fixed = TRUE)
  expect_match(conditionMessage(err), "has 1 row.", fixed = TRUE)
})

test_that("sample_na_loc() rejects a non-numeric obj by name", {
  # Regression: a logical matrix failed inside col_vars(), with the same
  # text but raised from there, and a character matrix reached C++ ("Not
  # compatible with requested type").
  err <- expect_error(
    sample_na_loc(matrix(c(TRUE, FALSE), 4, 5), n_cols = 1, n_rows = 1),
    "Assertion on 'obj' failed"
  )
  expect_identical(as.character(conditionCall(err)[[1]]), "sample_na_loc")
  expect_error(
    sample_na_loc(matrix(letters[1:20], 4), n_cols = 1, n_rows = 1),
    "Assertion on 'obj' failed"
  )
})

test_that("an integerish na_col_subset just under its integer is not truncated", {
  set.seed(26)
  m <- matrix(rnorm(10 * 30), nrow = 10)
  p <- 0.29 * 100
  expect_lt(p, 29)
  expect_true(checkmate::test_integerish(p))
  out <- sample_na_loc(m, n_cols = 1, n_rows = 1, na_col_subset = p, n_reps = 5)
  expect_identical(
    unique(vapply(out, function(r) r[1, "col"], integer(1))),
    29L
  )
})

# resolve_na_loc ----
test_that("resolve_na_loc() rounds near-integer positions", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:90), nrow = 30)

  # Regression: 0.29 * 100 is 28.999999999999996. It passes the integerish
  # check, and as.integer() truncated it to 28: the wrong cell was masked.
  p <- 0.29 * 100
  expect_identical(as.integer(p), 28L)
  expect_identical(unname(rna(obj, p)[[1]][1, ]), c(29L, 1L))
  expect_identical(
    unname(rna(obj, matrix(c(p, 2), ncol = 2))[[1]][1, ]),
    c(29L, 2L)
  )
})

test_that("resolve_na_loc() rejects repeated and missing matrix positions", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:50), nrow = 10)

  # Regression: the matrix form took (1, 1) twice, and tune_imp() scored the
  # cell twice. The vector form already refused the same repeat.
  dup <- matrix(c(1L, 1L, 2L, 1L, 1L, 2L), ncol = 2)
  expect_error(rna(obj, list(dup)), "positions must be unique")
  expect_error(rna(obj, list(c(1L, 1L, 12L))), "duplicated")
  # a repeat that only appears after rounding
  expect_error(rna(obj, list(c(0.29 * 100, 29))), "positions must be unique")
  fill0 <- function(obj) {
    obj[is.na(obj)] <- 0
    obj
  }
  expect_error(
    suppressMessages(tune_imp(
      obj,
      .f = fill0,
      na_loc = list(dup),
      .progress = FALSE
    )),
    "positions must be unique"
  )

  # Regression: an NA inside a matrix read as an out-of-bounds index
  # ("'na_loc[[1]] row indices' failed: Must be TRUE")
  expect_error(rna(obj, list(matrix(c(1, NA, 1, 1), 2))), "missing values")
})

test_that("resolve_na_loc() reads a matrix named col, row by name", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:50), nrow = 10)

  # Regression: cbind(col = 2, row = 5) was read by position, as row 2,
  # column 5.
  expect_identical(
    rna(obj, list(cbind(col = 2L, row = 5L)))[[1]],
    cbind(row = 5L, col = 2L)
  )
  # sample_na_loc() output, and unnamed matrices, are read as before
  expect_identical(
    rna(obj, list(cbind(row = 5L, col = 2L)))[[1]],
    cbind(row = 5L, col = 2L)
  )
  expect_identical(
    rna(obj, list(cbind(5L, 2L)))[[1]],
    cbind(row = 5L, col = 2L)
  )
})

test_that("resolve_na_loc() argument errors are cli conditions", {
  rna <- function(obj, na_loc) {
    resolve_na_loc(
      obj = obj,
      na_loc = na_loc,
      n_reps = 1L,
      num_na = NULL,
      n_cols = NULL,
      n_rows = 2L,
      rowmax = 0.9,
      colmax = 0.9,
      na_col_subset = NULL,
      max_attempts = 100
    )
  }
  obj <- matrix(as.numeric(1:50), nrow = 10)

  # Regression: both were base stop() calls
  err <- expect_error(rna(obj, "a"), "must be a 2-col matrix")
  expect_s3_class(err, "rlang_error")
  err <- expect_error(rna(obj, list(1:2, 1:3)), "same number of NA positions")
  expect_s3_class(err, "rlang_error")
})

# tune imp ----
test_that("tune_imp works", {
  slide_imp_par <- data.frame(
    window_size = c(100, 100),
    k = c(5, 10),
    overlap_size = c(10, 10),
    min_window_n = 20,
    method = "euclidean",
    post_imp = FALSE
  )
  set.seed(1234)
  obj <- sim_mat(50, 1000, perc_col_na = 0.5)$input
  expect_true(anyNA(obj))
  location <- 1:ncol(obj)
  # Check `slide_imp`
  expect_no_error({
    slide_imp_imp_res <- tune_imp(
      obj,
      slide_imp_par,
      .f = "slide_imp",
      location = location,
      n_reps = 1,
      num_na = 200
    )
  })

  # `slide_imp` requires parameters
  expect_error(
    {
      slide_imp_imp_res <- tune_imp(
        obj,
        .f = "slide_imp",
        location = location,
        n_reps = 1,
        num_na = 200
      )
    },
    regexp = "requires"
  )

  expect_true(
    all(
      vapply(
        slide_imp_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )

  # Check `knn_imp`
  knn_imp_par <- data.frame(
    k = c(5, 10),
    method = "euclidean",
    post_imp = TRUE
  )
  expect_no_error({
    knn_imp_res <- tune_imp(
      obj,
      knn_imp_par,
      .f = "knn_imp",
      n_reps = 1,
      num_na = 100
    )
  })

  expect_true(
    all(
      vapply(
        knn_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )

  # Check `pca_imp`
  pca_imp_par <- data.frame(ncp = 2, miniter = 2)
  expect_no_error({
    pca_imp_res <- tune_imp(
      obj,
      pca_imp_par,
      .f = "pca_imp",
      n_reps = 1,
      num_na = 100
    )
  })

  expect_true(
    all(
      vapply(
        pca_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )

  # Check custom function
  f1 <- function() {}
  custom_fun <- function(obj, value) {
    obj[is.na(obj)] <- value
    f1()
    return(obj)
  }
  custom_par <- data.frame(
    value = c(0, 1)
  )
  expect_no_error({
    custom_imp_res <- tune_imp(
      obj,
      custom_par,
      n_reps = 1,
      num_na = 100,
      .f = custom_fun
    )
  })

  expect_true(
    all(
      vapply(
        custom_imp_res$result,
        \(x) {
          class(x$estimate)
        },
        character(1)
      ) ==
        "numeric"
    )
  )
})

test_that("tune_imp runs the parallel (mirai) branch", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1234)
  obj <- sim_mat(50, 1000, perc_col_na = 0.5)$input
  knn_imp_par <- data.frame(k = c(5, 10), method = "euclidean", post_imp = TRUE)
  mirai::daemons(2, seed = 1234)
  withr::defer(mirai::daemons(0))
  expect_message(
    res <- tune_imp(obj, knn_imp_par, .f = "knn_imp", n_reps = 2, num_na = 100),
    regexp = "Running mode.*mirai"
  )
  expect_true(all(is.na(res$error)))
  expect_true(
    all(vapply(res$result, \(x) class(x$estimate), character(1)) == "numeric")
  )
})

test_that("tune_imp forces cores = 1 under active daemons", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1234)
  obj <- sim_mat(50, 1000, perc_col_na = 0.5)$input
  knn_imp_par <- data.frame(k = 5, method = "euclidean", post_imp = TRUE)
  mirai::daemons(2, seed = 1234)
  withr::defer(mirai::daemons(0))
  expect_message(
    output <- tune_imp(
      obj,
      knn_imp_par,
      .f = "knn_imp",
      n_reps = 1,
      num_na = 100,
      cores = 2
    ),
    regexp = "nested parallelism"
  )
})

test_that("tune_imp works when n_reps is a list of NA locations", {
  # Create a complete matrix (no NAs) for testing
  obj <- sim_mat(50, 200)$input
  obj[is.na(obj)] <- 0 # Fill any existing NAs

  # Create predefined NA location sets
  # Each set has 10 locations, all within matrix bounds
  set.seed(42)
  na_loc_list <- list(
    sample(1:length(obj), 10, replace = FALSE),
    sample(1:length(obj), 10, replace = FALSE),
    sample(1:length(obj), 10, replace = FALSE)
  )

  # Test with slide_imp
  slide_imp_par <- data.frame(
    window_size = 100,
    k = 5,
    overlap_size = 10,
    method = "euclidean",
    min_window_n = 10,
    post_imp = FALSE
  )

  location <- 1:ncol(obj)
  expect_no_error({
    slide_imp_res <- tune_imp(
      location = location,
      obj,
      slide_imp_par,
      .f = "slide_imp",
      na_loc = na_loc_list, # Using list instead of integer
    )
  })

  # Check that we get 3 results (one for each NA location set)
  expect_equal(nrow(slide_imp_res), 3)

  # Check that each result has the correct number of estimates (10 each)
  expect_true(
    all(vapply(slide_imp_res$result, function(x) nrow(x) == 10, logical(1)))
  )

  # Verify the truth values match the original matrix values at those locations
  for (i in 1:3) {
    truth_values <- slide_imp_res$result[[i]]$truth
    expected_truth <- obj[na_loc_list[[i]]]
    expect_equal(truth_values, expected_truth)
  }

  # Test with knn_imp
  knn_imp_par <- data.frame(
    k = c(5, 10),
    method = "euclidean",
    post_imp = FALSE
  )

  expect_no_error({
    knn_imp_res <- tune_imp(
      obj,
      knn_imp_par,
      .f = "knn_imp",
      na_loc = na_loc_list
    )
  })

  # Should have 2 parameters x 3 repetitions = 6 rows
  expect_equal(nrow(knn_imp_res), 6)

  # Check that results contain numeric estimates
  expect_true(
    all(vapply(
      knn_imp_res$result,
      function(x) {
        is.numeric(x$estimate) && is.numeric(x$truth)
      },
      logical(1)
    ))
  )

  # Test with custom function
  custom_fun <- function(obj, value) {
    obj[is.na(obj)] <- value
    return(obj)
  }

  custom_par <- data.frame(value = c(0.5, 1.5))

  expect_no_error({
    custom_res <- tune_imp(
      obj,
      custom_par,
      na_loc = na_loc_list,
      .f = custom_fun
    )
  })

  # Should have 2 parameters x 3 repetitions = 6 rows
  expect_equal(nrow(custom_res), 6)

  # Verify custom function fills with the specified values
  for (i in 1:nrow(custom_res)) {
    expected_value <- custom_res$value[i]
    estimates <- custom_res$result[[i]]$estimate
    expect_true(all(estimates == expected_value))
  }

  # Test with different length NA location sets
  varied_na_locs <- list(
    sample(1:length(obj), 5, replace = FALSE),
    sample(1:length(obj), 5, replace = FALSE)
  )

  location <- 1:ncol(obj)
  expect_no_error({
    varied_res <- tune_imp(
      obj,
      location = location,
      slide_imp_par,
      .f = "slide_imp",
      na_loc = varied_na_locs
    )
  })

  expect_equal(nrow(varied_res), 2)
  expect_equal(nrow(varied_res$result[[1]]), 5)
  expect_equal(nrow(varied_res$result[[2]]), 5)
})

test_that("tune_imp correctly uses provided NA locations from list", {
  # Create a simple matrix for easier verification
  set.seed(123)
  obj <- matrix(1:100, nrow = 10, ncol = 10)

  # Define specific NA locations
  na_locations <- list(
    c(1, 11, 21), # First column positions
    c(10, 20, 30), # Last position of first 3 rows
    c(50, 60, 70) # Middle positions
  )

  simple_imp <- function(obj, fill_value) {
    obj[is.na(obj)] <- fill_value
    return(obj)
  }

  params <- data.frame(fill_value = 42)

  result <- tune_imp(
    obj,
    params,
    na_loc = na_locations,
    .f = simple_imp
  )

  # Verify each repetition used the correct NA locations
  for (i in 1:3) {
    res <- result$result[[i]]

    # Check truth values match original matrix at specified locations
    expected_truth <- obj[na_locations[[i]]]
    expect_equal(res$truth, expected_truth)

    # Check all estimates are the fill value
    expect_true(all(res$estimate == 42))

    # Check we have the right number of values
    expect_equal(length(res$truth), length(na_locations[[i]]))
  }
})

test_that("tune_imp handles mixed linear and 2D positions in list", {
  set.seed(789)
  obj <- matrix(1:100, nrow = 10, ncol = 10)

  # mix of linear and 2D positions
  na_locations_mixed <- list(
    c(1, 11, 21), # linear
    matrix(c(10, 10, 10, 1, 2, 3), ncol = 2), # 2D, row 10, column 1, 2, 3
    c(45, 55, 65) # linear
  )

  simple_imp <- function(obj, fill_value) {
    obj[is.na(obj)] <- fill_value
    return(obj)
  }

  params <- data.frame(fill_value = 67)

  result <- tune_imp(
    obj,
    params,
    na_loc = na_locations_mixed,
    .f = simple_imp
  )

  expected_linear <- list(
    c(1, 11, 21),
    c(10, 20, 30),
    c(45, 55, 65)
  )

  for (i in 1:3) {
    res <- result$result[[i]]
    expected_truth <- obj[expected_linear[[i]]]
    expect_equal(res$truth, expected_truth)
    expect_true(all(res$estimate == 67))
  }
})

test_that("compute_metrics works with slideimp_tune and data.frame", {
  set.seed(123)
  obj <- matrix(1:100, nrow = 10, ncol = 10)

  simple_imp <- function(obj, mu) {
    miss <- is.na(obj)
    obj[miss] <- stats::rnorm(n = sum(miss), mean = mu)
    return(obj)
  }

  params <- data.frame(mu = 42)
  result_tune <- tune_imp(
    obj,
    params,
    n_reps = 2,
    num_na = 10,
    .f = simple_imp
  )

  # slideimp_tune object
  out_tune <- compute_metrics(result_tune)
  expect_s3_class(out_tune, "data.frame")
  expect_true(all(
    c(".metric", ".estimator", ".estimate", "n", "n_miss") %in% names(out_tune)
  ))
  expect_equal(sort(unique(out_tune$.metric)), c("mae", "rmse"))

  # Plain data.frame
  result_df <- as.data.frame(result_tune)
  class(result_df) <- "data.frame"
  out_df <- compute_metrics(result_df)
  expect_s3_class(out_df, "data.frame")
  expect_equal(out_df$.estimate, out_tune$.estimate)
})

test_that("compute_metrics correctly computes n and n_miss with NA estimates", {
  set.seed(456)
  obj <- matrix(1:100, nrow = 10, ncol = 10)

  simple_imp <- function(obj, mu) {
    miss <- is.na(obj)
    obj[miss] <- rnorm(n = sum(miss), mean = mu)
    return(obj)
  }

  params <- data.frame(mu = 42)
  result <- tune_imp(
    obj,
    params,
    n_reps = 2,
    num_na = 10,
    .f = simple_imp
  )

  # No NAs case: all estimates should be present
  out_clean <- compute_metrics(result)
  expect_true(all(out_clean$n == 10))
  expect_true(all(out_clean$n_miss == 0))

  # Inject NAs into the estimate column of each result element
  result$result[[1]]$estimate[c(1, 3)] <- NA
  result$result[[2]]$estimate[c(2, 5, 7)] <- NA

  out_na <- compute_metrics(
    result,
    metrics = c("mae", "rmse", "mape", "bias", "rsq", "rsq_trad")
  )

  # Rep 1: 10 rows, 2 missing
  rows_rep1 <- out_na[out_na$rep_id == 1, ]
  expect_true(all(rows_rep1$n == 10))
  expect_true(all(rows_rep1$n_miss == 2))

  # Rep 2: 10 rows, 3 missing
  rows_rep2 <- out_na[out_na$rep_id == 2, ]
  expect_true(all(rows_rep2$n == 10))
  expect_true(all(rows_rep2$n_miss == 3))

  # n and n_miss are consistent across metrics within the same rep
  for (r in unique(out_na$rep_id)) {
    subset <- out_na[out_na$rep_id == r, ]
    expect_length(unique(subset$n), 1)
    expect_length(unique(subset$n_miss), 1)
  }
})

test_that("compute_metrics.data.frame errors without required columns", {
  bad_df <- data.frame(x = 1:3)
  expect_error(compute_metrics(bad_df), "result")

  bad_result <- data.frame(result = I(list(data.frame(a = 1, b = 2))))
  expect_error(compute_metrics(bad_result), "truth.*estimate")
})

test_that("tune_imp works with custom function and list-column parameters", {
  set.seed(42)
  obj <- matrix(rnorm(200), nrow = 10, ncol = 20)

  # custom function that takes a vector of weights per column and fills NAs
  # with a weighted column mean
  weighted_fill <- function(obj, weights) {
    stopifnot(length(weights) == ncol(obj))
    for (j in seq_len(ncol(obj))) {
      col_mean <- mean(obj[, j], na.rm = TRUE)
      obj[is.na(obj[, j]), j] <- col_mean * weights[j]
    }
    return(obj)
  }

  # parameters with a list column: each row holds a different weight vector
  custom_par <- data.frame(
    weights = I(list(
      rep(1, 20),
      rep(0.5, 20),
      seq(0.1, 2, length.out = 20)
    ))
  )

  expect_no_error({
    res <- tune_imp(
      obj,
      custom_par,
      .f = weighted_fill,
      n_reps = 2,
      num_na = 15
    )
  })

  # Should have 3 param sets * 2 reps = 6 rows
  expect_equal(nrow(res), 6)
  expect_true("result" %in% names(res))
  expect_true(
    all(
      vapply(
        res$result,
        \(x) {
          is.numeric(x$estimate) && nrow(x) > 0
        },
        logical(1)
      )
    )
  )
  # The weights list column should be preserved in the output
  expect_true("weights" %in% names(res))
  expect_true(is.list(res$weights))
})

test_that("tune_imp works with custom function and NULL parameters", {
  set.seed(99)
  obj <- matrix(rnorm(150), nrow = 10, ncol = 15)

  # a function with only `obj` - fills NAs with 0
  zero_fill <- function(obj) {
    obj[is.na(obj)] <- 0
    return(obj)
  }

  expect_no_error({
    res <- tune_imp(
      obj,
      parameters = NULL,
      .f = zero_fill,
      n_reps = 3,
      num_na = 10
    )
  })

  # 1 param set * 3 reps = 3 rows
  expect_equal(nrow(res), 3)
  expect_true("result" %in% names(res))
  # placeholder column should be stripped
  expect_false(".placeholder" %in% names(res))
  expect_true(
    all(
      vapply(
        res$result,
        \(x) {
          is.numeric(x$estimate) && is.numeric(x$truth) && nrow(x) > 0
        },
        logical(1)
      )
    )
  )
})

test_that("tune_imp with NULL parameters and a function that has defaults", {
  set.seed(7)
  obj <- matrix(rnorm(100), nrow = 5, ncol = 20)

  fill_with_default <- function(obj, value = -999, scale = 1.0) {
    obj[is.na(obj)] <- value * scale
    return(obj)
  }

  # NULL parameters should run the function using its defaults
  expect_no_error({
    res_null <- tune_imp(
      obj,
      parameters = NULL,
      .f = fill_with_default,
      n_reps = 1,
      num_na = 10
    )
  })

  expect_equal(nrow(res_null), 1)
  # all imputed values should be -999 (the defaults)
  expect_true(all(res_null$result[[1]]$estimate == -999))

  # compare with explicit parameters to make sure NULL truly uses defaults
  explicit_par <- data.frame(value = -999, scale = 1.0)
  expect_no_error({
    res_explicit <- tune_imp(
      obj,
      parameters = explicit_par,
      .f = fill_with_default,
      n_reps = 1,
      num_na = 10
    )
  })

  expect_equal(res_null$result[[1]]$estimate, res_explicit$result[[1]]$estimate)
})

test_that("serial branch reuses `pre` correctly across param_sets within a rep", {
  set.seed(1234)
  obj <- matrix(rnorm(50 * 100), nrow = 50, ncol = 100) # clean, no pre-existing NAs

  seen <- list()
  capture_fun <- function(obj, tag) {
    seen[[length(seen) + 1L]] <<- which(is.na(obj), arr.ind = TRUE)
    obj[is.na(obj)] <- 0
    obj
  }

  params <- data.frame(tag = c("a", "b", "c")) # 3 param sets

  res <- tune_imp(
    obj,
    params,
    .f = capture_fun,
    n_reps = 2,
    num_na = 20,
    .progress = FALSE
  )

  # 1 probe call + 2 reps * 3 param_sets = 7
  expect_length(seen, 7L)

  rep1 <- seen[2:4]
  rep2 <- seen[5:7]

  # within-rep: NA mask identical across param_sets (pre is reused correctly)
  expect_equal(rep1[[1]], rep1[[2]])
  expect_equal(rep1[[2]], rep1[[3]])
  expect_equal(rep2[[1]], rep2[[2]])
  expect_equal(rep2[[2]], rep2[[3]])

  # across-rep: masks differ (pre is rebuilt with the new rep's NAs)
  expect_false(identical(rep1[[1]], rep2[[1]]))

  # truth_list is indexed by the correct rep_id.
  # sort both sides because na_loc ordering (from C++) != which() column-major order.
  for (i in seq_len(nrow(res))) {
    mask <- if (res$rep_id[i] == 1L) rep1[[1]] else rep2[[1]]
    expect_equal(sort(res$result[[i]]$truth), sort(obj[mask]))
  }
})

test_that("slide_imp flank mode is tuned correctly", {
  set.seed(1234)
  ncols <- 100
  obj <- sim_mat(n = 50, p = ncols, perc_col_na = 1)$input
  # we test imputing subset feature 5, 20, and 95
  testing_features <- paste0("feature", c(5, 20, 95))
  na_loc <- sample_na_loc(
    obj = obj,
    n_cols = 3,
    na_col_subset = testing_features,
    n_rows = 2,
    n_reps = 1
  )

  # we manually use slide imp flank mode by manual NA injection into obj
  slide_imp_pre <- obj
  truth <- slide_imp_pre[na_loc[[1]]]
  slide_imp_pre[na_loc[[1]]] <- NA
  location <- seq_len(ncols)

  ## slide_imp on the amputed object. Windows are
  #  # slideimp table: 3 x 5
  # start end window_n target  subset_local
  #     1  10       10      5 <integer [1]>
  #    15  25       11     20 <integer [1]>
  #    90 100       11     95 <integer [1]>

  slide_imp_knn <- slide_imp(
    obj = slide_imp_pre,
    location = location,
    k = 2,
    subset = testing_features,
    min_window_n = 3,
    window_size = 5,
    overlap_size = 0,
    flank = TRUE
  )

  slide_imp_pca <- slide_imp(
    obj = slide_imp_pre,
    location = location,
    ncp = 2,
    subset = testing_features,
    min_window_n = 3,
    window_size = 5,
    overlap_size = 0,
    flank = TRUE
  )

  ## now we setup tune_imp
  params <- data.frame(
    min_window_n = 3,
    window_size = 5,
    overlap_size = 0,
    flank = TRUE
  )
  knn_params <- params
  knn_params$k <- 2
  pca_params <- params
  pca_params$ncp <- 2

  tune_imp_knn <- tune_imp(
    obj,
    parameters = knn_params,
    na_loc = na_loc,
    .f = "slide_imp",
    location = location
  )

  tune_imp_pca <- tune_imp(
    obj,
    parameters = pca_params,
    na_loc = na_loc,
    .f = "slide_imp",
    location = location
  )

  expect_equal(
    tune_imp_knn$result[[1]]$estimate,
    slide_imp_knn[na_loc[[1]]]
  )
  expect_equal(
    tune_imp_knn$result[[1]]$truth,
    truth
  )
  expect_equal(
    tune_imp_pca$result[[1]]$estimate,
    slide_imp_pca[na_loc[[1]]]
  )
  expect_equal(
    tune_imp_pca$result[[1]]$truth,
    truth
  )
})

test_that("supplied na_loc positions must point at observed cells", {
  set.seed(1234)
  obj <- sim_mat(20, 20, perc_total_na = 0.2)$input
  obj[5, 5] <- NA
  obj[6, 6] <- NaN
  observed <- which(!is.na(obj), arr.ind = TRUE)

  fill0 <- function(obj, dummy = 1) {
    obj[is.na(obj)] <- 0
    obj
  }
  params <- data.frame(dummy = 1)

  # a position sitting on an existing hole has no truth behind it
  on_hole <- rbind(observed[1:3, ], c(5L, 5L))
  expect_error(
    tune_imp(
      obj,
      params,
      .f = fill0,
      na_loc = list(on_hole),
      .progress = FALSE
    ),
    "observed cells"
  )

  # NaN is a hole too
  on_nan <- rbind(observed[1:3, ], c(6L, 6L))
  expect_error(
    tune_imp(obj, params, .f = fill0, na_loc = list(on_nan), .progress = FALSE),
    "observed cells"
  )

  # linear positions go through the same check
  expect_error(
    tune_imp(
      obj,
      params,
      .f = fill0,
      na_loc = list(c(which(is.na(obj))[1], which(!is.na(obj))[1])),
      .progress = FALSE
    ),
    "observed cells"
  )

  # and fully observed positions are scored in full, with no silent shrink
  clean <- observed[1:6, ]
  res <- tune_imp(
    obj,
    params,
    .f = fill0,
    na_loc = list(clean),
    .progress = FALSE
  )
  d <- res$result[[1]]
  expect_equal(nrow(d), 6L)
  expect_false(anyNA(d$truth))
  expect_equal(sum(is.finite(d$truth - d$estimate)), 6L)
})

test_that("sample_na_loc positions are always observed cells", {
  set.seed(1234)
  obj <- sim_mat(20, 20, perc_total_na = 0.2)$input
  locs <- sample_na_loc(obj, n_cols = 5, n_rows = 2, n_reps = 3)

  for (pos in locs) {
    expect_false(anyNA(obj[pos]))
  }
})

test_that("sample_na_loc returns integer index matrices", {
  # sample_each_rep_cpp() hands back an arma::umat, which under 64-bit arma
  # indices arrives as doubles. Indices are bounded by dim(), so they are
  # coerced back to integer to match the documented return type.
  set.seed(1234)
  obj <- sim_mat(20, 30, perc_total_na = 0.2)$input
  locs <- sample_na_loc(obj, n_cols = 5, n_rows = 2, n_reps = 3)

  for (pos in locs) {
    expect_type(pos, "integer")
    expect_equal(colnames(pos), c("row", "col"))
    expect_length(obj[pos], nrow(pos))
  }
})

# tune_imp() driver ----
test_that("a grid built with expand.grid() tunes a built-in method", {
  # expand.grid() makes a character column a factor, and match.arg() in the
  # methods rejects a factor: every iteration used to fail
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  grid <- expand.grid(k = c(3, 5), method = c("euclidean", "manhattan"))
  expect_true(is.factor(grid$method))

  expect_no_warning(
    res <- suppressMessages(
      tune_imp(obj, grid, .f = "knn_imp", num_na = 10, .progress = FALSE)
    )
  )
  expect_identical(nrow(res), 4L)
  expect_true(all(is.na(res$error)))
  expect_true(all(vapply(res$result, nrow, integer(1)) == 10L))
  # the caller's column comes back as it was given
  expect_true(is.factor(res$method))
})

test_that("the output of a built-in method holds only the caller's parameters", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  res <- suppressMessages(
    tune_imp(
      obj,
      data.frame(k = 3),
      .f = "knn_imp",
      num_na = 10,
      .progress = FALSE
    )
  )
  expect_identical(names(res), c("k", "param_set", "rep_id", "result", "error"))

  # a `.progress` column is dropped, with a note, and the call still runs
  expect_message(
    res2 <- tune_imp(
      obj,
      data.frame(k = c(3, 3), .progress = c(TRUE, FALSE)),
      .f = "knn_imp",
      num_na = 10,
      .progress = FALSE
    ),
    "suppressed"
  )
  expect_identical(names(res2), names(res))
  expect_identical(nrow(res2), 1L)
  expect_true(all(is.na(res2$error)))
})

test_that("a misplaced k or ncp gets the message that names the other method", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  run <- function(parameters, .f) {
    suppressMessages(
      tune_imp(obj, parameters, .f = .f, num_na = 10, .progress = FALSE)
    )
  }
  expect_error(
    run(data.frame(k = 3, ncp = 2), "knn_imp"),
    "Did you mean .*pca_imp"
  )
  expect_error(
    run(data.frame(ncp = 2, k = 3), "pca_imp"),
    "Did you mean .*knn_imp"
  )
  # any other unknown name still gets the general message
  expect_error(
    run(data.frame(k = 3, subset = 1), "knn_imp"),
    "Unknown parameter"
  )
})

test_that("every tunable slide_imp() argument is accepted in parameters", {
  # what tune_imp() supplies itself, or what breaks the matrix return
  own <- c("obj", "location", "subset", "cores", "dry_run")
  tunable <- setdiff(names(formals(slide_imp)), own)
  expect_gt(length(tunable), 5L)
  expect_setequal(slide_tune_params(), tunable)
  expect_true(all(.param_registry$knn$allowed %in% slide_tune_params()))
  expect_true(all(.param_registry$pca$allowed %in% slide_tune_params()))
})

test_that("parameters cannot use a name the tuning output needs", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  f_n <- function(obj, n) {
    obj[is.na(obj)] <- n
    obj
  }
  err <- expect_error(
    tune_imp(
      obj,
      data.frame(n = c(3, 7)),
      .f = f_n,
      num_na = 10,
      .progress = FALSE
    ),
    "Name in use: \"n\""
  )
  expect_s3_class(err, "rlang_error")

  f_two <- function(obj, error, result) {
    obj[is.na(obj)] <- 0
    obj
  }
  expect_error(
    tune_imp(
      obj,
      data.frame(error = 1, result = 2),
      .f = f_two,
      num_na = 10,
      .progress = FALSE
    ),
    "Names in use: \"error\" and \"result\""
  )
})

test_that("a custom .f return is checked on every parameter set", {
  # the probe runs parameter set 1 only. A later set that returned another
  # shape was scored: a matrix index still finds cells in a taller matrix.
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  f_shape <- function(obj, mode) {
    obj[is.na(obj)] <- if (mode == 4) Inf else 0.5
    if (mode == 2) {
      return(rbind(obj, obj))
    }
    if (mode == 3) {
      return(as.vector(obj))
    }
    obj
  }
  expect_warning(
    res <- suppressMessages(
      tune_imp(
        obj,
        data.frame(mode = c(1, 2, 3, 4)),
        .f = f_shape,
        num_na = 10,
        .progress = FALSE
      )
    ),
    "may have failed"
  )
  expect_identical(res$mode, c(1, 2, 3, 4))
  expect_true(is.na(res$error[1]))
  expect_identical(nrow(res$result[[1]]), 10L)
  expect_match(res$error[2], "returned a 40 x 30 matrix, but `obj` is 20 x 30")
  expect_match(res$error[3], "did not return a numeric matrix")
  expect_match(res$error[4], "`Inf` for a masked cell")
  expect_identical(vapply(res$result[2:4], nrow, integer(1)), c(0L, 0L, 0L))
})

test_that("cores and metrics are validated by name", {
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  expect_error(
    suppressMessages(tune_imp(
      obj,
      data.frame(k = 3),
      .f = "knn_imp",
      num_na = 10,
      cores = NA,
      .progress = FALSE
    )),
    "'cores'"
  )
  res <- suppressMessages(
    tune_imp(
      obj,
      data.frame(k = 3),
      .f = "knn_imp",
      num_na = 10,
      .progress = FALSE
    )
  )
  expect_error(compute_metrics(res, metrics = character(0)), "'metrics'")
  expect_error(compute_metrics(res, metrics = c("mae", NA)), "'metrics'")
})

test_that("compute_metrics() errors are cli conditions", {
  err <- expect_error(compute_metrics(data.frame(x = 1)), "result")
  expect_s3_class(err, "rlang_error")

  bad_result <- data.frame(result = I(list(data.frame(a = 1))))
  err <- expect_error(compute_metrics(bad_result), "truth.*estimate")
  expect_s3_class(err, "rlang_error")

  ok <- data.frame(result = I(list(data.frame(truth = 1:3, estimate = 3:1))))
  err <- expect_error(
    compute_metrics(ok, metrics = c("mae", "nope", "zzz")),
    "Unknown metrics in .metrics.: nope, zzz"
  )
  expect_s3_class(err, "rlang_error")
})

test_that("mape does not drop a cell whose true value is 0", {
  # true value 0, estimate not 0: the ratio is Inf
  expect_identical(calc_mape(c(0, 1), c(0.5, 1)), Inf)
  # true value 0, estimate 0: 0 / 0. Dropping that cell gave 50, the error of
  # the other cell alone.
  expect_true(is.nan(calc_mape(c(0, 2), c(0, 1))))
  # a cell with no estimate is still dropped
  expect_equal(calc_mape(c(2, 4, 1), c(1, NA, 1)), 25)
  # no zero: as before
  expect_equal(calc_mape(c(2, 4), c(1, 5)), 37.5)
})

test_that("rsq is NA, without a warning, against a constant", {
  expect_no_warning(const_est <- calc_rsq(c(0.1, 0.5, 0.9), c(0.3, 0.3, 0.3)))
  expect_identical(const_est, NA_real_)
  expect_no_warning(const_truth <- calc_rsq(c(0.3, 0.3, 0.3), c(0.1, 0.5, 0.9)))
  expect_identical(const_truth, NA_real_)
  # and the usual value otherwise
  expect_equal(
    calc_rsq(c(1, 2, 3, 5), c(2, 4, 7, 8)),
    cor(c(1, 2, 3, 5), c(2, 4, 7, 8))^2
  )
})

test_that("mirai_daemons_set() answers without raising", {
  skip_if(mirai::daemons_set(), "daemons are set in this session")
  expect_false(mirai_daemons_set())
})

test_that("a task whose daemon stops is recorded as a failed iteration", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  # quit() stands in for a crash in compiled code or an out-of-memory kill.
  # The probe in the main session runs parameter set 1 only.
  die_on_two <- function(obj, v) {
    if (v == 2) {
      quit(save = "no", status = 1, runLast = FALSE)
    }
    obj[is.na(obj)] <- v
    obj
  }
  mirai::daemons(2)
  withr::defer(mirai::daemons(0))
  expect_true(mirai_daemons_set())

  expect_warning(
    res <- suppressMessages(
      tune_imp(
        obj,
        data.frame(v = c(1, 2)),
        .f = die_on_two,
        num_na = 10,
        .progress = FALSE
      )
    ),
    "1 iteration did not return from the mirai daemons"
  )
  expect_true(is.na(res$error[1]))
  expect_match(res$error[2], "The mirai task failed")
  expect_s3_class(res$result[[2]], "data.frame")
  expect_identical(nrow(res$result[[2]]), 0L)
  # the metrics no longer stop on that row
  cm <- compute_metrics(res)
  expect_identical(cm$n[cm$v == 2], c(0L, 0L))
  expect_identical(cm$n[cm$v == 1], c(10L, 10L))
})

test_that("mirai_error_text() reads both kinds of failed task", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  mirai::daemons(2)
  withr::defer(mirai::daemons(0))
  res <- mirai::mirai_map(1:2, function(i) {
    if (i == 1) {
      stop("boom in the task")
    }
    quit(save = "no", status = 1, runLast = FALSE)
  })[]
  expect_true(all(vapply(res, mirai::is_error_value, logical(1))))
  expect_identical(mirai_error_text(res[[1]]), "boom in the task")
  # a stopped daemon leaves an integer code, which has no conditionMessage()
  expect_false(inherits(res[[2]], "miraiError"))
  expect_match(mirai_error_text(res[[2]]), "errorValue")
})

test_that("under mirai a custom .f finds the packages attached on the daemon", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  mirai::daemons(1)
  withr::defer(mirai::daemons(0))

  # median() is in stats, not in base
  fill_median <- function(obj) {
    obj[is.na(obj)] <- median(obj, na.rm = TRUE)
    obj
  }
  res <- suppressMessages(
    tune_imp(
      obj,
      NULL,
      .f = fill_median,
      num_na = 10,
      n_reps = 2,
      .progress = FALSE
    )
  )
  expect_identical(nrow(res), 2L)
  expect_true(all(is.na(res$error)))

  # a variable of the defining session does not travel, and the error says so
  make_fill <- function(shift) {
    force(shift)
    function(obj) {
      obj[is.na(obj)] <- shift
      obj
    }
  }
  expect_warning(
    res <- suppressMessages(
      tune_imp(
        obj,
        NULL,
        .f = make_fill(0.5),
        num_na = 10,
        n_reps = 2,
        .progress = FALSE
      )
    ),
    "may have failed"
  )
  expect_match(res$error, "'shift' not found", all = TRUE)
})

test_that("under mirai a wrong return on a later set is recorded", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  f_shape <- function(obj, mode) {
    obj[is.na(obj)] <- 0.5
    if (mode == 2) {
      return(rbind(obj, obj))
    }
    obj
  }
  mirai::daemons(1)
  withr::defer(mirai::daemons(0))
  expect_warning(
    res <- suppressMessages(
      tune_imp(
        obj,
        data.frame(mode = c(1, 2)),
        .f = f_shape,
        num_na = 10,
        .progress = FALSE
      )
    ),
    "may have failed"
  )
  expect_true(is.na(res$error[1]))
  expect_match(res$error[2], "returned a 40 x 30 matrix")
  expect_identical(nrow(res$result[[2]]), 0L)
})

test_that(".progress = FALSE draws no progress bar on the mirai path", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  slow_fill <- function(obj) {
    Sys.sleep(0.3)
    obj[is.na(obj)] <- 0
    obj
  }
  withr::local_options(cli.progress_show_after = 0)
  mirai::daemons(1)
  withr::defer(mirai::daemons(0))

  # everything the call prints that is not one of its own status lines
  drawn <- function(flag) {
    txt <- utils::capture.output(
      out <- utils::capture.output(
        res <- tune_imp(
          obj,
          NULL,
          .f = slow_fill,
          num_na = 10,
          n_reps = 3,
          .progress = flag
        ),
        type = "output"
      ),
      type = "message"
    )
    txt <- c(out, txt)
    status <- "Tuning|Step [12]/2|Running mode|pin_blas|^\\s*$"
    length(txt[!grepl(status, txt)])
  }
  expect_identical(drawn(FALSE), 0L)
  # the detector can see a bar: without this the line above cannot fail
  expect_gt(drawn(TRUE), 0L)
  expect_identical(drawn(FALSE), 0L)
})
