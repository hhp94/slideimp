skip_on_cran()

# sample_na_loc ----
test_that("returns a list of length n_reps, each a 2-col row/col matrix", {
  m <- sim_mat(20, 10, perc_total_na = 0)$input
  out <- sample_na_loc(m, n_cols = 3, n_rows = 2, n_reps = 4)

  expect_type(out, "list")
  expect_length(out, 4)
  for (rep in out) {
    expect_identical(dim(rep), c(6L, 2L))
    expect_identical(colnames(rep), c("row", "col"))
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

test_that("a near-constant column is never drawn, as pca_imp() ignores it", {
  set.seed(3)
  m <- sim_mat(20, 6, perc_total_na = 0)$input
  m[, 3] <- 0.5 + 1e-9 * (seq_len(nrow(m)) %% 2 == 0)
  # the column sits in the gap the threshold closes: variance above zero,
  # below the double.eps that pca_imp() uses for eligibility
  v3 <- col_vars(m)[[3]]
  expect_gt(v3, 0)
  expect_lt(v3, .Machine$double.eps)

  out <- sample_na_loc(
    m,
    n_cols = 1,
    n_rows = 2,
    na_col_subset = 2:3,
    n_reps = 50
  )
  used <- vapply(out, function(r) unique(unname(r[, "col"])), integer(1))
  expect_true(all(used == 2L))

  expect_error(
    sample_na_loc(m, n_cols = 1, n_rows = 2, na_col_subset = 3),
    class = "slideimp_infeasible"
  )

  # Regression: the helper skipped only var <= 0, so a repetition could blank
  # column 3 alone. pca_imp() then refused it ("All columns with missing
  # values are ineligible") and the repetition came back as an error row.
  res <- suppressMessages(tune_imp(
    m,
    parameters = data.frame(ncp = 1),
    .f = "pca_imp",
    n_cols = 1,
    n_rows = 2,
    na_col_subset = 2:3,
    n_reps = 20,
    .progress = FALSE
  ))
  expect_identical(nrow(res), 20L)
  expect_true(all(is.na(res$error)))
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
    "holds infinite values, which are not supported",
    fixed = TRUE
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
    "3 attempts found no valid NA locations for repetition 1",
    fixed = TRUE
  )
})

# pre-existing NAs
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
  # budget, such rows holding one value, and such rows holding two. Every
  # fixture is sampled first; each property is then asserted once over all
  # of them, and a failure names the seeds that break it.
  rowmax <- 0.5
  colmax <- 0.75
  row_cap <- max_miss_count(6L, rowmax)
  col_cap <- max_miss_count(8L, colmax)
  seeds <- 1:150
  fixtures <- lapply(seeds, function(seed) {
    set.seed(seed)
    m <- matrix(sample(c(0.1, 0.2, 0.3, 0.4), 8 * 6, replace = TRUE), nrow = 8)
    m[sample.int(length(m), 12)] <- NA
    out <- tryCatch(
      sample_na_loc(
        m,
        n_cols = 3,
        n_rows = 2,
        rowmax = rowmax,
        colmax = colmax
      ),
      error = function(e) e
    )
    list(m = m, out = out)
  })
  no_bad_seeds <- function(s, ok, what) {
    bad <- s[!ok]
    expect(
      length(bad) == 0L,
      sprintf(
        "%s for %d seed(s): %s",
        what,
        length(bad),
        paste(utils::head(bad, 10L), collapse = ", ")
      )
    )
  }

  # a fixture is refused exactly when one of its own rows is already over
  # rowmax, and by that check: none is refused as infeasible
  refused <- vapply(fixtures, function(f) inherits(f$out, "error"), logical(1))
  over <- vapply(
    fixtures,
    function(f) any(rowSums(is.na(f$m)) > row_cap),
    logical(1)
  )
  no_bad_seeds(seeds, refused == over, "refused is not rows over rowmax")
  msgs <- vapply(
    fixtures[refused],
    function(f) conditionMessage(f$out),
    character(1)
  )
  no_bad_seeds(
    seeds[refused],
    grepl("before any NA is added", msgs, fixed = TRUE),
    "the refusal is not the rowmax pre-check"
  )

  done <- fixtures[!refused]
  done_seeds <- seeds[!refused]
  props <- vapply(
    done,
    function(f) {
      m <- f$m
      loc <- f$out[[1]]
      m2 <- apply_na(m, loc)
      cols <- unique(loc[, "col"])
      c(
        observed = !anyNA(m[loc]),
        unique = anyDuplicated(loc) == 0L,
        n_rows = identical(as.integer(table(loc[, "col"])), c(2L, 2L, 2L)),
        rowmax = max(rowSums(is.na(m2))) <= row_cap,
        colmax = max(colSums(is.na(m2))[cols]) <= col_cap,
        distinct = all(vapply(
          cols,
          function(j) length(unique(m2[!is.na(m2[, j]), j])) >= 2L,
          logical(1)
        ))
      )
    },
    logical(6)
  )
  for (p in rownames(props)) {
    no_bad_seeds(done_seeds, props[p, ], paste(p, "fails"))
  }

  # The branch each drawn column takes, replayed from the draw: the sampler
  # writes columns in the order it fills them, so the rows at their budget
  # when it reaches a column are known. Their observed values in that column
  # number 0, 1, or 2 (two or more).
  branch <- unlist(lapply(done, function(f) {
    m <- f$m
    loc <- f$out[[1]]
    room <- row_cap - rowSums(is.na(m))
    out <- integer(0)
    for (j in unique(loc[, "col"])) {
      held <- m[room <= 0L & !is.na(m[, j]), j]
      out <- c(out, min(length(unique(held)), 2L))
      rows <- loc[loc[, "col"] == j, "row"]
      room[rows] <- room[rows] - 1L
    }
    out
  }))

  # not vacuous: most fixtures are sampled, some are refused, and every
  # branch is reached
  expect_gt(length(done), 75L)
  expect_gt(sum(refused), 0L)
  expect_setequal(branch, 0:2)
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
  skip_if_not(capabilities("profmem"), "R was built without memory profiling")
  n <- 200L
  p <- 2000L
  set.seed(24)
  m <- matrix(rnorm(n * p), n)
  m[sample.int(length(m), 500)] <- NA

  # Regression: `!is.na(obj)` allocated two n x p logicals (1.6 MB each
  # here) to count missing values per row and column.
  f <- tempfile(fileext = ".out")
  withr::defer(unlink(f))
  Rprofmem(f, threshold = 4 * n * p)
  # switched off even if the call errors, not left on for the rest of the run
  withr::defer(Rprofmem(NULL))
  out <- sample_na_loc(m, n_cols = 5, n_rows = 2)
  Rprofmem(NULL)
  big <- grep("^[0-9]+ :", readLines(f), value = TRUE)
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

test_that("column names in na_col_subset need column names on obj", {
  m <- unname(sim_mat(20, 5, perc_total_na = 0)$input)
  err <- expect_error(
    sample_na_loc(m, n_cols = 1, n_rows = 1, na_col_subset = "feature1"),
    "`na_col_subset` names columns, but `obj` has no column names.",
    fixed = TRUE
  )
  expect_identical(class(err), c("rlang_error", "error", "condition"))
})

test_that("na_col_subset names not in obj are listed", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  err <- expect_error(
    sample_na_loc(
      m,
      n_cols = 1,
      n_rows = 1,
      na_col_subset = c("feature1", "nope")
    ),
    "`na_col_subset` names 1 column not in `obj`: nope.",
    fixed = TRUE
  )
  expect_identical(class(err), c("rlang_error", "error", "condition"))
  err <- expect_error(
    sample_na_loc(
      m,
      n_cols = 1,
      n_rows = 1,
      na_col_subset = c("nope", "feature1", "zip")
    ),
    "`na_col_subset` names 2 columns not in `obj`: nope, zip.",
    fixed = TRUE
  )
  expect_identical(class(err), c("rlang_error", "error", "condition"))
})

test_that("na_col_subset of another type is refused", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  for (bad in list(TRUE, list(1))) {
    err <- expect_error(
      sample_na_loc(m, n_cols = 1, n_rows = 1, na_col_subset = bad),
      "`na_col_subset` must be column names, column numbers, or `NULL`.",
      fixed = TRUE
    )
    expect_identical(class(err), c("rlang_error", "error", "condition"))
  }
})

test_that("rows already over rowmax are refused before sampling, by number", {
  m <- sim_mat(20, 5, perc_total_na = 0)$input
  # rowmax = 0.9 of 5 columns allows 4 NAs in a row; row 3 holds 5
  m[3, ] <- NA
  err <- expect_error(
    sample_na_loc(m, n_cols = 1, n_rows = 1),
    "1 row already exceeds `rowmax` before any NA is added.",
    fixed = TRUE
  )
  expect_identical(class(err), c("rlang_error", "error", "condition"))
  expect_match(conditionMessage(err), "Row numbers: 3.", fixed = TRUE)
  m[c(7, 9), ] <- NA
  err <- expect_error(
    sample_na_loc(m, n_cols = 1, n_rows = 1),
    "3 rows already exceed `rowmax` before any NA is added.",
    fixed = TRUE
  )
  expect_match(conditionMessage(err), "Row numbers: 3, 7, 9.", fixed = TRUE)
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

test_that("sample_na_loc returns integer index matrices", {
  # the documented return type, built by sample_each_rep_cpp() itself
  set.seed(1234)
  obj <- sim_mat(20, 30, perc_total_na = 0.2)$input
  locs <- sample_na_loc(obj, n_cols = 5, n_rows = 2, n_reps = 3)

  for (pos in locs) {
    expect_type(pos, "integer")
    expect_equal(colnames(pos), c("row", "col"))
    expect_length(obj[pos], nrow(pos))
  }
})

# sample_each_rep_cpp() ----
test_that("the sampler checks the indices it relies on", {
  # sample_each_rep_cpp() indexes `obj` and its budgets with these values,
  # so one out of range would read past the matrix rather than fail
  base <- list(
    obj = matrix(c(1, 2, 3, 4, 5, 6), 3, 2),
    pool_idx_in = 0:1,
    na_per_col = 1L,
    row_room = rep(1L, 3),
    pool_cap = c(1L, 1L),
    max_attempts = 5L
  )
  # well formed as written, so the failures below are the inputs
  set.seed(1)
  expect_identical(dim(do.call(sample_each_rep_cpp, base)), c(1L, 2L))

  bad <- base
  bad$pool_idx_in <- c(0L, 2L)
  expect_error(
    do.call(sample_each_rep_cpp, bad),
    "pool_idx_in' holds 2",
    fixed = TRUE
  )

  bad <- base
  bad$row_room <- c(1L, NA, 1L)
  expect_error(
    do.call(sample_each_rep_cpp, bad),
    "row_room' holds NA",
    fixed = TRUE
  )

  bad <- base
  bad$row_room <- 1L
  expect_error(
    do.call(sample_each_rep_cpp, bad),
    "row_room' has length 1",
    fixed = TRUE
  )

  bad <- base
  bad$pool_cap <- 1L
  expect_error(
    do.call(sample_each_rep_cpp, bad),
    "differ in length",
    fixed = TRUE
  )

  bad <- base
  bad$na_per_col <- -1L
  expect_error(
    do.call(sample_each_rep_cpp, bad),
    "na_per_col' holds -1",
    fixed = TRUE
  )

  bad <- base
  bad$max_attempts <- NA_integer_
  expect_error(do.call(sample_each_rep_cpp, bad), "max_attempts", fixed = TRUE)

  # a repeated column could fill two slots and return one cell twice
  bad <- base
  bad$pool_idx_in <- c(1L, 1L)
  expect_error(
    do.call(sample_each_rep_cpp, bad),
    "pool_idx_in' holds column index 1 more than once",
    fixed = TRUE
  )
})

test_that("the sampler's draws under a seed are pinned", {
  # Positions drawn under set.seed(1), unchanged since an earlier build of
  # the sampler (76d4c60). A change to the order or number of R::runif()
  # calls moves them. The pool is visited out of order (columns 2, 6, 3), so the
  # shuffle, the protected pair and the commit draws all feed the result.
  obj <- matrix(as.numeric(c(1:30, 30:1)), 10, 6)
  set.seed(1)
  pos <- sample_each_rep_cpp(
    obj,
    pool_idx_in = 0:5,
    na_per_col = c(2L, 2L, 3L),
    row_room = rep(2L, 10),
    pool_cap = rep(4L, 6),
    max_attempts = 10L
  )
  expect_identical(
    pos,
    matrix(c(8L, 3L, 9L, 2L, 4L, 9L, 5L, 2L, 2L, 6L, 6L, 3L, 3L, 3L), 7, 2)
  )
})
