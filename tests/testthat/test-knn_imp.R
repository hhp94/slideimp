test_that("`impute_knn_brute` calculate the missing location correctly", {
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.5, perc_col_na = 1)$input
  miss <- is.na(to_test)
  cmiss <- colSums(miss)
  miss_rate <- cmiss / nrow(to_test)
  # same preprocessing knn_imp() does before calling the low-level functions
  colmax <- 0.9
  eligible <- miss_rate <= colmax
  pre_imp_cols <- to_test[, eligible, drop = FALSE]
  pre_imp_miss <- miss[, eligible, drop = FALSE]

  local_has_miss <- which(cmiss[eligible] > 0L)
  grp_impute <- as.integer(local_has_miss - 1L)

  expected_local <- unname(which(pre_imp_miss, arr.ind = TRUE))

  imputed_brute <- impute_knn_brute(
    obj = pre_imp_cols,
    k = 5,
    grp_impute = grp_impute,
    grp_miss_no_imp = integer(0L),
    grp_complete = integer(0L),
    method = 0L,
    dist_pow = 1,
    cores = 1
  )

  # extract only the location columns (row, local_col_1based) that the C++ already returns
  idx_brute <- imputed_brute[, 1:2, drop = FALSE]

  # sort so the comparison is order-independent
  expected_sorted <- expected_local[
    order(expected_local[, 1], expected_local[, 2]),
    ,
    drop = FALSE
  ]
  brute_sorted <- idx_brute[
    order(idx_brute[, 1], idx_brute[, 2]),
    ,
    drop = FALSE
  ]

  expect_equal(brute_sorted, expected_sorted)
})

test_that("knn_imp fills every missing cell and touches nothing else", {
  set.seed(1234)
  obj <- sim_mat(50, 100)$input
  na_cells <- is.na(obj)
  expect_gt(sum(na_cells), 0L)

  r <- knn_imp(obj, k = 3, method = "euclidean")
  expect_false(anyNA(r[,]))
  # the observed cells come back bit for bit: they are copied, not recomputed
  expect_identical(r[,][!na_cells], obj[!na_cells])
})

# The comparison against impute::impute.knn lives in dev/probe_impute_knn.R,
# not here. `impute` is a Bioconductor package and is not in Suggests, so a
# `impute::` call anywhere under tests/ makes R CMD check report an unstated
# dependency even when it is guarded by skips. The script asserts the same
# thing the test would: knn_imp() agrees with impute.knn to 2.2e-16 absolute
# and 4.5e-16 relative on a 100 x 100 matrix, at k = 3 and k = 10. Run it with
# `Rscript dev/probe_impute_knn.R` after touching the neighbor search.

test_that("`subset` feature of `knn_imp` works with post_imp = FALSE/TRUE", {
  set.seed(1234)
  to_test <- sim_mat(20, 50, perc_total_na = 0.2, perc_col_na = 1)$input
  # Impute just 3 columns
  ## Check subset using numeric index
  r1 <- knn_imp(to_test, k = 3, post_imp = FALSE, subset = c(1, 3, 5))
  expect_true(!anyNA(r1[, c(1, 3, 5)]))
  expect_equal(is.na(r1[, -c(1, 3, 5)]), is.na(to_test[, -c(1, 3, 5)]))
  ## Check subset using character vector
  r2 <- knn_imp(
    to_test,
    k = 3,
    post_imp = FALSE,
    subset = paste0("feature", c(1, 3, 5))
  )
  expect_equal(r1, r2)

  # Test with post_imp = TRUE and a column requiring post imputation
  to_test_post <- to_test
  # Column 5 will be colMeans if post_imp is TRUE
  to_test_post[2:nrow(to_test_post), 5] <- NA
  r3 <- knn_imp(to_test_post, k = 3, post_imp = TRUE, subset = c(1, 3, 5))
  expect_true(!anyNA(r3[, c(1, 3, 5)]))
  # Expect that only the subset columns are imputed. The rests are untouched
  expect_equal(is.na(r3[, -c(1, 3, 5)]), is.na(to_test_post[, -c(1, 3, 5)]))
  # Verify post_imp on column 5
  col5_mean <- mean(to_test_post[, 5], na.rm = TRUE)
  expect_equal(unname(r3[, 5]), rep(col5_mean, nrow(to_test_post)))
  r4 <- knn_imp(
    to_test_post,
    k = 3,
    post_imp = TRUE,
    subset = paste0("feature", c(1, 3, 5))
  )
  expect_equal(r3, r4)
})

test_that("Behavior with extreme missing columns and rows", {
  set.seed(1234)
  to_test <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  # row 1 is all NA
  to_test[1, ] <- NA
  expect_no_error(knn_imp(to_test, k = 3, post_imp = FALSE))
  expect_true(!anyNA(knn_imp(to_test, k = 3, post_imp = TRUE)))

  to_test[, 1] <- NA
  expect_error(knn_imp(to_test, k = 3, post_imp = FALSE), "All NA/NaN")

  # not admissible
  mat <- matrix(NA, nrow = 20, ncol = 20)
  diag(mat) <- rnorm(20)
  expect_error(knn_imp(mat, k = 2), "exceeds usable columns")
})

test_that("colmax is inclusive: a column exactly at colmax is eligible", {
  set.seed(1234)
  obj <- matrix(rnorm(50), 10, 5)
  obj[1:5, 1] <- NA # miss_rate exactly 0.5

  # at the threshold: column 1 is imputed, same as pca_imp() treats it
  r <- knn_imp(obj, k = 2, colmax = 0.5, post_imp = FALSE)
  expect_false(anyNA(r[, 1]))
  expect_s3_class(r, "slideimp_results")

  # just below it: column 1 is the only column with missing values, so the
  # subset is empty and the second abort site fires
  expect_error(
    knn_imp(obj, k = 2, colmax = 0.49, post_imp = FALSE),
    class = "slideimp_infeasible"
  )
})

test_that("neighbors sharing no observed row with the target are dropped", {
  set.seed(1234)
  n <- 8
  target <- rnorm(n)
  target[5:8] <- NA
  dup <- target
  dup[5:8] <- rnorm(4) # identical to target on the observed rows
  disjoint <- rnorm(n)
  disjoint[1:4] <- NA # observed only where the target is missing
  obj <- cbind(target, dup, disjoint)

  # the disjoint column has infinite distance to the target and must not be
  # averaged in, whether or not distances are used as weights
  for (dp in c(0, 1)) {
    r <- knn_imp(
      obj,
      k = 2,
      dist_pow = dp,
      subset = 1,
      post_imp = FALSE,
      na_check = FALSE
    )
    got <- unname(r[5:8, 1])
    want <- dup[5:8]
    abs_err <- max(abs(got - want))
    rel_err <- max(abs(got - want) / abs(want))
    expect_lte(abs_err, 1e-12)
    expect_lte(rel_err, 1e-12)
    expect_gt(max(abs(got - unname(obj[5:8, 3]))), 0)
  }

  # no candidate overlaps the target at all: nothing to average, stays NA
  obj2 <- cbind(target, disjoint, disjoint2 = disjoint + 1)
  r2 <- knn_imp(
    obj2,
    k = 2,
    subset = 1,
    post_imp = FALSE,
    na_check = FALSE
  )
  expect_true(all(is.na(r2[5:8, 1])))
  expect_equal(unname(r2[1:4, 1]), target[1:4])
})

# Pure-R brute-force K-NN for a single target column, written from the Details
# section of ?knn_imp rather than from the kernel: the distance to a candidate
# is a mean over the rows where both are observed, a candidate sharing no
# observed row is not a neighbor, the k smallest distances win, and a neighbor
# at distance d is weighted (d_min / d)^dist_pow with d_min the smallest
# positive distance among them.
#
# `dist_fun` maps the vector of co-observed differences to the distance, so the
# same reference serves both metrics and the square-rooted euclidean distance.
knn_ref <- function(obj, target, miss_rows, k, dist_fun, dist_pow) {
  cand <- setdiff(seq_len(ncol(obj)), target)
  x <- obj[, target]
  obs <- !is.na(x)
  d <- vapply(
    cand,
    function(j) {
      both <- obs & !is.na(obj[, j])
      if (!any(both)) Inf else dist_fun(x[both] - obj[both, j])
    },
    numeric(1)
  )
  ord <- order(d)[seq_len(k)]
  ord <- ord[is.finite(d[ord])]
  dk <- d[ord]
  w <- if (dist_pow == 0) {
    rep(1, length(ord))
  } else {
    (min(dk[dk > 0]) / dk)^dist_pow
  }
  # The neighbors are chosen once for the column, but a neighbor that is
  # itself missing at the row being imputed drops out of that row's average,
  # numerator and denominator both.
  vapply(
    miss_rows,
    function(r) {
      v <- obj[r, cand[ord]]
      ok <- !is.na(v)
      if (!any(ok)) NA_real_ else sum(w[ok] * v[ok]) / sum(w[ok])
    },
    numeric(1)
  )
}

# Report max absolute and max relative error; assert on the worst case.
expect_close <- function(got, want, tol = 1e-12) {
  expect_gt(length(got), 0L)
  abs_err <- max(abs(got - want))
  big <- abs(want) > 1e-8
  expect_true(any(big))
  rel_err <- max(abs(got[big] - want[big]) / abs(want[big]))
  expect_lt(abs_err, tol)
  expect_lt(rel_err, tol)
}

test_that("the documented weighting rule reproduces the kernel", {
  set.seed(20260905)
  n_row <- 30
  k <- 3
  miss_rows <- c(3L, 11L, 24L)

  # Two candidate layouts, because they take different code paths. With every
  # candidate complete the kernel takes its group-3 route, which reads obj
  # unmasked and knows the co-observed count up front. With candidates that are
  # themselves missing, those go to group 2 and the masked route runs instead.
  complete_cand <- matrix(rnorm(n_row * 7L), nrow = n_row)
  complete_cand[miss_rows, 1] <- NA

  mixed_cand <- matrix(rnorm(n_row * 9L), nrow = n_row)
  mixed_cand[miss_rows, 1] <- NA
  mixed_cand[c(2L, 17L), 2] <- NA
  mixed_cand[c(5L, 11L, 28L), 3] <- NA
  mixed_cand[c(3L, 9L), 4] <- NA

  scenarios <- list(
    "complete candidates" = complete_cand,
    "mixed candidates" = mixed_cand
  )

  d_euclidean <- function(x) mean(x^2)
  d_manhattan <- function(x) mean(abs(x))
  d_rooted <- function(x) sqrt(mean(x^2))

  for (nm in names(scenarios)) {
    obj <- scenarios[[nm]]
    cand_miss <- colSums(is.na(obj))[-1]
    # prove the scenario is the one it claims to be, so neither branch of the
    # loop can pass by never reaching the path it is here to cover
    if (nm == "complete candidates") {
      expect_identical(sum(cand_miss > 0L), 0L)
    } else {
      expect_gt(sum(cand_miss > 0L), 0L)
      expect_gt(sum(cand_miss == 0L), 0L)
    }
    # more candidates than k, so the pruned replacement phase runs too
    expect_gt(ncol(obj) - 1L, k)

    for (dp in c(0, 1, 2)) {
      for (m in c("euclidean", "manhattan")) {
        got <- knn_imp(
          obj,
          k = k,
          method = m,
          dist_pow = dp,
          subset = 1,
          post_imp = FALSE,
          na_check = FALSE
        )
        want <- knn_ref(
          obj,
          1L,
          miss_rows,
          k,
          if (m == "euclidean") d_euclidean else d_manhattan,
          dp
        )
        expect_close(unname(got[miss_rows, 1]), want)
      }
    }

    # The euclidean distance is never square-rooted, so dist_pow acts as an
    # exponent of 2 * dist_pow on the square-rooted euclidean distance.
    for (dp in c(1, 2)) {
      got <- unname(
        knn_imp(
          obj,
          k = k,
          method = "euclidean",
          dist_pow = dp,
          subset = 1,
          post_imp = FALSE,
          na_check = FALSE
        )[miss_rows, 1]
      )
      expect_close(got, knn_ref(obj, 1L, miss_rows, k, d_rooted, 2 * dp))
      # and NOT the same as weighting the rooted distance by dist_pow itself,
      # which is what the help page used to claim.
      expect_gt(
        max(abs(got - knn_ref(obj, 1L, miss_rows, k, d_rooted, dp))),
        1e-6
      )
    }
  }
})

test_that("tied neighbors are settled by scan order, not by the library sort", {
  # Regression test. The fill phase of distance_vector_impl() used to order
  # its first k candidates with std::sort, which leaves equal distances in an
  # order of the library's choosing, and a closer candidate then evicts the
  # LAST entry - so the library decided which of several tied neighbors was
  # dropped. libstdc++ insertion-sorts up to 16 elements, which keeps ties in
  # order, so the bug showed only from k = 17, where it evicted tied
  # candidates from the front of the run. The rule is first-encountered wins,
  # in the kernel's scan order: complete columns, then columns with missing
  # values, each in column order. knn_ref() above breaks ties by column order
  # instead, and agrees with the kernel only because its rnorm data has none.
  #
  # Column 1 is the target, missing at row 1 only. The first k candidates in
  # scan order are the target plus a +-1 pattern on its observed rows, at
  # distance exactly 1 under both metrics; the n_closer after them use +-1/2
  # and each evicts one tied candidate. Candidate column j holds 2^(j - 2) at
  # row 1, so with dist_pow = 0 the imputed value times k is a sum of distinct
  # powers of two, one per chosen neighbor: the set can be read back from it.
  # Every sum here is of integers below 2^53, so both sides are exact.
  n_obs <- 8L
  n_closer <- 4L
  ks <- c(5L, 17L, 32L)
  # k > 16 is where libstdc++'s std::sort stops keeping ties in order
  expect_gt(max(ks), 16L)

  for (layout in c("complete", "masked", "mixed")) {
    for (k in ks) {
      m <- k + n_closer
      kind <- switch(
        layout,
        complete = rep("c", m),
        masked = rep("m", m),
        # complete and incomplete alternate, so scan order is not column order
        mixed = c(rep_len(c("c", "m"), k), rep("m", n_closer))
      )
      set.seed(k)
      t_obs <- sample(-5:5, n_obs, replace = TRUE)
      sgn <- matrix(sample(c(-1, 1), n_obs * m, replace = TRUE), n_obs, m)
      step <- rep(c(1, 0.5), c(k, n_closer))
      cand <- t_obs + sweep(sgn, 2, step, "*")
      for (j in which(kind == "m")) {
        cand[1L + j %% n_obs, j] <- NA # in a row the target observes
      }
      obj <- unname(rbind(c(NA, 2^(seq_len(m) - 1)), cbind(t_obs, cand)))

      # the case is the one it claims: every fill candidate is tied, every
      # candidate after the fill is strictly closer, and evictions happen
      gap <- abs(cand - t_obs)
      expect_true(all(gap[, step == 1] == 1, na.rm = TRUE))
      expect_true(all(gap[, step < 1] == 0.5, na.rm = TRUE))
      cols <- seq_len(m) + 1L
      scan <- c(cols[kind == "c"], cols[kind == "m"])
      tied <- scan[seq_len(k)]
      expect_identical(sort(tied), cols[step == 1])

      # the last n_closer tied candidates in scan order are the ones evicted
      closer <- cols[step < 1]
      want_nn <- sort(c(tied[seq_len(k - n_closer)], closer))
      want <- sum(obj[1L, want_nn]) / k
      # and any other eviction is detectable - here, evicting the first ones
      wrong_nn <- c(tied[-seq_len(n_closer)], closer)
      expect_gt(abs(sum(obj[1L, wrong_nn]) / k - want), 0.5)

      for (method in c("euclidean", "manhattan")) {
        r <- knn_imp(
          obj,
          k = k,
          method = method,
          dist_pow = 0,
          subset = 1,
          post_imp = FALSE,
          na_check = FALSE
        )
        got <- unname(r[1L, 1L])
        expect_close(got, want)
        s <- round(got * k)
        got_nn <- cols[floor(s / 2^(seq_len(m) - 1)) %% 2 == 1]
        expect_identical(got_nn, want_nn, info = paste(layout, k, method))
      }
    }
  }
})

test_that("the parallel path gives the same answer as the serial one", {
  # par_for() (src/par_for.h) runs the loop body inline at one thread and
  # through RcppThread's pool above that, so the two have to agree. Each
  # column is computed independently in either case, so the arithmetic per
  # column is the same sequence of operations: within one process that makes
  # identical() the right comparator.
  set.seed(31415)
  obj <- sim_mat(60, 40, perc_total_na = 0.15, perc_col_na = 1)$input
  expect_gt(sum(is.na(obj)), 0L)

  for (dp in c(0, 1)) {
    r1 <- knn_imp(obj, k = 5, dist_pow = dp, cores = 1, na_check = FALSE)
    r2 <- knn_imp(obj, k = 5, dist_pow = dp, cores = 4, na_check = FALSE)
    expect_identical(r1, r2)
  }

  expect_identical(col_vars(obj, cores = 1), col_vars(obj, cores = 4))
  expect_identical(
    mean_imp_col(obj, cores = 1),
    mean_imp_col(obj, cores = 4)
  )
})

test_that("the kernel checks the invariants it relies on", {
  # impute_knn_brute() is reachable only from inside the package, but the
  # tests call it directly and knn_imp() is not the only caller it has to
  # survive. These three inputs used to produce silent nonsense.
  set.seed(4242)
  obj <- matrix(rnorm(20 * 5), 20, 5)
  obj[3:6, 1] <- NA

  base <- list(
    obj = obj,
    k = 2,
    grp_impute = 0L,
    grp_miss_no_imp = integer(0),
    grp_complete = 1:4, # 0-based: columns 2 to 5
    method = 0L,
    dist_pow = 0,
    cores = 1L
  )
  # the call is well formed as written, so the failures below are the inputs
  expect_no_error(do.call(impute_knn_brute, base))

  # a column claimed by two groups would be its own neighbor at distance zero
  overlap <- base
  overlap$grp_complete <- 0:3
  expect_error(do.call(impute_knn_brute, overlap), "disjoint")

  # group 3 is read unmasked, so a missing value there is never excluded
  nan_complete <- base
  nan_complete$obj[2, 3] <- NA
  expect_error(do.call(impute_knn_brute, nan_complete), "fully observed")

  # a target with no observed row overlaps nothing, so it has no neighbors
  # and its cells stay missing rather than being filled from a NaN distance
  all_na <- base
  all_na$obj[, 1] <- NA
  r <- do.call(impute_knn_brute, all_na)
  expect_identical(nrow(r), nrow(obj))
  expect_true(all(is.na(r[, 3])))
})

test_that("post_imp fills columns K-NN could not finish", {
  # The two ways a subset column can still be missing after the kernel, which
  # are the two terms knn_imp() uses to decide whether mean_imp_col() has any
  # work: a column held out by colmax, and a column with no overlapping
  # candidate. Both must end up at their column mean.
  set.seed(606)
  n <- 12
  obj <- matrix(rnorm(n * 5), n, 5)
  obj[2:n, 2] <- NA # miss_rate 11/12 = 0.917, above the default colmax
  obj[1:6, 3] <- NA # ordinary missingness, K-NN handles it

  r <- knn_imp(obj, k = 2, post_imp = TRUE)
  expect_false(anyNA(r[,]))
  expect_equal(unname(r[2:n, 2]), rep(obj[1, 2], n - 1L))
  expect_false(anyNA(r[1:6, 3]))

  # a target that overlaps no candidate: K-NN leaves it, post_imp fills it
  target <- rnorm(n)
  target[7:n] <- NA
  disjoint <- rnorm(n)
  disjoint[1:6] <- NA
  obj2 <- cbind(target, disjoint, disjoint + 1, disjoint + 2)
  r2 <- knn_imp(obj2, k = 2, subset = 1, post_imp = TRUE)
  expect_false(anyNA(r2[, 1]))
  expect_equal(unname(r2[7:n, 1]), rep(mean(target[1:6]), n - 6L))
  # the untargeted columns keep their missingness
  expect_equal(is.na(r2[, -1]), is.na(obj2[, -1]))
})

test_that("early returns carry the slideimp_results class", {
  set.seed(99)
  obj <- matrix(rnorm(40), nrow = 10)
  obj[2, 4] <- NA

  # Empty subset: resolve_subset() reports it, and knn_imp must not repeat it.
  msgs <- testthat::capture_messages(
    r1 <- knn_imp(obj, k = 2, subset = integer(0))
  )
  expect_length(msgs, 1L)
  expect_match(msgs, "No features")
  expect_s3_class(r1, "slideimp_results")
  expect_identical(attr(r1, "imp_method"), "knn")
  expect_identical(attr(r1, "fallback"), FALSE)
  expect_equal(r1[,, drop = FALSE], obj)

  # Subset with nothing missing in it.
  msgs2 <- testthat::capture_messages(
    r2 <- knn_imp(obj, k = 2, subset = 1:3)
  )
  expect_length(msgs2, 1L)
  expect_match(msgs2, "No missing values in subset columns")
  expect_s3_class(r2, "slideimp_results")
  expect_equal(r2[,, drop = FALSE], obj)

  # na_check decides whether the NA accounting is done, on every return path.
  expect_false(is.null(attr(r2, "has_remaining_na")))
  r3 <- suppressMessages(knn_imp(obj, k = 2, subset = 1:3, na_check = FALSE))
  expect_null(attr(r3, "has_remaining_na"))

  # The class does not depend on which path the data takes.
  r_main <- knn_imp(obj, k = 2, subset = 4)
  expect_identical(class(r1), class(r_main))
  expect_identical(class(r2), class(r_main))
})

test_that("dist_pow is validated by name before it reaches C++", {
  set.seed(11)
  obj <- matrix(rnorm(40), 10, 4)
  obj[2, 1] <- NA

  # "1" used to pass, because "1" >= 0 is a string comparison in R, and the
  # failure surfaced at the .Call as a type error naming no argument at all.
  for (bad in list("1", TRUE, -1, Inf, c(1, 2), NA_real_, NULL)) {
    expect_error(knn_imp(obj, k = 2, dist_pow = bad), "dist_pow")
  }

  expect_no_error(knn_imp(obj, k = 2, dist_pow = 0L))
  expect_no_error(knn_imp(obj, k = 2, dist_pow = 2.5))
})

test_that("Throw on Inf", {
  set.seed(1234)
  to_test <- sim_mat(20, 20, perc_total_na = 0.2, perc_col_na = 1)$input
  to_test[1, 1] <- Inf
  expect_error(knn_imp(to_test, k = 3, post_imp = FALSE), "Infinite")
  to_test[1, 1] <- -Inf
  expect_error(knn_imp(to_test, k = 3, post_imp = FALSE), "Infinite")
})
