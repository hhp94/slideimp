skip_on_cran()

test_that("print.si_flank_plan() counts targets by status, and n truncates the table", {
  set.seed(7)
  y <- matrix(stats::rnorm(20 * 30), 20)
  dimnames(y) <- list(sprintf("s%02d", 1:20), sprintf("c%02d", 1:30))
  y[c(2, 6), 8] <- NA
  y[c(5, 11), 20] <- NA
  y[c(1, 3), 25] <- NA
  plan <- mipca_flank_plan(
    y,
    seq_len(30),
    window_size = 5,
    subset = c(3, 8, 20, 25),
    min_window_n = 5
  )
  tg <- plan$targets
  # c03 has no missing value, so both status counts are exercised
  n_imp <- sum(tg$status == "impute")
  n_cln <- sum(tg$status == "clean")
  expect_gt(n_imp, 0L)
  expect_gt(n_cln, 0L)
  expect_identical(n_imp + n_cln, nrow(tg))

  out <- capture.output(v <- withVisible(print(plan, n = 2L)))
  expect_false(v$visible)
  expect_identical(v$value, plan)
  expect_identical(
    out[1:2],
    c(
      sprintf(
        "si_flank_plan: %d targets on a %d x %d matrix, nothing fitted",
        nrow(tg),
        nrow(y),
        ncol(y)
      ),
      sprintf("  status: impute %d, clean %d", n_imp, n_cln)
    )
  )
  shown <- vapply(tg$name, function(s) any(grepl(s, out, fixed = TRUE)), NA)
  expect_identical(unname(shown), seq_len(nrow(tg)) <= 2L)
  expect_true(
    sprintf("  ... %d more rows in x$targets", nrow(tg) - 2L) %in% out
  )

  out_all <- capture.output(print(plan))
  expect_true(all(vapply(tg$name, function(s) any(grepl(s, out_all)), NA)))
  expect_false(any(grepl("more row", out_all, fixed = TRUE)))
})

test_that("mipca_flank_plan() reads only the columns a window can hold", {
  set.seed(11)
  X <- matrix(stats::rnorm(30 * 16), 30, 16)
  dimnames(X) <- list(sprintf("r%02d", 1:30), sprintf("c%02d", 1:16))
  X[1:2, c(2, 8, 12)] <- NA_real_
  plan_of <- function(Y) {
    mipca_flank_plan(
      Y,
      location = 1:16,
      window_size = 2,
      subset = c(2, 8, 12),
      min_window_n = 3,
      abb = c(3, 8, 15)
    )
  }
  # the spans of c02 (columns 1 to 4) and c12 (10 to 14), less the 'abb'
  # column c03 inside them, and the 'abb' target c08, whose span no window
  # needs
  read <- c(1L, 2L, 4L, 8L, 10:14)
  real_scan <- col_scan
  seen <- list()
  local_mocked_bindings(
    col_scan = function(obj, cols = NULL, ...) {
      seen[[length(seen) + 1L]] <<- cols
      real_scan(obj, cols = cols, ...)
    }
  )
  plan <- plan_of(X)
  expect_identical(seen, list(read))
  expect_identical(plan$targets$status, c("impute", "abb", "impute"))

  # an Inf in a column the plan does not read leaves the plan as it was, and
  # one in a column it reads is refused
  unread <- setdiff(1:16, read)
  expect_true(all(c(3L, 9L, 16L) %in% unread))
  for (j in unread) {
    Y <- X
    Y[5, j] <- Inf
    expect_identical(
      plan_of(Y)$targets,
      plan$targets,
      info = paste("column", j)
    )
  }
  for (j in read) {
    Y <- X
    Y[5, j] <- Inf
    expect_error(
      plan_of(Y),
      "'obj' holds infinite values, which are not supported",
      fixed = TRUE,
      info = j
    )
  }
})
