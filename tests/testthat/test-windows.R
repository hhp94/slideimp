skip_on_cran()

test_that("count_window_cover counts every window, duplicate ends included", {
  # windows 1-2 and 3-4 end on the same column
  start <- c(1L, 6L, 11L, 16L, 21L)
  end <- c(10L, 10L, 20L, 20L, 30L)
  expect_gt(anyDuplicated(end), 0L)

  expected <- integer(32L)
  for (i in seq_along(start)) {
    cols <- start[i]:end[i]
    expected[cols] <- expected[cols] + 1L
  }

  # Regression: the counts came from `delta[end + 1] <- delta[end + 1] - 1`,
  # which applies a repeated index once. Two windows ending on column 10
  # lowered the running count by 1 instead of 2, so every later column was
  # counted one too high. slide_imp() then halved singly covered columns and
  # left zeros in columns it should have restored.
  expect_identical(count_window_cover(start, end, 32L), expected)
  # columns past the last window are covered by nothing
  expect_identical(count_window_cover(start, end, 32L)[31:32], c(0L, 0L))
  # no windows at all
  expect_identical(count_window_cover(integer(0), integer(0), 4L), integer(4L))
})

test_that("find_windows emits no window nested inside the one before it", {
  # five clusters of 10 columns, 0.1 apart, each followed by a gap of 0.6
  loc <- unlist(lapply(0:4, function(g) g * 1.5 + seq(0, 0.9, by = 0.1)))
  w <- find_windows(loc, window = 1, overlap = 0.5)

  # Regression: the start found in the trailing overlap region was used even
  # when the next coordinate was a full window away, which produced [1, 10]
  # followed by [6, 10]: a window that adds no column and shares its end.
  expect_identical(w$start, c(1L, 11L, 21L, 31L, 41L))
  expect_identical(w$end, c(10L, 20L, 30L, 40L, 50L))

  # the property on irregular coordinates, at several overlaps
  set.seed(1)
  gappy <- sort(c(runif(200, 0, 10), runif(200, 30, 40), runif(50, 100, 101)))
  for (ov in c(0, 1, 2.5, 4.9)) {
    w <- find_windows(gappy, window = 5, overlap = ov)
    expect_true(all(diff(w$start) > 0))
    expect_true(all(diff(w$end) > 0))
    # every column is still inside some window
    covered <- sort(unique(sequence(w$end - w$start + 1L, from = w$start)))
    expect_identical(covered, seq_along(gappy))
    # and a positive overlap still produces overlapping windows, so the fix
    # did not turn every window disjoint
    if (ov > 0) {
      expect_true(any(w$start[-1L] <= w$end[-length(w$end)]))
    }
  }
})
