skip_if_not_manual()

# tune_imp() under mirai daemons, run on request ----
# The daemons load the INSTALLED slideimp, not the load_all() build.
test_that("tune_imp runs the parallel (mirai) branch", {
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

test_that("tune_imp says cores > 1 is ignored under active daemons", {
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
    regexp = "`cores` = 2 is ignored"
  )
})

test_that("a task whose daemon stops is recorded as a failed iteration", {
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
    "tuning iteration",
    fixed = TRUE
  )
  expect_match(res$error, "'shift' not found", all = TRUE)
})

test_that("under mirai a wrong return on a later set is recorded", {
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
    "tuning iteration",
    fixed = TRUE
  )
  expect_true(is.na(res$error[1]))
  expect_match(res$error[2], "returned a 40 x 30 matrix")
  expect_identical(nrow(res$result[[2]]), 0L)
})

test_that("under mirai every task masks a clean copy, and the daemons drop it", {
  skip_if_not_installed("withr")
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  # the estimate is the number of NA in the matrix the task was handed, so a
  # daemon copy that kept an earlier rep's mask would show in it
  fill_count <- function(obj, v) {
    obj[is.na(obj)] <- sum(is.na(obj)) + v
    obj
  }
  na_loc <- sample_na_loc(obj, num_na = 10, n_reps = 3)
  pars <- data.frame(v = c(0, 1))
  want <- suppressMessages(
    tune_imp(obj, pars, .f = fill_count, na_loc = na_loc, .progress = FALSE)
  )
  # what the sequential path hands each task: the input's own NA plus its 10
  n_in <- sum(is.na(obj))
  est <- unlist(lapply(want$result, `[[`, "estimate"))
  expect_identical(sort(unique(est)), n_in + 10 + c(0, 1))

  # one daemon, so all six tasks run on its one copy, in turn
  mirai::daemons(1)
  withr::defer(mirai::daemons(0))
  got <- suppressMessages(
    tune_imp(obj, pars, .f = fill_count, na_loc = na_loc, .progress = FALSE)
  )
  expect_true(all(is.na(got$error)))
  expect_identical(got$result, want$result)
  # the copy is gone once tune_imp() returns
  left <- unlist(mirai::everywhere(
    exists(".slideimp_tune_obj", envir = globalenv(), inherits = FALSE)
  )[])
  expect_identical(left, FALSE)
})

test_that("under mirai a daemon without the input records why", {
  skip_if_not_installed("withr")
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  # removing the daemon's copy stands in for a daemon that connected after
  # tune_imp() sent the input; the probe in the main session has none
  drop_copy <- function(obj) {
    suppressWarnings(rm(list = ".slideimp_tune_obj", envir = globalenv()))
    obj[is.na(obj)] <- 0
    obj
  }
  mirai::daemons(1)
  withr::defer(mirai::daemons(0))
  expect_warning(
    res <- suppressMessages(
      tune_imp(
        obj,
        NULL,
        .f = drop_copy,
        num_na = 10,
        n_reps = 2,
        .progress = FALSE
      )
    ),
    "tuning iteration",
    fixed = TRUE
  )
  expect_true(is.na(res$error[1]))
  expect_match(res$error[2], "did not reach this daemon", fixed = TRUE)
})

test_that("under mirai a hand-off that fails stops before any task", {
  skip_if_not_installed("withr")
  set.seed(1)
  obj <- sim_mat(20, 30)$input
  mirai::daemons(1)
  withr::defer(mirai::daemons(0))
  # what a daemon that cannot take its copy hands back
  failed <- mirai::mirai(stop("no room for obj"))[]
  expect_true(mirai::is_error_value(failed))
  ran <- FALSE
  local_mocked_bindings(
    everywhere = function(...) list(failed),
    mirai_map = function(...) {
      ran <<- TRUE
      stop("a task was dispatched")
    },
    .package = "mirai"
  )
  expect_error(
    suppressMessages(tune_imp(
      obj,
      data.frame(k = 3),
      .f = "knn_imp",
      num_na = 10,
      .progress = FALSE
    )),
    "'obj' did not reach 1 of the 1 mirai daemons. The first one reported: ",
    fixed = TRUE
  )
  expect_false(ran)
})

test_that(".progress = FALSE draws no progress bar on the mirai path", {
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
    status <- "Tuning|Step [12]/2|Running mode|^\\s*$"
    length(txt[!grepl(status, txt)])
  }
  expect_identical(drawn(FALSE), 0L)
  # the detector can see a bar: without this the line above cannot fail
  expect_gt(drawn(TRUE), 0L)
  expect_identical(drawn(FALSE), 0L)
})
