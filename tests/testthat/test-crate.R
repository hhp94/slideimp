# -----------------------------------------------------------------------------
# Reduced base-R re-implementation of carrier::crate()
#
# Derived from carrier: https://github.com/r-lib/carrier
# Copyright (c) Posit Software, PBC
# Original implementation by the carrier authors
# Licensed under the MIT License (see inst/LICENSE.note)
# -----------------------------------------------------------------------------

test_that(".crate() requires named `...` arguments", {
  crate <- getFromNamespace(".crate", "slideimp")

  expect_error(
    crate(function(x) identity(x), x = 1, y = 2, 3),
    "named"
  )

  expect_no_error(crate(function(x) identity(x), x = 1, y = 2))
  expect_no_error(crate(function(x) identity(x)))
})

test_that(".crate() requires a function", {
  crate <- getFromNamespace(".crate", "slideimp")

  expect_error(crate(1), "function")
})

test_that(".new_crate() creates crate objects", {
  new_crate <- getFromNamespace(".new_crate", "slideimp")

  fn <- new_crate(function() NULL)
  expect_s3_class(fn, "crate")
})

test_that(".crate() can supply data", {
  crate <- getFromNamespace(".crate", "slideimp")

  fn <- crate(function() toupper(foo), foo = "foo")

  expect_identical(fn(), "FOO")
})

test_that(".crate() isolates undeclared globals", {
  crate <- getFromNamespace(".crate", "slideimp")

  foobar <- "foobar"
  fn <- crate(function() toupper(foobar), foo = "foo")

  expect_error(fn(), "not found")
})

test_that(".crate() supports serialization roundtrip", {
  crate <- getFromNamespace(".crate", "slideimp")

  fn <- crate(function() toupper(foo), foo = "foo")
  out <- unserialize(serialize(fn, NULL))

  expect_equal(as.list(environment(fn)), as.list(environment(out)))
  expect_identical(fn(), out())
})

test_that("helper functions can be passed through .crate(...)", {
  crate <- getFromNamespace(".crate", "slideimp")

  really_do_it <- function() "foo"
  do_it <- function(x) really_do_it()

  environment(really_do_it) <- globalenv()
  environment(do_it) <- globalenv()

  fn <- crate(
    function(x) do_it(x),
    do_it = do_it,
    really_do_it = really_do_it
  )

  expect_identical(fn(NULL), "foo")
})

test_that("closures passed through .crate(...) lose undeclared parent locals", {
  crate <- getFromNamespace(".crate", "slideimp")

  outer <- function() {
    secret <- "bar"
    helper <- function() secret

    fn <- crate(function() helper(), helper = helper)
    fn()
  }

  expect_error(outer(), "not found")
})

test_that("closures passed through .crate(...) work when dependencies are declared", {
  crate <- getFromNamespace(".crate", "slideimp")

  outer <- function() {
    secret <- "bar"
    helper <- function() secret

    fn <- crate(
      function() helper(),
      helper = helper,
      secret = secret
    )

    fn()
  }

  expect_identical(outer(), "bar")
})

test_that("functions inside containers are not recursively re-homed", {
  crate <- getFromNamespace(".crate", "slideimp")

  outer <- function() {
    secret <- "bar"
    helper_list <- list(fn = function() secret)

    fn <- crate(
      function() helper_list$fn(),
      helper_list = helper_list
    )

    fn()
  }

  expect_identical(outer(), "bar")
})

test_that("namespace functions are not re-homed", {
  crate <- getFromNamespace(".crate", "slideimp")

  fn <- crate(
    function(x) pkg_fun(x),
    pkg_fun = stats::median
  )

  expect_equal(fn(c(3, 1, 2)), 2)

  pkg_fun <- get("pkg_fun", envir = environment(fn), inherits = FALSE)
  expect_true(isNamespace(environment(pkg_fun)))
})

test_that("unqualified non-base functions are not visible unless declared", {
  crate <- getFromNamespace(".crate", "slideimp")

  fn <- crate(function(x) median(x))

  expect_error(fn(1:3), "could not find function|not found")
})

test_that("existing crates passed through `...` are not re-crated", {
  crate <- getFromNamespace(".crate", "slideimp")

  fn2 <- crate(function() "ok")
  fn2_env <- environment(fn2)

  fn <- crate(function() fn2(), fn2 = fn2)

  expect_identical(fn(), "ok")
  expect_identical(environment(get("fn2", environment(fn))), fn2_env)
})

test_that(".crate() rejects duplicate `...` names", {
  crate <- getFromNamespace(".crate", "slideimp")
  expect_error(crate(function() NULL, x = 1, x = 2), "named")
})

test_that(".crate() honours a custom .parent_env", {
  crate <- getFromNamespace(".crate", "slideimp")
  e <- new.env()
  e$visible <- "yes"
  fn <- crate(function() visible, .parent_env = e)
  expect_identical(fn(), "yes")
})

test_that(".crate() works through mirai", {
  skip_if_not_manual()
  skip_if_not_installed("withr")
  crate <- getFromNamespace(".crate", "slideimp")
  mirai::daemons(1L)
  withr::defer(mirai::daemons(0L))
  helper <- function(x) x + offset
  fn <- crate(
    function(i) helper(i),
    helper = helper,
    offset = 10
  )
  out <- mirai::mirai_map(1:3, fn)[.progress = FALSE]
  expect_equal(unlist(out), 11:13)
})
