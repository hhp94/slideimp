skip_if_not_manual()

# -----------------------------------------------------------------------------
# Reduced base-R re-implementation of carrier::crate()
#
# Derived from carrier: https://github.com/r-lib/carrier
# Copyright (c) Posit Software, PBC
# Original implementation by the carrier authors
# Licensed under the MIT License (see inst/LICENSE.note)
# -----------------------------------------------------------------------------

test_that(".crate() works through mirai", {
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
