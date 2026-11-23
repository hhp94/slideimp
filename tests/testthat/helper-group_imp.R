# The text of every rlang_message `expr` raises, muffled, so that a test can
# pin one message's wording and class among the others a call reports. A
# message of another class is not collected, and reaches the caller.
rlang_messages <- function(expr) {
  msgs <- character()
  withCallingHandlers(
    expr,
    rlang_message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  msgs
}

# One matrix for the post_imp and mean-fill tests of group_imp(), in
# test-group_imp.R and test-group_imp-manual.R. Group "b" is columns 7:12, so
# its column over colmax, f9, is local column 3 and column 9 of `obj`.
post_imp_case <- function() {
  set.seed(1234)
  obj <- matrix(
    rnorm(20 * 12),
    20,
    12,
    dimnames = list(NULL, paste0("f", 1:12))
  )
  obj[1:3, 2] <- NA # group "a": ordinary missingness
  obj[2:20, 9] <- NA # group "b": miss_rate 0.95, above the default colmax
  obj[4:5, 11] <- NA # group "b": ordinary missingness
  list(
    obj = obj,
    group = data.frame(
      feature = colnames(obj),
      group = rep(c("a", "b"), c(6, 6))
    )
  )
}
