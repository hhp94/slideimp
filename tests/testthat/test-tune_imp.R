test_that("tune_imp .f = `SlideKnn`, `knn_imp` or function works", {
  data(khanmiss1)
  SlideKnn_par <- data.frame(
    n_feat = c(100, 100),
    k = c(5, 10),
    n_overlap = c(10, 10),
    method = "euclidean",
    post_imp = FALSE
  )
  set.seed(1234)
  # Tune `SlideKnn` function on a subset of khanmiss1
  obj <- t(khanmiss1)[1:30, sample.int(nrow(khanmiss1), size = 200)]
  expect_true(anyNA(obj))

  # Check `SlideKnn`
  expect_no_error(tune_imp(obj, SlideKnn_par, rep = 1, .f = "SlideKnn"))

  # Check `knn_imp`
  knn_imp_par <- data.frame(
    k = c(5, 10),
    method = "euclidean",
    post_imp = TRUE
  )
  expect_no_error(tune_imp(obj, knn_imp_par, rep = 1, .f = "knn_imp"))

  # Check custom function
  custom_fun <- function(obj, value) {
    obj[is.na(obj)] <- value
    return(obj)
  }
  custom_par <- data.frame(
    value = c(0, 1)
  )
  expect_no_error(tune_imp(obj, custom_par, rep = 1, .f = custom_fun))
})
