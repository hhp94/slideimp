#' Print ImputedMatrix
#'
#' @param x An ImputedMatrix
#' @param n Number of rows to print
#' @param m Number of cols to print
#' @param ... Not used
#'
#' @returns Invisible object of class ImputedMatrix
#' @export
#'
#' @examples
#' data(khanmiss1)
#' t_khanmiss1 <- t(khanmiss1)
#' result <- knn_imp(t_khanmiss1, k = 5)
#' print(result, n = 6, m = 6)
print.ImputedMatrix <- function(x, n = 5, m = 5, ...) {
  imp_method <- attr(x, "imp_method")

  # Cleaner header
  cat("ImputedMatrix (", toupper(imp_method), ")\n", sep = "")
  cat("Dimensions: ", nrow(x), " x ", ncol(x), "\n\n", sep = "")

  # Print subset
  subset_x <- x[seq_len(min(n, nrow(x))), seq_len(min(m, ncol(x))), drop = FALSE]
  print(subset_x, ...)

  # Better truncation message
  if (n < nrow(x) || m < ncol(x)) {
    cat("\n# Showing [1:", min(n, nrow(x)), ", 1:", min(m, ncol(x)),
      "] of full matrix\n",
      sep = ""
    )
  }

  invisible(x)
}
