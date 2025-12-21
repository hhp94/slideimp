## code to prepare `khanmiss1` dataset goes here
library(impute)
data(khanmiss)
# khanmiss has invalid names
khanmiss
khanmiss1 <- khanmiss[2:nrow(khanmiss), 3:ncol(khanmiss)]
khanmiss1 <- as.matrix(
  dplyr::mutate(
    khanmiss1,
    dplyr::across(.cols = dplyr::everything(), .fns = as.numeric)
  )
)
colnames(khanmiss1) <- names(khanmiss)[3:ncol(khanmiss)]
anyDuplicated(khanmiss$X[2:nrow(khanmiss)])
khanmiss$X <- c("", paste0("g", seq_len(nrow(khanmiss) - 1)))
rownames(khanmiss1) <- khanmiss$X[2:nrow(khanmiss)]
anyDuplicated(rownames(khanmiss1))
is.na(khanmiss1) |>
  rowSums() |>
  sort(decreasing = TRUE)
usethis::use_data(khanmiss1, overwrite = TRUE)
