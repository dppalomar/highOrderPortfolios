context("Checking embedded dataset")
#library(testthat)

data("X50")  # data in the package
# # to reduce required storage, only check sample mean and standard variance
# X50_vector <- as.vector(X50)
# dataset_features_check <- c(mean(X50_vector), sum(X50_vector), sd(X50_vector))
# save(dataset_features_check, file = "dataset_features_check.RData", version = 2)

test_that("the dataset used is the same", {
  X50_vector <- as.vector(X50)
  dataset_features <- c(mean(X50_vector), sum(X50_vector), sd(X50_vector))
  load("dataset_features_check.RData")
  expect_equivalent(dataset_features, dataset_features_check)
})
