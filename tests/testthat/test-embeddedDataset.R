context("Checking embedded dataset")

data("X50")  # data in the package
# # to reduce required storage, only check sample mean and standard variance
# dataset_check <- as.vector(X50)
# dataset_features_check <- c(mean(dataset_check), sum(dataset_check), sd(dataset_check))
# save(dataset_features_check, file = "dataset_features_check.RData", version = 2)

test_that("the dataset used is the same", {
  dataset_check <- as.vector(X50)
  dataset_features <- c(mean(dataset_check), sum(dataset_check), sd(dataset_check))
  load("dataset_features_check.RData")
  expect_equivalent(dataset_features, dataset_features_check)
})
