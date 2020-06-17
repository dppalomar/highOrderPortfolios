context("Checking moments parameter estimation and calculation")

data("X50")  # data in the package
N <- ncol(X50)
w <- rep(1/N, N)

test_that("parameter estimation and moments conincide with the precomputed ones", {
  # params <- estMomParams(X50)
  # moments_check <- evalMoms(w = w, mom_params = params)
  # save(moments_check, file = "moments_check.RData", version = 2)
  params <- estMomParams(X50)
  moments <- evalMoms(w = w, mom_params = params)
  load("moments_check.RData")
  expect_equivalent(moments, moments_check)
})