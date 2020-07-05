context("Checking moments parameter estimation and calculation")

data("X50")  # data in the package
N <- ncol(X50)
w <- rep(1/N, N)

test_that("parameter estimation and moments conincide with the precomputed ones", {
  # params <- estMomParams(X50)
  # moments_check <- evalMoms(w = w, mom_params = params)
  # save(moments_check, file = "moments_check.RData", version = 2)
  X_moments <- estimate_moments(X50)
  moments <- eval_portfolio_moments(w = w, X_moments = X_moments)
  load("moments_check.RData")
  expect_equivalent(moments, moments_check)
})