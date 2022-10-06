context("Checking sample moments parameter estimation and calculation")
#library(testthat)

data("X50")  # data in the package
N <- ncol(X50)

test_that("sample parameter estimation and portfolio moments conincide with the precomputed ones", {
  # X_moments <- estimate_moments(X50)
  # moments_check <- eval_portfolio_moments(w = rep(1/N, N), X_moments = X_moments)
  # save(moments_check, file = "moments_check.RData", version = 2)
  X_moments <- estimate_sample_moments(X50)
  moments <- eval_portfolio_moments(w = rep(1/N, N), X_statistics = X_moments)
  load("moments_check.RData")
  expect_equivalent(moments, moments_check)
  
  X_moments_adjusted <- estimate_sample_moments(X50, adjust_magnitude = TRUE)
  moments_adjusted <- eval_portfolio_moments(w = rep(1/N, N), X_statistics = X_moments_adjusted)
  expect_equivalent(moments_adjusted, c(1, 1, -1, 1))
})
