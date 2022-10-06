context("Checking sample moments parameter estimation and calculation")
#library(testthat)

data("X50")  # data in the package
N <- ncol(X50)



test_that("skew t parameter estimation via function \"fit_ghMST()\" and portfolio moments conincide with the precomputed ones", {
  ghMST_model <- estimate_skew_t(X50, nu_lb = 9, PXEM = FALSE, max_iter = 10000)
  load("fit_ghMST_check.RData")
  expect_equal(ghMST_model[c("mu", "gamma", "scatter", "nu")],
               ghMST_model_check[c("mu", "gamma", "scatter", "nu")],
               tolerance = 1e-5)
  
  
  
  # compare parametric vs nonparametric
  moments_param <- eval_portfolio_moments(w = rep(1/N, N), X_statistics = ghMST_model)
  X_moments <- estimate_sample_moments(X50)
  moments_nonparam <- eval_portfolio_moments(w = rep(1/N, N), X_statistics = X_moments)
  expect_equal(moments_nonparam, moments_param, tolerance = 1e-5)
})
