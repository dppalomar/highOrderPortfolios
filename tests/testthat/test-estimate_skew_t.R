context("Checking sample moments parameter estimation and calculation")
#library(testthat)

data("X50")  # data in the package
N <- ncol(X50)
 
test_that("skew t parameter estimation via function \"estimate_skew_t()\" and portfolio moments conincide with the precomputed ones", {
  # ghMST_model_check <- estimate_skew_t(X50, nu_lb = 9, PXEM = FALSE, max_iter = 10000)
  # save(ghMST_model_check, file = "estimate_skew_t.RData", version = 2)
  
  ghMST_model <- estimate_skew_t(X50, nu_lb = 9, PXEM = FALSE, max_iter = 10000)
  load("estimate_skew_t.RData")
  expect_equal(ghMST_model[c("mu", "gamma", "scatter", "nu")],
               ghMST_model_check[c("mu", "gamma", "scatter", "nu")],
               tolerance = 1e-5)
    
  # compare parametric vs nonparametric
  moments_param <- eval_portfolio_moments(w = rep(1/N, N), X_statistics = ghMST_model)
  X_moments <- estimate_sample_moments(X50)
  moments_nonparam <- eval_portfolio_moments(w = rep(1/N, N), X_statistics = X_moments) 
  
  expect_equal(moments_nonparam, moments_param, tolerance = 1e-3)
})
