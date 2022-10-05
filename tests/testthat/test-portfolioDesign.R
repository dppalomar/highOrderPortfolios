context("Checking high-order portfolio design")
#library(testthat)

data("X50")  # data in the package
N <- ncol(X50)


test_that("MVSK portfolio design", {
  X_moments <- estimate_sample_moments(X50)
  xi <- 10
  lmd <- c(1, xi/2, xi*(xi+1)/6, xi*(xi+1)*(xi+2)/24)
  
  # Method "Q-MVSK"
  # sol_MVSK_check <- design_MVSK_portfolio_via_sample_moments(lmd = lmd, X_moments = X_moments, method = "Q-MVSK")
  # save(sol_MVSK_check, file = "sol_MVSK_QMVSK_check.RData", version = 2)
  sol_MVSK <- design_MVSK_portfolio_via_sample_moments(lmd = lmd, X_moments = X_moments)
  load("sol_MVSK_QMVSK_check.RData")
  expect_equivalent(sol_MVSK[-2], sol_MVSK_check[-2])
  

  # Method "MM"
  # sol_MVSK_check <- design_MVSK_portfolio_via_sample_moments(lmd = lmd, X_moments = X_moments, method = "MM")
  # save(sol_MVSK_check, file = "sol_MVSK_MM_check.RData", version = 2)
  sol_MVSK <- design_MVSK_portfolio_via_sample_moments(lmd = lmd, X_moments = X_moments, method = "MM")
  load("sol_MVSK_MM_check.RData")
  expect_equivalent(sol_MVSK[-2], sol_MVSK_check[-2])
  
  
  # Method "DC"
  # sol_MVSK_check <- design_MVSK_portfolio_via_sample_moments(lmd = lmd, X_moments = X_moments, method = "DC")
  # save(sol_MVSK_check, file = "sol_MVSK_DC_check.RData", version = 2)
  sol_MVSK <- design_MVSK_portfolio_via_sample_moments(lmd = lmd, X_moments = X_moments, method = "DC")
  load("sol_MVSK_DC_check.RData")
  expect_equivalent(sol_MVSK[-2], sol_MVSK_check[-2])  
})





test_that("MVSK tilting portfolio design", {
  X_moments <- estimate_sample_moments(X50, adjust_magnitude = TRUE)
  w0 <- rep(1/N, N)
  w0_moments <- eval_portfolio_moments(w = w0, X_statistics = X_moments)
  d <- abs(w0_moments) 
  kappa <- 0.3 * sqrt(w0 %*% X_moments$Sgm %*% w0)
  
  
  # Method "L-MVSKT"
  # sol_tilting_check <- design_MVSKtilting_portfolio_via_sample_moments(method = "L-MVSKT", d = d, X_moments = X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)
  # save(sol_tilting_check, file = "sol_MVSKtilting_LMVSKT_check.RData", version = 2)
  sol_tilting <- design_MVSKtilting_portfolio_via_sample_moments(method = "L-MVSKT", d = d, X_moments = X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)
  load("sol_MVSKtilting_LMVSKT_check.RData")
  expect_equivalent(sol_tilting[-3], sol_tilting_check[-3])

  
  # Method "Q-MVSKT"
  # sol_tilting_check <- design_MVSKtilting_portfolio_via_sample_moments(method = "Q-MVSKT", d = d, X_moments = X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)
  # save(sol_tilting_check, file = "sol_MVSKtilting_QMVSKT_check.RData", version = 2)
  sol_tilting <- design_MVSKtilting_portfolio_via_sample_moments(method = "Q-MVSKT", d = d, X_moments = X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)
  load("sol_MVSKtilting_QMVSKT_check.RData")
  expect_equivalent(sol_tilting[-3], sol_tilting_check[-3])
})
