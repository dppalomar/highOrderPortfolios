context("Checking high-order portfolio design")

data("X50")  # data in the package
N <- ncol(X50)


test_that("MVSK portfolio design", {
  X_moments <- estimate_moments(X50)
  xi <- 10
  lmd <- c(1, xi/2, xi*(xi+1)/6, xi*(xi+1)*(xi+2)/24)
  # sol_MVSK_check <- design_MVSK_portfolio(lmd = lmd, X_moments = X_moments)[-2]
  # save(sol_MVSK_check, file = "sol_MVSK_check.RData", version = 2)
  sol_MVSK <- design_MVSK_portfolio(lmd = lmd, X_moments = X_moments)[-2]
  load("sol_MVSK_check.RData")
  expect_equivalent(sol_MVSK, sol_MVSK_check)
})


test_that("MVSK tilting portfolio design", {
  X_moments <- estimate_moments(X50, adjust_magnitude = TRUE)
  w0 <- rep(1/N, N)
  w0_moments <- eval_portfolio_moments(w = w0, X_moments = X_moments)
  d <- abs(w0_moments) 
  kappa <- sqrt(w0%*%X_moments$Sgm%*%w0) * 0.3
  # sol_tilting_check <- design_MVSKtilting_portfolio(d = d, X_moments = X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)[-3]
  # save(sol_tilting_check, file = "sol_tilting_check.RData", version = 2)
  sol_tilting <- design_MVSKtilting_portfolio(d = d, X_moments = X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)[-3]
  load("sol_tilting_check.RData")
  expect_equivalent(sol_tilting, sol_tilting_check)
})