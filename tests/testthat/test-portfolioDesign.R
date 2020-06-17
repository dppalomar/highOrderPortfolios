context("Checking high-order portfolio design")

data("X50")  # data in the package
N <- ncol(X50)


test_that("MVSK portfolio design", {
  params <- estMomParams(X50)
  xi <- 10
  lmd <- c(1, xi/2, xi*(xi+1)/6, xi*(xi+1)*(xi+2)/24)
  # sol_MVSK_check <- MVSK(lmd = lmd, mom_params = params)[-2]
  # save(sol_MVSK_check, file = "sol_MVSK_check.RData", version = 2)
  sol_MVSK <- MVSK(lmd = lmd, mom_params = params)[-2]
  load("sol_MVSK_check.RData")
  expect_equivalent(sol_MVSK, sol_MVSK_check)
})


test_that("MVSK tilting portfolio design", {
  params <- estMomParams(X50, align_order = TRUE)
  w0 <- rep(1/N, N)
  moms0 <- evalMoms(w = w0, mom_params = params)
  d <- abs(moms0) 
  kappa <- sqrt(w0%*%params$Sgm%*%w0) * 0.3
  # sol_tilting_check <- MVSKtilting(d = d, mom_params = params, w_init = w0, w0 = w0, moms0 = moms0, kappa = kappa)[-3]
  # save(sol_tilting_check, file = "sol_tilting_check.RData", version = 2)
  sol_tilting <- MVSKtilting(d = d, mom_params = params, w_init = w0, w0 = w0, moms0 = moms0, kappa = kappa)[-3]
  load("sol_tilting_check.RData")
  expect_equivalent(sol_tilting, sol_tilting_check)
})