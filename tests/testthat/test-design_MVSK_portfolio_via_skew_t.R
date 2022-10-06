test_that("Checking high-order design via parametric approach", {
  X <- X50[, 1:10]
  w0 <- rep(1/10, 10)
  X_parameters <- estimate_skew_t(X)

  lambda <- c(1, 4, 10, 20)
  load("MVSK_skew_t.RData")

  L_MVSK_result <- design_MVSK_portfolio_via_skew_t(lambda, X_parameters, method = "L-MVSK", maxiter = 10000, ftol = 1e-10, wtol = 1e-10)
  expect_equal(L_MVSK_result$w, w_check, tolerance = 1e-4)

  Q_MVSK_result <- design_MVSK_portfolio_via_skew_t(lambda, X_parameters, method = "Q-MVSK", maxiter = 10000, ftol = 1e-10, wtol = 1e-10)
  expect_equal(Q_MVSK_result$w, w_check, tolerance = 1e-4)

  DC_result <- design_MVSK_portfolio_via_skew_t(lambda, X_parameters, method = "DC", maxiter = 10000, ftol = 1e-10, wtol = 1e-10)
  expect_equal(DC_result$w, w_check, tolerance = 1e-4)

  PGD_result <- design_MVSK_portfolio_via_skew_t(lambda, X_parameters, method = "PGD", maxiter = 10000, ftol = 1e-10, wtol = 1e-10, initial_eta = 100)
  expect_equal(PGD_result$w, w_check, tolerance = 1e-4)

  RFPA_result <- design_MVSK_portfolio_via_skew_t(lambda, X_parameters, method = "RFPA", maxiter = 10000, ftol = 1e-10, wtol = 1e-10, tau = 100)
  expect_equal(RFPA_result$w, w_check, tolerance = 1e-4)
})
