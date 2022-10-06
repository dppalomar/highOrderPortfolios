test_that("Function \"fit_ghMST()\"", {
  ghMST_model <- estimate_skew_t(X50, nu_lb = 9, PXEM = FALSE, max_iter = 10000)
  load("fit_ghMST_check.RData")
  expect_equal(ghMST_model[c("mu", "gamma", "scatter", "nu")],
               ghMST_model_check[c("mu", "gamma", "scatter", "nu")],
               tolerance = 1e-5)
})
