rm(list = ls())
require(highOrderPortfolios)
require(magrittr)

data(X50)
N <- ncol(X50)

# estimate moment parameters
mom_params <- estimate_moments(X50, adjust_magnitude = TRUE)
w0 <- rep(1/N, N)
moms0 <- eval_portfolio_moments(w = w0, X_moments = mom_params)
d <- abs(moms0) 

# scale the moments to align the order of magnitude
kappa <- sqrt(w0%*%mom_params$Sgm%*%w0) * 0.3

# portfolio optimization
sol_L_MVSKT <- design_MVSKtilting_portfolio(d = d, tau_w = 30, tau_delta = 5, zeta = 1e-2, ftol = 1e-6, wtol = 1e-6, maxiter = 5e2, X_moments = mom_params, w_init = w0, w0 = w0, w0_moments = moms0, kappa = kappa, method = "L-MVSKT")
sol_Q_MVSKT <- design_MVSKtilting_portfolio(d = d, tau_w = 1e-6, tau_delta = 1e-6, zeta = 1e-10, ftol = 1e-6, wtol = 1e-6, X_moments = mom_params, w_init = w0, w0 = w0, w0_moments = moms0, kappa = kappa, method = "Q-MVSKT")
t1 <- proc.time()[3]
sol_pkg <- mvskPortfolios::mvskPortfolio(m1 = mom_params$mu, M2 = mom_params$Sgm, M3 = mom_params$Phi, M4 = mom_params$Psi, w0 = w0, g = d, href = "TEvol", kappa = kappa)
time_pkg <- as.numeric(proc.time()[3] - t1)

plot(sol_L_MVSKT$cpu_time, sol_L_MVSKT$objs, type = "b")
points(sol_Q_MVSKT$cpu_time, sol_Q_MVSKT$objs, type = "h")
abline(h = -sol_pkg$delta, v = time_pkg)

rbind(sol_L_MVSKT$w, sol_Q_MVSKT$w, sol_pkg$w) %>% barplot(., beside = TRUE)
# w <- as.vector(sol_pkg$w)
# (w-w0)%*%mom_params$Sgm%*%(w-w0)





library(highOrderPortfolios)
data(X50)

# estimate moments
X_moments <- estimate_moments(X50)

# decide problem setting
w0 <- rep(1/50, 50)
w0_moments <- eval_portfolio_moments(w0, X_moments)
d <- abs(w0_moments) 
kappa <- sqrt(w0%*%X_moments$Sgm%*%w0) * 0.3

# portfolio optimization
sol <- design_MVSKtilting_portfolio(d, X_moments, w_init = w0, w0 = w0, w0_moments = w0_moments, kappa = kappa)
sol

save(sol, file = "R_buildignore/mvsk_tilting_sanity.RData")


