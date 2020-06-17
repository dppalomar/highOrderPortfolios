rm(list = ls())
require(highOrderPortfolios)
require(magrittr)

data(X50)
N <- ncol(X50)

# estimate moment parameters
mom_params <- estMomParams(X50, align_order = TRUE)
w0 <- rep(1/N, N)
moms0 <- evalMoms(w = w0, mom_params = mom_params)
d <- abs(moms0) 

# scale the moments to align the order of magnitude
kappa <- sqrt(w0%*%mom_params$Sgm%*%w0) * 0.3

# portfolio optimization
sol_L_MVSKT <- MVSKtilting(d = d, tau_w = 30, tau_delta = 5, zeta = 1e-2, ftol = 1e-6, wtol = 1e-6, maxiter = 5e2, mom_params = mom_params, w_init = w0, w0 = w0, moms0 = moms0, kappa = kappa, method = "L-MVSKT")
sol_Q_MVSKT <- MVSKtilting(d = d, tau_w = 1e-6, tau_delta = 1e-6, zeta = 1e-10, ftol = 1e-6, wtol = 1e-6, mom_params = mom_params, w_init = w0, w0 = w0, moms0 = moms0, kappa = kappa, method = "Q-MVSKT")
t1 <- proc.time()[3]
sol_pkg <- mvskPortfolios::mvskPortfolio(m1 = mom_params$mu, M2 = mom_params$Sgm, M3 = mom_params$Phi, M4 = mom_params$Psi, w0 = w0, g = d, href = "TEvol", kappa = kappa)
time_pkg <- as.numeric(proc.time()[3] - t1)

plot(sol_L_MVSKT$times, sol_L_MVSKT$objs, type = "b")
points(sol_Q_MVSKT$times, sol_Q_MVSKT$objs, type = "h")
abline(h = -sol_pkg$delta, v = time_pkg)

rbind(sol_L_MVSKT$w, sol_Q_MVSKT$w, sol_pkg$w) %>% barplot(., beside = TRUE)
# w <- as.vector(sol_pkg$w)
# (w-w0)%*%mom_params$Sgm%*%(w-w0)
