rm(list = ls())
require(highOrderPortfolios)
require(magrittr)

data(X50)

# estimate moment parameters
mom_params <- estMomParams(X50)

# lambda, moment weights
xi <- 10
lmd <- c(1, xi/2, xi*(xi+1)/6, xi*(xi+1)*(xi+2)/24)

sol_nloptr_MVSK <- highOrderPortfolios:::.MVSKnloptr(lmd = lmd, mom_params = mom_params, ftol_rel = 1e-6, xtol_rel = 1e-6, stopval = NULL)
stopval <- sol_nloptr_MVSK$obj

# portfolio optimization
sol_SCA_MVSK <- MVSK(lmd = lmd, mom_params = mom_params, ftol = 1e-6, wtol = 1e-6, method = "Q-MVSK", stopval = stopval)
sol_MM_MVSK  <- MVSK(lmd = lmd, mom_params = mom_params, method = "MM", maxiter = 2e2)
sol_DC_MVSK  <- MVSK(lmd = lmd, mom_params = mom_params, method = "DC", ftol = -1, wtol = -1, maxiter = 2e2)

require(ggplot2)
objs <- rbind(
  data.frame("times" = sol_DC_MVSK$times,  "objective" = sol_DC_MVSK$objs,  "alg" = "DC"),
  data.frame("times" = sol_MM_MVSK$times,  "objective" = sol_MM_MVSK$objs,  "alg" = "MM"),
  data.frame("times" = sol_SCA_MVSK$times, "objective" = sol_SCA_MVSK$objs, "alg" = "Q-MVSK"),
  data.frame("times" = sol_nloptr_MVSK$time, "objective" = sol_nloptr_MVSK$obj, "alg" = "nloptr")
)
ggplot(objs, aes(x = times, y = objective, color = alg, shape = alg)) +
  geom_hline(yintercept = objs[objs$alg == "nloptr", ]$objective, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = sol_nloptr_MVSK$time, linetype = "dashed", color = "gray30") +
  geom_point(size = 2) +
  geom_line() +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #               labels = trans_format("log10", math_format(10^.x))) +
  scale_shape_manual(values = c(25, 24, 23, 19)) +
  labs(title = NULL, x = "CPU time (seconds)", y = "Objective") + 
  theme(legend.position = c(0.8, 0.4),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 8))
