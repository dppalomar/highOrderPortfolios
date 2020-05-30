library(portfolioBacktest)  # for dataset10

X50 <- diff(log(dataset10[[1]]$adjusted))[-1]

save(X50, file = "data/X50.RData", version = 2)
