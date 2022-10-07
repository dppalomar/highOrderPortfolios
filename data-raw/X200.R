library(portfolioBacktest)
SP500 <- stockDataDownload(stock_symbols = SP500_symbols,
                           from = "2017-01-01", to = "2021-12-31")
Return <- diff(log(SP500$close), lag=1)[-1]

N <- 100
set.seed(2200)
data_matrix <- Return
data_matrix <- data_matrix[ , colSums(is.na(data_matrix)) == 0]
data_matrix <- data_matrix[ , sample(ncol(data_matrix), N)]

X200 <- data_matrix[1:1000, ]

save(X200, file = "data/X200.RData", version = 2)