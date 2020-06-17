require(magrittr)  # for dataset10

N <- 50

# load data
set.seed(3)
load("data-raw/stockdata_from_2009-01-01_to_2018-12-31_(70514247d845bb494d2269e5d5f9521e).RData")
X50 <- portfolioBacktest::stockDataResample(stockdata, N_sample = N, T_sample = N*5+1, num_datasets = 1)[[1]]$adjusted %>% log() %>% diff(., na.pad = FALSE) %>% as.matrix()


colnames(X50) <- rownames(X50) <- NULL

save(X50, file = "data/X50.RData", version = 2)