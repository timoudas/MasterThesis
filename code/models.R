library(tsDyn)
library(vars)
source('stationarity_tests.R')

var = VAR(model_data, lag=4, type="const")
tvar <- TVAR(model_data, lag=4, nthresh=1, thDelay=2, trim = 0.10, plot=TRUE)
AIC(var)
AIC(tvar)


