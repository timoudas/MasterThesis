library(tsDyn)
library(vars)
library(aod)
source('stationarity_tests.R')

#Modeldata
model_data <- window(ts.intersect(BNP, CPI1,  diff(INT1), CLIF_qts), start=c(1993,2))
colnames(model_data) <-c("BNP", 'CPI', "DINT", "CLIF")

#AIC(4), BIC(1)
VARselect(model_data, lag.max = 10, type = 'const')$selection


TVAR.LRtest(model_data, lag = 4, thDelay = 3, trim = 0.15, mTh = 4, nboot=50, plot=TRUE)


tvar <- TVAR(model_data, lag=4, nthresh=1, thDelay=3, trim = 0.15, mTh = 4, plot=TRUE, trace = TRUE)

irf_th1_l2_L <- irf(tvar, regime = "Bup")
irf_th1_l2_H <- irf(tvar, regime = "Bdown")

wald.test(Sigma = vcov(tvar), b = coef(tvar)$Bup, Terms = 1:2)


plot(model_data[,'CLIF'])
c(AIC(tvar), BIC(tvar), logLik(tvar))

AIC(var)
AIC(tvar)


