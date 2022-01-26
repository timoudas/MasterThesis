#Augmented dicky-fuller testing
source('data_import.R')
source("adf.R")
library(urca)

"
ADF-test: Test to check if series has a stochasic trend
Null: Series has a stochastic trend
Alt: Series is stationary
Con: Low power, when (1-psi)=theta, theta close to zero, i.e drift is high

KPSS-test: Test to check is series is stationatry
Null: Series is stationary
Alt: Series is not stationary
Usage: Use to validate ADF test if drift is high, and ADF null is not rejected
"

##### IN LEVELS #####
#CLIF INDEX
ADF.test(sample[, 'CLIF']) #STATIONARY
#BNP
ADF.test(sample[, 'BNP']) #STATIONARY
#INFLATION
ADF.test(sample[, 'CPI']) #STATIONARY
#INTEREST RATE
ADF.test(sample[, 'INT']) #UNIT ROOT
summary(ur.kpss(sample[, 'INT'])) #NOT STATIONARY

##### FIRST DIFFERENCE #####
ADF.test(diff(sample[, 'INT'])) #STATIONARY
summary(ur.kpss(diff(sample[, 'INT']))) #STATIONARY


