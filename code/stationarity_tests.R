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
ADF.test(DATA_CLIF[, 'CLIF']) #STATIONARY
#BNP
ADF.test(DATA_CLIF[, 'BNP']) #STATIONARY
#INFLATION
ADF.test(DATA_CLIF[, 'CPI']) #STATIONARY
#INTEREST RATE
ADF.test(DATA_CLIF[, 'INT']) #UNIT ROOT
summary(ur.kpss(DATA_CLIF[, 'INT'])) #NOT STATIONARY

##### FIRST DIFFERENCE #####
ADF.test(diff(DATA_CLIF[, 'INT'])) #STATIONARY
summary(ur.kpss(diff(DATA_CLIF[, 'INT']))) #STATIONARY



##### IN LEVELS #####
#RIKSBANKEN FSI INDEX
plot(DATA_SFSI[, 'SFSI'])
ADF.test(DATA_SFSI[, 'SFSI']) #NOT STATIONARY
summary(ur.kpss(DATA_SFSI[, 'SFSI'])) #STATIONARY
#BNP
ADF.test(DATA_SFSI[, 'BNP']) #STATIONARY
#INFLATION
ADF.test(DATA_SFSI[, 'CPI']) #STATIONARY
#INTEREST RATE
ADF.test(DATA_SFSI[, 'INT']) #UNIT ROOT
summary(ur.kpss(DATA_SFSI[, 'INT'])) #NOT STATIONARY

##### FIRST DIFFERENCE #####
ADF.test(diff(DATA_SFSI[, 'INT'])) #STATIONARY
summary(ur.kpss(diff(DATA_SFSI[, 'INT']))) #STATIONARY

##### IN LEVELS #####
#ECB FSI INDEX
ADF.test(DATA_ECB[, 'ECBFSI']) #NOT STATIONARY
summary(ur.kpss(DATA_ECB[, 'ECBFSI'])) #STATIONARY
#BNP
ADF.test(DATA_ECB[, 'BNP']) #STATIONARY
#INFLATION
ADF.test(DATA_ECB[, 'CPI']) #STATIONARY
#INTEREST RATE
ADF.test(DATA_ECB[, 'INT']) #UNIT ROOT
summary(ur.kpss(DATA_ECB[, 'INT'])) #NOT STATIONARY

##### FIRST DIFFERENCE #####
ADF.test(diff(DATA_ECB[, 'INT'])) #STATIONARY
summary(ur.kpss(diff(DATA_ECB[, 'INT']))) #STATIONARY

model_data.SFI <- ts.intersect(DATA_SFSI[, 'SFSI'], DATA_SFSI[, 'BNP'], DATA_SFSI[, 'CPI'], diff(DATA_SFSI[, 'INT']))
colnames(model_data.SFI) <- c('SFSI', 'BNP', 'CPI', 'DINT')
                        

model_data.ECB <- ts.intersect(DATA_ECB[, 'ECBFSI'], DATA_ECB[, 'BNP'], DATA_ECB[, 'CPI'], diff(DATA_ECB[, 'INT']))
colnames(model_data.ECB) <- c('ECBFSI', 'BNP', 'CPI', 'DINT')

model_data.CLIFS <- ts.intersect(DATA_CLIF[, 'CLIF'], DATA_CLIF[, 'BNP'], DATA_CLIF[, 'CPI'], diff(DATA_CLIF[, 'INT']))
colnames(model_data.CLIFS) <- c('CLIFS', 'BNP', 'CPI', 'DINT')


plot(model_data.SFI[,1], ylab='SFSI')


summary(model_data.SFI)
plot(model_data.SFI[,'SFSI'])
hist(model_data.SFI[,'SFSI'], breaks=40)

