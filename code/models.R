# TVAR MODEL STEPS
# MODEL SPECIFICATION: CLIF, BNP, INF, INT
# CHECK FOR STATIONARITY - DIFF NON STATIONARY VARIABLES
# OPTIMAL LAG SELECTION: AIC/BIC for VAR-MODEL
# CHECK FOR AUTOCORRELATION - IN VAR MODEL / IF NONE CONTINUE / ELSE ADD LAG
# THRESHOLD TEST - TSAR TEST / CHECKS AGAINST VAR(p) MODEL
"The lag order p and the threshold lag d need to be determined a priori,
which in case of p is achieved by applying the normal information criteria in the
linear VAR estimation. For the choice of d we will rely on economic reasoning."
# FIND THRESHOLD - GRID SEARCH
# ESTIMATE MODEL
# GIRF
#
#
#
#
#

library(tsDyn)
library(vars)
source('stationarity_tests.R')
library(remotes)
suppressMessages(library(nonlinearTseries))
source("https://raw.githubusercontent.com/MatthieuStigler/tsDyn_GIRF/master/GIRF2")
#install_github("angusmoore/tvarGIRF", ref = "v0.1.2")
#library(tvarGIRF)
dev.off()
set.seed(1234)

norm.thresh <- list()
norm.thresh$CLIF <- (model_data.CLIFS[, 'CLIFS'] - mean(model_data.CLIFS[, 'CLIFS'])) / sd(model_data.CLIFS[, 'CLIFS'])
norm.thresh$SFSI <- (model_data.SFI[, 'SFSI'] - mean(model_data.SFI[, 'SFSI'])) / sd(model_data.SFI[, 'SFSI'])
norm.thresh$ECB <- (model_data.ECB[, 'ECBFSI'] - mean(model_data.ECB[, 'ECBFSI'])) / sd(model_data.ECB[, 'ECBFSI'])

model_data.CLIFS[, 'CLIFS'] <-  (model_data.CLIFS[, 'CLIFS'] - mean(model_data.CLIFS[, 'CLIFS'])) / sd(model_data.CLIFS[, 'CLIFS'])
model_data.SFI[, 'SFSI'] <- (model_data.SFI[, 'SFSI'] - mean(model_data.SFI[, 'SFSI'])) / sd(model_data.SFI[, 'SFSI'])
model_data.ECB[, 'ECBFSI'] <- (model_data.ECB[, 'ECBFSI'] - mean(model_data.ECB[, 'ECBFSI'])) / sd(model_data.ECB[, 'ECBFSI'])




#LAG SELECTION FOR BASELINE VAR MODEL
VARselect(model_data.CLIFS, lag.max = 10, type = 'const')$selection #AIC(4), BIC(1)
VARselect(model_data.SFI, lag.max = 10, type = 'const')$selection #AIC(4), BIC(1)
VARselect(model_data.ECB, lag.max = 10, type = 'const')$selection #AIC(10), BIC(1)

#ESTIMATE VAR MODEL
var.clif <- VAR(model_data.CLIFS, p=1, type='const')
var.sfi <- VAR(model_data.SFI, p=1, type='const')
var.ecb <- VAR(model_data.ECB, p=1, type='const')

#TEST FOR AUTOCORRELATION
serial.test(var.clif) #NO AUTOCORR
serial.test(var.sfi) #NO AUTOCORR
serial.test(var.ecb) #NO AUTOCORR
#ROOTS
roots(var.clif) #STABLE ROOTS
roots(var.sfi) #STABLE ROOTS
roots(var.ecb) #STABLE ROOTS

#PLOT THESH_VARIABLES
plot(model_data.CLIFS[, 'CLIFS'])
plot(model_data.SFI[, 'SFSI'])
plot(model_data.ECB[, 'ECBFSI'])


#TSAY TEST CLIF
tsayTest(model_data.CLIFS[, 'CLIFS'], 1) #DOES NOT FOLLOW AN AR(1)
tsayTest(model_data.CLIFS[, 'CLIFS'], 2) #FOLLOW AN AR(2)
tsayTest(model_data.CLIFS[, 'CLIFS'], 3) #DOES NOT FOLLOW AN AR(1)

#TSAY TEST SFSI
model_data.SFI
tsayTest(model_data.SFI[, 'SFSI'], 1) #FOLLOWS AN AR(1)
tsayTest(model_data.SFI[, 'SFSI'], 2) #FOLLOWS AN AR(2)
tsayTest(model_data.SFI[, 'SFSI'], 3) #DOES NOT FOLLOW AN AR(3)

#TSAY TEST ECB
tsayTest(model_data.ECB[, 'ECBFSI'], 1) #FOLLOW AN AR(1)
tsayTest(model_data.ECB[, 'ECBFSI'], 2) #FOLLOW AN AR(2)
tsayTest(model_data.ECB[, 'ECBFSI'], 3) #FOLLOW AN AR(3)

#THRESHOLD TEST
thresholdTest(model_data.SFI[, 'SFSI'], p = 1, d = 1, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.SFI[, 'SFSI'], p = 1, d = 2, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.SFI[, 'SFSI'], p = 1, d = 3, lower.percent = 0.15, upper.percent = 0.85) #NS

thresholdTest(model_data.SFI[, 'SFSI'], p = 4, d = 1, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.SFI[, 'SFSI'], p = 4, d = 2, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.SFI[, 'SFSI'], p = 4, d = 3, lower.percent = 0.15, upper.percent = 0.85) #S

#CLIFS INDEX
thresholdTest(model_data.CLIFS[, 'CLIFS'], p = 1, d = 1, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.CLIFS[, 'CLIFS'], p = 1, d = 2, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.CLIFS[, 'CLIFS'], p = 1, d = 3, lower.percent = 0.15, upper.percent = 0.85) #NS

thresholdTest(model_data.CLIFS[, 'CLIFS'], p = 4, d = 1, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.CLIFS[, 'CLIFS'], p = 4, d = 2, lower.percent = 0.15, upper.percent = 0.85) #NS
thresholdTest(model_data.CLIFS[, 'CLIFS'], p = 4, d = 3, lower.percent = 0.15, upper.percent = 0.85) #NS

#ECB INDEX
thresholdTest(model_data.ECB[, 'ECBFSI'], p = 1, d = 1, lower.percent = 0.15, upper.percent = 0.85) #NS


###LIKELIHOOD RATIO TEST

#CLIFS INDEX
TVAR.LRtest(model_data.CLIFS, lag = 1, thDelay = 1, model = "TAR", trim = 0.15, mTh = 1, nboot=400, plot=TRUE)

#ECB
TVAR.LRtest(model_data.ECB, lag = 1, thDelay = 1, model = "TAR", trim = 0.15, mTh = 1, nboot=200, plot=TRUE)

#SFSI SIGNIFICATLY DIFFERENT FROM LINEAR VAR
TVAR.LRtest(model_data.SFI, lag = 1, thDelay = 1, model = "TAR", trim = 0.15, mTh = 1, nboot=500, plot=TRUE)


tvar.sfsi <- TVAR(model_data.SFI, lag=3, thDelay = 2, nthresh=1, trim = 0.15, mTh = 1, plot=TRUE, trace = TRUE)

md <- model_data.SFI[,2:4]
th <- model_data.SFI[, 1]

tvar <- TVAR(md, lag=1, nthresh=1, thDelay=1, thVar = th, trim = 0.15, plot=TRUE, trace = TRUE)

e <- resid(tvar)
cov_matrix <- t(e) %*% e
R <- chol(cov_matrix)  ## upper tri
L <-t(R)  ## lower tri
L

tvar <- TVAR(model_data.SFI, lag=1, nthresh=1, thDelay=1, trim = 0.15, mTh = 1, plot=TRUE, trace = TRUE)



source("tar.R")
t <- unclass(model_data.SFI[, 'SFSI'])
summary(t > 0.2694938)
out1 <- tar(diff(t), 4, omit = 1, r1 = 0.15, r2 = 0.85, rep = 1000, 1, 1)


plot(model_data.SFI[, 'SFSI'], at=c(c(1995, 2):c(2022, 1)), labels=c(1995,2):c(2022, 1))
abline(h=0.2694938, col="blue")



e <- resid(tvar.sfsi)
cov_matrix <- t(e) %*% e
R <- chol(cov_matrix)  ## upper tri
L <-t(R)  ## lower tri
L[2,]

resGIRF <- GIRF(tvar.sfsi)
print(resGIRF)

data(zeroyld)

#par(mar = rep(2, 4))

B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),
         c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
colnames(B) <- paste(rep(c("Const", "Lag_1_var1", "Lag_1_var2"), 2), c("Low", "High"), sep="_")
TVAR.sim(B=B,nthresh=1,n=500, mTh=1, Thresh=5, starting=matrix(c(5.2, 5.5), nrow=1))
B




