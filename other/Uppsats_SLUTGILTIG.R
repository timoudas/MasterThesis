library(quantmod)
library(urca)
library(vars)
library(seasonal)
library(xtable)
library(tseries)
library(sjPlot)
library(Hmisc)
library(psych)

##ADF-test function
globalVariables(c("ur.df"))
ADF.test <- function(data){
  ADF <- function(type, data) {
    result1 <- ur.df(data,
                     type = type,
                     lags = 3*frequency(data),
                     selectlags = "AIC")
    result2 <- cbind(t(result1@teststat),
                     result1@cval,
                     coefficients(result1@testreg)["z.lag.1", 1], result1@lags)
    round(result2, 2)
  }
  types   <- c("trend", "drift", "none")
  result3 <- apply(t(types), 2, ADF, data)
  cat(rep("#", 20),'\n')
  cat(rep("#", 20),'\n')
  cat("Augmented Dickey--Fuller test\n")
  cat(rep("#", 20),'\n')
  cat("type:", "  trend ", "drift ", "none\n")
  cat("AR1:   ", result3[[1]][1, 5] + 1, " ", result3[[2]][1, 5] + 1, " ", result3[[3]][5] + 1, "\n")
  cat("lags:  ", result3[[1]][1, 6], "   ", result3[[2]][1, 6], "   ", result3[[3]][6], "\n")
  cat(rep("#", 20),'\n')
  result5 <- rbind(result3[[1]][c(1,3),1:4],
                   result3[[2]][1:2,1:4],
                   result3[[3]][1:4])
  rownames(result5)[5] <- "tau1" 
  result5
}

#Constant Price Gross Domestic Product in Sweden, Seasonally Adjusted
#Percent Change from Year Ago, Quarterly, Seasonally Adjusted
getSymbols("SWEGDPRQPSMEI", src = "FRED")
BNP <- ts(as.ts(SWEGDPRQPSMEI), start = c(1961, 1), frequency = 4)

#Unemployment Rate: Aged 15-64: All Persons for Sweden
#Percent, Quarterly, Seasonally Adjusted
getSymbols("LRUN64TTSEQ156S", src = "FRED")
U <- ts(as.ts(LRUN64TTSEQ156S), start = c(2001, 1), frequency = 4)

 
#Consumer Price Index: Total All Items for Sweden Growth Rate Previous Period,
#Not Seasonally Adjusted
getSymbols("CPALTT01SEQ657N", src = "FRED")
CPI <- ts(as.ts(CPALTT01SEQ657N), start = c(1960, 1), frequency = 4)

#Global price of Brent Crude
#U.S. Dollars per Barrel, Not Seasonally Adjusted
getSymbols("POILWTIUSDQ", src = "FRED")
OIL <- ts(as.ts(POILWTIUSDQ), start = c(1990, 1), frequency = 4)
OIL <- log(OIL)


#3-Month or 90-day Rates and Yields: Treasury Securities for Sweden
#Percent,Not Seasonally Adjusted
getSymbols("IR3TTS01SEQ156N", src = "FRED")
INT <- ts(as.ts(IR3TTS01SEQ156N), start = c(1982, 1), frequency = 4)

##SEASON AJDUSt##

# seasonal adjust time series using X11. 
CPI1 <- final(seas(as.ts((CPI),freq=4)))
OIL1 <- final(seas(as.ts((OIL),freq=4)))


sweden.data <- ts.intersect(OIL, U, BNP, CPI, INT)
#OIL & CPI S.ADJ
sweden.data1 <- ts.intersect(OIL1, U, BNP, CPI1, INT)
plot(sweden.data)
plot(sweden.data1)
# save(okun, file = "okun.rda")
# Data nedladdade 2018-12-11.



var.data <- window(sweden.data, start=c(2001, 1), end=c(2019, 1))
var.data1 <- window(sweden.data1, start=c(2001, 1), end=c(2019, 1))

newtS<-as.ts(var.data1)
describe(newtS)


par(oma=c(0,0,0,0))
plot(var.data1[,"OIL1"], 
     main="Log Oil", 
     ylab="")
plot(var.data1[,"U"], 
     main="Arbetslöshet", 
     ylab="")
plot(var.data1[,"CPI1"],
     main="Inflations takt", 
     ylab="")
plot(var.data1[,"BNP"],
     main="BNP-tillväxt", 
     ylab="")
plot(var.data1[,"INT"], 
     main="3-månaders ränta", 
     ylab="")



var.diff <- ts.intersect(diff(var.data[, "OIL"]), diff(var.data[, "U"]), var.data[,"BNP"],var.data[,"CPI"], var.data[,"INT"])
colnames(var.diff) <-c("DOIL", "DU","BNP","CPI","INT")
var.diff #OIL, Unemployment diffade

var.diff1 <- ts.intersect(diff(var.data1[, "OIL1"]), diff(var.data[, "U"]), diff(var.data[,"BNP"]), var.data1[,"CPI1"], var.data[,"INT"])
colnames(var.diff1) <-c("DOIL", "DU","DBNP","CPI","INT")
var.diff1 #OIL, Unemployment, BNP
plot(var.diff1)

############################ SEASON ADJUSTED ###############################
var.sa <- ts.intersect(diff(var.data1[, "OIL1"]), diff(var.data1[, "U"]), var.data1[,"BNP"],var.data1[,"CPI1"], var.data1[,"INT"])
colnames(var.sa) <-c("DOIL", "DU","BNP","CPI","INT")
var.sa #OIL, Unemployment diffade
plot(var.sa)
var.sa
plot(var.sa[,"DOIL"])
plot(var.sa[,"DU"])
plot(var.sa[,"CPI"])
plot(var.sa[,"BNP"])
plot(var.sa[,"INT"])
var.sa

var.sa1 <- ts.intersect(diff(var.data1[, "OIL1"]), diff(var.data1[, "U"]), var.data1[,"BNP"],var.data1[,"CPI1"], diff(var.data1[,"INT"]))
colnames(var.sa1) <-c("DOIL", "DU","BNP","CPI","DINT")
var.sa1 #OIL, Unemployment, Interest-retas diffade


ADF.test(var.data1[,"OIL1"])
ADF.test(var.data1[,"U"])
ADF.test(var.data1[,"BNP"])
ADF.test(var.data1[,"CPI1"])
ADF.test(var.data1[,"INT"])
summary(ur.kpss(var.data1[,"OIL1"]))
summary(ur.kpss(var.data1[,"U"]))
summary(ur.kpss(var.data1[,"BNP"]))
summary(ur.kpss(var.data1[,"CPI1"]))
summary(ur.kpss(var.data1[,"INT"]))
ADF.test(diff(var.data1[,"OIL1"]))
ADF.test(diff(var.data1[,"U"]))
ADF.test(diff(var.data1[,"BNP"]))
summary(ur.kpss(diff(var.data1[,"OIL1"])))
summary(ur.kpss(diff(var.data1[,"U"])))
summary(ur.kpss(diff(var.data1[,"BNP"])))
summary(ur.kpss(diff(var.data1[,"INT"])))



ADF.test(var.data1[,"BNP"]) 
ADF.test(var.data1[,"U"])
ADF.test(var.data1[,"CPI1"])
ADF.test(var.data1[,"INT"])


summary(ur.kpss(var.diff[,"INT"]))
summary(ur.kpss(var.data1[,"BNP"]))
summary(ur.kpss(var.data1[,"CPI1"]))

VARselect(var.sa, lag.max = 10, type="const")$selection
############################################################################
###################### ICKE SÄSONGSJUSTEDE VARS ############################

#VAR(1)
model.1<-VAR(var.diff, type = c("const"), p=1 ,ic = c("AIC"))
model.11<-VAR(var.diff1, type = c("const"), p=1 ,ic = c("AIC"))
AIC(model.1)
AIC(model.11)
roots(model.1)
roots(model.11)

serial.test(model.1)
serial.test(model.11)
normality.test(model.1)$jb.mul$JB
normality.test(model.11)$jb.mul$JB



#Most of the interpretation is already in the comments to the code. 
#First, a VAR(1) model is estimated. It is tested for autocorrelation in errors using a portmanteau test. 
#The null hypothesis of no autocorrelation is rejected since the 𝑝-value of 0.01996 is lower than the significance level of 0.05.

#VAR2
model.2<-VAR(var.diff, type = c("const"), p=2 ,ic = c("AIC"))
model.22<-VAR(var.diff1, type = c("const"), p=2 ,ic = c("AIC"))
AIC(model.2)
AIC(model.22)
roots(model.2)
roots(model.22)

serial.test(model.2)
serial.test(model.22)
normality.test(model.2)$jb.mul$JB
normality.test(model.22)$jb.mul$JB



############################################################
##################### SÄSONGS JUSTERADE VARS ##############
#VAR(1)
model.sa1<-VAR(var.sa, type = c("const"), p=1 ,ic = c("AIC"))
model.sa11<-VAR(var.sa1, type = c("const"), p=1 ,ic = c("AIC"))
AIC(model.sa1)
AIC(model.sa11)
roots(model.sa1)
roots(model.sa11)

serial.test(model.sa1, type = "PT.asymptotic")
serial.test(model.sa11)
normality.test(model.sa1)$jb.mul$JB
normality.test(model.sa11)$jb.mul$JB
#Most of the interpretation is already in the comments to the code. 
#First, a VAR(1) model is estimated. It is tested for autocorrelation in errors using a portmanteau test. 
#The null hypothesis of no autocorrelation is rejected since the 𝑝-value of 0.01996 is lower than the significance level of 0.05.

#VAR2
model.sa2<-VAR(var.sa, type = c("const"), p=2 ,ic = c("AIC"))
model.sa22<-VAR(var.sa1, type = c("const"), p=2 ,ic = c("AIC"))
AIC(model.sa2)
AIC(model.sa22)
roots(model.sa2)
roots(model.sa22)
serial.test(model.sa2)
serial.test(model.sa22)
normality.test(model.sa2)$jb.mul$JB
normality.test(model.sa22)$jb.mul$JB
kk
#########################################
#########################################
#########################################
##############SENASTE, HELA##############
VARselect(var.diff1, type = "const", lag.max = 10)
model.sa11<-VAR(var.diff1, type = c("const"), p=1)
model.sa1122<-VAR(var.diff1, type = c("const"), p=2)
BIC(model.sa11)
serial.test(model.sa11)
normality.test(model.sa11)$jb.mul$JB
BIC(model.sa1122)
serial.test(model.sa1122)
normality.test(model.sa1122)$jb.mul$JB
res.sa11<-residuals(model.sa11)
DOIL.min44 <- res.sa11[, "DOIL"] == min(res.sa11[, "DOIL"])
DU.max44 <- res.sa11[,"DU"] == max(res.sa11[,"DU"])
DU.min44 <- res.sa11[,"DU"] == min(res.sa11[,"DU"])
DOIL.min44 #30-th obs
DU.max44 #16-obs
DU.min44 #17-obs
var.diff1.extDD <- cbind(DOIL.min44, DU.max44, DU.min44)
var.diff1.DD <- window(var.diff1, start=c(2001, 3))
VARselect(var.diff1.DD, type = "const", exogen = var.diff1.extDD)
model.sano11<-VAR(var.diff1.DD, type= "const", p=1, exogen = var.diff1.extDD)
model.sano22<-VAR(var.diff1.DD, type= "const", p=2, exogen = var.diff1.extDD)
model.sano33<-VAR(var.diff1.DD, type= "const", p=3, exogen = var.diff1.extDD)
roots(model.sano11)
roots(model.sano22)
BIC(model.sano11)
BIC(model.sano22)
BIC(model.sano33)
normality.test(model.sano11)$jb.mul$JB
normality.test(model.sano22)$jb.mul$JB
normality.test(model.sano33)$jb.mul$JB
serial.test(model.sano11)
serial.test(model.sano22)
serial.test(model.sano33)
res.sano11<-residuals(model.sano11)

oir.bnp111 <- irf(model.sano11, impulse = "DOIL", response = "DBNP",
               n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.du111 <- irf(model.sano11, impulse = "DOIL", response = "DU",
              n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.cpi111 <- irf(model.sano11, impulse = "DOIL", response = "CPI",
               n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.int111 <- irf(model.sano11, impulse = "DOIL", response = "INT",
               n.ahead = 30, ortho = TRUE, runs = 5000, seed = 12345)
#########################################
#########################################
#########################################
plot(oir.bnp111, main="BNP Response from Impulse in OIL")
plot(oir.du111, main="U Response from Impulse in OIL")
plot(oir.cpi111, main="CPI Response from Impulse in OIL")
plot(oir.int111, main="INT Response from Impulse in OIL")

qqnorm(res.sa11[,"DOIL"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="OIL") 
qqline(res.sa11[,"DOIL"])

qqnorm(res.sa11[,"DU"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="Unemployment") 
qqline(res.sa11[,"DU"])

qqnorm(res.sa11[,"DBNP"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="BNP") 
qqline(res.sa11[,"DBNP"])

qqnorm(res.sa11[,"CPI"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="CPI") 
qqline(res.sa11[,"CPI"])

qqnorm(res.sa11[,"INT"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="INT") 
qqline(res.sa11[,"INT"])

#########################################
#########################################
#########################################
##################################################

qqnorm(res.sano11[,"DOIL"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="OIL") 
qqline(res.sano11[,"DOIL"])

qqnorm(res.sano11[,"DU"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="Unemployment") 
qqline(res.sano11[,"DU"])

qqnorm(res.sano11[,"DBNP"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="BNP") 
qqline(res.sano11[,"DBNP"])

qqnorm(res.sano11[,"CPI"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="CPI") 
qqline(res.sano11[,"CPI"])

qqnorm(res.sano11[,"INT"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="INT") 
qqline(res.sano11[,"INT"])


#############################

###### AIC SCORES ##########
AIC(model.1) #283.805
AIC(model.11) #287.3576
AIC(model.sa1) #228.3256 # 1
BIC(model.sa1)
AIC(model.sa11) #230.3381 # 2
AIC(model.2) #277.0217
AIC(model.22) #277.8302
AIC(model.sa2) #240.7151 # 3
BIC(model.sa2)
AIC(model.sa22) #241.9763 #4

############################
###### NORMALITY ##########

res.sa1<-residuals(model.sa1)
res.sa2<-residuals(model.sano2)
res.sa1
DOIL.min4 <- res.sa1[, "DOIL"] == min(res.sa1[, "DOIL"])
DOIL.min4 #30-th obs
DU.max4 #16-obs
DU.max4 <- res.sa1[,"DU"] == max(res.sa1[,"DU"])
#BNP.min4 <- res.sa1[, "BNP"] == min(res.sa1[, "BNP"])
#INT.min4 <- res.sa1[, "INT"] == min(res.sa1[, "INT"])
var.sa.extDD <- cbind(DOIL.min4, DU.max4)
var.sa.DD <- window(var.sa, start=c(2001, 3))
VARselect(var.sa.DD, type = "const", exogen = var.sa.extDD)


VARselect(var.sa.DD, lag.max = 10, type = c("const"),
          exogen = var.sa.extDD)
model.sano1<-VAR(var.sa.DD, type= "const", p=1, exogen = var.sa.extDD)
model.sano2<-VAR(var.sa.DD, type= "const", p=2, exogen = var.sa.extDD)
model.sano3<-VAR(var.sa.DD, type= "const", p=3, exogen = var.sa.extDD)
model.sano4<-VAR(var.sa.DD, type= "const", p=4, exogen = var.sa.extDD)
model.sano5<-VAR(var.sa.DD, type= "const", p=5, exogen = var.sa.extDD)
AIC(model.sano3)

model.sano1
model.sano2


AIC(model.sano1)
AIC(model.sano2)
AIC(model.sano3)
AIC(model.sano4)
AIC(model.sano5)

BIC(model.sano1)
BIC(model.sano2)
BIC(model.sano3)
BIC(model.sano4)
BIC(model.sano5)


roots(model.sano1)
roots(model.sano2)

serial.test(model.sano1)
serial.test(model.sano2)
serial.test(model.sano3)
serial.test(model.sano4)
serial.test(model.sano5)

normality.test(model.sano1)$jb.mul$JB
normality.test(model.sano2)$jb.mul$JB
normality.test(model.sano3)$jb.mul$JB
normality.test(model.sano4)$jb.mul$JB
normality.test(model.sano5)$jb.mul$JB

qqnorm(res.sa1[,"DOIL"], 
             ylab="Standardized Residuals", 
             xlab="Normal Scores", 
            main="OIL") 
qqline(res.sa1[,"DOIL"])

qqnorm(res.sa1[,"DU"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="Unemployment") 
qqline(res.sa2[,"DU"])

qqnorm(res.sa1[,"BNP"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="BNP") 
qqline(res.sa2[,"BNP"])

qqnorm(res.sa1[,"CPI"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="CPI") 
qqline(res.sa1[,"CPI"])

qqnorm(res.sa1[,"INT"], 
       ylab="Standardized Residuals", 
       xlab="Normal Scores", 
       main="INT") 
qqline(res.sa2[,"INT"])

############################
############################################################
############### IRF-FUNCTIONS ##############################

oir.bnp <- irf(model.sano2, impulse = "DOIL", response = "BNP",
           n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.du <- irf(model.sano2, impulse = "DOIL", response = "DU",
           n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.cpi <- irf(model.sano2, impulse = "DOIL", response = "CPI",
           n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.int <- irf(model.sano2, impulse = "DOIL", response = "INT",
           n.ahead = 30, ortho = TRUE, runs = 5000, seed = 12345)

plot(oir.bnp, main="BNP Response from Impulse in OIL")
plot(oir.du, main="U Response from Impulse in OIL")
plot(oir.cpi, main="CPI Response from Impulse in OIL")
plot(oir.int, main="INT Response from Impulse in OIL")

###############################################################################

all <- irf(model.sa2, impulse = "DOIL", response = c("BNP", "DU", "CPI", "INT"), 
n.ahead = 15, ortho = TRUE, runs = 5000, seed = 12345)

all1 <- irf(model.22, impulse = "DOIL", response = c("BNP", "DU", "CPI", "DINT"), 
           n.ahead = 15, ortho = TRUE, runs = 500, seed = 12345)

plot(all,main="Impulse Response functions")
plot(all1)

############ EXTRA####################

# Calculate summary statistics
model_summary <- summary(model.sa2)
model_summary1 <- summary(model.aic1)
model_summary2 <- summary(model.aic2)

# Obtain variance-covariance matrix
model_summary$covres
chol <- t(chol(model_summary$covres))
print(xtable(chol), type="html")
write.table(chol, file = "olstab.txt", sep = ",", quote = FALSE, row.names = F)

t(chol(model_summary1$covres))
t(chol(model_summary2$covres))

oir <- irf(model.aic, impulse = "oil", response = "un",
           n.ahead = 12, ortho = TRUE, runs = 5000, seed = 12345)
plot(oir)

roots(model.aic2)





###############################################################
### EXPORT GRAPHICS ###
pdf("plots_VAR(1)_FINAL.pdf")
par(mfrow=c(2,2))
plot(oir.bnp111, main="DBNP Response from Impulse in DOIL")
plot(oir.du111, main="DU Response from Impulse in DOIL")
plot(oir.cpi111, main="CPI Response from Impulse in DOIL")
plot(oir.int111, main="INT Response from Impulse in DOIL")
par(mfrow=c(1,1))
graphics.off()

######################
### EXPORT TABLES ###

sink("IRF.txt")
oir.bnp111
oir.du111
oir.cpi111
oir.int111
sink()



pdf("data_plots.pdf")
par(mfrow=c(3,2))
plot(var.data1[,"OIL1"], 
     main="Log Oil", 
     ylab="",
     xlab="År")
plot(var.data1[,"U"], 
     main="Arbetslöshet", 
     ylab="",
     xlab="År")
plot(var.data1[,"CPI1"],
     main="Inflations tkt", 
     ylab="",
     xlab="År")
plot(var.data1[,"BNP"],
     main="BNP-tillväxt", 
     ylab="",
     xlab="År")
plot(var.data1[,"INT"], 
     main="3-månaders rnta", 
     ylab="",
     xlab="År")
par(mfrow=c(1,1))
graphics.off()

##############################
oir.bnp4 <- irf(model.sano3, impulse = "DOIL", response = "BNP",
               n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.du4 <- irf(model.sano3, impulse = "DOIL", response = "DU",
              n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.cpi4 <- irf(model.sano3, impulse = "DOIL", response = "CPI",
               n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.int4 <- irf(model.sano3, impulse = "DOIL", response = "INT",
               n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)

plot(oir.bnp4, main="BNP Response from Impulse in OIL")
plot(oir.du4, main="U Response from Impulse in OIL")
plot(oir.cpi4, main="CPI Response from Impulse in OIL")
plot(oir.int4, main="INT Response from Impulse in OIL")

oir.bnp1 <- irf(model.sano1, impulse = "DOIL", response = "BNP",
                n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.du1 <- irf(model.sano1, impulse = "DOIL", response = "DU",
               n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.cpi1 <- irf(model.sano1, impulse = "DOIL", response = "CPI",
                n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)
oir.int1 <- irf(model.sano1, impulse = "DOIL", response = "INT",
                n.ahead = 30, ortho = TRUE, runs = 500, seed = 12345)

plot(oir.bnp1, main="BNP Response from Impulse in OIL")
plot(oir.du1, main="U Response from Impulse in OIL")
plot(oir.cpi1, main="CPI Response from Impulse in OIL")
plot(oir.int1, main="INT Response from Impulse in OIL")

