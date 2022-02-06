library(quantmod)
library(tseries)
library(seasonal)
library("readxl")


#Constant Price Gross Domestic Product in Sweden, Seasonally Adjusted
#Percent Change from Year Ago, Quarterly, Seasonally Adjusted
getSymbols("SWEGDPRQPSMEI", src = "FRED")
BNP <- ts(as.ts(SWEGDPRQPSMEI), start = c(1961, 1), frequency = 4)
plot(BNP)

#Consumer Price Index: Total All Items for Sweden Growth Rate Previous Period,
#Not Seasonally Adjusted
getSymbols("CPALTT01SEQ657N", src = "FRED")
CPI <- ts(as.ts(CPALTT01SEQ657N), start = c(1960, 1), frequency = 4)
# seasonal adjust time series using X11. 
CPI <- final(seas(as.ts((CPI),freq=4)))

#3-Month or 90-day Rates and Yields: Treasury Securities for Sweden
#Percent,Not Seasonally Adjusted
getSymbols("IR3TTS01SEQ156N", src = "FRED")
INT <- ts(as.ts(IR3TTS01SEQ156N), start = c(1982, 1), frequency = 4)
# seasonal adjust time series using X11. 
INT <- final(seas(as.ts((INT),freq=4)))

#CLIF INDEX
CLIF_INDEX <- read.csv('../ext_data/CLIFS.csv', header = FALSE)
CLIF_INDEX <-CLIF_INDEX[nrow(CLIF_INDEX):1,]
CLIF_INDEX <- ts(CLIF_INDEX[,2], start = c(1970, 1), frequency = 12)
CLIF_INDEX = aggregate(CLIF_INDEX, nfrequency = 4) #SUM TO QUATERLY


SFSI_INDEX <- read_excel("../ext_data/stressindex.xlsx", sheet = "SVE")
SFSI_INDEX <- xts(SFSI_INDEX[,-1], order.by=as.Date(SFSI_INDEX$Datum, "%Y-%m-%d"))
SFSI_INDEX <- apply.quarterly(SFSI_INDEX, sum) #SUM TO QUATERLY
SFSI_INDEX <- ts(SFSI_INDEX$Stressindex, start = c(1995, 3, 31), frequency = 4)

ECB_INDEX <- read_excel("../ext_data/stressindex.xlsx", sheet = "ECB")
ECB_INDEX <- xts(ECB_INDEX[,-1], order.by=as.Date(ECB_INDEX$Datum, "%Y-%m-%d"))
ECB_INDEX <- apply.quarterly(ECB_INDEX, sum) #SUM TO QUATERLY
ECB_INDEX <- ts(ECB_INDEX$Stressindex, start = c(1999, 3, 26), frequency = 4)


#DATA IN LEVELS WITH RIKSBANK FIN STRESS INDEX
DATA_SFSI <- ts.intersect(SFSI_INDEX, BNP, CPI,  INT)
colnames(DATA_SFSI) <-c("SFSI", "BNP", 'CPI', "INT")

#DATA IN LEVELS WITH RIKSBANK FIN STRESS INDEX
DATA_ECB <- ts.intersect(ECB_INDEX, BNP, CPI,  INT)
colnames(DATA_ECB) <-c("ECBFSI", "BNP", 'CPI', "INT")

#Data IN LEVELS WITH CLIFS INDEX
DATA_CLIF <- ts.intersect(CLIF_INDEX, BNP, CPI,  INT)
colnames(DATA_CLIF) <-c("CLIF", "BNP", 'CPI', "INT")
#Start sample at 1993 to not have to account for paradigm shift
#Sweden let currency float autumn 1992
DATA_CLIF <- window(DATA_CLIF, start=c(1993,1))

