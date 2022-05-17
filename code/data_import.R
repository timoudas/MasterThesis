library(quantmod)
library(tseries)
library(seasonal)
library("readxl")
library(ggplot2)
library(xts)


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


DATA_DSFSI <- ts.intersect(diff(SFSI_INDEX), diff(BNP), CPI, INT)
colnames(DATA_DSFSI) <-c("DSFSI", "DBNP", 'CPI', "INT")
write.csv(DATA_DSFSI, file="data2.csv")


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


SFSI_INDEX_p <- read_excel("../ext_data/stressindex.xlsx", sheet = "SVE")
SFSI_INDEX_p$Datum <- as.Date(SFSI_INDEX_p$Datum) 



ggplot(data=SFSI_INDEX_p, aes(x=Datum, y=Stressindex)) +
    geom_line(color="darkgreen", size=0.5) +
    ylab("FSI level") +
    xlab('Date')


BANK_INDEX <- read.csv2("bank_index.csv", sep=';', header = TRUE)

BANK_INDEX <- xts(BANK_INDEX, order.by=as.Date(BANK_INDEX$Datum, "%Y-%m-%d"))
BANK_INDEX$Datum <- NULL
storage.mode(BANK_INDEX) <- "numeric"
BANK_INDEX <- apply.quarterly(BANK_INDEX, mean) #SUM TO QUATERLY
BANK_INDEX <- ts(BANK_INDEX$Stängning, start = c(2000, 3, 31), frequency = 4)
BANK_INDEX <- final(seas(as.ts((BANK_INDEX),freq=4)))
SD <- rollapply(data = diff(BANK_INDEX), width=4,FUN=sd) 

Sd_data <- data.frame(Y=as.matrix(SD), date=as.Date(SD))
write.csv(Sd_data, file = 'bank_vol.csv')
