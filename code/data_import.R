library(quantmod)
library(tseries)
library(seasonal)

#Constant Price Gross Domestic Product in Sweden, Seasonally Adjusted
#Percent Change from Year Ago, Quarterly, Seasonally Adjusted
getSymbols("SWEGDPRQPSMEI", src = "FRED")
BNP <- ts(as.ts(SWEGDPRQPSMEI), start = c(1961, 1), frequency = 4)


#Consumer Price Index: Total All Items for Sweden Growth Rate Previous Period,
#Not Seasonally Adjusted
getSymbols("CPALTT01SEQ657N", src = "FRED")
CPI <- ts(as.ts(CPALTT01SEQ657N), start = c(1960, 1), frequency = 4)
# seasonal adjust time series using X11. 
CPI1 <- final(seas(as.ts((CPI),freq=4)))

#3-Month or 90-day Rates and Yields: Treasury Securities for Sweden
#Percent,Not Seasonally Adjusted
getSymbols("IR3TTS01SEQ156N", src = "FRED")
INT <- ts(as.ts(IR3TTS01SEQ156N), start = c(1982, 1), frequency = 4)
# seasonal adjust time series using X11. 
INT1 <- final(seas(as.ts((INT),freq=4)))

#CLIF INDEX
CLIF_INDEX <- read.csv('../ext_data/CLIFS.csv', header = FALSE)
CLIF_INDEX <-CLIF_INDEX[nrow(CLIF_INDEX):1,]
CLIF_mts <- ts(CLIF_INDEX[,2], start = c(1970, 1), frequency = 12)
CLIF_qts = aggregate(CLIF_mts, nfrequency = 4)


#Data in levels
data <- ts.intersect(CLIF_qts, BNP, CPI1,  INT1)
colnames(data) <-c("CLIF", "BNP", 'CPI', "INT")
#Start sample at 1993 to not have to account for paradigm shift
#Sweden let currency float autumn 1992
sample <- window(data, start=c(1993,1))

