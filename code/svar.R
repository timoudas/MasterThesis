source('stationarity_tests.R')
library(vars)


VARselect(model_data.SFI, lag.max = 10, type='const')$selection

model <- VAR(model_data.SFI, type='const', lag=1)

FEIR <- irf(model, impulse = "SFSI", response = "BNP", n.ahead = 20, ci = 0.90, type = "oir", runs=1000)
par(mfrow=c(1,1))
plot(FEIR,  
     main = "",
     col = 2, 
)

