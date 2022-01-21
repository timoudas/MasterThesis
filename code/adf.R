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