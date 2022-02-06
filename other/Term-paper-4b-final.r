  ###############################################################################
  ################################################################################
  ######                                                       ###################
  ######                  TERMPAPER - AUTUMN 2021              ###################
  ######                                                       ###################
  ###### By: Petter Söderlind                                  ###################
  ###### Professor: Annika Alexius                             ###################
  ###### Interest rates and house prices                       ###################
  ######                                                       ###################
  ######                                                       ###################
  ################################################################################
  ################################################################################
  ################# Clearing the environment #####################################
  rm(list = ls())
  graphics.off()
  ################ Installing necessary packages #################################
  # install.packages("pxweb")
  # install.packages("readxl")
  # install.packages("bvartools", dependencies = TRUE)
  # install.packages("pracma")
  # install.packages("lmtest")
  # install.packages("LaplacesDemon", dependencies = TRUE)
  # install.packages("minnesota_prior{bvartools}", force = TRUE)
  # install.packages("forecast")
  # install.packages("tseries")
  ############### loading packages ###############################################
  library(bvartools)
  library(readxl)
  library(pxweb)
  library(pracma)
  library(lmtest)
  library(urca)
  library(forecast)
  library(tseries)
  ############## Setting working directory #######################################
  # setwd("/Users/pettersoderlind/Documents/R/BayesianTEST")
  ################################################################################
  #############                                  #################################
  #############            DOWNLOADING DATA      #################################
  #############                                  #################################
  ################################################################################
  CPI_query_list <- 
    list("ContentsCode"=c("000004VU"),
         "Tid"=c("1990M01","1990M02","1990M03","1990M04","1990M05","1990M06","1990M07","1990M08","1990M09","1990M10","1990M11","1990M12","1991M01","1991M02","1991M03","1991M04","1991M05","1991M06","1991M07","1991M08","1991M09","1991M10","1991M11","1991M12","1992M01","1992M02","1992M03","1992M04","1992M05","1992M06","1992M07","1992M08","1992M09","1992M10","1992M11","1992M12","1993M01","1993M02","1993M03","1993M04","1993M05","1993M06","1993M07","1993M08","1993M09","1993M10","1993M11","1993M12","1994M01","1994M02","1994M03","1994M04","1994M05","1994M06","1994M07","1994M08","1994M09","1994M10","1994M11","1994M12","1995M01","1995M02","1995M03","1995M04","1995M05","1995M06","1995M07","1995M08","1995M09","1995M10","1995M11","1995M12","1996M01","1996M02","1996M03","1996M04","1996M05","1996M06","1996M07","1996M08","1996M09","1996M10","1996M11","1996M12","1997M01","1997M02","1997M03","1997M04","1997M05","1997M06","1997M07","1997M08","1997M09","1997M10","1997M11","1997M12","1998M01","1998M02","1998M03","1998M04","1998M05","1998M06","1998M07","1998M08","1998M09","1998M10","1998M11","1998M12","1999M01","1999M02","1999M03","1999M04","1999M05","1999M06","1999M07","1999M08","1999M09","1999M10","1999M11","1999M12","2000M01","2000M02","2000M03","2000M04","2000M05","2000M06","2000M07","2000M08","2000M09","2000M10","2000M11","2000M12","2001M01","2001M02","2001M03","2001M04","2001M05","2001M06","2001M07","2001M08","2001M09","2001M10","2001M11","2001M12","2002M01","2002M02","2002M03","2002M04","2002M05","2002M06","2002M07","2002M08","2002M09","2002M10","2002M11","2002M12","2003M01","2003M02","2003M03","2003M04","2003M05","2003M06","2003M07","2003M08","2003M09","2003M10","2003M11","2003M12","2004M01","2004M02","2004M03","2004M04","2004M05","2004M06","2004M07","2004M08","2004M09","2004M10","2004M11","2004M12","2005M01","2005M02","2005M03","2005M04","2005M05","2005M06","2005M07","2005M08","2005M09","2005M10","2005M11","2005M12","2006M01","2006M02","2006M03","2006M04","2006M05","2006M06","2006M07","2006M08","2006M09","2006M10","2006M11","2006M12","2007M01","2007M02","2007M03","2007M04","2007M05","2007M06","2007M07","2007M08","2007M09","2007M10","2007M11","2007M12","2008M01","2008M02","2008M03","2008M04","2008M05","2008M06","2008M07","2008M08","2008M09","2008M10","2008M11","2008M12","2009M01","2009M02","2009M03","2009M04","2009M05","2009M06","2009M07","2009M08","2009M09","2009M10","2009M11","2009M12","2010M01","2010M02","2010M03","2010M04","2010M05","2010M06","2010M07","2010M08","2010M09","2010M10","2010M11","2010M12","2011M01","2011M02","2011M03","2011M04","2011M05","2011M06","2011M07","2011M08","2011M09","2011M10","2011M11","2011M12","2012M01","2012M02","2012M03","2012M04","2012M05","2012M06","2012M07","2012M08","2012M09","2012M10","2012M11","2012M12","2013M01","2013M02","2013M03","2013M04","2013M05","2013M06","2013M07","2013M08","2013M09","2013M10","2013M11","2013M12","2014M01","2014M02","2014M03","2014M04","2014M05","2014M06","2014M07","2014M08","2014M09","2014M10","2014M11","2014M12","2015M01","2015M02","2015M03","2015M04","2015M05","2015M06","2015M07","2015M08","2015M09","2015M10","2015M11","2015M12","2016M01","2016M02","2016M03","2016M04","2016M05","2016M06","2016M07","2016M08","2016M09","2016M10","2016M11","2016M12","2017M01","2017M02","2017M03","2017M04","2017M05","2017M06","2017M07","2017M08","2017M09","2017M10","2017M11","2017M12","2018M01","2018M02","2018M03","2018M04","2018M05","2018M06","2018M07","2018M08","2018M09","2018M10","2018M11","2018M12","2019M01","2019M02","2019M03","2019M04","2019M05","2019M06","2019M07","2019M08","2019M09","2019M10","2019M11","2019M12","2020M01","2020M02","2020M03","2020M04","2020M05","2020M06","2020M07","2020M08","2020M09","2020M10","2020M11","2020M12","2021M01","2021M02","2021M03","2021M04","2021M05","2021M06","2021M07","2021M08","2021M09"))
  CPI_data <- 
    pxweb_get(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/PR/PR0101/PR0101A/KPItotM",
              query = CPI_query_list) # Download data 
  ############## Convert to data.frame  ##########################################
  CPI_data_frame <- as.data.frame(CPI_data, column.name.type = "text", variable.value.type = "text")
  CPImonthly <- ts(CPI_data_frame[, -c(1)], start = c(1990,1), frequency = 12) ## Monthly data
  CPI <- aggregate(CPImonthly, start = c(1990, 1), nfrequency = 4) ## Quarterly data
  ################################################################################
  GDP_query_list <- 
    list("Anvandningstyp"=c("BNPM"),
         "ContentsCode"=c("NR0103BV"),
         "Tid"=c("1990K1","1990K2","1990K3","1990K4","1991K1","1991K2","1991K3","1991K4","1992K1","1992K2","1992K3","1992K4","1993K1","1993K2","1993K3","1993K4","1994K1","1994K2","1994K3","1994K4","1995K1","1995K2","1995K3","1995K4","1996K1","1996K2","1996K3","1996K4","1997K1","1997K2","1997K3","1997K4","1998K1","1998K2","1998K3","1998K4","1999K1","1999K2","1999K3","1999K4","2000K1","2000K2","2000K3","2000K4","2001K1","2001K2","2001K3","2001K4","2002K1","2002K2","2002K3","2002K4","2003K1","2003K2","2003K3","2003K4","2004K1","2004K2","2004K3","2004K4","2005K1","2005K2","2005K3","2005K4","2006K1","2006K2","2006K3","2006K4","2007K1","2007K2","2007K3","2007K4","2008K1","2008K2","2008K3","2008K4","2009K1","2009K2","2009K3","2009K4","2010K1","2010K2","2010K3","2010K4","2011K1","2011K2","2011K3","2011K4","2012K1","2012K2","2012K3","2012K4","2013K1","2013K2","2013K3","2013K4","2014K1","2014K2","2014K3","2014K4","2015K1","2015K2","2015K3","2015K4","2016K1","2016K2","2016K3","2016K4","2017K1","2017K2","2017K3","2017K4","2018K1","2018K2","2018K3","2018K4","2019K1","2019K2","2019K3","2019K4","2020K1","2020K2","2020K3","2020K4","2021K1","2021K2","2021K3"))
  
  GDP_data <- 
    pxweb_get(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/NR/NR0103/NR0103A/NR0103ENS2010T01Kv",
              query = GDP_query_list) # Download data 
  
  ############## Convert to data.frame  ##########################################
  GDP_data_frame <- as.data.frame(GDP_data, column.name.type = "text", variable.value.type = "text")
  GDP <- ts(GDP_data_frame[, -c(1:2)], start = c(1990,1), frequency = 4) ## Quarterly
  ################################################################################
  HP_query_list <- 
    list("ContentsCode"=c("BO0501K6"),
         "Tid"=c("1990K1","1990K2","1990K3","1990K4","1991K1","1991K2","1991K3","1991K4","1992K1","1992K2","1992K3","1992K4","1993K1","1993K2","1993K3","1993K4","1994K1","1994K2","1994K3","1994K4","1995K1","1995K2","1995K3","1995K4","1996K1","1996K2","1996K3","1996K4","1997K1","1997K2","1997K3","1997K4","1998K1","1998K2","1998K3","1998K4","1999K1","1999K2","1999K3","1999K4","2000K1","2000K2","2000K3","2000K4","2001K1","2001K2","2001K3","2001K4","2002K1","2002K2","2002K3","2002K4","2003K1","2003K2","2003K3","2003K4","2004K1","2004K2","2004K3","2004K4","2005K1","2005K2","2005K3","2005K4","2006K1","2006K2","2006K3","2006K4","2007K1","2007K2","2007K3","2007K4","2008K1","2008K2","2008K3","2008K4","2009K1","2009K2","2009K3","2009K4","2010K1","2010K2","2010K3","2010K4","2011K1","2011K2","2011K3","2011K4","2012K1","2012K2","2012K3","2012K4","2013K1","2013K2","2013K3","2013K4","2014K1","2014K2","2014K3","2014K4","2015K1","2015K2","2015K3","2015K4","2016K1","2016K2","2016K3","2016K4","2017K1","2017K2","2017K3","2017K4","2018K1","2018K2","2018K3","2018K4","2019K1","2019K2","2019K3","2019K4","2020K1","2020K2","2020K3","2020K4","2021K1","2021K2","2021K3"))
  
  HP_data <- 
    pxweb_get(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/BO/BO0501/BO0501A/FastpiFritidshusKv",
              query = HP_query_list)# Download data 
  HP_data_frame <- as.data.frame(HP_data, column.name.type = "text", variable.value.type = "text")
  HP <- ts(HP_data_frame[, -c(1)], start = c(1990,1), frequency = 4) ## Quarterly
  ################################################################################
  # PXWEB query 
  PI_query_list <- 
    list("ContentsCode"=c("PR0101C5"),
         "Tid"=c("1990M01","1990M02","1990M03","1990M04","1990M05","1990M06","1990M07","1990M08","1990M09","1990M10","1990M11","1990M12","1991M01","1991M02","1991M03","1991M04","1991M05","1991M06","1991M07","1991M08","1991M09","1991M10","1991M11","1991M12","1992M01","1992M02","1992M03","1992M04","1992M05","1992M06","1992M07","1992M08","1992M09","1992M10","1992M11","1992M12","1993M01","1993M02","1993M03","1993M04","1993M05","1993M06","1993M07","1993M08","1993M09","1993M10","1993M11","1993M12","1994M01","1994M02","1994M03","1994M04","1994M05","1994M06","1994M07","1994M08","1994M09","1994M10","1994M11","1994M12","1995M01","1995M02","1995M03","1995M04","1995M05","1995M06","1995M07","1995M08","1995M09","1995M10","1995M11","1995M12","1996M01","1996M02","1996M03","1996M04","1996M05","1996M06","1996M07","1996M08","1996M09","1996M10","1996M11","1996M12","1997M01","1997M02","1997M03","1997M04","1997M05","1997M06","1997M07","1997M08","1997M09","1997M10","1997M11","1997M12","1998M01","1998M02","1998M03","1998M04","1998M05","1998M06","1998M07","1998M08","1998M09","1998M10","1998M11","1998M12","1999M01","1999M02","1999M03","1999M04","1999M05","1999M06","1999M07","1999M08","1999M09","1999M10","1999M11","1999M12","2000M01","2000M02","2000M03","2000M04","2000M05","2000M06","2000M07","2000M08","2000M09","2000M10","2000M11","2000M12","2001M01","2001M02","2001M03","2001M04","2001M05","2001M06","2001M07","2001M08","2001M09","2001M10","2001M11","2001M12","2002M01","2002M02","2002M03","2002M04","2002M05","2002M06","2002M07","2002M08","2002M09","2002M10","2002M11","2002M12","2003M01","2003M02","2003M03","2003M04","2003M05","2003M06","2003M07","2003M08","2003M09","2003M10","2003M11","2003M12","2004M01","2004M02","2004M03","2004M04","2004M05","2004M06","2004M07","2004M08","2004M09","2004M10","2004M11","2004M12","2005M01","2005M02","2005M03","2005M04","2005M05","2005M06","2005M07","2005M08","2005M09","2005M10","2005M11","2005M12","2006M01","2006M02","2006M03","2006M04","2006M05","2006M06","2006M07","2006M08","2006M09","2006M10","2006M11","2006M12","2007M01","2007M02","2007M03","2007M04","2007M05","2007M06","2007M07","2007M08","2007M09","2007M10","2007M11","2007M12","2008M01","2008M02","2008M03","2008M04","2008M05","2008M06","2008M07","2008M08","2008M09","2008M10","2008M11","2008M12","2009M01","2009M02","2009M03","2009M04","2009M05","2009M06","2009M07","2009M08","2009M09","2009M10","2009M11","2009M12","2010M01","2010M02","2010M03","2010M04","2010M05","2010M06","2010M07","2010M08","2010M09","2010M10","2010M11","2010M12","2011M01","2011M02","2011M03","2011M04","2011M05","2011M06","2011M07","2011M08","2011M09","2011M10","2011M11","2011M12","2012M01","2012M02","2012M03","2012M04","2012M05","2012M06","2012M07","2012M08","2012M09","2012M10","2012M11","2012M12","2013M01","2013M02","2013M03","2013M04","2013M05","2013M06","2013M07","2013M08","2013M09","2013M10","2013M11","2013M12","2014M01","2014M02","2014M03","2014M04","2014M05","2014M06","2014M07","2014M08","2014M09","2014M10","2014M11","2014M12","2015M01","2015M02","2015M03","2015M04","2015M05","2015M06","2015M07","2015M08","2015M09","2015M10","2015M11","2015M12","2016M01","2016M02","2016M03","2016M04","2016M05","2016M06","2016M07","2016M08","2016M09","2016M10","2016M11","2016M12","2017M01","2017M02","2017M03","2017M04","2017M05","2017M06","2017M07","2017M08","2017M09","2017M10","2017M11","2017M12","2018M01","2018M02","2018M03","2018M04","2018M05","2018M06","2018M07","2018M08","2018M09","2018M10","2018M11","2018M12","2019M01","2019M02","2019M03","2019M04","2019M05","2019M06","2019M07","2019M08","2019M09","2019M10","2019M11","2019M12","2020M01","2020M02","2020M03","2020M04","2020M05","2020M06","2020M07","2020M08","2020M09","2020M10","2020M11","2020M12","2021M01","2021M02","2021M03","2021M04","2021M05","2021M06","2021M07","2021M08","2021M09"))
  
  # Download data 
  PI_data <- 
    pxweb_get(url = "http://api.scb.se/OV0104/v1/doris/en/ssd/PR/PR0101/PR0101A/KPIFFMP",
              query = PI_query_list)
  
  # Convert to data.frame 
  PI_data_frame <- as.data.frame(PI_data, column.name.type = "text", variable.value.type = "text")
  ############## Convert to data.frame  ##########################################
  monthly <- ts(PI_data_frame[ -c(1) ], start = 1990, frequency = 12) #quarterly inflation
  quarterly <- aggregate(monthly, nfrequency = 4)
  PI <- ts(quarterly, start = c(1990), frequency = 4) #quarterly inflation(remeber to check if this is correct!)
  ################################################################################
  r_data_frame <- read_excel("ADVANCED RESULT_2021-12-14_13_49.xlsx" ) # Get interest rate data
  R1 <- ts(r_data_frame[ -c(1:4, 112:118), ]$`...2`, start = c(1995,1), frequency = 4) # Quarterly
  R <- ts(as.numeric(R1), start = c(1995,1), frequency = 4) #quarterly interest rates
  ################################################################################
  #############                                  #################################
  #############     FINALIZING DATA-SET          #################################
  #############                                  #################################
  ################################################################################
  realHP <- diff(log((HP)/(CPI)))
  realHP <- realHP * 100
  realGDP <- diff(log((GDP)/(CPI)))
  realGDP <- realGDP * 100
  realR <- R-PI ## Getting the real interest rates, taylor rule = R* = R - PI 
  ##We are interested in the real interest rate NOT the first difference.There is a trend here
  # Is it deterministic or not? De-trending 
  myGDP <- ts(realGDP[-c(1:20)], start = 1995, frequency = 4)#Creating the subset of interest
  myHP <- ts(realHP[-c(1:20)], start = 1995, frequency = 4)
  myR <- ts(realR[-c(1:20)], start = 1995, frequency = 4)
  trend=lm(realR~c(1:length(realR))) 
  myR =ts(residuals(trend),start = c(1995,1), frequency = 4)# Detrended real interest rates.
  priorGDP <- ts(realGDP[(1:20)], start = 1990, frequency = 4)#Creating data for constructing prior values.
  priorHP <- ts(realHP[(1:20)], start = 1990, frequency = 4)
  priorR <- ts(realR[(1:20)], start = 1990, frequency = 4)
  ################################################################################
  #############                                  #################################
  #############     FINAL DATA-SET               #################################
  #############                                  #################################
  ################################################################################
  fulldata <- ts.intersect(realGDP, realR, realHP) # One data set with all variables
  mydata <- ts.intersect(myGDP, myR,myHP) # Subset of interest
  priordata <- ts.intersect(priorGDP, priorR,priorHP) # prior-data
  rawdata <- ts.intersect(GDP,R,HP,CPI, PI) #raw data
  
  plot(rawdata, main = "Raw data")
  plot(priordata, main = "prior data")
  plot(mydata, main = "the data") #plotting data 1995Q1-2021Q3
  ################################################################################
  #############                                  #################################
  #############     STATIONARITY TEST            #################################
  #############                                  #################################
  ################################################################################
  # We now look at whether our variables are stationary. From first glance we
  # Observe that these variables are not stationary in the raw data, using the stationarity
  # rule of thumb the trimmed and transformed data-set, "mydata" looks stationary.
  # using "ndiffs" to confirm stationary, no further differenciation should be needed for 
  # stationarity in the data.
  ndiffs(myHP)
  ndiffs(myGDP)
  ndiffs(myR)
  ################################################################################
  #############                                  #################################
  #############  FREQUENTIST ESTIMATION          #################################
  #############  (FOR PRIORS/BENCHMARK)          #################################
  #############                                  #################################
  ################################################################################
  # To get a benchmark, the data used is for the years 1995-2000. 
  data <- gen_var(priordata, p = 3, deterministic = "const", 
                  iterations = 10000, burnin = 5000)
  y <- t(data$data$Y)
  x <- t(data$data$Z)
  A_freq <- tcrossprod(y, x) %*% solve(tcrossprod(x)) # Calculate estimates
  round(A_freq, 3) # Round estimates and print
  u_freq <- y - A_freq %*% x
  u_sigma_freq <- tcrossprod(u_freq) / (ncol(y) - nrow(x))
  u_sigma_freq <- round(u_sigma_freq, 2)
  ################################################################################
  #############                                  #################################
  #############     BVAR ANALYSIS                #################################
  #############                                  #################################
  ################################################################################
  # The gen_var function produces the inputs y and x for the estimator, 
  # where y is a matrix of dependent variables and x is the matrix of regressors for the model
  data <- gen_var(mydata, p = 3, deterministic = "const", 
                  iterations = 30000, burnin = 10000)
  y <- t(data$data$Y)
  x <- t(data$data$Z)
  ############## Bayesian estimation ######################################
  # Reset random number generator for reproducibility
  set.seed(1234567)
  iter <- 30000 # Number of iterations of the Gibbs sampler
  burnin <- 10000 # Number of burn-in draws
  store <- iter - burnin
  tt <- ncol(y) # Number of observations
  k <- nrow(y) # Number of endogenous variables
  m <- k * nrow(x) # Number of estimated coefficients
  #Set priors
  a_mu_prior <- matrix(0, m) # Vector of prior parameter means
  a_v_i_prior <- diag(1, m) # Inverse of the prior covariance matrix(ie the precision)
  u_sigma_df_prior <- 6 # Prior degrees of freedom
  u_sigma_scale_prior <- diag(1, k) # Prior covariance matrix
  u_sigma_df_post <- tt + u_sigma_df_prior # Posterior degrees of freedom
  # Initial values
  u_sigma_i <- solve(u_sigma_freq)
  # Data containers for posterior draws
  draws_a <- matrix(NA, m, store)
  draws_sigma <- matrix(NA, k * k, store)
  # Start Gibbs sampler
  for (draw in 1:iter) {
    
    a <- post_normal(y, x, u_sigma_i, a_mu_prior, a_v_i_prior) # Draw conditional mean parameters
    
    # Draw variance-covariance matrix'
    
    u <- y - matrix(a, k) %*% x # Obtain residuals
    u_sigma_scale_post <- solve(u_sigma_scale_prior + tcrossprod(u))
    u_sigma_i <- matrix(rWishart(1, u_sigma_df_post, u_sigma_scale_post)[,, 1], k)
    u_sigma <- solve(u_sigma_i) # Invert Sigma_i to obtain Sigma
    
    
    
    if (draw > burnin) {
      draws_a[, draw - burnin] <- a
      draws_sigma[, draw - burnin] <- u_sigma
    }
  }# Store draws
  
  
  ## After we run the Gibbs sampler we obtain point estimates for the coefficient matrix = the mean of the posterior draws. 
  
  A <- rowMeans(draws_a) # Obtain means for every row
  A <- matrix(A, k) # Transform mean vector into a matrix
  A <- round(A, 3) # Round values
  dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]]) # Rename matrix dimensions
  
  A # Print
  
  Sigma <- rowMeans(draws_sigma) # Obtain means for every row
  Sigma <- matrix(Sigma, k) # Transform mean vector into a matrix
  Sigma <- round(Sigma, 2) # Round values
  dimnames(Sigma) <- list(dimnames(y)[[1]], dimnames(y)[[1]]) # Rename matrix dimensions
  
  Sigma # Print
  
  bvar_est <- bvar(y = data$data$Y, x = data$data$Z, A = draws_a[1:18,],
                   C = draws_a[19:21, ], Sigma = draws_sigma)
  bvar_est1 <- thin_posterior(bvar_est, thin = 15)
  summary(bvar_est1)
  
  
  
  ## Forecast error impulse response
  
  FEIR <- irf(bvar_est1, impulse = "myR", response = "myHP", n.ahead = 8, ci = 0.90)
  FEIR1 <- irf(bvar_est, impulse = "myR", response = "myGDP", n.ahead = 8, ci = 0.90)
  FEIR2 <- irf(bvar_est, impulse = "myGDP", response = "myHP", n.ahead = 8, ci = 0.90)
  FEIR3 <- irf(bvar_est, impulse = "myGDP", response = "myR", n.ahead = 8, ci = 0.90)
  FEIR4 <- irf(bvar_est, impulse = "myHP", response = "myGDP", n.ahead = 8, ci = 0.90)
  FEIR5 <- irf(bvar_est, impulse = "myHP", response = "myR", n.ahead = 8, ci = 0.90)
  
  par(mfrow=c(2,3))
  plot(FEIR,  
       main = "HP response from shock to R",
       col = 2, 
  )
  
  plot(FEIR1,  
       main = "GDP response from shock to R",
       col = 2, 
  )
  
  plot(FEIR2,  
       main = "HP response from shock to GDP",
       col = 2, 
  )
  
  plot(FEIR3,  
       main = "R response from shock to GDP",
       col = 2, 
  )
  
  plot(FEIR4,  
       main = "GDP response from shock to HP",
       col = 2, 
  )
  
  plot(FEIR5,  
       main = "R response from shock to HP",
       col = 2, 
  )
  
  
  ## Variance decomposition
  par(mfrow=c(1,1))
  decomp<- fevd(bvar_est, response = "myHP")
  
  plot(decomp, main = "Variance decomposition: House prices")
  

  
  
  
  
