library(quantmod)
library(urca)
library(vars)
library(tseries)
library(seasonal)
library(BVAR)
library(mfbvar)
library(dplyr)
library(ggridges)
library(ggplot2)
library(parallel)
library(vars)

if(!exists("foo", mode="function")) source("adf.R")

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
# seasonal adjust time series using X11. 
CPI1 <- final(seas(as.ts((CPI),freq=4)))

#3-Month or 90-day Rates and Yields: Treasury Securities for Sweden
#Percent,Not Seasonally Adjusted
getSymbols("IR3TTS01SEQ156N", src = "FRED")
INT <- ts(as.ts(IR3TTS01SEQ156N), start = c(1982, 1), frequency = 4)


#CLIF INDEX
CLIF_INDEX <- read.csv('CLIFS.csv', header = FALSE)
CLIF_INDEX <-CLIF_INDEX[nrow(CLIF_INDEX):1,]
CLIF_mts <- ts(CLIF_INDEX[,2], start = c(1970, 1), frequency = 12)
CLIF_qts = aggregate(CLIF_mts, nfrequency = 4)

t.test(window(CLIF_mts), start=c(2001, 1))
t.test(window(CLIF_qts), start=c(2001, 1))



#Data in levels
data <- ts.intersect(CLIF_qts, U, BNP, CPI1,  INT)
colnames(data) <-c("CLIF", "U", "BNP", 'CPI', "INT")
data

#Stationary data
sdata <- ts.intersect(CLIF_qts, diff(U), BNP, CPI1, INT)
colnames(sdata) <-c("CLIF", "DU", "BNP", 'CPI', "INT")
CLIF_qts


ADF.test(data[, 'BNP']) #Stationaty
ADF.test(data[, 'CLIF']) #Stationary
ADF.test(data[, 'U']) #Unitroot
ADF.test(data[, 'CPI']) #Stationary
ADF.test(data[, 'INT']) #Unitroot
summary(ur.kpss(data[,"INT"])) #Stationary

#We have to create a list with ts objects to be able to us mfbvar
#Start by setting the window,
list_data <- list(CLIF_mts, diff(U), BNP, CPI1, INT)
variables <- c("CLIF", "DU", "BNP", "CPI", "INT")

names(list_data) <- variables
list_data <- mapply(window, x = list_data,
                  start = list(c(2001, 1), c(2001, 1), c(2001, 1), c(2001, 1), c(2001, 1)))



prior_means <- list()
prior_means$clif <- t.test(sdata[, 'CLIF'])
prior_means$du <- t.test(sdata[, 'DU'])
prior_means$bnp <- t.test(sdata[, 'BNP'])
prior_means$cpi <- t.test(sdata[, 'CPI'])
prior_means$int <- t.test(sdata[, 'INT'])


#0.3664805 0.4673713
prior_means$clif
prio
#-0.02102656  0.13177813
prior_means$du
#1.386510 2.632993
prior_means$bnp
#0.2262972 0.4121828
prior_means$cpi
#0.9164065 1.6419063
prior_means$int 



prior <- set_prior(Y = list_data, n_lags = 5, n_reps = 1000, freq = c('m', 'q', 'q', 'q', 'q'))


prior_intervals <- matrix(c(0.1, 0.3, #CLIF INDEX
                               -1, 1, #Rate U
                               1, 3, #BNP
                               0, 1, #Inflation rate
                               0, 2), ncol=2, byrow=TRUE)


moments <- interval_to_moments(prior_intervals)


prior <- update_prior(prior,
                      d = "intercept",
                      prior_psi_mean = moments$prior_psi_mean,
                      prior_psi_Omega = moments$prior_psi_Omega)



prior <- update_prior(prior, n_fcst = 24)

summary(prior)


mod_ss_iw <- estimate_mfbvar(prior, prior = "ss", variance = "iw")
mod_ssng_iw <- estimate_mfbvar(prior, prior = "ssng", variance = "iw")

mod_ss_csv <- estimate_mfbvar(prior, prior = "ss", variance = "csv")
mod_ss_fsv <- estimate_mfbvar(prior, prior = "ss", variance = "fsv", n_fac = 1)

irf(mod_ss_iw, conf_bands, n_thin = 1L)


predict(mod_ssng_iw, pred_bands = 0.8)
plot(mod_ss_iw, nrow_facet = 3)
plot(mod_ssng_iw,  nrow_facet = 5)


pred_iw <- predict(mod_ss_iw, pred_bands = NULL)
pred_csv <- predict(mod_ss_csv, pred_bands = NULL)
pred_fsv <- predict(mod_ss_fsv, pred_bands = NULL)
pred_df <- bind_rows("Inverse Wishart" = pred_iw,
                        "Common stochastic volatility" = pred_csv,
                        "Factor stochastic volatility" = pred_fsv,
                        .id = "Variance") %>% filter(variable == "BNP")

ggplot(pred_df, aes(y = factor(fcst_date), x = fcst, fill = Variance)) +
  ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha = 0.5) +
  labs(x = "Swedish GDP Growth",
         y = "Date of Forecast") +
  coord_cartesian(xlim = c(-5, 15)) +
  theme_minimal() +
  scale_fill_brewer(palette = "YlGnBu")

const_vol <- median(sqrt(mod_ss_iw$Sigma[5, 5, ]))

varplot(mod_ss_fsv, variables = "BNP") +
  geom_hline(yintercept = const_vol ,
               color = "red", linetype = "dashed") +
  coord_cartesian(ylim = c(0, 20))

varplot(mod_ss_csv, variables = "BNP") +
  geom_hline(yintercept = const_vol,
               color = "red", linetype = "dashed") +
  coord_cartesian(ylim = c(0, 20))



par_fun <- function(lambda1, prior) {
  set.seed(2019)
  mod_par <- estimate_mfbvar(prior, prior = "ss", variance = "iw",
                               lambda1 = lambda1, lambda3 = 1)
  mdd(mod_par)
  }

cl <- makeCluster(4)
clusterEvalQ(cl, library("mfbvar"))
lambda1_seq <- seq(0.05, 1, by = 0.05)
result <- parSapply(cl, lambda1_seq, par_fun, prior = prior)
stopCluster(cl)

plot_df <- tibble(lambda1 = lambda1_seq, mdd = result)

ggplot(plot_df, aes(x = lambda1, y = mdd)) +
  geom_line() +
  geom_point(data = filter(plot_df, mdd == max(mdd))) +
  labs(y = "Marginal data density (log)",
         x = bquote(lambda[1])) +
  theme_minimal()

conf_bands(sdata)
