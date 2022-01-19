library(urca)
library(vars)
library(tseries)
library(vars)


source('data_prep.R')
if(!exists("foo", mode="function")) source("adf.R")

# Estimate VAR
var_aic <- VAR(data, type = "const", lag.max = 6, ic = "AIC")

# Lag order suggested by AIC
var_aic$p


# Estimate
vec <- ca.jo(data, ecdet = "none", type = "trace",
             K = 4, spec = "transitory", season = 4)

summary(vec)

# Beta
round(vec@V, 2)

# Alpha, error correction terms
round(vec@W, 2)

#Convert to var to do IRF
vec_var <- vec2var(vec, r = 1)



calc_irf <- function(model, model_name, shock_var, resp_vars){
  for(i in 1:length(resp_vars)){
    ir <- irf(model, n.ahead = 20, impulse = shock_var, response = resp_vars[i],
        sign.constr=sign.constr, runs = 500, ortho=FALSE)
    plot(ir)
  }
  
}

shock_vars <- c('U', 'BNP', 'CPI', 'INT')
calc_irf(vec_var, 'vecm', 'CLIF', shock_vars)

# Obtain IRF
ir <- irf(vec_var, n.ahead = 20, impulse = "CLIF", response = "U",
          sign.constr=sign.constr, runs = 500, ortho = TRUE)

# Plot
plot(ir)


