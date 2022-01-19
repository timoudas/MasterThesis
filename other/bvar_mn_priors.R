library(urca)
library(vars)
library(tseries)
library(seasonal)
library(BVAR)
library(vars)

source('data_prep.R')

bv_lambda(mode = 0.2, sd = 0.4, min = 0.0001, max = 5)
bv_alpha(mode = 2, sd = 0.25, min = 1, max = 3)
bv_psi(scale = 0.004, shape = 0.004, mode = "auto", min = "auto", max = "auto")

bv_mn(
  lambda = bv_lambda(),
  alpha = bv_alpha(),
  psi = bv_psi(),
  var = 10000000,
  b = 1
)


bvar_est <- bvar(data = sdata, 
          lags = 5,
          n_draw = 10000L,
          n_burn = 5000L,
          n_thin = 1L,
          priors = bv_priors(hyper='auto', mn=bv_mn()),
          mh = bv_mh(),
          fcast = NULL,
          irf = NULL,
          verbose = TRUE)

OIR <- irf(bvar_est, horizon=20)

plot(OIR, main = "Orthogonalised Impulse Response", xlab = "Period", ylab = "Response", )
