library(bvartools)
source('stationarity_tests.R')
library(vars)

#GDP -> CPI -> DINT -> SFSI
md <- model_data.SFI[,c(2, 3, 4, 1)]

priordata<- window(md, start=c(1995, 4), end=c(2000,3)) # Subset of interest
mydata <- window(md, start=c(2000, 4)) # prior-data


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


bvar_est <- bvar(y = data$data$Y, x = data$data$Z, A = draws_a[1:48,],
                 C = draws_a[49:52, ], Sigma = draws_sigma)


bvar_est1 <- thin_posterior(bvar_est, thin = 15)

summary(bvar_est1)                                     


FEIR <- irf(bvar_est1, impulse = "SFSI", response = "BNP", n.ahead = 20, ci = 0.90, type = "oir")
FEIR1 <- irf(bvar_est1, impulse = "SFSI", response = "CPI", n.ahead = 20, ci = 0.90, type = "oir")
FEIR2 <- irf(bvar_est1, impulse = "SFSI", response = "DINT", n.ahead = 20, ci = 0.90, type = "oir")
FEIR3 <- irf(bvar_est1, impulse = "SFSI", response = "SFSI", n.ahead = 20, ci = 0.90, type = "oir")



par(mfrow=c(1,1))
plot(FEIR3,  
     main = "",
     col = 2, 
)
