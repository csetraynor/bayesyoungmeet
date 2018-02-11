## create initial estimates

init_exp <- function() list(
  tau_s_raw = abs(rnorm(1)),
  tau_raw = abs(rnorm(M)),
  beta_raw = rnorm(M)
)
## Specify the variables for which you want history and density plots
parametersToPlot <- c("beta_raw")

## Additional variables to monitor
otherRVs <- c( "log_lik", "yhat_uncens")

parameters <- c(parametersToPlot, otherRVs)

################################################################################################
# run Stan

library(rstan)
nchains <- 1
niter <- 1000
nwarmup <- 500


fit_exp_clinical_bg <- stan(file = "exponential_clinical.stan", 
                            data = data, 
                            init = init_exp, 
                            chains = nchains, 
                            iter = niter,
                            warmup = nwarmup,
                            control = list(stepsize = 0.01, adapt_delta = 0.99),
                            cores = min(nchains, parallel::detectCores())) 

